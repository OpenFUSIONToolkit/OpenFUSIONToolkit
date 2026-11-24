!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file xmhd_2d.F90
!
!> Solve time-dependent MHD equations in 2D with lagrange basis set
!!
!! Solves the single-temperature resistive MHD equations (continuity, momentum, temperature, and
!! induction for the poloidal flux and out-of-plane field), or an incompressible variant, in planar
!! or axisymmetric geometry using Lagrange elements for all fields. The full system of equations is
!! defined by the nonlinear residual (@ref xmhd_2d::nlfun_apply "nlfun_apply") and is documented in
!! \ref doc_mug_xmhd_2d_eqns.
!---------------------------------------------------------------------------
MODULE xmhd_2d
#if !defined(XMHD_RST_LEN)
#define XMHD_RST_LEN 5
#endif
#define XMHD_ITCACHE 10
USE oft_base
USE oft_io, ONLY: hdf5_read, hdf5_write, oft_file_exist, &
  hdf5_field_exist, oft_bin_file, xdmf_plot_file
USE oft_quadrature
USE oft_mesh_type, ONLY: oft_bmesh, cell_is_curved
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
!
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_local_mat, oft_vector_ptr, &
  vector_extrapolate, oft_graph, oft_graph_ptr, OFT_ACCESS_ACCUMULATE
USE oft_solver_utils, ONLY: create_solver_xml, create_diag_pre
USE oft_deriv_matrices, ONLY: oft_noop_matrix, oft_mf_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre, create_native_pre
!
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_fem_comp_type
USE fem_utils, ONLY: fem_dirichlet_diag, fem_dirichlet_vec, bfem_map_flag
USE oft_lag_basis, ONLY: oft_lag_setup,oft_scalar_bfem, oft_blag_eval, oft_blag_geval, oft_2D_lagrange_cast, &
  oft_blag_npos
USE oft_blag_operators, ONLY: oft_blag_vproject,oft_blag_project, oft_blag_getmop, oft_lag_bginterp
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: mu0, elec_charge, proton_mass
USE oft_gs, ONLY: gs_epsilon
USE oft_la_utils, ONLY: create_matrix
IMPLICIT NONE
#include "local.h"
#if !defined(TDIFF_RST_LEN)
#define TDIFF_RST_LEN 5
#endif
PRIVATE
REAL(r8), PARAMETER :: k_boltz = elec_charge !< Boltzmann constant (elec_charge, so temperatures are in eV)

!------------------------------------------------------------------------------
!> Object to compute left-hand side of system for nonlinear solves
!------------------------------------------------------------------------------
TYPE, extends(oft_noop_matrix) :: xmhd_2d_nlfun
  REAL(r8) :: dt = -1.d0 !< time step
  REAL(r8) :: diag_vals(5) = 0.d0 !< diagnostic values
  TYPE(oft_xmhd_2d_sim), pointer :: parent_sim => NULL() !< pointer to parent simulation object for access to parameters
  CONTAINS
  !> Apply the matrix
  PROCEDURE :: apply_real => nlfun_apply
END TYPE xmhd_2d_nlfun

!------------------------------------------------------------------------------
!> Object to compute left-hand side of system for linear solves
!------------------------------------------------------------------------------
TYPE, extends(oft_noop_matrix) :: xmhd_2d_mfun
  REAL (r8) :: diag_vals(5) = 0.d0 !< diagnostic values
  TYPE(oft_xmhd_2d_sim), pointer :: parent_sim => NULL() !< pointer to parent simulation object for access to parameters
  CONTAINS
  !> Apply the matrix
  PROCEDURE :: apply_real => mfun_apply
END TYPE xmhd_2d_mfun
!------------------------------------------------------------------------------
!> Fields, geometry, and material properties at a quadrature point, shared by
!! the residual (@ref nlfun_apply), the mass function (@ref mfun_apply), and the
!! approximate Jacobian
!!
!! Geometric factors are set so that each term can be written once for
!! both geometries: `rmul` is the volume weighting (R in cylindrical geometry)
!! and `rden` is the regularized radius (R + @ref oft_gs::gs_epsilon) used where
!! terms are divided by R. Both are one in Cartesian geometry.
!------------------------------------------------------------------------------
TYPE :: xmhd_2d_qpt
  REAL(r8) :: n !< Density
  REAL(r8) :: dn(3) !< Density gradient
  REAL(r8) :: vel(3) !< Velocity
  REAL(r8) :: dvel(3,3) !< Velocity gradient (dvel(i,j) = d vel_i / d x_j)
  REAL(r8) :: div_vel !< Velocity divergence
  REAL(r8) :: T !< Temperature (pressure if incompressible)
  REAL(r8) :: dT(3) !< Temperature (pressure if incompressible) gradient
  REAL(r8) :: psi !< Poloidal flux, including the vacuum offset
  REAL(r8) :: dpsi(3) !< Poloidal flux gradient, including the vacuum offset
  REAL(r8) :: by !< Out-of-plane magnetic field component
  REAL(r8) :: dby(3) !< Out-of-plane magnetic field gradient
  REAL(r8) :: btmp(3) !< Magnetic field
  REAL(r8) :: rmul !< Volume weighting factor (R in cylindrical geometry, otherwise 1)
  REAL(r8) :: rden !< Radial divisor (R + gs_epsilon in cylindrical geometry, otherwise 1)
  REAL(r8) :: int_factor !< Quadrature weight times cell Jacobian determinant
  REAL(r8) :: chi !< Thermal diffusivity
  REAL(r8) :: m_i !< Ion mass
  REAL(r8) :: eta(2) !< Electrical resistivity (in plane, out of plane)
  REAL(r8) :: nu !< Kinematic viscosity
  REAL(r8) :: gamma !< Adiabatic index
  REAL(r8) :: D_diff !< Density diffusivity
END TYPE xmhd_2d_qpt
!------------------------------------------------------------------------------
!> Main 2D MHD simulation object
!------------------------------------------------------------------------------
TYPE, public :: oft_xmhd_2d_sim
  LOGICAL :: mfnk = .FALSE. !< Use matrix free method?
  LOGICAL :: cyl_flag = .FALSE. !Use cylindrical coordinates
  LOGICAL :: incomp = .FALSE. !Run incompressible MHD?
  LOGICAL :: linear = .FALSE. !Run sim. in linear mode
  LOGICAL :: timestep_cn = .FALSE. ! Use Crank-Nicolson timestep
  INTEGER(i4) :: nsteps = -1 !< # of timesteps
  INTEGER(i4) :: rst_base = 0 !< starting index of restart files
  INTEGER(i4) :: rst_freq = 10 !< frequency to generate restart files
  INTEGER(i4) :: pre_freq = 1 !< frequency to update preconditioner
  INTEGER(i4) :: ittarget = 40 !< Needs docs
  REAL(r8) :: dt = -1.d0 !< timestep
  REAL(r8) :: jac_dt = -1.d0 !< timestep for computing approximate jacobian
  REAL(r8) :: t = 0.d0 !< current time
  REAL(r8), ALLOCATABLE :: chi(:) !< thermal diffusivity, by region
  REAL(r8), ALLOCATABLE :: nu(:) !< kinematic viscosity [m^2/s], by region
  REAL(r8), ALLOCATABLE :: gamma(:) !< adiabatic index, by region, by region
  REAL(r8), ALLOCATABLE :: D_diff(:) !< diffusivity, by region
  REAL(r8) :: den_scale = 1.d19 !< scale factor used to normalize density
  REAL(r8), ALLOCATABLE :: m_i(:) !< ion mass, by region
  REAL(r8), ALLOCATABLE :: eta(:, :) !< electrical resistivity, in units of mu0, by region. First column is in plane, second column is out of plane (for anisotropic resistivity)
  LOGICAL, ALLOCATABLE :: ignore_rmask(:) !< Mask for regions to ignore when populating physics (defaults to false)
  REAL(r8) :: lin_tol = 1.d-8 !< absolute tolerance for linear solver
  REAL(r8) :: nl_tol = 1.d-5 !< absolute tolerance for nonlinear solver
  !> Uniform, static background magnetic field (only a Z component is allowed in cylindrical geometry)
  !!
  !! @note On a non-periodic domain this is equivalent to adding the flux \f$ B_{0,z} x - B_{0,x} z \f$
  !! (planar) or \f$ B_{0,Z} R^2/2 \f$ (cylindrical) to `psi_vac`, and a uniform \f$ B_{0,y} \f$ to the
  !! initial out-of-plane field. On a periodic domain a uniform in-plane field cannot be represented
  !! by `psi_vac`, because its flux is not periodic, so `B_0` is required there.
  REAL (r8) :: B_0(3) = 0.d0
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: n_bc => NULL() !< n BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velx_bc => NULL() !< vel BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: vely_bc => NULL() !< vel BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velz_bc => NULL() !< vel BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: T_bc => NULL() !< T BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: psi_bc => NULL() !< psi BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: by_bc => NULL() !< by BC flag
  TYPE(oft_fem_comp_type), POINTER :: fe_rep => NULL() !< Finite element representation for solution field
  TYPE(xdmf_plot_file) :: xdmf_plot
  CLASS(oft_vector), POINTER :: u => NULL() !< current solution vector
  CLASS(oft_vector), POINTER :: psi_vac => NULL() !< vacuum poloidal flux offset added to the psi field (matches psi/oft_blagrange space)
  CLASS(oft_vector), POINTER :: u0 => NULL() !< equilibrium solution vector (used for linearization if linear = True)
  CLASS(oft_matrix), POINTER :: jacobian => NULL() !< approximate jacobian matrix
  TYPE(oft_mf_matrix), POINTER :: mf_mat => NULL() !< Matrix free operator
  TYPE(xmhd_2d_nlfun), POINTER :: nlfun => NULL() !< Nonlinear function
  TYPE(xmhd_2d_mfun), POINTER :: mfun => NULL() !< mass matrix function
  TYPE(xml_node):: xml_root !< XML root element
  TYPE(xml_node):: xml_pre_def !< XML element for preconditioner definition
  contains
  !> Setup
  PROCEDURE :: setup => setup
  !> Set default boundary conditions
  PROCEDURE :: setup_bc => setup_bc
  !> Run simulation
  PROCEDURE :: run_simulation => run_simulation
  !> Run simulation in linear mode
  PROCEDURE :: run_lin_simulation => run_lin_simulation
  !> Save restart file
  PROCEDURE :: rst_save => rst_save
  !> Load restart file
  PROCEDURE :: rst_load => rst_load
END TYPE oft_xmhd_2d_sim
TYPE(oft_xmhd_2d_sim), POINTER :: current_sim => NULL() !Active 2D MHD simulation
!------------------------------------------------------------------------------
!> Interface for approximate Jacobian construction
!------------------------------------------------------------------------------
ABSTRACT INTERFACE
  SUBROUTINE xmhd_2d_jac_builder(self,a)
    IMPORT oft_xmhd_2d_sim, oft_vector
    CLASS(oft_xmhd_2d_sim), INTENT(INOUT) :: self !< Simulation object
    CLASS(oft_vector), INTENT(INOUT) :: a !< Solution for computing Jacobian
  END SUBROUTINE xmhd_2d_jac_builder
END INTERFACE
!> Jacobian builder for the active simulation (compressible or incompressible), set by the run drivers
PROCEDURE(xmhd_2d_jac_builder), POINTER :: build_jacobian => NULL()
!
CLASS(multigrid_mesh), POINTER :: mg_mesh => NULL() !< Multigrid mesh object
CLASS(oft_bmesh), POINTER, PUBLIC :: mesh => NULL() !< Current mesh level
TYPE(oft_ml_fem_type), TARGET, PUBLIC :: ML_oft_blagrange !< Multilevel finite element representation
CLASS(oft_scalar_bfem), POINTER, PUBLIC :: oft_blagrange => NULL() !< Lagrange finite element representation
TYPE(oft_ml_fem_type), TARGET, PUBLIC :: ML_oft_blagrange_p !< Multilevel finite element representation for pressure (if incompressible)
CLASS(oft_scalar_bfem), POINTER, PUBLIC :: oft_blagrange_p => NULL() !< Lagrange finite element representation for pressure (if incompressible)
PUBLIC xmhd_2d_plot, build_approx_jacobian_comp, build_approx_jacobian_incomp
CONTAINS

!---------------------------------------------------------------------------
!Main driver for 2D MHD nonlinear time advance
!---------------------------------------------------------------------------
subroutine run_simulation(self)
class(oft_xmhd_2d_sim), target, intent(inout) :: self
type(oft_nksolver) :: nksolver
!---Solver objects
class(oft_vector), pointer :: u,v,up
type(oft_native_gmres_solver), target :: solver
type(oft_timer) :: mytimer
!---History file fields
TYPE(oft_bin_file) :: hist_file
integer(i4) :: hist_i4(3)
real(4) :: hist_r4(6)
!---Timestep control
integer(i4) :: itcount(XMHD_ITCACHE)
real(r8) :: dtin,dthist(XMHD_ITCACHE)
!---
!---Extrapolation fields
integer(i4), parameter :: maxextrap=2
integer(i4) :: nextrap
real(r8), allocatable, dimension(:) :: extrapt
type(oft_vector_ptr), allocatable, dimension(:) :: extrap_fields
!---
character(LEN=TDIFF_RST_LEN) :: rst_char
integer(i4) :: i,j,io_stat,rst_tmp,npre
real(r8) :: n_avg, u_avg(3), T_avg, psi_avg, by_avg,elapsed_time
real(r8), pointer :: plot_vals(:),plot_vec(:,:)
current_sim=>self
CALL check_settings(self)
!---------------------------------------------------------------------------
! Create solver fields
!---------------------------------------------------------------------------
call self%fe_rep%vec_create(u)
call self%fe_rep%vec_create(up)
call self%fe_rep%vec_create(v)
!---

ALLOCATE(extrap_fields(maxextrap),extrapt(maxextrap))
DO i=1,maxextrap
  CALL self%fe_rep%vec_create(extrap_fields(i)%f)
  extrapt(i)=0.d0
END DO
nextrap=0
!---Set initial time and load initial solution vector from driver
self%t=0.d0
CALL u%add(0.d0,1.d0,self%u)
!---Create initial conditions restart file
104 FORMAT (I TDIFF_RST_LEN.TDIFF_RST_LEN)
WRITE(rst_char,104)0
CALL self%rst_save(u, self%t, self%dt, 'xmhd2d_'//rst_char//'.rst', 'U')
CALL save_background(self, 'xmhd2d_'//rst_char//'.rst')

!---Fill timestep control cache
itcount=self%ittarget
dthist=self%dt
dtin=self%dt
!---------------------------------------------------------------------------
! Setup linear solver
!---------------------------------------------------------------------------
IF(self%incomp)THEN
  build_jacobian=>build_approx_jacobian_incomp
ELSE
  build_jacobian=>build_approx_jacobian_comp
END IF
self%jac_dt=self%dt
IF(self%timestep_cn)self%jac_dt=self%dt/2.d0
CALL build_jacobian(self,u)
IF(self%mfnk)THEN
  ALLOCATE(self%mf_mat)
  CALL up%set(1.d0)
  CALL self%mf_mat%setup(u,self%nlfun,up)
  self%mf_mat%b0=1.d-5
  solver%A=>self%mf_mat
ELSE
  solver%A=>self%jacobian
END IF
solver%its=self%ittarget*4
solver%atol=self%lin_tol
solver%itplot=1
solver%nrits=40
solver%pm=oft_env%pm
NULLIFY(solver%pre)
IF(self%xml_pre_def%associated())THEN
  CALL create_solver_xml(solver%pre,self%xml_pre_def)
ELSE IF(self%incomp)THEN
  ! If incompressible, default to a complete LU preconditioner
  CALL create_native_pre(solver%pre,"lu")
ELSE
  CALL create_diag_pre(solver%pre)
END IF
solver%pre%A=>self%jacobian

!---------------------------------------------------------------------------
! Setup nonlinear solver
!---------------------------------------------------------------------------
nksolver%A=>self%nlfun
nksolver%J_inv=>solver
nksolver%its=30
nksolver%atol=self%nl_tol
nksolver%rtol=1.d-20 ! Disable relative tolerance
IF(self%mfnk)THEN
  nksolver%J_update=>mfnk_update
  nksolver%up_freq=1
ELSE
  nksolver%J_update=>update_jacobian
  nksolver%up_freq=4
END IF

!---------------------------------------------------------------------------
! Setup history file
!---------------------------------------------------------------------------
IF(oft_env%head_proc)THEN
  CALL hist_file%setup('xmhd_2d.hist', desc="History file for non-linear 2D xMHD run")
  CALL hist_file%add_field('ts',   'i4', desc="Time step index")
  CALL hist_file%add_field('lits', 'i4', desc="Linear iteration count")
  CALL hist_file%add_field('nlits','i4', desc="Non-linear iteration count")
  CALL hist_file%add_field('time', 'r4', desc="Simulation time [s]")
  CALL hist_file%add_field('men',  'r4', desc="Magnetic energy [J*(2*mu0)]")
  CALL hist_file%add_field('ven',  'r4', desc="Kinetic energy [J]")
  CALL hist_file%add_field('ne',   'r4', desc="Average density [m^-3]")
  CALL hist_file%add_field('ti',   'r4', desc="Average temperature [eV]")
  CALL hist_file%add_field('stime','r4', desc="Walltime [s]")
  CALL hist_file%write_header
  CALL hist_file%open ! Open history file
END IF

!---------------------------------------------------------------------------
! Begin time stepping
!---------------------------------------------------------------------------
npre=0
DO i=1,self%nsteps
  IF(oft_env%head_proc)CALL mytimer%tick()
  IF(self%timestep_cn)THEN
    self%nlfun%dt=-self%dt/2.0
  ELSE
    self%nlfun%dt=0.d0
  END IF
  CALL self%nlfun%apply(u,v)
  IF(self%timestep_cn)THEN
    self%nlfun%dt=self%dt/2.0
    self%jac_dt=self%dt/2.d0
  ELSE
    self%nlfun%dt=self%dt
    self%jac_dt=self%dt
  END IF
  npre = npre + 1
  IF((.NOT.self%mfnk).OR.MOD(npre,self%pre_freq)==0)THEN
    CALL update_jacobian(u)
    CALL solver%pre%update(.TRUE.)
  END IF
  !
  DO j=maxextrap,2,-1
    CALL extrap_fields(j)%f%add(0.d0,1.d0,extrap_fields(j-1)%f)
    extrapt(j)=extrapt(j-1)
  END DO
  IF(nextrap<maxextrap)nextrap=nextrap+1
  CALL extrap_fields(1)%f%add(0.d0,1.d0,u)
  extrapt(1)=self%t
  IF(i>maxextrap)CALL vector_extrapolate(extrapt,extrap_fields,nextrap,self%t+self%dt,u)
  CALL nksolver%apply(u,v)
  IF(nksolver%cits<0)CALL oft_abort("Nonlinear solve failed","run_simulation",__FILE__)
  !---------------------------------------------------------------------------
  ! Write out initial solution progress
  !---------------------------------------------------------------------------
  n_avg=self%nlfun%diag_vals(3)/self%nlfun%diag_vals(5)
  T_avg = self%nlfun%diag_vals(4)/ self%nlfun%diag_vals(5)
  IF(oft_env%head_proc)THEN
    elapsed_time=mytimer%tock()
    hist_i4=[self%rst_base+i-1,nksolver%lits,nksolver%nlits]
    hist_r4=REAL([self%t,self%nlfun%diag_vals(1), self%nlfun%diag_vals(2), n_avg, T_avg,elapsed_time],4)
103 FORMAT(' Timestep',I8,ES14.6,2X,I4,2X,I4,F12.3,ES12.2)
    WRITE(*,103)self%rst_base+i,self%t+self%dt,nksolver%lits,nksolver%nlits,elapsed_time,self%dt
    IF(oft_debug_print(1))WRITE(*,*)
    CALL hist_file%write(data_i4=hist_i4, data_r4=hist_r4)
  END IF
  !---------------------------------------------------------------------------
  ! Update timestep and save solution
  !---------------------------------------------------------------------------
  self%t=self%t+self%dt
  CALL self%u%add(0.d0,1.d0,u)
  IF(MOD(i,self%rst_freq)==0)THEN
    IF(oft_env%head_proc)CALL mytimer%tick
    !---Create restart file
    WRITE(rst_char,104)self%rst_base+i
    READ(rst_char,104,IOSTAT=io_stat)rst_tmp
    IF((io_stat/=0).OR.(rst_tmp/=self%rst_base+i))CALL oft_abort("Step count exceeds format width", "run_simulation", __FILE__)
    CALL self%rst_save(u, self%t, self%dt, 'xmhd2d_'//rst_char//'.rst', 'U')
    IF(oft_env%head_proc)THEN
      elapsed_time=mytimer%tock()
      WRITE(*,'(2X,A,F12.3)')'I/O Time = ',elapsed_time
      CALL hist_file%flush
    END IF
    !---
  END IF
  !---
  itcount(MOD(i,XMHD_ITCACHE)+1)=nksolver%lits
  dthist(MOD(i,XMHD_ITCACHE)+1)=self%dt
  self%dt=self%ittarget*SUM(dthist)/SUM(itcount)
  IF(self%dt>dtin)self%dt=dtin
  IF(self%dt<1.d-3*dtin)CALL oft_abort('Time step dropped too low!','run_simulation',__FILE__)
  IF(ABS(1.d0-self%dt/dthist(MOD(i,XMHD_ITCACHE)+1))>0.1d0)npre=-1
END DO
CALL hist_file%close()
CALL nksolver%delete()
CALL solver%delete()
CALL u%delete()
CALL up%delete()
CALL v%delete()
DEALLOCATE(u,up,v)
!---extrapolation fields
IF(maxextrap>0)THEN
  DO i=1,maxextrap
    CALL extrap_fields(i)%f%delete
    DEALLOCATE(extrap_fields(i)%f)
  END DO
  DEALLOCATE(extrap_fields,extrapt)
END IF
end subroutine run_simulation
!---------------------------------------------------------------------------
!> Main driver for 2D MHD linear time advance
!---------------------------------------------------------------------------
SUBROUTINE run_lin_simulation(self)
CLASS(oft_xmhd_2d_sim), TARGET, INTENT(INOUT) :: self
!---Solver objects
CLASS(oft_vector), POINTER :: u, v, up,un1,un2
TYPE(oft_native_gmres_solver), TARGET :: solver
TYPE(OFT_TIMER) :: mytimer
!---History file fields
TYPE(oft_bin_file) :: hist_file
INTEGER(i4) :: hist_i4(3)
REAL(4) :: hist_r4(6)
!---Extrapolation fields
integer(i4), parameter :: maxextrap=2
integer(i4) :: nextrap
real(r8), allocatable, dimension(:) :: extrapt
type(oft_vector_ptr), allocatable, dimension(:) :: extrap_fields
!---
character(LEN=TDIFF_RST_LEN) :: rst_char
integer(i4) :: i,j,io_stat,rst_tmp,npre
real(r8) :: n_avg, u_avg(3), T_avg, psi_avg, by_avg,elapsed_time,diag_save(5)
real(r8), pointer :: plot_vals(:),plot_vec(:,:)
current_sim=>self
!---------------------------------------------------------------------------
! Set boundary conditions not alreadys set
!---------------------------------------------------------------------------
CALL self%setup_bc()
!---------------------------------------------------------------------------
! Create solver fields
!---------------------------------------------------------------------------
IF(self%cyl_flag)CALL oft_abort( &
        "Cylindrical mode not currently supported for linear solve","run_lin_simulation",__FILE__)
CALL check_settings(self)

call self%fe_rep%vec_create(u)
call self%fe_rep%vec_create(up)
call self%fe_rep%vec_create(v)
!---

ALLOCATE(extrap_fields(maxextrap),extrapt(maxextrap))
DO i=1,maxextrap
  CALL self%fe_rep%vec_create(extrap_fields(i)%f)
  extrapt(i)=0.d0
END DO
nextrap=0

self%t=0.d0
CALL u%add(0.d0,1.d0,self%u)
!---Create initial conditions restart file
104 FORMAT (I TDIFF_RST_LEN.TDIFF_RST_LEN)
WRITE(rst_char,104)0
CALL self%rst_save(self%u0, self%t, self%dt, 'xmhd2d_eq.rst', 'U')
CALL self%rst_save(u, self%t, self%dt, 'xmhd2d_'//rst_char//'.rst', 'U')
CALL save_background(self, 'xmhd2d_'//rst_char//'.rst')

ALLOCATE(self%mfun)
self%mfun%parent_sim=>self

! Construct the linear advance matrix with equilibrium fields
IF(self%incomp)THEN
  build_jacobian=>build_approx_jacobian_incomp
ELSE
  build_jacobian=>build_approx_jacobian_comp
END IF
self%jac_dt=self%dt
CALL build_jacobian(self,self%u0)
!---------------------------------------------------------------------------
! Setup linear solver
!---------------------------------------------------------------------------
solver%A=>self%jacobian
solver%its=400
solver%atol=self%lin_tol
solver%itplot=1
solver%nrits=20
solver%pm=oft_env%pm
NULLIFY(solver%pre)
IF(self%xml_pre_def%associated())THEN
  CALL create_solver_xml(solver%pre,self%xml_pre_def)
ELSE IF(self%incomp)THEN
  ! If incompressible, default to a complete LU preconditioner
  CALL create_native_pre(solver%pre,"lu")
ELSE
  CALL create_diag_pre(solver%pre)
END IF
solver%pre%A=>self%jacobian

!---------------------------------------------------------------------------
! Setup history file
!---------------------------------------------------------------------------
IF(oft_env%head_proc)THEN
  CALL hist_file%setup('xmhd_2d.hist', desc="History file for linear 2D xMHD run")
  CALL hist_file%add_field('ts',   'i4', desc="Time step index")
  CALL hist_file%add_field('lits', 'i4', desc="Linear iteration count")
  CALL hist_file%add_field('nlits','i4', desc="Non-linear iteration count")
  CALL hist_file%add_field('time', 'r4', desc="Simulation time [s]")
  CALL hist_file%add_field('men',  'r4', desc="Magnetic energy [J*(2*mu0)]")
  CALL hist_file%add_field('ven',  'r4', desc="Kinetic energy [J]")
  CALL hist_file%add_field('ne',   'r4', desc="Average density [m^-3]")
  CALL hist_file%add_field('ti',   'r4', desc="Average temperature [eV]")
  CALL hist_file%add_field('stime','r4', desc="Walltime [s]")
  CALL hist_file%write_header
  CALL hist_file%open ! Open history file
END IF
!------------------------------------------------------------------------------
! Begin time stepping
!------------------------------------------------------------------------------
IF(self%timestep_cn)THEN
  CALL self%fe_rep%vec_create(un1)
  CALL self%fe_rep%vec_create(un2)
  CALL un2%add(0.d0,1.d0,u)
  CALL un1%add(0.d0,1.d0,u)
END IF
DO i=1,self%nsteps
  IF(oft_env%head_proc)CALL mytimer%tick
  DO j=maxextrap,2,-1
    CALL extrap_fields(j)%f%add(0.d0,1.d0,extrap_fields(j-1)%f)
    extrapt(j)=extrapt(j-1)
  END DO
  IF(nextrap<maxextrap)nextrap=nextrap+1
  CALL extrap_fields(1)%f%add(0.d0,1.d0,u)
  extrapt(1)=self%t
  IF(self%timestep_cn)THEN
    CALL un2%add(0.d0,1.d0,un1)
    CALL un1%add(0.d0,1.d0,u)
  END IF
  IF(self%timestep_cn.AND.(i>1))THEN
    IF(i==2)THEN
      self%jac_dt=self%dt*2.d0/3.d0
      CALL build_jacobian(self,self%u0)
      CALL solver%pre%update(.TRUE.)
    END IF
    CALL un2%add(-1.d0,4.d0,un1)
    CALL un2%scale(1.d0/3.d0)
    CALL self%mfun%apply(u,v) ! This double call and diag swap should be replaced with a better approach
    diag_save=self%mfun%diag_vals
    CALL self%mfun%apply(un2,v)
    self%mfun%diag_vals=diag_save
  ELSE
    CALL self%mfun%apply(u,v)
  END IF
  IF(i>maxextrap)CALL vector_extrapolate(extrapt,extrap_fields,nextrap,self%t+self%dt,u)
  CALL solver%apply(u,v)
  !---------------------------------------------------------------------------
  ! Write out initial solution progress
  !---------------------------------------------------------------------------
  n_avg=self%mfun%diag_vals(3)/self%mfun%diag_vals(5)
  T_avg = self%mfun%diag_vals(4)/ self%mfun%diag_vals(5)
  IF(oft_env%head_proc)THEN
    elapsed_time=mytimer%tock()
    hist_i4=(/self%rst_base+i-1,solver%cits,0/)
    hist_r4=REAL([self%t,self%mfun%diag_vals(1), self%mfun%diag_vals(2), n_avg, T_avg,elapsed_time],4)
103 FORMAT(' Timestep',I8,ES14.6,2X,I4,F12.3,ES12.2)
    WRITE(*,103)self%rst_base+i,self%t,solver%cits,elapsed_time,self%dt
    IF(oft_debug_print(1))WRITE(*,*)
    CALL hist_file%write(data_i4=hist_i4, data_r4=hist_r4)
  END IF
  !---------------------------------------------------------------------------
  ! Update timestep and save solution
  !---------------------------------------------------------------------------
  self%t=self%t+self%dt
  CALL self%u%add(0.d0,1.d0,u)
  IF(MOD(i,self%rst_freq)==0)THEN
    IF(oft_env%head_proc)CALL mytimer%tick
    !---Create restart file
    WRITE(rst_char,104)self%rst_base+i
    READ(rst_char,104,IOSTAT=io_stat)rst_tmp
    IF((io_stat/=0).OR.(rst_tmp/=self%rst_base+i))CALL oft_abort("Step count exceeds format width", "run_simulation", __FILE__)
    CALL self%rst_save(u, self%t, self%dt, 'xmhd2d_'//rst_char//'.rst', 'U')
    IF(oft_env%head_proc)THEN
      elapsed_time=mytimer%tock()
      WRITE(*,'(2X,A,F12.3)')'I/O Time = ',elapsed_time
      CALL hist_file%flush
    END IF
    !---
  END IF
END DO
CALL hist_file%close()
CALL solver%delete()
CALL u%delete()
CALL up%delete()
CALL v%delete()
DEALLOCATE(u,up,v)
IF(self%timestep_cn)THEN
  CALL un1%delete()
  CALL un2%delete()
  DEALLOCATE(un1,un2)
END IF
!---extrapolation fields
IF(maxextrap>0)THEN
  DO i=1,maxextrap
    CALL extrap_fields(i)%f%delete
    DEALLOCATE(extrap_fields(i)%f)
  END DO
  DEALLOCATE(extrap_fields,extrapt)
END IF
END SUBROUTINE run_lin_simulation
!---------------------------------------------------------------------------
!> Check that required material properties are set and valid
!!
!! Material properties default to -1 and must be set for every region that is not
!! ignored (`ignore_rmask`). Properties that are not used by the chosen model
!! (eg. \f$ \chi \f$ and \f$ \gamma \f$ when incompressible) are not checked.
!---------------------------------------------------------------------------
SUBROUTINE check_settings(self)
CLASS(oft_xmhd_2d_sim), INTENT(in) :: self !< Simulation object
INTEGER(i4) :: i
DO i=1,SIZE(self%ignore_rmask)
  IF(self%ignore_rmask(i))CYCLE
  IF(self%m_i(i)<=0.d0)CALL oft_abort("Ion mass must be positive","check_settings",__FILE__)
  IF(self%nu(i)<0.d0)CALL oft_abort("Viscosity not set or negative","check_settings",__FILE__)
  IF(self%eta(i,1)<0.d0)CALL oft_abort("In-plane resistivity not set or negative","check_settings",__FILE__)
  IF(self%eta(i,2)<=0.d0)CALL oft_abort("Out-of-plane resistivity must be positive","check_settings",__FILE__)
  IF(.NOT.self%incomp)THEN
    IF(self%D_diff(i)<0.d0)CALL oft_abort("Particle diffusivity not set or negative","check_settings",__FILE__)
    IF(self%chi(i)<0.d0)CALL oft_abort("Thermal diffusivity not set or negative","check_settings",__FILE__)
    IF(self%gamma(i)<=1.d0)CALL oft_abort("Adiabatic index must be greater than 1","check_settings",__FILE__)
  END IF
END DO
END SUBROUTINE check_settings
!---------------------------------------------------------------------------
!> Compute the mass matrix for RHS
!---------------------------------------------------------------------------
SUBROUTINE mfun_apply(self,a,b)
class(xmhd_2d_mfun), intent(inout) :: self !< mass function object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
type(oft_quad_type), pointer :: quad
LOGICAL :: incomp
INTEGER(i4) :: i
REAL(r8) :: diag_vals(5)
REAL(r8), POINTER, DIMENSION(:) :: n_weights,T_weights,psi_weights,by_weights,psi_vac_weights
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: T_res,n_res,psi_res,by_res,velx_res,vely_res,velz_res
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights
quad=>oft_blagrange%quad
NULLIFY(n_res, velx_res, vely_res, velz_res, T_res, psi_res, by_res)

!---Get weights from solution vector
CALL get_field_weights(self%parent_sim,a,n_weights,vel_weights,T_weights,psi_weights,by_weights,psi_vac_weights)
incomp = self%parent_sim%incomp

!---Zero result and get storage array
CALL b%set(0.d0)
CALL b%borrow_local(n_res,OFT_ACCESS_ACCUMULATE,1)
CALL b%borrow_local(velx_res,OFT_ACCESS_ACCUMULATE,2)
CALL b%borrow_local(vely_res,OFT_ACCESS_ACCUMULATE,3)
CALL b%borrow_local(velz_res,OFT_ACCESS_ACCUMULATE,4)
CALL b%borrow_local(T_res,OFT_ACCESS_ACCUMULATE,5)
CALL b%borrow_local(psi_res,OFT_ACCESS_ACCUMULATE,6)
CALL b%borrow_local(by_res,OFT_ACCESS_ACCUMULATE,7)
diag_vals=0.d0

BLOCK
INTEGER(i4) :: m, jr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs
LOGICAL :: curved
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals,T_loc,n_loc,psi_loc,by_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: vel_loc,basis_grads,res_loc
REAL(r8) :: jac_det,jac_mat(3,4),coords(3),int_factor,rmul,rden
TYPE(xmhd_2d_qpt) :: pt
!$omp parallel private(m,jr,cell_dofs,curved,basis_vals,T_loc,n_loc,psi_loc,by_loc,vel_loc, &
!$omp basis_grads,res_loc,jac_det,jac_mat,coords,int_factor,rmul,rden,pt) reduction(+:diag_vals)
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(T_loc(oft_blagrange%nce),n_loc(oft_blagrange%nce),psi_loc(oft_blagrange%nce), &
  by_loc(oft_blagrange%nce),vel_loc(3,oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce),res_loc(oft_blagrange%nce,7))
!$omp do ordered
DO i=1,mesh%nc
  IF(self%parent_sim%ignore_rmask(mesh%reg(i)))CYCLE ! Skip cells in ignored regions
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  res_loc = 0.d0 ! Zero local (cell) contribution to function
  n_loc = n_weights(cell_dofs)
  vel_loc = vel_weights(:, cell_dofs)
  !---Field 5 is pressure when incompressible -> no RHS term
  !-- Set weights to 0 so it is not reconstructed
  IF(incomp)THEN
    T_loc = 0.d0
  ELSE
    T_loc = T_weights(cell_dofs)
  END IF
  !---The vacuum flux (psi_vac) is part of the equilibrium only, so it is not added to the perturbation
  psi_loc = psi_weights(cell_dofs)
  by_loc = by_weights(cell_dofs)
  !---Set material properties
  pt%m_i = self%parent_sim%m_i(mesh%reg(i))
  pt%eta = self%parent_sim%eta(mesh%reg(i), :)
  pt%gamma = self%parent_sim%gamma(mesh%reg(i))
  !---------------------------------------------------------------------------
  ! Quadrature Loop
  !---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    CALL eval_basis(oft_blagrange,i,quad%pts(:,m),jac_mat,basis_vals,basis_grads)
    coords = mesh%log2phys(i,quad%pts(:,m))
    int_factor = jac_det*quad%wts(m)
    CALL eval_fields(self%parent_sim,pt,coords,basis_vals,basis_grads,n_loc,vel_loc,psi_loc,by_loc, &
      basis_vals,basis_grads,T_loc)
    rmul = pt%rmul; rden = pt%rden

    diag_vals(1) = diag_vals(1) + DOT_PRODUCT(pt%btmp,pt%btmp)*int_factor*rmul
    diag_vals(2) = diag_vals(2) + 0.5d0*pt%m_i*pt%n*DOT_PRODUCT(pt%vel,pt%vel)*int_factor*rmul
    diag_vals(3) = diag_vals(3) + pt%n*int_factor*rmul
    diag_vals(4) = diag_vals(4) + pt%T*int_factor*rmul
    diag_vals(5) = diag_vals(5) + int_factor*rmul !total volume
    !---Compute local function contributions
    DO jr=1,oft_blagrange%nce
      !---Diffusion
      res_loc(jr,1) = res_loc(jr, 1) &
        + basis_vals(jr)*pt%n*int_factor*rmul
      !---Momentum
      res_loc(jr, 2:4) = res_loc(jr, 2:4) &
      + basis_vals(jr)*pt%vel*int_factor*rmul
      !---Temperature (compressible only)
      IF(.NOT.incomp) res_loc(jr,5) = res_loc(jr, 5) &
        + basis_vals(jr)*pt%T*int_factor*rmul/(pt%gamma-1)
      !---Induction
      res_loc(jr, 6) = res_loc(jr, 6) &
        + basis_vals(jr)*pt%psi*int_factor/(pt%eta(2)*rden)
      res_loc(jr, 7) = res_loc(jr, 7) &
      + basis_vals(jr)*pt%by*int_factor/rden
    END DO
  END DO
  !---Add local values to full vector
  !$omp ordered
  DO jr=1,oft_blagrange%nce
    !$omp atomic
    n_res(cell_dofs(jr)) = n_res(cell_dofs(jr)) + res_loc(jr,1)/self%parent_sim%den_scale
    !$omp atomic
    velx_res(cell_dofs(jr)) = velx_res(cell_dofs(jr)) + res_loc(jr,2)
    !$omp atomic
    vely_res(cell_dofs(jr)) = vely_res(cell_dofs(jr)) + res_loc(jr,3)
    !$omp atomic
    velz_res(cell_dofs(jr)) = velz_res(cell_dofs(jr)) + res_loc(jr,4)
    IF(.NOT.incomp)THEN !---Skip field 5 (pressure)for incompressible
      !$omp atomic
      T_res(cell_dofs(jr)) = T_res(cell_dofs(jr)) + res_loc(jr,5)
    END IF
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr,6)
    !$omp atomic
    by_res(cell_dofs(jr)) = by_res(cell_dofs(jr)) + res_loc(jr,7)
  END DO
  !$omp end ordered
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads,n_loc,T_loc,vel_loc,psi_loc,by_loc,cell_dofs,res_loc)
!$omp end parallel
END BLOCK
IF(oft_debug_print(2))write(*,'(4X,A)')'Applying BCs'
CALL fem_dirichlet_vec(oft_blagrange,n_weights,n_res,self%parent_sim%n_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(1, :),velx_res,self%parent_sim%velx_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(2, :),vely_res,self%parent_sim%vely_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(3, :),velz_res,self%parent_sim%velz_bc)
IF(incomp)THEN
  !---Field 5 is the pressure, which uses the lower-order pressure space
  CALL fem_dirichlet_vec(oft_blagrange_p,T_weights,T_res,self%parent_sim%T_bc)
ELSE
  CALL fem_dirichlet_vec(oft_blagrange,T_weights,T_res,self%parent_sim%T_bc)
END IF
CALL fem_dirichlet_vec(oft_blagrange,psi_weights,psi_res,self%parent_sim%psi_bc)
CALL fem_dirichlet_vec(oft_blagrange,by_weights,by_res,self%parent_sim%by_bc)

!---Put results into full vector
CALL b%release_local(n_res,OFT_ACCESS_ACCUMULATE,.FALSE.,1)
CALL b%release_local(velx_res,OFT_ACCESS_ACCUMULATE,.FALSE.,2)
CALL b%release_local(vely_res,OFT_ACCESS_ACCUMULATE,.FALSE.,3)
CALL b%release_local(velz_res,OFT_ACCESS_ACCUMULATE,.FALSE.,4)
CALL b%release_local(T_res,OFT_ACCESS_ACCUMULATE,.FALSE.,5)
CALL b%release_local(psi_res,OFT_ACCESS_ACCUMULATE,.FALSE.,6)
CALL b%release_local(by_res,OFT_ACCESS_ACCUMULATE,.TRUE.,7) ! Completes accumulation (sum across seams)
self%diag_vals=oft_mpi_sum(diag_vals,5)
!---Cleanup remaining storage
DEALLOCATE(n_weights,vel_weights, T_weights, psi_weights, by_weights)
DEALLOCATE(psi_vac_weights)
end subroutine mfun_apply
!---------------------------------------------------------------------------
!> Compute the NL error function, where we are solving F(x) = 0
!!
!! b = F(a)
!---------------------------------------------------------------------------
subroutine nlfun_apply(self,a,b)
class(xmhd_2d_nlfun), intent(inout) :: self !< NL function object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
type(oft_quad_type), pointer :: quad
CLASS(oft_scalar_bfem), POINTER :: T_fe
LOGICAL :: cyl_flag, incomp
INTEGER(i4) :: i
REAL(r8) :: m_i=proton_mass
REAL(r8) :: chi, eta(2), nu, D_diff, gamma, diag_vals(5), B_0(3)
REAL(r8), POINTER, DIMENSION(:) :: n_weights,T_weights,psi_weights,by_weights,psi_vac_weights
REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:) :: T_res,n_res,psi_res,by_res,velx_res,vely_res,velz_res
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights
quad=>oft_blagrange%quad
NULLIFY(n_res, velx_res, vely_res, velz_res, T_res, psi_res, by_res)

!---Get weights from solution vector
CALL get_field_weights(self%parent_sim,a,n_weights,vel_weights,T_weights,psi_weights,by_weights,psi_vac_weights)

B_0 = self%parent_sim%B_0
cyl_flag = self%parent_sim%cyl_flag
!---Check that B_0 is only vertical in cylindrical coordinates
IF(cyl_flag.AND.ANY(B_0(1:2)/=0.d0))CALL oft_abort( &
  "Only vertical B_0 (B_0(3)) is supported in cylindrical coordinates","nlfun_apply",__FILE__)
incomp = self%parent_sim%incomp
!---Field 5 is pressure on its own space when incompressible, otherwise temperature
IF(incomp)THEN
  T_fe=>oft_blagrange_p
ELSE
  T_fe=>oft_blagrange
END IF

!---Zero result and get storage array
CALL b%set(0.d0)
CALL b%borrow_local(n_res,OFT_ACCESS_ACCUMULATE,1)
CALL b%borrow_local(velx_res,OFT_ACCESS_ACCUMULATE,2)
CALL b%borrow_local(vely_res,OFT_ACCESS_ACCUMULATE,3)
CALL b%borrow_local(velz_res,OFT_ACCESS_ACCUMULATE,4)
CALL b%borrow_local(T_res,OFT_ACCESS_ACCUMULATE,5)
CALL b%borrow_local(psi_res,OFT_ACCESS_ACCUMULATE,6)
CALL b%borrow_local(by_res,OFT_ACCESS_ACCUMULATE,7)
diag_vals=0.d0

BLOCK
LOGICAL :: curved
INTEGER(i4) :: k,m,jr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs, cell_dofs_T
REAL(r8) :: n,vel(3),T,psi,by,dT(3),dn(3),dpsi(3),dpsi_0(3),dby(3),rmul,rden
REAL(r8) :: dvel(3,3),div_vel,jac_mat(3,4),jac_det,int_factor,btmp(3),tmp1(3),coords(3)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals,basis_vals_T,T_loc,n_loc,psi_loc,by_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: vel_loc,basis_grads,basis_grads_T,res_loc
TYPE(xmhd_2d_qpt) :: pt
!$omp parallel private(k,m,jr,curved,coords,cell_dofs,cell_dofs_T,basis_vals,basis_grads,basis_vals_T, &
!$omp basis_grads_T,T_loc,n_loc,psi_loc,by_loc,vel_loc,res_loc,jac_mat,jac_det,int_factor,pt, &
!$omp T,n,psi,by,vel,dT,dn,dpsi,dpsi_0,dby,dvel,div_vel,btmp,tmp1,rmul,rden, &
!$omp chi, m_i, eta, nu, gamma, D_diff) reduction(+:diag_vals)
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(basis_vals_T(T_fe%nce),basis_grads_T(3,T_fe%nce))
ALLOCATE(n_loc(oft_blagrange%nce),psi_loc(oft_blagrange%nce),by_loc(oft_blagrange%nce), &
  vel_loc(3,oft_blagrange%nce),T_loc(T_fe%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce),cell_dofs_T(T_fe%nce),res_loc(oft_blagrange%nce,7))
!$omp do ordered
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  IF(self%parent_sim%ignore_rmask(mesh%reg(i)))CYCLE ! Skip cells in ignored regions
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  IF(incomp)THEN
    call T_fe%ncdofs(i,cell_dofs_T) ! Get global index of local field-5 DOFs
  ELSE
    cell_dofs_T=cell_dofs
  END IF
  res_loc = 0.d0 ! Zero local (cell) contribution to function
  T_loc = T_weights(cell_dofs_T)
  vel_loc = vel_weights(:, cell_dofs)
  n_loc = n_weights(cell_dofs)
  psi_loc = psi_weights(cell_dofs) + psi_vac_weights(cell_dofs)
  by_loc = by_weights(cell_dofs)

  !Set material properties
  chi = self%parent_sim%chi(mesh%reg(i))
  m_i = self%parent_sim%m_i(mesh%reg(i))
  eta = self%parent_sim%eta(mesh%reg(i), :)
  nu = self%parent_sim%nu(mesh%reg(i))
  gamma = self%parent_sim%gamma(mesh%reg(i))
  D_diff = self%parent_sim%D_diff(mesh%reg(i))
  !---------------------------------------------------------------------------
  ! Quadrature Loop
  !---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    CALL eval_basis(oft_blagrange,i,quad%pts(:,m),jac_mat,basis_vals,basis_grads)
    IF(incomp)THEN
      CALL eval_basis(T_fe,i,quad%pts(:,m),jac_mat,basis_vals_T,basis_grads_T)
    ELSE
      basis_vals_T=basis_vals
      basis_grads_T=basis_grads
    END IF
    !--Extract spatial coordinates at current point
    coords = mesh%log2phys(i,quad%pts(:,m))
    int_factor = jac_det*quad%wts(m)
    !---Reconstruct values of solution fields and geometric factors
    !Note: in cylindrical coordinates, dvel is not the actual gradient of the velocity vector,
    ! but the gradient of each component, in a matrix
    CALL eval_fields(self%parent_sim,pt,coords,basis_vals,basis_grads,n_loc,vel_loc,psi_loc,by_loc, &
      basis_vals_T,basis_grads_T,T_loc)
    n = pt%n; dn = pt%dn; vel = pt%vel; dvel = pt%dvel; div_vel = pt%div_vel
    T = pt%T; dT = pt%dT; psi = pt%psi; dpsi = pt%dpsi; by = pt%by; dby = pt%dby
    btmp = pt%btmp; rmul = pt%rmul; rden = pt%rden
    ! dpsi_0 excludes the vacuum flux psi_vac (gradient of the solved psi only)
    dpsi_0 = 0.d0
    DO jr=1,oft_blagrange%nce
      dpsi_0 = dpsi_0 + psi_weights(cell_dofs(jr))*basis_grads(:,jr)
    END DO
    diag_vals(1) = diag_vals(1) + DOT_PRODUCT(btmp,btmp)*int_factor*rmul
    diag_vals(2) = diag_vals(2) + 0.5d0*m_i*n*DOT_PRODUCT(vel,vel)*int_factor*rmul
    diag_vals(3) = diag_vals(3) + n*int_factor*rmul
    diag_vals(4) = diag_vals(4) + T*int_factor*rmul
    diag_vals(5) = diag_vals(5) + int_factor*rmul !total volume
    !---Compute local function contributions
    IF (incomp) THEN
      ! The divergence constraint is algebraic, so it is only applied at the new time level
      ! (dt>0). Omitting it from explicit (RHS) evaluations (dt<=0) enforces div(V^{n+1})=0
      ! instead of div(V^{n+1})=div(V^n).
      IF(self%dt>0.d0)THEN
        DO jr=1,T_fe%nce
          !---Continuity
          res_loc(jr,5) = res_loc(jr, 5) &
            + basis_vals_T(jr)*div_vel*int_factor*rmul
        END DO
      END IF
    ELSE
      DO jr=1,oft_blagrange%nce
        !---Temperature
        IF (cyl_flag) THEN
          res_loc(jr,5) = res_loc(jr, 5) &
            + basis_vals(jr)*T*int_factor*rmul/(gamma-1) &
            + self%dt*basis_vals(jr)*DOT_PRODUCT(vel, dT)*int_factor*rmul/(gamma-1) &
            + self%dt*basis_vals(jr)*T*div_vel*int_factor*rmul &
            + self%dt*chi*DOT_PRODUCT(dT, basis_grads(:,jr))*int_factor*rmul &
            - self%dt*chi*basis_vals(jr)*DOT_PRODUCT(dn, dT)*int_factor*rmul/n
        ELSE
          res_loc(jr,5) = res_loc(jr, 5) &
            + basis_vals(jr)*T*int_factor/(gamma-1) &
            + self%dt*basis_vals(jr)*DOT_PRODUCT(vel, dT)*int_factor/(gamma-1) &
            + self%dt*basis_vals(jr)*T*div_vel*int_factor & 
            + self%dt*chi*DOT_PRODUCT(dT, basis_grads(:,jr))*int_factor &
            - self%dt*chi*basis_vals(jr)*DOT_PRODUCT(dn, dT)*int_factor/n 
        END IF
      END DO
    END IF

    DO jr=1,oft_blagrange%nce
      !---Diffusion
      IF(cyl_flag)THEN
        res_loc(jr,1) = res_loc(jr, 1) &
          + basis_vals(jr)*n*int_factor*rmul &
          +self%dt*basis_vals(jr)*DOT_PRODUCT(dn, vel)*int_factor*rmul &
          +self%dt*basis_vals(jr)*n*div_vel*int_factor*rmul &
          +self%dt*D_diff*DOT_PRODUCT(dn,basis_grads(:,jr))*int_factor*rmul
      ELSE
        res_loc(jr,1) = res_loc(jr, 1) &
          + basis_vals(jr)*n*int_factor &
          +self%dt*basis_vals(jr)*DOT_PRODUCT(dn, vel)*int_factor &
          +self%dt*basis_vals(jr)*n*div_vel*int_factor &
          +self%dt*D_diff*DOT_PRODUCT(dn,basis_grads(:,jr))*int_factor
      END IF
      !---Momentum
      IF(cyl_flag)THEN
        res_loc(jr, 2:4) = res_loc(jr, 2:4) &
          + basis_vals(jr)*vel*int_factor*rmul &
          + self%dt*DOT_PRODUCT(btmp,basis_grads(:,jr))*btmp*int_factor*rmul/(mu0*m_i*n) &
          - self%dt*DOT_PRODUCT(btmp,btmp)*basis_grads(:,jr)*int_factor*rmul/(2*mu0*m_i*n)
        IF (incomp) THEN
          res_loc(jr, 2:4) = res_loc(jr, 2:4) &
          -self%dt*T*basis_grads(:,jr)*int_factor*rmul/(m_i*n)
          res_loc(jr,2) = res_loc(jr,2) &
          -self%dt*basis_vals(jr)*T*int_factor/(m_i*n)
        ELSE
          res_loc(jr, 2:4) = res_loc(jr, 2:4) &
          - self%dt*basis_vals(jr)*DOT_PRODUCT(dn,btmp)*btmp*int_factor*rmul/(mu0*m_i*n**2) &
          + self%dt*basis_vals(jr)*DOT_PRODUCT(btmp,btmp)*dn*int_factor*rmul/(2*mu0*m_i*n**2) &
          + self%dt*basis_vals(jr)*2.d0*k_boltz*dT*int_factor*rmul/m_i &
          + self%dt*basis_vals(jr)*2.d0*k_boltz*T*dn*int_factor*rmul/(m_i*n)
        END IF
        DO k=1,3
          res_loc(jr,k+1) = res_loc(jr, k+1) &
            + basis_vals(jr)*self%dt*DOT_PRODUCT(vel,dvel(k,:))*int_factor*rmul &
            + nu*self%dt*DOT_PRODUCT(basis_grads(:,jr),dvel(k,:))*int_factor*rmul &
            - nu*basis_vals(jr)*self%dt*DOT_PRODUCT(dn,dvel(k,:))*int_factor*rmul/n
        END DO
        res_loc(jr,2) = res_loc(jr,2) &
        + self%dt*basis_vals(jr)*(btmp(2)**2-btmp(1)**2-btmp(3)**2)*int_factor/(2.d0*mu0*m_i*n) &
        + self%dt*basis_vals(jr)*nu*vel(1)*int_factor/rden  &
        - self%dt*basis_vals(jr)*vel(2)**2*int_factor
        res_loc(jr,3) = res_loc(jr,3) &
        - self%dt*basis_vals(jr)*btmp(1)*btmp(2)*int_factor/(mu0*m_i*n) & !!nonzero
        + self%dt*basis_vals(jr)*nu*vel(2)*int_factor/rden &
        + self%dt*basis_vals(jr)*vel(1)*vel(2)*int_factor
      ELSE
        res_loc(jr, 2:4) = res_loc(jr, 2:4) &
          + basis_vals(jr)*vel*int_factor&
          + self%dt*DOT_PRODUCT(btmp,basis_grads(:,jr))*btmp*int_factor/(mu0*m_i*n) & 
          - self%dt*DOT_PRODUCT(btmp,btmp)*basis_grads(:,jr)*int_factor/(2*mu0*m_i*n)
        IF(incomp) THEN
          res_loc(jr, 2:4) = res_loc(jr, 2:4) &
          -self%dt*T*basis_grads(:,jr)*int_factor/(m_i*n)
        ELSE
          res_loc(jr,2:4) = res_loc(jr,2:4) &
          - self%dt*basis_vals(jr)*DOT_PRODUCT(dn,btmp)*btmp*int_factor/(mu0*m_i*n**2) &
          + self%dt*basis_vals(jr)*DOT_PRODUCT(btmp,btmp)*dn*int_factor/(2*mu0*m_i*n**2) &
          + self%dt*basis_vals(jr)*2.d0*k_boltz*dT*int_factor/m_i &
          + self%dt*basis_vals(jr)*2.d0*k_boltz*T*dn*int_factor/(m_i*n)
        END IF
        DO k=1,3
          res_loc(jr,k+1) = res_loc(jr, k+1) &
            + basis_vals(jr)*self%dt*DOT_PRODUCT(vel,dvel(k,:))*int_factor &
            + nu*self%dt*DOT_PRODUCT(basis_grads(:,jr),dvel(k,:))*int_factor &
            - nu*basis_vals(jr)*self%dt*DOT_PRODUCT(dn,dvel(k,:))*int_factor/n
        END DO
      END IF
      !---Psi
      tmp1 = cross_product(B_0,vel)
      IF(cyl_flag)THEN
        res_loc(jr, 6) = res_loc(jr, 6) &
        + basis_vals(jr)*psi*int_factor/(eta(2)*rden) &
        + basis_vals(jr)*self%dt*DOT_PRODUCT(vel, dpsi)*int_factor/(eta(2)*rden) &
        + basis_vals(jr)*self%dt*tmp1(2)*int_factor/eta(2) &
        + self%dt*DOT_PRODUCT(basis_grads(:,jr), dpsi_0)*int_factor/rden
      ELSE
        res_loc(jr, 6) = res_loc(jr, 6) &
        + basis_vals(jr)*psi*int_factor/eta(2) &
        + basis_vals(jr)*self%dt*DOT_PRODUCT(vel, dpsi)*int_factor/eta(2) &
        + basis_vals(jr)*self%dt*tmp1(2)*int_factor/eta(2) &
        + self%dt*DOT_PRODUCT(basis_grads(:,jr), dpsi_0)*int_factor
      END IF
      ! --By
      tmp1 = cross_product(dpsi,dvel(2, :))
      IF(cyl_flag)THEN
        res_loc(jr,7) = res_loc(jr,7) &
        + basis_vals(jr)*by*int_factor/rden &
        + self%dt*basis_vals(jr)*tmp1(2)*int_factor/rden &
        - self%dt*basis_vals(jr)*vel(2)*dpsi(3)*int_factor/rden**2 &
        + self%dt*basis_vals(jr)*dvel(1,1)*by*int_factor/rden &
        + self%dt*basis_vals(jr)*dvel(3,3)*by*int_factor/rden &
        + self%dt*basis_vals(jr)*DOT_PRODUCT(vel, dby)*int_factor/rden &
        - self%dt*basis_vals(jr)*vel(1)*by*int_factor/rden**2 &
        + self%dt*eta(1)*DOT_PRODUCT(basis_grads(:,jr), dby)*int_factor/rden &
        ! Background field: -[curl(v x B_0)]_y, only B_0z allowed)
        - self%dt*basis_vals(jr)*B_0(3)*dvel(2,3)*int_factor
      ELSE
        res_loc(jr, 7) = res_loc(jr, 7) &
        + basis_vals(jr)*by*int_factor &
        + basis_vals(jr)*self%dt*tmp1(2)*int_factor &
        + basis_vals(jr)*self%dt*DOT_PRODUCT(vel, dby)*int_factor &
        + basis_vals(jr)*self%dt*by*div_vel*int_factor &
        + self%dt*eta(1)*DOT_PRODUCT(basis_grads(:,jr), dby)*int_factor &
        ! Background field: -[curl(v x B_0)]_y 
        - self%dt*basis_vals(jr)*DOT_PRODUCT(B_0, dvel(2,:))*int_factor &
        + self%dt*basis_vals(jr)*B_0(2)*div_vel*int_factor
      END IF
    END DO
  END DO
  !---Add local values to full vector
  !$omp ordered
  DO jr=1,T_fe%nce
    !$omp atomic
    T_res(cell_dofs_T(jr)) = T_res(cell_dofs_T(jr)) + res_loc(jr,5)
  END DO
  !$omp end ordered
  !$omp ordered
  DO jr=1,oft_blagrange%nce
    !$omp atomic
    n_res(cell_dofs(jr)) = n_res(cell_dofs(jr)) + res_loc(jr,1)/self%parent_sim%den_scale
    !$omp atomic
    velx_res(cell_dofs(jr)) = velx_res(cell_dofs(jr)) + res_loc(jr,2)
    !$omp atomic
    vely_res(cell_dofs(jr)) = vely_res(cell_dofs(jr)) + res_loc(jr,3)
    !$omp atomic
    velz_res(cell_dofs(jr)) = velz_res(cell_dofs(jr)) + res_loc(jr,4)
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr,6)
    !$omp atomic
    by_res(cell_dofs(jr)) = by_res(cell_dofs(jr)) + res_loc(jr,7)
  END DO
  !$omp end ordered
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads,basis_vals_T,basis_grads_T,n_loc,T_loc,vel_loc,psi_loc,by_loc, &
  cell_dofs,cell_dofs_T,res_loc)
!$omp end parallel
END BLOCK
IF(oft_debug_print(2))write(*,'(4X,A)')'Applying BCs'
CALL fem_dirichlet_vec(T_fe,T_weights,T_res,self%parent_sim%T_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(1, :),velx_res,self%parent_sim%velx_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(2, :),vely_res,self%parent_sim%vely_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(3, :),velz_res,self%parent_sim%velz_bc)
CALL fem_dirichlet_vec(oft_blagrange,n_weights,n_res,self%parent_sim%n_bc)
CALL fem_dirichlet_vec(oft_blagrange,psi_weights,psi_res,self%parent_sim%psi_bc)
CALL fem_dirichlet_vec(oft_blagrange,by_weights,by_res,self%parent_sim%by_bc)
!---Put results into full vector
CALL b%release_local(n_res,OFT_ACCESS_ACCUMULATE,.FALSE.,1)
CALL b%release_local(velx_res,OFT_ACCESS_ACCUMULATE,.FALSE.,2)
CALL b%release_local(vely_res,OFT_ACCESS_ACCUMULATE,.FALSE.,3)
CALL b%release_local(velz_res,OFT_ACCESS_ACCUMULATE,.FALSE.,4)
CALL b%release_local(T_res,OFT_ACCESS_ACCUMULATE,.FALSE.,5)
CALL b%release_local(psi_res,OFT_ACCESS_ACCUMULATE,.FALSE.,6)
CALL b%release_local(by_res,OFT_ACCESS_ACCUMULATE,.TRUE.,7) ! Completes accumulation (sum across seams)
self%diag_vals=oft_mpi_sum(diag_vals,5)
!---Cleanup remaining storage
DEALLOCATE(n_weights,vel_weights, T_weights, psi_weights, by_weights)
DEALLOCATE(psi_vac_weights)
END SUBROUTINE nlfun_apply
!---------------------------------------------------------------------------
!> Compute the approximate Jacobian matrix for the compressible system
!!
!! Field 5 is temperature, represented on the same Lagrange space as the
!! density, velocity, and magnetic fields.
!---------------------------------------------------------------------------
SUBROUTINE build_approx_jacobian_comp(self,a)
CLASS(oft_xmhd_2d_sim), INTENT(INOUT) :: self !< Simulation object
CLASS(oft_vector), INTENT(INOUT) :: a !< Solution for computing Jacobian
CALL build_approx_jacobian(self,a,oft_blagrange,.FALSE.)
END SUBROUTINE build_approx_jacobian_comp
!---------------------------------------------------------------------------
!> Compute the approximate Jacobian matrix for the incompressible system
!!
!! Field 5 is pressure, represented on a separate Lagrange space
!! (@ref oft_blagrange_p) and acting as a Lagrange multiplier for the
!! divergence-free constraint on velocity.
!---------------------------------------------------------------------------
SUBROUTINE build_approx_jacobian_incomp(self,a)
CLASS(oft_xmhd_2d_sim), INTENT(INOUT) :: self !< Simulation object
CLASS(oft_vector), INTENT(INOUT) :: a !< Solution for computing Jacobian
CALL build_approx_jacobian(self,a,oft_blagrange_p,.TRUE.)
END SUBROUTINE build_approx_jacobian_incomp
!---------------------------------------------------------------------------
!> Assemble the approximate Jacobian matrix for either system
!!
!! Loops over cells and quadrature points, adding the terms shared by both
!! systems (@ref jac_add_shared) and then the field-5 couplings of the chosen
!! system (@ref jac_add_comp or @ref jac_add_incomp). Field 5 (temperature or
!! pressure) is represented on `T_fe`; in the compressible system that is
!! @ref oft_blagrange itself, so its basis and DOFs are copied rather than
!! evaluated a second time.
!---------------------------------------------------------------------------
SUBROUTINE build_approx_jacobian(self,a,T_fe,incomp)
CLASS(oft_xmhd_2d_sim), INTENT(INOUT) :: self !< Simulation object
CLASS(oft_vector), INTENT(INOUT) :: a !< Solution for computing Jacobian
CLASS(oft_scalar_bfem), INTENT(INOUT) :: T_fe !< Lagrange representation for field 5
LOGICAL, INTENT(IN) :: incomp !< Field 5 is pressure (incompressible) rather than temperature
INTEGER(i4) :: i
REAL(r8) :: dt_fac
REAL(r8), POINTER, DIMENSION(:) :: n_weights,T_weights,psi_weights,by_weights,psi_vac_weights
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights
INTEGER(KIND=omp_lock_kind), ALLOCATABLE, DIMENSION(:) :: tlocks
TYPE(oft_quad_type), POINTER :: quad
quad=>oft_blagrange%quad
CALL self%jacobian%zero
CALL get_field_weights(self,a,n_weights,vel_weights,T_weights,psi_weights,by_weights,psi_vac_weights)
dt_fac=self%jac_dt
!---Setup thread locks
ALLOCATE(tlocks(self%fe_rep%nfields))
DO i=1,self%fe_rep%nfields
  CALL omp_init_lock(tlocks(i))
END DO
BLOCK
LOGICAL :: curved
INTEGER(i4) :: jr,m
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs,cell_dofs_T
REAL(r8) :: jac_mat(3,4),jac_det,coords(3)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals,basis_vals_T,n_loc,T_loc,psi_loc,by_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads,basis_grads_T,vel_loc
TYPE(xmhd_2d_qpt) :: pt
TYPE(oft_1d_int), ALLOCATABLE, DIMENSION(:) :: iloc
TYPE(oft_local_mat), ALLOCATABLE, DIMENSION(:,:) :: jac_loc
!$omp parallel private(jr,m,curved,cell_dofs,cell_dofs_T,jac_mat,jac_det,coords,basis_vals,basis_grads, &
!$omp basis_vals_T,basis_grads_T,n_loc,vel_loc,T_loc,psi_loc,by_loc,pt,iloc,jac_loc)
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(basis_vals_T(T_fe%nce),basis_grads_T(3,T_fe%nce))
ALLOCATE(n_loc(oft_blagrange%nce),vel_loc(3,oft_blagrange%nce),T_loc(T_fe%nce), &
  psi_loc(oft_blagrange%nce),by_loc(oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce),cell_dofs_T(T_fe%nce))
ALLOCATE(jac_loc(self%fe_rep%nfields,self%fe_rep%nfields),iloc(self%fe_rep%nfields))
DO jr=1,self%fe_rep%nfields
  iloc(jr)%v=>cell_dofs
END DO
iloc(5)%v=>cell_dofs_T
CALL self%fe_rep%mat_setup_local(jac_loc)
!$omp do ordered
DO i=1,mesh%nc
  IF(self%ignore_rmask(mesh%reg(i)))CYCLE ! Skip cells in ignored regions
  curved=cell_is_curved(mesh,i) ! Straight cell test
  CALL oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  IF(incomp)THEN
    CALL T_fe%ncdofs(i,cell_dofs_T) ! Get global index of local field-5 DOFs
  ELSE
    cell_dofs_T=cell_dofs
  END IF
  CALL self%fe_rep%mat_zero_local(jac_loc) ! Zero local (cell) contribution to matrix
  n_loc=n_weights(cell_dofs)
  vel_loc=vel_weights(:,cell_dofs)
  T_loc=T_weights(cell_dofs_T)
  psi_loc=psi_weights(cell_dofs)+psi_vac_weights(cell_dofs)
  by_loc=by_weights(cell_dofs)
  !---Set material properties
  pt%chi=self%chi(mesh%reg(i))
  pt%m_i=self%m_i(mesh%reg(i))
  pt%eta=self%eta(mesh%reg(i),:)
  pt%nu=self%nu(mesh%reg(i))
  pt%gamma=self%gamma(mesh%reg(i))
  pt%D_diff=self%D_diff(mesh%reg(i))
  !---Quadrature loop
  DO m=1,quad%np
    IF(curved.OR.(m==1))CALL mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    CALL eval_basis(oft_blagrange,i,quad%pts(:,m),jac_mat,basis_vals,basis_grads)
    IF(incomp)THEN
      CALL eval_basis(T_fe,i,quad%pts(:,m),jac_mat,basis_vals_T,basis_grads_T)
    ELSE
      basis_vals_T=basis_vals
      basis_grads_T=basis_grads
    END IF
    coords=mesh%log2phys(i,quad%pts(:,m))
    pt%int_factor=jac_det*quad%wts(m)
    CALL eval_fields(self,pt,coords,basis_vals,basis_grads,n_loc,vel_loc,psi_loc,by_loc, &
      basis_vals_T,basis_grads_T,T_loc)
    CALL jac_add_shared(self,jac_loc,pt,basis_vals,basis_grads,dt_fac)
    IF(incomp)THEN
      CALL jac_add_incomp(self,jac_loc,pt,basis_vals,basis_grads,basis_vals_T,dt_fac)
    ELSE
      CALL jac_add_comp(self,jac_loc,pt,basis_vals,basis_grads,dt_fac)
    END IF
  END DO
  CALL jac_finish_local(self,jac_loc,cell_dofs,cell_dofs_T)
  !---Add local contributions to matrix
  !$omp ordered
  CALL self%fe_rep%mat_add_local(self%jacobian,jac_loc,iloc,tlocks)
  !$omp end ordered
END DO
!---Cleanup thread-local storage
CALL self%fe_rep%mat_destroy_local(jac_loc)
DEALLOCATE(basis_vals,basis_grads,basis_vals_T,basis_grads_T,n_loc,vel_loc,T_loc,psi_loc,by_loc, &
  cell_dofs,cell_dofs_T,jac_loc,iloc)
!$omp end parallel
END BLOCK
!---Destroy thread locks
DO i=1,self%fe_rep%nfields
  CALL omp_destroy_lock(tlocks(i))
END DO
DEALLOCATE(tlocks)
CALL jac_finalize(self,T_fe)
DEALLOCATE(n_weights,vel_weights,T_weights,psi_weights,by_weights,psi_vac_weights)
END SUBROUTINE build_approx_jacobian
!---------------------------------------------------------------------------
!> Get local weights of all fields from a solution vector
!!
!! The poloidal flux weights are those of the solved field alone; the vacuum
!! offset (@ref oft_xmhd_2d_sim::psi_vac) is returned separately, so callers
!! that need the total flux add the two.
!---------------------------------------------------------------------------
SUBROUTINE get_field_weights(self,a,n_weights,vel_weights,T_weights,psi_weights,by_weights,psi_vac_weights)
CLASS(oft_xmhd_2d_sim), INTENT(INOUT) :: self !< Simulation object
CLASS(oft_vector), INTENT(INOUT) :: a !< Solution vector
REAL(r8), POINTER, DIMENSION(:), INTENT(OUT) :: n_weights !< Density weights
REAL(r8), POINTER, DIMENSION(:,:), INTENT(OUT) :: vel_weights !< Velocity weights [3,ne]
REAL(r8), POINTER, DIMENSION(:), INTENT(OUT) :: T_weights !< Temperature (or pressure) weights
REAL(r8), POINTER, DIMENSION(:), INTENT(OUT) :: psi_weights !< Poloidal flux weights (solved field only)
REAL(r8), POINTER, DIMENSION(:), INTENT(OUT) :: by_weights !< Out-of-plane field weights
REAL(r8), POINTER, DIMENSION(:), INTENT(OUT) :: psi_vac_weights !< Vacuum poloidal flux offset weights
REAL(r8), POINTER, DIMENSION(:) :: vtmp
NULLIFY(n_weights,vel_weights,T_weights,psi_weights,by_weights,psi_vac_weights)
CALL a%get_local(n_weights,1)
ALLOCATE(vel_weights(3,oft_blagrange%ne))
vtmp=>vel_weights(1,:)
CALL a%get_local(vtmp,2)
vtmp=>vel_weights(2,:)
CALL a%get_local(vtmp,3)
vtmp=>vel_weights(3,:)
CALL a%get_local(vtmp,4)
CALL a%get_local(T_weights,5)
CALL a%get_local(psi_weights,6)
CALL a%get_local(by_weights,7)
CALL self%psi_vac%get_local(psi_vac_weights)
END SUBROUTINE get_field_weights
!---------------------------------------------------------------------------
!> Evaluate a 2D Lagrange basis and its gradients at a point in a cell
!!
!! Gradients are returned in (x, y, z) ordering with the in-plane components in
!! x and z and the out-of-plane (y) component set to zero.
!---------------------------------------------------------------------------
SUBROUTINE eval_basis(fe,cell,f,jac_mat,basis_vals,basis_grads)
CLASS(oft_scalar_bfem), INTENT(IN) :: fe !< Lagrange representation
INTEGER(i4), INTENT(IN) :: cell !< Cell for evaluation
REAL(r8), INTENT(IN) :: f(:) !< Position in cell in logical space
REAL(r8), INTENT(IN) :: jac_mat(3,4) !< Cell Jacobian matrix at point (f)
REAL(r8), INTENT(OUT) :: basis_vals(:) !< Basis values [nce]
REAL(r8), INTENT(OUT) :: basis_grads(:,:) !< Basis gradients [3,nce]
INTEGER(i4) :: jr
DO jr=1,fe%nce
  CALL oft_blag_eval(fe,cell,jr,f,basis_vals(jr))
  CALL oft_blag_geval(fe,cell,jr,f,basis_grads(:,jr),jac_mat)
END DO
basis_grads(3,:)=basis_grads(2,:)
basis_grads(2,:)=0.d0
END SUBROUTINE eval_basis
!---------------------------------------------------------------------------
!> Reconstruct all fields and the geometric factors at a quadrature point
!!
!! Field 5 (temperature, or pressure if incompressible) may live on its own
!! Lagrange space, so it takes its own basis.
!---------------------------------------------------------------------------
SUBROUTINE eval_fields(self,pt,coords,basis_vals,basis_grads,n_loc,vel_loc,psi_loc,by_loc, &
  basis_vals_T,basis_grads_T,T_loc)
CLASS(oft_xmhd_2d_sim), INTENT(IN) :: self !< Simulation object
TYPE(xmhd_2d_qpt), INTENT(INOUT) :: pt !< Quadrature point state
REAL(r8), INTENT(IN) :: coords(3) !< Physical coordinates of point
REAL(r8), INTENT(IN) :: basis_vals(:) !< Basis values [nce]
REAL(r8), INTENT(IN) :: basis_grads(:,:) !< Basis gradients [3,nce]
REAL(r8), INTENT(IN) :: n_loc(:) !< Local density weights [nce]
REAL(r8), INTENT(IN) :: vel_loc(:,:) !< Local velocity weights [3,nce]
REAL(r8), INTENT(IN) :: psi_loc(:) !< Local poloidal flux weights, including vacuum offset [nce]
REAL(r8), INTENT(IN) :: by_loc(:) !< Local out-of-plane field weights [nce]
REAL(r8), INTENT(IN) :: basis_vals_T(:) !< Field-5 basis values [nce_T]
REAL(r8), INTENT(IN) :: basis_grads_T(:,:) !< Field-5 basis gradients [3,nce_T]
REAL(r8), INTENT(IN) :: T_loc(:) !< Local field-5 weights [nce_T]
INTEGER(i4) :: jr
pt%n=0.d0; pt%dn=0.d0; pt%vel=0.d0; pt%dvel=0.d0
pt%psi=0.d0; pt%dpsi=0.d0; pt%by=0.d0; pt%dby=0.d0
DO jr=1,SIZE(basis_vals)
  pt%vel=pt%vel+vel_loc(:,jr)*basis_vals(jr)
  pt%n=pt%n+n_loc(jr)*basis_vals(jr)
  pt%psi=pt%psi+psi_loc(jr)*basis_vals(jr)
  pt%by=pt%by+by_loc(jr)*basis_vals(jr)
  pt%dvel(:,1)=pt%dvel(:,1)+vel_loc(:,jr)*basis_grads(1,jr)
  pt%dvel(:,3)=pt%dvel(:,3)+vel_loc(:,jr)*basis_grads(3,jr)
  pt%dn=pt%dn+n_loc(jr)*basis_grads(:,jr)
  pt%dpsi=pt%dpsi+psi_loc(jr)*basis_grads(:,jr)
  pt%dby=pt%dby+by_loc(jr)*basis_grads(:,jr)
END DO
pt%n=pt%n*self%den_scale
pt%dn=pt%dn*self%den_scale
pt%T=0.d0; pt%dT=0.d0
DO jr=1,SIZE(basis_vals_T)
  pt%T=pt%T+T_loc(jr)*basis_vals_T(jr)
  pt%dT=pt%dT+T_loc(jr)*basis_grads_T(:,jr)
END DO
!---Geometric factors and velocity divergence
IF(self%cyl_flag)THEN
  pt%rmul=coords(1)
  pt%rden=coords(1)+gs_epsilon
  pt%div_vel=pt%dvel(1,1)+pt%vel(1)/pt%rden+pt%dvel(3,3)
ELSE
  pt%rmul=1.d0
  pt%rden=1.d0
  pt%div_vel=pt%dvel(1,1)+pt%dvel(3,3)
END IF
pt%btmp=cross_product(pt%dpsi/pt%rden,[0.d0,1.d0,0.d0])+pt%by*[0.d0,1.d0,0.d0]/pt%rden+self%B_0
END SUBROUTINE eval_fields
!---------------------------------------------------------------------------
!> Add Jacobian contributions common to the compressible and incompressible systems
!!
!! Covers the density, momentum, and induction equations, excluding couplings
!! to field 5 (temperature or pressure). Terms are written once for both
!! geometries using the factors in @ref xmhd_2d_qpt (unity in Cartesian
!! geometry); the extra hoop-stress and curvature terms in cylindrical geometry
!! follow each block.
!---------------------------------------------------------------------------
SUBROUTINE jac_add_shared(self,jac_loc,pt,basis_vals,basis_grads,dt_fac)
CLASS(oft_xmhd_2d_sim), INTENT(IN) :: self !< Simulation object
TYPE(oft_local_mat), INTENT(INOUT) :: jac_loc(:,:) !< Local (cell) Jacobian
TYPE(xmhd_2d_qpt), INTENT(IN) :: pt !< Quadrature point state
REAL(r8), INTENT(IN) :: basis_vals(:) !< Basis values [nce]
REAL(r8), INTENT(IN) :: basis_grads(:,:) !< Basis gradients [3,nce]
REAL(r8), INTENT(IN) :: dt_fac !< Time step factor
LOGICAL :: cyl_flag,linear,incomp
INTEGER(i4) :: jr,jc,k,l,nce
REAL(r8) :: n,dn(3),vel(3),dvel(3,3),div_vel,T,dT(3),dpsi(3),by,dby(3),btmp(3),B_0(3),dn_c(3)
REAL(r8) :: rmul,rden,int_factor,m_i,nu,eta(2),D_diff
REAL(r8) :: db_psi(3),db_by(3),tmp2(3),b0_cross(3,3),psi_cross(3),dvel_cross(3)
!---Local copies (avoids reloads through possible aliasing with jac_loc pointers)
cyl_flag=self%cyl_flag; linear=self%linear; incomp=self%incomp; B_0=self%B_0
n=pt%n; dn=pt%dn; vel=pt%vel; dvel=pt%dvel; div_vel=pt%div_vel; T=pt%T; dT=pt%dT
dpsi=pt%dpsi; by=pt%by; dby=pt%dby; btmp=pt%btmp; rmul=pt%rmul; rden=pt%rden
int_factor=pt%int_factor; m_i=pt%m_i; nu=pt%nu; eta=pt%eta; D_diff=pt%D_diff
nce=SIZE(basis_vals)
!---The grad(n) terms of the Maxwell stress are only included in the compressible model
IF(incomp)THEN
  dn_c=0.d0
ELSE
  dn_c=dn
END IF
!---B_0 x e_l, which does not depend on the basis functions
DO l=1,3
  tmp2 = 0.d0
  tmp2(l) = 1.d0
  b0_cross(:,l) = cross_product(B_0,tmp2)
END DO
!---Terms common to every geometry and mode. jc is the outer loop so that quantities depending
!   only on the column basis function are computed once per column, and jr runs along memory.
DO jc=1,nce
  db_psi = cross_product(basis_grads(:,jc),[0.d0,1.d0/rden,0.d0]) ! Perturbed B from delta_psi
  db_by = [0.d0,basis_vals(jc)/rden,0.d0] ! Perturbed B from delta_by
  psi_cross = cross_product(dpsi,basis_grads(:,jc))
  dvel_cross = cross_product(basis_grads(:,jc),dvel(2,:))
  DO jr=1,nce
    !---Continuity
    ! n, n
    jac_loc(1,1)%m(jr,jc) = jac_loc(1,1)%m(jr,jc) &
      + basis_vals(jr)*basis_vals(jc)*int_factor*rmul &
      + dt_fac*basis_vals(jr)*DOT_PRODUCT(basis_grads(:,jc),vel)*int_factor*rmul & ! u dot grad(delta_n)
      + dt_fac*basis_vals(jr)*basis_vals(jc)*div_vel*int_factor*rmul & ! delta_n*div(u)
      + dt_fac*D_diff*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*int_factor*rmul
    ! n, vel
    DO l=1,3
      jac_loc(1,l+1)%m(jr,jc) = jac_loc(1,l+1)%m(jr,jc) &
        + basis_vals(jr)*dt_fac*n*basis_grads(l,jc)*int_factor*rmul & ! n*div(delta_u)
        + basis_vals(jr)*dt_fac*basis_vals(jc)*dn(l)*int_factor*rmul ! delta_u dot grad(n)
    END DO
    !---Momentum
    ! vel, vel
    DO k=1,3
      jac_loc(k+1,k+1)%m(jr,jc) = jac_loc(k+1,k+1)%m(jr,jc) &
        + basis_vals(jr)*basis_vals(jc)*int_factor*rmul &
        + dt_fac*basis_vals(jr)*DOT_PRODUCT(vel,basis_grads(:,jc))*int_factor*rmul &
        + dt_fac*nu*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*int_factor*rmul &
        - dt_fac*nu*basis_vals(jr)*DOT_PRODUCT(dn,basis_grads(:,jc))*int_factor*rmul/n
      DO l=1,3
        jac_loc(k+1,l+1)%m(jr,jc) = jac_loc(k+1,l+1)%m(jr,jc) &
          + dt_fac*basis_vals(jr)*basis_vals(jc)*dvel(k,l)*int_factor*rmul
      END DO
    END DO
    ! vel, psi
    DO l=1,3
      jac_loc(l+1,6)%m(jr,jc) = jac_loc(l+1,6)%m(jr,jc) &
        + dt_fac*DOT_PRODUCT(basis_grads(:,jr),db_psi)*btmp(l)*int_factor*rmul/(m_i*n*mu0) &
        + dt_fac*DOT_PRODUCT(basis_grads(:,jr),btmp)*db_psi(l)*int_factor*rmul/(m_i*n*mu0) &
        - dt_fac*basis_grads(l,jr)*DOT_PRODUCT(db_psi,btmp)*int_factor*rmul/(m_i*n*mu0) &
        - dt_fac*basis_vals(jr)*DOT_PRODUCT(dn_c,db_psi)*btmp(l)*int_factor*rmul/(m_i*n**2*mu0) &
        - dt_fac*basis_vals(jr)*DOT_PRODUCT(dn_c,btmp)*db_psi(l)*int_factor*rmul/(m_i*n**2*mu0) &
        + dt_fac*basis_vals(jr)*dn_c(l)*DOT_PRODUCT(db_psi,btmp)*int_factor*rmul/(m_i*n**2*mu0)
    END DO
    ! vel, by
    DO l=1,3
      jac_loc(l+1,7)%m(jr,jc) = jac_loc(l+1,7)%m(jr,jc) &
        + dt_fac*DOT_PRODUCT(basis_grads(:,jr),db_by)*btmp(l)*int_factor*rmul/(m_i*n*mu0) &
        + dt_fac*DOT_PRODUCT(basis_grads(:,jr),btmp)*db_by(l)*int_factor*rmul/(m_i*n*mu0) &
        - dt_fac*basis_grads(l,jr)*DOT_PRODUCT(db_by,btmp)*int_factor*rmul/(m_i*n*mu0) &
        - dt_fac*basis_vals(jr)*DOT_PRODUCT(dn_c,db_by)*btmp(l)*int_factor*rmul/(m_i*n**2*mu0) &
        - dt_fac*basis_vals(jr)*DOT_PRODUCT(dn_c,btmp)*db_by(l)*int_factor*rmul/(m_i*n**2*mu0) &
        + dt_fac*basis_vals(jr)*dn_c(l)*DOT_PRODUCT(db_by,btmp)*int_factor*rmul/(m_i*n**2*mu0)
    END DO
    !---Induction
    ! psi, vel
    DO l=1,3
      jac_loc(6,l+1)%m(jr,jc) = jac_loc(6,l+1)%m(jr,jc) &
        + basis_vals(jr)*dt_fac*basis_vals(jc)*dpsi(l)*int_factor/(eta(2)*rden) &
        + basis_vals(jr)*dt_fac*basis_vals(jc)*b0_cross(2,l)*int_factor/eta(2)
    END DO
    ! psi, psi
    jac_loc(6,6)%m(jr,jc) = jac_loc(6,6)%m(jr,jc) &
      + basis_vals(jr)*basis_vals(jc)*int_factor/(eta(2)*rden) &
      + dt_fac*basis_vals(jr)*DOT_PRODUCT(vel,basis_grads(:,jc))*int_factor/(eta(2)*rden) &
      + dt_fac*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*int_factor/rden
    ! by, vel
    DO l=1,3
      jac_loc(7,l+1)%m(jr,jc) = jac_loc(7,l+1)%m(jr,jc) &
        + basis_vals(jr)*dt_fac*basis_vals(jc)*dby(l)*int_factor/rden &
        + basis_vals(jr)*dt_fac*by*basis_grads(l,jc)*int_factor/rden
    END DO
    jac_loc(7,3)%m(jr,jc) = jac_loc(7,3)%m(jr,jc) &
      + basis_vals(jr)*dt_fac*psi_cross(2)*int_factor/rden
    ! by, psi
    jac_loc(7,6)%m(jr,jc) = jac_loc(7,6)%m(jr,jc) &
      + basis_vals(jr)*dt_fac*dvel_cross(2)*int_factor/rden
    ! by, by
    jac_loc(7,7)%m(jr,jc) = jac_loc(7,7)%m(jr,jc) &
      + basis_vals(jr)*basis_vals(jc)*int_factor/rden &
      + basis_vals(jr)*dt_fac*DOT_PRODUCT(basis_grads(:,jc),vel)*int_factor/rden &
      + basis_vals(jr)*dt_fac*basis_vals(jc)*(dvel(1,1)+dvel(3,3))*int_factor/rden &
      + dt_fac*eta(1)*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*int_factor/rden
  END DO
END DO
!---Momentum, vel-n coupling, which takes a different form in linear mode
IF(linear)THEN
  DO jc=1,nce
    DO jr=1,nce
      DO l=1,3
        jac_loc(l+1,1)%m(jr,jc) = jac_loc(l+1,1)%m(jr,jc) &
          + dt_fac*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*dT(l)*int_factor*rmul/(m_i*n) &
          + dt_fac*basis_vals(jr)*2.d0*k_boltz*T*basis_grads(l,jc)*int_factor*rmul/(m_i*n) &
          + basis_vals(jr)*basis_vals(jc)*dt_fac*DOT_PRODUCT(vel,dvel(l,:))*int_factor*rmul/n
      END DO
    END DO
  END DO
ELSE
  DO jc=1,nce
    DO jr=1,nce
      DO l=1,3
        jac_loc(l+1,1)%m(jr,jc) = jac_loc(l+1,1)%m(jr,jc) &
          ! Maxwell stress (1/n factor)
          - dt_fac*basis_vals(jc)*DOT_PRODUCT(basis_grads(:,jr),btmp)*btmp(l)*int_factor*rmul/(mu0*m_i*n**2.d0) &
          + dt_fac*basis_vals(jc)*basis_grads(l,jr)*DOT_PRODUCT(btmp,btmp)*int_factor*rmul/(2.d0*mu0*m_i*n**2.d0) &
          ! velocity dissipation terms
          - dt_fac*nu*basis_vals(jr)*DOT_PRODUCT(basis_grads(:,jc),dvel(l,:))*int_factor*rmul/n &
          + dt_fac*nu*basis_vals(jr)*basis_vals(jc)*DOT_PRODUCT(dn,dvel(l,:))*int_factor*rmul/n**2
        IF(incomp)THEN
          ! Pressure (1/n factor)
          jac_loc(l+1,1)%m(jr,jc) = jac_loc(l+1,1)%m(jr,jc) &
            + dt_fac*basis_vals(jc)*T*basis_grads(l,jr)*int_factor*rmul/(m_i*n**2)
        ELSE
          jac_loc(l+1,1)%m(jr,jc) = jac_loc(l+1,1)%m(jr,jc) &
            ! grad T terms
            + dt_fac*basis_vals(jr)*2.d0*k_boltz*T*basis_grads(l,jc)*int_factor*rmul/(m_i*n) &
            - dt_fac*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*dn(l)*T*int_factor*rmul/(m_i*n**2.d0) &
            ! Maxwell stress (grad(n) terms)
            - dt_fac*basis_vals(jr)*DOT_PRODUCT(basis_grads(:,jc),btmp)*btmp(l)*int_factor*rmul/(mu0*m_i*n**2.d0) &
            + dt_fac*basis_vals(jr)*2.d0*basis_vals(jc)*DOT_PRODUCT(dn,btmp)*btmp(l)*int_factor*rmul/(mu0*m_i*n**3.d0) &
            + dt_fac*basis_vals(jr)*basis_grads(l,jc)*DOT_PRODUCT(btmp,btmp)*int_factor*rmul/(2.d0*mu0*m_i*n**2.d0) &
            - dt_fac*basis_vals(jr)*basis_vals(jc)*dn(l)*DOT_PRODUCT(btmp,btmp)*int_factor*rmul/(mu0*m_i*n**3.d0)
        END IF
      END DO
    END DO
  END DO
END IF
IF(cyl_flag)THEN
  !---Extra terms in cylindrical geometry (hoop stress, curvature and 1/R factors)
  DO jc=1,nce
    db_psi = cross_product(basis_grads(:,jc),[0.d0,1.d0/rden,0.d0]) ! Perturbed B from delta_psi
    db_by = [0.d0,basis_vals(jc)/rden,0.d0] ! Perturbed B from delta_by
    DO jr=1,nce
      ! n, vel
      jac_loc(1,2)%m(jr,jc) = jac_loc(1,2)%m(jr,jc) &
        + basis_vals(jr)*dt_fac*n*basis_vals(jc)*int_factor
      ! vel, n
      jac_loc(2,1)%m(jr,jc) = jac_loc(2,1)%m(jr,jc) &
        - dt_fac*basis_vals(jr)*basis_vals(jc)*(btmp(2)**2-btmp(1)**2-btmp(3)**2)*int_factor/(2.d0*mu0*m_i*n**2)
      ! Curvature term of the incompressible pressure gradient (integrated by parts)
      IF(incomp)jac_loc(2,1)%m(jr,jc) = jac_loc(2,1)%m(jr,jc) &
        + dt_fac*basis_vals(jr)*basis_vals(jc)*T*int_factor/(m_i*n**2)
      jac_loc(3,1)%m(jr,jc) = jac_loc(3,1)%m(jr,jc) &
        + dt_fac*basis_vals(jr)*basis_vals(jc)*(btmp(1)*btmp(2))*int_factor/(mu0*m_i*n**2)
      ! vel, vel
      jac_loc(2,2)%m(jr,jc) = jac_loc(2,2)%m(jr,jc) &
        + dt_fac*nu*basis_vals(jr)*basis_vals(jc)*int_factor/rden
      jac_loc(2,3)%m(jr,jc) = jac_loc(2,3)%m(jr,jc) &
        - dt_fac*basis_vals(jr)*2.d0*vel(2)*basis_vals(jc)*int_factor
      jac_loc(3,2)%m(jr,jc) = jac_loc(3,2)%m(jr,jc) &
        + dt_fac*basis_vals(jr)*basis_vals(jc)*vel(2)*int_factor
      jac_loc(3,3)%m(jr,jc) = jac_loc(3,3)%m(jr,jc) &
        + dt_fac*basis_vals(jr)*basis_vals(jc)*vel(1)*int_factor &
        + dt_fac*nu*basis_vals(jr)*basis_vals(jc)*int_factor/rden
      ! vel, psi
      jac_loc(2,6)%m(jr,jc) = jac_loc(2,6)%m(jr,jc) &
        + dt_fac*basis_vals(jr)*(btmp(2)*db_psi(2)-btmp(1)*db_psi(1)-btmp(3)*db_psi(3))*int_factor/(m_i*n*mu0)
      jac_loc(3,6)%m(jr,jc) = jac_loc(3,6)%m(jr,jc) &
        - dt_fac*basis_vals(jr)*(btmp(1)*db_psi(2)+btmp(2)*db_psi(1))*int_factor/(m_i*n*mu0)
      ! vel, by
      jac_loc(2,7)%m(jr,jc) = jac_loc(2,7)%m(jr,jc) &
        + dt_fac*basis_vals(jr)*(btmp(2)*db_by(2)-btmp(1)*db_by(1)-btmp(3)*db_by(3))*int_factor/(m_i*n*mu0)
      jac_loc(3,7)%m(jr,jc) = jac_loc(3,7)%m(jr,jc) &
        - dt_fac*basis_vals(jr)*(btmp(1)*db_by(2)+btmp(2)*db_by(1))*int_factor/(m_i*n*mu0)
      ! by, vel
      jac_loc(7,2)%m(jr,jc) = jac_loc(7,2)%m(jr,jc) &
        - basis_vals(jr)*dt_fac*by*basis_vals(jc)*int_factor/rden**2
      jac_loc(7,3)%m(jr,jc) = jac_loc(7,3)%m(jr,jc) &
        - basis_vals(jr)*dt_fac*basis_vals(jc)*dpsi(3)*int_factor/rden**2 &
        ! Background field (only B_0z allowed)
        - basis_vals(jr)*dt_fac*B_0(3)*basis_grads(3,jc)*int_factor
      ! by, psi
      jac_loc(7,6)%m(jr,jc) = jac_loc(7,6)%m(jr,jc) &
        - basis_vals(jr)*dt_fac*vel(2)*basis_grads(3,jc)*int_factor/rden**2
      ! by, by
      jac_loc(7,7)%m(jr,jc) = jac_loc(7,7)%m(jr,jc) &
        - dt_fac*basis_vals(jr)*vel(1)*basis_vals(jc)*int_factor/rden**2
    END DO
  END DO
ELSE
  !---Background field terms in Cartesian geometry
  DO jc=1,nce
    DO jr=1,nce
      ! by, vel: B_0y*div(delta v)
      DO l=1,3
        jac_loc(7,l+1)%m(jr,jc) = jac_loc(7,l+1)%m(jr,jc) &
          + basis_vals(jr)*dt_fac*B_0(2)*basis_grads(l,jc)*int_factor
      END DO
      ! by, vel: -B_0.grad(delta v_y)
      jac_loc(7,3)%m(jr,jc) = jac_loc(7,3)%m(jr,jc) &
        - basis_vals(jr)*dt_fac*DOT_PRODUCT(B_0,basis_grads(:,jc))*int_factor
    END DO
  END DO
END IF
END SUBROUTINE jac_add_shared
!---------------------------------------------------------------------------
!> Add Jacobian contributions specific to the compressible system (temperature couplings)
!---------------------------------------------------------------------------
SUBROUTINE jac_add_comp(self,jac_loc,pt,basis_vals,basis_grads,dt_fac)
CLASS(oft_xmhd_2d_sim), INTENT(IN) :: self !< Simulation object
TYPE(oft_local_mat), INTENT(INOUT) :: jac_loc(:,:) !< Local (cell) Jacobian
TYPE(xmhd_2d_qpt), INTENT(IN) :: pt !< Quadrature point state
REAL(r8), INTENT(IN) :: basis_vals(:) !< Basis values [nce]
REAL(r8), INTENT(IN) :: basis_grads(:,:) !< Basis gradients [3,nce]
REAL(r8), INTENT(IN) :: dt_fac !< Time step factor
LOGICAL :: cyl_flag,linear
INTEGER(i4) :: jr,jc,l,nce
REAL(r8) :: n,dn(3),vel(3),div_vel,T,dT(3),rmul,int_factor,m_i,chi,gamma
!---Local copies (avoids reloads through possible aliasing with jac_loc pointers)
cyl_flag=self%cyl_flag; linear=self%linear
n=pt%n; dn=pt%dn; vel=pt%vel; div_vel=pt%div_vel; T=pt%T; dT=pt%dT
rmul=pt%rmul; int_factor=pt%int_factor; m_i=pt%m_i; chi=pt%chi; gamma=pt%gamma
nce=SIZE(basis_vals)
!---Terms common to every geometry and mode
DO jc=1,nce
  DO jr=1,nce
    ! vel, T
    DO l=1,3
      jac_loc(l+1,5)%m(jr,jc) = jac_loc(l+1,5)%m(jr,jc) &
        + dt_fac*basis_vals(jr)*2.d0*k_boltz*dn(l)*basis_vals(jc)*int_factor*rmul/(m_i*n) &
        + dt_fac*basis_vals(jr)*2.d0*k_boltz*basis_grads(l,jc)*int_factor*rmul/m_i
    END DO
    ! T, vel
    DO l=1,3
      jac_loc(5,l+1)%m(jr,jc) = jac_loc(5,l+1)%m(jr,jc) &
        + basis_vals(jr)*dt_fac*basis_vals(jc)*dT(l)*int_factor*rmul/(gamma-1) & ! delta_u dot grad(T)
        + basis_vals(jr)*dt_fac*T*basis_grads(l,jc)*int_factor*rmul ! T*div(delta_u)
    END DO
    ! T, T
    jac_loc(5,5)%m(jr,jc) = jac_loc(5,5)%m(jr,jc) &
      + basis_vals(jr)*basis_vals(jc)*int_factor*rmul/(gamma-1) & ! delta_T
      + basis_vals(jr)*dt_fac*DOT_PRODUCT(vel,basis_grads(:,jc))*int_factor*rmul/(gamma-1) & ! u dot grad(delta_T)
      + basis_vals(jr)*dt_fac*basis_vals(jc)*div_vel*int_factor*rmul & ! delta_T*div(u)
      + dt_fac*chi*DOT_PRODUCT(basis_grads(:,jc),basis_grads(:,jr))*int_factor*rmul & ! Conduction
      - dt_fac*basis_vals(jr)*chi*DOT_PRODUCT(dn,basis_grads(:,jc))*int_factor*rmul/n
  END DO
END DO
!---T, n coupling, which takes a different form in linear mode
IF(linear)THEN
  DO jc=1,nce
    DO jr=1,nce
      ! The 1/n dependence of the conduction term is evaluated with the steady-state
      ! temperature equation for the equilibrium
      jac_loc(5,1)%m(jr,jc) = jac_loc(5,1)%m(jr,jc) &
        + dt_fac*basis_vals(jr)*basis_vals(jc)*DOT_PRODUCT(vel,dT)*int_factor*rmul/(n*(gamma-1)) &
        + dt_fac*basis_vals(jr)*basis_vals(jc)*T*div_vel*int_factor*rmul/n &
        + dt_fac*basis_vals(jc)*chi*DOT_PRODUCT(basis_grads(:,jr),dT)*int_factor*rmul/n &
        - dt_fac*basis_vals(jr)*chi*basis_vals(jc)*DOT_PRODUCT(dn,dT)*int_factor*rmul/(n**2)
    END DO
  END DO
ELSE
  DO jc=1,nce
    DO jr=1,nce
      jac_loc(5,1)%m(jr,jc) = jac_loc(5,1)%m(jr,jc) &
        - dt_fac*chi*basis_vals(jr)*DOT_PRODUCT(basis_grads(:,jc),dT)*int_factor*rmul/n & ! grad(delta_n)
        + dt_fac*chi*basis_vals(jr)*basis_vals(jc)*DOT_PRODUCT(dn,dT)*int_factor*rmul/(n**2) ! delta_n
    END DO
  END DO
END IF
!---Extra T, vel term in cylindrical geometry
IF(cyl_flag)THEN
  DO jc=1,nce
    DO jr=1,nce
      jac_loc(5,2)%m(jr,jc) = jac_loc(5,2)%m(jr,jc) &
        + basis_vals(jr)*dt_fac*T*basis_vals(jc)*int_factor
    END DO
  END DO
END IF
END SUBROUTINE jac_add_comp
!---------------------------------------------------------------------------
!> Add Jacobian contributions specific to the incompressible system (pressure couplings)
!---------------------------------------------------------------------------
SUBROUTINE jac_add_incomp(self,jac_loc,pt,basis_vals,basis_grads,basis_vals_p,dt_fac)
CLASS(oft_xmhd_2d_sim), INTENT(IN) :: self !< Simulation object
TYPE(oft_local_mat), INTENT(INOUT) :: jac_loc(:,:) !< Local (cell) Jacobian
TYPE(xmhd_2d_qpt), INTENT(IN) :: pt !< Quadrature point state
REAL(r8), INTENT(IN) :: basis_vals(:) !< Basis values [nce]
REAL(r8), INTENT(IN) :: basis_grads(:,:) !< Basis gradients [3,nce]
REAL(r8), INTENT(IN) :: basis_vals_p(:) !< Pressure basis values [nce_p]
REAL(r8), INTENT(IN) :: dt_fac !< Time step factor
LOGICAL :: cyl_flag
INTEGER(i4) :: jr,jc,l
REAL(r8) :: n,rmul,int_factor,m_i
!---Local copies (avoids reloads through possible aliasing with jac_loc pointers)
cyl_flag=self%cyl_flag
n=pt%n; rmul=pt%rmul; int_factor=pt%int_factor; m_i=pt%m_i
! p, vel (divergence constraint)
DO jr=1,SIZE(basis_vals_p)
  DO jc=1,SIZE(basis_vals)
    DO l=1,3
      jac_loc(5,l+1)%m(jr,jc) = jac_loc(5,l+1)%m(jr,jc) &
        + basis_vals_p(jr)*basis_grads(l,jc)*int_factor*rmul
    END DO
    IF(cyl_flag)THEN
      jac_loc(5,2)%m(jr,jc) = jac_loc(5,2)%m(jr,jc) &
        + basis_vals_p(jr)*basis_vals(jc)*int_factor
    END IF
  END DO
END DO
! vel, p (pressure gradient)
DO jr=1,SIZE(basis_vals)
  DO jc=1,SIZE(basis_vals_p)
    DO l=1,3
      jac_loc(l+1,5)%m(jr,jc) = jac_loc(l+1,5)%m(jr,jc) &
        - dt_fac*basis_vals_p(jc)*basis_grads(l,jr)*int_factor*rmul/(m_i*n)
    END DO
    IF(cyl_flag)THEN
      jac_loc(2,5)%m(jr,jc) = jac_loc(2,5)%m(jr,jc) &
        - dt_fac*basis_vals_p(jc)*basis_vals(jr)*int_factor/(m_i*n)
    END IF
  END DO
END DO
END SUBROUTINE jac_add_incomp
!---------------------------------------------------------------------------
!> Apply density scaling and boundary conditions to a local (cell) Jacobian
!---------------------------------------------------------------------------
SUBROUTINE jac_finish_local(self,jac_loc,cell_dofs,T_dofs)
CLASS(oft_xmhd_2d_sim), INTENT(INOUT) :: self !< Simulation object
TYPE(oft_local_mat), INTENT(INOUT) :: jac_loc(:,:) !< Local (cell) Jacobian
INTEGER(i4), INTENT(IN) :: cell_dofs(:) !< Global indices of cell DOFs
INTEGER(i4), INTENT(IN) :: T_dofs(:) !< Global indices of cell DOFs for field 5 (temperature or pressure)
INTEGER(i4) :: jr
DO jr=1,7
  jac_loc(1,jr)%m=jac_loc(1,jr)%m/self%den_scale
  jac_loc(jr,1)%m=jac_loc(jr,1)%m*self%den_scale
END DO
CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%T_bc(T_dofs),5)
CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%n_bc(cell_dofs),1)
CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%velx_bc(cell_dofs),2)
CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%vely_bc(cell_dofs),3)
CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%velz_bc(cell_dofs),4)
CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%psi_bc(cell_dofs),6)
CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%by_bc(cell_dofs),7)
END SUBROUTINE jac_finish_local
!---------------------------------------------------------------------------
!> Set Dirichlet rows and assemble the global Jacobian
!---------------------------------------------------------------------------
SUBROUTINE jac_finalize(self,T_fe)
CLASS(oft_xmhd_2d_sim), INTENT(INOUT) :: self !< Simulation object
CLASS(oft_scalar_bfem), INTENT(INOUT) :: T_fe !< Lagrange representation for field 5 (temperature or pressure)
CLASS(oft_vector), POINTER :: tmp
IF(oft_debug_print(2))WRITE(*,'(4X,A)')'Setting BCs'
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%n_bc,1)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%velx_bc,2)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%vely_bc,3)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%velz_bc,4)
CALL fem_dirichlet_diag(T_fe,self%jacobian,self%T_bc,5)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%psi_bc,6)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%by_bc,7)
!---Assemble
CALL self%fe_rep%vec_create(tmp)
CALL self%jacobian%assemble(tmp)
CALL tmp%delete
DEALLOCATE(tmp)
END SUBROUTINE jac_finalize

!---------------------------------------------------------------------------
!> Update matrix-free Jacobian on all levels with new solution
!---------------------------------------------------------------------------
subroutine mfnk_update(uin)
class(oft_vector), target, intent(inout) :: uin !< Current field
IF(oft_debug_print(1))write(*,*)'Updating 2D MUG MF-Jacobian'
CALL current_sim%mf_mat%update(uin)
END SUBROUTINE mfnk_update
!---------------------------------------------------------------------------
!> Update Jacobian matrices on all levels with new solution
!---------------------------------------------------------------------------
subroutine update_jacobian(uin)
class(oft_vector), target, intent(inout) :: uin !< Current solution
IF(oft_debug_print(1))write(*,*)'Updating 2D MUG approximate Jacobian'
CALL build_jacobian(current_sim,uin)
END SUBROUTINE update_jacobian
!---------------------------------------------------------------------------
!> Setup composite FE representation and ML environment
!---------------------------------------------------------------------------
subroutine setup(self,mg_mesh_in, order, fe_rep_in)
class(oft_xmhd_2d_sim), intent(inout), target :: self
CLASS(multigrid_mesh), TARGET, intent(in) :: mg_mesh_in
integer(i4), intent(in) :: order
CLASS(oft_scalar_bfem), TARGET,optional, intent(in) :: fe_rep_in
integer(i4) :: i,j, ierr,io_unit, cond_ind, coil_ind, type, order_p
LOGICAL, ALLOCATABLE :: vert_flag(:),edge_flag(:), boundary_flag(:)
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs

mg_mesh=>mg_mesh_in
mesh=>mg_mesh_in%smesh

IF(ASSOCIATED(self%fe_rep))CALL oft_abort("Setup can only be called once","setup",__FILE__)
IF(ASSOCIATED(oft_blagrange))CALL oft_abort("FE space already built","setup",__FILE__)

!---Look for XML defintion elements
IF(ASSOCIATED(oft_env%xml))THEN
  CALL xml_get_element(oft_env%xml,"xmhd2d",self%xml_root,ierr)
  IF(ierr==0)THEN
    !---Look for pre node
    CALL xml_get_element(self%xml_root,"pre",self%xml_pre_def,ierr)
  END IF
END IF

!---Setup FE representation
IF (PRESENT(fe_rep_in)) THEN
  IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Using user-provided FE representation'
  oft_blagrange => fe_rep_in
ELSE
  IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Building lagrange FE space'
  CALL oft_lag_setup(mg_mesh,order,ML_blag_obj=ML_oft_blagrange,minlev=-1)
  IF(.NOT.oft_2D_lagrange_cast(oft_blagrange,ML_oft_blagrange%current_level))CALL oft_abort("Invalid lagrange FE object","setup",__FILE__)
END IF

IF (self%incomp) THEN
  order_p = order-1
  IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Building pressure lagrange FE space'
  CALL oft_lag_setup(mg_mesh,order_p, ML_blag_obj=ML_oft_blagrange_p,minlev=-1)
  IF(.NOT.oft_2D_lagrange_cast(oft_blagrange_p,ML_oft_blagrange_p%current_level))CALL oft_abort("Invalid lagrange FE object","setup",__FILE__)
END IF

!---Build composite FE definition for solution field
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Creating FE type'
ALLOCATE(self%fe_rep)
self%fe_rep%nfields=7
ALLOCATE(self%fe_rep%fields(self%fe_rep%nfields))
ALLOCATE(self%fe_rep%field_tags(self%fe_rep%nfields))
self%fe_rep%fields(1)%fe=>oft_blagrange
self%fe_rep%field_tags(1)='n'
self%fe_rep%fields(2)%fe=>oft_blagrange
self%fe_rep%field_tags(2)='velx'
self%fe_rep%fields(3)%fe=>oft_blagrange
self%fe_rep%field_tags(3)='vely'
self%fe_rep%fields(4)%fe=>oft_blagrange
self%fe_rep%field_tags(4)='velz'
IF(self%incomp)THEN
  self%fe_rep%fields(5)%fe=>oft_blagrange_p
  self%fe_rep%fields(5)%fe%type = 2
  self%fe_rep%field_tags(5)='p'
ELSE
  self%fe_rep%fields(5)%fe=>oft_blagrange
  self%fe_rep%field_tags(5)='T'
END IF
self%fe_rep%fields(6)%fe=>oft_blagrange
self%fe_rep%field_tags(6)='psi'
self%fe_rep%fields(7)%fe=>oft_blagrange
self%fe_rep%field_tags(7)='by'
!---Create solution vector
CALL self%fe_rep%vec_create(self%u)
CALL self%fe_rep%vec_create(self%u0)
!---Create vacuum poloidal flux offset for the psi field (initialized to zero)
CALL oft_blagrange%vec_create(self%psi_vac)
CALL self%psi_vac%set(0.d0)
!---Create Jacobian matrix
CALL self%fe_rep%mat_create(self%jacobian)

!-- Set material properties not already set
IF (.NOT. ALLOCATED(self%nu)) THEN
  ALLOCATE(self%nu(mesh%nreg))
  self%nu = -1.d0
END IF
IF (.NOT. ALLOCATED(self%D_diff)) THEN
  ALLOCATE(self%D_diff(mesh%nreg))
  self%D_diff = -1.d0
END IF
IF (.NOT. ALLOCATED(self%chi)) THEN
  ALLOCATE(self%chi(mesh%nreg))
  self%chi = -1.d0
END IF
IF (.NOT. ALLOCATED(self%m_i)) THEN
  ALLOCATE(self%m_i(mesh%nreg))
  self%m_i = proton_mass
END IF
IF (.NOT. ALLOCATED(self%eta)) THEN
  ALLOCATE(self%eta(mesh%nreg, 2))
  self%eta = -1.d0
END IF
IF (.NOT. ALLOCATED(self%gamma)) THEN
  ALLOCATE(self%gamma(mesh%nreg))
  self%gamma = -1.d0
END IF
IF (.NOT. ALLOCATED(self%ignore_rmask)) THEN
  ALLOCATE(self%ignore_rmask(mesh%nreg))
  self%ignore_rmask = .FALSE.
END IF
! Set boundary conditions not alreadys set
CALL self%setup_bc()

!---------------------------------------------------------------------------
! Allocate nonlinear function object
!---------------------------------------------------------------------------
ALLOCATE(self%nlfun)
self%nlfun%parent_sim=>self
end subroutine setup

!---------------------------------------------------------------------------
!> Set any boundary conditions not already set with default values
!> Note that default BCs on pressure 
!---------------------------------------------------------------------------
subroutine setup_bc(self)
class(oft_xmhd_2d_sim), intent(inout) :: self
INTEGER(i4) :: i
IF(.NOT.ASSOCIATED(self%n_bc))THEN
  !If incompressible the density is not evolved
  IF(self%incomp)THEN
    ALLOCATE(self%n_bc(oft_blagrange%ne))
    self%n_bc=.TRUE.
  ELSE
    self%n_bc=>oft_blagrange%global%gbe
  END IF
END IF
IF(.NOT.ASSOCIATED(self%velx_bc))self%velx_bc=>oft_blagrange%global%gbe
IF(.NOT.ASSOCIATED(self%vely_bc))self%vely_bc=>oft_blagrange%global%gbe
IF(.NOT.ASSOCIATED(self%velz_bc))self%velz_bc=>oft_blagrange%global%gbe
IF (self%incomp) THEN
  !Default pressure BC is freezing one node to fix gauge
  !Only works for single MHD region simulations, must be overwritten if several MHD regions exist
  IF(.NOT.ASSOCIATED(self%T_bc))THEN
    ALLOCATE(self%T_bc(oft_blagrange_p%ne))
    self%T_bc = .FALSE.
    DO i=1,oft_blagrange_p%map%nslice
      IF(oft_blagrange_p%global%le(oft_blagrange_p%map%slice(i))==1)THEN
        self%T_bc(oft_blagrange_p%map%slice(i)) = .TRUE.
        EXIT
      END IF
    END DO
  END IF
ELSE
  IF(.NOT.ASSOCIATED(self%T_bc))self%T_bc=>oft_blagrange%global%gbe
END IF
IF(.NOT.ASSOCIATED(self%psi_bc))self%psi_bc=>oft_blagrange%global%gbe
IF(.NOT.ASSOCIATED(self%by_bc))self%by_bc=>oft_blagrange%global%gbe
end subroutine setup_bc

!---------------------------------------------------------------------------
!> Save xMHD solution state to a restart file
!---------------------------------------------------------------------------
subroutine rst_save(self,u,t,dt,filename,path)
class(oft_xmhd_2d_sim), intent(inout) :: self
class(oft_vector), pointer, intent(inout) :: u !< Solution to save
real(r8), intent(in) :: t !< Current solution time
real(r8), intent(in) :: dt !< Current timestep
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to store solution vector in file
DEBUG_STACK_PUSH
CALL self%fe_rep%vec_save(u,filename,path)
IF(oft_env%head_proc)THEN
  CALL hdf5_write(t,filename,'t')
  CALL hdf5_write(dt,filename,'dt')
END IF
DEBUG_STACK_POP
end subroutine rst_save
!---------------------------------------------------------------------------
!> Save the fixed background fields (vacuum flux `psi_vac` and background field `B_0`)
!! to an existing restart file
!!
!! The run drivers add these to the step-0 restart file, so that @ref xmhd_2d_plot can
!! include them in the plotted totals.
!---------------------------------------------------------------------------
subroutine save_background(self,filename)
class(oft_xmhd_2d_sim), intent(inout) :: self
character(LEN=*), intent(in) :: filename !< Name of restart file (must already exist)
DEBUG_STACK_PUSH
CALL oft_blagrange%vec_save(self%psi_vac,filename,'psi_vac',append=.TRUE.)
IF(oft_env%head_proc)CALL hdf5_write(self%B_0,filename,'B_0')
DEBUG_STACK_POP
end subroutine save_background

!---------------------------------------------------------------------------
!> Plot the current solution state using XDMF/HDF5 output
!! Runtime options are set in the main input file using the group `xmhd_plot_options`.
!! **Option group:** `xmhd_plot_options`
!! |  Option    |  Description  | Type [dim] |
!! |------------|---------------|------------|
!! | `rst_start=0`      | First restart file to read | int |
!! | `rst_end=2000`     | Last restart file to read | int |
!!
!! The plotted poloidal flux `psi` and magnetic field (`B`, or `B*R` in cylindrical geometry)
!! are the total values, including the vacuum flux (`psi_vac`) and the background field (`B_0`).
!! Both are read from the step-0 restart file (`xmhd2d_00000.rst`), where the run drivers save
!! them. If they are not found there, a warning is printed and they are not included.
!------------------------------------------------------------------------------
subroutine xmhd_2d_plot(self)
class(oft_xmhd_2d_sim), intent(inout) :: self
class(oft_vector), pointer :: ux,uy,uz,v_lag
type(oft_lag_bginterp) :: grad_psi
class(oft_matrix), pointer :: lmop => NULL()
real(r8), pointer :: plot_vals(:),plot_vec(:,:),plot_u0(:),psi_bg(:)
real(r8) :: b0_plot(3)
integer(i4) :: bg_err
logical :: bg_found
class(oft_vector), pointer :: psi_vac_plot
!---Solver objects
CLASS(oft_solver), POINTER :: lminv => NULL()
class(oft_vector), pointer :: u,v,up, u0
class(oft_matrix), pointer :: lap_mat => NULL()
INTEGER(i4):: rst_start=0
INTEGER(i4):: rst_end=2000
INTEGER(i4) :: rst_cur, rst_tmp, ierr, io_stat, io_unit
CHARACTER(LEN=OFT_PATH_SLEN) :: file_tmp
LOGICAL :: rst_exist
real(r8) :: t
TYPE(xdmf_plot_file) :: xdmf_plot, xdmf_plot_p
character(LEN=XMHD_RST_LEN) :: rst_char
namelist/xmhd_plot_options/rst_start,rst_end
!---------------------------------------------------------------------------
! Read plotting options
!---------------------------------------------------------------------------
open(NEWUNIT=io_unit,FILE=oft_env%ifile)
read(io_unit,xmhd_plot_options,IOSTAT=ierr)
close(io_unit)
!---------------------------------------------------------------------------
! Create solver fields
!---------------------------------------------------------------------------
call self%fe_rep%vec_create(u)
call self%fe_rep%vec_create(up)
call self%fe_rep%vec_create(v)
call self%fe_rep%vec_create(u0)
!---
call oft_blagrange%vec_create(grad_psi%u)
call oft_blagrange%vec_create(ux)
call oft_blagrange%vec_create(uy)
call oft_blagrange%vec_create(uz)
call oft_blagrange%vec_create(v_lag)
NULLIFY(lmop)
call oft_blag_getmop(oft_blagrange,lmop)
CALL create_cg_solver(lminv)
lminv%A=>lmop
lminv%its=-2
CALL create_diag_pre(lminv%pre)
ALLOCATE(plot_vec(3,v_lag%n))
NULLIFY(plot_vals,plot_u0,psi_bg)
CALL grad_psi%setup(oft_blagrange)
!---------------------------------------------------------------------------
! Background fields (psi_vac and B_0) saved in the step-0 restart file
!---------------------------------------------------------------------------
100 FORMAT (I XMHD_RST_LEN.XMHD_RST_LEN)
WRITE(rst_char,100)0
file_tmp='xmhd2d_'//rst_char//'.rst'
CALL oft_blagrange%vec_create(psi_vac_plot)
CALL psi_vac_plot%set(0.d0)
b0_plot=0.d0
bg_found=oft_file_exist(TRIM(file_tmp))
CALL oft_mpi_barrier(ierr)
IF(bg_found)THEN
  CALL oft_blagrange%vec_load(psi_vac_plot,TRIM(file_tmp),'psi_vac',err_flag=bg_err)
  bg_found=(bg_err==0)
END IF
IF(bg_found)CALL hdf5_read(b0_plot,TRIM(file_tmp),'B_0',success=bg_found)
IF(.NOT.bg_found)THEN
  CALL psi_vac_plot%set(0.d0)
  b0_plot=0.d0
  IF(oft_env%head_proc)CALL oft_warn('Background fields (psi_vac, B_0) not found in "'//TRIM(file_tmp)// &
    '"; if used, they are not included in the plotted psi and B')
END IF
!---------------------------------------------------------------------------
! Background flux added to the plotted psi: the vacuum flux plus the flux of the
! in-plane background field (B_0z*x - B_0x*z planar, B_0Z*R^2/2 cylindrical)
!---------------------------------------------------------------------------
CALL psi_vac_plot%get_local(psi_bg)
CALL psi_vac_plot%delete
DEALLOCATE(psi_vac_plot)
BLOCK
INTEGER(i4) :: i,j
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs
REAL(r8) :: f(3),pt(3)
ALLOCATE(cell_dofs(oft_blagrange%nce))
DO i=1,mesh%nc
  CALL oft_blagrange%ncdofs(i,cell_dofs)
  DO j=1,oft_blagrange%nce
    CALL oft_blag_npos(oft_blagrange,i,j,f)
    pt=mesh%log2phys(i,f)
    IF(self%cyl_flag)THEN
      plot_vec(1,cell_dofs(j))=b0_plot(3)*pt(1)**2/2.d0
    ELSE
      plot_vec(1,cell_dofs(j))=b0_plot(3)*pt(1)-b0_plot(1)*pt(2)
    END IF
  END DO
END DO
DEALLOCATE(cell_dofs)
END BLOCK
psi_bg=psi_bg+plot_vec(1,:)

!---------------------------------------------------------------------------
! Pass 1: main fields (n, V, psi, B, and T if compressible) on the full-order mesh
!---------------------------------------------------------------------------
CALL xdmf_plot%setup("xmhd_2d")
CALL mesh%setup_io(xdmf_plot,oft_blagrange%order)

!---------------------------------------------------------------------------
! If linear, extract equilibrium fields
!---------------------------------------------------------------------------
IF(self%linear)THEN
  file_tmp='xmhd2d_eq.rst'
  rst_exist=oft_file_exist(TRIM(file_tmp))
  CALL oft_mpi_barrier(ierr)
  IF(.NOT.rst_exist)THEN
    CALL oft_abort('Equilibrium restart file does not exist.','xmhd_plot',__FILE__)
  ELSE
    CALL hdf5_read(t,file_tmp,'t')
    CALL self%rst_load(u0,file_tmp, 'U')
  END IF
END IF

!---------------------------------------------------------------------------
! Loop through rst files
!---------------------------------------------------------------------------
rst_cur = rst_start
DO
  IF(rst_cur > rst_end) EXIT
  ! Determine file name
  WRITE(rst_char,100)rst_cur
  READ(rst_char,100,IOSTAT=io_stat)rst_tmp
  IF((io_stat/=0).OR.(rst_tmp/=rst_cur))CALL oft_abort("Step count exceeds format width", "xmhd_plot", __FILE__)
  file_tmp='xmhd2d_'//rst_char//'.rst'
  !Read file
  rst_exist=oft_file_exist(TRIM(file_tmp))
  CALL oft_mpi_barrier(ierr)
  IF(.NOT.rst_exist)THEN
    EXIT !CALL oft_abort('Restart file does not exist.','xmhd_plot',__FILE__)
  ELSE
    CALL hdf5_read(t,file_tmp,'t')
    CALL self%rst_load(u,file_tmp, 'U')
  END IF
  IF(self%linear)THEN
    CALL u%add(1.d0, 1.d0, u0)
  END IF
  !Plot data
  CALL xdmf_plot%add_timestep(t)

  !Plot density
  CALL u%get_local(plot_vals,1)
  plot_vals = plot_vals*self%den_scale
  CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot,'n')
  !Plot velocity
  CALL u%get_local(plot_vals,2)
  plot_vec(1,:)=plot_vals
  CALL u%get_local(plot_vals,3)
  plot_vec(3,:)=plot_vals
  CALL u%get_local(plot_vals,4)
  plot_vec(2,:)=plot_vals
  CALL mesh%save_vertex_vector(plot_vec,xdmf_plot,'V')
  !Plot temperature (compressible only)
  IF(.NOT.self%incomp)THEN
    NULLIFY(plot_vals)
    CALL u%get_local(plot_vals,5)
    CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot,'T')
  END IF
  !Plot total psi (including vacuum flux and background field)
  NULLIFY(plot_vals)
  CALL u%get_local(plot_vals,6)
  plot_vals=plot_vals+psi_bg
  CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot,'psi')
  !Plot B
  CALL grad_psi%u%restore_local(plot_vals)
  CALL grad_psi%setup(oft_blagrange)
  CALL oft_blag_vproject(oft_blagrange,grad_psi,ux,uy,uz)
  CALL v_lag%set(0.d0)
  CALL lminv%apply(v_lag,ux)
  CALL ux%add(0.d0,1.d0,v_lag)
  CALL v_lag%set(0.d0)
  CALL lminv%apply(v_lag,uy)
  CALL uy%add(0.d0,1.d0,v_lag)
  CALL uy%get_local(plot_vals)
  plot_vec(1,:)=-plot_vals
  CALL ux%get_local(plot_vals)
  plot_vec(2,:)=plot_vals
  CALL u%get_local(plot_vals,7)
  plot_vec(3,:)=plot_vals
  ! Out-of-plane background field (only valid in planar geometry)
  IF(.NOT.self%cyl_flag)plot_vec(3,:)=plot_vec(3,:)+b0_plot(2)
  IF(self%cyl_flag)THEN
    CALL mesh%save_vertex_vector(plot_vec,xdmf_plot,'B*R')
  ELSE
    CALL mesh%save_vertex_vector(plot_vec,xdmf_plot,'B')
  END IF

  !Move to next file
  rst_cur = rst_cur + self%rst_freq
END DO
!---------------------------------------------------------------------------
! Pass 2: incompressible pressure (order-1)
!---------------------------------------------------------------------------
IF(self%incomp)THEN
  CALL xdmf_plot_p%setup("xmhd_2d_p", "pressure/")
  CALL mesh%setup_io(xdmf_plot_p,oft_blagrange%order-1)
  rst_cur = rst_start
  DO
    IF(rst_cur > rst_end) EXIT
    WRITE(rst_char,100)rst_cur
    file_tmp='xmhd2d_'//rst_char//'.rst'
    rst_exist=oft_file_exist(TRIM(file_tmp))
    CALL oft_mpi_barrier(ierr)
    IF(.NOT.rst_exist) EXIT
    CALL hdf5_read(t,file_tmp,'t')
    CALL self%rst_load(u,file_tmp,'U')
    IF(self%linear)CALL u%add(1.d0, 1.d0, u0)
    CALL xdmf_plot_p%add_timestep(t)
    NULLIFY(plot_vals)
    CALL u%get_local(plot_vals,5)
    CALL mesh%save_vertex_scalar(plot_vals,xdmf_plot_p,'p')
    rst_cur = rst_cur + self%rst_freq
  END DO
END IF
CALL u%delete
CALL up%delete
CALL v%delete
CALL u0%delete
CALL ux%delete
CALL uy%delete
CALL uz%delete
DEALLOCATE(u,up,v,u0,ux,uy,uz)
IF(ASSOCIATED(plot_vals))DEALLOCATE(plot_vals)
IF(ASSOCIATED(plot_vec))DEALLOCATE(plot_vec)
IF(ASSOCIATED(plot_u0))DEALLOCATE(plot_u0)
IF(ASSOCIATED(psi_bg))DEALLOCATE(psi_bg)
end subroutine xmhd_2d_plot
!---------------------------------------------------------------------------
!> Load xMHD solution state from a restart file
!---------------------------------------------------------------------------
subroutine rst_load(self,u,filename,path,t,dt)
class(oft_xmhd_2d_sim), intent(inout) :: self
class(oft_vector), pointer, intent(inout) :: u !< Solution to load
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to store solution vector in file
real(r8), optional, intent(out) :: t !< Time of loaded solution
real(r8), optional, intent(out) :: dt !< Timestep at loaded time
DEBUG_STACK_PUSH
CALL self%fe_rep%vec_load(u,filename,path)
IF(PRESENT(t))CALL hdf5_read(t,filename,'t')
IF(PRESENT(dt))CALL hdf5_read(dt,filename,'dt')
DEBUG_STACK_POP
end subroutine rst_load

END MODULE xmhd_2d