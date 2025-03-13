!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_grad_shaf.F90
!
!> Grad-Shafranov implementation for TokaMaker
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_physics
!---------------------------------------------------------------------------
MODULE oft_gs
USE oft_base
USE oft_sort, ONLY: sort_matrix, sort_array
USE oft_io, ONLY: hdf5_field_exist, hdf5_read, hdf5_write, &
  xdmf_plot_file
  ! hdf5_rst_field_exist
USE oft_quadrature, ONLY: oft_quad_type
USE oft_gauss_quadrature, ONLY: set_quad_1d
USE oft_tet_quadrature, ONLY: set_quad_2d
USE oft_mesh_type, ONLY: oft_bmesh, bmesh_findcell, cell_is_curved
USE oft_mesh_local_util, ONLY: mesh_local_findedge
!---
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, oft_matrix, oft_graph, oft_graph_ptr
USE oft_la_utils, ONLY: create_matrix, graph_add_dense_blocks
USE oft_solver_base, ONLY: oft_solver, oft_eigsolver, oft_solver_bc
USE oft_solver_utils, ONLY: create_diag_pre, create_cg_solver
USE oft_native_solvers, ONLY: oft_native_cg_eigsolver
USE oft_lu, ONLY: oft_lusolver, lapack_matinv
!
USE fem_base, ONLY: oft_ml_fem_type, oft_afem_type
USE fem_utils, ONLY: bfem_interp, bfem_map_flag
USE oft_lag_basis, ONLY: oft_blag_eval, oft_blag_geval, &
  oft_blag_npos, oft_blag_d2eval, oft_scalar_bfem
USE oft_blag_operators, ONLY: oft_blag_project, oft_lag_brinterp, oft_lag_bginterp, &
  oft_lag_bg2interp, oft_blag_zerob, oft_blag_zerogrnd, oft_blag_getmop, &
  oft_blag_vproject
!---
USE fem_utils, ONLY: fem_interp
USE mhd_utils, ONLY: mu0
USE axi_green, ONLY: axi_coil_set, decay_eigenmodes, grad_green, green
USE tracing_2d, ONLY: active_tracer, tracinginv_fs, set_tracer, cylinv_interp
IMPLICIT NONE
#include "local.h"
INTEGER(4), PARAMETER :: max_xpoints = 20
!---------------------------------------------------------------------------
! CLASS oft_scalar_torus
!---------------------------------------------------------------------------
!> Interpolation class for an axisymmetric gaussian source in toroidal geometry
!!
!! In toroidal coordinates defined by the class the scalar field is defined
!! as
!! \f[ S = e^{-\frac{r^2}{\lambda}} \f]
!---------------------------------------------------------------------------
type, extends(bfem_interp) :: circular_curr
  real(r8) :: x0(2) = [1.d0,0.d0] !< Circle center
  real(r8) :: a = 0.d0 !< Minor radius
  real(r8) :: delta = 0.d0 !< Triangularity
  real(r8) :: kappa = 1.d0 !< Elongation
contains
  !> Reconstruct field
  procedure :: interp => circle_interp
end type circular_curr
!---------------------------------------------------------------------------
! CLASS flux_func
!---------------------------------------------------------------------------
!> Abstract flux function prototype
!---------------------------------------------------------------------------
TYPE, ABSTRACT :: flux_func
  INTEGER(i4) :: ncofs = 0 !< Number of free coefficients
  REAL(r8) :: f_offset = 0.d0
  REAL(r8) :: plasma_bounds(2) = (/-1.d99,1.d99/)
CONTAINS
  !>
  PROCEDURE(flux_func_eval), DEFERRED :: f
  !>
  PROCEDURE(flux_func_eval), DEFERRED :: fp
  PROCEDURE :: fpp => dummy_fpp
  !>
  PROCEDURE(flux_func_update), DEFERRED :: update
  !>
  PROCEDURE(flux_cofs_set), DEFERRED :: set_cofs
  !>
  PROCEDURE(flux_cofs_get), DEFERRED :: get_cofs
END TYPE flux_func
!---------------------------------------------------------------------------
! TYPE coil_region
!---------------------------------------------------------------------------
!> Internal coil region structure
!---------------------------------------------------------------------------
TYPE :: coil_region
  INTEGER(i4) :: nc = 0
  INTEGER(i4) :: id = 0
  INTEGER(i4), POINTER, DIMENSION(:) :: lc => NULL()
  REAL(r8) :: curr = 0.d0
  REAL(r8) :: vcont_gain = 0.d0
  REAL(r8) :: area = 0.d0
END TYPE coil_region
!---------------------------------------------------------------------------
! TYPE cond_region
!---------------------------------------------------------------------------
!> Internal wall region structure
!---------------------------------------------------------------------------
TYPE :: cond_region
  INTEGER(i4) :: nc = 0
  INTEGER(i4) :: nc_quad = 0
  INTEGER(i4) :: neigs = 0
  INTEGER(i4) :: id = 0
  INTEGER(i4) :: pair = -1
  ! LOGICAL :: limiter = .TRUE.
  LOGICAL :: continuous = .TRUE.
  REAL(r8) :: coverage = 1.d0
  REAL(r8) :: extent(2) = 0.d0
  REAL(r8) :: eta = -1.d0
  LOGICAL, POINTER, DIMENSION(:) :: fixed => NULL()
  INTEGER(i4), POINTER, DIMENSION(:) :: mtype => NULL()
  INTEGER(i4), POINTER, DIMENSION(:) :: mind => NULL()
  INTEGER(i4), POINTER, DIMENSION(:) :: eig_map => NULL()
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lc => NULL()
  REAL(r8), POINTER, DIMENSION(:) :: weights => NULL()
  REAL(r8), POINTER, DIMENSION(:) :: pair_signs => NULL()
  REAL(r8), POINTER, DIMENSION(:) :: fit_scales => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: cond_vals => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: rc => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: cond_curr => NULL()
  REAL(r8), POINTER, DIMENSION(:,:,:,:) :: corr_3d => NULL()
  CLASS(oft_vector_ptr), POINTER, DIMENSION(:) :: psi_eig => NULL()
END TYPE cond_region
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
TYPE :: gs_region_info
  INTEGER(i4) :: nnonaxi = 0
  INTEGER(i4) :: block_max = 0
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: reg_map => NULL()
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: dense_flag => NULL()
  INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:,:) :: node_mark => NULL()
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: nonaxi_vals => NULL()
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: noaxi_nodes => NULL()
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: bc_nodes => NULL()
END TYPE gs_region_info
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
type, extends(oft_solver_bc) :: oft_gs_zerob
  logical, pointer, dimension(:) :: node_flag => NULL()
  CLASS(oft_scalar_bfem), POINTER :: fe_rep => NULL() !< FE representation
contains
  procedure :: apply => zerob_apply
  procedure :: delete => zerob_delete
end type oft_gs_zerob
!---------------------------------------------------------------------------
! CLASS gs_eq
!---------------------------------------------------------------------------
!> Grad-Shafranov equilibrium object
!---------------------------------------------------------------------------
TYPE :: gs_eq
  INTEGER(i4) :: ierr = 0
  INTEGER(i4) :: maxits = 30
  INTEGER(i4) :: mode = 0
  INTEGER(i4) :: ninner = 4
  INTEGER(i4) :: nR0_ramp = 6
  INTEGER(i4) :: nlcfs = 100
  INTEGER(i4) :: nx_points = 0
  INTEGER(i4) :: ncoils = 0
  INTEGER(i4) :: ncoils_ext = 0
  INTEGER(i4) :: ncoil_regs = 0
  INTEGER(i4) :: nregularize = 0
  INTEGER(i4) :: nlimiter_pts = 0
  INTEGER(i4) :: nlimiter_nds = 0
  INTEGER(i4) :: ncond_regs = 0
  INTEGER(i4) :: ncond_eigs = 0
  INTEGER(i4) :: bc_nrhs = 0
  INTEGER(i4) :: isoflux_ntargets = 0
  INTEGER(i4) :: saddle_ntargets = 0
  INTEGER(i4) :: flux_ntargets = 0
  INTEGER(i4) :: nlim_con = 0
  INTEGER(i4) :: lim_nloops = 0
  REAL(r8) :: eps = 1.d-12
  REAL(r8) :: rmax = 0.d0
  REAL(r8) :: urf = .2d0
  REAL(r8) :: psimin = 0.d0
  REAL(r8) :: psiscale = 1.d0
  REAL(r8) :: psimax = 1.d0
  REAL(r8) :: alam = 1.d0
  REAL(r8) :: pnorm = 1.d0
  REAL(r8) :: rmin = 0.d0
  REAL(r8) :: dt = -1.d0
  REAL(r8) :: dt_last = -1.d0
  REAL(r8) :: Itor_target = -1.d0
  REAL(r8) :: estore_target = -1.d0
  REAL(r8) :: pax_target = -1.d0
  REAL(r8) :: Ip_ratio_target = -1.d99
  REAL(r8) :: R0_target = -1.d0
  REAL(r8) :: V0_target = -1.d99
  REAL(r8) :: nl_tol = 1.d-8
  REAL(r8) :: plasma_bounds(2) = (/-1.d99,1.d99/)
  REAL(r8) :: spatial_bounds(2,2) = RESHAPE((/-1.d99,1.d99,-1.d99,1.d99/),(/2,2/))
  REAL(r8) :: lim_zmax = 1.d99
  REAL(r8) :: lim_area = -1.d0
  REAL(r8) :: o_point(2) = (/-1.d0,1.d99/)
  REAL(r8) :: lim_point(2) = (/-1.d0,1.d99/)
  REAL(r8) :: x_points(2,max_xpoints) = 0.d0
  REAL(r8) :: x_vecs(2,max_xpoints) = 0.d0
  REAL(r8) :: vcontrol_val = 0.d0
  REAL(r8) :: vcont_I_gain = 1.d-8
  REAL(r8) :: vcont_D_gain = 0.d0
  REAL(r8) :: timing(4) = 0.d0
  REAL(r8) :: isoflux_grad_wt_lim = -1.d0
  LOGICAL, POINTER, DIMENSION(:) :: fe_flag => NULL()
  LOGICAL, POINTER, DIMENSION(:) :: saddle_pmask => NULL()
  LOGICAL, POINTER, DIMENSION(:) :: saddle_cmask => NULL()
  LOGICAL, POINTER, DIMENSION(:) :: saddle_rmask => NULL()
  INTEGER(i4), POINTER, DIMENSION(:) :: limiter_nds => NULL()
  INTEGER(i4), POINTER, DIMENSION(:) :: bc_rhs_list => NULL()
  INTEGER(i4), POINTER, DIMENSION(:) :: olbp => NULL() !< Oriented list of boundary points
  INTEGER(i4), POINTER, DIMENSION(:) :: lim_con => NULL()
  INTEGER(i4), POINTER, DIMENSION(:) :: lim_ptr => NULL()
  REAL(r8), POINTER, DIMENSION(:) :: cond_weights => NULL()
  REAL(r8), POINTER, DIMENSION(:) :: coil_reg_targets => NULL()
  REAL(r8), POINTER, DIMENSION(:) :: coil_currs => NULL()
  REAL(r8), POINTER, DIMENSION(:) :: coil_vcont => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: rlimiter_nds => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: limiter_pts => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: bc_mat => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: bc_coil_mat => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: bc_lmat => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: bc_bmat => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: rlcfs => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: isoflux_targets => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: saddle_targets => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: flux_targets => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: coil_reg_mat => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: coil_bounds => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: coil_nturns => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: Lcoils => NULL()
  LOGICAL :: free = .FALSE.
  LOGICAL :: compute_chi = .FALSE.
  LOGICAL :: plot_step = .TRUE.
  LOGICAL :: plot_final = .TRUE.
  LOGICAL :: boundary_limiter = .TRUE.
  LOGICAL :: diverted = .FALSE.
  LOGICAL :: has_plasma = .TRUE.
  LOGICAL :: full_domain = .FALSE.
#if defined( HAVE_SUPERLU ) || defined( HAVE_UMFPACK ) || defined( HAVE_MKL )
  LOGICAL :: use_lu = .TRUE.
#else
  LOGICAL :: use_lu = .FALSE.
#endif
  LOGICAL :: assym = .FALSE.
  ! LOGICAL :: active_x = .FALSE.
  LOGICAL :: limited_only = .FALSE.
  LOGICAL :: save_visit = .TRUE.
  LOGICAL :: isoflux_grad_weight = .TRUE.
  CHARACTER(LEN=OFT_PATH_SLEN) :: coil_file = 'none'
  CHARACTER(LEN=OFT_PATH_SLEN) :: limiter_file = 'none'
  TYPE(xdmf_plot_file) :: xdmf
  TYPE(oft_lusolver) :: lu_solver
  TYPE(oft_lusolver) :: lu_solver_dt
  ! CLASS(oft_solver), POINTER :: solver => NULL()
  ! CLASS(oft_solver), POINTER :: mop_solver => NULL()
  TYPE(axi_coil_set), POINTER, DIMENSION(:) :: coils_ext => NULL()
  TYPE(coil_region), POINTER, DIMENSION(:) :: coil_regions => NULL()
  TYPE(cond_region), POINTER, DIMENSION(:) :: cond_regions => NULL()
  TYPE(gs_region_info) :: region_info
  CLASS(oft_vector), POINTER :: psi => NULL() !<
  CLASS(oft_vector), POINTER :: chi => NULL() !<
  CLASS(oft_vector), POINTER :: u_hom => NULL() !<
  CLASS(oft_vector_ptr), POINTER, DIMENSION(:) :: psi_coil => NULL()
  ! CLASS(oft_vector_ptr), POINTER, DIMENSION(:) :: psi_cond => NULL()
  CLASS(oft_vector), POINTER :: psi_dt => NULL() !<
  CLASS(oft_matrix), POINTER :: dels => NULL()
  CLASS(oft_matrix), POINTER :: dels_dt => NULL()
  CLASS(oft_matrix), POINTER :: dels_full => NULL()
  CLASS(oft_matrix), POINTER :: mrop => NULL()
  CLASS(oft_matrix), POINTER :: mop => NULL()
  CLASS(flux_func), POINTER :: I => NULL() !< F*F' flux function
  CLASS(flux_func), POINTER :: P => NULL() !< Pressure flux function
  CLASS(flux_func), POINTER :: eta => NULL() !< Resistivity flux function
  CLASS(flux_func), POINTER :: I_NI => NULL() !< Non-inductive F*F' flux function
  CLASS(oft_bmesh), POINTER :: mesh => NULL() !< Mesh
  CLASS(oft_scalar_bfem), POINTER :: fe_rep => NULL()
  TYPE(oft_ml_fem_type), POINTER :: ML_fe_rep => NULL()
  TYPE(oft_blag_zerob), POINTER :: zerob_bc => NULL()
  TYPE(oft_blag_zerogrnd), POINTER :: zerogrnd_bc => NULL()
  TYPE(oft_gs_zerob), POINTER :: gs_zerob_bc => NULL()
  PROCEDURE(region_eta_set), NOPASS, POINTER :: set_eta => NULL()
CONTAINS
  !
  PROCEDURE :: setup => gs_setup
  !
  PROCEDURE :: init => gs_init
  !
  PROCEDURE :: init_psi => gs_init_psi
  !
  PROCEDURE :: load_coils => gs_load_coils
  !
  PROCEDURE :: load_limiters => gs_load_limiters
  !
  PROCEDURE :: solve => gs_solve
  !
  PROCEDURE :: lin_solve => gs_lin_solve
  !
  PROCEDURE :: vac_solve => gs_vac_solve
  !
  PROCEDURE :: get_chi => gs_get_chi
  !
  PROCEDURE :: itor => gs_itor
  !
  PROCEDURE :: beta => gs_beta
  !
  PROCEDURE :: delete => gs_destroy
END TYPE gs_eq
!---------------------------------------------------------------------------
!> Interpolate G-S profiles at a specific point in space
!!
!! @extends fem_base::bfem_interp
!---------------------------------------------------------------------------
type, extends(bfem_interp) :: gs_prof_interp
  INTEGER(i4) :: mode = 0
  class(gs_eq), pointer :: gs => NULL() !< Field for interpolation
  type(oft_lag_brinterp), pointer :: psi_eval => NULL()
  type(oft_lag_bginterp), pointer :: psi_geval => NULL()
contains
  procedure :: setup => gs_prof_interp_setup
  procedure :: delete => gs_prof_interp_delete
  !> Reconstruct a Lagrange scalar field
  procedure :: interp => gs_prof_interp_apply
end type gs_prof_interp
!---------------------------------------------------------------------------
!> Interpolate a Lagrange field.
!!
!! @extends fem_base::fem_interp
!---------------------------------------------------------------------------
type, extends(gs_prof_interp) :: gs_b_interp
contains
  !> Reconstruct a Lagrange scalar field
  procedure :: interp => gs_b_interp_apply
end type gs_b_interp
!---------------------------------------------------------------------------
! CLASS gsinv_interp
!---------------------------------------------------------------------------
!> Need docs
!---------------------------------------------------------------------------
type, extends(cylinv_interp) :: gsinv_interp
  LOGICAL :: compute_geom = .FALSE.
  real(8), pointer, dimension(:) :: uvals => NULL()
  class(oft_vector), pointer :: u => NULL() !< Field for interpolation
  class(oft_scalar_bfem), pointer :: lag_rep => NULL() !< Lagrange FE representation
contains
  procedure :: setup => gsinv_setup
  !> Reconstruct the gradient of a Lagrange scalar field
  procedure :: interp => gsinv_apply
  !>
  procedure :: delete => gsinv_destroy
end type gsinv_interp
!---
abstract interface
!---------------------------------------------------------------------------
! FUNCTION flux_func_eval
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
   function flux_func_eval(self,psi) result(b)
      import flux_func, r8
      class(flux_func), intent(inout) :: self
      real(r8), intent(in) :: psi
      real(r8) :: b
   end function flux_func_eval
!---------------------------------------------------------------------------
! SUBROUTINE flux_func_update
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
   subroutine flux_func_update(self,gseq)
      import flux_func, gs_eq
      class(flux_func), intent(inout) :: self
      class(gs_eq), intent(inout) :: gseq
   end subroutine flux_func_update
!---------------------------------------------------------------------------
! FUNCTION flux_cofs_set
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
   function flux_cofs_set(self,c) result(ierr)
      import flux_func, r8, i4
      class(flux_func), intent(inout) :: self
      real(r8), intent(in) :: c(:)
      integer(i4) :: ierr
   end function flux_cofs_set
!---------------------------------------------------------------------------
! SUBROUTINE flux_cofs_get
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
   subroutine flux_cofs_get(self,c)
      import flux_func, r8
      class(flux_func), intent(inout) :: self
      real(r8), intent(out) :: c(:)
   end subroutine flux_cofs_get
!---------------------------------------------------------------------------
! FUNCTION region_eta_set
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
   function region_eta_set(rc,id) result(eta)
      import r8, i4
      real(8), intent(in) :: rc(2)
      integer(4), intent(in) :: id
      real(8) :: eta
   end function region_eta_set
end interface
! logical, allocatable, dimension(:) :: node_flag
real(r8), PARAMETER :: gs_epsilon = 1.d-12
!
integer(i4) :: cell_active = 0
real(r8) :: pt_con_active(2) = [0.d0,0.d0]
real(r8) :: vec_con_active(2) = [0.d0,0.d0]
real(r8) :: psi_target_active = 0.d0
real(r8) :: qp_int_tol = 1.d-12
type(oft_lag_brinterp), pointer :: psi_eval_active => NULL()
type(oft_lag_bginterp), pointer :: psi_geval_active => NULL()
type(oft_lag_bg2interp), pointer :: psi_g2eval_active => NULL()
!$omp threadprivate(cell_active,pt_con_active,vec_con_active,psi_target_active)
!$omp threadprivate(psi_eval_active,psi_geval_active,psi_g2eval_active)
contains
!
function dummy_fpp(self,psi) result(b)
class(flux_func), intent(inout) :: self
real(r8), intent(in) :: psi
real(r8) :: b
b=0.d0
end function dummy_fpp
!---------------------------------------------------------------------------
!> Needs Docs
!!
!! @param[in,out] self G-S object
!---------------------------------------------------------------------------
subroutine gs_setup(self,ML_lag_2d)
class(gs_eq), intent(inout) :: self
class(oft_ml_fem_type), target, intent(inout) :: ML_lag_2d
SELECT TYPE(this=>ML_lag_2d%current_level)
  CLASS IS(oft_scalar_bfem)
    self%ML_fe_rep=>ML_lag_2d
    self%fe_rep=>this
  CLASS DEFAULT
    CALL oft_abort("Invalid FE space","gs_setup",__FILE__)
END SELECT
self%mesh=>self%fe_rep%mesh
ALLOCATE(self%zerob_bc)
self%zerob_bc%ML_lag_rep=>self%ML_fe_rep
ALLOCATE(self%zerogrnd_bc)
self%zerogrnd_bc%ML_lag_rep=>self%ML_fe_rep
end subroutine gs_setup
!---------------------------------------------------------------------------
!> Needs Docs
!!
!! @param[in,out] self G-S object
!---------------------------------------------------------------------------
subroutine gs_load_coils(self,ignore_inmesh)
class(gs_eq), intent(inout) :: self
logical, optional, intent(in) :: ignore_inmesh
!---XML solver fields
integer(i4) :: nread
TYPE(xml_node), POINTER :: doc,group_node,coil_set,coil,tmaker_group
TYPE(xml_nodelist) :: coil_sets,coils
!---
INTEGER(i4) :: i,j,ierr,cell
REAL(r8) :: f(3)
REAL(r8), PARAMETER :: tol=1.d-8
LOGICAL :: check_inmesh
check_inmesh=.TRUE.
IF(PRESENT(ignore_inmesh))check_inmesh=.NOT.ignore_inmesh
IF(TRIM(self%coil_file)=='none')RETURN
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Loading external coils:'
CALL oft_increase_indent
WRITE(*,'(3A)')oft_indent,'coil_file = ',TRIM(self%coil_file)
doc=>xml_parseFile(TRIM(self%coil_file),iostat=ierr)
CALL xml_get_element(doc,"tokamaker",tmaker_group,ierr)
CALL xml_get_element(tmaker_group,"coils",group_node,ierr)
!---Count coil sets
CALL xml_get_element(group_node,"coil_set",coil_sets,ierr)
self%ncoils_ext=coil_sets%n
ALLOCATE(self%coils_ext(self%ncoils_ext))
!---Setup coil sets
DO i=1,self%ncoils_ext
  coil_set=>coil_sets%nodes(i)%this
  !---
  CALL xml_extractDataAttribute(coil_set,"current",self%coils_ext(i)%curr,iostat=ierr)
  !---
  CALL xml_get_element(coil_set,"coil",coil_sets,ierr)
  self%coils_ext(i)%ncoils=coils%n
  ALLOCATE(self%coils_ext(i)%pt(2,self%coils_ext(i)%ncoils))
  ALLOCATE(self%coils_ext(i)%scale(self%coils_ext(i)%ncoils))
  self%coils_ext(i)%scale=1.d0
  DO j=1,self%coils_ext(i)%ncoils
    coil=>coils%nodes(j)%this
    CALL xml_extractDataContent(coil,self%coils_ext(i)%pt(:,j),num=nread,iostat=ierr)
    cell=0
    CALL bmesh_findcell(self%fe_rep%mesh,cell,self%coils_ext(i)%pt(:,j),f)
    IF((MAXVAL(f)<1.d0+tol).AND.(MINVAL(f)>-tol).AND.check_inmesh)THEN
      WRITE(*,*)'BAD COIL Found: ',i,self%coils_ext(i)%pt(:,j)
      CALL oft_abort('External coil in mesh','gs_load_coils',__FILE__)
    END IF
    !---Get polarity
    IF(xml_hasAttribute(coil,"scale"))CALL xml_extractDataAttribute(coil,"scale",self%coils_ext(i)%scale(j),num=nread,iostat=ierr)
  END DO
  IF(ASSOCIATED(coils%nodes))DEALLOCATE(coils%nodes)
END DO
IF(ASSOCIATED(coil_sets%nodes))DEALLOCATE(coil_sets%nodes)
!---
IF(oft_debug_print(2))THEN
  WRITE(*,*)
  WRITE(*,'(2A)')oft_indent,'Coils set definitions'
  WRITE(*,'(2A)')oft_indent,'========================='
  WRITE(*,'(2A,I4,A)')oft_indent,'Found ',self%ncoils_ext,' coil sets'
  CALL oft_increase_indent
  DO i=1,self%ncoils_ext
    WRITE(*,'(2A,ES11.3)')oft_indent,'Current [A] : ',self%coils_ext(i)%curr
    CALL oft_increase_indent
    DO j=1,self%coils_ext(i)%ncoils
      WRITE(*,'(2A,2ES11.3)')oft_indent,'Position  : ',self%coils_ext(i)%pt(:,j)
    END DO
    CALL oft_decrease_indent
    WRITE(*,*)
  END DO
  CALL oft_decrease_indent
  WRITE(*,'(2A)')oft_indent,'========================='
  WRITE(*,*)
ELSE
  WRITE(*,'(2A,I4,A)')oft_indent,'Found ',self%ncoils_ext,' coil sets'
END IF
!---Normalize currents
DO i=1,self%ncoils_ext
  self%coils_ext(i)%curr=self%coils_ext(i)%curr*mu0
END DO
CALL oft_decrease_indent
CALL gs_load_regions(self)
end subroutine gs_load_coils
!---------------------------------------------------------------------------
!> Needs Docs
!!
!! @param[in,out] self G-S object
!---------------------------------------------------------------------------
subroutine gs_load_regions(self)
class(gs_eq), intent(inout) :: self
!---XML solver fields
integer(4) :: nread
TYPE(xml_node), POINTER :: doc,region,field,tmaker_group
TYPE(xml_nodelist) :: regions,fields
!---
INTEGER(4) :: i,j,ierr,id,nregions,nreg_defs,nfields
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: region_flag,region_map
CHARACTER(LEN=10) :: reg_type
!---
IF(TRIM(self%coil_file)=='none')RETURN
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Loading internal coil and wall regions:'
CALL oft_increase_indent
WRITE(*,'(3A)')oft_indent,'coil_file = ',TRIM(self%coil_file)
doc=>xml_parseFile(TRIM(self%coil_file),iostat=ierr)
CALL xml_get_element(doc,"tokamaker",tmaker_group,ierr)
!---Count coil regions
CALL xml_get_element(tmaker_group,"region",regions,ierr)
nreg_defs=regions%n
nregions=MAXVAL(self%fe_rep%mesh%reg)
ALLOCATE(region_map(nreg_defs),region_flag(nreg_defs))
region_flag=0
region_map=0
!---
self%ncoil_regs=0
self%ncond_regs=0
DO i=1,nreg_defs
  region=>regions%nodes(i)%this
  CALL xml_extractDataAttribute(region,"id",id,num=nread,iostat=ierr)
  CALL xml_extractDataAttribute(region,"type",reg_type,iostat=ierr)
  IF(id<=0.OR.id>nregions)CALL oft_abort("Invalid region ID.","gs_load_regions",__FILE__)
  region_map(i)=id
  SELECT CASE(TRIM(reg_type))
    CASE("plasma")
      region_flag(i)=1
    CASE("wall")
      region_flag(i)=2
      self%ncond_regs=self%ncond_regs+1
    CASE("coil")
      region_flag(i)=3
      self%ncoil_regs=self%ncoil_regs+1
    CASE DEFAULT
      CALL oft_abort("Unknown region type.","gs_load_regions",__FILE__)
  END SELECT
END DO
!
ALLOCATE(self%cond_regions(self%ncond_regs))
ALLOCATE(self%coil_regions(self%ncoil_regs))
self%ncond_regs=0
self%ncoil_regs=0
!---Setup coil sets
DO i=1,nreg_defs
  region=>regions%nodes(i)%this
  SELECT CASE(region_flag(i))
    CASE(2)
      self%ncond_regs=self%ncond_regs+1
      self%cond_regions(self%ncond_regs)%id=region_map(i)
      !---
      CALL xml_get_element(region,"neigs",field,ierr)
      IF(ierr==0)THEN
        CALL xml_extractDataContent(field,self%cond_regions(self%ncond_regs)%neigs, &
             num=nread,iostat=ierr)
      END IF
      !---
      CALL xml_get_element(region,"eta",field,ierr)
      IF(ierr==0)THEN
        CALL xml_extractDataContent(field,self%cond_regions(self%ncond_regs)%eta, &
             num=nread,iostat=ierr)
      END IF
      !---
      IF(self%cond_regions(self%ncond_regs)%neigs>0)THEN
        ALLOCATE(self%cond_regions(self%ncond_regs)%fixed(self%cond_regions(self%ncond_regs)%neigs))
        self%cond_regions(self%ncond_regs)%fixed=.FALSE.
        CALL xml_get_element(region,"fixed",field,ierr)
        IF(ierr==0)THEN
          CALL xml_extractDataContent(field,self%cond_regions(self%ncond_regs)%fixed, &
               num=nread,iostat=ierr)
        END IF
        !
        ALLOCATE(self%cond_regions(self%ncond_regs)%weights(self%cond_regions(self%ncond_regs)%neigs))
        self%cond_regions(self%ncond_regs)%weights=1.d-5
        CALL xml_get_element(region,"weights",field,ierr)
        IF(ierr==0)THEN
          CALL xml_extractDataContent(field,self%cond_regions(self%ncond_regs)%weights, &
               num=nread,iostat=ierr)
        END IF
        !
        ALLOCATE(self%cond_regions(self%ncond_regs)%mtype(self%cond_regions(self%ncond_regs)%neigs))
        self%cond_regions(self%ncond_regs)%mtype=1
        CALL xml_get_element(region,"mtype",field,ierr)
        IF(ierr==0)THEN
          CALL xml_extractDataContent(field,self%cond_regions(self%ncond_regs)%mtype, &
               num=nread,iostat=ierr)
        END IF
        !
        ALLOCATE(self%cond_regions(self%ncond_regs)%mind(self%cond_regions(self%ncond_regs)%neigs))
        self%cond_regions(self%ncond_regs)%mind=(/(j,j=1,self%cond_regions(self%ncond_regs)%neigs)/)
        CALL xml_get_element(region,"mind",field,ierr)
        IF(ierr==0)THEN
          CALL xml_extractDataContent(field,self%cond_regions(self%ncond_regs)%mind, &
               num=nread,iostat=ierr)
        END IF
        !
        CALL xml_get_element(region,"pair",field,ierr)
        IF(ierr==0)THEN
          CALL xml_extractDataContent(field,self%cond_regions(self%ncond_regs)%pair, &
               num=nread,iostat=ierr)
        END IF
        !
        ALLOCATE(self%cond_regions(self%ncond_regs)%fit_scales(self%cond_regions(self%ncond_regs)%neigs))
        self%cond_regions(self%ncond_regs)%fit_scales = ABS(1.d0/self%cond_regions(self%ncond_regs)%weights)
        CALL xml_get_element(region,"fit_scales",field,ierr)
        IF(ierr==0)THEN
          CALL xml_extractDataContent(field,self%cond_regions(self%ncond_regs)%fit_scales, &
               num=nread,iostat=ierr)
        END IF
      END IF
      !---
      CALL xml_get_element(region,"continuous",field,ierr)
      IF(ierr==0)THEN
        CALL xml_extractDataContent(field,self%cond_regions(self%ncond_regs)%continuous, &
             num=nread,iostat=ierr)
        IF(.NOT.self%cond_regions(self%ncond_regs)%continuous)THEN
          CALL xml_get_element(region,"extent",field,ierr)
          IF(ierr==0)THEN
            CALL xml_extractDataContent(field,self%cond_regions(self%ncond_regs)%extent, &
                 num=nread,iostat=ierr)
          ELSE
            CALL oft_abort("No extents for non-continuous region","gs_load_regions",__FILE__)
          END IF
          !---Get toroidal coverage
          CALL xml_get_element(region,"coverage",field,ierr)
          IF(ierr==0)THEN
            CALL xml_extractDataContent(field,self%cond_regions(self%ncond_regs)%coverage, &
                 num=nread,iostat=ierr)
          END IF
        END IF
      END IF
      ! !---
      ! CALL xml_get_element(region,"limiter",field,ierr)
      ! IF(ierr==0)THEN
      !   CALL xml_extractDataContent(field,self%cond_regions(self%ncond_regs)%limiter, &
      !        num=nread,iostat=ierr)
      ! END IF
    CASE(3)
      self%ncoil_regs=self%ncoil_regs+1
      self%coil_regions(self%ncoil_regs)%id=region_map(i)
      !---
      CALL xml_get_element(region,"current",field,ierr)
      IF(ierr==0)THEN
        CALL xml_extractDataContent(field,self%coil_regions(self%ncoil_regs)%curr, &
             num=nread,iostat=ierr)
        self%coil_regions(self%ncoil_regs)%curr=self%coil_regions(self%ncoil_regs)%curr*mu0
      END IF
      !---
      CALL xml_get_element(region,"vcont_gain",field,ierr)
      IF(ierr==0)THEN
        CALL xml_extractDataContent(field,self%coil_regions(self%ncoil_regs)%vcont_gain, &
             num=nread,iostat=ierr)
      END IF
  END SELECT
END DO
IF(ASSOCIATED(regions%nodes))DEALLOCATE(regions%nodes)
!---
WRITE(*,'(2A,I4,A)')oft_indent,'Found ',self%ncond_regs,' conducting regions'
WRITE(*,'(2A,I4,A)')oft_indent,'Found ',self%ncoil_regs,' coil regions'
WRITE(*,*)
CALL oft_decrease_indent
end subroutine gs_load_regions
!---------------------------------------------------------------------------
! SUBROUTINE gs_load_limiters
!---------------------------------------------------------------------------
!> Needs Docs
!!
!! @param[in,out] self G-S object
!---------------------------------------------------------------------------
subroutine gs_load_limiters(self)
class(gs_eq), intent(inout) :: self
!---
INTEGER(4) :: i,io_unit,iostat
IF(TRIM(self%limiter_file)=='none')RETURN
IF(oft_debug_print(1))WRITE(*,'(2A,I4,A)')oft_indent,'Loading limiters'
OPEN(NEWUNIT=io_unit,FILE=TRIM(self%limiter_file))
READ(io_unit,*)self%nlimiter_pts
ALLOCATE(self%limiter_pts(2,self%nlimiter_pts))
DO i=1,self%nlimiter_pts
  READ(io_unit,*,IOSTAT=iostat)self%limiter_pts(:,i)
  IF(iostat<0)CALL oft_abort('EOF reached while reading limiter points', &
    'gs_load_limiters',__FILE__)
END DO
CLOSE(io_unit)
!---
IF(oft_debug_print(1))WRITE(*,'(2A,2X,I4,A)')oft_indent,'Found ', &
  self%nlimiter_pts,' limiter points'
end subroutine gs_load_limiters
!---------------------------------------------------------------------------
! SUBROUTINE gs_setup_walls
!---------------------------------------------------------------------------
!> Needs Docs
!!
!! @param[in,out] self G-S object
!---------------------------------------------------------------------------
subroutine gs_setup_walls(self,skip_load,make_plot)
class(gs_eq), intent(inout) :: self
logical, optional, intent(in) :: skip_load,make_plot
INTEGER(4) :: i,j,k,l,m,nphi_3d,cell
INTEGER(4), ALLOCATABLE :: eflag(:),j_lag(:),emark(:,:),pmark(:)
REAL(8) :: f(3),rcenter(2),vol,gop(3,3),pt(3)
REAL(8), ALLOCATABLE :: eigs(:)
LOGICAL :: file_exists,do_load,do_plot
CHARACTER(LEN=2) :: num_str,num_str2
CLASS(oft_bmesh), POINTER :: smesh
do_load=.TRUE.
do_plot=.TRUE.
IF(PRESENT(skip_load))do_load=skip_load
IF(PRESENT(make_plot))do_plot=make_plot
!---
IF(oft_debug_print(2))THEN
  WRITE(*,'(2A)')oft_indent,'Setup internal regions:'
  CALL oft_increase_indent
END IF
!---
smesh=>self%fe_rep%mesh
ALLOCATE(pmark(smesh%np))
pmark=0
DO i=1,smesh%np
  IF(smesh%bp(i))CYCLE
  pmark(i)=smesh%reg(smesh%lpc(smesh%kpc(i)))
  DO j=smesh%kpc(i),smesh%kpc(i+1)-1
    IF(pmark(i)/=smesh%reg(smesh%lpc(j)))THEN
      pmark(i)=-1
      EXIT
    END IF
  END DO
END DO
!---
DO i=1,self%ncoil_regs
  self%coil_regions(i)%nc=0
  DO j=1,smesh%nc
    IF(smesh%reg(j)==self%coil_regions(i)%id)self%coil_regions(i)%nc=self%coil_regions(i)%nc+1
  END DO
  IF(self%coil_regions(i)%nc==0)CYCLE
  ALLOCATE(self%coil_regions(i)%lc(self%coil_regions(i)%nc))
  self%coil_regions(i)%nc=0
  self%coil_regions(i)%lc=0
  self%coil_regions(i)%area=0.d0
  f=1.d0/3.d0
  DO j=1,smesh%nc
    IF(smesh%reg(j)==self%coil_regions(i)%id)THEN
      self%coil_regions(i)%nc=self%coil_regions(i)%nc+1
      self%coil_regions(i)%lc(self%coil_regions(i)%nc)=j
      CALL smesh%jacobian(j,f,gop,vol)
      self%coil_regions(i)%area=self%coil_regions(i)%area+vol
    END IF
  END DO
  self%coil_regions(i)%curr=self%coil_regions(i)%curr/self%coil_regions(i)%area
  self%coil_regions(i)%vcont_gain = &
    self%coil_regions(i)%vcont_gain/self%coil_regions(i)%area
END DO
!---
IF(do_load.AND.hdf5_field_exist('wall_eig.rst', 'ngrid_3d'))THEN
  CALL hdf5_read(vol,'wall_eig.rst','ngrid_3d')
  nphi_3d=INT(vol,4)
ELSE
  nphi_3d=0
END IF
self%ncond_eigs=0
DO i=1,self%ncond_regs
  !---
  self%cond_regions(i)%nc=0
  DO j=1,smesh%nc
    IF(smesh%reg(j)==self%cond_regions(i)%id)self%cond_regions(i)%nc=self%cond_regions(i)%nc+1
  END DO
  IF(self%cond_regions(i)%nc==0)CYCLE
  IF(MOD(self%cond_regions(i)%nc,4)/=0.OR.(self%cond_regions(i)%neigs==0))THEN
    self%cond_regions(i)%nc_quad=self%cond_regions(i)%nc
    ALLOCATE(self%cond_regions(i)%lc(1,self%cond_regions(i)%nc_quad))
    ALLOCATE(self%cond_regions(i)%cond_vals(self%cond_regions(i)%nc_quad,self%cond_regions(i)%neigs))
    ALLOCATE(self%cond_regions(i)%cond_curr(self%cond_regions(i)%nc_quad-1,self%cond_regions(i)%neigs))
    ALLOCATE(self%cond_regions(i)%rc(2,self%cond_regions(i)%nc_quad))
    self%cond_regions(i)%nc=0
    self%cond_regions(i)%lc=0
    self%cond_regions(i)%cond_curr=0.d0
    self%cond_regions(i)%cond_vals=0.d0
    self%cond_regions(i)%rc=0.d0
    DO j=1,smesh%nc
      IF(smesh%reg(j)==self%cond_regions(i)%id)THEN
        self%cond_regions(i)%nc=self%cond_regions(i)%nc+1
        self%cond_regions(i)%lc(1,self%cond_regions(i)%nc)=j
        self%cond_regions(i)%rc(:,self%cond_regions(i)%nc)= &
          (smesh%r(1:2,smesh%lc(1,j)) + smesh%r(1:2,smesh%lc(2,j)) &
          + smesh%r(1:2,smesh%lc(3,j)))/3.d0
      END IF
    END DO
    CYCLE
  END IF
  !CALL oft_abort('BAD COND region','gs_setup_walls',__FILE__)
  self%cond_regions(i)%nc_quad=self%cond_regions(i)%nc/4
  !---
  ALLOCATE(self%cond_regions(i)%lc(4,self%cond_regions(i)%nc_quad))
  ALLOCATE(self%cond_regions(i)%cond_vals(self%cond_regions(i)%nc_quad,self%cond_regions(i)%neigs))
  ALLOCATE(self%cond_regions(i)%cond_curr(self%cond_regions(i)%nc_quad-1,self%cond_regions(i)%neigs))
  ALLOCATE(self%cond_regions(i)%rc(2,self%cond_regions(i)%nc_quad))
  self%cond_regions(i)%nc=0
  self%cond_regions(i)%lc=0
  self%cond_regions(i)%cond_curr=0.d0
  self%cond_regions(i)%cond_vals=0.d0
  f=1.d0/3.d0
  k=0
  DO j=1,smesh%np
    IF(pmark(j)==self%cond_regions(i)%id)THEN
      k=k+1
      self%cond_regions(i)%rc(:,k) = smesh%r(1:2,j)
      ! WRITE(*,*)j
      IF(smesh%kpc(j+1)-smesh%kpc(j)>4)THEN
        WRITE(*,*)j,pmark(j),smesh%kpc(j+1)-smesh%kpc(j)
        DO l=smesh%kpc(j),smesh%kpc(j+1)-1
          WRITE(*,*)smesh%lpc(l),smesh%reg(smesh%lpc(l))
        END DO
        CALL oft_abort("Bad conductor","Bad",__FILE__)
      END IF
      DO l=smesh%kpc(j),smesh%kpc(j+1)-1
        self%cond_regions(i)%nc=self%cond_regions(i)%nc+1
        self%cond_regions(i)%lc(MOD(self%cond_regions(i)%nc-1,4)+1,CEILING(self%cond_regions(i)%nc/4.d0))=smesh%lpc(l)
      END DO
    END IF
  END DO
  ! DO j=1,smesh%nc
  !   IF(smesh%reg(j)==self%cond_regions(i)%id)THEN
  !     self%cond_regions(i)%nc=self%cond_regions(i)%nc+1
  !     self%cond_regions(i)%lc(MOD(self%cond_regions(i)%nc-1,4)+1,CEILING(self%cond_regions(i)%nc/4.d0))=j
  !     !---
  !     IF(MOD(self%cond_regions(i)%nc,4)==0)THEN
  !       k=k+1
  !       self%cond_regions(i)%rc(:,k) = smesh%r(1:2,smesh%lc(1,j))
  !     END IF
  !   END IF
  ! END DO
  IF(k/=self%cond_regions(i)%nc_quad)CALL oft_abort('BAD COUNTS!','gs_setup_walls',__FILE__)
  rcenter=SUM(self%cond_regions(i)%rc,DIM=2)/REAL(self%cond_regions(i)%nc_quad,8)
  IF(self%cond_regions(i)%neigs>0.AND.do_load)THEN
    WRITE(num_str,'(I2.2)')i
    INQUIRE(FILE='wall_eig.rst',EXIST=file_exists)
    IF(.NOT.file_exists)CALL oft_abort("No eddy current shape file","gs_setup_walls",__FILE__)
    CALL hdf5_read(self%cond_regions(i)%cond_vals, 'wall_eig.rst', 'eig_'//num_str)
    WRITE(*,*)MINVAL(self%cond_regions(i)%cond_vals)
    WRITE(*,*)MAXVAL(self%cond_regions(i)%cond_vals)
    self%cond_regions(i)%cond_vals=self%cond_regions(i)%cond_vals/4.d0
    !
    IF(nphi_3d>0)THEN
      ALLOCATE(self%cond_regions(i)%corr_3d(3,nphi_3d,smesh%np+smesh%ne,self%cond_regions(i)%neigs))
      self%cond_regions(i)%corr_3d=0.d0
      DO j=1,self%cond_regions(i)%neigs
        WRITE(num_str2,'(I2.2)')j
        IF(hdf5_field_exist('wall_eig.rst', 'rz_corr_r'//num_str//'_'//num_str2))THEN
          CALL hdf5_read(self%cond_regions(i)%corr_3d(1,:,:,j), 'wall_eig.rst', 'rz_corr_r'//num_str//'_'//num_str2)
          CALL hdf5_read(self%cond_regions(i)%corr_3d(2,:,:,j), 'wall_eig.rst', 'rz_corr_t'//num_str//'_'//num_str2)
          CALL hdf5_read(self%cond_regions(i)%corr_3d(3,:,:,j), 'wall_eig.rst', 'rz_corr_z'//num_str//'_'//num_str2)
        END IF
      END DO
    END IF
  END IF
  !---
  DO j=1,self%cond_regions(i)%neigs
    !---Compute poloidal flows
    self%cond_regions(i)%cond_curr(1,j)=self%cond_regions(i)%cond_vals(1,j)
    DO k=2,self%cond_regions(i)%nc_quad-1
      self%cond_regions(i)%cond_curr(k,j)=self%cond_regions(i)%cond_curr(k-1,j) &
      + self%cond_regions(i)%cond_vals(k,j)
    END DO
  END DO
  IF(self%cond_regions(i)%pair<0)THEN
    ALLOCATE(self%cond_regions(i)%eig_map(self%cond_regions(i)%neigs))
    self%cond_regions(i)%eig_map=(/(j,j=1,self%cond_regions(i)%neigs)/)+self%ncond_eigs
    self%ncond_eigs=self%ncond_eigs+self%cond_regions(i)%neigs
  ELSE
    ALLOCATE(self%cond_regions(i)%pair_signs(self%cond_regions(i)%neigs))
    self%cond_regions(i)%pair_signs = self%cond_regions(i)%weights
  END IF
  !---
  IF(oft_debug_print(2))THEN
    WRITE(*,'(2A,I8,A)')oft_indent,'Found ',self%cond_regions(i)%nc,' conductor cells'
    WRITE(*,'(2A,2ES11.3)')oft_indent,'Rcenter = ',rcenter
    WRITE(*,*)
  END IF
END DO
DO i=1,self%ncond_regs
  IF(self%cond_regions(i)%pair>0)THEN
    self%cond_regions(i)%eig_map=>self%cond_regions(self%cond_regions(i)%pair)%eig_map
  END IF
END DO
ALLOCATE(self%cond_weights(self%ncond_eigs))
CALL gs_get_cond_weights(self,self%cond_weights,.FALSE.)
CALL gs_set_cond_weights(self,self%cond_weights,.FALSE.)
!---Setup limiters
ALLOCATE(eflag(self%fe_rep%ne),j_lag(self%fe_rep%nce))
eflag=0
! DO j=1,self%ncond_regs
!   IF(self%cond_regions(j)%limiter)THEN
!     DO k=1,self%cond_regions(j)%nc_quad
!       DO m=1,4
!         i=self%cond_regions(j)%lc(m,k)
!         CALL self%fe_rep%ncdofs(i,j_lag)
!         DO l=1,self%fe_rep%nce
!           eflag(j_lag(l))=1
!         END DO
!       END DO
!     END DO
!   END IF
! END DO
!---Add nodes on boundary of region 1
ALLOCATE(emark(2,self%fe_rep%ne))
emark=0
DO i=1,smesh%nc
  IF(MAXVAL(ABS(smesh%r(2,smesh%lc(:,i))))>self%lim_zmax)CYCLE
  CALL self%fe_rep%ncdofs(i,j_lag)
  IF(smesh%reg(i)==1)THEN
    emark(1,j_lag)=1
  ELSE
    emark(2,j_lag)=1
  END IF
END DO
DO i=1,self%fe_rep%ne
  IF(emark(1,i)==1.AND.emark(2,i)==1)eflag(i)=1
  IF(emark(1,i)==1.AND.self%fe_rep%be(i))eflag(i)=1
END DO
!---
self%nlimiter_nds=SUM(eflag)
ALLOCATE(self%limiter_nds(self%nlimiter_nds),self%rlimiter_nds(2,self%nlimiter_nds))
self%nlimiter_nds=0
DO i=1,self%fe_rep%ne
  IF(eflag(i)==1)THEN
    self%nlimiter_nds=self%nlimiter_nds+1
    self%limiter_nds(self%nlimiter_nds)=i
  END IF
END DO
!
DO i=1,self%nlimiter_nds
  ! Get position
  cell=self%fe_rep%lec(self%fe_rep%kec(self%limiter_nds(i)))
  CALL self%fe_rep%ncdofs(cell,j_lag)
  DO l=1,self%fe_rep%nce
    IF(j_lag(l)==self%limiter_nds(i))EXIT
  END DO
  IF(l>self%fe_rep%nce)CALL oft_abort("BAD", "gs_setup_walls", __FILE__)
  CALL oft_blag_npos(self%fe_rep,cell,l,f)
  pt=smesh%log2phys(cell,f)
  self%rlimiter_nds(:,i)=pt(1:2)
END DO
DEALLOCATE(eflag,j_lag)
IF(oft_debug_print(2))THEN!.AND.((self%ncoil_regs>0).OR.(self%ncond_regs>0)))THEN
  WRITE(*,'(2A,I8,A)')oft_indent,'Found ',self%nlimiter_nds,' material limiter nodes'
  WRITE(*,*)
END IF
!---Mark non-continuous regions
CALL set_noncontinuous
!---
IF(do_plot)THEN
  ALLOCATE(eigs(smesh%nc))
  ! DO j=1,smesh%nc
  !   eigs(j)=REAL(smesh%reg(j),8)
  ! END DO
  ! CALL smesh%save_cell_scalar(eigs,self%xdmf,'Regions')
  !---
  eigs=0.d0
  DO j=1,self%ncond_regs
    IF(self%cond_regions(j)%neigs>0)THEN
      DO k=1,self%cond_regions(j)%nc_quad
        DO l=1,4
          i=self%cond_regions(j)%lc(l,k)
          eigs(i)=eigs(i) &
          + DOT_PRODUCT(self%cond_regions(j)%weights,self%cond_regions(j)%cond_vals(k,:))
        END DO
      END DO
    END IF
  END DO
  ! CALL smesh%save_cell_scalar(eigs,self%xdmf,'Eig')
  !---
  eigs=0.d0
  DO j=1,self%ncond_regs
    IF(self%cond_regions(j)%neigs>0)THEN
      DO k=1,self%cond_regions(j)%nc_quad-1
        DO l=1,4
          i=self%cond_regions(j)%lc(l,k)
          eigs(i)=eigs(i) + self%cond_regions(j)%cond_curr(k,1)
        END DO
      END DO
    END IF
  END DO
  ! CALL smesh%save_cell_scalar(eigs,self%xdmf,'Curr')
  DEALLOCATE(eigs)
END IF
IF(oft_debug_print(2))CALL oft_decrease_indent
contains
!
subroutine set_noncontinuous
integer(4) :: i,m,jr,jc
integer(4), allocatable, dimension(:) :: j
self%region_info%nnonaxi=0
DO jr=1,self%ncond_regs
    IF(self%cond_regions(jr)%continuous)CYCLE
    self%region_info%nnonaxi=self%region_info%nnonaxi+1
END DO
ALLOCATE(self%region_info%reg_map(smesh%nreg))
self%region_info%reg_map=0
IF(self%region_info%nnonaxi>0)THEN
  ALLOCATE(self%region_info%node_mark(smesh%nreg,self%fe_rep%ne+1))
  self%region_info%node_mark=0
  allocate(j(self%fe_rep%nce))
  DO i=1,self%fe_rep%mesh%nc
      call self%fe_rep%ncdofs(i,j)
      DO jc=1,self%fe_rep%nce
          IF(self%region_info%node_mark(smesh%reg(i),j(jc))/=0)CYCLE
          self%region_info%node_mark(smesh%reg(i),self%fe_rep%ne+1)=self%region_info%node_mark(smesh%reg(i),self%fe_rep%ne+1)+1
          self%region_info%node_mark(smesh%reg(i),j(jc))=self%region_info%node_mark(smesh%reg(i),self%fe_rep%ne+1)
      END DO
  END DO
  deallocate(j)
  ALLOCATE(self%region_info%noaxi_nodes(self%region_info%nnonaxi))
  self%region_info%block_max=0
  m=0
  DO jr=1,self%ncond_regs
      IF(self%cond_regions(jr)%continuous)CYCLE
      i=self%cond_regions(jr)%id
      m=m+1
      self%region_info%reg_map(i)=m
      self%region_info%noaxi_nodes(m)%n=self%region_info%node_mark(i,self%fe_rep%ne+1)
      ALLOCATE(self%region_info%noaxi_nodes(m)%v(self%region_info%noaxi_nodes(m)%n))
      self%region_info%noaxi_nodes(m)%v=0
      DO jc=1,self%fe_rep%ne
          IF(self%region_info%node_mark(i,jc)>0)self%region_info%noaxi_nodes(m)%v(self%region_info%node_mark(i,jc)) = jc
      END DO
      self%region_info%block_max=MAX(self%region_info%block_max,self%region_info%noaxi_nodes(m)%n)
  END DO
END IF
end subroutine set_noncontinuous
end subroutine gs_setup_walls
! !---------------------------------------------------------------------------
! ! SUBROUTINE gs_setup_cflag
! !---------------------------------------------------------------------------
! !> Need docs
! !---------------------------------------------------------------------------
! subroutine gs_setup_cflag(self)
! class(gs_eq), intent(inout) :: self
! integer(4) :: i,j,k,l
! real(8) :: pt(2)
! logical :: in_bounds
! !---Set flag
! ALLOCATE(self%cflag(smesh%nc))
! self%cflag=3
! DO i=1,smesh%nc
!   DO j=1,3
!     pt=smesh%r(1:2,smesh%lc(j,i))
!     in_bounds=.TRUE.
!     in_bounds=in_bounds.AND.(pt(1)>=self%spatial_bounds(1,1).AND.pt(1)<=self%spatial_bounds(2,1))
!     in_bounds=in_bounds.AND.(pt(2)>=self%spatial_bounds(1,2).AND.pt(2)<=self%spatial_bounds(2,2))
!     IF(in_bounds)self%cflag(i)=0
!   END DO
! END DO
! DO j=1,self%ncond_regs
!   DO k=1,self%cond_regions(j)%nc_quad
!     DO l=1,4
!       i=self%cond_regions(j)%lc(l,k)
!       self%cflag(i)=1
!     END DO
!   END DO
! END DO
! DO j=1,self%ncoil_regs
!   DO k=1,self%coil_regions(j)%nc
!     i=self%coil_regions(j)%lc(k)
!     self%cflag(i)=2
!   END DO
! END DO
! i=COUNT(self%cflag==0)
! j=COUNT(self%cflag==3)
! k=smesh%nc-i-j
! IF(oft_debug_print(1))THEN
!   WRITE(*,'(A,2X,A)')oft_indent,'Marking cells:'
!   WRITE(*,'(A,4X,A,I8)')oft_indent,'Plasma = ',i
!   WRITE(*,'(A,4X,A,I8)')oft_indent,'Coil   = ',k
!   WRITE(*,'(A,4X,A,I8)')oft_indent,'Cond   = ',j
!   WRITE(*,'(A,4X,A,I8)')oft_indent,'Total  = ',smesh%nc
! END IF
! end subroutine gs_setup_cflag
!---------------------------------------------------------------------------
!> Initialize Grad-Shafranov solution with the Taylor state
!!
!! @param[in,out] self G-S object
!---------------------------------------------------------------------------
subroutine gs_init(self)
class(gs_eq), intent(inout) :: self
! logical, optional, intent(in) :: compute
! real(8), optional, intent(in) :: r0(2),a,kappa,delta
type(oft_native_cg_eigsolver) :: eigsolver
class(oft_vector), pointer :: tmp_vec,tmp_vec2
integer(4), pointer, dimension(:) :: cdofs
real(r8), pointer, dimension(:) :: psi_vals
type(oft_lag_brinterp) :: psi_eval 
integer(4) :: i,j,k,mind,nCon,ierr
integer(4), allocatable :: cells(:)
real(r8) :: itor,curr,f(3),goptmp(3,4),pol_val(1),v,pt(2),theta
real(r8), allocatable :: err_mat(:,:),rhs(:),err_inv(:,:),currs(:)
character(LEN=3) :: coil_tag
character(LEN=2) :: cond_tag,eig_tag
CLASS(oft_bmesh), POINTER :: smesh
! logical :: do_compute
! do_compute=.TRUE.
! IF(PRESENT(compute))do_compute=compute
smesh=>self%fe_rep%mesh
! self%ML_fe_rep=>ML_oft_blagrange
! self%fe_rep=>oft_blagrange
! ALLOCATE(self%zerob_bc)
! self%zerob_bc%ML_lag_rep=>self%ML_fe_rep
! ALLOCATE(self%zerogrnd_bc)
! self%zerogrnd_bc%ML_lag_rep=>self%ML_fe_rep
!---Get Vector
call self%fe_rep%vec_create(self%psi)
call self%psi%set(0.d0)
call self%fe_rep%vec_create(self%chi)
call self%chi%set(0.d0)
IF(self%free)THEN
  call self%fe_rep%vec_create(self%u_hom)
  call self%u_hom%set(0.d0)
END IF
self%rmax=-1.d0 !MAXVAL(smesh%r(1,:))
DO i=1,smesh%nc
  IF(smesh%reg(i)/=1)CYCLE
  DO j=1,3
    self%rmax=MAX(self%rmax,smesh%r(1,smesh%lc(j,i)))
  END DO
END DO
!---Build operators
IF(.NOT.ASSOCIATED(self%dels))THEN
  IF(self%free)THEN
    CALL compute_bcmat(self)
    CALL build_dels(self%dels,self,"free")
  ELSE
    CALL build_dels(self%dels,self,"zerob")
  END IF
END IF
IF(.NOT.ASSOCIATED(self%mrop))CALL build_mrop(self%fe_rep,self%mrop,"none")
IF(.NOT.ASSOCIATED(self%mop))CALL oft_blag_getmop(self%fe_rep,self%mop,"none")
!---Setup boundary conditions
ALLOCATE(self%gs_zerob_bc)
self%gs_zerob_bc%fe_rep=>self%fe_rep
ALLOCATE(self%gs_zerob_bc%node_flag(self%fe_rep%ne),cdofs(self%fe_rep%nce))
IF(.NOT.ASSOCIATED(self%saddle_rmask))THEN
  ALLOCATE(self%saddle_rmask(smesh%nreg))
  self%saddle_rmask=.TRUE.
  self%saddle_rmask(1)=.FALSE.
END IF
ALLOCATE(self%saddle_cmask(smesh%nc),self%saddle_pmask(smesh%np))
self%spatial_bounds(:,1)=[1.d99,-1.d99]
self%spatial_bounds(:,2)=[1.d99,-1.d99]
self%gs_zerob_bc%node_flag=self%fe_rep%be
self%saddle_pmask=.TRUE.
self%saddle_cmask=.TRUE.
DO i=1,smesh%nc
  IF(smesh%reg(i)>1)THEN
    CALL self%fe_rep%ncdofs(i,cdofs)
    DO j=1,self%fe_rep%nce
      self%gs_zerob_bc%node_flag(cdofs(j))=.TRUE.
    END DO
    IF(.NOT.self%saddle_rmask(smesh%reg(i)))THEN
      self%saddle_cmask(i)=.FALSE.
      DO j=1,3
        self%saddle_pmask(smesh%lc(j,i))=.FALSE.
      END DO
    END IF
  ELSE
    self%saddle_cmask(i)=.FALSE.
    DO j=1,3
      self%saddle_pmask(smesh%lc(j,i))=.FALSE.
      self%spatial_bounds(1,1)=MIN(self%spatial_bounds(1,1),smesh%r(1,smesh%lc(j,i)))
      self%spatial_bounds(2,1)=MAX(self%spatial_bounds(2,1),smesh%r(1,smesh%lc(j,i)))
      self%spatial_bounds(1,2)=MIN(self%spatial_bounds(1,2),smesh%r(2,smesh%lc(j,i)))
      self%spatial_bounds(2,2)=MAX(self%spatial_bounds(2,2),smesh%r(2,smesh%lc(j,i)))
    END DO
  END IF
END DO
DEALLOCATE(cdofs)
!
DO i=1,smesh%np
  IF(self%saddle_pmask(i))CYCLE
  DO j=smesh%kpc(i),smesh%kpc(i+1)-1
    self%saddle_cmask(smesh%lpc(j))=.FALSE.
  END DO
END DO
self%saddle_pmask=.FALSE.
DO i=1,smesh%nc
  IF(self%saddle_cmask(i))self%saddle_pmask(smesh%lc(:,i))=.TRUE.
END DO
self%saddle_pmask=self%saddle_pmask.OR.smesh%bp
CALL get_limiter
self%lim_area=0.d0
DO i=1,smesh%nc
  IF(smesh%reg(i)==1)self%lim_area=self%lim_area+smesh%ca(i)
END DO
!
NULLIFY(tmp_vec,psi_vals)
call self%psi%new(tmp_vec)
!---Compute coil fields
IF(self%ncoils==0)THEN
  self%ncoils=self%ncoil_regs+self%ncoils_ext
  ALLOCATE(self%coil_nturns(smesh%nreg+self%ncoils_ext,self%ncoils))
  ALLOCATE(self%coil_currs(self%ncoils),self%coil_vcont(self%ncoils))
  self%coil_nturns=0.d0
  self%coil_vcont=0.d0
  self%coil_currs=0.d0
  DO i=1,self%ncoil_regs
    self%coil_nturns(self%coil_regions(i)%id,i)=1.d0
    self%coil_currs(i)=self%coil_regions(i)%curr*self%coil_regions(i)%area
    self%coil_vcont(self%coil_regions(i)%id)=self%coil_regions(i)%vcont_gain
  END DO
  DO i=1,self%ncoils_ext
    self%coil_nturns(smesh%nreg+i,self%ncoil_regs+i)=1.d0
    self%coil_currs(self%ncoil_regs+i)=self%coils_ext(i)%curr
  END DO
ELSE
  DO i=1,self%ncoil_regs
    self%coil_nturns(self%coil_regions(i)%id,:) = &
      self%coil_nturns(self%coil_regions(i)%id,:)/self%coil_regions(i)%area ! Normalize turns by coil area
  END DO
END IF
ALLOCATE(self%psi_coil(self%ncoils),self%Lcoils(self%ncoils+1,self%ncoils+1))
self%Lcoils=0.d0
DO i=1,self%ncoils
  CALL self%psi%new(self%psi_coil(i)%f)
  CALL gs_coil_source(self,i,tmp_vec)
  CALL self%zerob_bc%apply(tmp_vec)
  CALL gs_vacuum_solve(self,self%psi_coil(i)%f,tmp_vec)
  CALL self%psi_coil(i)%f%get_local(psi_vals)
  WRITE(coil_tag,'(I3.3)')i
  IF(self%save_visit)CALL smesh%save_vertex_scalar(psi_vals,self%xdmf,'Psi_coil'//coil_tag)
END DO
DO i=1,self%ncoils
  DO j=i,self%ncoils
    CALL gs_coil_mutual(self,i,self%psi_coil(j)%f,self%Lcoils(i,j))
    self%Lcoils(i,j)=self%Lcoils(i,j)
    IF(j>i)self%Lcoils(j,i)=self%Lcoils(i,j)
  END DO
END DO
! ALLOCATE(self%psi_cond(self%ncond_eigs))
call self%psi%new(tmp_vec2)
DO i=1,self%ncond_regs
  ALLOCATE(self%cond_regions(i)%psi_eig(self%cond_regions(i)%neigs))
  DO j=1,self%cond_regions(i)%neigs
    CALL self%psi%new(self%cond_regions(i)%psi_eig(j)%f)
    CALL gs_cond_source(self,i,j,tmp_vec)
    CALL self%zerob_bc%apply(tmp_vec)
    CALL gs_vacuum_solve(self,self%cond_regions(i)%psi_eig(j)%f,tmp_vec)
    CALL self%cond_regions(i)%psi_eig(j)%f%get_local(psi_vals)
    WRITE(cond_tag,'(I2.2)')i
    WRITE(eig_tag,'(I2.2)')j
    IF(self%save_visit)CALL smesh%save_vertex_scalar(psi_vals,self%xdmf,'Psi_cond'//cond_tag//"_"//eig_tag)
  END DO
END DO
CALL tmp_vec2%delete()
DEALLOCATE(tmp_vec2)
! IF(self%isoflux_ntargets>0)CALL gs_fit_isoflux(self,ierr)
! IF(do_compute)THEN
!   !---Add coil/conductor fields to IC
!   DO i=1,self%ncoil_regs
!     curr = self%coil_regions(i)%curr & 
!       + self%coil_regions(i)%vcont_gain*self%vcontrol_val
!     CALL self%psi%add(1.d0,curr,self%psi_coil(i)%f)
!   END DO
!   DO i=1,self%ncond_regs
!     DO j=1,self%cond_regions(i)%neigs
!       ! k=self%cond_regions(i)%eig_map(j)
!       CALL self%psi%add(1.d0,self%cond_regions(i)%weights(j), &
!         self%cond_regions(i)%psi_eig(j)%f)
!     END DO
!   END DO
! END IF
! !
! CALL self%psi%get_local(psi_vals)
! IF(self%save_visit)CALL smesh%save_vertex_scalar(psi_vals,self%xdmf,'Psi_init')
CALL tmp_vec%delete()
DEALLOCATE(tmp_vec)
contains
subroutine get_limiter()
integer(4) :: i,ii,istart,iloop,j,k,l,orient(2),nmax
integer(4), allocatable :: tmp_ptr(:)
real(8) :: val
logical, allocatable :: eflag(:)
allocate(eflag(smesh%ne))
eflag=.FALSE.
DO i=1,smesh%nc
  IF(smesh%reg(i)/=1)CYCLE
  DO j=1,3
    IF(smesh%lcc(j,i)==0)THEN
      eflag(ABS(smesh%lce(j,i)))=.TRUE.
    ELSE
      IF(smesh%reg(smesh%lcc(j,i))/=1)eflag(ABS(smesh%lce(j,i)))=.TRUE.
    END IF
  END DO
END DO
self%nlim_con=COUNT(eflag)
ALLOCATE(self%lim_con(self%nlim_con),tmp_ptr(self%nlim_con))
tmp_ptr=0
self%lim_con=0
istart=1
self%lim_nloops=0
nmax=0
DO WHILE(.TRUE.)
  DO i=1,smesh%ne
    IF(eflag(i))EXIT
  END DO
  IF(i>smesh%ne)EXIT
  self%lim_nloops=self%lim_nloops+1
  tmp_ptr(self%lim_nloops)=istart
  self%lim_con(istart)=smesh%le(1,i)
  eflag(i)=.FALSE. ! Unmark edge after use
  !
  DO ii=istart+1,self%nlim_con
    DO j=smesh%kpe(self%lim_con(ii-1)),smesh%kpe(self%lim_con(ii-1)+1)-1
      i=smesh%lpe(j)
      IF(.NOT.eflag(i))CYCLE
      k=SUM(smesh%le(:,i))-self%lim_con(ii-1)
      IF(ii>2)THEN
        IF(k==self%lim_con(ii-2))CYCLE
      END IF
      self%lim_con(ii)=k
      eflag(i)=.FALSE. ! Unmark edge after use
      EXIT
    END DO
    IF(j>=smesh%kpe(self%lim_con(ii-1)+1))EXIT
  END DO
  istart=ii
END DO
IF(MINVAL(self%lim_con)==0)CALL oft_abort('Zero value detected in lim_con','gs_init::get_limiter',__FILE__)
ALLOCATE(self%lim_ptr(self%lim_nloops+1))
self%lim_ptr(self%lim_nloops+1)=self%nlim_con+1
self%lim_ptr(1:self%lim_nloops)=tmp_ptr(1:self%lim_nloops)
deallocate(eflag,tmp_ptr)
! !---Check orientation
! orient=0.d0
! DO i=1,smesh%nbp
!     IF(i<smesh%nbp/2)THEN
!       orient(1)=orient(1)+smesh%r(2,olbp(i))
!     ELSE
!       orient(2)=orient(2)+smesh%r(2,olbp(i))
!     END IF
! END DO
! !---Reverse list if necessary
! IF(orient(2)>orient(1))THEN
!   DO i=2,smesh%nbp/2+1
!     k=olbp(i)
!     olbp(i)=olbp(smesh%nbp+2-i)
!     olbp(smesh%nbp+2-i)=k
!   END DO
! END IF
end subroutine get_limiter
end subroutine gs_init
!---------------------------------------------------------------------------
! SUBROUTINE gs_init_psi
!---------------------------------------------------------------------------
!> Initialize Grad-Shafranov solution with the Taylor state
!!
!! @param[in,out] self G-S object
!---------------------------------------------------------------------------
subroutine gs_init_psi(self,ierr,r0,a,kappa,delta,curr_source)
class(gs_eq), intent(inout) :: self
integer(4), intent(out) :: ierr
real(8), optional, intent(in) :: r0(2),a,kappa,delta,curr_source(:)
type(oft_native_cg_eigsolver) :: eigsolver
class(oft_vector), pointer :: tmp_vec,tmp_vec2
integer(4), pointer, dimension(:) :: cdofs
real(r8), pointer, dimension(:) :: psi_vals
type(circular_curr) :: circle_init
type(oft_lag_brinterp) :: psi_eval 
integer(4) :: i,j,k,mind,nCon
integer(4), allocatable :: cells(:)
real(r8) :: itor,curr,f(3),goptmp(3,4),pol_val(1),v,pt(2),theta
real(r8), allocatable :: err_mat(:,:),rhs(:),err_inv(:,:),currs(:)
character(LEN=3) :: coil_tag
character(LEN=2) :: cond_tag,eig_tag
!
ierr=0
NULLIFY(tmp_vec,psi_vals)
call self%psi%new(tmp_vec)
IF(PRESENT(r0))THEN
  IF(r0(1)>0.d0)THEN
    circle_init%x0=r0
    circle_init%a=a
    circle_init%kappa=kappa
    circle_init%delta=delta
    circle_init%mesh=>self%fe_rep%mesh
    CALL gs_gen_source(self,circle_init,tmp_vec)
  ELSE
    CALL self%psi%set(1.d0)
    CALL self%mop%apply(self%psi,tmp_vec)
  END IF
  CALL self%psi%set(0.d0)
  CALL gs_vacuum_solve(self,self%psi,tmp_vec)
  itor=self%itor()
  IF(self%Itor_target>0.d0)CALL self%psi%scale(self%Itor_target/itor)
ELSE IF(PRESENT(curr_source))THEN
  CALL self%psi%set(0.d0)
  CALL tmp_vec%restore_local(curr_source)
  CALL gs_vacuum_solve(self,self%psi,tmp_vec)
  itor=self%itor()
  IF(self%Itor_target>0.d0)CALL self%psi%scale(self%Itor_target/itor)
ELSE
  !---Setup Solver
  eigsolver%A=>self%dels
  eigsolver%M=>self%mrop
  eigsolver%bc=>self%gs_zerob_bc
  eigsolver%its=-2
  !---Setup Preconditioner
  ! if(oft_blagrange_level==1)then
    CALL create_diag_pre(eigsolver%pre)
  ! else
  !   call lag_mlprecon_dels(eigsolver%pre,level=oft_blagrange_level)
  ! end if
  !---Compute eigenstate
  CALL self%psi%get_local(psi_vals)
  CALL oft_random_number(psi_vals,self%psi%n)
  CALL self%psi%restore_local(psi_vals)
  WRITE(*,'(2A)')oft_indent,'Computing force-free eigenmode'
  CALL oft_increase_indent
  call eigsolver%apply(self%psi,self%alam)
  self%alam=sqrt(self%alam)
  WRITE(*,'(2A,ES11.3)')oft_indent,'Lambda = ',self%alam
  CALL oft_decrease_indent
  !---
  IF(.NOT.self%free)THEN
    CALL self%psi%get_local(psi_vals)
    mind=maxloc(abs(psi_vals),DIM=1)
    call self%psi%scale(1.d0/psi_vals(mind))
    self%psimax=1.d0
  END IF
  IF(self%Itor_target>0.d0)THEN
    itor=self%itor()
    CALL self%psi%scale(self%Itor_target/itor)
  END IF
  !---Cleanup
  CALL eigsolver%pre%delete
  CALL eigsolver%delete
  DEALLOCATE(psi_vals)
END IF
CALL tmp_vec%delete()
DEALLOCATE(tmp_vec)
!
CALL self%psi%get_local(psi_vals)
self%plasma_bounds=[0.d0,MAXVAL(ABS(psi_vals))]
DEALLOCATE(psi_vals)
! CALL tmp_vec2%delete()
! DEALLOCATE(tmp_vec2)
IF(self%isoflux_ntargets+self%flux_ntargets+self%saddle_ntargets>0)THEN
  CALL gs_fit_isoflux(self,self%psi,ierr)
  IF(ierr/=0)THEN
    ierr=-7
    RETURN
  END IF
END IF
!---Add coil/conductor fields to IC
DO i=1,self%ncoils
  ! curr = self%coil_regions(i)%curr & 
  !   + self%coil_regions(i)%vcont_gain*self%vcontrol_val
  CALL self%psi%add(1.d0,self%coil_currs(i),self%psi_coil(i)%f)
END DO
! CALL self%psi%add(1.d0,self%vcontrol_val,self%psi_coil(self%ncoils+1)%f)
DO i=1,self%ncond_regs
  DO j=1,self%cond_regions(i)%neigs
    ! k=self%cond_regions(i)%eig_map(j)
    CALL self%psi%add(1.d0,self%cond_regions(i)%weights(j), &
      self%cond_regions(i)%psi_eig(j)%f)
  END DO
END DO
!
self%plasma_bounds=[-1.d99,1.d99]
CALL gs_update_bounds(self)
CALL self%I%update(self)
CALL self%p%update(self)
!
IF(self%save_visit)THEN
  CALL self%psi%get_local(psi_vals)
  CALL self%fe_rep%mesh%save_vertex_scalar(psi_vals,self%xdmf,'Psi_init')
  DEALLOCATE(psi_vals)
END IF
! CALL tmp_vec%delete()
! DEALLOCATE(tmp_vec)
end subroutine gs_init_psi
!---------------------------------------------------------------------------
!> Zero a surface Lagrange scalar field at all edge nodes
!!
!! @param[in,out] a Field to be zeroed
!---------------------------------------------------------------------------
subroutine zerob_apply(self,a)
class(oft_gs_zerob), intent(inout) :: self
class(oft_vector), intent(inout) :: a
integer(i4) :: i,j
real(r8), pointer, dimension(:) :: vloc
DEBUG_STACK_PUSH
NULLIFY(vloc)
!---Cast to vector type
NULLIFY(vloc)
call a%get_local(vloc)
!---Zero boundary values
!$omp parallel do
do i=1,self%fe_rep%ne
  IF(self%node_flag(i))vloc(i)=0.d0
end do
!---
call a%restore_local(vloc)
deallocate(vloc)
DEBUG_STACK_POP
end subroutine zerob_apply
!---------------------------------------------------------------------------
!> Zero a surface Lagrange scalar field at all edge nodes
!!
!! @param[in,out] a Field to be zeroed
!---------------------------------------------------------------------------
subroutine zerob_delete(self)
class(oft_gs_zerob), intent(inout) :: self
NULLIFY(self%fe_rep)
IF(ASSOCIATED(self%node_flag))DEALLOCATE(self%node_flag)
end subroutine zerob_delete
!---------------------------------------------------------------------------
!> Evaluate torus source
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [1]
!---------------------------------------------------------------------------
subroutine circle_interp(self,cell,f,gop,val)
class(circular_curr), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,3)
real(r8), intent(out) :: val(:)
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: r(2),lam,coord_tmp(2),pt(3)
!---MINPACK variables
real(8) :: ftol,xtol,gtol,epsfcn,factor,cofs(1),error(1),gpsitmp(2)
real(8), allocatable, dimension(:) :: diag,wa1,wa2,wa3,wa4,qtf
real(8), allocatable, dimension(:,:) :: fjac
integer(4) :: maxfev,mode,nprint,info,nfev,ldfjac,ncons,ncofs
integer(4), allocatable, dimension(:) :: ipvt
!---Get coordinates
pt=self%mesh%log2phys(cell,f)
r = pt(1:2) - self%x0

ncons=2
ncofs=2
allocate(diag(ncofs),fjac(ncons,ncofs))
allocate(qtf(ncofs),wa1(ncofs),wa2(ncofs))
allocate(wa3(ncofs),wa4(ncons))
allocate(ipvt(ncofs))
mode = 1
factor = 1.d0
maxfev = 100
ftol = 1.d-9
xtol = 1.d-8
gtol = 1.d-8
epsfcn = 1.d-4
nprint = 0
ldfjac = ncons
pt_con_active=r
coord_tmp=[SQRT(SUM(r**2)),ATAN2(r(2),r(1))]
call lmdif(circ_error,ncons,ncofs,coord_tmp,gpsitmp, &
            ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
            nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
deallocate(diag,fjac,qtf,wa1,wa2)
deallocate(wa3,wa4,ipvt)
!
lam=self%a*10.d0
val=(1.d0-TANH((coord_tmp(1)-self%a)/lam))/2.d0
CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE psimax_error
!------------------------------------------------------------------------------
!>
!------------------------------------------------------------------------------
SUBROUTINE circ_error(m,n,cofs,err,iflag)
integer(4), intent(in) :: m,n
real(8), intent(in) :: cofs(n)
real(8), intent(out) :: err(m)
integer(4), intent(inout) :: iflag
real(8) :: pt_eval(2)
pt_eval = cofs(1)*[COS(cofs(2)+ASIN(self%delta)*SIN(cofs(2))), &
  self%kappa*SIN(cofs(2))]
err = pt_eval-pt_con_active
end subroutine circ_error
end subroutine circle_interp
!---------------------------------------------------------------------------
! SUBROUTINE gs_vacuum_solve
!---------------------------------------------------------------------------
!> Compute Grad-Shafranov solution for current flux definitions
!!
!! @param[in,out] self G-S object
!! @param[out] ierr Error flag
!---------------------------------------------------------------------------
subroutine gs_vacuum_solve(self,pol_flux,source,ierr)
class(gs_eq), intent(inout) :: self
class(oft_vector), intent(inout) :: pol_flux,source
integer(4), optional, intent(out) :: ierr
class(oft_vector), pointer :: rhs,u_hom,psi_bc
logical :: pm_save
IF(self%use_lu)THEN
  IF(.NOT.ASSOCIATED(self%lu_solver%A))THEN
    self%lu_solver%A=>self%dels
    ALLOCATE(self%lu_solver%sec_rhs(self%psi%n,2))
  END IF
ELSE
  CALL oft_abort("LU solver required for GS solve","gs_solve",__FILE__)
END IF
!---Create worker vectors
CALL pol_flux%new(rhs)
!---Solve
CALL rhs%add(0.d0,1.d0,source)
CALL pol_flux%set(0.d0)
CALL self%zerob_bc%apply(rhs)
!---Solve linear system
pm_save=oft_env%pm; oft_env%pm=.FALSE.
CALL self%lu_solver%apply(pol_flux,rhs)
oft_env%pm=pm_save
!---Cleanup
CALL rhs%delete()
DEALLOCATE(rhs)
END SUBROUTINE gs_vacuum_solve
!---------------------------------------------------------------------------
! SUBROUTINE gs_gen_source
!---------------------------------------------------------------------------
!> Needs Docs
!!
!! @param[in,out] self G-S object
!! @param[in,out] a Psi field
!! @param[in,out] b Source field
!---------------------------------------------------------------------------
subroutine gs_gen_source(self,source_fun,b)
class(gs_eq), intent(inout) :: self
CLASS(bfem_interp), intent(inout) :: source_fun
CLASS(oft_vector), intent(inout) :: b
real(r8), pointer, dimension(:) :: btmp
real(8) :: psitmp,goptmp(3,3),det,pt(3),v,ffp(3),t1,source_tmp(1)
real(8), allocatable :: rhs_loc(:),cond_fac(:),rop(:),vcache(:)
integer(4) :: j,m,l,k
integer(4), allocatable :: j_lag(:)
logical :: curved
t1=omp_get_wtime()
!---
NULLIFY(btmp)
! IF(.NOT.ASSOCIATED(self%cflag))CALL gs_setup_cflag(self)
call b%set(0.d0)
CALL b%get_local(btmp)
!---
!!$omp parallel private(j,rhs_loc,j_lag,ffp,curved,goptmp,v,m,det,pt,psitmp,l,rop,source_tmp)
allocate(rhs_loc(self%fe_rep%nce))
allocate(rop(self%fe_rep%nce))
allocate(j_lag(self%fe_rep%nce))
!!$omp do schedule(static,1)
DO j=1,self%fe_rep%mesh%nc
  IF(self%fe_rep%mesh%reg(j)/=1)CYCLE
  call self%fe_rep%ncdofs(j,j_lag)
  rhs_loc=0.d0
  curved=cell_is_curved(self%fe_rep%mesh,j)
  do m=1,self%fe_rep%quad%np
    if(curved.OR.(m==1))call self%fe_rep%mesh%jacobian(j,self%fe_rep%quad%pts(:,m),goptmp,v)
    det=v*self%fe_rep%quad%wts(m)
    DO l=1,self%fe_rep%nce
      CALL oft_blag_eval(self%fe_rep,j,l,self%fe_rep%quad%pts(:,m),rop(l))
    END DO
    CALL source_fun%interp(j,self%fe_rep%quad%pts(:,m),goptmp,source_tmp)
    !$omp simd
    do l=1,self%fe_rep%nce
      rhs_loc(l)=rhs_loc(l)+rop(l)*source_tmp(1)*det
    end do
  end do
  !---Get local to global DOF mapping
  do l=1,self%fe_rep%nce
    m = j_lag(l)
    !!$omp atomic
    btmp(m)=btmp(m)+rhs_loc(l)
  end do
end do
deallocate(rhs_loc,j_lag,rop)
!!$omp end parallel
! DEALLOCATE(cond_fac)
CALL b%restore_local(btmp,add=.TRUE.)
DEALLOCATE(btmp)
! self%timing(2)=self%timing(2)+(omp_get_wtime()-t1)
end subroutine gs_gen_source
!---------------------------------------------------------------------------
! SUBROUTINE gs_coil_source
!---------------------------------------------------------------------------
!> Calculates coil current source
!---------------------------------------------------------------------------
subroutine gs_coil_source(self,iCoil,b)
class(gs_eq), intent(inout) :: self !< G-S Object
integer(4), intent(in) :: iCoil !< Coil index
CLASS(oft_vector), intent(inout) :: b !< Coil current source
real(r8), pointer, dimension(:) :: btmp
real(8) :: psitmp,goptmp(3,3),det,pt(3),v,ffp(3),t1,nturns
real(8), allocatable :: rhs_loc(:),cond_fac(:),rop(:),vcache(:)
integer(4) :: j,m,l,k
integer(4), allocatable :: j_lag(:)
logical :: curved
t1=omp_get_wtime()
!---
NULLIFY(btmp)
! IF(.NOT.ASSOCIATED(self%cflag))CALL gs_setup_cflag(self)
call b%set(0.d0)
CALL b%get_local(btmp)
!---
!$omp parallel private(rhs_loc,j_lag,ffp,curved,goptmp,v,m,det,pt,psitmp,l,rop,nturns)
allocate(rhs_loc(self%fe_rep%nce))
allocate(rop(self%fe_rep%nce))
allocate(j_lag(self%fe_rep%nce))
!$omp do schedule(static,1)
DO j=1,self%fe_rep%mesh%nc
  nturns=self%coil_nturns(self%fe_rep%mesh%reg(j),iCoil)
  IF(ABS(nturns)<1.d-10)CYCLE
  call self%fe_rep%ncdofs(j,j_lag)
  rhs_loc=0.d0
  curved=cell_is_curved(self%fe_rep%mesh,j)
  do m=1,self%fe_rep%quad%np
    if(curved.OR.(m==1))call self%fe_rep%mesh%jacobian(j,self%fe_rep%quad%pts(:,m),goptmp,v)
    det=v*self%fe_rep%quad%wts(m)
    DO l=1,self%fe_rep%nce
      CALL oft_blag_eval(self%fe_rep,j,l,self%fe_rep%quad%pts(:,m),rop(l))
    END DO
    !$omp simd
    do l=1,self%fe_rep%nce
      rhs_loc(l)=rhs_loc(l)+rop(l)*det
    end do
  end do
  !---Get local to global DOF mapping
  do l=1,self%fe_rep%nce
    m = j_lag(l)
    !$omp atomic
    btmp(m)=btmp(m)+rhs_loc(l)*nturns
  end do
end do
deallocate(rhs_loc,j_lag,rop)
!$omp end parallel
! DEALLOCATE(cond_fac)
CALL b%restore_local(btmp,add=.TRUE.)
DEALLOCATE(btmp)
! self%timing(2)=self%timing(2)+(omp_get_wtime()-t1)
end subroutine gs_coil_source
!---------------------------------------------------------------------------
!> Calculates field contribution due to coil with non-uniform current distribution
!---------------------------------------------------------------------------
subroutine gs_coil_source_distributed(self,iCoil,b,curr_dist)
class(gs_eq), intent(inout) :: self !< G-S object
integer(4), intent(in) :: iCoil !< Coil index
CLASS(oft_vector), intent(inout) :: b !< Coil current source
REAL(8), POINTER, DIMENSION(:), intent(in) :: curr_dist !< Relative current density
real(r8), pointer, dimension(:) :: btmp
real(8) :: psitmp,goptmp(3,3),det,pt(3),v,t1,nturns
real(8), allocatable :: rhs_loc(:),rop(:)
integer(4) :: j,m,l,k
integer(4), allocatable :: j_lag(:)
class(oft_bmesh), pointer :: smesh
! t1=omp_get_wtime()
smesh=>self%fe_rep%mesh
!---
NULLIFY(btmp)
CALL b%set(0.d0)
CALL b%get_local(btmp)
!$omp parallel private(j,rhs_loc,j_lag,goptmp,v,m,det,pt,psitmp,l,rop,nturns)
allocate(rhs_loc(self%fe_rep%nce))
allocate(rop(self%fe_rep%nce))
allocate(j_lag(self%fe_rep%nce))
!$omp do schedule(static,1)
DO j=1,smesh%nc
  nturns=self%coil_nturns(smesh%reg(j),iCoil)
  IF(ABS(nturns)<1.d-10)CYCLE
  call self%fe_rep%ncdofs(j,j_lag)
  rhs_loc=0.d0
  do m=1,self%fe_rep%quad%np
    call smesh%jacobian(j,self%fe_rep%quad%pts(:,m),goptmp,v)
    det=v*self%fe_rep%quad%wts(m)
    DO l=1,self%fe_rep%nce
      CALL oft_blag_eval(self%fe_rep,j,l,self%fe_rep%quad%pts(:,m),rop(l)) 
    END DO
    !$omp simd
    do l=1,self%fe_rep%nce
      rhs_loc(l)=rhs_loc(l)+rop(l)*det*curr_dist(j_lag(l))
    end do
  end do
  !---Get local to global DOF mapping
  do l=1,self%fe_rep%nce
    m = j_lag(l)
    !$omp atomic
    btmp(m)=btmp(m)+rhs_loc(l)*nturns
  end do
END DO
deallocate(rhs_loc,j_lag,rop)
!$omp end parallel
CALL b%restore_local(btmp,add=.TRUE.)
DEALLOCATE(btmp)
end subroutine gs_coil_source_distributed
!---------------------------------------------------------------------------
! SUBROUTINE gs_cond_source
!---------------------------------------------------------------------------
!> Needs Docs
!!
!! @param[in,out] self G-S object
!! @param[in,out] a Psi field
!! @param[in,out] b Source field
!---------------------------------------------------------------------------
subroutine gs_cond_source(self,iCond,iMode,b)
class(gs_eq), intent(inout) :: self
integer(4), intent(in) :: iCond,iMode
CLASS(oft_vector), intent(inout) :: b
real(r8), pointer, dimension(:) :: btmp
real(8) :: psitmp,goptmp(3,3),det,pt(3),v,ffp(3),t1
real(8), allocatable :: rhs_loc(:),cond_fac(:),rop(:),vcache(:)
integer(4) :: j,m,l,k,kk
integer(4), allocatable :: j_lag(:)
logical :: curved
t1=omp_get_wtime()
!---
NULLIFY(btmp)
! IF(.NOT.ASSOCIATED(self%cflag))CALL gs_setup_cflag(self)
call b%set(0.d0)
CALL b%get_local(btmp)
!---
!$omp parallel private(j,kk,rhs_loc,j_lag,ffp,curved,goptmp,v,m,det,pt,psitmp,l,rop,vcache)
allocate(rhs_loc(self%fe_rep%nce))
allocate(rop(self%fe_rep%nce))
allocate(j_lag(self%fe_rep%nce))
!$omp do schedule(static,1)
DO k=1,self%cond_regions(iCond)%nc_quad
DO kk=1,4
  j=self%cond_regions(iCond)%lc(kk,k)
  psitmp=self%cond_regions(iCond)%cond_vals(k,iMode)/self%fe_rep%mesh%ca(j)
  call self%fe_rep%ncdofs(j,j_lag)
  rhs_loc=0.d0
  curved=cell_is_curved(self%fe_rep%mesh,j)
  do m=1,self%fe_rep%quad%np
    if(curved.OR.(m==1))call self%fe_rep%mesh%jacobian(j,self%fe_rep%quad%pts(:,m),goptmp,v)
    det=v*self%fe_rep%quad%wts(m)
    DO l=1,self%fe_rep%nce
      CALL oft_blag_eval(self%fe_rep,j,l,self%fe_rep%quad%pts(:,m),rop(l))
    END DO
    !$omp simd
    do l=1,self%fe_rep%nce
      rhs_loc(l)=rhs_loc(l)+rop(l)*psitmp*det
    end do
  end do
  !---Get local to global DOF mapping
  do l=1,self%fe_rep%nce
    m = j_lag(l)
    !$omp atomic
    btmp(m)=btmp(m)+rhs_loc(l)
  end do
END DO
end do
deallocate(rhs_loc,j_lag,rop)
!$omp end parallel
! DEALLOCATE(cond_fac)
CALL b%restore_local(btmp,add=.TRUE.)
DEALLOCATE(btmp)
! self%timing(2)=self%timing(2)+(omp_get_wtime()-t1)
end subroutine gs_cond_source
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gs_wall_source(self,dpsi_dt,b)
class(gs_eq), intent(inout) :: self !< G-S object
CLASS(oft_vector), intent(inout) :: dpsi_dt,b
real(r8), pointer, dimension(:) :: btmp
real(8) :: psitmp,goptmp(3,3),det,pt(3),v,ffp(3),t1
real(8), allocatable :: rhs_loc(:),cond_fac(:),rop(:),vcache(:),eta_reg(:),reg_source(:)
real(8), pointer :: psi_vals(:)
integer(4) :: j,m,l,k,kk,iCond
integer(4), allocatable :: j_lag(:)
logical :: curved
! t1=omp_get_wtime()
!---
NULLIFY(btmp,psi_vals)
call b%set(0.d0)
CALL b%get_local(btmp)
CALL dpsi_dt%get_local(psi_vals)
!
ALLOCATE(eta_reg(self%fe_rep%mesh%nreg),reg_source(self%fe_rep%mesh%nreg))
reg_source=0.d0
eta_reg=-1.d0
DO l=1,self%ncond_regs
  j=self%cond_regions(l)%id
  eta_reg(j)=self%cond_regions(l)%eta
END DO
IF(ASSOCIATED(self%region_info%nonaxi_vals))THEN
  DO l=1,self%fe_rep%mesh%nreg
    reg_source(l)=DOT_PRODUCT(psi_vals,self%region_info%nonaxi_vals(:,l))
  END DO
END IF
!---
!$omp parallel private(rhs_loc,j_lag,curved,goptmp,v,m,det,pt,psitmp,l,rop)
allocate(rhs_loc(self%fe_rep%nce))
allocate(rop(self%fe_rep%nce))
allocate(j_lag(self%fe_rep%nce))
!$omp do schedule(static,1)
DO j=1,self%fe_rep%mesh%nc
  IF(eta_reg(self%fe_rep%mesh%reg(j))<0.d0)CYCLE
  call self%fe_rep%ncdofs(j,j_lag)
  rhs_loc=0.d0
  curved=cell_is_curved(self%fe_rep%mesh,j)
  do m=1,self%fe_rep%quad%np
    if(curved.OR.(m==1))call self%fe_rep%mesh%jacobian(j,self%fe_rep%quad%pts(:,m),goptmp,v)
    det=v*self%fe_rep%quad%wts(m)
    pt=self%fe_rep%mesh%log2phys(j,self%fe_rep%quad%pts(:,m))
    psitmp=0.d0
    DO l=1,self%fe_rep%nce
      CALL oft_blag_eval(self%fe_rep,j,l,self%fe_rep%quad%pts(:,m),rop(l))
      psitmp=psitmp+psi_vals(j_lag(l))*rop(l)
    END DO
    psitmp=psitmp/eta_reg(self%fe_rep%mesh%reg(j))/(pt(1)+gs_epsilon)
    !$omp simd
    do l=1,self%fe_rep%nce
      rhs_loc(l)=rhs_loc(l)+rop(l)*(psitmp+reg_source(self%fe_rep%mesh%reg(j)))*det
    end do
  end do
  !---Get local to global DOF mapping
  do l=1,self%fe_rep%nce
    m = j_lag(l)
    !$omp atomic
    btmp(m)=btmp(m)+rhs_loc(l)
  end do
END DO
deallocate(rhs_loc,j_lag,rop)
!$omp end parallel
CALL b%restore_local(btmp,add=.TRUE.)
DEALLOCATE(btmp,psi_vals,eta_reg,reg_source)
! self%timing(2)=self%timing(2)+(omp_get_wtime()-t1)
end subroutine gs_wall_source
!---------------------------------------------------------------------------
!> Compute inductance between coil and given poloidal flux
!---------------------------------------------------------------------------
subroutine gs_coil_mutual(self,iCoil,b,mutual)
class(gs_eq), intent(inout) :: self !< G-S object
integer(4), intent(in) :: iCoil !< Coil index
CLASS(oft_vector), intent(inout) :: b !< \f$ \psi \f$ for mutual calculation

real(8), intent(out) :: mutual !< Mutual inductance \f$ \int I_C \psi dV / I_C \f$
real(r8), pointer, dimension(:) :: btmp
real(8) :: goptmp(3,3),det,v,t1,psi_tmp,nturns
real(8), allocatable :: rhs_loc(:),cond_fac(:),rop(:)
integer(4) :: j,m,l,k
integer(4), allocatable :: j_lag(:)
logical :: curved
!---
NULLIFY(btmp)
CALL b%get_local(btmp)
!---
mutual=0.d0
!$omp parallel private(j,j_lag,curved,goptmp,v,m,det,psi_tmp,l,rop,nturns) reduction(+:mutual)
allocate(rop(self%fe_rep%nce))
allocate(j_lag(self%fe_rep%nce))
!$omp do schedule(static,1)
DO j=1,self%fe_rep%mesh%nc
  nturns=self%coil_nturns(self%fe_rep%mesh%reg(j),iCoil)
  IF(ABS(nturns)<1.d-10)CYCLE
  call self%fe_rep%ncdofs(j,j_lag)
  curved=cell_is_curved(self%fe_rep%mesh,j)
  do m=1,self%fe_rep%quad%np
    if(curved.OR.(m==1))call self%fe_rep%mesh%jacobian(j,self%fe_rep%quad%pts(:,m),goptmp,v)
    det=v*self%fe_rep%quad%wts(m)
    psi_tmp=0.d0
    DO l=1,self%fe_rep%nce
      CALL oft_blag_eval(self%fe_rep,j,l,self%fe_rep%quad%pts(:,m),rop(l))
      psi_tmp=psi_tmp+btmp(j_lag(l))*rop(l)
    END DO
    mutual = mutual + psi_tmp*det*nturns
  end do
end do
deallocate(j_lag,rop)
!$omp end parallel
mutual=mu0*2.d0*pi*mutual
DEALLOCATE(btmp)
end subroutine gs_coil_mutual
!---------------------------------------------------------------------------
!> Compute inductance between a coil with non-uniform current distribution and given poloidal flux
!---------------------------------------------------------------------------
subroutine gs_coil_mutual_distributed(self, iCoil, b, curr_dist, mutual)
class(gs_eq), intent(inout) :: self !< G-S object
integer(4), intent(in) :: iCoil !< Coil index
CLASS(oft_vector), intent(inout) :: b !< \f$ \psi \f$ for mutual calculation
REAL(8), POINTER, DIMENSION(:), intent(in) :: curr_dist !< Relative current density
real(8), intent(out) :: mutual !< Mutual inductance \f$ \int I_C \psi dV / I_C \f$
real(r8), pointer, dimension(:) :: btmp
real(8) :: goptmp(3,3),det,v,t1,psi_tmp,nturns,j_phi
real(8), allocatable :: rhs_loc(:),cond_fac(:),rop(:)
integer(4) :: j,m,l,k
integer(4), allocatable :: j_lag(:)
logical :: curved
class(oft_bmesh), pointer :: smesh
smesh=>self%fe_rep%mesh
!---
NULLIFY(btmp)
CALL b%get_local(btmp)
!---
mutual=0.d0
!$omp parallel private(j,j_lag,curved,goptmp,v,m,det,l,rop,nturns,j_phi,psi_tmp) reduction(+:mutual)
allocate(rop(self%fe_rep%nce))
allocate(j_lag(self%fe_rep%nce))
!$omp do schedule(static,1)
DO j=1,smesh%nc
  nturns=self%coil_nturns(smesh%reg(j),iCoil)
  IF(ABS(nturns)<1.d-10)CYCLE
  call self%fe_rep%ncdofs(j,j_lag)
  do m=1,self%fe_rep%quad%np
    call smesh%jacobian(j,self%fe_rep%quad%pts(:,m),goptmp,v)
    det=v*self%fe_rep%quad%wts(m)
    psi_tmp=0.d0
    j_phi=0.d0
    DO l=1,self%fe_rep%nce
      CALL oft_blag_eval(self%fe_rep,j,l,self%fe_rep%quad%pts(:,m),rop(l))
      psi_tmp=psi_tmp+btmp(j_lag(l))*rop(l)
      j_phi=j_phi+curr_dist(j_lag(l))*rop(l)
    END DO
    mutual = mutual + psi_tmp*det*nturns*j_phi
  end do
end do
deallocate(j_lag,rop)
!$omp end parallel
mutual=mu0*2.d0*pi*mutual
DEALLOCATE(btmp)
end subroutine gs_coil_mutual_distributed
!---------------------------------------------------------------------------
!> Compute inductance between plasma current and given poloidal flux
!---------------------------------------------------------------------------
subroutine gs_plasma_mutual(self,b,mutual,itor)
class(gs_eq), intent(inout) :: self !< G-S solver object
CLASS(oft_vector), intent(inout) :: b !< \f$ \psi \f$ for mutual calculation
real(8), intent(out) :: mutual !< Mutual inductance \f$ \int J_p \psi dV / I_p \f$
real(8), intent(out) :: itor !< Plasma toroidal current
real(r8), pointer, dimension(:) :: btmp
real(8) :: psitmp(1),goptmp(3,3),det,pt(3),v,t1,b_tmp,itor_loc
real(8), allocatable :: rhs_loc(:),cond_fac(:),rop(:)
integer(4) :: j,m,l,k
integer(4), allocatable :: j_lag(:)
logical :: curved
type(oft_lag_brinterp), target :: psi_eval
! t1=omp_get_wtime()
!---
NULLIFY(btmp)
CALL b%get_local(btmp)
psi_eval%u=>self%psi
CALL psi_eval%setup(self%fe_rep)
!---
mutual=0.d0
itor=0.d0
!$omp parallel private(j,j_lag,curved,goptmp,v,m,det,pt,psitmp,b_tmp,l,rop,itor_loc) reduction(+:mutual) &
!$omp reduction(+:itor)
allocate(rop(self%fe_rep%nce))
allocate(j_lag(self%fe_rep%nce))
!$omp do schedule(static,1)
DO j=1,self%fe_rep%mesh%nc
  IF(self%fe_rep%mesh%reg(j)/=1)CYCLE
  call self%fe_rep%ncdofs(j,j_lag)
  curved=cell_is_curved(self%fe_rep%mesh,j)
  do m=1,self%fe_rep%quad%np
    if(curved.OR.(m==1))call self%fe_rep%mesh%jacobian(j,self%fe_rep%quad%pts(:,m),goptmp,v)
    call psi_eval%interp(j,self%fe_rep%quad%pts(:,m),goptmp,psitmp)
    pt=self%fe_rep%mesh%log2phys(j,self%fe_rep%quad%pts(:,m))
    !---Compute Magnetic Field
    IF(gs_test_bounds(self,pt).AND.psitmp(1)>self%plasma_bounds(1))THEN
      IF(self%mode==0)THEN
        itor_loc = (self%pnorm*pt(1)*self%P%Fp(psitmp(1)) &
        + (self%alam**2)*self%I%Fp(psitmp(1))*(self%I%f(psitmp(1))+self%I%f_offset/self%alam)/(pt(1)+gs_epsilon))
      ELSE
        itor_loc = (self%pnorm*pt(1)*self%P%Fp(psitmp(1)) &
        + .5d0*self%alam*self%I%Fp(psitmp(1))/(pt(1)+gs_epsilon))
      END IF
      b_tmp=0.d0
      DO l=1,self%fe_rep%nce
        CALL oft_blag_eval(self%fe_rep,j,l,self%fe_rep%quad%pts(:,m),rop(l))
        b_tmp=b_tmp+btmp(j_lag(l))*rop(l)
      END DO
      det = v*self%fe_rep%quad%wts(m)
      itor = itor + itor_loc*det
      mutual = mutual + b_tmp*itor_loc*det
    END IF
  end do
end do
deallocate(j_lag,rop)
!$omp end parallel
mutual=mu0*2.d0*pi*mutual/itor
CALL psi_eval%delete()
DEALLOCATE(btmp)
! self%timing(2)=self%timing(2)+(omp_get_wtime()-t1)
end subroutine gs_plasma_mutual
!---------------------------------------------------------------------------
! SUBROUTINE gs_fit_isoflux
!---------------------------------------------------------------------------
!> Compute Grad-Shafranov solution for current flux definitions
!!
!! @param[in,out] self G-S object
!! @param[out] ierr Error flag
!---------------------------------------------------------------------------
subroutine gs_fit_isoflux(self,psi_full,ierr)
class(gs_eq), intent(inout) :: self
class(oft_vector), target, intent(inout) :: psi_full
integer(4), intent(out) :: ierr
type(oft_lag_brinterp) :: psi_eval
type(oft_lag_bginterp) :: psi_geval,psi_geval2
integer(4) :: i,j,k,mind,nCon,io_unit,coffset,roffset
integer(4), allocatable :: cells(:)
real(r8) :: itor,curr,f(3),goptmp(3,4),pol_val(1),v,pt(2),theta,gpsi(3),wt_max,wt_min
real(r8), allocatable :: err_mat(:,:),rhs(:),err_inv(:,:),currs(:),wt_tmp(:)
logical :: pm_save
! !---Read coil constraint matrix
! IF(.NOT.ASSOCIATED(self%coil_reg_mat))THEN
!   ALLOCATE(self%coil_reg_mat(self%ncoil_regs+1,self%ncoil_regs+1))
!   OPEN(NEWUNIT=io_unit,FILE='coil_reg_mat.dat')
!   DO i=1,self%ncoil_regs+1
!     READ(io_unit,*)self%coil_reg_mat(i,:)
!   END DO
!   CLOSE(io_unit)
! END IF
!
nCon = self%isoflux_ntargets+self%flux_ntargets+2*self%saddle_ntargets+self%nregularize
ALLOCATE(err_mat(nCon,self%ncoils+1),err_inv(self%ncoils+1,self%ncoils+1))
ALLOCATE(rhs(nCon),currs(self%ncoils+1))
ALLOCATE(cells(self%isoflux_ntargets+self%flux_ntargets+self%saddle_ntargets))
err_mat=0.d0
rhs=0.d0
cells=-1
psi_eval%u=>self%psi
CALL psi_eval%setup(self%fe_rep)
CALL psi_geval%shared_setup(psi_eval)
ALLOCATE(wt_tmp(self%isoflux_ntargets))
wt_tmp=1.d0
IF(self%isoflux_grad_wt_lim>0.d0)THEN
  psi_geval2%u=>psi_full
  CALL psi_geval2%setup(self%fe_rep)
END IF
!---Get RHS
DO j=1,self%isoflux_ntargets
  CALL bmesh_findcell(self%fe_rep%mesh,cells(j),self%isoflux_targets(1:2,j),f)
  CALL psi_eval%interp(cells(j),f,goptmp,pol_val)
  rhs(j)=pol_val(1)!*self%isoflux_targets(3,j)
  IF(self%isoflux_grad_wt_lim>0.d0)THEN
    CALL self%fe_rep%mesh%jacobian(cells(j),f,goptmp,v)
    CALL psi_geval2%interp(cells(j),f,goptmp,gpsi)
    wt_tmp(j)=SQRT(SUM(gpsi(1:2)**2))
  END IF
END DO
roffset=self%isoflux_ntargets
coffset=self%isoflux_ntargets
DO j=1,self%flux_ntargets
  CALL bmesh_findcell(self%fe_rep%mesh,cells(coffset+j),self%flux_targets(1:2,j),f)
  CALL psi_eval%interp(cells(coffset+j),f,goptmp,pol_val)
  rhs(roffset+j)=pol_val(1)-self%flux_targets(3,j)!*self%isoflux_targets(3,j)
  ! IF(self%isoflux_grad_wt_lim>0.d0)THEN
  !   CALL self%fe_rep%mesh%jacobian(cells(coffset+j),f,goptmp,v)
  !   CALL psi_geval2%interp(cells(coffset+j),f,goptmp,gpsi)
  !   wt_tmp(roffset+j)=SQRT(SUM(gpsi(1:2)**2))
  ! END IF
END DO
roffset=roffset+self%flux_ntargets
coffset=coffset+self%flux_ntargets
DO j=1,self%saddle_ntargets
  CALL bmesh_findcell(self%fe_rep%mesh,cells(coffset+j),self%saddle_targets(1:2,j),f)
  CALL self%fe_rep%mesh%jacobian(cells(coffset+j),f,goptmp,v)
  CALL psi_geval%interp(cells(coffset+j),f,goptmp,gpsi)
  rhs(roffset+2*(j-1)+1)=gpsi(1)*self%saddle_targets(3,j)
  rhs(roffset+2*(j-1)+2)=gpsi(2)*self%saddle_targets(3,j)
END DO
!---Build L-S Matrix
DO i=1,self%ncoils
  psi_eval%u=>self%psi_coil(i)%f
  CALL psi_eval%setup(self%fe_rep)
  psi_geval%u=>self%psi_coil(i)%f
  CALL psi_geval%setup(self%fe_rep)
  DO j=1,self%isoflux_ntargets
    CALL bmesh_findcell(self%fe_rep%mesh,cells(j),self%isoflux_targets(1:2,j),f)
    ! CALL self%fe_rep%mesh%jacobian(cells(j),f,goptmp,v)
    CALL psi_eval%interp(cells(j),f,goptmp,pol_val)
    err_mat(j,i)=pol_val(1)!*self%isoflux_targets(3,j)
  END DO
  roffset=self%isoflux_ntargets
  coffset=self%isoflux_ntargets
  DO j=1,self%flux_ntargets
    CALL bmesh_findcell(self%fe_rep%mesh,cells(coffset+j),self%flux_targets(1:2,j),f)
    ! CALL self%fe_rep%mesh%jacobian(cells(j),f,goptmp,v)
    CALL psi_eval%interp(cells(coffset+j),f,goptmp,pol_val)
    err_mat(roffset+j,i)=pol_val(1)
  END DO
  roffset=roffset+self%flux_ntargets
  coffset=coffset+self%flux_ntargets
  DO j=1,self%saddle_ntargets
    CALL bmesh_findcell(self%fe_rep%mesh,cells(coffset+j),self%saddle_targets(1:2,j),f)
    CALL self%fe_rep%mesh%jacobian(cells(coffset+j),f,goptmp,v)
    CALL psi_geval%interp(cells(coffset+j),f,goptmp,gpsi)
    err_mat(roffset+2*(j-1)+1,i)=gpsi(1)*self%saddle_targets(3,j)
    err_mat(roffset+2*(j-1)+2,i)=gpsi(2)*self%saddle_targets(3,j)
  END DO
END DO
DO i=1,self%ncoils
  err_mat(:,self%ncoils+1)=err_mat(:,self%ncoils+1) &
    + self%coil_vcont(i)*err_mat(:,i)
END DO
!---Enforce difference of fluxes
wt_max=MAXVAL(wt_tmp)
wt_min=self%isoflux_grad_wt_lim*wt_max
DO j=1,self%isoflux_ntargets-1
  err_mat(j+1,:)=(err_mat(j+1,:)-err_mat(1,:))*self%isoflux_targets(3,j+1)*wt_max/MAX(wt_min,wt_tmp(j+1))
  rhs(j+1)=(rhs(j+1)-rhs(1))*self%isoflux_targets(3,j+1)*wt_max/MAX(wt_min,wt_tmp(j+1))
END DO
DEALLOCATE(wt_tmp)
!---Apply weights to flux constraints
roffset=self%isoflux_ntargets
DO j=1,self%flux_ntargets
  err_mat(roffset+j,:)=err_mat(roffset+j,:)*self%flux_targets(4,j)!*wt_max/MAX(wt_min,wt_tmp(j+1))
  rhs(roffset+j)=rhs(roffset+j)*self%flux_targets(4,j)!*wt_max/MAX(wt_min,wt_tmp(j+1))
END DO
!---Coil regularization
roffset=self%isoflux_ntargets+self%flux_ntargets+2*self%saddle_ntargets
DO i=1,self%nregularize
  DO j=1,self%ncoils
    err_mat(roffset+i,j)=self%coil_reg_mat(i,j) !*self%coil_regions(j)%area
  END DO
  err_mat(roffset+i,self%ncoils+1)=self%coil_reg_mat(i,self%ncoils+1)
  rhs(roffset+i)=-self%coil_reg_targets(i)
END DO
!---Solve L-S system
IF(ASSOCIATED(self%coil_bounds))THEN
BLOCK
INTEGER(4) :: nsetp
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: index
REAL(8) :: rnorm
REAL(8), ALLOCATABLE, DIMENSION(:) :: w
ALLOCATE(w(self%ncoils+1),index(self%ncoils+1))
  CALL bvls(ncon-1,self%ncoils+1,err_mat(2:nCon,:),rhs(2:nCon), &
    self%coil_bounds,currs,rnorm,nsetp,w,index,ierr)
  ! WRITE(*,*)ierr,currs
DEALLOCATE(w,index)
END BLOCK
ELSE
  err_inv=MATMUL(TRANSPOSE(err_mat(2:nCon,:)),err_mat(2:nCon,:))
  pm_save=oft_env%pm; oft_env%pm=.FALSE.
  CALL lapack_matinv(self%ncoils+1,err_inv,ierr)
  oft_env%pm=pm_save
  currs=MATMUL(err_inv,MATMUL(TRANSPOSE(err_mat(2:nCon,:)),rhs(2:nCon)))
END IF
!---Add coil/conductor fields to IC
self%vcontrol_val=-currs(self%ncoils+1)
self%coil_currs=-currs(1:self%ncoils)
DEALLOCATE(err_mat,err_inv,rhs,currs,cells)
CALL psi_eval%delete
CALL psi_geval%delete
end subroutine gs_fit_isoflux
!---------------------------------------------------------------------------
! SUBROUTINE gs_fit_walleigs
!---------------------------------------------------------------------------
!> Compute Grad-Shafranov solution for current flux definitions
!!
!! @param[in,out] self G-S object
!! @param[out] ierr Error flag
!---------------------------------------------------------------------------
subroutine gs_fit_walleigs(self,ierr)
class(gs_eq), intent(inout) :: self
integer(4), intent(out) :: ierr
type(oft_lag_brinterp) :: psi_eval
type(oft_lag_bginterp) :: psi_geval
integer(4) :: i,j,k,kk,mind,nCon,io_unit
integer(4), allocatable :: cells(:)
real(r8) :: itor,curr,f(3),goptmp(3,4),pol_val(1),v,pt(2),theta,gpsi(3)
real(r8), allocatable :: err_mat(:,:),rhs(:),err_inv(:,:),currs(:)
logical :: pm_save
! !---Read coil constraint matrix
! IF(.NOT.ASSOCIATED(self%coil_reg_mat))THEN
!   ALLOCATE(self%coil_reg_mat(self%ncoil_regs+1,self%ncoil_regs+1))
!   OPEN(NEWUNIT=io_unit,FILE='coil_reg_mat.dat')
!   DO i=1,self%ncoil_regs+1
!     READ(io_unit,*)self%coil_reg_mat(i,:)
!   END DO
!   CLOSE(io_unit)
! END IF
!
nCon = self%flux_ntargets+self%ncond_eigs
ALLOCATE(err_mat(nCon,self%ncond_eigs),err_inv(self%ncond_eigs,self%ncond_eigs))
ALLOCATE(rhs(nCon),currs(self%ncond_eigs),cells(self%flux_ntargets))
err_mat=0.d0
rhs=0.d0
cells=-1
psi_eval%u=>self%psi
CALL psi_eval%setup(self%fe_rep)
! CALL psi_geval%shared_setup(psi_eval)
!---Get RHS
DO j=1,self%flux_ntargets
  CALL bmesh_findcell(self%fe_rep%mesh,cells(j),self%flux_targets(1:2,j),f)
  ! CALL self%fe_rep%mesh%jacobian(cells(j),f,goptmp,v)
  CALL psi_eval%interp(cells(j),f,goptmp,pol_val)
  rhs(j)=self%flux_targets(3,j)-pol_val(1)
END DO
!---Build L-S Matrix
kk=0
DO i=1,self%ncond_regs
  DO k=1,self%cond_regions(i)%neigs
    kk=kk+1
    err_mat(self%flux_ntargets+kk,kk)=1.E-5 ! Coil regularization
    psi_eval%u=>self%cond_regions(i)%psi_eig(k)%f
    CALL psi_eval%setup(self%fe_rep)
    DO j=1,self%flux_ntargets
      CALL bmesh_findcell(self%fe_rep%mesh,cells(j),self%flux_targets(1:2,j),f)
      ! CALL self%fe_rep%mesh%jacobian(cells(j),f,goptmp,v)
      CALL psi_eval%interp(cells(j),f,goptmp,pol_val)
      err_mat(j,kk)=pol_val(1)
    END DO
    END DO
  END DO
!---Solve L-S system
err_inv=MATMUL(TRANSPOSE(err_mat),err_mat)
pm_save=oft_env%pm; oft_env%pm=.FALSE.
CALL lapack_matinv(self%ncond_eigs,err_inv,ierr)
oft_env%pm=pm_save
CALL gs_get_cond_weights(self,currs,.FALSE.)
! WRITE(*,*)'Before',currs
currs=MATMUL(err_inv,MATMUL(TRANSPOSE(err_mat),rhs))
! WRITE(*,*)'After',currs
!---Add coil/conductor fields to IC
CALL gs_set_cond_weights(self,currs,.FALSE.)
DEALLOCATE(err_mat,err_inv,rhs,currs,cells)
CALL psi_eval%delete
CALL psi_geval%delete
end subroutine gs_fit_walleigs
!---------------------------------------------------------------------------
! SUBROUTINE gs_solve
!---------------------------------------------------------------------------
!> Compute Grad-Shafranov solution for current flux definitions
!!
!! @param[in,out] self G-S object
!! @param[out] ierr Error flag
!---------------------------------------------------------------------------
subroutine gs_solve(self,ierr)
class(gs_eq), intent(inout) :: self
integer(4), optional, intent(out) :: ierr
class(oft_vector), pointer :: rhs,rhs_bc,psip,psiin,psi_bc,psi_eddy,psi_dt
class(oft_vector), pointer :: tmp_vec,psi_alam,psi_press,psi_vac,psi_vcont
real(r8), pointer, DIMENSION(:) :: vals_tmp
type(oft_lag_bginterp), target :: psi_geval
real(8) :: goptmp(3,3),pt(2),v,pmax,pmin,dpnorm,curr
real(8) :: opoint(2),R0_in,f(3),V0_in,V0_tmp,estored
REAL(8) :: nl_res,psimax,alamin,alamp,itor,pnorm0,pnormp,itor_alam,itor_press
REAL(8) :: R0_tmp,R0_hist(2),gpsi0(3),gpsi1(3),gpsi2(3),t0,t1
REAL(8) :: param_mat(3,3),param_vec(3),param_rhs(3)
integer(4) :: i,ii,j,k,error_flag,cell,ierr_loc
integer(4), save :: eq_count = 0
CHARACTER(LEN=40) :: err_reason
logical :: pm_save,fail_test
!---
error_flag=0
IF(self%use_lu)THEN
  IF(.NOT.ASSOCIATED(self%lu_solver%A))THEN
    self%lu_solver%A=>self%dels
    ALLOCATE(self%lu_solver%sec_rhs(self%psi%n,2))
  END IF
ELSE
  CALL oft_abort("LU solver required for GS solve","gs_solve",__FILE__)
END IF
!---
ALLOCATE(vals_tmp(self%psi%n))
CALL self%psi%new(tmp_vec)
CALL self%psi%new(psip)
CALL self%psi%new(psiin)
CALL psiin%add(0.d0,1.d0,self%psi)
CALL psip%add(0.d0,1.d0,self%psi)
alamin=self%alam
alamp=self%alam
!---
CALL self%psi%new(rhs)
CALL self%psi%new(rhs_bc)
CALL self%psi%new(psi_bc)
CALL self%psi%new(psi_vac)
CALL self%psi%new(psi_vcont)
CALL self%psi%new(psi_eddy)
CALL self%psi%new(psi_dt)
CALL self%psi%new(psi_alam)
CALL self%psi%new(psi_press)
t0=omp_get_wtime()
IF(.NOT.self%free)THEN
  CALL rhs_bc%add(0.d0,1.d0,self%psi)
  CALL psi_bc%add(0.d0,1.d0,rhs_bc)
  CALL self%zerob_bc%apply(psi_bc)
  CALL rhs_bc%add(1.d0,-1.d0,psi_bc)
  CALL psi_bc%set(0.d0)
  pm_save=oft_env%pm; oft_env%pm=.FALSE.
  CALL self%lu_solver%apply(psi_bc,rhs_bc)
  oft_env%pm=pm_save
  CALL rhs_bc%set(0.d0)
END IF
!---Update flux functions
self%o_point(1)=-1.d0
CALL gs_update_bounds(self)
CALL self%I%update(self)
CALL self%p%update(self)
!---Get J_phi source term
CALL gs_source(self,self%psi,rhs,psi_alam,psi_press,itor_alam,itor_press,estored)
IF(self%dt>0.d0)THEN
  IF(self%dt/=self%dt_last)THEN
    CALL build_dels(self%dels_dt,self,"free",self%dt)
    self%dt_last=self%dt
    self%lu_solver_dt%refactor=.TRUE.
  END IF
  self%lu_solver_dt%A=>self%dels_dt
END IF
!---Update vacuum field part
CALL psi_vac%set(0.d0)
CALL psi_vac%add(1.d0,1.d0,psi_bc)
DO j=1,self%ncoils
  ! curr = self%coil_regions(j)%curr
  CALL psi_vac%add(1.d0,self%coil_currs(j),self%psi_coil(j)%f)
END DO
!
CALL psi_eddy%set(0.d0)
DO j=1,self%ncond_regs
  DO k=1,self%cond_regions(j)%neigs
    ! ii=self%cond_regions(j)%eig_map(k)
    ! CALL psi_vac%add(1.d0,self%cond_regions(j)%weights(k), &
    !   self%cond_regions(j)%psi_eig(k)%f)
    CALL psi_eddy%add(1.d0,self%cond_regions(j)%weights(k), &
      self%cond_regions(j)%psi_eig(k)%f)
  END DO
END DO
CALL psi_vac%add(1.d0,1.d0,psi_eddy)
!
CALL psi_vcont%set(0.d0)
DO j=1,self%ncoils
  ! curr = self%coil_regions(j)%vcont_gain
  CALL psi_vcont%add(1.d0,self%coil_vcont(j),self%psi_coil(j)%f)
END DO
!---Save input solution
IF(self%save_visit.AND.self%plot_final.AND.(eq_count==0))THEN
  CALL self%xdmf%add_timestep(REAL(eq_count,8))
  CALL self%psi%get_local(vals_tmp)
  IF(self%plasma_bounds(1)<-1.d98)THEN
    CALL self%fe_rep%mesh%save_vertex_scalar(vals_tmp,self%xdmf,'Psi')
  ELSE
    CALL self%fe_rep%mesh%save_vertex_scalar(vals_tmp-self%plasma_bounds(1),self%xdmf,'Psi')
  END IF
  CALL psi_vac%get_local(vals_tmp)
  CALL self%fe_rep%mesh%save_vertex_scalar(vals_tmp,self%xdmf,'Psi_vac')
  CALL psi_eddy%get_local(vals_tmp)
  CALL self%fe_rep%mesh%save_vertex_scalar(vals_tmp,self%xdmf,'Psi_eddy')
  CALL psi_vcont%get_local(vals_tmp)
  vals_tmp=vals_tmp*self%vcontrol_val
  CALL self%fe_rep%mesh%save_vertex_scalar(vals_tmp,self%xdmf,'Psi_vcont')
END IF
!---
nl_res=1.d99
pnorm0 = self%pnorm
pnormp = self%pnorm
R0_in = self%o_point(1)
V0_in = self%o_point(2)
cell=0
IF(self%R0_target>0.d0)THEN
  IF((self%estore_target>0.d0).OR.(self%pax_target>0.d0).OR.(self%Ip_ratio_target>-1.d98))THEN
    CALL oft_warn("Conflicting pressure targets specified, ignoring R0_target")
    self%R0_target=-1.d0
  END IF
END IF
IF((self%estore_target>0.d0).AND.((self%pax_target>0.d0).OR.(self%Ip_ratio_target>-1.d98)))THEN
  CALL oft_warn("Conflicting pressure targets specified, ignoring estore_target")
  self%estore_target=-1.d0
END IF
IF((self%pax_target>0.d0).AND.(self%Ip_ratio_target>-1.d98))THEN
  CALL oft_warn("Conflicting pressure targets specified, ignoring pax_target")
  self%pax_target=-1.d0
END IF
IF(self%V0_target>-1.d98)THEN
  IF(self%isoflux_ntargets>0.OR.self%flux_ntargets>0)THEN
    CALL oft_warn("V0_target and isoflux_targets specified, ignoring V0_target")
    self%V0_target=-1.d99
    self%vcontrol_val=0.d0
  END IF
END IF
!---
IF(oft_env%pm)THEN
  WRITE(*,'(2A)')oft_indent,'Starting non-linear GS solver'
  CALL oft_increase_indent
END IF
DO i=1,self%maxits
  !---Ramp R0 target
  R0_tmp=(i-1)*(self%R0_target-R0_in)/REAL(self%nR0_ramp,8) + R0_in
  V0_tmp=(i-1)*(self%V0_target-V0_in)/REAL(self%nR0_ramp,8) + V0_in
  IF(i>self%nR0_ramp)R0_tmp=self%R0_target
  IF(i>self%nR0_ramp)V0_tmp=self%V0_target
  !---
  CALL psip%add(0.d0,1.d0,self%psi)

  !---Compute toroidal flux contribution
  CALL rhs%add(0.d0,1.d0,psi_alam)
  CALL self%zerob_bc%apply(rhs)
  !---Solve linear system
  CALL rhs%get_local(vals_tmp)
  self%lu_solver%sec_rhs(:,1) = vals_tmp
  !---Compute pressure contribution
  CALL rhs%add(0.d0,1.d0,psi_press)
  CALL self%zerob_bc%apply(rhs)
  pm_save=oft_env%pm; oft_env%pm=.FALSE.
  t1=omp_get_wtime()
  self%lu_solver%nrhs=2
  CALL psi_press%set(0.d0)
  CALL self%lu_solver%apply(psi_press,rhs)
  CALL psi_alam%restore_local(self%lu_solver%sec_rhs(:,1))
  CALL psi_alam%scale(1.d0/self%alam)
  self%lu_solver%nrhs=1
  self%timing(3)=self%timing(3)+(omp_get_wtime()-t1)
  oft_env%pm=pm_save

  param_mat=0.d0
  param_rhs=0.d0
  IF(self%Itor_target>0.d0)THEN
    ! itor_alam=self%itor(psi_alam)
    ! itor_press=self%itor(psi_press)
    param_mat(1,:)=[itor_alam,itor_press,0.d0]
    param_rhs(1)=self%Itor_target
  ELSE
    param_mat(1,1)=1.d0
    param_rhs(1)=self%alam
  END IF

  !---Get desired O-point location for linear fit
  cell=0
  pt=self%o_point
  IF(self%R0_target>0.d0)pt(1)=R0_tmp
  IF(self%V0_target>-1.d98)pt(2)=V0_tmp
  CALL bmesh_findcell(self%fe_rep%mesh,cell,pt,f)
  CALL self%fe_rep%mesh%jacobian(cell,f,goptmp,v)

  !---Add row for radial control (beta)
  IF(self%R0_target>0.d0)THEN
    !
    psi_geval%u=>psi_vac
    CALL psi_geval%setup(self%fe_rep)
    CALL psi_geval%interp(cell,f,goptmp,gpsi0)
    param_rhs(2)=-gpsi0(1)
    !
    psi_geval%u=>psi_vcont
    CALL psi_geval%setup(self%fe_rep)
    CALL psi_geval%interp(cell,f,goptmp,gpsi0)
    psi_geval%u=>psi_alam
    CALL psi_geval%setup(self%fe_rep)
    CALL psi_geval%interp(cell,f,goptmp,gpsi1)
    psi_geval%u=>psi_press
    CALL psi_geval%setup(self%fe_rep)
    CALL psi_geval%interp(cell,f,goptmp,gpsi2)
    param_mat(2,:)=[gpsi1(1),gpsi2(1),gpsi0(1)]
  ELSE IF(self%estore_target>0.d0)THEN
    param_rhs(2)=self%estore_target
    param_mat(2,2)=estored*3.d0/2.d0
  ELSE IF(self%pax_target>0.d0)THEN
    param_rhs(2)=self%pax_target
    param_mat(2,2)=self%P%f(self%plasma_bounds(2))
  ELSE IF(self%Ip_ratio_target>-1.d98)THEN
    param_rhs(2)=0.d0
    param_mat(2,:)=[itor_alam,-itor_press*self%Ip_ratio_target,0.d0]
  ELSE
    param_mat(2,2)=1.d0
    param_rhs(2)=self%pnorm
  END IF

  !---Add row for vertical control
  IF(self%V0_target>-1.d98)THEN
    ! 
    CALL bmesh_findcell(self%fe_rep%mesh,cell,pt,f)
    CALL self%fe_rep%mesh%jacobian(cell,f,goptmp,v)
    psi_geval%u=>psi_vac
    CALL psi_geval%setup(self%fe_rep)
    CALL psi_geval%interp(cell,f,goptmp,gpsi0)
    param_rhs(3)=-gpsi0(2)
    ! 
    psi_geval%u=>psi_vcont
    CALL psi_geval%setup(self%fe_rep)
    CALL psi_geval%interp(cell,f,goptmp,gpsi0)
    psi_geval%u=>psi_alam
    CALL psi_geval%setup(self%fe_rep)
    CALL psi_geval%interp(cell,f,goptmp,gpsi1)
    psi_geval%u=>psi_press
    CALL psi_geval%setup(self%fe_rep)
    CALL psi_geval%interp(cell,f,goptmp,gpsi2)
    param_mat(3,:)=[gpsi1(2),gpsi2(2),gpsi0(2)]
  ELSE
    param_mat(3,3)=1.d0
    param_rhs(3)=0.d0
  END IF

  ! Solve for parameters
  pm_save=oft_env%pm; oft_env%pm=.FALSE.
  CALL lapack_matinv(3,param_mat,ierr_loc)
  oft_env%pm=pm_save
  IF(ierr_loc/=0)THEN
    error_flag=-6
    EXIT
  END IF
  param_vec=MATMUL(param_mat,param_rhs)
  self%alam=param_vec(1)
  self%pnorm=param_vec(2)
  self%vcontrol_val=param_vec(3)

  ! Create plasma poloidal flux
  CALL self%psi%set(0.d0)
  CALL self%psi%add(0.d0,self%alam,psi_alam,self%pnorm,psi_press)

  ! Fit coils if operating in isoflux mode
  IF(self%isoflux_ntargets+self%flux_ntargets+self%saddle_ntargets>0)THEN
    CALL self%psi%add(1.d0,1.d0,psi_eddy)
    CALL self%psi%add(1.d0,1.d0,psi_bc)
    CALL gs_fit_isoflux(self,psip,ierr_loc)
    IF(ierr_loc/=0)THEN
      error_flag=-7
      EXIT
    END IF
    CALL self%psi%add(1.d0,-1.d0,psi_eddy)
    CALL self%psi%add(1.d0,-1.d0,psi_bc)
    !---Update vacuum field part
    CALL psi_vac%set(0.d0)
    CALL psi_vac%add(1.d0,1.d0,psi_bc)
    DO j=1,self%ncoils
      ! curr = self%coil_regions(j)%curr
      CALL psi_vac%add(1.d0,self%coil_currs(j),self%psi_coil(j)%f)
    END DO
    ! CALL psi_eddy%set(0.d0)
    ! DO j=1,self%ncond_regs
    !   DO k=1,self%cond_regions(j)%neigs
    !     ! ii=self%cond_regions(j)%eig_map(k)
    !     CALL psi_vac%add(1.d0,self%cond_regions(j)%weights(k), &
    !       self%cond_regions(j)%psi_eig(k)%f)
    !     CALL psi_eddy%add(1.d0,self%cond_regions(j)%weights(k), &
    !       self%cond_regions(j)%psi_eig(k)%f)
    !   END DO
    ! END DO
  END IF

  ! Fit wall eigenmodes if fluxes are available
  ! IF(self%flux_ntargets>0)THEN
  !   CALL self%psi%add(1.d0,1.d0,psi_vac)
  !   CALL gs_fit_walleigs(self,ierr_loc)
  !   IF(ierr_loc/=0)THEN
  !     error_flag=-8
  !     EXIT
  !   END IF
  !   CALL self%psi%add(1.d0,-1.d0,psi_vac)
  !   !---Update vacuum field part
  !   CALL psi_eddy%set(0.d0)
  !   DO j=1,self%ncond_regs
  !     DO k=1,self%cond_regions(j)%neigs
  !       ! ii=self%cond_regions(j)%eig_map(k)
  !       ! CALL psi_vac%add(1.d0,self%cond_regions(j)%weights(k), &
  !       !   self%cond_regions(j)%psi_eig(k)%f)
  !       CALL psi_eddy%add(1.d0,self%cond_regions(j)%weights(k), &
  !         self%cond_regions(j)%psi_eig(k)%f)
  !     END DO
  !   END DO
  !   CALL psi_vac%add(1.d0,1.d0,psi_eddy)
  ! ELSE IF(self%isoflux_ntargets>0)THEN
  IF(self%isoflux_ntargets+self%flux_ntargets+self%saddle_ntargets>0)THEN
    !---Update vacuum field part
    CALL psi_eddy%set(0.d0)
    DO j=1,self%ncond_regs
      DO k=1,self%cond_regions(j)%neigs
        ! ii=self%cond_regions(j)%eig_map(k)
        ! CALL psi_vac%add(1.d0,self%cond_regions(j)%weights(k), &
        !   self%cond_regions(j)%psi_eig(k)%f)
        CALL psi_eddy%add(1.d0,self%cond_regions(j)%weights(k), &
          self%cond_regions(j)%psi_eig(k)%f)
      END DO
    END DO
    CALL psi_vac%add(1.d0,1.d0,psi_eddy)
  END IF

  ! Add vacuum fields to solution
  CALL self%psi%add(1.d0,1.d0,psi_vac)
  CALL self%psi%add(1.d0,self%vcontrol_val,psi_vcont)

  ! Compute passive eddy currents
  IF(self%dt>0.d0)THEN
    CALL psi_dt%set(0.d0)
    CALL psi_dt%add(0.d0,-1.d0/self%dt,self%psi,1.d0/self%dt,self%psi_dt)
    CALL gs_wall_source(self,psi_dt,tmp_vec)
    CALL self%zerob_bc%apply(tmp_vec)
    pm_save=oft_env%pm; oft_env%pm=.FALSE.
    CALL self%lu_solver_dt%apply(psi_dt,tmp_vec)
    oft_env%pm=pm_save
    CALL psi_eddy%add(1.d0,1.d0,psi_dt)
    CALL psi_vac%add(1.d0,1.d0,psi_dt)
    CALL self%psi%add(1.d0,1.d0,psi_dt)
  END IF

  !---Check for scale issues
  CALL self%psi%get_local(vals_tmp)
  psimax=maxval(ABS(vals_tmp))
  fail_test=psimax<gs_epsilon
  fail_test=fail_test.OR.((self%plasma_bounds(2) < self%plasma_bounds(1)).AND.self%has_plasma)
  fail_test=fail_test.OR.((self%o_point(1) < self%rmin).AND.self%has_plasma)
  IF(fail_test)THEN
    !WRITE(*,*)psimax,self%plasma_bounds,self%o_point
    IF(psimax<gs_epsilon)error_flag=-2
    IF((self%plasma_bounds(2) < self%plasma_bounds(1)).AND.self%has_plasma)error_flag=-3
    IF((self%o_point(1) < self%rmin).AND.self%has_plasma)error_flag=-4
    EXIT
  END IF
  !---Under-relax solution
  ! IF(MOD(i,self%ninner)/=0.OR.(.NOT.self%free))
  CALL self%psi%add(1.d0-self%urf,self%urf,psip)
  !---Update flux scale for free and fixed boundary
  ! CALL self%psi%add(1.d0,-1.d0,psi_vac)
  ! CALL self%psi%add(1.d0,-self%vcontrol_val,psi_vcont)
  ! CALL self%zerob_bc%apply(self%psi)
  IF(self%free)THEN
    ! CALL self%psi%add(1.d0,1.d0,psi_bc)
    ! CALL self%psi%add(1.d0,1.d0,psi_vac)
    ! CALL self%psi%add(1.d0,self%vcontrol_val,psi_vcont)
  ELSE
    CALL self%psi%add(1.d0,-1.d0,psi_vac)
    CALL self%psi%add(1.d0,-self%vcontrol_val,psi_vcont)
    CALL self%zerob_bc%apply(self%psi)
    IF(self%I%f_offset==0.d0)THEN
      CALL self%psi%get_local(vals_tmp)
      self%psimax=MAXVAL(vals_tmp)
      CALL self%psi%scale(1.d0/self%psimax)
      self%alam=SQRT((self%alam**2)/self%psimax)
      self%psimax=1.d0
    ! ELSE IF(self%Itor_target>0.d0)THEN
    !   itor=self%itor()
    !   IF(itor<=0.d0)THEN
    !     error_flag=-5
    !     EXIT
    !   END IF
    !   alamp=itor/self%Itor_target
    !   CALL self%psi%scale(1.d0/alamp)
    !   self%alam=SQRT((self%alam**2)/alamp)
    END IF
    CALL self%psi%add(1.d0,1.d0,psi_vac)
    CALL self%psi%add(1.d0,self%vcontrol_val,psi_vcont)
    alamp=self%alam
  END IF
  !---Update flux functions
  CALL gs_update_bounds(self)
  CALL self%I%update(self)
  CALL self%p%update(self)
  !---Output
  IF(self%save_visit.AND.self%plot_step)THEN
    eq_count=eq_count+1
    CALL self%xdmf%add_timestep(REAL(eq_count,8))
    CALL self%psi%get_local(vals_tmp)
    IF(self%plasma_bounds(1)<-1.d98)THEN
      CALL self%fe_rep%mesh%save_vertex_scalar(vals_tmp,self%xdmf,'Psi')
    ELSE
      CALL self%fe_rep%mesh%save_vertex_scalar(vals_tmp-self%plasma_bounds(1),self%xdmf,'Psi')
    END IF
    CALL psi_vac%get_local(vals_tmp)
    CALL self%fe_rep%mesh%save_vertex_scalar(vals_tmp,self%xdmf,'Psi_vac')
    CALL psi_eddy%get_local(vals_tmp)
    CALL self%fe_rep%mesh%save_vertex_scalar(vals_tmp,self%xdmf,'Psi_eddy')
    CALL psi_vcont%get_local(vals_tmp)
    vals_tmp=vals_tmp*self%vcontrol_val
    CALL self%fe_rep%mesh%save_vertex_scalar(vals_tmp,self%xdmf,'Psi_vcont')
  END IF
  ! !---Under-relax pressure in R0 control mode
  ! IF(self%R0_target>0.d0.AND.MOD(i,self%ninner)==0)THEN
  !   self%pnorm=(self%pnorm+pnormp)/2.d0
  !   pnormp=self%pnorm
  ! END IF
  !---Update vacuum field part
  CALL gs_source(self,self%psi,rhs,psi_alam,psi_press,itor_alam,itor_press,estored)
  !---Compute error in NL function
  CALL tmp_vec%set(0.d0)
  CALL tmp_vec%add(0.d0,1.d0,self%psi,-1.d0,psi_vac)
  CALL tmp_vec%add(1.d0,-self%vcontrol_val,psi_vcont)
  CALL self%dels%apply(tmp_vec,psip)
  CALL psip%add(1.d0,-1.d0,rhs)
  IF(.NOT.self%free)CALL self%zerob_bc%apply(psip)
  nl_res=psip%dot(psip)
  !---Output progress
  IF(oft_env%pm)WRITE(*,'(A,I4,6ES12.4)')oft_indent,i,self%alam,self%pnorm, &
    SQRT(nl_res),self%o_point(1),self%o_point(2),self%vcontrol_val/mu0
  !---Check if converged
  IF((self%R0_target>0.d0).AND.(i<self%nR0_ramp))CYCLE
  IF(SQRT(nl_res)<self%nl_tol)EXIT !.AND.i>3*self%ninner)EXIT
end do
IF(oft_env%pm)CALL oft_decrease_indent
IF(i>self%maxits)error_flag=-1
!---Output
IF(self%save_visit.AND.self%plot_final)THEN
  eq_count=eq_count+1
  CALL self%xdmf%add_timestep(REAL(eq_count,8))
  CALL self%psi%get_local(vals_tmp)
  IF(self%plasma_bounds(1)<-1.d98)THEN
    CALL self%fe_rep%mesh%save_vertex_scalar(vals_tmp,self%xdmf,'Psi')
  ELSE
    CALL self%fe_rep%mesh%save_vertex_scalar(vals_tmp-self%plasma_bounds(1),self%xdmf,'Psi')
  END IF
  CALL psi_vac%get_local(vals_tmp)
  CALL self%fe_rep%mesh%save_vertex_scalar(vals_tmp,self%xdmf,'Psi_vac')
  CALL psi_eddy%get_local(vals_tmp)
  CALL self%fe_rep%mesh%save_vertex_scalar(vals_tmp,self%xdmf,'Psi_eddy')
  CALL psi_vcont%get_local(vals_tmp)
  vals_tmp=vals_tmp*self%vcontrol_val
  CALL self%fe_rep%mesh%save_vertex_scalar(vals_tmp,self%xdmf,'Psi_vcont')
END IF
self%timing(1)=self%timing(1)+(omp_get_wtime()-t0)
IF(oft_env%pm)THEN
  WRITE(*,*)'Timing:',self%timing(1)
  WRITE(*,*)'  Source:  ',self%timing(2)
  WRITE(*,*)'  Solve:   ',self%timing(3)
  WRITE(*,*)'  Boundary:',self%timing(4)
  WRITE(*,*)'  Other:   ',self%timing(1)-SUM(self%timing(2:4))
END IF
!---
CALL rhs%delete
CALL psip%delete
CALL psiin%delete
CALL rhs_bc%delete
CALL psi_bc%delete
CALL psi_vac%delete
CALL psi_vcont%delete
CALL psi_eddy%delete
CALL psi_alam%delete
CALL psi_press%delete
DEALLOCATE(rhs,psip,psiin)
DEALLOCATE(rhs_bc,psi_bc,psi_alam,psi_press,psi_vac,psi_vcont)
DEALLOCATE(vals_tmp,psi_eddy)
!---
IF(self%compute_chi)CALL self%get_chi
self%ierr=error_flag
IF(PRESENT(ierr))THEN
  ierr=error_flag
ELSE
  IF(error_flag<0)THEN
    err_reason=gs_err_reason(error_flag)
    WRITE(*,'(3A)')oft_indent,'Equilibrium solve Failed: ',TRIM(err_reason)
  END IF
END IF
end subroutine gs_solve
!---------------------------------------------------------------------------
! SUBROUTINE gs_lin_solve
!---------------------------------------------------------------------------
!> Compute Grad-Shafranov solution for current flux definitions
!!
!! @param[in,out] self G-S object
!! @param[out] ierr Error flag
!---------------------------------------------------------------------------
subroutine gs_lin_solve(self,adjust_r0,ierr)
class(gs_eq), intent(inout) :: self
logical, intent(in) :: adjust_r0
integer(4), optional, intent(out) :: ierr
class(oft_vector), pointer :: rhs,rhs_bc,psip,psiin,psi_bc,psi_alam,psi_press
class(oft_vector), pointer :: psi_vac,psi_vcont
real(r8), pointer, dimension(:) :: vals_tmp
type(oft_lag_bginterp), target :: psi_geval
real(8) :: goptmp(3,3),pt(2),v,pmax,pmin,dpnorm
real(8) :: opoint(2),f(3),t0,t1,curr,itor_alam,itor_press,estored
REAL(8) :: psimax,alamin,alamp,pnorm0,R0_hist(2),gpsi0(3),gpsi1(3),gpsi2(3)
REAL(8) :: param_mat(3,3),param_vec(3),param_rhs(3)
integer(4) :: ii,j,k,error_flag,cell,ierr_loc
CHARACTER(LEN=40) :: err_reason
logical :: pm_save
t0=omp_get_wtime()
!---
error_flag=0
!---
ALLOCATE(vals_tmp(self%psi%n))
CALL self%psi%new(psip)
CALL self%psi%new(psiin)
CALL psiin%add(0.d0,1.d0,self%psi)
CALL psip%add(0.d0,1.d0,self%psi)
alamin=self%alam
alamp=self%alam
!---
CALL self%psi%new(rhs)
CALL self%psi%new(rhs_bc)
CALL self%psi%new(psi_bc)
CALL self%psi%new(psi_vac)
CALL self%psi%new(psi_vcont)
CALL self%psi%new(psi_alam)
CALL self%psi%new(psi_press)
IF(oft_env%pm)THEN
  WRITE(*,'(2A)')oft_indent,'Starting Linearized GS solver'
  CALL oft_increase_indent
END IF
!---Update flux functions
CALL gs_update_bounds(self)
CALL self%I%update(self)
CALL self%p%update(self)
!---Get J_phi source term
CALL gs_source(self,self%psi,rhs,psi_alam,psi_press,itor_alam,itor_press,estored)
CALL psi_alam%scale(1.d0/self%alam)
!---Update vacuum field part
CALL psi_vac%set(0.d0)
DO j=1,self%ncoils
  ! curr = self%coil_regions(j)%curr
  CALL psi_vac%add(1.d0,self%coil_currs(j),self%psi_coil(j)%f)
END DO
DO j=1,self%ncond_regs
  DO k=1,self%cond_regions(j)%neigs
    ii=self%cond_regions(j)%eig_map(k)
    CALL psi_vac%add(1.d0,self%cond_regions(j)%weights(k), &
      self%cond_regions(j)%psi_eig(k)%f)
  END DO
END DO
CALL psi_vcont%set(0.d0)
DO j=1,self%ncoils
  ! curr = self%coil_regions(j)%vcont_gain
  CALL psi_vcont%add(1.d0,self%coil_vcont(j),self%psi_coil(j)%f)
END DO

param_mat=0.d0
param_rhs=0.d0
IF(self%Itor_target>0.d0)THEN
  itor_alam=self%itor(psi_alam)
  itor_press=self%itor(psi_press)
  param_mat(1,:)=[itor_alam,itor_press,0.d0]
  param_rhs(1)=self%Itor_target
ELSE
  param_mat(1,1)=1.d0
  param_rhs(1)=self%alam
END IF

!---Get desired O-point location for linear fit
cell=0
pt=self%o_point
IF(self%R0_target>0.d0)pt(1)=self%R0_target
IF(self%V0_target>-1.d98)pt(2)=self%V0_target
CALL bmesh_findcell(self%fe_rep%mesh,cell,pt,f)
CALL self%fe_rep%mesh%jacobian(cell,f,goptmp,v)

!---Add row for radial control (beta)
IF((self%R0_target>0.d0).AND.adjust_r0)THEN
  !
  psi_geval%u=>psi_vac
  CALL psi_geval%setup(self%fe_rep)
  CALL psi_geval%interp(cell,f,goptmp,gpsi0)
  param_rhs(2)=-gpsi0(1)
  !
  psi_geval%u=>psi_vcont
  CALL psi_geval%setup(self%fe_rep)
  CALL psi_geval%interp(cell,f,goptmp,gpsi0)
  psi_geval%u=>psi_alam
  CALL psi_geval%setup(self%fe_rep)
  CALL psi_geval%interp(cell,f,goptmp,gpsi1)
  psi_geval%u=>psi_press
  CALL psi_geval%setup(self%fe_rep)
  CALL psi_geval%interp(cell,f,goptmp,gpsi2)
  param_mat(2,:)=[gpsi1(1),gpsi2(1),gpsi0(1)]
ELSE
  param_mat(2,2)=1.d0
  param_rhs(2)=self%pnorm
END IF

!---Add row for vertical control
IF((self%V0_target>-1.d98).AND.adjust_r0)THEN
  ! 
  CALL bmesh_findcell(self%fe_rep%mesh,cell,pt,f)
  CALL self%fe_rep%mesh%jacobian(cell,f,goptmp,v)
  psi_geval%u=>psi_vac
  CALL psi_geval%setup(self%fe_rep)
  CALL psi_geval%interp(cell,f,goptmp,gpsi0)
  param_rhs(3)=-gpsi0(2)
  ! 
  psi_geval%u=>psi_vcont
  CALL psi_geval%setup(self%fe_rep)
  CALL psi_geval%interp(cell,f,goptmp,gpsi0)
  psi_geval%u=>psi_alam
  CALL psi_geval%setup(self%fe_rep)
  CALL psi_geval%interp(cell,f,goptmp,gpsi1)
  psi_geval%u=>psi_press
  CALL psi_geval%setup(self%fe_rep)
  CALL psi_geval%interp(cell,f,goptmp,gpsi2)
  param_mat(3,:)=[gpsi1(2),gpsi2(2),gpsi0(2)]
ELSE
  param_mat(3,3)=1.d0
  param_rhs(3)=0.d0
END IF

! Solve for parameters
pm_save=oft_env%pm; oft_env%pm=.FALSE.
CALL lapack_matinv(3,param_mat,ierr_loc)
oft_env%pm=pm_save
param_vec=MATMUL(param_mat,param_rhs)
self%alam=param_vec(1)
self%pnorm=param_vec(2)
self%vcontrol_val=param_vec(3)

! Create plasma poloidal flux
CALL self%psi%set(0.d0)
CALL self%psi%add(0.d0,self%alam,psi_alam,self%pnorm,psi_press)
CALL self%psi%add(1.d0,1.d0,psi_vac)
CALL self%psi%add(1.d0,self%vcontrol_val,psi_vcont)


! !---Compute initial guess based on zeroing r-gradient
! IF(self%R0_target>0.d0.AND.adjust_r0)THEN
!   !---Compute toroidal flux contribution
!   CALL rhs%add(0.d0,1.d0,psi_alam)
!   IF(self%free)THEN ! Set BC for dirichlet flux
!     CALL psi_bc%set(0.d0)
!     CALL gs_set_bc(self,self%u_hom,psi_bc)
!     CALL self%zerob_bc%apply(rhs)
!     CALL rhs%add(1.d0,1.d0,psi_bc)
!   ELSE
!     CALL self%zerob_bc%apply(rhs)
!   END IF
!   !---Solve linear system
!   CALL rhs%get_local(vals_tmp)
!   self%lu_solver%sec_rhs(:,1) = vals_tmp
!   !---Compute pressure contribution
!   CALL rhs%add(0.d0,1.d0,psi_press)
!   CALL self%zerob_bc%apply(rhs)
!   pm_save=oft_env%pm; oft_env%pm=.FALSE.
!   t1=omp_get_wtime()
!   self%lu_solver%nrhs=2
!   CALL self%lu_solver%apply(psi_press,rhs)
!   CALL psi_alam%restore_local(self%lu_solver%sec_rhs(:,1))
!   self%lu_solver%nrhs=1
!   self%timing(3)=self%timing(3)+(omp_get_wtime()-t1)
!   oft_env%pm=pm_save
!   !---Update pressure with fixed R0 using bissection
!   pmin=-1.d99
!   pmax=-1.d99
!   pnorm0=self%pnorm
!   dpnorm=.1d0
!   R0_hist=0.d0
!   R0_tmp=self%R0_target
!   opoint=self%o_point
!   DO j=1,30
!     !---Compute initial guess based on zeroing r-gradient
!     IF(j==1)THEN
!       pt=(/R0_tmp,self%o_point(2)/)
!       cell=0
!       CALL bmesh_findcell(self%fe_rep%mesh,cell,pt,f)
!       CALL self%fe_rep%mesh%jacobian(cell,f,goptmp,v)
!       !
!       psi_geval%u=>psi_alam
!       CALL psi_geval%setup(self%fe_rep)
!       CALL psi_geval%interp(cell,f,goptmp,gpsi1)
!       psi_geval%u=>psi_press
!       CALL psi_geval%setup(self%fe_rep)
!       CALL psi_geval%interp(cell,f,goptmp,gpsi2)
!       IF(gpsi2(1)/=0.d0)self%pnorm=-gpsi1(1)/gpsi2(1)
!       opoint=pt
!     END IF
!     !
!     CALL self%psi%set(0.d0)
!     CALL self%psi%add(0.d0,1.d0,psi_alam,self%pnorm,psi_press)
!     IF(self%R0_target<0.d0)EXIT
!     !---Compute R0
!     CALL gs_psimax(self,psimax,opoint(1),opoint(2))
!     IF(ABS(opoint(1)-R0_tmp)<1.d-4)EXIT
!     !---Compute next guess for pnorm
!     IF(opoint(1)>R0_tmp)THEN
!       R0_hist(2)=opoint(1)
!       pmax=self%pnorm
!     END IF
!     IF(opoint(1)<R0_tmp)THEN
!       R0_hist(1)=opoint(1)
!       pmin=self%pnorm
!     END IF
!     IF(pmin>-1.d98.AND.pmax>-1.d98)THEN
!       self%pnorm = (R0_tmp-R0_hist(1))*(pmax-pmin)/(R0_hist(2)-R0_hist(1)) + pmin
!       CYCLE
!     END IF
!     IF(opoint(1)>R0_tmp)self%pnorm = self%pnorm - dpnorm
!     IF(opoint(1)<R0_tmp)self%pnorm = self%pnorm + dpnorm
!   END DO
! ELSE
!   !---Compute full contribution
!   CALL rhs%add(0.d0,1.d0,psi_alam,self%pnorm,psi_press)
!   IF(self%free)THEN ! Set BC for dirichlet flux
!     CALL psi_bc%set(0.d0)
!     CALL gs_set_bc(self,self%u_hom,psi_bc)
!     IF(self%use_lu)THEN
!       CALL self%zerob_bc%apply(rhs)
!       CALL rhs%add(1.d0,1.d0,psi_bc)
!     ! ELSE
!     !   CALL self%solver%a%apply(psi_bc,rhs_bc)
!     !   CALL rhs%add(1.d0,-1.d0,rhs_bc)
!     END IF
!   ELSE
!     CALL self%zerob_bc%apply(rhs)
!   END IF
!   !---Solve linear system
!   t1=omp_get_wtime()
!   IF(self%use_lu)THEN
!     CALL self%lu_solver%apply(self%psi,rhs)
!   ! ELSE
!   !   pm_save=oft_env%pm; oft_env%pm=.FALSE.
!   !   CALL self%psi%add(0.d0,1.d0,psi_bc)
!   !   CALL self%solver%apply(self%psi,rhs)
!   !   oft_env%pm=pm_save
!   END IF
!   self%timing(3)=self%timing(3)+(omp_get_wtime()-t1)
! END IF
IF(oft_env%pm)CALL oft_decrease_indent
! !---Update flux functions
! CALL gs_update_bounds(self)
! CALL self%I%update(self)
! CALL self%p%update(self)
!---
CALL rhs%delete
CALL psip%delete
CALL psiin%delete
CALL rhs_bc%delete
CALL psi_bc%delete
CALL psi_vac%delete
CALL psi_vcont%delete
CALL psi_alam%delete
CALL psi_press%delete
DEALLOCATE(rhs,psip,psiin)
DEALLOCATE(rhs_bc,psi_bc,psi_alam,psi_press,psi_vac,psi_vcont)
DEALLOCATE(vals_tmp)
!---
IF(self%compute_chi)CALL self%get_chi
self%ierr=error_flag
self%timing(1)=self%timing(1)+(omp_get_wtime()-t0)
IF(PRESENT(ierr))THEN
  ierr=error_flag
ELSE
  IF(error_flag<0)THEN
    err_reason=gs_err_reason(error_flag)
    WRITE(*,'(3A)')oft_indent,'Equilibrium solve Failed: ',TRIM(err_reason)
  END IF
END IF
end subroutine gs_lin_solve
!---------------------------------------------------------------------------
!> Compute Grad-Shafranov solution for vacuum (no plasma)
!---------------------------------------------------------------------------
subroutine gs_vac_solve(self,psi_sol,rhs_source,ierr)
class(gs_eq), intent(inout) :: self !< G-S object
class(oft_vector), intent(inout) :: psi_sol !< Input: BCs for \f$ \psi \f$, Output: solution
CLASS(bfem_interp), optional, intent(inout) :: rhs_source !< Specified current source (optional)
integer(4), optional, intent(out) :: ierr !< Error flag
class(oft_vector), pointer :: rhs_bc,psi_bc,psi_eddy,psi_dt,tmp_vec,psi_vac,psi_vcont
integer(4) :: j,k,error_flag
REAL(8) :: psimax
logical :: pm_save
!---
ierr=0
IF(self%use_lu)THEN
  IF(.NOT.ASSOCIATED(self%lu_solver%A))THEN
    self%lu_solver%A=>self%dels
    ALLOCATE(self%lu_solver%sec_rhs(self%psi%n,2))
  END IF
ELSE
  CALL oft_abort("LU solver required for GS solve","gs_solve",__FILE__)
END IF
!
! self%o_point(1)=-1.d0
CALL self%psi%new(psi_vac)
CALL self%psi%new(psi_vcont)
CALL self%psi%new(psi_eddy)
IF(self%dt>0.d0)THEN
  IF(self%dt/=self%dt_last)THEN
    CALL build_dels(self%dels_dt,self,"free",self%dt)
    self%dt_last=self%dt
    self%lu_solver_dt%refactor=.TRUE.
  END IF
  self%lu_solver_dt%A=>self%dels_dt
  CALL self%psi%new(psi_dt)
  CALL self%psi%new(tmp_vec)
ELSE
  IF(PRESENT(rhs_source))THEN
    CALL self%psi%new(psi_dt)
    CALL self%psi%new(tmp_vec)
  END IF
END IF
!---Update vacuum field part
CALL psi_vac%set(0.d0)
DO j=1,self%ncoils
  ! curr = self%coil_regions(j)%curr
  CALL psi_vac%add(1.d0,self%coil_currs(j),self%psi_coil(j)%f)
END DO
CALL psi_eddy%set(0.d0)
DO j=1,self%ncond_regs
  DO k=1,self%cond_regions(j)%neigs
    ! ii=self%cond_regions(j)%eig_map(k)
    ! CALL psi_vac%add(1.d0,self%cond_regions(j)%weights(k), &
    !   self%cond_regions(j)%psi_eig(k)%f)
    CALL psi_eddy%add(1.d0,self%cond_regions(j)%weights(k), &
    self%cond_regions(j)%psi_eig(k)%f)
  END DO
END DO
CALL psi_vac%add(1.d0,1.d0,psi_eddy)
!
CALL psi_vcont%set(0.d0)
DO j=1,self%ncoils
  ! curr = self%coil_regions(j)%vcont_gain
  CALL psi_vcont%add(1.d0,self%coil_vcont(j),self%psi_coil(j)%f)
END DO
IF((self%dt>0.d0).AND.oft_env%pm)THEN
  WRITE(*,'(2A)')oft_indent,'Starting vacuum GS solver'
  CALL oft_increase_indent
END IF
! Compute inhomogeneous part
psimax=psi_sol%dot(psi_sol)
IF(psimax>1.d-14)THEN
  CALL self%psi%new(rhs_bc)
  CALL self%psi%new(psi_bc)
  CALL rhs_bc%add(0.d0,1.d0,psi_sol)
  CALL psi_bc%add(0.d0,1.d0,rhs_bc)
  CALL self%zerob_bc%apply(psi_bc)
  CALL rhs_bc%add(1.d0,-1.d0,psi_bc)
  CALL psi_bc%set(0.d0)
  pm_save=oft_env%pm; oft_env%pm=.FALSE.
  CALL self%lu_solver%apply(psi_bc,rhs_bc)
  oft_env%pm=pm_save
  CALL psi_vac%add(1.d0,1.d0,psi_bc)
  CALL rhs_bc%delete()
  CALL psi_bc%delete()
  DEALLOCATE(rhs_bc,psi_bc)
END IF
!
CALL psi_sol%set(0.d0)
CALL psi_sol%add(0.d0,1.d0,psi_vac)
CALL psi_sol%add(1.d0,self%vcontrol_val,psi_vcont)
IF(self%dt>0.d0)THEN
  CALL tmp_vec%set(0.d0)
  CALL tmp_vec%add(0.d0,-1.d0/self%dt,psi_sol,1.d0/self%dt,self%psi_dt)
  CALL gs_wall_source(self,tmp_vec,psi_dt)
  IF(PRESENT(rhs_source))THEN
    CALL gs_gen_source(self,rhs_source,tmp_vec)
    CALL psi_dt%add(1.d0,1.d0,tmp_vec)
  END IF
  CALL self%zerob_bc%apply(psi_dt)
  pm_save=oft_env%pm; oft_env%pm=.FALSE.
  CALL self%lu_solver_dt%apply(tmp_vec,psi_dt)
  oft_env%pm=pm_save
  CALL psi_sol%add(1.d0,1.d0,tmp_vec)
ELSE
  IF(PRESENT(rhs_source))THEN
    CALL tmp_vec%set(0.d0)
    CALL gs_gen_source(self,rhs_source,psi_dt)
    CALL self%zerob_bc%apply(psi_dt)
    pm_save=oft_env%pm; oft_env%pm=.FALSE.
    CALL self%lu_solver%apply(tmp_vec,psi_dt)
    oft_env%pm=pm_save
    CALL psi_sol%add(1.d0,1.d0,tmp_vec)
  END IF
END IF
psimax=psi_sol%dot(psi_sol)
IF(oft_env%pm)WRITE(*,'(A,I4,1ES12.4)')oft_indent,1,psimax
IF((self%dt>0.d0).AND.oft_env%pm)THEN
  CALL oft_decrease_indent
END IF
CALL psi_vac%delete
CALL psi_vcont%delete
CALL psi_eddy%delete
IF(self%dt>0.d0)THEN
  CALL psi_dt%delete
  CALL tmp_vec%delete
  DEALLOCATE(psi_dt,tmp_vec)
END IF
DEALLOCATE(psi_vac,psi_vcont,psi_eddy)
end subroutine gs_vac_solve
!---------------------------------------------------------------------------
! FUNCTION gs_err_reason
!---------------------------------------------------------------------------
!> Compute Grad-Shafranov solution for current flux definitions
!!
!! @param[in,out] self G-S object
!! @param[out] ierr Error flag
!---------------------------------------------------------------------------
function gs_err_reason(ierr) result(err_reason)
integer(4), intent(in) :: ierr
CHARACTER(LEN=40) :: err_reason
SELECT CASE(ierr)
  CASE(-1)
    err_reason='Exceeded "maxits"'
  CASE(-2)
    err_reason='Total poloidal flux is zero'
  CASE(-3)
    err_reason='Closed flux volume lost'
  CASE(-4)
    err_reason='Axis dropped below "rmin"'
  CASE(-5)
    err_reason='Toroidal current droppped too low'
  CASE(-6)
    err_reason='Matrix solve failed for targets'
  CASE(-7)
    err_reason='Isoflux fitting failed'
  CASE(-8)
    err_reason='Wall eigenmode flux loop fitting failed'
  CASE DEFAULT
    err_reason='Unkown reason'
END SELECT
end function gs_err_reason
!---------------------------------------------------------------------------
! SUBROUTINE gs_fixed_vflux
!---------------------------------------------------------------------------
!> Compute required vacuum flux for fixed boundary equilibrium
!!
!! @param[in,out] self G-S object
!---------------------------------------------------------------------------
subroutine gs_fixed_vflux(self,pts,fluxes)
class(gs_eq), intent(inout) :: self
real(8), pointer, intent(inout) :: pts(:,:),fluxes(:)
class(oft_vector), pointer :: rhs,psi_fixed,psi_dummy
real(r8), pointer, DIMENSION(:) :: vals_tmp
real(8) :: itor_alam,itor_press,estored
integer(4) :: i,io_unit
logical :: pm_save
CLASS(oft_matrix), POINTER :: dels_free
TYPE(oft_lusolver) :: lu_solver
!---
IF(self%free)CALL oft_abort("Equilibrium is free-boundary","gs_fixed_vflux",__FILE__)
WRITE(*,'(2A)')oft_indent,"Computing fixed boundary vacuum flux"
!---
NULLIFY(dels_free)
CALL compute_bcmat(self)
CALL build_dels(dels_free,self,"free")
lu_solver%A=>dels_free
!---Compute boundary term
CALL self%psi%new(rhs)
CALL self%psi%new(psi_fixed)
CALL self%psi%new(psi_dummy)
CALL gs_source(self,self%psi,rhs,psi_fixed,psi_dummy,itor_alam,itor_press,estored)
CALL psi_fixed%set(0.d0)
CALL self%zerob_bc%apply(rhs)
CALL lu_solver%apply(psi_fixed,rhs)
!---Write out error at boundary points
NULLIFY(vals_tmp)
CALL psi_fixed%get_local(vals_tmp)
! OPEN(NEWUNIT=io_unit,FILE='fixed_vflux.dat')
ALLOCATE(pts(2,self%fe_rep%mesh%nbp),fluxes(self%fe_rep%mesh%nbp))
IF(.NOT.ASSOCIATED(self%olbp))CALL get_olbp(self%mesh,self%olbp)
DO i=1,self%fe_rep%mesh%nbp
  pts(:,i)=self%fe_rep%mesh%r(1:2,self%olbp(i))
  fluxes(i)=-vals_tmp(self%olbp(i))*self%psiscale
  ! WRITE(io_unit,*)self%fe_rep%mesh%r(1:2,self%fe_rep%mesh%lbp(i)),-vals_tmp(self%fe_rep%mesh%lbp(i))*self%psiscale
END DO
! CLOSE(io_unit)
!---
CALL lu_solver%delete()
CALL dels_free%delete()
DEALLOCATE(dels_free)
!---Destroy temporary vectors
CALL rhs%delete
CALL psi_fixed%delete
CALL psi_dummy%delete
DEALLOCATE(rhs,psi_fixed,psi_dummy)
DEALLOCATE(vals_tmp)
end subroutine gs_fixed_vflux
!---------------------------------------------------------------------------
! SUBROUTINE gs_get_cond_source
!---------------------------------------------------------------------------
!> Needs Docs
!!
!! @param[in,out] self G-S object
!! @param[in,out] a Psi field
!! @param[in,out] b Source field
!---------------------------------------------------------------------------
subroutine gs_get_cond_source(self,cond_fac)
CLASS(gs_eq), intent(inout) :: self
REAL(8), intent(inout) :: cond_fac(:)
real(8) :: curr
integer(4) :: i,j,l,k
!---
cond_fac=0.d0
DO j=1,self%ncond_regs
  IF(self%cond_regions(j)%neigs>0)THEN
    DO k=1,self%cond_regions(j)%nc_quad
      curr = SUM(self%cond_regions(j)%weights*self%cond_regions(j)%cond_vals(k,:))
      DO l=1,4
        i=self%cond_regions(j)%lc(l,k)
        cond_fac(i)=cond_fac(i) + curr/self%fe_rep%mesh%ca(i)
      END DO
    END DO
  END IF
END DO
END SUBROUTINE gs_get_cond_source
!---------------------------------------------------------------------------
! SUBROUTINE gs_source
!---------------------------------------------------------------------------
!> Needs Docs
!!
!! @param[in,out] self G-S object
!! @param[in,out] a Psi field
!! @param[in,out] b Source field
!---------------------------------------------------------------------------
subroutine gs_source(self,a,b,b2,b3,itor_alam,itor_press,estore)
class(gs_eq), intent(inout) :: self
class(oft_vector), TARGET, intent(inout) :: a
CLASS(oft_vector), intent(inout) :: b
CLASS(oft_vector), intent(inout) :: b2
CLASS(oft_vector), intent(inout) :: b3
REAL(8), INTENT(out) :: itor_alam,itor_press,estore
real(r8), pointer, dimension(:) :: atmp,btmp,b2tmp,b3tmp
real(8) :: psitmp,goptmp(3,3),det,pt(3),v,ffp(3),t1
real(8), allocatable :: rhs_loc(:,:),cond_fac(:),rop(:),vcache(:)
integer(4) :: j,m,l
integer(4), allocatable :: j_lag(:)
logical :: curved
t1=omp_get_wtime()
!---
NULLIFY(atmp,btmp,b2tmp,b3tmp)
call b%set(0.d0)
call b2%set(0.d0)
call b3%set(0.d0)
IF(.NOT.self%has_plasma)RETURN
CALL b%get_local(btmp)
CALL b2%get_local(b2tmp)
CALL b3%get_local(b3tmp)
CALL a%get_local(atmp)
!---
itor_alam=0.d0
itor_press=0.d0
estore=0.d0
!$omp parallel private(rhs_loc,j_lag,ffp,curved,goptmp,v,m,det,pt,psitmp,l,rop,vcache) &
!$omp reduction(+:itor_alam) reduction(+:itor_press) reduction(+:estore)
allocate(rhs_loc(self%fe_rep%nce,3))
allocate(rop(self%fe_rep%nce),vcache(self%fe_rep%nce))
allocate(j_lag(self%fe_rep%nce))
!$omp do schedule(static,1)
do j=1,self%fe_rep%mesh%nc
  ! IF(self%cflag(j)==3)CYCLE ! Vacuum region (no source)
  IF(self%fe_rep%mesh%reg(j)/=1)CYCLE ! Only compute in plasma region
  call self%fe_rep%ncdofs(j,j_lag)
  rhs_loc=0.d0
  DO l=1,self%fe_rep%nce
    vcache(l) = atmp(j_lag(l))
  END DO
  curved=cell_is_curved(self%fe_rep%mesh,j)
  do m=1,self%fe_rep%quad%np
    if(curved.OR.(m==1))call self%fe_rep%mesh%jacobian(j,self%fe_rep%quad%pts(:,m),goptmp,v)
    det=v*self%fe_rep%quad%wts(m)
    DO l=1,self%fe_rep%nce
      CALL oft_blag_eval(self%fe_rep,j,l,self%fe_rep%quad%pts(:,m),rop(l))
    END DO
    ffp=0.d0
    pt=self%fe_rep%mesh%log2phys(j,self%fe_rep%quad%pts(:,m))
    IF(gs_test_bounds(self,pt))THEN
      psitmp=0.d0
      !$omp simd reduction(+:psitmp)
      DO l=1,self%fe_rep%nce
        psitmp=psitmp+vcache(l)*rop(l)
      END DO
      IF(psitmp>self%plasma_bounds(1))THEN
        IF(self%mode==0)THEN
          ffp(1:2)=((self%alam**2)*self%I%f(psitmp)+self%alam*self%I%f_offset)*self%I%fp(psitmp)
          itor_alam = itor_alam + self%I%Fp(psitmp)*(self%I%f(psitmp)+self%I%f_offset)/(pt(1)+gs_epsilon)
        ELSE
          ffp(1:2)=0.5d0*self%alam*self%I%fp(psitmp)
          itor_alam = itor_alam + 0.5d0*self%I%Fp(psitmp)/(pt(1)+gs_epsilon)*v*self%fe_rep%quad%wts(m)
        END IF
        ffp((/1,3/)) = ffp((/1,3/)) + (/self%pnorm,1.d0/)*self%P%fp(psitmp)*(pt(1)**2)
        !
        estore = estore + (self%P%F(psitmp))*v*self%fe_rep%quad%wts(m)*pt(1)
        itor_press = itor_press + pt(1)*self%P%Fp(psitmp)*v*self%fe_rep%quad%wts(m)
      END IF
    END IF
    pt(1) = MAX(pt(1),gs_epsilon)
    ffp = ffp*det/pt(1)
    !$omp simd
    do l=1,self%fe_rep%nce
      rhs_loc(l,:)=rhs_loc(l,:)+rop(l)*ffp
    end do
  end do
  !---Get local to global DOF mapping
  do l=1,self%fe_rep%nce
    m = j_lag(l)
    !$omp atomic
    btmp(m)=btmp(m)+rhs_loc(l,1)
    !$omp atomic
    b2tmp(m)=b2tmp(m)+rhs_loc(l,2)
    !$omp atomic
    b3tmp(m)=b3tmp(m)+rhs_loc(l,3)
  end do
end do
deallocate(rhs_loc,j_lag,rop,vcache)
!$omp end parallel
CALL b%restore_local(btmp,add=.TRUE.)
CALL b2%restore_local(b2tmp,add=.TRUE.)
CALL b3%restore_local(b3tmp,add=.TRUE.)
DEALLOCATE(atmp,btmp,b2tmp,b3tmp)
estore = estore*2.d0*pi*self%psiscale
itor_alam = itor_alam*self%psiscale
itor_press = itor_press*self%psiscale
self%timing(2)=self%timing(2)+(omp_get_wtime()-t1)
end subroutine gs_source
!---------------------------------------------------------------------------
! SUBROUTINE gs_get_chi
!---------------------------------------------------------------------------
!> Compute flux potential from Grad-Shafranov solution
!!
!! @param[in,out] self G-S object
!---------------------------------------------------------------------------
subroutine gs_get_chi(self)
class(gs_eq), intent(inout) :: self
class(oft_solver), POINTER :: solver
type(oft_lag_brinterp) :: psi_interp
class(oft_vector), pointer :: psihat,rhs
real(r8), pointer :: vals_tmp(:)
class(oft_matrix), POINTER :: dels_grnd
real(8) :: psitmp(1),goptmp(3,3),det,pt(3),rop,v,f
real(8), allocatable :: rhs_loc(:)
integer(4) :: j,m,l
integer(4), allocatable :: j_lag(:)
logical :: curved
!---
NULLIFY(dels_grnd,vals_tmp)
CALL build_dels(dels_grnd,self,"grnd")
!---Setup Solver
CALL create_cg_solver(solver)
solver%A=>dels_grnd
solver%bc=>self%zerogrnd_bc
solver%its=-2
CALL create_diag_pre(solver%pre)
!---
call self%psi%new(psihat)
call psihat%add(0.d0,1.d0/self%psimax,self%psi)
!---
call self%psi%new(rhs)
call rhs%set(0.d0)
CALL rhs%get_local(vals_tmp)
!---
psi_interp%u=>psihat
CALL psi_interp%setup(self%fe_rep)
!---
!$omp parallel private(rhs_loc,j_lag,f,curved,goptmp,v,m,det,pt,psitmp,l,rop)
allocate(rhs_loc(self%fe_rep%nce))
allocate(j_lag(self%fe_rep%nce))
!$omp do
do j=1,self%fe_rep%mesh%nc
  rhs_loc=0.d0
  curved=cell_is_curved(self%fe_rep%mesh,j)
  do m=1,self%fe_rep%quad%np
    if(curved.OR.(m==1))call self%fe_rep%mesh%jacobian(j,self%fe_rep%quad%pts(:,m),goptmp,v)
    det=v*self%fe_rep%quad%wts(m)
    pt=self%fe_rep%mesh%log2phys(j,self%fe_rep%quad%pts(:,m))
    call psi_interp%interp(j,self%fe_rep%quad%pts(:,m),goptmp,psitmp)
    IF(self%mode==0)THEN
      f=self%alam*self%I%f(psitmp(1))+self%I%f_offset
    ELSE
      f=SQRT(self%alam*self%I%F(psitmp(1)) + self%I%f_offset**2)
    END IF
    do l=1,self%fe_rep%nce
      call oft_blag_eval(self%fe_rep,j,l,self%fe_rep%quad%pts(:,m),rop)
      rhs_loc(l)=rhs_loc(l)+rop*f*det/(pt(1)+gs_epsilon)
    end do
  end do
  !---Get local to global DOF mapping
  call self%fe_rep%ncdofs(j,j_lag)
  do l=1,self%fe_rep%nce
    m=j_lag(l)
    !$omp atomic
    vals_tmp(m)=vals_tmp(m)+rhs_loc(l)
  end do
end do
!$omp end do
deallocate(rhs_loc,j_lag)
!$omp end parallel
CALL rhs%restore_local(vals_tmp,add=.TRUE.)
call self%zerob_bc%apply(rhs)
call solver%apply(self%chi,rhs)
!---
call rhs%delete
call psihat%delete
call dels_grnd%delete
DEALLOCATE(rhs,psihat,dels_grnd)
end subroutine gs_get_chi
!---------------------------------------------------------------------------
! SUBROUTINE gs_delete
!---------------------------------------------------------------------------
!> Cleanup Grad-Shafranov object internal storage
!!
!! @param[in,out] self G-S object
!---------------------------------------------------------------------------
subroutine gs_delete(self)
CLASS(gs_eq), INTENT(inout) :: self
CALL self%psi%delete
NULLIFY(self%I,self%P)
end subroutine gs_delete
!---------------------------------------------------------------------------
! FUNCTION gs_itor
!---------------------------------------------------------------------------
!> Compute toroidal current for Grad-Shafranov equilibrium
!!
!! @param[in,out] self G-S object
!! @result itor Toroidal current
!---------------------------------------------------------------------------
function gs_itor(self,psi_vec) result(itor)
class(gs_eq), intent(inout) :: self
class(oft_vector), optional, intent(inout) :: psi_vec
real(8):: itor
class(oft_vector), pointer :: x
real(r8), pointer, dimension(:) :: vals_tmp
call self%psi%new(x)
!---Compute total current
IF(PRESENT(psi_vec))THEN
  call self%dels%apply(psi_vec,x)
ELSE
  call self%dels%apply(self%psi,x)
END IF
call self%gs_zerob_bc%apply(x)
NULLIFY(vals_tmp)
call x%get_local(vals_tmp)
itor=sum(vals_tmp)*self%psiscale
call x%delete()
DEALLOCATE(x,vals_tmp)
end function gs_itor
!---------------------------------------------------------------------------
! SUBROUTINE gs_itor_nl
!---------------------------------------------------------------------------
!> Compute toroidal current for Grad-Shafranov equilibrium
!!
!! @param[in,out] self G-S object
!! @param[out] itor Toroidal current
!! @param[out] centroid Current centroid (optional) [2]
!---------------------------------------------------------------------------
subroutine gs_itor_nl(self,itor,centroid)
class(gs_eq), intent(inout) :: self
real(8), intent(out) :: itor
real(8), optional, intent(out) :: centroid(2)
type(oft_lag_brinterp), target :: psi_eval
real(8) :: itor_loc,goptmp(3,3),v,psitmp(1)
real(8) :: pt(3),curr_cent(2)
integer(4) :: i,m
IF(.NOT.self%has_plasma)THEN
  IF(PRESENT(centroid))centroid = 0.d0
  itor=0.d0
  RETURN
END IF
!---
psi_eval%u=>self%psi
CALL psi_eval%setup(self%fe_rep)
!---
itor=0.d0
curr_cent=0.d0
do i=1,self%fe_rep%mesh%nc
  IF(self%fe_rep%mesh%reg(i)/=1)CYCLE
  do m=1,self%fe_rep%quad%np
    call self%fe_rep%mesh%jacobian(i,self%fe_rep%quad%pts(:,m),goptmp,v)
    call psi_eval%interp(i,self%fe_rep%quad%pts(:,m),goptmp,psitmp)
    ! IF(psitmp(1)<self%plasma_bounds(1).OR.psitmp(1)>self%plasma_bounds(2))CYCLE
    pt=self%fe_rep%mesh%log2phys(i,self%fe_rep%quad%pts(:,m))
    !---Compute Magnetic Field
    IF(gs_test_bounds(self,pt).AND.psitmp(1)>self%plasma_bounds(1))THEN
      IF(self%mode==0)THEN
        itor_loc = (self%pnorm*pt(1)*self%P%Fp(psitmp(1)) &
        + (self%alam**2)*self%I%Fp(psitmp(1))*(self%I%f(psitmp(1))+self%I%f_offset/self%alam)/(pt(1)+gs_epsilon))
      ELSE
        itor_loc = (self%pnorm*pt(1)*self%P%Fp(psitmp(1)) &
        + .5d0*self%alam*self%I%Fp(psitmp(1))/(pt(1)+gs_epsilon))
      END IF
      itor = itor + itor_loc*v*self%fe_rep%quad%wts(m)
      curr_cent = curr_cent + itor_loc*pt(1:2)*v*self%fe_rep%quad%wts(m)
    END IF
  end do
end do
IF(PRESENT(centroid))centroid = curr_cent/itor
itor=itor*self%psiscale
end subroutine gs_itor_nl
!---------------------------------------------------------------------------
! SUBROUTINE gs_update_bounds
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gs_update_bounds(self,track_opoint)
class(gs_eq), intent(inout) :: self
logical, optional, intent(in) :: track_opoint
type(oft_lag_brinterp) :: psi_interp
integer(4) :: i,j,cell,ierr,ind_sort(max_xpoints),trace_err,itmp
REAL(8) :: old_bounds(2),f(3),goptmp(3,3),v,psitmp(1),alt_max
REAL(8) :: alt_r,alt_z,zmin,zmax,pttmp(2),vtmp
REAL(8) :: x_point(2,max_xpoints),t1,t2,x_psi(max_xpoints),x_psi_sort(max_xpoints),x_tmp
REAL(8) :: x_point_old(2,max_xpoints),x_vec_old(2,max_xpoints)
REAL(8), POINTER :: psi_vals(:)
REAL(8), PARAMETER :: tol=1.d-8, bounds_exp=2.5d-2
logical, allocatable :: x_masked(:)
t1=omp_get_wtime()
trace_err=0
old_bounds=self%plasma_bounds
IF(.NOT.ASSOCIATED(self%rlcfs))ALLOCATE(self%rlcfs(self%nlcfs,3))
!---Check boundary as limiter
NULLIFY(psi_vals)
CALL self%psi%get_local(psi_vals)
! CALL vector_cast(psiv,self%psi)

! t1=omp_get_wtime()
IF(PRESENT(track_opoint))THEN
  IF(.NOT.track_opoint)self%o_point(1)=-1.d0
ELSE
  self%o_point(1)=-1.d0
END IF
CALL gs_analyze_saddles(self, self%o_point, self%plasma_bounds(2), x_point, x_psi)
! t2=omp_get_wtime()
! WRITE(*,*)'Analyze',t2-t1
! alt_r=-1.d0
! CALL gs_psimax(self,alt_max,alt_r,alt_z)
! WRITE(*,*)'CHK:'
! WRITE(*,*)'  ',self%plasma_bounds(2),alt_max
! WRITE(*,*)'  ',self%o_point(1),alt_r
! WRITE(*,*)'  ',self%o_point(2),alt_z
IF(self%o_point(1)<0.d0)THEN
  self%plasma_bounds(1)=MAXVAL(psi_vals) !MINVAL(psiv%v)
  self%plasma_bounds(2)=MAXVAL(psi_vals)
  self%nx_points=0
  self%o_point=[0.5d0, 0.d0]
  DEALLOCATE(psi_vals)
  RETURN
END IF

self%plasma_bounds(1)=-1.d99
! CALL reset_lcfs

!---Look for X-points
self%nx_points=0
self%diverted=.FALSE.
x_psi_sort=1.d99
zmin=self%spatial_bounds(1,2)
zmax=self%spatial_bounds(2,2)
DO i=1,max_xpoints
  ! IF(ABS(x_point(2,i))<0.2d0)CYCLE
  IF(x_point(1,i)>0.d0.AND.SQRT(SUM((x_point(:,i)-self%o_point)**2))>5.d-2)THEN
    IF(x_psi(i)>self%plasma_bounds(1))THEN
      self%plasma_bounds(1)=x_psi(i)
      self%diverted=.TRUE.
    END IF
    IF(self%o_point(2)>x_point(2,i))zmin=MAX(zmin,x_point(2,i))
    IF(self%o_point(2)<x_point(2,i))zmax=MIN(zmax,x_point(2,i))
    ! IF(oft_debug_print(1))THEN
    !   WRITE(*,'(2A,5ES11.3)')oft_indent,'  X-point:',x_psi(i),x_point(:,i),self%o_point
    ! END IF
    self%nx_points=self%nx_points+1
    ! self%x_points(:,self%nx_points)=x_point(:,i)
    ind_sort(self%nx_points)=i
    x_psi_sort(self%nx_points)=x_psi(i)
  END IF
END DO

self%x_points(1,:)=-1.d0
IF(self%nx_points>0)THEN
  CALL sort_array(x_psi_sort,ind_sort,self%nx_points)
  DO i=1,self%nx_points
    self%x_points(:,i)=x_point(:,ind_sort(i))
    self%x_vecs(:,i)=self%o_point-self%x_points(:,i)
    IF(oft_debug_print(1))THEN
      WRITE(*,'(2A,5ES11.3)')oft_indent,'  X-point:',x_psi_sort(i),self%x_points(:,i),self%x_vecs(:,i)
    END IF
  END DO
  ALLOCATE(x_masked(self%nx_points))
  x_masked=.TRUE.
  DO i=1,self%nx_points
    x_vec_old(:,i)=self%x_vecs(:,i)
    x_point_old(:,i)=self%x_points(:,i)
    DO j=1,self%nx_points
      IF(i==j)CYCLE
      x_masked(i)=x_masked(i).AND.(DOT_PRODUCT(self%x_points(:,i)-self%x_points(:,j),self%x_vecs(:,j))>0.d0)
    END DO
    IF(oft_debug_print(1))THEN
      WRITE(*,'(2A,5ES11.3,L)')oft_indent,'  X-point:',x_psi_sort(i),self%x_points(:,i),self%x_vecs(:,i),.NOT.x_masked(i)
    END IF
  END DO
  j=0
  DO i=1,self%nx_points
    IF(x_masked(i))CYCLE
    j=j+1
    self%x_points(:,j)=x_point_old(:,i)
    self%x_vecs(:,j)=x_vec_old(:,i)
  END DO
  DO i=1,self%nx_points
    IF(.NOT.x_masked(i))CYCLE
    j=j+1
    self%x_points(:,j)=x_point_old(:,i)
    self%x_vecs(:,j)=x_vec_old(:,i)
    self%plasma_bounds(1)=x_psi_sort(i)
  END DO
  self%nx_points=j
  DEALLOCATE(x_masked)
END IF

psi_interp%u=>self%psi
CALL psi_interp%setup(self%fe_rep)

!---Get plasma boundary contour
IF(self%plasma_bounds(1)>-1.d98)THEN
  ! CALL get_lcfs
  self%rlcfs(1,1)=zmin
  self%rlcfs(self%nlcfs,1)=zmax
END IF
! t1=omp_get_wtime()
! WRITE(*,*)'LCFS',t1-t2

! self%plasma_bounds(1)=-1.d99
v=self%plasma_bounds(1)
! WRITE(*,*)'CHK',v
!$omp parallel private(i,j,cell,f,psitmp,vtmp,pttmp)
! reduction(max:v)
vtmp=self%plasma_bounds(1)
!$omp do
DO i=1,self%nlimiter_nds
  j=self%limiter_nds(i)
  IF(psi_vals(j)>vtmp.AND.gs_test_bounds(self,self%rlimiter_nds(:,i)))THEN
    vtmp=psi_vals(j)
    pttmp=self%rlimiter_nds(:,i)
    ! WRITE(*,*)i,self%rlimiter_nds(:,i),vtmp
  END IF
END DO
!$omp enddo nowait

!---Check limiters
cell=0
! psi_interp%u=>self%psi
!$omp do
DO i=1,self%nlimiter_pts
  IF(.NOT.gs_test_bounds(self,self%limiter_pts(:,i)))CYCLE
  CALL bmesh_findcell(self%fe_rep%mesh,cell,self%limiter_pts(:,i),f)
  IF((MAXVAL(f)>1.d0+tol).OR.(MINVAL(f)<-tol))CYCLE
  CALL psi_interp%interp(cell,f,goptmp,psitmp)
  IF(psitmp(1)>vtmp)THEN
    vtmp=psitmp(1)
    pttmp=self%limiter_pts(:,i)
  END IF
  ! v=MAX(v,psitmp(1))
END DO

!$omp critical
IF(vtmp>v)THEN
  v=vtmp
  self%lim_point=pttmp
END IF
!$omp end critical
!$omp end parallel
DEALLOCATE(psi_vals)
! t2=omp_get_wtime()
! WRITE(*,*)'Limiter',t2-t1

IF(v>self%plasma_bounds(1).OR.self%nx_points==0)THEN
  self%plasma_bounds(1)=v
  self%diverted=.FALSE.
  ! self%nx_points=0
  ! CALL get_lcfs
  IF(self%plasma_bounds(1)<-1.d98)self%plasma_bounds(1)=v
ELSE
  self%lim_point=self%x_points(:,self%nx_points)
  IF(oft_debug_print(1))THEN
    WRITE(*,'(2A,5ES11.3)')oft_indent,'  Active X-point:',self%plasma_bounds(1), &
      self%x_points(:,self%nx_points),self%o_point
  END IF
END IF
CALL psi_interp%delete()
IF(trace_err/=0)CALL oft_warn("gs_update_bounds: Trace did not complete")
! t1=omp_get_wtime()
! WRITE(*,*)'Finalize',t1-t2
IF(self%full_domain)THEN
  self%plasma_bounds=[-1.d99,1.d99]
  self%diverted=.FALSE.
  self%nx_points=0
END IF
!
IF(.NOT.self%has_plasma)RETURN
IF(oft_debug_print(1).AND.oft_env%pm)THEN
  WRITE(*,'(2A,4ES11.3)')oft_indent,'New plasma bounds',self%plasma_bounds,self%o_point
END IF
self%timing(4)=self%timing(4)+(omp_get_wtime()-t1)
end subroutine gs_update_bounds
!---------------------------------------------------------------------------
! SUBROUTINE gs_analyze_saddles
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gs_analyze_saddles(self, o_point, o_psi, x_point, x_psi)
class(gs_eq), intent(inout) :: self
real(8), intent(inout) :: o_point(2)
real(8), intent(out) :: o_psi,x_point(2,max_xpoints),x_psi(max_xpoints)
integer(4), PARAMETER :: npts = 10, max_unique = 20
integer(4) :: i,j,m,n_unique,stype,stypes(max_unique),cell,nx_points
integer(4), allocatable :: ncuts(:)
real(8) :: saddle_loc(2),saddle_psi,unique_saddles(3,max_unique),ptmp(2),f(3),loc_vals(3),psi_scale_len
real(8) :: region(2,2) = RESHAPE([-1.d99,1.d99,-1.d99,1.d99], [2,2])
type(oft_lag_brinterp), target :: psi_eval
type(oft_lag_bginterp), target :: psi_geval
type(oft_lag_bg2interp), target :: psi_g2eval
CLASS(oft_bmesh), POINTER :: smesh
!
psi_eval%u=>self%psi
psi_eval_active=>psi_eval
CALL psi_eval%setup(self%fe_rep)
CALL psi_geval%shared_setup(psi_eval)
CALL psi_g2eval%shared_setup(psi_eval)
psi_geval_active=>psi_geval
psi_g2eval_active=>psi_g2eval
!
smesh=>self%fe_rep%mesh
ALLOCATE(ncuts(smesh%np))
ncuts=0
!$omp parallel do simd private(loc_vals)
DO i=1,smesh%nc
  ! IF(smesh%reg(i)/=1)CYCLE
  IF(self%saddle_cmask(i))CYCLE
  loc_vals=psi_eval%vals(smesh%lc(:,i))
  IF((loc_vals(1)-loc_vals(2))*(loc_vals(3)-loc_vals(1))>0)THEN
    !$omp atomic
    ncuts(smesh%lc(1,i))=ncuts(smesh%lc(1,i))+1
  END IF
  IF((loc_vals(2)-loc_vals(3))*(loc_vals(1)-loc_vals(2))>0)THEN
    !$omp atomic
    ncuts(smesh%lc(2,i))=ncuts(smesh%lc(2,i))+1
  END IF
  IF((loc_vals(3)-loc_vals(1))*(loc_vals(2)-loc_vals(3))>0)THEN
    !$omp atomic
    ncuts(smesh%lc(3,i))=ncuts(smesh%lc(3,i))+1
  END IF
END DO
!$omp do simd
DO i=1,smesh%np
  IF(self%saddle_pmask(i))ncuts(i)=-1
END DO
!
psi_scale_len = ABS(self%plasma_bounds(2)-self%plasma_bounds(1))*5.d0/(SQRT(self%lim_area))
unique_saddles=-1.d99
o_psi=-1.d99
n_unique=0
!
DO i=1,smesh%np
  IF((ncuts(i)==0).OR.(ncuts(i)>3))THEN
    IF((ncuts(i)==0).AND.(o_point(1)>0.d0))THEN
      saddle_loc=o_point
    ELSE
      saddle_loc=smesh%r(1:2,i)
    END IF
    ! IF(ALL(smesh%reg(smesh%lpc(smesh%kpc(i):smesh%kpc(i+1)-1))/=1))CYCLE
    IF(self%fe_rep%order>1)THEN
      CALL gs_find_saddle(self,psi_scale_len,saddle_psi,saddle_loc,stype)
    ELSE
      saddle_psi=psi_eval%vals(i)
      IF(ncuts(i)==0)stype=1
      IF(ncuts(i)>3)stype=3
    END IF
    IF(stype<0)CYCLE
    IF(saddle_psi>-1.d98)THEN
      DO m=1,n_unique
        IF(SQRT(SUM((saddle_loc-unique_saddles(1:2,m))**2))<2.d-2)EXIT
      END DO
      IF(m>n_unique.AND.n_unique<max_unique)THEN
        n_unique = n_unique + 1
        unique_saddles(1:2,n_unique) = saddle_loc
        unique_saddles(3,n_unique) = saddle_psi
        stypes(n_unique) = stype
        !
        cell=0
        CALL bmesh_findcell(smesh,cell,saddle_loc,f)
        IF(smesh%reg(cell)/=1)CYCLE
        o_psi = MAX(o_psi,saddle_psi)
      END IF
    END IF
  END IF
END DO
DEALLOCATE(ncuts)
CALL psi_eval%delete
CALL psi_geval%delete
CALL psi_g2eval%delete
IF(oft_debug_print(2))THEN
  WRITE(*,*)
  WRITE(*,*)'Saddle points',n_unique
END IF
o_point(1)=-1.d0
x_point(1,:)=-1.d0
nx_points=0
DO m=1,n_unique
  ptmp=unique_saddles(1:2,m)
  cell=0
  CALL bmesh_findcell(smesh,cell,unique_saddles(1:2,m),f)
  IF(oft_debug_print(2))WRITE(*,*)stypes(m),unique_saddles(:,m),smesh%reg(cell)
  IF(self%saddle_rmask(smesh%reg(cell)))CYCLE
  IF(ABS(o_psi-unique_saddles(3,m))<1.d-8)THEN
    IF(smesh%reg(cell)/=1)CYCLE
    o_point=unique_saddles(1:2,m)
    CYCLE
  END IF
  !
  nx_points=nx_points+1
  IF(nx_points>max_xpoints)CALL oft_abort("Found too many X-points","gs_analyze_saddles",__FILE__)
  x_point(:,nx_points)=unique_saddles(1:2,m)
  x_psi(nx_points)=unique_saddles(3,m)
END DO
IF(oft_debug_print(2))WRITE(*,*)
end subroutine gs_analyze_saddles
!---------------------------------------------------------------------------
! SUBROUTINE gs_find_saddle
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gs_find_saddle(self,psi_scale_len,psi_x,pt,stype)
class(gs_eq), intent(inout) :: self
real(8), intent(in) :: psi_scale_len
real(8), intent(inout) :: psi_x
real(8), intent(inout) :: pt(2)
integer(4), intent(out) :: stype
integer(4) :: i
real(8) :: ptmp(2),gpsitmp(2),f(3),goptmp(3,3),v,mag_min,psi_circ(36),pt_circ(2)
!---MINPACK variables
real(8) :: ftol,xtol,gtol,epsfcn,factor,dx
real(8), allocatable, dimension(:) :: diag,wa1,wa2,wa3,wa4,qtf
real(8), allocatable, dimension(:,:) :: fjac
integer(4) :: maxfev,mode,nprint,info,nfev,ldfjac,ncons,ncofs,njev
integer(4), allocatable, dimension(:) :: ipvt
!---Determine initial guess from cell search
! psi_eval%u=>self%psi
! psi_geval%u=>self%psi
stype=-1
f=1.d0/3.d0
mag_min=1.d99
psi_x=-1.d99
goptmp=1.d0
cell_active=0
CALL bmesh_findcell(self%fe_rep%mesh,cell_active,pt,f)
IF((cell_active==0).OR.(minval(f)<-1.d-3).OR.(maxval(f)>1.d0+1.d-3))RETURN
!---Use MINPACK to find maximum (zero gradient)
ncons=2
ncofs=2
allocate(diag(ncofs),fjac(ncons,ncofs))
allocate(qtf(ncofs),wa1(ncofs),wa2(ncofs))
allocate(wa3(ncofs),wa4(ncons))
allocate(ipvt(ncofs))
mode = 1
factor = 1.d0
maxfev = 100
ftol = 1.d-9
xtol = 1.d-8
gtol = 1.d-8
epsfcn = SQRT(self%fe_rep%mesh%ca(cell_active)*2.d0)/REAL(self%fe_rep%order,8)*0.04d0 !5.d-4
nprint = 0
ldfjac = ncons
ptmp=pt
call lmdif(psimax_error,ncons,ncofs,ptmp,gpsitmp, &
            ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
            nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
! call lmder(psimax_error_grad,ncons,ncofs,ptmp,gpsitmp,fjac,ldfjac, &
!             ftol,xtol,gtol,maxfev,diag,mode,factor,nprint,info,nfev,njev, &
!             ipvt,qtf,wa1,wa2,wa3,wa4)
deallocate(diag,fjac,qtf,wa1,wa2)
deallocate(wa3,wa4,ipvt)
!---Get axis values
CALL bmesh_findcell(self%fe_rep%mesh,cell_active,ptmp,f)
IF((cell_active==0).OR.(minval(f)<-1.d-3).OR.(maxval(f)>1.d0+1.d-3))THEN
  ! CALL psi_eval%delete()
  RETURN
END IF
IF(self%saddle_rmask(self%fe_rep%mesh%reg(cell_active)))RETURN ! Dont allow saddles outside of allowable regions
IF(SQRT(SUM(gpsitmp**2))>psi_scale_len)RETURN
call psi_eval_active%interp(cell_active,f,goptmp,gpsitmp(1:1))
psi_x=gpsitmp(1)
pt=ptmp
stype=1
end subroutine gs_find_saddle
!---------------------------------------------------------------------------
! FUNCTION gs_test_bounds
!---------------------------------------------------------------------------
!> Test whether a point is inside the LCFS
!!
!! @param[in,out] self G-S object
!! @param[in] pt Location to test in/out of plasma
!---------------------------------------------------------------------------
function gs_test_bounds(self,pt) result(in_bounds)
class(gs_eq), intent(inout) :: self
real(8), intent(in) :: pt(2)
integer(4) :: i
real(8) :: rmin,rmax
logical :: in_bounds
in_bounds=.TRUE.
in_bounds=in_bounds.AND.(pt(1)>=self%spatial_bounds(1,1).AND.pt(1)<=self%spatial_bounds(2,1))
in_bounds=in_bounds.AND.(pt(2)>=self%spatial_bounds(1,2).AND.pt(2)<=self%spatial_bounds(2,2))
DO i=1,self%nx_points
  in_bounds=in_bounds.AND.(DOT_PRODUCT(pt-self%x_points(:,i),self%x_vecs(:,i))>0.d0)
END DO
end function gs_test_bounds
!---------------------------------------------------------------------------
! SUBROUTINE gs_save_fields
!---------------------------------------------------------------------------
!> Compute magnetic fields from Grad-Shafranov equilibrium
!!
!! @param[in,out] self G-S object
!! @param[in] pts Sampling locations [2,npts]
!! @param[in] npts Number of points to sample
!! @param[in] filename Output file for field data
!---------------------------------------------------------------------------
subroutine gs_save_fields(self,pts,npts,filename)
class(gs_eq), intent(inout) :: self
real(8), intent(in) :: pts(2,npts)
integer(4), intent(in) :: npts
character(LEN=*), intent(in) :: filename
type(oft_lag_brinterp), target :: psi_eval
type(oft_lag_bginterp), target :: psi_geval
real(8) :: v,psitmp(1),gpsitmp(3),f(3),goptmp(3,3),B(5),pttmp(3)
integer(4) :: i,cell,io_unit
!---
psi_eval%u=>self%psi
CALL psi_eval%setup(self%fe_rep)
CALL psi_geval%shared_setup(psi_eval)
!---
cell=0
OPEN(NEWUNIT=io_unit,FILE=TRIM(filename))
WRITE(io_unit,'(A)')"# R, Z, Br, Bt, Bz, Psi-Psi_a, P"
!---
DO i=1,npts
  cell=0
  pttmp(1:2)=pts(:,i)
  CALL bmesh_findcell(self%fe_rep%mesh,cell,pttmp,f)
  CALL self%fe_rep%mesh%jacobian(cell,f,goptmp,v)
  CALL psi_eval%interp(cell,f,goptmp,psitmp)
  CALL psi_geval%interp(cell,f,goptmp,gpsitmp)
  !---
  B((/1,3/))=(/-gpsitmp(2),gpsitmp(1)/)/pts(1,i)
  IF(self%mode==0)THEN
    B(2)=self%alam*(self%I%f(psitmp(1))+self%I%f_offset/self%alam)/pts(1,i)
  ELSE
    B(2)=SQRT(self%alam*self%I%f(psitmp(1)) + self%I%f_offset**2)/pts(1,i)
  END IF
  B=B*self%psiscale
  !
  B(5)=self%pnorm*self%P%f(psitmp(1))*self%psiscale/mu0
  IF(self%plasma_bounds(1)<-1.d98)THEN
    B(4)=psitmp(1)
  ELSE
    B(4)=psitmp(1)-self%plasma_bounds(1)
  END IF
  !---
  WRITE(io_unit,'(7E18.10)')pts(:,i),B
END DO
CLOSE(io_unit)
end subroutine gs_save_fields
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE psi2pt_error(m,n,cofs,err,iflag)
integer(4), intent(in) :: m,n
real(8), intent(in) :: cofs(n)
real(8), intent(out) :: err(m)
integer(4), intent(inout) :: iflag
real(8) :: f(3),goptmp(3,3),psitmp(1),pt(2)
real(8), parameter :: tol=1.d-10
!---
pt=cofs(1)*vec_con_active + pt_con_active
call bmesh_findcell(psi_eval_active%mesh,cell_active,pt,f)
IF(cell_active==0)THEN
  err(1)=psi_target_active
  RETURN
END IF
call psi_eval_active%interp(cell_active,f,goptmp,psitmp)
err(1)=psitmp(1)-psi_target_active
end subroutine psi2pt_error
!---------------------------------------------------------------------------
!> Find position of psi along a vector search direction
!---------------------------------------------------------------------------
subroutine gs_psi2pt(self,psi_target,pt,pt_con,vec,psi_int)
class(gs_eq), intent(inout) :: self
real(8), intent(in) :: psi_target !< Target \f$ \psi \f$ value to find
real(8), intent(inout) :: pt(2) !< Guess location (input); Closest point found (output) [2]
real(8), intent(in) :: pt_con(2) !< Location defining origin of search path
real(8), intent(in) :: vec(2) !< Vector defining direction of search path
type(oft_lag_brinterp), target, optional, intent(inout) :: psi_int !< Interpolation object (created internally if not passed)
type(oft_lag_brinterp), target :: psi_eval
!---MINPACK variables
real(8) :: ftol,xtol,gtol,epsfcn,factor,cofs(1),error(1)
real(8), allocatable, dimension(:) :: diag,wa1,wa2,wa3,wa4,qtf
real(8), allocatable, dimension(:,:) :: fjac
integer(4) :: maxfev,mode,nprint,info,nfev,ldfjac,ncons,ncofs
integer(4), allocatable, dimension(:) :: ipvt
!---
IF(PRESENT(psi_int))THEN
  psi_eval_active=>psi_int
ELSE
  psi_eval%u=>self%psi
  CALL psi_eval%setup(self%fe_rep)
  psi_eval_active=>psi_eval
END IF
psi_target_active=psi_target
cell_active=0
! z_con_active=z
! rax_con_active=self%o_point(1)
pt_con_active=pt_con
vec_con_active=vec
!---Use MINPACK to find maximum (zero gradient)
ncons=1
ncofs=1
allocate(diag(ncofs),fjac(ncons,ncofs))
allocate(qtf(ncofs),wa1(ncofs),wa2(ncofs))
allocate(wa3(ncofs),wa4(ncons))
allocate(ipvt(ncofs))
mode = 1
factor = 1.d0
maxfev = 100
ftol = 1.d-9
xtol = 1.d-8
gtol = 1.d-8
epsfcn = 1.d-4
nprint = 0
ldfjac = ncons
cofs(1)=DOT_PRODUCT(pt-pt_con_active,vec_con_active)
call lmdif(psi2pt_error,ncons,ncofs,cofs,error, &
              ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
              nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
deallocate(diag,fjac,qtf,wa1,wa2)
deallocate(wa3,wa4,ipvt)
IF(.NOT.PRESENT(psi_int))CALL psi_eval%delete()
!---Save back result
pt=pt_con_active+cofs(1)*vec_con_active
end subroutine gs_psi2pt
!---------------------------------------------------------------------------
! SUBROUTINE gs_psi2r
!---------------------------------------------------------------------------
!> Find position of psi along a radial chord
!!
!! @param[in,out] self G-S object
!! @param[in] psi_target
!! @param[in,out] r
!! @param[in] z
!---------------------------------------------------------------------------
subroutine gs_psi2r(self,psi_target,pt,psi_int)
class(gs_eq), intent(inout) :: self
real(8), intent(in) :: psi_target
real(8), intent(inout) :: pt(2)
type(oft_lag_brinterp), target, optional, intent(inout) :: psi_int
real(8) :: vec(2)
vec=[1.d0,0.d0]
CALL gs_psi2pt(self,psi_target,pt,[self%o_point(1),pt(2)],vec,psi_int)
end subroutine gs_psi2r
!---------------------------------------------------------------------------
! FUNCTION gs_beta
!---------------------------------------------------------------------------
!> Compute plasma Beta
!!
!! @param[in] beta_mr Minor radius for optional calculations
!! @result Beta by several different metrics
!---------------------------------------------------------------------------
function gs_beta(self,beta_mr) result(beta)
class(gs_eq), intent(inout) :: self
real(8), optional, intent(in) :: beta_mr
type(oft_lag_brinterp), target :: psi_eval
type(oft_lag_bginterp), target :: gpsi_eval
real(8) :: goptmp(3,3),v,gpsitmp(3),psitmp(1),pt(3),vol,beta(3)
real(8) :: Pvol,Psq,Bsq,Bpsq,Itor,B(3),Bmax
integer(4) :: i,m
!---
CALL gs_itor_nl(self,Itor)
Itor=Itor/self%psiscale
!---
psi_eval%u=>self%psi
CALL psi_eval%setup(self%fe_rep)
CALL gpsi_eval%shared_setup(psi_eval)
!---
Psq=0.d0
Pvol=0.d0
Bsq=0.d0
Bpsq=0.d0
vol=0.d0
Bmax=0.d0
do i=1,self%fe_rep%mesh%nc
  IF(self%fe_rep%mesh%reg(i)/=1)CYCLE
  do m=1,self%fe_rep%quad%np
    call self%fe_rep%mesh%jacobian(i,self%fe_rep%quad%pts(:,m),goptmp,v)
    call psi_eval%interp(i,self%fe_rep%quad%pts(:,m),goptmp,psitmp)
    ! IF(psitmp(1)<self%plasma_bounds(1).OR.psitmp(1)>self%plasma_bounds(2))CYCLE
    pt=self%fe_rep%mesh%log2phys(i,self%fe_rep%quad%pts(:,m))
    IF(gs_test_bounds(self,pt).AND.(psitmp(1)>self%plasma_bounds(1)))THEN
      call gpsi_eval%interp(i,self%fe_rep%quad%pts(:,m),goptmp,gpsitmp)
      !---Compute Magnetic Field
      B(1) = -gpsitmp(2)/(pt(1)+gs_epsilon)
      B(2) = gpsitmp(1)/(pt(1)+gs_epsilon)
      IF(self%mode==0)THEN
        B(3) = self%alam*(self%I%F(psitmp(1))+self%I%f_offset/self%alam)/(pt(1)+gs_epsilon)
      ELSE
        B(3) = (SQRT(self%alam*self%I%F(psitmp(1)) + self%I%f_offset**2) &
        + self%I%f_offset*(1.d0-SIGN(1.d0,self%I%f_offset)))/(pt(1)+gs_epsilon)
      END IF
      !---Update integrand
      Pvol = Pvol + (self%pnorm*self%P%F(psitmp(1)))*v*self%fe_rep%quad%wts(m)*pt(1)
      Psq = Psq + ((self%pnorm*self%P%F(psitmp(1)))**2)*v*self%fe_rep%quad%wts(m)*pt(1)
      Bsq = Bsq + DOT_PRODUCT(B,B)*v*self%fe_rep%quad%wts(m)*pt(1)
      Bpsq = Bpsq + DOT_PRODUCT(B(1:2),B(1:2))*v*self%fe_rep%quad%wts(m)*pt(1)
      vol = vol + v*self%fe_rep%quad%wts(m)*pt(1)
      !---Update max beta
      Bmax = MAX(Bmax,2.d0*self%pnorm*self%P%F(psitmp(1))/DOT_PRODUCT(B,B))
    END IF
  end do
end do
!---
beta(1)=2.d0*Pvol/Bsq
beta(2)=Bmax
beta(3)=2.d0*Pvol/Bpsq
!---
IF(PRESENT(beta_mr))THEN
  beta(1)=2.d0*SQRT(Psq/vol)/((Itor/(2.d0*pi*beta_mr))**2)
  beta(2)=2.d0*(Pvol/vol)/((Itor/(2.d0*pi*beta_mr))**2)
  IF(self%I%f_offset/=0.d0)THEN
    B(3)=self%I%f_offset/self%o_point(1)
    beta(3)=2.d0*(Pvol/vol)/(B(3)**2)
  ELSE
    B(3)=0.d0
  END IF
END IF
end function gs_beta
!---------------------------------------------------------------------------
! FUNCTION gs_estored
!---------------------------------------------------------------------------
!> Compute plasma stored energy
!!
!! @result Plasma stored energy \f$ \int P dV \f$ [J]
!---------------------------------------------------------------------------
function gs_estored(self) result(wstored)
class(gs_eq), intent(inout) :: self
type(oft_lag_brinterp), target :: psi_eval
real(8) :: goptmp(3,3),v,psitmp(1),pt(3),wstored
integer(4) :: i,m
logical, parameter :: curved = .FALSE.
!---
psi_eval%u=>self%psi
CALL psi_eval%setup(self%fe_rep)
wstored=0.d0
!$omp parallel do private(m,goptmp,v,psitmp,pt) reduction(+:wstored)
do i=1,self%fe_rep%mesh%nc
  IF(self%fe_rep%mesh%reg(i)/=1)CYCLE
  ! Fetch whether curved or not
  do m=1,self%fe_rep%quad%np
    pt=self%fe_rep%mesh%log2phys(i,self%fe_rep%quad%pts(:,m))
    IF(gs_test_bounds(self,pt))THEN
      if(curved.OR.(m==1))call self%fe_rep%mesh%jacobian(i,self%fe_rep%quad%pts(:,m),goptmp,v)
      call psi_eval%interp(i,self%fe_rep%quad%pts(:,m),goptmp,psitmp)
      IF(psitmp(1)>self%plasma_bounds(1))THEN
        !---Update integrand
        wstored = wstored + (self%pnorm*self%P%F(psitmp(1)))*v*self%fe_rep%quad%wts(m)*pt(1)
      END IF
    END IF
  end do
end do
wstored = wstored*2.d0*pi*self%psiscale
end function gs_estored
!---------------------------------------------------------------------------
! FUNCTION gs_dflux
!---------------------------------------------------------------------------
!> Compute diamagentic flux
!!
!! @result Toroidal flux increment due to plasma
!---------------------------------------------------------------------------
function gs_dflux(self) result(dflux)
class(gs_eq), intent(inout) :: self
type(oft_lag_brinterp), target :: psi_eval
real(8) :: goptmp(3,3),v,psitmp(1),pt(3)
real(8) :: Btor,dflux
integer(4) :: i,m
!---
psi_eval%u=>self%psi
CALL psi_eval%setup(self%fe_rep)
dflux=0.d0
!$omp parallel do private(m,pt,goptmp,v,psitmp,Btor) reduction(+:dflux)
do i=1,self%fe_rep%mesh%nc
  IF(self%fe_rep%mesh%reg(i)/=1)CYCLE
  do m=1,self%fe_rep%quad%np
    call self%fe_rep%mesh%jacobian(i,self%fe_rep%quad%pts(:,m),goptmp,v)
    call psi_eval%interp(i,self%fe_rep%quad%pts(:,m),goptmp,psitmp)
    ! IF(psitmp(1)<self%plasma_bounds(1).OR.psitmp(1)>self%plasma_bounds(2))CYCLE
    pt=self%fe_rep%mesh%log2phys(i,self%fe_rep%quad%pts(:,m))
    IF(gs_test_bounds(self,pt).AND.(psitmp(1)>self%plasma_bounds(1)))THEN
      pt(1)=MAX(pt(1),gs_epsilon)
      !---Compute differential toroidal Field
      IF(self%mode==0)THEN
        Btor = self%alam*(self%I%F(psitmp(1)))/pt(1)
      ELSE
        Btor = (SIGN(1.d0,self%I%f_offset)*SQRT(self%alam*self%I%F(psitmp(1)) + self%I%f_offset**2) &
        - self%I%f_offset)/pt(1)
      END IF
      !---Update integrand
      dflux = dflux + Btor*v*self%fe_rep%quad%wts(m)
    END IF
  end do
end do
dflux=dflux*self%psiscale
end function gs_dflux
!---------------------------------------------------------------------------
! FUNCTION gs_tflux
!---------------------------------------------------------------------------
!> Compute total enclosed toroidal flux
!!
!! @result Toroidal flux enclosed by LCFS
!---------------------------------------------------------------------------
function gs_tflux(self) result(tflux)
class(gs_eq), intent(inout) :: self
type(oft_lag_brinterp), target :: psi_eval
real(8) :: goptmp(3,3),v,psitmp(1),pt(3)
real(8) :: Btor,tflux
integer(4) :: i,m
!---
psi_eval%u=>self%psi
CALL psi_eval%setup(self%fe_rep)
tflux=0.d0
!$omp parallel do private(m,pt,goptmp,v,psitmp,Btor) reduction(+:tflux)
do i=1,self%fe_rep%mesh%nc
  IF(self%fe_rep%mesh%reg(i)/=1)CYCLE
  do m=1,self%fe_rep%quad%np
    call self%fe_rep%mesh%jacobian(i,self%fe_rep%quad%pts(:,m),goptmp,v)
    call psi_eval%interp(i,self%fe_rep%quad%pts(:,m),goptmp,psitmp)
    ! IF(psitmp(1)<self%plasma_bounds(1).OR.psitmp(1)>self%plasma_bounds(2))CYCLE
    pt=self%fe_rep%mesh%log2phys(i,self%fe_rep%quad%pts(:,m))
    IF(gs_test_bounds(self,pt).AND.(psitmp(1)>self%plasma_bounds(1)))THEN
      !---Compute total toroidal Field
      IF(self%mode==0)THEN
        Btor = (self%alam*(self%I%F(psitmp(1))) + self%I%f_offset)/(pt(1)+gs_epsilon)
      ELSE
        Btor = (SIGN(1.d0,self%I%f_offset)*SQRT(self%alam*self%I%F(psitmp(1)) + self%I%f_offset**2))/(pt(1)+gs_epsilon)
      END IF
      !---Update integrand
      tflux = tflux + Btor*v*self%fe_rep%quad%wts(m)
    END IF
  end do
end do
tflux=tflux*self%psiscale
end function gs_tflux
!---------------------------------------------------------------------------
! SUBROUTINE gs_li
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gs_li(self,psi_lim,li)
class(gs_eq), intent(inout) :: self
real(8), intent(in) :: psi_lim
real(8), intent(inout) :: li
type(oft_lag_brinterp), target :: psi_eval
type(oft_lag_bginterp), target :: gpsi_eval
real(8) :: goptmp(3,3),v,psitmp(1),gpsitmp(3),B(3),It,pt(3),Bmsq
integer(4) :: i,m
!---
psi_eval%u=>self%psi
CALL psi_eval%setup(self%fe_rep)
gpsi_eval%u=>self%psi
CALL gpsi_eval%setup(self%fe_rep)
!---
li=0.d0
It=0.d0
do i=1,self%fe_rep%mesh%nc
  do m=1,self%fe_rep%quad%np
    call self%fe_rep%mesh%jacobian(i,self%fe_rep%quad%pts(:,m),goptmp,v)
    call psi_eval%interp(i,self%fe_rep%quad%pts(:,m),goptmp,psitmp)
    call gpsi_eval%interp(i,self%fe_rep%quad%pts(:,m),goptmp,gpsitmp)
    pt=self%fe_rep%mesh%log2phys(i,self%fe_rep%quad%pts(:,m))
    if(psitmp(1)<psi_lim)cycle
    !---Compute Magnetic Field
    B(1) = -gpsitmp(2)/(pt(1)+gs_epsilon)
    B(2) = gpsitmp(1)/(pt(1)+gs_epsilon)
    B(3) = self%alam*self%I%F(psitmp(1))/(pt(1)+gs_epsilon)
    Bmsq = SUM(B(1:2)**2)
    !---Update integrand
    li = li + Bmsq*v*self%fe_rep%quad%wts(m)*pt(1)
  end do
end do

end subroutine gs_li
!---------------------------------------------------------------------------
! SUBROUTINE gs_epar
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gs_epar(self,psi_lim,epar)
class(gs_eq), intent(inout) :: self
real(8), intent(in) :: psi_lim
real(8), intent(inout) :: epar
type(oft_lag_brinterp), target :: psi_eval
type(oft_lag_bginterp), target :: gpsi_eval,gchi_eval
real(8) :: goptmp(3,3),v,psitmp(1),gchitmp(3),gpsitmp(3),B(3),A(3),pt(3)
integer(4) :: i,m
!---
psi_eval%u=>self%psi
CALL psi_eval%setup(self%fe_rep)
gpsi_eval%u=>self%psi
CALL gpsi_eval%setup(self%fe_rep)
gchi_eval%u=>self%chi
CALL gchi_eval%setup(self%fe_rep)
!---
epar=0.d0
do i=1,self%fe_rep%mesh%nc
  do m=1,self%fe_rep%quad%np
    call self%fe_rep%mesh%jacobian(i,self%fe_rep%quad%pts(:,m),goptmp,v)
    call psi_eval%interp(i,self%fe_rep%quad%pts(:,m),goptmp,psitmp)
    call gpsi_eval%interp(i,self%fe_rep%quad%pts(:,m),goptmp,gpsitmp)
    call gchi_eval%interp(i,self%fe_rep%quad%pts(:,m),goptmp,gchitmp)
    pt=self%fe_rep%mesh%log2phys(i,self%fe_rep%quad%pts(:,m))
    if(psitmp(1)<psi_lim)cycle
    !---Compute vector potential
    A(1) = -gchitmp(2)/(pt(1)+gs_epsilon)
    A(2) = gchitmp(1)/(pt(1)+gs_epsilon)
    A(3) = psitmp(1)/(pt(1)+gs_epsilon)
    !---Compute Magnetic Field
    B(1) = -gpsitmp(2)/(pt(1)+gs_epsilon)
    B(2) = gpsitmp(1)/(pt(1)+gs_epsilon)
    B(3) = self%alam*self%I%F(psitmp(1))/(pt(1)+gs_epsilon)
    !---Update integrand
    epar = epar + (dot_product(A,B)/SQRT(SUM(B**2)))*v*self%fe_rep%quad%wts(m)*pt(1)
  end do
end do

end subroutine gs_epar
!---------------------------------------------------------------------------
! SUBROUTINE gs_helicity
!---------------------------------------------------------------------------
!> Compute the magnetic energy and helicity of a fixed boundary equilibrium
!!
!! @note Helicity computed by this subroutine is only valid for equilibria
!! with no normal field on the boundary.
!!
!! @param[out] ener Total magnetic energy
!! @param[out] helic Total magnetic helicity
!---------------------------------------------------------------------------
subroutine gs_helicity(self,ener,helic)
class(gs_eq), intent(inout) :: self
real(8), intent(out) :: ener,helic
type(oft_lag_brinterp), target :: psi_eval
type(oft_lag_bginterp), target :: gpsi_eval,gchi_eval
real(8) :: goptmp(3,3),v,psitmp(1),gchitmp(3),gpsitmp(3),B(3),A(3),pt(3)
integer(4) :: i,m
logical :: pm_save
!---
pm_save=oft_env%pm; oft_env%pm=.FALSE.
CALL self%get_chi
oft_env%pm=pm_save
!---
psi_eval%u=>self%psi
CALL psi_eval%setup(self%fe_rep)
gpsi_eval%u=>self%psi
CALL gpsi_eval%setup(self%fe_rep)
gchi_eval%u=>self%chi
CALL gchi_eval%setup(self%fe_rep)
!---
ener=0.d0
helic=0.d0
!$omp parallel do private(m,goptmp,v,psitmp,gpsitmp,gchitmp,pt,A,B) reduction(+:ener) reduction(+:helic)
do i=1,self%fe_rep%mesh%nc
  do m=1,self%fe_rep%quad%np
    call self%fe_rep%mesh%jacobian(i,self%fe_rep%quad%pts(:,m),goptmp,v)
    call psi_eval%interp(i,self%fe_rep%quad%pts(:,m),goptmp,psitmp)
    call gpsi_eval%interp(i,self%fe_rep%quad%pts(:,m),goptmp,gpsitmp)
    call gchi_eval%interp(i,self%fe_rep%quad%pts(:,m),goptmp,gchitmp)
    pt=self%fe_rep%mesh%log2phys(i,self%fe_rep%quad%pts(:,m))
    !---Compute vector potential
    A(1) = -gchitmp(2)/(pt(1)+gs_epsilon)
    A(2) = gchitmp(1)/(pt(1)+gs_epsilon)
    A(3) = psitmp(1)/(pt(1)+gs_epsilon)
    !---Compute Magnetic Field
    B(1) = -gpsitmp(2)/(pt(1)+gs_epsilon)
    B(2) = gpsitmp(1)/(pt(1)+gs_epsilon)
    B(3) = self%alam*self%I%F(psitmp(1))/(pt(1)+gs_epsilon)
    !---Update integrand
    helic = helic + DOT_PRODUCT(A,B)*v*self%fe_rep%quad%wts(m)*pt(1)
    ener = ener + DOT_PRODUCT(B,B)*v*self%fe_rep%quad%wts(m)*pt(1)
  end do
end do
end subroutine gs_helicity
!------------------------------------------------------------------------------
! SUBROUTINE psimax_error
!------------------------------------------------------------------------------
!>
!------------------------------------------------------------------------------
SUBROUTINE psimax_error(m,n,cofs,err,iflag)
integer(4), intent(in) :: m,n
real(8), intent(in) :: cofs(n)
real(8), intent(out) :: err(m)
integer(4), intent(inout) :: iflag
real(8) :: f(3),goptmp(3,3),v,err_tmp(3)
!---
call bmesh_findcell(psi_geval_active%mesh,cell_active,cofs,f)
IF(cell_active==0)THEN
  err(1:2)=0.d0
  RETURN
END IF
call psi_geval_active%mesh%jacobian(cell_active,f,goptmp,v)
call psi_geval_active%interp(cell_active,f,goptmp,err_tmp)
err(1:2)=err_tmp(1:2)
end subroutine psimax_error
!------------------------------------------------------------------------------
! SUBROUTINE psimax_error
!------------------------------------------------------------------------------
!>
!------------------------------------------------------------------------------
SUBROUTINE psimax_error_grad(m,n,cofs,err,jac_mat,ldjac_mat,iflag)
integer(4), intent(in) :: m,n,ldjac_mat
real(8), intent(in) :: cofs(n)
real(8), intent(out) :: err(m),jac_mat(ldjac_mat,n)
integer(4), intent(in) :: iflag
real(8) :: f(3),goptmp(3,3),v,d2_tmp(6),err_tmp(3)
!---
IF(iflag==1)THEN
  call bmesh_findcell(psi_geval_active%mesh,cell_active,cofs,f)
  call psi_geval_active%mesh%jacobian(cell_active,f,goptmp,v)
  call psi_geval_active%interp(cell_active,f,goptmp,err_tmp)
  err(1:2)=err_tmp(1:2)
ELSE
  call psi_g2eval_active%interp(cell_active,f,goptmp,d2_tmp)
  jac_mat(1,1)=d2_tmp(1)
  jac_mat(2,1)=d2_tmp(2)
  jac_mat(1,2)=d2_tmp(2)
  jac_mat(2,2)=d2_tmp(4)
END IF
end subroutine psimax_error_grad
!---------------------------------------------------------------------------
! SUBROUTINE gs_psimax
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gs_psimax(self,psi_max,r,z)
class(gs_eq), intent(inout) :: self
real(8), intent(inout) :: psi_max
real(8), intent(inout) :: r
real(8), intent(inout) :: z
type(oft_lag_brinterp), target :: psi_eval
type(oft_lag_bginterp), target :: psi_geval
type(oft_lag_bg2interp), target :: psi_g2eval
integer(4) :: i
real(8) :: ptmp(3),gpsitmp(2),f(3),goptmp(3,3),v,rb(2),zb(2)
!---MINPACK variables
real(8) :: ftol,xtol,gtol,epsfcn,factor
real(8), allocatable, dimension(:) :: diag,wa1,wa2,wa3,wa4,qtf
real(8), allocatable, dimension(:,:) :: fjac
integer(4) :: maxfev,mode,nprint,info,nfev,ldfjac,ncons,ncofs,njev
integer(4), allocatable, dimension(:) :: ipvt
! rb=self%spatial_bounds(:,1)
! zb=self%spatial_bounds(:,2)
!---Determine initial guess from cell search
psi_eval%u=>self%psi
CALL psi_eval%setup(self%fe_rep)
f=1.d0/3.d0
psi_max=-1.d99
goptmp=1.d0
cell_active=0
IF(r<0.d0)THEN
  DO i=1,self%fe_rep%mesh%nc
    ptmp = self%fe_rep%mesh%log2phys(i,f)
    ! IF(ptmp(1)<rb(1).OR.ptmp(1)>rb(2))CYCLE
    ! IF(ptmp(2)<zb(1).OR.ptmp(2)>zb(2))CYCLE
    IF(self%fe_rep%mesh%reg(i)/=1)CYCLE
    CALL psi_eval%interp(i,f,goptmp,gpsitmp(1:1))
    IF(gpsitmp(1)>psi_max)THEN
      psi_max=gpsitmp(1)
      r=ptmp(1)
      z=ptmp(2)
      cell_active=i
    END IF
  END DO
END IF
!RETURN
!---Setup gradient interpolation
psi_geval%u=>self%psi
CALL psi_geval%setup(self%fe_rep)
psi_geval_active=>psi_geval
psi_g2eval%u=>self%psi
CALL psi_g2eval%setup(self%fe_rep)
psi_g2eval_active=>psi_g2eval
!---Use MINPACK to find maximum (zero gradient)
ncons=2
ncofs=2
allocate(diag(ncofs),fjac(ncons,ncofs))
allocate(qtf(ncofs),wa1(ncofs),wa2(ncofs))
allocate(wa3(ncofs),wa4(ncons))
allocate(ipvt(ncofs))
mode = 1
factor = 1.d0
maxfev = 100
ftol = 1.d-9
xtol = 1.d-8
gtol = 1.d-8
epsfcn = 1.d-4
nprint = 0
ldfjac = ncons
ptmp(1)=r
ptmp(2)=z
call lmder(psimax_error_grad,ncons,ncofs,ptmp,gpsitmp,fjac,ldfjac, &
           ftol,xtol,gtol,maxfev,diag,mode,factor,nprint,info,nfev,njev, &
           ipvt,qtf,wa1,wa2,wa3,wa4)
deallocate(diag,fjac,qtf,wa1,wa2)
deallocate(wa3,wa4,ipvt)
!---Get axis values
call bmesh_findcell(self%fe_rep%mesh,cell_active,ptmp,f)
call self%fe_rep%mesh%jacobian(cell_active,f,goptmp,v)
call psi_eval%interp(cell_active,f,goptmp,gpsitmp(1:1))
psi_max=gpsitmp(1)
r=ptmp(1)
z=ptmp(2)
call psi_eval%delete()
call psi_geval%delete()
call psi_g2eval%delete()
end subroutine gs_psimax
!---------------------------------------------------------------------------
! SUBROUTINE gs_find_xpoint
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gs_find_xpoint(self,psi_x,r,z)
class(gs_eq), intent(inout) :: self
real(8), intent(inout) :: psi_x
real(8), intent(inout) :: r
real(8), intent(inout) :: z
type(oft_lag_brinterp), target :: psi_eval
type(oft_lag_bginterp), target :: psi_geval
type(oft_lag_bg2interp), target :: psi_g2eval
integer(4) :: i
real(8) :: ptmp(3),gpsitmp(3),f(3),goptmp(3,3),v,mag_min
!---MINPACK variables
real(8) :: ftol,xtol,gtol,epsfcn,factor
real(8), allocatable, dimension(:) :: diag,wa1,wa2,wa3,wa4,qtf
real(8), allocatable, dimension(:,:) :: fjac
integer(4) :: maxfev,mode,nprint,info,nfev,ldfjac,ncons,ncofs,njev
integer(4), allocatable, dimension(:) :: ipvt
!---Determine initial guess from cell search
psi_eval%u=>self%psi
CALL psi_eval%setup(self%fe_rep)
psi_geval%u=>self%psi
CALL psi_geval%setup(self%fe_rep)
f=1.d0/3.d0
mag_min=1.d99
psi_x=-1.d99
goptmp=1.d0
cell_active=0
DO i=1,self%fe_rep%mesh%nc
  ptmp = self%fe_rep%mesh%log2phys(i,f)
  call self%fe_rep%mesh%jacobian(i,f,goptmp,v)
  CALL psi_geval%interp(i,f,goptmp,gpsitmp)
  IF(magnitude(gpsitmp(1:2))<mag_min.AND.SQRT(SUM((ptmp(1:2)-self%o_point)**2))>2.d-2)THEN
    r=ptmp(1)
    z=ptmp(2)
    cell_active=i
    mag_min=magnitude(gpsitmp(1:2))
  END IF
END DO
IF(cell_active==0)RETURN
!---Setup gradient interpolation
psi_geval%u=>self%psi
psi_geval_active=>psi_geval
psi_g2eval%u=>self%psi
CALL psi_g2eval%setup(self%fe_rep)
psi_g2eval_active=>psi_g2eval
!---Use MINPACK to find maximum (zero gradient)
ncons=2
ncofs=2
allocate(diag(ncofs),fjac(ncons,ncofs))
allocate(qtf(ncofs),wa1(ncofs),wa2(ncofs))
allocate(wa3(ncofs),wa4(ncons))
allocate(ipvt(ncofs))
mode = 1
factor = 1.d0
maxfev = 100
ftol = 1.d-9
xtol = 1.d-8
gtol = 1.d-8
epsfcn = 1.d-4
nprint = 0
ldfjac = ncons
ptmp(1)=r
ptmp(2)=z
call lmder(psimax_error_grad,ncons,ncofs,ptmp,gpsitmp,fjac,ldfjac, &
           ftol,xtol,gtol,maxfev,diag,mode,factor,nprint,info,nfev,njev, &
           ipvt,qtf,wa1,wa2,wa3,wa4)
deallocate(diag,fjac,qtf,wa1,wa2)
deallocate(wa3,wa4,ipvt)
!---Get axis values
call psi_eval%interp(cell_active,f,goptmp,gpsitmp(1:1))
psi_x=gpsitmp(1)
r=ptmp(1)
z=ptmp(2)
end subroutine gs_find_xpoint
!---------------------------------------------------------------------------
!> Get q profile for equilibrium
!---------------------------------------------------------------------------
subroutine gs_get_qprof(gseq,nr,psi_q,prof,dl,rbounds,zbounds,ravgs)
class(gs_eq), intent(inout) :: gseq !< G-S object
integer(4), intent(in) :: nr !< Number of flux surfaces to sample
real(8), intent(in) :: psi_q(nr) !< Locations to sample in normalized flux
real(8), intent(out) :: prof(nr) !< q value at each sampling location
real(8), intent(out) :: dl !< Arc length of surface `psi_q(1)`
real(8), intent(out) :: rbounds(2,2) !< Radial bounds of surface `psi_q(1)`
real(8), intent(out) :: zbounds(2,2) !< Vertical bounds of surface `psi_q(1)`
real(8), optional, intent(out) :: ravgs(nr,2) !< Flux surface averages <R> and <1/R>
real(8) :: psi_surf,rmax,x1,x2,raxis,zaxis,fpol,qpsi
real(8) :: pt(3),pt_last(3),f(3),psi_tmp(1),gop(3,3)
type(oft_lag_brinterp), target :: psi_int
real(8), pointer :: ptout(:,:)
real(8), parameter :: tol=1.d-10
integer(4) :: i,j,cell
type(gsinv_interp), pointer :: field
CHARACTER(LEN=OFT_ERROR_SLEN) :: error_str
!---
raxis=gseq%o_point(1)
zaxis=gseq%o_point(2)
x1=0.d0; x2=1.d0
IF(gseq%plasma_bounds(1)>-1.d98)THEN
  x1=gseq%plasma_bounds(1); x2=gseq%plasma_bounds(2)
!   x1 = x1 + (x2-x1)*2.d-2
!   x2 = x2 + (x1-x2)*2.d-2
END IF
! IF(.NOT.gseq%free)x1 = x1 + (x2-x1)*2.d-2
psi_int%u=>gseq%psi
CALL psi_int%setup(gseq%fe_rep)
!---Find Rmax along Zaxis
rmax=raxis
cell=0
DO j=1,100
  pt=[(gseq%rmax-raxis)*j/REAL(100,8)+raxis,zaxis,0.d0]
  CALL bmesh_findcell(gseq%fe_rep%mesh,cell,pt,f)
  IF( (MAXVAL(f)>1.d0+tol) .OR. (MINVAL(f)<-tol) )EXIT
  CALL psi_int%interp(cell,f,gop,psi_tmp)
  IF( psi_tmp(1) < x1)EXIT
  rmax=pt(1)
END DO
pt_last=[(.1d0*rmax+.9d0*raxis),zaxis,0.d0]
!---
IF(oft_debug_print(1))THEN
  WRITE(*,'(2A)')oft_indent,'Axis Position:'
  CALL oft_increase_indent
  WRITE(*,'(2A,ES11.3)')oft_indent,'R    = ',raxis
  WRITE(*,'(2A,ES11.3)')oft_indent,'Z    = ',zaxis
  WRITE(*,'(2A,ES11.3)')oft_indent,'Rmax = ',rmax
  CALL oft_decrease_indent
END IF
!---Trace
call set_tracer(1)
!$omp parallel private(psi_surf,pt,ptout,fpol,qpsi,field) firstprivate(pt_last)
ALLOCATE(field)
field%u=>gseq%psi
CALL field%setup(gseq%fe_rep)
IF(PRESENT(ravgs))THEN
  field%compute_geom=.TRUE.
  active_tracer%neq=5
ELSE
  field%compute_geom=.FALSE.
  active_tracer%neq=3
END IF
active_tracer%B=>field
active_tracer%maxsteps=8e4
active_tracer%raxis=raxis
active_tracer%zaxis=zaxis
active_tracer%inv=.TRUE.
ALLOCATE(ptout(3,active_tracer%maxsteps+1))
!$omp do schedule(dynamic,1)
do j=1,nr
  !---------------------------------------------------------------------------
  ! Trace contour
  !---------------------------------------------------------------------------
  ! psi_surf=(x2-x1)*((j-1)/REAL(nr,8))
  ! psi_surf=x2 - psi_surf
  psi_surf=psi_q(j)*(x2-x1) + x1
  IF(gseq%diverted.AND.psi_q(j)<0.02d0)THEN ! Use higher tracing tolerance near divertor
    active_tracer%tol=1.d-10
  ELSE
    active_tracer%tol=1.d-8
  END IF
  !
  pt=pt_last
  !$omp critical
  CALL gs_psi2r(gseq,psi_surf,pt,psi_int)
  !$omp end critical
  CALL tracinginv_fs(gseq%fe_rep%mesh,pt(1:2),ptout)
  pt_last=pt
  !---Skip point if trace fails
  if(active_tracer%status/=1)THEN
    WRITE(error_str,"(A,F10.4)")"gs_get_qprof: Trace did not complete at psi = ",1.d0-psi_q(j)
    CALL oft_warn(error_str)
    CYCLE
  end if
  IF(j==1)THEN
    dl = 0.d0
    rbounds(:,1)=ptout(2:3,1); rbounds(:,2)=ptout(2:3,1)
    zbounds(:,1)=ptout(2:3,1); zbounds(:,2)=ptout(2:3,1)
    DO i=2,active_tracer%nsteps
      dl = dl + SQRT(SUM((ptout(2:3,i)-ptout(2:3,i-1))**2))
      IF(ptout(2,i)<rbounds(1,1))THEN
        rbounds(:,1)=ptout(2:3,i)
      ELSE IF(ptout(2,i)>rbounds(1,2))THEN
        rbounds(:,2)=ptout(2:3,i)
      END IF
      IF(ptout(3,i)<zbounds(2,1))THEN
        zbounds(:,1)=ptout(2:3,i)
      ELSE IF(ptout(3,i)>zbounds(2,2))THEN
        zbounds(:,2)=ptout(2:3,i)
      END IF
    END DO
    IF(active_tracer%status/=1)dl=-1.d0
  END IF
  !---Get flux variables
  IF(gseq%mode==0)THEN
    fpol=gseq%alam*gseq%I%f(psi_surf)+gseq%I%f_offset
  ELSE
    fpol=SQRT(gseq%alam*gseq%I%f(psi_surf) + gseq%I%f_offset**2) &
    + gseq%I%f_offset*(1.d0-SIGN(1.d0,gseq%I%f_offset))
  END IF
  !---Safety Factor (q)
  qpsi=fpol*active_tracer%v(3)/(2*pi)
  prof(j)=qpsi
  IF(PRESENT(ravgs))THEN
    ravgs(j,1)=active_tracer%v(4)/active_tracer%v(2)
    ravgs(j,2)=active_tracer%v(5)/active_tracer%v(2)
  END IF
end do
CALL active_tracer%delete
DEALLOCATE(ptout)
CALL field%delete()
DEALLOCATE(field)
!$omp end parallel
CALL psi_int%delete()
end subroutine gs_get_qprof
!---------------------------------------------------------------------------
!>
!---------------------------------------------------------------------------
subroutine gs_trace_surf(gseq,psi_in,points,npoints)
class(gs_eq), intent(inout) :: gseq
real(8), intent(in) :: psi_in
real(8), pointer, dimension(:,:), intent(out) :: points
integer(4), intent(out) :: npoints
real(8) :: rmax,x1,x2,raxis,zaxis,fpol,qpsi,pmin,psi_surf
real(8) :: pt(3),pt_last(3),f(3),psi_tmp(1),gop(3,3)
type(oft_lag_brinterp) :: psi_int
real(8), pointer :: ptout(:,:)
real(8), parameter :: tol=1.d-10
integer(4) :: i,j,cell
type(gsinv_interp), target :: field
!---
raxis=gseq%o_point(1)
zaxis=gseq%o_point(2)
x1=0.d0; x2=1.d0
IF(gseq%plasma_bounds(1)>-1.d98)THEN
  x1=gseq%plasma_bounds(1); x2=gseq%plasma_bounds(2)
END IF
psi_surf = x1 + (x2-x1)*psi_in
psi_int%u=>gseq%psi
CALL psi_int%setup(gseq%fe_rep)
!---Find Rmax along Zaxis
rmax=raxis
cell=0
pmin=1.d99
DO j=1,100
  pt=[(gseq%rmax-raxis)*j/REAL(100,8)+raxis,zaxis,0.d0]
  CALL bmesh_findcell(gseq%fe_rep%mesh,cell,pt,f)
  IF( (MAXVAL(f)>1.d0+tol) .OR. (MINVAL(f)<-tol) )EXIT
  CALL psi_int%interp(cell,f,gop,psi_tmp)
  IF( ABS(psi_tmp(1)-psi_surf)<pmin)THEN
    pmin = ABS(psi_tmp(1)-psi_surf)
    rmax = pt(1)
  END IF
END DO
pt_last=[rmax,zaxis,0.d0]
! !---
! IF(oft_debug_print(1))THEN
!   WRITE(*,'(2A)')oft_indent,'Axis Position:'
!   CALL oft_increase_indent
!   WRITE(*,'(2A,ES11.3)')oft_indent,'R    = ',raxis
!   WRITE(*,'(2A,ES11.3)')oft_indent,'Z    = ',zaxis
!   WRITE(*,'(2A,ES11.3)')oft_indent,'Rmax = ',rmax
!   CALL oft_decrease_indent
! END IF
!---Trace
call set_tracer(1)
!!$omp parallel private(j,psi_surf,pt,ptout,fpol,qpsi,field) firstprivate(pt_last)
field%u=>gseq%psi
CALL field%setup(gseq%fe_rep)
active_tracer%neq=3
active_tracer%B=>field
active_tracer%maxsteps=8e4
IF(gseq%diverted.AND.psi_in<0.02d0)THEN ! Use higher tracing tolerance near divertor
  active_tracer%tol=1.d-10
ELSE
  active_tracer%tol=1.d-8
END IF
active_tracer%raxis=raxis
active_tracer%zaxis=zaxis
active_tracer%inv=.TRUE.
ALLOCATE(ptout(3,active_tracer%maxsteps+1))
pt=pt_last
!!$omp critical
CALL gs_psi2r(gseq,psi_surf,pt)
!!$omp end critical
CALL tracinginv_fs(gseq%fe_rep%mesh,pt(1:2),ptout)
!---Skip point if trace fails
if(active_tracer%status/=1)THEN
  ! CALL oft_warn("gs_trace_surf: Trace did not complete")
  npoints=-1
else
  ALLOCATE(points(2,active_tracer%nsteps))
  points=ptout(2:3,1:active_tracer%nsteps)
  npoints=active_tracer%nsteps
end if
CALL active_tracer%delete
DEALLOCATE(ptout)
!!$omp end parallel
CALL psi_int%delete
end subroutine gs_trace_surf
!---------------------------------------------------------------------------
! SUBROUTINE gs_set_cond_weights
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gs_set_cond_weights(self,vals,skip_fixed)
class(gs_eq), intent(inout) :: self
real(8), intent(in) :: vals(:)
logical, intent(in) :: skip_fixed
integer(4) :: i,j,k,kk
!---
k=0
kk=0
DO i=1,self%ncond_regs
  IF(self%cond_regions(i)%pair<0)THEN
    DO j=1,self%cond_regions(i)%neigs
      k=k+1
      IF(self%cond_regions(i)%fixed(j).AND.skip_fixed)CYCLE
      kk=kk+1
      self%cond_weights(k)=vals(kk)
      self%cond_regions(i)%weights(j)=self%cond_weights(k)
    END DO
  ELSE
    DO j=1,self%cond_regions(i)%neigs
      self%cond_regions(i)%weights(j)=self%cond_regions(i)%pair_signs(j)* &
        self%cond_regions(self%cond_regions(i)%pair)%weights(j)
    END DO
  END IF
END DO
end subroutine gs_set_cond_weights
!---------------------------------------------------------------------------
! SUBROUTINE gs_get_cond_weights
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gs_get_cond_weights(self,vals,skip_fixed)
class(gs_eq), intent(inout) :: self
real(8), intent(inout) :: vals(:)
logical, intent(in) :: skip_fixed
integer(4) :: i,j,k,kk
!---
k=0
kk=0
DO i=1,self%ncond_regs
  IF(self%cond_regions(i)%pair<0)THEN
    DO j=1,self%cond_regions(i)%neigs
      k=k+1
      self%cond_weights(k)=self%cond_regions(i)%weights(j)
      IF(self%cond_regions(i)%fixed(j).AND.skip_fixed)CYCLE
      kk=kk+1
      vals(kk)=self%cond_weights(k)
    END DO
  END IF
END DO
end subroutine gs_get_cond_weights
!---------------------------------------------------------------------------
! SUBROUTINE gs_get_cond_scales
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gs_get_cond_scales(self,vals,skip_fixed)
class(gs_eq), intent(inout) :: self
real(8), intent(inout) :: vals(:)
logical, intent(in) :: skip_fixed
integer(4) :: i,j,k
k=0
DO i=1,self%ncond_regs
  IF(self%cond_regions(i)%pair<0)THEN
    DO j=1,self%cond_regions(i)%neigs
      IF(self%cond_regions(i)%fixed(j).AND.skip_fixed)CYCLE
      k=k+1
      vals(k) = self%cond_regions(i)%fit_scales(j)
    END DO
  END IF
END DO
end subroutine gs_get_cond_scales
!---------------------------------------------------------------------------
! FUNCTION gs_eta_spitzer
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
function gs_eta_spitzer(T,lambda) result(eta)
real(8), intent(in) :: T
real(8), optional, intent(in) :: lambda
real(8) :: eta,lam
lam=EXP(11.d0)
IF(PRESENT(lambda))lam=lambda
!---
eta = 1.65d-9*LOG(lam)/(T**1.5d0) ! Spitzer 1.65e-9 * Log(Lambda)/(T^3/2) [T in KeV]
end function gs_eta_spitzer
!
subroutine gs_prof_interp_setup(self,gs)
class(gs_prof_interp), intent(inout) :: self
class(gs_eq), target, intent(inout) :: gs
self%gs=>gs
self%mesh=>gs%mesh
ALLOCATE(self%psi_eval,self%psi_geval)
self%psi_eval%u=>self%gs%psi
CALL self%psi_eval%setup(self%gs%fe_rep)
CALL self%psi_geval%shared_setup(self%psi_eval)
end subroutine gs_prof_interp_setup
!
subroutine gs_prof_interp_delete(self)
class(gs_prof_interp), intent(inout) :: self
CALL self%psi_eval%delete()
CALL self%psi_geval%delete()
DEALLOCATE(self%psi_eval,self%psi_geval)
NULLIFY(self%gs,self%mesh)
end subroutine gs_prof_interp_delete
!---------------------------------------------------------------------------
! SUBROUTINE gs_rinterp
!---------------------------------------------------------------------------
!> Reconstruct a Grad-Shafranov field
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [1]
!---------------------------------------------------------------------------
subroutine gs_prof_interp_apply(self,cell,f,gop,val)
class(gs_prof_interp), intent(inout) :: self
integer(4), intent(in) :: cell
real(8), intent(in) :: f(:)
real(8), intent(in) :: gop(3,3)
real(8), intent(out) :: val(:)
real(8) :: psitmp(1),gpsitmp(3),pt(3)
logical :: in_plasma
!---
pt=self%gs%fe_rep%mesh%log2phys(cell,f)
in_plasma=.TRUE.
IF(gs_test_bounds(self%gs,pt).AND.(self%gs%fe_rep%mesh%reg(cell)==1))THEN
  in_plasma=.TRUE.
ELSE
  in_plasma=.FALSE.
END IF
!---
SELECT CASE(self%mode)
  CASE(1)
    CALL self%psi_eval%interp(cell,f,gop,psitmp)
    val(1)=self%gs%psiscale*psitmp(1)
  CASE(2)
    pt=self%gs%fe_rep%mesh%log2phys(cell,f)
    CALL self%psi_eval%interp(cell,f,gop,psitmp)
    IF(in_plasma.AND.(psitmp(1)>self%gs%plasma_bounds(1)))THEN
      IF(self%gs%mode==0)THEN
        val(1)=self%gs%psiscale*self%gs%alam*(self%gs%I%f(psitmp(1))+self%gs%I%f_offset/self%gs%alam)
      ELSE
        val(1)=self%gs%psiscale*SQRT(self%gs%alam*self%gs%I%f(psitmp(1)) + self%gs%I%f_offset**2)
      END IF
    ELSE
      val(1)=self%gs%psiscale*self%gs%I%f_offset
    END IF
  CASE(3)
    CALL self%psi_eval%interp(cell,f,gop,psitmp)
    IF(in_plasma.AND.(psitmp(1)>self%gs%plasma_bounds(1)))THEN
      val(1)=(self%gs%psiscale**2)*self%gs%pnorm*self%gs%P%F(psitmp(1))/mu0
    ELSE
      val(1)=0.d0
    END IF
  CASE DEFAULT
    CALL oft_abort('Unkown field mode','gs_prof_interp_apply',__FILE__)
END SELECT
end subroutine gs_prof_interp_apply
!---------------------------------------------------------------------------
!> Reconstruct a Grad-Shafranov field
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [1]
!---------------------------------------------------------------------------
subroutine gs_b_interp_apply(self,cell,f,gop,val)
class(gs_b_interp), intent(inout) :: self
integer(4), intent(in) :: cell
real(8), intent(in) :: f(:)
real(8), intent(in) :: gop(3,3)
real(8), intent(out) :: val(:)
real(8) :: psitmp(1),gpsitmp(3),pt(3)
logical :: in_plasma
pt=self%gs%fe_rep%mesh%log2phys(cell,f)
in_plasma=.TRUE.
IF(gs_test_bounds(self%gs,pt).AND.(self%gs%fe_rep%mesh%reg(cell)==1))THEN
  in_plasma=.TRUE.
ELSE
  in_plasma=.FALSE.
END IF
! Sample fields
CALL self%psi_eval%interp(cell,f,gop,psitmp)
CALL self%psi_geval%interp(cell,f,gop,gpsitmp)
! Evaluate B-field
val(1)=-self%gs%psiscale*gpsitmp(2)/(pt(1)+gs_epsilon)
val(3)=self%gs%psiscale*gpsitmp(1)/(pt(1)+gs_epsilon)
IF(in_plasma.AND.(psitmp(1)>self%gs%plasma_bounds(1)))THEN
  IF(self%gs%mode==0)THEN
    val(2)=self%gs%psiscale*self%gs%alam*(self%gs%I%f(psitmp(1))+self%gs%I%f_offset/self%gs%alam)/(pt(1)+gs_epsilon)
  ELSE
    val(2)=self%gs%psiscale*SQRT(self%gs%alam*self%gs%I%f(psitmp(1)) + self%gs%I%f_offset**2)/(pt(1)+gs_epsilon)
  END IF
ELSE
  val(2)=self%gs%psiscale*self%gs%I%f_offset/(pt(1)+gs_epsilon)
END IF
end subroutine gs_b_interp_apply
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gsinv_setup(self,lag_rep)
class(gsinv_interp), intent(inout) :: self
class(oft_afem_type), target, intent(inout) :: lag_rep
SELECT TYPE(lag_rep)
CLASS IS(oft_scalar_bfem)
  self%lag_rep=>lag_rep
CLASS DEFAULT
  CALL oft_abort("Incorrect FE type","gsinv_setup",__FILE__)
END SELECT
self%mesh=>self%lag_rep%mesh
NULLIFY(self%uvals)
CALL self%u%get_local(self%uvals)
end subroutine gsinv_setup
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gsinv_destroy(self)
class(gsinv_interp), intent(inout) :: self
IF(ASSOCIATED(self%uvals))DEALLOCATE(self%uvals)
end subroutine gsinv_destroy
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gsinv_apply(self,cell,f,gop,val)
class(gsinv_interp), intent(inout) :: self
integer(4), intent(in) :: cell
real(8), intent(in) :: f(:)
real(8), intent(in) :: gop(3,3)
real(8), intent(out) :: val(:)
integer(4), allocatable :: j(:)
integer(4) :: jc
real(8) :: rop(3),d2op(6),pt(3),grad(3),tmp
real(8) :: s,c
!---Get dofs
allocate(j(self%lag_rep%nce))
call self%lag_rep%ncdofs(cell,j)
!---Reconstruct gradient
grad=0.d0
do jc=1,self%lag_rep%nce
  call oft_blag_geval(self%lag_rep,cell,jc,f,rop,gop)
  grad=grad+self%uvals(j(jc))*rop
end do
!---Get radial position
pt=self%mesh%log2phys(cell,f)
!---
s=SIN(self%t)
c=COS(self%t)
!---Position
val(1)=(self%rho*(grad(1)*s-grad(2)*c))/(grad(1)*c+grad(2)*s)
val(2)=pt(1)*SQRT((self%rho**2+val(1)**2)/SUM(grad**2))
!---Safety Factor
val(3)=val(2)/pt(1)**2
IF(self%compute_geom)THEN
  val(4)=pt(1)*val(2)
  val(5)=val(2)/pt(1)
END IF
! val(3:8)=val(3:8)
deallocate(j)
end subroutine gsinv_apply
!---------------------------------------------------------------------------
! SUBROUTINE gs_save_fgrid
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gs_save_fgrid(self,filename)
class(gs_eq), target, intent(inout) :: self
character(LEN=*), optional, intent(in) :: filename
class(oft_vector), pointer :: a,br,bt,bz
real(r8), pointer :: vals_tmp(:)
class(oft_solver), pointer :: solver
type(gs_b_interp) :: Bfield
type(gs_prof_interp) :: field
integer(4) :: i,io_unit
real(8), allocatable :: Fout(:,:)
logical :: pm_save
!---
NULLIFY(vals_tmp)
ALLOCATE(Fout(5,self%fe_rep%ne))
CALL self%psi%new(a)
CALL self%psi%new(br)
CALL self%psi%new(bt)
CALL self%psi%new(bz)
! Bfield%gs=>self
! field%gs=>self
CALL Bfield%setup(self)
CALL field%setup(self)
!---Setup Solver
CALL create_cg_solver(solver)
solver%A=>self%mop
solver%its=-2
CALL create_diag_pre(solver%pre)
pm_save=oft_env%pm; oft_env%pm=.FALSE.
!---Project B-field
CALL oft_blag_vproject(self%fe_rep,Bfield,br,bt,bz)
CALL a%set(0.d0)
CALL solver%apply(a,br)
CALL a%get_local(vals_tmp)
Fout(1,:)=vals_tmp
CALL a%set(0.d0)
CALL solver%apply(a,bt)
CALL a%get_local(vals_tmp)
Fout(2,:)=vals_tmp
CALL a%set(0.d0)
CALL solver%apply(a,bz)
CALL a%get_local(vals_tmp)
Fout(3,:)=vals_tmp
!---Project pressure
field%mode=3
CALL oft_blag_project(self%fe_rep,field,br)
CALL a%set(0.d0)
CALL solver%apply(a,br)
CALL a%get_local(vals_tmp)
Fout(4,:)=vals_tmp
!---Get poloidal flux
CALL self%psi%get_local(vals_tmp)
Fout(5,:)=vals_tmp
!---Output
IF(PRESENT(filename))THEN
  OPEN(NEWUNIT=io_unit,FILE=TRIM(filename))
  WRITE(io_unit,*)self%fe_rep%order
  DO i=1,self%fe_rep%ne
    WRITE(io_unit,'(5E18.10)')Fout(:,i)
  END DO
  CLOSE(io_unit)
END IF
!---Save fields to plot
CALL self%fe_rep%mesh%save_vertex_scalar(Fout(1,:),self%xdmf,'Br')
CALL self%fe_rep%mesh%save_vertex_scalar(Fout(2,:),self%xdmf,'Bt')
CALL self%fe_rep%mesh%save_vertex_scalar(Fout(3,:),self%xdmf,'Bz')
CALL self%fe_rep%mesh%save_vertex_scalar(Fout(4,:),self%xdmf,'P')
!---Clean up
oft_env%pm=pm_save
CALL a%delete
CALL br%delete
CALL bt%delete
CALL bz%delete
DEALLOCATE(a,br,bt,bz,Fout,vals_tmp)
CALL Bfield%delete()
CALL field%delete()
CALL solver%pre%delete
DEALLOCATE(solver%pre)
CALL solver%delete
DEALLOCATE(solver)
end subroutine gs_save_fgrid
!---------------------------------------------------------------------------
! SUBROUTINE gs_save_prof
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine gs_save_prof(self,filename,mpsi_sample)
class(gs_eq), target, intent(inout) :: self
character(LEN=*), intent(in) :: filename
INTEGER(4), OPTIONAL, intent(in) :: mpsi_sample
integer(4) :: i,m,io_unit
real(8) :: x1,x2,outtmp(5),r
m=100
IF(PRESENT(mpsi_sample))m=mpsi_sample
!---
x1=0.d0; x2=1.d0
IF(self%plasma_bounds(1)>-1.d98)THEN
  x1=self%plasma_bounds(1); x2=self%plasma_bounds(2)
END IF
100 FORMAT (5E18.10)
OPEN(NEWUNIT=io_unit,FILE=TRIM(filename))
DO i=0,m
  r=i*(x2-x1)/m + x1
  IF(i==0)r = r + (x2-x1)*1.d-10
  outtmp(1)=self%psiscale*r
  IF(self%mode==0)THEN
    outtmp(2)=self%alam*self%I%fp(r)
    outtmp(3)=self%psiscale*self%alam*self%I%f(r) + self%I%f_offset
  ELSE
    outtmp(3)=SQRT(self%psiscale*self%alam*self%I%f(r) + self%I%f_offset**2)
    outtmp(2)=self%alam*self%I%fp(r)/(2.d0*outtmp(3))
  END IF
  outtmp(4)=self%psiscale*self%pnorm*self%P%fp(r)
  outtmp(5)=self%psiscale*self%psiscale*self%pnorm*self%P%f(r)
  WRITE(io_unit,100)outtmp
END DO
CLOSE(io_unit)
end subroutine gs_save_prof
!---------------------------------------------------------------------------
! SUBROUTINE: build_mrop
!---------------------------------------------------------------------------
!> Construct mass matrix for a boundary Lagrange scalar representation
!!
!! Supported boundary conditions
!! - `'none'` Full matrix
!! - `'zerob'` Dirichlet for all boundary DOF
!!
!! @param[in,out] mat Matrix object
!! @param[in] bc Boundary condition
!---------------------------------------------------------------------------
subroutine build_mrop(fe_rep,mat,bc)
type(oft_scalar_bfem), intent(inout) :: fe_rep
class(oft_matrix), pointer, intent(inout) :: mat
character(LEN=*), intent(in) :: bc
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time,pt(3)
real(r8), allocatable :: rop(:),mop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_lag_vec
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing Boundary LAG::MOP'
  CALL mytimer%tick()
END IF
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL fe_rep%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,rop,det,mop,curved,goptmp,m,vol,jc,jr,pt)
allocate(j(fe_rep%nce)) ! Local DOF and matrix indices
allocate(rop(fe_rep%nce)) ! Reconstructed gradient operator
allocate(mop(fe_rep%nce,fe_rep%nce)) ! Local laplacian matrix
!$omp do
! schedule(dynamic,1) ordered
do i=1,fe_rep%mesh%nc
  !---Get local reconstructed operators
  mop=0.d0
  do m=1,fe_rep%quad%np ! Loop over quadrature points
    call fe_rep%mesh%jacobian(i,fe_rep%quad%pts(:,m),goptmp,vol)
    det=vol*fe_rep%quad%wts(m)
    pt=fe_rep%mesh%log2phys(i,fe_rep%quad%pts(:,m))
    do jc=1,fe_rep%nce ! Loop over degrees of freedom
      call oft_blag_eval(fe_rep,i,jc,fe_rep%quad%pts(:,m),rop(jc))
    end do
    !---Compute local matrix contributions
    do jr=1,fe_rep%nce
      do jc=1,fe_rep%nce
        mop(jr,jc) = mop(jr,jc) + rop(jr)*rop(jc)*det/(pt(1)+gs_epsilon)
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call fe_rep%ncdofs(i,j)
  !---Add local values to global matrix
  !!$omp ordered
  call mat%atomic_add_values(j,j,mop,fe_rep%nce,fe_rep%nce)
  !!$omp end ordered
end do
deallocate(j,rop,mop)
!$omp end parallel
CALL fe_rep%vec_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine build_mrop
!---------------------------------------------------------------------------
! SUBROUTINE: build_dels
!---------------------------------------------------------------------------
!> Construct laplacian matrix for Lagrange scalar representation
!!
!! Supported boundary conditions
!! - `'none'` Full matrix
!! - `'zerob'` Dirichlet for all boundary DOF
!! - `'grnd'`  Dirichlet for only groundin point
!!
!! @param[in,out] mat Matrix object
!! @param[in] bc Boundary condition
!---------------------------------------------------------------------------
subroutine build_dels(mat,self,bc,dt,scale)
class(oft_matrix), pointer, intent(inout) :: mat
class(gs_eq), intent(inout) :: self
character(LEN=*), intent(in) :: bc
real(8), optional, intent(in) :: dt
real(8), optional, intent(in) :: scale
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time,pt(3),dt_in,main_scale
real(r8), allocatable :: rop(:),gop(:,:),lop(:,:),eta_reg(:)
logical :: curved
integer(i4) :: nnonaxi
integer(i4), allocatable :: dense_flag(:)
real(r8), pointer, dimension(:) :: nonaxi_tmp
real(r8), pointer, dimension(:,:) :: nonaxi_vals
type(oft_1d_int), pointer, dimension(:) :: bc_nodes
CLASS(oft_vector), POINTER :: oft_lag_vec
TYPE(oft_graph_ptr) :: graphs(1,1)
TYPE(oft_graph), TARGET :: graph1,graph2
type(oft_timer) :: mytimer
CLASS(oft_bmesh), POINTER :: smesh
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing Boundary LAG::LOP'
  CALL mytimer%tick()
END IF
dt_in=-1.d0
IF(PRESENT(dt))dt_in=dt
main_scale=1.d0
IF(PRESENT(scale))main_scale=scale
smesh=>self%fe_rep%mesh
!
nnonaxi=0
IF(dt_in>0.d0)nnonaxi=self%region_info%nnonaxi
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  !---
  graph1%nr=self%fe_rep%ne
  graph1%nrg=self%fe_rep%global%ne
  graph1%nc=self%fe_rep%ne
  graph1%ncg=self%fe_rep%global%ne
  graph1%nnz=self%fe_rep%nee
  graph1%kr=>self%fe_rep%kee
  graph1%lc=>self%fe_rep%lee
  !---Add dense blocks for non-continuous regions
  IF(nnonaxi>0)THEN
    !---Add dense blocks
    ALLOCATE(dense_flag(self%fe_rep%ne))
    dense_flag=0
    DO m=1,self%region_info%nnonaxi
      dense_flag(self%region_info%noaxi_nodes(m)%v)=m
    END DO
    CALL graph_add_dense_blocks(graph1,graph2,dense_flag,self%region_info%noaxi_nodes)
    NULLIFY(graph1%kr,graph1%lc)
    graph1%nnz=graph2%nnz
    graph1%kr=>graph2%kr
    graph1%lc=>graph2%lc
    DEALLOCATE(dense_flag)
  END IF
  !---Add dense block for boundary
  IF(TRIM(bc)=="free")THEN
    ! CALL gs_mat_create(mat)
    ALLOCATE(bc_nodes(1))
    bc_nodes(1)%n=self%fe_rep%nbe
    bc_nodes(1)%v=>self%fe_rep%lbe
    ALLOCATE(dense_flag(self%fe_rep%ne))
    dense_flag=0
    dense_flag(bc_nodes(1)%v)=1
    !---Add dense blocks
    CALL graph_add_dense_blocks(graph1,graph2,dense_flag,bc_nodes)
    NULLIFY(graph1%kr,graph1%lc)
    graph1%nnz=graph2%nnz
    graph1%kr=>graph2%kr
    graph1%lc=>graph2%lc
    DEALLOCATE(dense_flag)
  END IF
  !---Create matrix
  graphs(1,1)%g=>graph1
  CALL self%fe_rep%vec_create(oft_lag_vec)
  CALL create_matrix(mat,graphs,oft_lag_vec,oft_lag_vec)
  CALL oft_lag_vec%delete
  DEALLOCATE(oft_lag_vec)
  NULLIFY(graphs(1,1)%g)
ELSE
  CALL mat%zero
END IF
ALLOCATE(eta_reg(smesh%nreg))
eta_reg=-1.d0
IF(dt_in>0.d0)THEN
  DO i=1,self%ncond_regs
    jr=self%cond_regions(i)%id
    eta_reg(jr)=self%cond_regions(i)%eta
  END DO
END IF
!---------------------------------------------------------------------------
! Operator integration
!---------------------------------------------------------------------------
IF(nnonaxi>0)THEN
  ALLOCATE(nonaxi_vals(self%region_info%block_max+1,self%region_info%nnonaxi))
  nonaxi_vals=0.d0
END IF
!$omp parallel private(j,rop,gop,det,lop,curved,goptmp,m,vol,jc,jr,pt,nonaxi_tmp)
allocate(j(self%fe_rep%nce)) ! Local DOF and matrix indices
allocate(rop(self%fe_rep%nce),gop(3,self%fe_rep%nce)) ! Reconstructed gradient operator
allocate(lop(self%fe_rep%nce,self%fe_rep%nce)) ! Local laplacian matrix
IF(nnonaxi>0)allocate(nonaxi_tmp(self%fe_rep%nce))
!$omp do schedule(dynamic,1) ordered
do i=1,self%fe_rep%mesh%nc
  !---Get local reconstructed operators
  lop=0.d0
  IF(nnonaxi>0)nonaxi_tmp=0.d0
  do m=1,self%fe_rep%quad%np ! Loop over quadrature points
    call self%fe_rep%mesh%jacobian(i,self%fe_rep%quad%pts(:,m),goptmp,vol)
    det=vol*self%fe_rep%quad%wts(m)
    pt=smesh%log2phys(i,self%fe_rep%quad%pts(:,m))
    do jc=1,self%fe_rep%nce ! Loop over degrees of freedom
      call oft_blag_eval(self%fe_rep,i,jc,self%fe_rep%quad%pts(:,m),rop(jc))
      call oft_blag_geval(self%fe_rep,i,jc,self%fe_rep%quad%pts(:,m),gop(:,jc),goptmp)
    end do
    !---Compute local matrix contributions
    do jr=1,self%fe_rep%nce
      do jc=1,self%fe_rep%nce
        lop(jr,jc) = lop(jr,jc) + DOT_PRODUCT(gop(1:2,jr),gop(1:2,jc))*det/(pt(1)+gs_epsilon)
      end do
    end do
    IF(dt_in>0.d0.AND.eta_reg(smesh%reg(i))>0.d0)THEN
      do jr=1,self%fe_rep%nce
        do jc=1,self%fe_rep%nce
          lop(jr,jc) = lop(jr,jc) + rop(jr)*rop(jc)*det/(pt(1)+gs_epsilon)/(dt_in*eta_reg(smesh%reg(i)))
        end do
        IF(nnonaxi>0.AND.self%region_info%reg_map(smesh%reg(i))>0)THEN
          nonaxi_tmp(jr)=nonaxi_tmp(jr) + rop(jr)*det/(pt(1)+gs_epsilon)/(dt_in*eta_reg(smesh%reg(i)))
        END IF
      end do
    END IF
    IF(nnonaxi>0.AND.self%region_info%reg_map(smesh%reg(i))>0)THEN
      !$omp atomic
      nonaxi_vals(self%region_info%block_max+1,self%region_info%reg_map(smesh%reg(i))) = &
        nonaxi_vals(self%region_info%block_max+1,self%region_info%reg_map(smesh%reg(i))) + det
    END IF
  end do
  !---Get local to global DOF mapping
  call self%fe_rep%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,self%fe_rep%nce
        IF(self%fe_rep%be(j(jr)))lop(jr,:)=0.d0
      END DO
    CASE("free")
      DO jr=1,self%fe_rep%nce
        IF(self%fe_rep%be(j(jr)))lop(jr,:)=0.d0
      END DO
  END SELECT
  !---Add local values to global matrix
  lop=lop*main_scale
  !$omp ordered
  call mat%atomic_add_values(j,j,lop,self%fe_rep%nce,self%fe_rep%nce)
  IF(nnonaxi>0.AND.self%region_info%reg_map(smesh%reg(i))>0)THEN
    DO jr=1,self%fe_rep%nce
      !$omp atomic
      nonaxi_vals(self%region_info%node_mark(smesh%reg(i),j(jr)),self%region_info%reg_map(smesh%reg(i))) = &
        nonaxi_vals(self%region_info%node_mark(smesh%reg(i),j(jr)),self%region_info%reg_map(smesh%reg(i))) &
        + nonaxi_tmp(jr)
    END DO
  END IF
  !$omp end ordered
end do
IF(nnonaxi>0)THEN
  !$omp do schedule(dynamic,1) ordered
  do i=1,self%fe_rep%mesh%nc
    IF(self%region_info%reg_map(smesh%reg(i))==0)CYCLE
    !---Get local to global DOF mapping
    call self%fe_rep%ncdofs(i,j)
    !---Get local reconstructed operators
    lop(:,1)=0.d0
    do m=1,self%fe_rep%quad%np ! Loop over quadrature points
      call self%fe_rep%mesh%jacobian(i,self%fe_rep%quad%pts(:,m),goptmp,vol)
      det=vol*self%fe_rep%quad%wts(m)
      do jc=1,self%fe_rep%nce ! Loop over degrees of freedom
        call oft_blag_eval(self%fe_rep,i,jc,self%fe_rep%quad%pts(:,m),rop(jc))
        lop(jc,1) = lop(jc,1) + rop(jc)*det
      end do
    end do
    !---Add local values to global matrix
    m=self%region_info%reg_map(smesh%reg(i))
    lop(:,1)=-lop(:,1)/nonaxi_vals(self%region_info%block_max+1,m)
    lop(:,1)=lop(:,1)*main_scale
    !$omp ordered
    do jc=1,self%fe_rep%nce
      call mat%atomic_add_values(j(jc:jc),self%region_info%noaxi_nodes(m)%v, &
        lop(jc,1)*nonaxi_vals(:,m),1,self%region_info%noaxi_nodes(m)%n)
    end do
    !$omp end ordered
  end do
END IF
deallocate(j,rop,gop,lop)
IF(nnonaxi>0)deallocate(nonaxi_tmp)
!$omp end parallel
ALLOCATE(lop(1,1),j(1))
IF(nnonaxi>0)THEN
  IF(.NOT.ASSOCIATED(self%region_info%nonaxi_vals))ALLOCATE(self%region_info%nonaxi_vals(self%fe_rep%ne,smesh%nreg))
  self%region_info%nonaxi_vals=0.d0
  DO i=1,smesh%nreg
    IF(self%region_info%reg_map(i)==0)CYCLE
    DO jc=1,self%region_info%noaxi_nodes(self%region_info%reg_map(i))%n
      self%region_info%nonaxi_vals(self%region_info%noaxi_nodes(self%region_info%reg_map(i))%v(jc),i) = &
        -nonaxi_vals(jc,self%region_info%reg_map(i))*dt_in/nonaxi_vals(self%region_info%block_max+1,self%region_info%reg_map(i))
    END DO
  END DO
  DEALLOCATE(nonaxi_vals)
END IF
!---Set diagonal entries for dirichlet rows
SELECT CASE(TRIM(bc))
  CASE("zerob")
    lop(1,1)=1.d0
    DO i=1,self%fe_rep%nbe
      IF(.NOT.self%fe_rep%linkage%leo(i))CYCLE
      j=self%fe_rep%lbe(i)
      call mat%add_values(j,j,lop,1,1)
    END DO
  CASE("free")
    CALL set_bcmat(self,mat)
END SELECT
DEALLOCATE(j,lop)
CALL self%fe_rep%vec_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec,eta_reg)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine build_dels
!---------------------------------------------------------------------------
!> Initialize Grad-Shafranov solution with the Taylor state
!!
!! @param[in,out] self G-S object
!---------------------------------------------------------------------------
subroutine gs_destroy(self)
class(gs_eq), intent(inout) :: self
integer(i4) :: i,j
IF(ASSOCIATED(self%gs_zerob_bc%node_flag))DEALLOCATE(self%gs_zerob_bc%node_flag)
!
! IF(ASSOCIATED(self%cflag))DEALLOCATE(self%cflag)
IF(ASSOCIATED(self%lim_con))DEALLOCATE(self%lim_con)
IF(ASSOCIATED(self%lim_ptr))DEALLOCATE(self%lim_ptr)
IF(ASSOCIATED(self%limiter_nds))DEALLOCATE(self%limiter_nds)
IF(ASSOCIATED(self%bc_rhs_list))DEALLOCATE(self%bc_rhs_list)
IF(ASSOCIATED(self%olbp))DEALLOCATE(self%olbp)
!
IF(ASSOCIATED(self%rlimiter_nds))DEALLOCATE(self%rlimiter_nds)
IF(ASSOCIATED(self%limiter_pts))DEALLOCATE(self%limiter_pts)
IF(ASSOCIATED(self%cond_weights))DEALLOCATE(self%cond_weights)
IF(ASSOCIATED(self%bc_mat))DEALLOCATE(self%bc_mat)
IF(ASSOCIATED(self%bc_coil_mat))DEALLOCATE(self%bc_coil_mat)
IF(ASSOCIATED(self%rlcfs))DEALLOCATE(self%rlcfs)
IF(ASSOCIATED(self%isoflux_targets))DEALLOCATE(self%isoflux_targets)
IF(ASSOCIATED(self%saddle_targets))DEALLOCATE(self%saddle_targets)
IF(ASSOCIATED(self%coil_reg_mat))DEALLOCATE(self%coil_reg_mat)
IF(ASSOCIATED(self%coil_reg_targets))DEALLOCATE(self%coil_reg_targets)
IF(ASSOCIATED(self%coil_bounds))DEALLOCATE(self%coil_bounds)
IF(ASSOCIATED(self%saddle_cmask))DEALLOCATE(self%saddle_cmask)
IF(ASSOCIATED(self%saddle_pmask))DEALLOCATE(self%saddle_pmask)
IF(ASSOCIATED(self%Lcoils))DEALLOCATE(self%Lcoils)
!---
CALL self%lu_solver%delete()
CALL self%lu_solver_dt%delete()
!---
IF(self%ncoils_ext>0)THEN
  DO i=1,self%ncoils_ext
    IF(ASSOCIATED(self%coils_ext(i)%scale))DEALLOCATE(self%coils_ext(i)%scale)
    IF(ASSOCIATED(self%coils_ext(i)%pt))DEALLOCATE(self%coils_ext(i)%pt)
  END DO
  DEALLOCATE(self%coils_ext)
  self%ncoils_ext=0
END IF
!---
IF(self%ncoil_regs>0)THEN
  DO i=1,self%ncoil_regs
    IF(ASSOCIATED(self%coil_regions(i)%lc))DEALLOCATE(self%coil_regions(i)%lc)
  END DO
  DEALLOCATE(self%coil_regions)
  self%ncoil_regs=0
END IF
!---
IF(self%ncoils>0)THEN
  DO i=1,self%ncoils
    IF(ASSOCIATED(self%psi_coil(i)%f))CALL self%psi_coil(i)%f%delete()
  END DO
  DEALLOCATE(self%psi_coil,self%coil_currs,self%coil_vcont,self%coil_nturns)
  self%ncoils=0
END IF
!---
IF(self%ncond_regs>0)THEN
  DO i=1,self%ncond_regs
    IF(ASSOCIATED(self%cond_regions(i)%fixed))DEALLOCATE(self%cond_regions(i)%fixed)
    IF(ASSOCIATED(self%cond_regions(i)%mtype))DEALLOCATE(self%cond_regions(i)%mtype)
    IF(ASSOCIATED(self%cond_regions(i)%mind))DEALLOCATE(self%cond_regions(i)%mind)
    IF(ASSOCIATED(self%cond_regions(i)%eig_map))DEALLOCATE(self%cond_regions(i)%eig_map)
    IF(ASSOCIATED(self%cond_regions(i)%lc))DEALLOCATE(self%cond_regions(i)%lc)
    IF(ASSOCIATED(self%cond_regions(i)%weights))DEALLOCATE(self%cond_regions(i)%weights)
    IF(ASSOCIATED(self%cond_regions(i)%pair_signs))DEALLOCATE(self%cond_regions(i)%pair_signs)
    IF(ASSOCIATED(self%cond_regions(i)%fit_scales))DEALLOCATE(self%cond_regions(i)%fit_scales)
    IF(ASSOCIATED(self%cond_regions(i)%cond_vals))DEALLOCATE(self%cond_regions(i)%cond_vals)
    IF(ASSOCIATED(self%cond_regions(i)%rc))DEALLOCATE(self%cond_regions(i)%rc)
    IF(ASSOCIATED(self%cond_regions(i)%cond_curr))DEALLOCATE(self%cond_regions(i)%cond_curr)
    IF(ASSOCIATED(self%cond_regions(i)%corr_3d))DEALLOCATE(self%cond_regions(i)%corr_3d)
    IF(self%cond_regions(i)%neigs>0)THEN
      DO j=1,self%cond_regions(i)%neigs
        IF(ASSOCIATED(self%cond_regions(i)%psi_eig(j)%f))CALL self%cond_regions(i)%psi_eig(j)%f%delete()
      END DO
      DEALLOCATE(self%cond_regions(i)%psi_eig)
      self%cond_regions(i)%neigs=0
    END IF
  END DO
  DEALLOCATE(self%cond_regions)
  self%ncond_regs=0
END IF
!---
IF(ASSOCIATED(self%psi))THEN
  CALL self%psi%delete()
  DEALLOCATE(self%psi)
END IF
IF(ASSOCIATED(self%chi))THEN
  CALL self%chi%delete()
  DEALLOCATE(self%chi)
END IF
IF(ASSOCIATED(self%u_hom))THEN
  CALL self%u_hom%delete()
  DEALLOCATE(self%u_hom)
END IF
!---
CALL self%dels%delete()
CALL self%mrop%delete()
CALL self%mop%delete()
DEALLOCATE(self%dels,self%mrop,self%mop)
IF(ASSOCIATED(self%dels_dt))THEN
  CALL self%dels_dt%delete()
  DEALLOCATE(self%dels_dt)
END IF
! Destory I and P in the future
NULLIFY(self%I,self%P)
NULLIFY(self%set_eta)
end subroutine gs_destroy
!---------------------------------------------------------------------------
! SUBROUTINE compute_bcmat
!---------------------------------------------------------------------------
!> Needs Docs
!!
!! @param[in,out] self G-S object
!! @param[in,out] a Psi field
!! @param[in,out] b Source field
!---------------------------------------------------------------------------
subroutine compute_bcmat(self)
class(gs_eq), intent(inout) :: self
!---
integer(4) :: i,j,m,l,jr,jc,k,io_unit,nrhs,ierr,i_inds(1),j_inds(1)
integer(4), allocatable :: elist(:,:),marker(:),bemap(:),el1(:),el2(:),eflag(:)
integer(4), allocatable :: loc_map1(:),loc_map2(:)
real(8) :: f(3),pt(3),goptmp(3,3),gop(3),v,dl(2),dn(2),pt_int(3),pt2(3)
real(8), allocatable :: ltmp(:),vflux_mat(:,:)
real(8), allocatable :: rop1(:),rop2(:),gop1(:,:),gop2(:,:), massmat(:,:)
logical :: file_exists
integer(4) :: cell1,ed1,kk
real(8) :: pts1(2,2),pts2(2,2),dl1(2),dl2(2),dl1_mag,dl2_mag,val,f1(3),f2(3)
real(8) :: goptmp1(3,3),goptmp2(3,3),one_val(1,1),dn1(2),dn2(2),offset,grad_tmp(2)
logical, allocatable, dimension(:) :: vert_flag,edge_flag
type(oft_quad_type) :: quad,quad_hp,sing_quad
CLASS(oft_bmesh), POINTER :: smesh
!
integer(4), parameter :: qp_div_lim = 15
integer(4) :: neval,last,iwork(qp_div_lim),jc_active,nfail
real(8) :: abserr,work(5*qp_div_lim)
character(len=6) :: nfail_str
integer(4), save :: cell2,jc_int,ed2
real(8), save :: pt1(3)
!$omp threadprivate(pt1,cell2,jc_int,ed2)
IF(ASSOCIATED(self%bc_lmat))RETURN
WRITE(*,*)'Computing flux BC matrix '
CALL set_quad_1d(quad,self%fe_rep%order+2)
CALL set_quad_1d(quad_hp,4*self%fe_rep%order+2)
smesh=>self%fe_rep%mesh
!
sing_quad%dim=1
sing_quad%np=7
sing_quad%order=6
allocate(sing_quad%pts(1,sing_quad%np),sing_quad%wts(sing_quad%np))
sing_quad%pts(1,:)=[0.175965211846577428056264284949d-2,0.244696507125133674276453373497d-1, &
  0.106748056858788954180259781083d0,0.275807641295917383077859512057d0, &
  0.517855142151833716158668961982d0,0.771815485362384900274646869494d0, &
  0.952841340581090558994306588503d0]
sing_quad%wts=[0.663266631902570511783904989051d-2,0.457997079784753341255767348120d-1, &
  0.123840208071318194550489564922d0,0.212101926023811930107914875456d0, &
  0.261390645672007725646580606859d0,0.231636180290909384318815526104d0, &
  0.118598665644451726132783641957d0]
!
allocate(self%olbp(smesh%nbp+1))
CALL get_olbp(self%mesh,self%olbp)
!---Compute oriented edge list
ALLOCATE(elist(2,smesh%nbe))
DO i=1,smesh%nbe
  j=ABS(mesh_local_findedge(smesh,[self%olbp(i),self%olbp(i+1)]))
  elist(2,i)=smesh%lec(smesh%kec(j))
  DO m=1,3
    IF(j==ABS(smesh%lce(m,elist(2,i))))THEN
      elist(1,i)=m
      EXIT
    END IF
  END DO
  IF(m>3)CALL oft_abort("bad","",__FILE__)
END DO
!---Count rhs elements
allocate(marker(self%fe_rep%ne))
marker=0
ALLOCATE(el1(self%fe_rep%nce))
self%bc_nrhs=0
DO i=1,smesh%nbc
  j=smesh%lbc(i)
  call self%fe_rep%ncdofs(j,el1)
  DO m=1,self%fe_rep%nce
    IF(marker(el1(m))==0)THEN
      self%bc_nrhs=self%bc_nrhs+1
      marker(el1(m))=self%bc_nrhs
    END IF
  END DO
END DO
DEALLOCATE(el1)
!---
allocate(self%bc_rhs_list(self%bc_nrhs))
allocate(self%bc_lmat(self%fe_rep%nbe,self%fe_rep%nbe))
allocate(self%bc_bmat(self%fe_rep%nbe,self%bc_nrhs))
allocate(vflux_mat(self%bc_nrhs,self%fe_rep%nbe))
allocate(massmat(self%fe_rep%nbe,self%fe_rep%nbe))
self%bc_lmat=0.d0
self%bc_bmat=0.d0
vflux_mat=0.d0
massmat=0.d0
!
DO i=1,self%fe_rep%ne
  IF(marker(i)/=0)self%bc_rhs_list(marker(i))=i
END DO
allocate(bemap(self%fe_rep%ne))
bemap=0
DO i=1,self%fe_rep%nbe
  bemap(self%fe_rep%lbe(i))=i
END DO
!---Compute boundary current to volume flux projection matrix
nfail=0
!$omp parallel private(j,jr,jc,k,kk,rop1,gop1,loc_map1,cell1,el1,f1,ed1,dl1,dn1,dl1_mag,pts1,val, &
!$omp rop2,gop2,loc_map2,el2,f2,dl2,dn2,dl2_mag,pt2,pts2,work,neval,ierr,iwork,last,ltmp,goptmp1) &
!$omp reduction(+:nfail)
ALLOCATE(el1(self%fe_rep%nce),loc_map1(self%fe_rep%nce))
ALLOCATE(rop1(self%fe_rep%nce),gop1(3,self%fe_rep%nce))
ALLOCATE(el2(self%fe_rep%nce),loc_map2(self%fe_rep%nce))
ALLOCATE(rop2(self%fe_rep%nce),gop2(3,self%fe_rep%nce))
!$omp do schedule(static,10)
DO i=1,self%bc_nrhs
  IF(self%fe_rep%be(self%bc_rhs_list(i)))CYCLE
  cell1 = self%fe_rep%lec(self%fe_rep%kec(self%bc_rhs_list(i)))
  CALL self%fe_rep%ncdofs(cell1,el1)
  DO k=1,self%fe_rep%nce
    IF(el1(k)==self%bc_rhs_list(i))EXIT
  END DO
  CALL oft_blag_npos(self%fe_rep,cell1,k,f1)
  CALL oft_blag_eval(self%fe_rep,cell1,k,f1,rop1(1))
  pt1=smesh%log2phys(cell1,f1)
  DO j=1,smesh%nbe
    cell2=elist(2,j)
    ed2=elist(1,j)
    pts2(:,1)=smesh%r(1:2,smesh%lc(smesh%cell_ed(1,ed2),cell2))
    pts2(:,2)=smesh%r(1:2,smesh%lc(smesh%cell_ed(2,ed2),cell2))
    IF((pts2(1,1)<1.d-8).AND.(pts2(1,2)<1.d-8))CYCLE
    CALL self%fe_rep%ncdofs(cell2,el2)
    loc_map2=bemap(el2)
    !
    dl2=pts2(:,1)-pts2(:,2)
    dl2_mag=SQRT(SUM(dl2**2))
    DO jc=1,self%fe_rep%nce
      IF(loc_map2(jc)==0)CYCLE
      jc_int=jc
      CALL dqagse(integrand1,0.d0,1.d0,qp_int_tol,1.d2*qp_int_tol,qp_div_lim,val,abserr,neval,ierr, &
        work(1),work(qp_div_lim+1),work(2*qp_div_lim+1),work(3*qp_div_lim+1),iwork,last)
      ! ierr=-1
      IF(ierr/=0)THEN
        nfail=nfail+1
        val = 0.d0
        DO kk=1,quad_hp%np
          val=val + integrand1(quad_hp%pts(1,kk))*quad_hp%wts(kk)
        END DO
      END IF
      vflux_mat(i,loc_map2(jc))=vflux_mat(i,loc_map2(jc)) &
        + val*dl2_mag
    END DO
  END DO
END DO
!$omp end do nowait
!---Compute inductance and mass matrices for boundary
ALLOCATE(ltmp(self%fe_rep%nbe))
!$omp do schedule(static,10)
DO i=1,smesh%nbe
  cell1=elist(2,i)
  ed1=elist(1,i)
  pts1(:,1)=smesh%r(1:2,smesh%lc(smesh%cell_ed(1,ed1),cell1))
  pts1(:,2)=smesh%r(1:2,smesh%lc(smesh%cell_ed(2,ed1),cell1))
  IF((pts1(1,1)<1.d-8).AND.(pts1(1,2)<1.d-8))CYCLE
  dl1=pts1(:,1)-pts1(:,2)
  IF(self%olbp(i)==smesh%lc(smesh%cell_ed(2,ed1),cell1))dl1=-dl1
  dl1_mag=SQRT(SUM(dl1**2))
  dn1=[-dl1(2),dl1(1)]
  CALL self%fe_rep%ncdofs(cell1,el1)
  loc_map1=bemap(el1)
  DO k=1,quad%np
    f1 = 0.d0
    f1(smesh%cell_ed(1,ed1))=quad%pts(1,k)
    f1(smesh%cell_ed(2,ed1))=1.d0 - quad%pts(1,k)
    pt1=smesh%log2phys(cell1,f1)
    CALL smesh%jacobian(cell1,f1,goptmp1,v)
    DO jc=1,self%fe_rep%nce
      CALL oft_blag_eval(self%fe_rep,cell1,jc,f1,rop1(jc))
      CALL oft_blag_geval(self%fe_rep,cell1,jc,f1,gop1(:,jc),goptmp1)
    END DO
    !---
    DO jc=1,self%fe_rep%nce
      DO jr=1,self%fe_rep%nce
        IF(loc_map1(jr)==0)CYCLE
        !$omp atomic
        self%bc_bmat(loc_map1(jr),marker(el1(jc)))=self%bc_bmat(loc_map1(jr),marker(el1(jc))) &
          + rop1(jr)*DOT_PRODUCT(gop1(1:2,jc),dn1)*quad%wts(k)/pt1(1)
      END DO
    END DO
    !---Mass matrix
    DO jc=1,self%fe_rep%nce
      IF(loc_map1(jc)==0)CYCLE
      DO jr=1,self%fe_rep%nce
        IF(loc_map1(jr)==0)CYCLE
        !$omp atomic
        massmat(loc_map1(jr),loc_map1(jc))=massmat(loc_map1(jr),loc_map1(jc)) &
          + rop1(jr)*rop1(jc)*quad%wts(k)*dl1_mag
      END DO
    END DO
  END DO
  DO k=1,quad_hp%np
    f1 = 0.d0
    f1(smesh%cell_ed(1,ed1))=quad_hp%pts(1,k)
    f1(smesh%cell_ed(2,ed1))=1.d0 - quad_hp%pts(1,k)
    pt1=smesh%log2phys(cell1,f1)
    DO jc=1,self%fe_rep%nce
      CALL oft_blag_eval(self%fe_rep,cell1,jc,f1,rop1(jc))
    END DO
    !
    ltmp=0.d0
    DO j=1,smesh%nbe
      cell2=elist(2,j)
      ed2=elist(1,j)
      pts2(:,1)=smesh%r(1:2,smesh%lc(smesh%cell_ed(1,ed2),cell2))
      pts2(:,2)=smesh%r(1:2,smesh%lc(smesh%cell_ed(2,ed2),cell2))
      IF((pts2(1,1)<1.d-8).AND.(pts2(1,2)<1.d-8))CYCLE
      CALL self%fe_rep%ncdofs(cell2,el2)
      loc_map2=bemap(el2)
      IF(j==i)THEN
        ! Segment 1
        dl2=pt1(1:2)-pts2(:,2)
        dl2_mag=SQRT(SUM(dl2**2))
        dn2=[-dl2(2),dl2(1)]
        IF(self%olbp(j)==smesh%lc(smesh%cell_ed(2,ed2),cell2))dn2=-dn2
        DO jc=1,self%fe_rep%nce
          IF(loc_map2(jc)==0)CYCLE
          DO kk=1,sing_quad%np
            f2 = f1*(1.d0 - sing_quad%pts(1,kk))
            f2(smesh%cell_ed(2,ed2))=f2(smesh%cell_ed(2,ed2))+sing_quad%pts(1,kk)
            CALL oft_blag_eval(self%fe_rep,cell2,jc,f2,rop2(jc))
            pt2=smesh%log2phys(cell2,f2)
            ! Lmat \int phi^T \int phi * green * dl2 * dl1
            val=green(pt2(1),pt2(2),pt1(1),pt1(2))
            ltmp(loc_map2(jc))=ltmp(loc_map2(jc)) + rop2(jc)*val*sing_quad%wts(kk)*dl2_mag
          END DO
        END DO
        ! Segment 2
        dl2=pts2(:,1)-pt1(1:2)
        dl2_mag=SQRT(SUM(dl2**2))
        dn2=[-dl2(2),dl2(1)]
        IF(self%olbp(j)==smesh%lc(smesh%cell_ed(2,ed2),cell2))dn2=-dn2
        DO jc=1,self%fe_rep%nce
          IF(loc_map2(jc)==0)CYCLE
          DO kk=1,sing_quad%np
            f2 = f1*(1.d0 - sing_quad%pts(1,kk))
            f2(smesh%cell_ed(1,ed2))=f2(smesh%cell_ed(1,ed2))+sing_quad%pts(1,kk)
            CALL oft_blag_eval(self%fe_rep,cell2,jc,f2,rop2(jc))
            pt2=smesh%log2phys(cell2,f2)
            ! Lmat \int phi^T \int phi * green * dl2 * dl1
            val=green(pt2(1),pt2(2),pt1(1),pt1(2))
            ltmp(loc_map2(jc))=ltmp(loc_map2(jc)) + rop2(jc)*val*sing_quad%wts(kk)*dl2_mag
          END DO
        END DO
      ELSE
        dl2=pts2(:,1)-pts2(:,2)
        dl2_mag=SQRT(SUM(dl2**2))
        dn2=[-dl2(2),dl2(1)]
        IF(self%olbp(j)==smesh%lc(smesh%cell_ed(2,ed2),cell2))dn2=-dn2
        DO jc=1,self%fe_rep%nce
          IF(loc_map2(jc)==0)CYCLE
          jc_int=jc
          CALL dqagse(integrand2,0.d0,1.d0,qp_int_tol,1.d2*qp_int_tol,qp_div_lim,val,abserr,neval,ierr, &
            work(1),work(qp_div_lim+1),work(2*qp_div_lim+1),work(3*qp_div_lim+1),iwork,last)
          ! ierr=-1
          IF(ierr/=0)THEN
            nfail=nfail+1
            val=0.d0
            DO kk=1,quad_hp%np
              val=val + integrand2(quad_hp%pts(1,kk))*quad_hp%wts(kk)
            END DO
          END IF
          ltmp(loc_map2(jc))=ltmp(loc_map2(jc)) + val*dl2_mag
        END DO
      END IF
    END DO
    !---Accumulate to matrices
    DO jc=1,self%fe_rep%nce
      IF(loc_map1(jc)==0)CYCLE
      DO jr=1,self%fe_rep%nbe
        !$omp atomic
        self%bc_lmat(jr,loc_map1(jc))=self%bc_lmat(jr,loc_map1(jc)) &
          + ltmp(jr)*rop1(jc)*quad_hp%wts(k)*dl1_mag
      END DO
    END DO
  END DO
END DO
DEALLOCATE(el1,el2,rop1,rop2,gop1,gop2,loc_map1,loc_map2,ltmp)
!$omp end parallel
IF(nfail>0)THEN
  WRITE(nfail_str,'(I6)')nfail
  CALL oft_warn("QUADPACK integration failed "//TRIM(nfail_str)//" times in vacuum BC calculation.")
END IF
!
ALLOCATE(vert_flag(smesh%np),edge_flag(smesh%ne),self%fe_flag(self%fe_rep%ne))
vert_flag=smesh%r(1,:)<1.d-8
edge_flag=.FALSE.
DO i=1,smesh%ne
  IF(vert_flag(smesh%le(1,i)).AND.vert_flag(smesh%le(2,i)))edge_flag(i)=.TRUE.
END DO
CALL bfem_map_flag(self%fe_rep,vert_flag,edge_flag,self%fe_flag)
DEALLOCATE(vert_flag,edge_flag)
!
DO i=1,self%fe_rep%nbe
  IF(self%fe_flag(self%fe_rep%lbe(i)))THEN
    self%bc_lmat(i,:)=0.d0; self%bc_lmat(:,i)=0.d0; self%bc_lmat(i,i)=1.d0
    massmat(i,:)=0.d0; massmat(:,i)=0.d0; massmat(i,i)=1.d0
  END IF
END DO
!
CALL lapack_matinv(self%fe_rep%nbe,self%bc_lmat,ierr)
IF(ierr>0)WRITE(*,*)'ERR1',ierr,smesh%nbp
!
vflux_mat = MATMUL(vflux_mat,MATMUL(self%bc_lmat,massmat))
self%bc_lmat = MATMUL(massmat,MATMUL(self%bc_lmat,massmat))
DO i=1,self%bc_nrhs
  IF(self%fe_rep%be(self%bc_rhs_list(i)))THEN
    vflux_mat(i,:)=0.d0; vflux_mat(i,bemap(self%bc_rhs_list(i)))=1.d0
  END IF
END DO
!--- Final BC is I(psi_b) = B_t(psi) - B_t(psi_b) [B_t with B_n = 0]
self%bc_lmat=self%bc_lmat - MATMUL(self%bc_bmat,vflux_mat)

DO i=1,self%fe_rep%nbe
  IF(self%fe_flag(self%fe_rep%lbe(i)))THEN
    self%bc_lmat(i,:)=0.d0
    self%bc_bmat(i,:)=0.d0
  END IF
END DO
CALL quad%delete()
CALL sing_quad%delete()
DEALLOCATE(massmat,marker,bemap,vflux_mat)
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Complete'
contains
function integrand1(x) result(itegrand)
real(8), intent(in) :: x
real(8) :: itegrand,pt2(3),rop2(1),val,f2(3)
f2 = 0.d0
f2(smesh%cell_ed(1,ed2))=x
f2(smesh%cell_ed(2,ed2))=1.d0 - x
CALL oft_blag_eval(self%fe_rep,cell2,jc_int,f2,rop2(1))
pt2=smesh%log2phys(cell2,f2)
! Lmat \int phi^T \int phi * green * dl2 * dl1
val=green(pt2(1),pt2(2),pt1(1),pt1(2))
itegrand=rop2(1)*val
end function integrand1
!
function integrand2(x) result(itegrand)
real(8), intent(in) :: x
real(8) :: itegrand,pt2(3),rop2(1),val,f2(3)
f2 = 0.d0
f2(smesh%cell_ed(1,ed2))=x
f2(smesh%cell_ed(2,ed2))=1.d0 - x
CALL oft_blag_eval(self%fe_rep,cell2,jc_int,f2,rop2(1))
pt2=smesh%log2phys(cell2,f2)
! Lmat \int phi^T \int phi * green * dl2 * dl1
val=green(pt1(1),pt1(2),pt2(1),pt2(2))
itegrand=rop2(1)*val
end function integrand2
end subroutine compute_bcmat
!
subroutine get_olbp(smesh,olbp)
class(oft_bmesh), intent(inout) :: smesh
integer(4), intent(out) :: olbp(:)
integer(4) :: i,ii,j,k,l
real(8) :: val,orient(2)
olbp=0
j=0
val=1.d99
DO i=1,smesh%nbp
  IF(ABS(smesh%r(1,smesh%lbp(i)))+ABS(smesh%r(2,smesh%lbp(i)))<val)THEN
    val=ABS(smesh%r(1,smesh%lbp(i)))+ABS(smesh%r(2,smesh%lbp(i)))
    j=smesh%lbp(i)
  END IF
END DO
olbp(1)=j
olbp(smesh%nbp+1)=j
DO j=smesh%kpc(olbp(1)),smesh%kpc(olbp(1)+1)-1
  i=smesh%lpc(j)
  IF(.NOT.smesh%bc(i))CYCLE
  DO l=1,3
    k=ABS(smesh%lce(l,i))
    IF(.NOT.smesh%be(k))CYCLE
    IF(smesh%le(1,k)==olbp(1))THEN
      olbp(2)=smesh%le(2,k)
      EXIT
    END IF
    IF(smesh%le(2,k)==olbp(1))THEN
      olbp(2)=smesh%le(1,k)
      EXIT
    END IF
  END DO
END DO
!---
DO ii=3,smesh%nbp
  DO j=smesh%kpc(olbp(ii-1)),smesh%kpc(olbp(ii-1)+1)-1
    i=smesh%lpc(j)
    IF(.NOT.smesh%bc(i))CYCLE
    DO l=1,3
      k=ABS(smesh%lce(l,i))
      IF(.NOT.smesh%be(k))CYCLE
      IF(ANY(smesh%le(:,k)==olbp(ii-1)).AND.ALL(smesh%le(:,k)/=olbp(ii-2)))THEN
        IF(smesh%le(1,k)==olbp(ii-1))THEN
          olbp(ii)=smesh%le(2,k)
        ELSE
          olbp(ii)=smesh%le(1,k)
        END IF
      END IF
    END DO
  END DO
END DO
IF(MINVAL(olbp)==0)CALL oft_abort('Zero value detected in OLBP','get_olbp',__FILE__)
!---Check orientation
orient=0.d0
DO i=1,smesh%nbp
    IF(i<smesh%nbp/2)THEN
      orient(1)=orient(1)+smesh%r(2,olbp(i))
    ELSE
      orient(2)=orient(2)+smesh%r(2,olbp(i))
    END IF
END DO
!---Reverse list if necessary
IF(orient(2)>orient(1))THEN
  DO i=2,smesh%nbp/2+1
    k=olbp(i)
    olbp(i)=olbp(smesh%nbp+2-i)
    olbp(smesh%nbp+2-i)=k
  END DO
END IF
end subroutine get_olbp
!---------------------------------------------------------------------------
! SUBROUTINE set_bcmat
!---------------------------------------------------------------------------
!> Needs Docs
!!
!! @param[in,out] self G-S object
!! @param[in,out] a Psi field
!! @param[in,out] b Source field
!---------------------------------------------------------------------------
subroutine set_bcmat(self,mat)
class(gs_eq), intent(inout) :: self
class(oft_matrix), intent(inout) :: mat
integer(4) :: i,j,i_inds(1),j_inds(1)
real(8) :: one_val(1,1)
!---Add to matrix
! | A_ii A_ib |
! | M_bi M_bb + M*L^-1*M |
one_val=1.d0
DO i=1,self%fe_rep%nbe
  i_inds=self%fe_rep%lbe(i)
  IF(self%fe_flag(self%fe_rep%lbe(i)))THEN
    CALL mat%add_values(i_inds,i_inds,one_val,1,1)
  ELSE
    DO j=1,self%bc_nrhs
      j_inds=self%bc_rhs_list(j)
      IF(ABS(self%bc_bmat(i,j))<1.d-20)CYCLE
      CALL mat%add_values(i_inds,j_inds,self%bc_bmat(i:i,j:j),1,1)
    END DO
    DO j=1,self%fe_rep%nbe
      j_inds=self%fe_rep%lbe(j)
      CALL mat%add_values(i_inds,j_inds,self%bc_lmat(i:i,j:j),1,1)
    END DO
  END IF
END DO
end subroutine set_bcmat
!
subroutine gs_mat_create(fe_rep,new)
class(oft_scalar_bfem), intent(inout) :: fe_rep
CLASS(oft_matrix), POINTER, INTENT(out) :: new
CLASS(oft_vector), POINTER :: tmp_vec
INTEGER(4) :: i,j,nr
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: ltmp
TYPE(oft_graph_ptr) :: graphs(1,1)
DEBUG_STACK_PUSH
CALL fe_rep%vec_create(tmp_vec)
ALLOCATE(graphs(1,1)%g)
graphs(1,1)%g%nr=fe_rep%ne
graphs(1,1)%g%nrg=fe_rep%global%ne
graphs(1,1)%g%nc=fe_rep%ne
graphs(1,1)%g%ncg=fe_rep%global%ne
graphs(1,1)%g%nnz=0
ALLOCATE(graphs(1,1)%g%kr(fe_rep%ne+1))
graphs(1,1)%g%kr=0
ALLOCATE(ltmp(fe_rep%ne))
DO i=1,fe_rep%ne
  IF(fe_rep%be(i))THEN
    ltmp=2*fe_rep%ne
    ltmp(1:fe_rep%nbe)=fe_rep%lbe
    DO j=fe_rep%kee(i),fe_rep%kee(i+1)-1
      ltmp(j-fe_rep%kee(i)+fe_rep%nbe+1)=fe_rep%lee(j)
    END DO
    CALL sort_array(ltmp,fe_rep%nbe+fe_rep%kee(i+1)-fe_rep%kee(i))
    graphs(1,1)%g%kr(i)=1
    DO j=2,fe_rep%nbe+fe_rep%kee(i+1)-fe_rep%kee(i)
      IF(ltmp(j)>ltmp(j-1))graphs(1,1)%g%kr(i)=graphs(1,1)%g%kr(i)+1
    END DO
    ! graphs(1,1)%g%kr(i)=fe_rep%ne
  ELSE
    graphs(1,1)%g%kr(i)=fe_rep%kee(i+1)-fe_rep%kee(i)
  END IF
END DO
graphs(1,1)%g%nnz=SUM(graphs(1,1)%g%kr)
graphs(1,1)%g%kr(fe_rep%ne+1)=graphs(1,1)%g%nnz+1
DO i=fe_rep%ne,1,-1
  graphs(1,1)%g%kr(i)=graphs(1,1)%g%kr(i+1)-graphs(1,1)%g%kr(i)
END DO
IF(graphs(1,1)%g%kr(1)/=1)CALL oft_abort('Bad new graph setup','gs_mat_create', &
__FILE__)
ALLOCATE(graphs(1,1)%g%lc(graphs(1,1)%g%nnz))
DO i=1,fe_rep%ne
  IF(fe_rep%be(i))THEN
    ltmp=2*fe_rep%ne
    ltmp(1:fe_rep%nbe)=fe_rep%lbe
    DO j=fe_rep%kee(i),fe_rep%kee(i+1)-1
      ltmp(j-fe_rep%kee(i)+fe_rep%nbe+1)=fe_rep%lee(j)
    END DO
    CALL sort_array(ltmp,fe_rep%nbe+fe_rep%kee(i+1)-fe_rep%kee(i))
    nr=0
    graphs(1,1)%g%lc(graphs(1,1)%g%kr(i))=ltmp(1)
    DO j=2,fe_rep%nbe+fe_rep%kee(i+1)-fe_rep%kee(i)
      IF(ltmp(j)>ltmp(j-1))THEN
        nr=nr+1
        graphs(1,1)%g%lc(graphs(1,1)%g%kr(i)+nr)=ltmp(j)
      END IF
    END DO
    ! DO j=1,fe_rep%ne
    !   graphs(1,1)%g%lc(graphs(1,1)%g%kr(i)+j-1)=j
    ! END DO
  ELSE
    DO j=fe_rep%kee(i),fe_rep%kee(i+1)-1
      graphs(1,1)%g%lc(graphs(1,1)%g%kr(i)+j-fe_rep%kee(i))=fe_rep%lee(j)
    END DO
  END IF
END DO
!---
CALL create_matrix(new,graphs,tmp_vec,tmp_vec)
CALL tmp_vec%delete
DEALLOCATE(ltmp)
DEALLOCATE(graphs(1,1)%g,tmp_vec)
DEBUG_STACK_POP
end subroutine gs_mat_create
end module oft_gs
