!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!> @file hcurl_grad_operators.F90
!
!> FE operators for full H(Curl) + Grad(H^1) space
!! - Operator construction
!!   - MOP: mass matrix \f$ \int \left( u^T \cdot v \right) dV \f$
!! - Interpolation classes
!! - Explicit operators
!!   - Gradient
!!   - Divergence
!!   - GradT
!!   - MC gradient
!!   - Divergence cleaning
!! - Boundary conditions
!! - Multi-Grid setup and operators
!! - Default preconditioner setup
!!   - Optimal smoother calculation
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_hcurl
!------------------------------------------------------------------------------
MODULE oft_hcurl_grad_operators
USE oft_base
USE oft_sort, ONLY: sort_array
USE oft_quadrature
USE oft_mesh_type, ONLY: oft_mesh, oft_bmesh, cell_is_curved
USE oft_trimesh_type, ONLY: oft_trimesh
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, oft_matrix, oft_matrix_ptr, &
  oft_graph_ptr, oft_graph
USE oft_solver_base, ONLY: oft_solver, oft_solver_bc
USE oft_deriv_matrices, ONLY: create_diagmatrix
#ifdef HAVE_ARPACK
USE oft_arpack, ONLY: oft_irlm_eigsolver
#endif
USE oft_la_utils, ONLY: create_matrix, combine_matrices
USE oft_solver_utils, ONLY: create_mlpre, create_cg_solver, create_diag_pre, &
  create_native_pre
USE fem_base, ONLY: oft_afem_type, oft_fem_type, fem_max_levels, oft_ml_fem_type
USE fem_utils, ONLY: fem_interp
USE fem_composite, ONLY: oft_ml_fem_comp_type, oft_fem_comp_type, oft_ml_fe_comp_vecspace
USE oft_lag_basis, ONLY: oft_lag_geval_all
USE oft_h1_basis, ONLY: oft_h1_geval_all, oft_h1_d2eval, oft_h1_fem, &
  oft_3D_h1_cast
USE oft_h1_operators, ONLY: oft_h1_zerob, oft_h1_zerogrnd, oft_h1_getlop, oft_h1_gop
USE oft_hcurl_basis, ONLY: oft_hcurl_eval_all, &
  oft_hcurl_ceval_all, oft_hcurl_get_cgops, oft_hcurl_fem, oft_3D_hcurl_cast
USE oft_hcurl_operators, ONLY: oft_hcurl_rop, oft_hcurl_cop
IMPLICIT NONE
#include "local.h"
!------------------------------------------------------------------------------
!> Interpolate a H(Curl) + Grad(H^1) vector field
!------------------------------------------------------------------------------
type, extends(fem_interp) :: oft_hcurl_grad_rinterp
  class(oft_vector), pointer :: u => NULL() !< Field to interpolate
  integer(i4), pointer, dimension(:) :: cache_cell => NULL() !< Needs Docs
  real(r8), pointer, dimension(:) :: grad_vals => NULL() !< Local gradient values
  real(r8), pointer, dimension(:) :: curl_vals => NULL() !< Local curl values
  real(r8), pointer, dimension(:,:) :: cache_grad => NULL() !< Needs Docs
  real(r8), pointer, dimension(:,:) :: cache_curl => NULL() !< Needs Docs
  class(oft_h1_fem), pointer :: hgrad_rep => NULL() !< Grad(H^1) FE representation
  class(oft_hcurl_fem), pointer :: hcurl_rep => NULL() !< H(Curl) FE representation
contains
  !> Retrieve local values for interpolation
  generic :: setup => setup1, setup2
  procedure :: setup1 => rinterp_setup1
  procedure :: setup2 => rinterp_setup2
  !> Reconstruct field
  procedure :: interp => rinterp_apply
  !> Destroy temporary internal storage
  procedure :: delete => rinterp_delete
end type oft_hcurl_grad_rinterp
!------------------------------------------------------------------------------
!> Interpolate \f$ \nabla \times \f$ of a H(Curl) + Grad(H^1) vector field
!------------------------------------------------------------------------------
type, extends(oft_hcurl_grad_rinterp) :: oft_hcurl_grad_cinterp
contains
  !> Reconstruct field
  procedure :: interp => cinterp_apply
end type oft_hcurl_grad_cinterp
!------------------------------------------------------------------------------
!> Interpolate \f$ \nabla \times \f$ of a H(Curl) + Grad(H^1) vector field
!------------------------------------------------------------------------------
type, extends(oft_hcurl_grad_rinterp) :: oft_hcurl_grad_dinterp
contains
  !> Reconstruct field
  procedure :: interp => dinterp_apply
end type oft_hcurl_grad_dinterp
!------------------------------------------------------------------------------
!> Clean the divergence from a H(Curl) + Grad(H^1) vector field
!!
!! Divergence is removed by adding a gradient field, such that \f$ f = f +
!! \nabla \phi \f$, where \f$ \nabla^2 \phi = - \nabla \cdot f \f$. Cleaning
!! may also be applied to a field which is pre-multiplied by the mass matrix
!! in which case \f$ \nabla^2 \phi = - \nabla^T f \f$ and \f$ f = f + M \nabla
!! \phi \f$. The mass matrix version is applied if `mop` is associated with the
!! corresponding mass matrix.
!------------------------------------------------------------------------------
type, extends(oft_solver_bc) :: oft_hcurl_grad_divout
  integer(i4) :: count = 0 !< Number of times apply has been called
  integer(i4) :: app_freq = 1 !< Frequency to apply solver
  logical :: keep_boundary = .FALSE. !< Flag for keeping boundary gradients
  logical :: pm = .FALSE. !< Flag for solver convergence monitor
  logical :: internal_solver = .FALSE. !< Solver was constructed internally?
  class(oft_solver), pointer :: solver => NULL() !< Solver object for H^1::LOP operator
  class(oft_solver_bc), pointer :: bc => NULL() !< Boundary condition
  class(oft_matrix), pointer :: mop => NULL() !< Mass matrix, applies divoutm if associated
  class(oft_vector), pointer :: bnorm => NULL() !< Normal field source on boundary
  class(oft_ml_fem_comp_type), pointer :: ML_hcurl_full => NULL()
  class(oft_ml_fem_type), pointer :: ML_grad => NULL()
  class(oft_ml_fem_type), pointer :: ML_curl => NULL()
contains
  !> Setup matrix and solver with default
  procedure :: setup => divout_setup
  !> Clean divergence from field
  procedure :: apply => divout_apply
  !> Clean-up internal storage
  procedure :: delete => divout_delete
end type oft_hcurl_grad_divout
!------------------------------------------------------------------------------
!> Orthogonalize a H(Curl) + Grad(H^1) vector field by zeroing the gradient subspace
!------------------------------------------------------------------------------
type, extends(oft_solver_bc) :: oft_hcurl_grad_gzero
contains
  !> Perform orthoganlization
  procedure :: apply => zerograd_apply
  !> Clean-up internal variables
  procedure :: delete => zerograd_delete
end type oft_hcurl_grad_gzero
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(oft_solver_bc) :: oft_hcurl_grad_zerob
  class(oft_ml_fem_comp_type), pointer :: ML_hcurl_grad_rep => NULL() !< FE representation
contains
  procedure :: apply => zerob_apply
  procedure :: delete => zerob_delete
end type oft_hcurl_grad_zerob
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(oft_hcurl_grad_zerob) :: oft_hcurl_grad_czerob
contains
  procedure :: apply => curl_zerob_apply
end type oft_hcurl_grad_czerob
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(oft_hcurl_grad_zerob) :: oft_hcurl_grad_gzerop
contains
  procedure :: apply => grad_zerop_apply
end type oft_hcurl_grad_gzerop
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(oft_hcurl_grad_zerob) :: oft_hcurl_grad_zeroi
contains
  procedure :: apply => zeroi_apply
end type oft_hcurl_grad_zeroi
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, PUBLIC, extends(oft_ml_fe_comp_vecspace) :: oft_ml_hcurl_grad_vecspace
contains
  !> Needs docs
  PROCEDURE :: interp => ml_vecspace_interp
  !> Needs docs
  PROCEDURE :: inject => ml_vecspace_inject
end type oft_ml_hcurl_grad_vecspace
!---Pre options
integer(i4), private :: nu_mop(fem_max_levels)=0 !< Needs Docs
real(r8), private :: df_mop(fem_max_levels)=-1.d99 !< Needs Docs
contains
!------------------------------------------------------------------------------
!> Read-in options for the basic H(Curl) + Grad(H^1) space ML preconditioners
!------------------------------------------------------------------------------
subroutine hcurl_grad_mloptions
integer(i4) :: ierr,io_unit
namelist/hcurl_grad_op_options/df_mop,nu_mop
DEBUG_STACK_PUSH
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,hcurl_grad_op_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr>0)THEN
  CALL oft_abort('Error reading "hcurl_grad_op_options" group in input file','hcurl_grad_mloptions',__FILE__)
END IF
IF(df_mop(1)<-1.d90)THEN
  IF(oft_env%head_proc)THEN
    CALL oft_warn("No H(Curl) + Grad(H^1) MG smoother settings found:")
    CALL oft_warn("  Using default values, which may result in convergence failure.")
  END IF
  nu_mop=2
  df_mop=.2d0
END IF
DEBUG_STACK_POP
end subroutine hcurl_grad_mloptions
!------------------------------------------------------------------------------
!> Setup interpolator for H(Curl) + Grad(H^1) vector fields
!!
!! Fetches local representation used for interpolation from vector object
!------------------------------------------------------------------------------
subroutine rinterp_setup1(self,hcurl_grad_rep)
class(oft_hcurl_grad_rinterp), intent(inout) :: self
class(oft_fem_comp_type), target, intent(inout) :: hcurl_grad_rep
IF(ASSOCIATED(self%parent))THEN
  SELECT TYPE(this=>self%parent)
    CLASS IS(oft_hcurl_grad_rinterp)
      self%curl_vals=>this%curl_vals
      self%grad_vals=>this%grad_vals
      self%hgrad_rep=>this%hgrad_rep
      self%hcurl_rep=>this%hcurl_rep
      self%mesh=>this%mesh
      self%cache_cell=>this%cache_cell
      self%cache_grad=>this%cache_grad
      self%cache_curl=>this%cache_curl
    CLASS DEFAULT
      CALL oft_abort('Parent interpolator must be of same class','rinterp_setup',__FILE__)
  END SELECT
ELSE
  !---Get local slice
  CALL self%u%get_local(self%curl_vals,1)
  CALL self%u%get_local(self%grad_vals,2)
  SELECT TYPE(this=>hcurl_grad_rep%fields(1)%fe)
    CLASS IS(oft_hcurl_fem)
      self%hcurl_rep=>this
      self%mesh=>this%mesh
    CLASS DEFAULT
      CALL oft_abort("Invalid HCurl space","rinterp_setup",__FILE__)
  END SELECT
  SELECT TYPE(this=>hcurl_grad_rep%fields(2)%fe)
    CLASS IS(oft_h1_fem)
      self%hgrad_rep=>this
    CLASS DEFAULT
      CALL oft_abort("Invalid HGrad space","rinterp_setup",__FILE__)
  END SELECT
  IF(.NOT.ASSOCIATED(self%cache_cell))THEN
    ALLOCATE(self%cache_cell(0:oft_env%nthreads-1))
    ALLOCATE(self%cache_grad(self%hgrad_rep%nce,0:oft_env%nthreads-1))
    ALLOCATE(self%cache_curl(self%hcurl_rep%nce,0:oft_env%nthreads-1))
  END IF
  self%cache_cell=-1
  self%cache_grad=0.d0
  self%cache_curl=0.d0
END IF
end subroutine rinterp_setup1
!------------------------------------------------------------------------------
!> Setup interpolator for H(Curl) + Grad(H^1) vector fields
!!
!! Fetches local representation used for interpolation from vector object
!------------------------------------------------------------------------------
subroutine rinterp_setup2(self,hcurl_rep,hgrad_rep)
class(oft_hcurl_grad_rinterp), intent(inout) :: self
class(oft_afem_type), target, intent(inout) :: hcurl_rep
class(oft_afem_type), target, intent(inout) :: hgrad_rep
IF(ASSOCIATED(self%parent))THEN
  SELECT TYPE(this=>self%parent)
    CLASS IS(oft_hcurl_grad_rinterp)
      self%curl_vals=>this%curl_vals
      self%grad_vals=>this%grad_vals
      self%hgrad_rep=>this%hgrad_rep
      self%hcurl_rep=>this%hcurl_rep
      self%mesh=>this%mesh
      self%cache_cell=>this%cache_cell
      self%cache_grad=>this%cache_grad
      self%cache_curl=>this%cache_curl
    CLASS DEFAULT
      CALL oft_abort('Parent interpolator must be of same class','rinterp_setup',__FILE__)
  END SELECT
ELSE
  !---Get local slice
  CALL self%u%get_local(self%curl_vals,1)
  CALL self%u%get_local(self%grad_vals,2)
  SELECT TYPE(hcurl_rep)
    CLASS IS(oft_hcurl_fem)
      self%hcurl_rep=>hcurl_rep
      self%mesh=>hcurl_rep%mesh
    CLASS DEFAULT
      CALL oft_abort("Invalid HCurl space","rinterp_setup",__FILE__)
  END SELECT
  SELECT TYPE(hgrad_rep)
    CLASS IS(oft_h1_fem)
      self%hgrad_rep=>hgrad_rep
    CLASS DEFAULT
      CALL oft_abort("Invalid HGrad space","rinterp_setup",__FILE__)
  END SELECT
  IF(.NOT.ASSOCIATED(self%cache_cell))THEN
    ALLOCATE(self%cache_cell(0:oft_env%nthreads-1))
    ALLOCATE(self%cache_grad(self%hgrad_rep%nce,0:oft_env%nthreads-1))
    ALLOCATE(self%cache_curl(self%hcurl_rep%nce,0:oft_env%nthreads-1))
  END IF
  self%cache_cell=-1
  self%cache_grad=0.d0
  self%cache_curl=0.d0
END IF
end subroutine rinterp_setup2
!------------------------------------------------------------------------------
!> Destroy temporary internal storage
!------------------------------------------------------------------------------
subroutine rinterp_delete(self)
class(oft_hcurl_grad_rinterp), intent(inout) :: self
IF(ASSOCIATED(self%parent))THEN
  NULLIFY(self%curl_vals,self%grad_vals)
  NULLIFY(self%cache_cell,self%cache_grad,self%cache_curl)
  NULLIFY(self%hcurl_rep,self%hgrad_rep,self%u)
  NULLIFY(self%parent)
ELSE
!---Deallocate local storage
IF(ASSOCIATED(self%curl_vals))THEN
  DEALLOCATE(self%curl_vals,self%grad_vals)
  DEALLOCATE(self%cache_cell,self%cache_grad,self%cache_curl)
END IF
NULLIFY(self%hcurl_rep,self%hgrad_rep,self%u)
END IF
end subroutine rinterp_delete
!------------------------------------------------------------------------------
!> Reconstruct a H(Curl) + Grad(H^1) vector field
!------------------------------------------------------------------------------
subroutine rinterp_apply(self,cell,f,gop,val)
class(oft_hcurl_grad_rinterp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [4]
real(r8), intent(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
real(r8), intent(out) :: val(:) !< Reconstructed field at f [3]
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: rop(3)
real(r8), allocatable :: hgrad_rop(:,:)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%grad_vals))CALL oft_abort('Setup has not been called!','rinterp_apply',__FILE__)
!---Pull cache
IF(self%cache_cell(oft_tid)/=cell)THEN
  allocate(j(self%hcurl_rep%nce))
  call self%hcurl_rep%ncdofs(cell,j)
  do jc=1,self%hcurl_rep%nce
    self%cache_curl(jc,oft_tid)=self%curl_vals(j(jc))
  end do
  deallocate(j)
  allocate(j(self%hgrad_rep%nce))
  call self%hgrad_rep%ncdofs(cell,j)
  do jc=1,self%hgrad_rep%nce
    self%cache_grad(jc,oft_tid)=self%grad_vals(j(jc))
  end do
  deallocate(j)
  self%cache_cell(oft_tid)=cell
END IF
!---Reconstruct field
val=0.d0
IF(ASSOCIATED(oft_hcurl_rop))THEN
  do jc=1,self%hcurl_rep%nce
    val=val+self%cache_curl(jc,oft_tid)*oft_hcurl_rop(:,jc)
  end do
ELSE
  ALLOCATE(hgrad_rop(3,self%hcurl_rep%nce))
  CALL oft_hcurl_eval_all(self%hcurl_rep,cell,f,hgrad_rop,gop)
  do jc=1,self%hcurl_rep%nce
    val=val+self%cache_curl(jc,oft_tid)*hgrad_rop(:,jc)
  end do
  DEALLOCATE(hgrad_rop)
END IF
IF(ASSOCIATED(oft_h1_gop))THEN
  do jc=1,self%hgrad_rep%nce
    val=val+self%cache_grad(jc,oft_tid)*oft_h1_gop(:,jc)
  end do
ELSE
  ALLOCATE(hgrad_rop(3,self%hgrad_rep%nce))
  CALL oft_h1_geval_all(self%hgrad_rep,cell,f,hgrad_rop,gop)
  do jc=1,self%hgrad_rep%nce
    val=val+self%cache_grad(jc,oft_tid)*hgrad_rop(:,jc)
  end do
  DEALLOCATE(hgrad_rop)
END IF
DEBUG_STACK_POP
end subroutine rinterp_apply
!------------------------------------------------------------------------------
!> Reconstruct \f$ \nabla \times \f$ of a H(Curl) + Grad(H^1) vector field
!------------------------------------------------------------------------------
subroutine cinterp_apply(self,cell,f,gop,val)
class(oft_hcurl_grad_cinterp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [4]
real(r8), intent(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
real(r8), intent(out) :: val(:) !< Reconstructed field at f [3]
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: cop(3),cgop(3,6)
real(r8), allocatable :: hcurl_cop(:,:)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%grad_vals))CALL oft_abort('Setup has not been called!','cinterp_apply',__FILE__)
!---Pull cache
IF(self%cache_cell(oft_tid)/=cell)THEN
  allocate(j(self%hcurl_rep%nce))
  call self%hcurl_rep%ncdofs(cell,j)
  do jc=1,self%hcurl_rep%nce
    self%cache_curl(jc,oft_tid)=self%curl_vals(j(jc))
  end do
  deallocate(j)
  allocate(j(self%hgrad_rep%nce))
  call self%hgrad_rep%ncdofs(cell,j)
  do jc=1,self%hgrad_rep%nce
    self%cache_grad(jc,oft_tid)=self%grad_vals(j(jc))
  end do
  deallocate(j)
  self%cache_cell(oft_tid)=cell
END IF
!---Reconstruct field
val=0.d0
IF(ASSOCIATED(oft_hcurl_cop))THEN
  do jc=1,self%hcurl_rep%nce
    val=val+self%cache_curl(jc,oft_tid)*oft_hcurl_cop(:,jc)
  end do
ELSE
  ALLOCATE(hcurl_cop(3,self%hcurl_rep%nce))
  CALL oft_hcurl_get_cgops(gop,cgop)
  CALL oft_hcurl_ceval_all(self%hcurl_rep,cell,f,hcurl_cop,cgop)
  do jc=1,self%hcurl_rep%nce
    val=val+self%cache_curl(jc,oft_tid)*hcurl_cop(:,jc)
  end do
  DEALLOCATE(hcurl_cop)
END IF
DEBUG_STACK_POP
end subroutine cinterp_apply
!------------------------------------------------------------------------------
!> Reconstruct \f$ \nabla \cdot \f$ of a H(Curl) + Grad(H^1) vector field
!------------------------------------------------------------------------------
subroutine dinterp_apply(self,cell,f,gop,val)
class(oft_hcurl_grad_dinterp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [4]
real(r8), intent(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
real(r8), intent(out) :: val(:) !< Reconstructed field at f [1]
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: dop(6),g2op(6,10),Kmat(10,3)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%grad_vals))CALL oft_abort('Setup has not been called!','dinterp_apply',__FILE__)
!---Pull cache
IF(self%cache_cell(oft_tid)/=cell)THEN
  allocate(j(self%hcurl_rep%nce))
  call self%hcurl_rep%ncdofs(cell,j)
  do jc=1,self%hcurl_rep%nce
    self%cache_curl(jc,oft_tid)=self%curl_vals(j(jc))
  end do
  deallocate(j)
  !
  allocate(j(self%hgrad_rep%nce))
  call self%hgrad_rep%ncdofs(cell,j)
  do jc=1,self%hgrad_rep%nce
    self%cache_grad(jc,oft_tid)=self%grad_vals(j(jc))
  end do
  deallocate(j)
  self%cache_cell(oft_tid)=cell
END IF
!---Reconstruct field
val=0.d0
CALL self%hcurl_rep%mesh%hessian(cell,f,g2op,Kmat)
do jc=1,self%hgrad_rep%nce
  call oft_h1_d2eval(self%hgrad_rep,cell,jc,f,dop,g2op)
  val=val+self%cache_grad(jc,oft_tid)*(dop(1)+dop(4)+dop(6))
end do
DEBUG_STACK_POP
end subroutine dinterp_apply
!------------------------------------------------------------------------------
!> Zero a H(Curl) + Grad(H^1) vector field at all boundary nodes
!------------------------------------------------------------------------------
subroutine zerob_apply(self,a)
class(oft_hcurl_grad_zerob), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to be zeroed
real(r8), pointer, dimension(:) :: agrad,acurl
integer(i4) :: i,j
class(oft_afem_type), pointer :: hcurl_rep,hgrad_rep
DEBUG_STACK_PUSH
!---Cast to vector type
NULLIFY(agrad,acurl)
CALL a%get_local(acurl,1)
CALL a%get_local(agrad,2)
hcurl_rep=>self%ML_hcurl_grad_rep%current_level%fields(1)%fe
hgrad_rep=>self%ML_hcurl_grad_rep%current_level%fields(2)%fe
! Apply operator
do i=1,hcurl_rep%nbe
  j=hcurl_rep%lbe(i)
  if(hcurl_rep%global%gbe(j))acurl(j)=0.d0
end do
! Apply operator
do i=1,hgrad_rep%nbe
  j=hgrad_rep%lbe(i)
  if(hgrad_rep%global%gbe(j))agrad(j)=0.d0
end do
CALL a%restore_local(acurl,1)
CALL a%restore_local(agrad,2)
DEALLOCATE(acurl,agrad)
DEBUG_STACK_POP
end subroutine zerob_apply
!------------------------------------------------------------------------------
!> Zero a H(Curl) + Grad(H^1) vector field at all boundary nodes
!------------------------------------------------------------------------------
subroutine zerob_delete(self)
class(oft_hcurl_grad_zerob), intent(inout) :: self
NULLIFY(self%ML_hcurl_grad_rep)
end subroutine zerob_delete
!------------------------------------------------------------------------------
!> Zero the curl components of a H(Curl) + Grad(H^1) vector field at all boundary nodes
!------------------------------------------------------------------------------
subroutine curl_zerob_apply(self,a)
class(oft_hcurl_grad_czerob), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to be zeroed
real(r8), pointer, dimension(:) :: acurl
integer(i4) :: i,j
class(oft_afem_type), pointer :: hcurl_rep
DEBUG_STACK_PUSH
SELECT TYPE(this=>self%ML_hcurl_grad_rep%current_level%fields(1)%fe)
CLASS IS(oft_fem_type)
  !---Cast to vector type
  NULLIFY(acurl)
  CALL a%get_local(acurl,1)
  ! Apply operator
  do i=1,this%nbe
    j=this%lbe(i)
    if(this%global%gbe(j))acurl(j)=0.d0
  end do
  CALL a%restore_local(acurl,1)
  DEALLOCATE(acurl)
CLASS DEFAULT
  CALL oft_abort("Invalid fe object","curl_zerob_apply",__FILE__)
END SELECT
DEBUG_STACK_POP
end subroutine curl_zerob_apply
!------------------------------------------------------------------------------
!> Zero the gradient components of a H(Curl) + Grad(H^1) vector field at all boundary nodes
!------------------------------------------------------------------------------
subroutine grad_zerop_apply(self,a)
class(oft_hcurl_grad_gzerop), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to be zeroed
real(r8), pointer, dimension(:) :: agrad
integer(i4) :: i,j
DEBUG_STACK_PUSH
SELECT TYPE(this=>self%ML_hcurl_grad_rep%current_level%fields(2)%fe)
CLASS IS(oft_fem_type)
  !---Cast to vector type
  NULLIFY(agrad)
  CALL a%get_local(agrad,2)
  ! Apply operator
  do i=1,this%mesh%np
    agrad(i)=0.d0
  end do
  CALL a%restore_local(agrad,2)
  DEALLOCATE(agrad)
CLASS DEFAULT
  CALL oft_abort("Invalid fe object","grad_zerop_apply",__FILE__)
END SELECT
DEBUG_STACK_POP
end subroutine grad_zerop_apply
!------------------------------------------------------------------------------
!> Zero a H(Curl) + Grad(H^1) vector field at all interior nodes
!------------------------------------------------------------------------------
subroutine zeroi_apply(self,a)
class(oft_hcurl_grad_zeroi), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to be zeroed
real(r8), pointer, dimension(:) :: agrad,acurl
integer(i4) :: i
class(oft_afem_type), pointer :: hcurl_rep,hgrad_rep
DEBUG_STACK_PUSH
!---Cast to vector type
NULLIFY(agrad,acurl)
CALL a%get_local(acurl,1)
CALL a%get_local(agrad,2)
hcurl_rep=>self%ML_hcurl_grad_rep%current_level%fields(1)%fe
hgrad_rep=>self%ML_hcurl_grad_rep%current_level%fields(2)%fe
! Apply operator
DO i=1,hcurl_rep%ne
  IF(hcurl_rep%global%gbe(i))CYCLE
  acurl(i)=0.d0
END DO
! Apply operator
DO i=1,hgrad_rep%ne
  IF(hgrad_rep%global%gbe(i))CYCLE
  agrad(i)=0.d0
END DO
CALL a%restore_local(acurl,1)
CALL a%restore_local(agrad,2)
DEALLOCATE(acurl,agrad)
DEBUG_STACK_POP
end subroutine zeroi_apply
!------------------------------------------------------------------------------
!> Compute the 0-th order gradient due to a jump plane
!!
!! The jump is represented as a circular surface defined by a center possition
!! and surface normal. This method requires that the mesh contain a
!! matching internal surface, such that no edge crosses the jump plane.
!!
!! @note The radius of the surface is represented by \f$ \left| hcpv \right| \f$
!------------------------------------------------------------------------------
subroutine hcurl_grad_mc(mesh,a,hcpc,hcpv,new_tol)
class(oft_mesh), intent(inout) :: mesh
class(oft_vector), intent(inout) :: a !< Jump field
real(r8), intent(in) :: hcpc(3) !< Jump plane center possition [3]
real(r8), intent(in) :: hcpv(3) !< Jump plane normal vector [3]
real(r8), optional, intent(in) :: new_tol
real(r8), pointer, dimension(:) :: acurl
integer(i4) :: i,j,k,l
real(r8) :: reg
real(r8) :: r1(3),r2(3),r1dv,r2dv,r1cv,tol
DEBUG_STACK_PUSH
tol=1.d-8
IF(PRESENT(new_tol))tol=new_tol
!---Get local values
NULLIFY(acurl)
CALL a%get_local(acurl,1)
!$omp parallel do private(r1,r2,r1dv,r2dv,r1cv)
do i=1,mesh%ne
  acurl(i)=0.d0
  r1=mesh%r(:,mesh%le(1,i))-hcpc(:)
  r2=mesh%r(:,mesh%le(2,i))-hcpc(:)
  r1dv=DOT_PRODUCT(r1,hcpv(:))
  r2dv=DOT_PRODUCT(r2,hcpv(:))
  if (XOR(ABS(r1dv)<tol,ABS(r2dv)<tol)) then
    if (ABS(r2dv)<tol) r1=r2
    r1cv=SUM(cross_product(r1,hcpv)**2)
    if(r1cv<1.d0)acurl(i)=1.d0*SIGN(.5d0,r2dv-r1dv)*SIGN(INT(1,8),mesh%global%le(i))
  end if
end do
CALL a%restore_local(acurl,1)
DEALLOCATE(acurl)
DEBUG_STACK_POP
end subroutine hcurl_grad_mc
!------------------------------------------------------------------------------
!> Compute the 0-th order gradient due to a jump plane
!!
!! The jump is represented as a circular surface defined by a center possition
!! and surface normal. This method requires that the mesh contain a
!! matching internal surface, such that no edge crosses the jump plane.
!!
!! @note The radius of the surface is represented by \f$ \left| hcpv \right| \f$
!------------------------------------------------------------------------------
subroutine hcurl_grad_bmc(mesh,a,hcpc,hcpv,new_tol)
class(oft_mesh), intent(inout) :: mesh
class(oft_vector), intent(inout) :: a !< Jump field
real(r8), intent(in) :: hcpc(3) !< Jump plane center possition [3]
real(r8), intent(in) :: hcpv(3) !< Jump plane normal vector [3]
real(r8), optional, intent(in) :: new_tol
real(r8), pointer, dimension(:) :: acurl
integer :: i,j,k,l,ed,mark
real(r8) :: reg
real(r8) :: r1(3),r2(3),r1dv,r2dv,r1cv,tol
DEBUG_STACK_PUSH
tol=1.d-8
IF(PRESENT(new_tol))tol=new_tol
!---Get local values
NULLIFY(acurl)
CALL a%get_local(acurl,1)
!$omp parallel do private(r1,r2,r1dv,r2dv,r1cv)
DO i=1,mesh%ne
  acurl(i)=0.d0
  IF(.NOT.mesh%global%gbe(i))CYCLE
  r1=mesh%r(:,mesh%le(1,i))-hcpc(:)
  r2=mesh%r(:,mesh%le(2,i))-hcpc(:)
  r1dv=DOT_PRODUCT(r1,hcpv(:))
  r2dv=DOT_PRODUCT(r2,hcpv(:))
  IF(ABS(r1dv)<=tol.AND.ABS(r2dv)<=tol)CALL oft_abort('Bad edge found','hcurl_grad_bmc',__FILE__)
  IF((r1dv>=tol.AND.r2dv<=-tol).OR.(r1dv<=-tol.AND.r2dv>=tol))THEN
    r1 = (r1+r2)/2.d0
    r1cv=SUM(cross_product(r1,hcpv)**2)
    IF(r1cv<1.d0)acurl(i)=SIGN(1.d0,r2dv-r1dv)*SIGN(INT(1,8),mesh%global%le(i))
  END IF
END DO
CALL a%restore_local(acurl,1)
DEALLOCATE(acurl)
DEBUG_STACK_POP
end subroutine hcurl_grad_bmc
!------------------------------------------------------------------------------
!> Add the gradient of a H^1 scalar field to a H(Curl) + Grad(H^1) vector field
!!
!! @note By default the 0-th order gradient subspace is represented on the
!! H(Curl) DOF, use the `keep_boundary` flag otherwise
!------------------------------------------------------------------------------
subroutine hcurl_grad_grad(hcurl_grad_fe,a,b,keep_boundary)
class(oft_fem_comp_type), intent(inout) :: hcurl_grad_fe
class(oft_vector), intent(inout) :: a !< Scalar field
class(oft_vector), intent(inout) :: b !< Vector field for gradient
logical, OPTIONAL, INTENT(in) :: keep_boundary !< Flag to keep 0-th order boundary component (optional)
real(r8), pointer, dimension(:) :: aloc
real(r8), pointer, dimension(:) :: bgrad,bcurl
integer(i4), allocatable :: emap(:)
integer(i4) :: i,j,k,l
real(r8) :: reg
LOGICAL :: zero_boundary
class(oft_mesh), pointer :: mesh
CLASS(oft_hcurl_fem), POINTER :: curl_rep
CLASS(oft_h1_fem), POINTER :: grad_rep
DEBUG_STACK_PUSH
NULLIFY(aloc,bcurl,bgrad)
zero_boundary=.FALSE.
IF(PRESENT(keep_boundary))zero_boundary=keep_boundary
IF(.NOT.oft_3D_hcurl_cast(curl_rep,hcurl_grad_fe%fields(1)%fe))CALL oft_abort("Incorrect HCurl FE type","hcurl_grad_grad",__FILE__)
IF(.NOT.oft_3D_h1_cast(grad_rep,hcurl_grad_fe%fields(2)%fe))CALL oft_abort("Incorrect HCurl FE type","hcurl_grad_grad",__FILE__)
mesh=>curl_rep%mesh
!---Get local values
CALL a%get_local(aloc)
!---Zero output and get aliases
CALL b%get_local(bcurl,1)
CALL b%get_local(bgrad,2)
!---Diagonal part
IF(zero_boundary)THEN
  DO i=1,mesh%np
    IF(mesh%global%gbp(i))THEN
      bgrad(i)=bgrad(i)+aloc(i)
      aloc(i)=0.d0
    END IF
  END DO
END IF
DO i=mesh%np+1,grad_rep%ne
  bgrad(i)=bgrad(i)+aloc(i)
END DO
!---Off-diagonal part
DO i=1,mesh%ne
  ! bcurl(i)=bcurl(i) + &
  ! (aloc(mesh%le(2,i))-aloc(mesh%le(1,i)))*SIGN(1_i8,mesh%global%le(i))
  bcurl((i-1)*curl_rep%gstruct(2)+1)=bcurl((i-1)*curl_rep%gstruct(2)+1) + &
    (aloc(mesh%le(2,i))-aloc(mesh%le(1,i)))*SIGN(1_i8,mesh%global%le(i))
END DO
!---
CALL b%restore_local(bcurl,1,wait=.TRUE.)
CALL b%restore_local(bgrad,2)
DEALLOCATE(aloc,bcurl,bgrad)
DEBUG_STACK_POP
end subroutine hcurl_grad_grad
!------------------------------------------------------------------------------
!> Apply the transposed gradient operator to a H(Curl) + Grad(H^1) vector field
!------------------------------------------------------------------------------
subroutine hcurl_grad_gradtp(h1_fe,a,b)
class(oft_afem_type), intent(inout) :: h1_fe
class(oft_vector), intent(inout) :: a !< Input field
class(oft_vector), intent(inout) :: b !< \f$ G^{T} a \f$
real(r8), pointer, dimension(:) :: agrad,acurl
real(r8), pointer, dimension(:) :: bloc
integer(i4), allocatable, dimension(:) :: emap
integer(i4) :: i,j,k,l
real(r8) :: reg
CLASS(oft_h1_fem), POINTER :: grad_rep
CLASS(oft_mesh), POINTER :: mesh
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_h1_cast(grad_rep,h1_fe))CALL oft_abort("Incorrect Grad FE type","hcurl_grad_gradtp",__FILE__)
mesh=>grad_rep%mesh
NULLIFY(acurl,agrad,bloc)
!---Cast local values
CALL a%get_local(acurl,1)
CALL a%get_local(agrad,2)
!---Zero output and get aliases
call b%set(0.d0)
CALL b%get_local(bloc)
!---Get boundary map
allocate(emap(mesh%ne))
CALL get_inverse_map(mesh%lbe,mesh%nbe,emap,mesh%ne)
!---Compute off-diagonal GT
!$omp parallel do private(j,k,l,reg)
do i=1,mesh%np
  bloc(i)=0.d0
  do j=mesh%kpe(i),mesh%kpe(i+1)-1
    k=mesh%lpe(j)
    if(.NOT.mesh%fullmesh)then
      if(emap(k)>0)then
        ! if(.NOT.mesh%linkage%leo(emap(k)))cycle
        if(.NOT.mesh%estitch%leo(emap(k)))cycle
      end if
    end if
    if(mesh%global%le(k)<0)then
      bloc(i)=bloc(i)+acurl(k)*(1-2*min(1,abs(mesh%le(1,k)-i)))
    else
      bloc(i)=bloc(i)-acurl(k)*(1-2*min(1,abs(mesh%le(1,k)-i)))
    end if
  enddo
enddo
deallocate(emap)
!---Add diagonal part of GT
!$omp parallel do
do i=1,grad_rep%ne
  bloc(i)=bloc(i)+agrad(i)
end do
CALL b%restore_local(bloc,add=.TRUE.)
DEALLOCATE(bloc,acurl,agrad)
DEBUG_STACK_POP
end subroutine hcurl_grad_gradtp
!------------------------------------------------------------------------------
!> Apply the divergence operator to a H(Curl) + Grad(H^1) field
!------------------------------------------------------------------------------
subroutine hcurl_grad_div(hcurl_grad_fe,a,b)
class(oft_fem_comp_type), intent(inout) :: hcurl_grad_fe
class(oft_vector), intent(inout) :: a !< Needs docs
class(oft_vector), intent(inout) :: b !< Needs docs
integer(i4) :: i,m,jr
integer(i4), allocatable :: j_curl(:),j_grad(:)
real(r8) :: vol,det,goptmp(3,4),aloc(3)
real(r8), pointer, contiguous, dimension(:) :: ac_loc,ag_loc,bloc
real(r8), allocatable, dimension(:) :: ac_tmp,ag_tmp,btmp
real(r8), allocatable, dimension(:,:) :: rop_curl,rop_grad
logical :: curved
CLASS(oft_hcurl_fem), POINTER :: curl_rep
CLASS(oft_h1_fem), POINTER :: grad_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_hcurl_cast(curl_rep,hcurl_grad_fe%fields(1)%fe))CALL oft_abort("Incorrect Curl FE type","hcurl_grad_div",__FILE__)
IF(.NOT.oft_3D_h1_cast(grad_rep,hcurl_grad_fe%fields(2)%fe))CALL oft_abort("Incorrect Grad FE type","hcurl_grad_div",__FILE__)
!------------------------------------------------------------------------------
! Allocate Laplacian Op
!------------------------------------------------------------------------------
NULLIFY(ac_loc,ag_loc,bloc)
!---Get local values
CALL a%get_local(ac_loc,1)
CALL a%get_local(ag_loc,2)
!---Zero output and get aliases
CALL b%set(0.d0)
CALL b%get_local(bloc)
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j_curl,j_grad,rop_curl,rop_grad,det,btmp,ac_tmp,ag_tmp, &
!$omp aloc,curved,goptmp,vol,m,jr)
allocate(j_curl(curl_rep%nce),j_grad(grad_rep%nce))
allocate(rop_grad(3,grad_rep%nce),rop_curl(3,curl_rep%nce))
allocate(ac_tmp(curl_rep%nce),ag_tmp(grad_rep%nce))
allocate(btmp(grad_rep%nce))
!$omp do schedule(guided)
do i=1,grad_rep%mesh%nc
  !---Get local to global DOF mapping
  call curl_rep%ncdofs(i,j_curl)
  call grad_rep%ncdofs(i,j_grad)
  curved=cell_is_curved(grad_rep%mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  do jr=1,curl_rep%nce
    ac_tmp(jr)=ac_loc(j_curl(jr))
  end do
  do jr=1,grad_rep%nce
    ag_tmp(jr)=ag_loc(j_grad(jr))
  end do
  btmp=0.d0
  do m=1,grad_rep%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call grad_rep%mesh%jacobian(i,grad_rep%quad%pts(:,m),goptmp,vol)
    det=vol*grad_rep%quad%wts(m)
    call oft_hcurl_eval_all(curl_rep,i,grad_rep%quad%pts(:,m),rop_curl,goptmp)
    call oft_h1_geval_all(grad_rep,i,grad_rep%quad%pts(:,m),rop_grad,goptmp)
    !---Compute local operator contribution
    aloc = 0.d0
    do jr=1,curl_rep%nce
      aloc = aloc + rop_curl(:,jr)*ac_tmp(jr)
    end do
    do jr=1,grad_rep%nce
      aloc = aloc + rop_grad(:,jr)*ag_tmp(jr)
    end do
    do jr=1,grad_rep%nce
      btmp(jr)=btmp(jr)+DOT_PRODUCT(rop_grad(:,jr),aloc)*det
    end do
  end do
  !---Add local values to global vector
  do jr=1,grad_rep%nce
    !$omp atomic
    bloc(j_grad(jr))=bloc(j_grad(jr))+btmp(jr)
  end do
end do
deallocate(j_curl,j_grad)
deallocate(rop_curl,rop_grad,btmp,ac_tmp,ag_tmp)
!$omp end parallel
CALL b%restore_local(bloc,add=.TRUE.)
DEALLOCATE(ac_loc,ag_loc,bloc)
DEBUG_STACK_POP
end subroutine hcurl_grad_div
!------------------------------------------------------------------------------
!> Apply the curl transpose operator to a H(Curl) + Grad(H^1) field
!------------------------------------------------------------------------------
subroutine hcurl_grad_curltp(hcurl_grad_fe,a,b)
class(oft_fem_comp_type), intent(inout) :: hcurl_grad_fe
class(oft_vector), intent(inout) :: a !< Needs docs
class(oft_vector), intent(inout) :: b !< Needs docs
real(r8), pointer, dimension(:) :: agrad,acurl
real(r8), pointer, dimension(:) :: bcurl
integer(i4) :: i,jr,jc,m
integer(i4), allocatable :: j_curl(:),j_grad(:)
real(r8) :: goptmp(3,4),cgop(3,6)
real(r8) :: v,f(4),det,vol
logical :: curved
real(r8), allocatable :: rop_curl(:,:),cop_curl(:,:),rop_grad(:,:),btmp(:)
CLASS(oft_hcurl_fem), POINTER :: curl_rep
CLASS(oft_h1_fem), POINTER :: grad_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_hcurl_cast(curl_rep,hcurl_grad_fe%fields(1)%fe))CALL oft_abort("Incorrect Curl FE type","hcurl_grad_curltp",__FILE__)
IF(.NOT.oft_3D_h1_cast(grad_rep,hcurl_grad_fe%fields(2)%fe))CALL oft_abort("Incorrect Grad FE type","hcurl_grad_curltp",__FILE__)
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
NULLIFY(acurl,agrad,bcurl)
!---Get local values
CALL a%get_local(acurl,1)
CALL a%get_local(agrad,2)
!---Zero output and get aliases
CALL b%set(0.d0)
CALL b%get_local(bcurl,1)
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j_curl,j_grad,rop_curl,cop_curl,rop_grad,det,btmp,curved,goptmp,cgop,m,v,jc,jr,f)
allocate(j_curl(curl_rep%nce)) ! Local DOF and matrix indices
allocate(j_grad(grad_rep%nce)) ! Local DOF and matrix indices
allocate(rop_grad(3,grad_rep%nce)) ! Reconstructed gradient operator
allocate(rop_curl(3,curl_rep%nce)) ! Reconstructed gradient operator
allocate(cop_curl(3,curl_rep%nce)) ! Reconstructed gradient operator
allocate(btmp(curl_rep%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,grad_rep%mesh%nc
  !---Get local to global DOF mapping
  call curl_rep%ncdofs(i,j_curl)
  call grad_rep%ncdofs(i,j_grad)
  curved=cell_is_curved(grad_rep%mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  btmp=0.d0
  do m=1,curl_rep%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)then
      call grad_rep%mesh%jacobian(i,curl_rep%quad%pts(:,m),goptmp,v)
      CALL oft_hcurl_get_cgops(goptmp,cgop)
    end if
    det=v*curl_rep%quad%wts(m)
    call oft_hcurl_eval_all(curl_rep,i,curl_rep%quad%pts(:,m),rop_curl,goptmp)
    call oft_hcurl_ceval_all(curl_rep,i,curl_rep%quad%pts(:,m),cop_curl,cgop)
    call oft_h1_geval_all(grad_rep,i,curl_rep%quad%pts(:,m),rop_grad,goptmp)
    !---Compute local operator contribution
    do jr=1,curl_rep%nce
      do jc=1,curl_rep%nce
        btmp(jr)=btmp(jr)+DOT_PRODUCT(cop_curl(:,jr),rop_curl(:,jc))*acurl(j_curl(jc))*det
      end do
      do jc=1,grad_rep%nce
        btmp(jr)=btmp(jr)+DOT_PRODUCT(cop_curl(:,jr),rop_grad(:,jc))*agrad(j_grad(jc))*det
      end do
    end do
  end do
  !---Add local values to global vector
  do jr=1,curl_rep%nce
    !$omp atomic
    bcurl(j_curl(jr))=bcurl(j_curl(jr))+btmp(jr)
  end do
end do
deallocate(j_curl,j_grad)
deallocate(rop_curl,cop_curl,rop_grad,btmp)
!$omp end parallel
CALL b%restore_local(bcurl,1,add=.TRUE.)
DEALLOCATE(acurl,agrad,bcurl)
DEBUG_STACK_POP
end subroutine hcurl_grad_curltp
!------------------------------------------------------------------------------
!> Construct mass matrix for a H(Curl) + Grad(H^1) representation
!!
!! Supported boundary conditions
!! - `'none'` Full matrix
!! - `'zerob'` Dirichlet for all boundary DOF
!------------------------------------------------------------------------------
subroutine hcurl_grad_getmop(hcurl_grad_rep,mat,bc)
class(oft_fem_comp_type), intent(inout) :: hcurl_grad_rep
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix object
character(LEN=*), intent(in) :: bc !< Boundary condition
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j_curl(:),j_grad(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time
real(r8), allocatable, dimension(:,:) :: rop_curl,rop_grad
real(r8), allocatable, dimension(:,:) :: mop11,mop12,mop21,mop22
logical :: curved
CLASS(oft_vector), POINTER :: oft_hcurl_grad_vec
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
CLASS(oft_h1_fem), POINTER :: hgrad_rep
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing H(Curl) + Grad(H^1) mass matrix'
  CALL mytimer%tick()
END IF
!---
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,hcurl_grad_rep%fields(1)%fe))CALL oft_abort("Incorrect Curl FE type","hcurl_grad_getmop",__FILE__)
IF(.NOT.oft_3D_h1_cast(hgrad_rep,hcurl_grad_rep%fields(2)%fe))CALL oft_abort("Incorrect Grad FE type","hcurl_grad_getmop",__FILE__)
!------------------------------------------------------------------------------
! Allocate matrix
!------------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL hcurl_grad_rep%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!------------------------------------------------------------------------------
! Operator integration
!------------------------------------------------------------------------------
!$omp parallel private(j_curl,j_grad,rop_curl,rop_grad,det,mop11,mop12,mop21,mop22, &
!$omp curved,goptmp,m,vol,jc,jr)
allocate(j_curl(hcurl_rep%nce)) ! Local DOF and matrix indices
allocate(j_grad(hgrad_rep%nce)) ! Local DOF and matrix indices
allocate(rop_grad(3,hgrad_rep%nce)) ! Reconstructed gradient operator
allocate(rop_curl(3,hcurl_rep%nce)) ! Reconstructed gradient operator
allocate(mop11(hcurl_rep%nce,hcurl_rep%nce)) ! Local laplacian matrix
allocate(mop12(hcurl_rep%nce,hgrad_rep%nce)) ! Local laplacian matrix
allocate(mop21(hgrad_rep%nce,hcurl_rep%nce)) ! Local laplacian matrix
allocate(mop22(hgrad_rep%nce,hgrad_rep%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,hgrad_rep%mesh%nc
  curved=cell_is_curved(hgrad_rep%mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  mop11=0.d0; mop12=0.d0
  mop21=0.d0; mop22=0.d0
  do m=1,hcurl_rep%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call hgrad_rep%mesh%jacobian(i,hcurl_rep%quad%pts(:,m),goptmp,vol)
    det=vol*hcurl_rep%quad%wts(m)
    call oft_hcurl_eval_all(hcurl_rep,i,hcurl_rep%quad%pts(:,m),rop_curl,goptmp)
    call oft_h1_geval_all(hgrad_rep,i,hcurl_rep%quad%pts(:,m),rop_grad,goptmp)
    !---Compute local matrix contributions
    do jr=1,hcurl_rep%nce
      do jc=1,hcurl_rep%nce
        mop11(jr,jc) = mop11(jr,jc) + DOT_PRODUCT(rop_curl(:,jr),rop_curl(:,jc))*det
      end do
      do jc=1,hgrad_rep%nce
        mop12(jr,jc) = mop12(jr,jc) + DOT_PRODUCT(rop_curl(:,jr),rop_grad(:,jc))*det
      end do
    end do
    do jr=1,hgrad_rep%nce
      do jc=1,hgrad_rep%nce
        mop22(jr,jc) = mop22(jr,jc) + DOT_PRODUCT(rop_grad(:,jr),rop_grad(:,jc))*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call hcurl_rep%ncdofs(i,j_curl)
  call hgrad_rep%ncdofs(i,j_grad)
  mop21=TRANSPOSE(mop12)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("none")
      mop21(1:hgrad_rep%mesh%cell_np,:)=0.d0
      mop22(1:hgrad_rep%mesh%cell_np,:)=0.d0
    CASE("zerob")
      DO jr=1,hcurl_rep%nce
        IF(hcurl_rep%global%gbe(j_curl(jr)))THEN
          mop11(jr,:)=0.d0
          mop12(jr,:)=0.d0
        END IF
      END DO
      mop21(1:hgrad_rep%mesh%cell_np,:)=0.d0
      mop22(1:hgrad_rep%mesh%cell_np,:)=0.d0
      DO jr=hgrad_rep%mesh%cell_np,hgrad_rep%nce
        IF(hgrad_rep%global%gbe(j_grad(jr)))THEN
          mop12(jr,:)=0.d0
          mop22(jr,:)=0.d0
        END IF
      END DO
  END SELECT
  !---Add local values to global matrix
  call mat%atomic_add_values(j_curl,j_curl,mop11,hcurl_rep%nce,hcurl_rep%nce,1,1)
  call mat%atomic_add_values(j_curl,j_grad,mop12,hcurl_rep%nce,hgrad_rep%nce,1,2)
  call mat%atomic_add_values(j_grad,j_curl,mop21,hgrad_rep%nce,hcurl_rep%nce,2,1)
  call mat%atomic_add_values(j_grad,j_grad,mop22,hgrad_rep%nce,hgrad_rep%nce,2,2)
end do
deallocate(j_curl,j_grad)
deallocate(rop_curl,rop_grad)
deallocate(mop11,mop12,mop21,mop22)
!$omp end parallel
ALLOCATE(mop11(1,1),j_curl(1))
!---Set diagonal entries for dirichlet rows
SELECT CASE(TRIM(bc))
  CASE("none")
    mop11(1,1)=1.d0
    DO i=1,hgrad_rep%mesh%np
      IF(hgrad_rep%mesh%bp(i))CYCLE
      j_curl=i
      call mat%add_values(j_curl,j_curl,mop11,1,1,2,2)
    END DO
    DO i=1,hgrad_rep%mesh%nbp
      jr=hgrad_rep%mesh%lbp(i)
      ! IF(.NOT.hgrad_rep%mesh%linkage%lpo(i))CYCLE
      IF(.NOT.hgrad_rep%mesh%pstitch%leo(i))CYCLE
      j_curl=jr
      call mat%add_values(j_curl,j_curl,mop11,1,1,2,2)
    END DO
  CASE("zerob")
    mop11(1,1)=1.d0
    DO i=1,hcurl_rep%nbe
      jr=hcurl_rep%lbe(i)
      IF(.NOT.hcurl_rep%global%gbe(jr))CYCLE
      IF(.NOT.hcurl_rep%linkage%leo(i))CYCLE
      j_curl=jr
      call mat%add_values(j_curl,j_curl,mop11,1,1,1,1)
    END DO
    DO i=1,hgrad_rep%mesh%np
      IF(.NOT.hgrad_rep%linkage%be(i))CYCLE
      j_curl=i
      call mat%add_values(j_curl,j_curl,mop11,1,1,1,1)
    END DO
    DO i=1,hgrad_rep%nbe
      jr=hgrad_rep%lbe(i)
      IF(.NOT.hgrad_rep%global%gbe(jr))CYCLE
      IF(.NOT.hgrad_rep%linkage%leo(i))CYCLE
      j_curl=jr
      call mat%add_values(j_curl,j_curl,mop11,1,1,1,1)
    END DO
END SELECT
DEALLOCATE(j_curl,mop11)
CALL hcurl_grad_rep%vec_create(oft_hcurl_grad_vec)
CALL mat%assemble(oft_hcurl_grad_vec)
CALL oft_hcurl_grad_vec%delete
DEALLOCATE(oft_hcurl_grad_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine hcurl_grad_getmop
!------------------------------------------------------------------------------
!> Project a vector field onto a H(Curl) + Grad(H^1) basis
!!
!! @note This subroutine only performs the integration of the field with
!! test functions for a H(Curl) + Grad(H^1) basis. To retrieve the correct projection the
!! result must be multiplied by the inverse of the mass matrix
!------------------------------------------------------------------------------
subroutine oft_hcurl_grad_project(hcurl_grad_rep,field,x)
class(oft_fem_comp_type), intent(inout) :: hcurl_grad_rep
class(fem_interp), intent(inout) :: field !< Vector field for projection
class(oft_vector), intent(inout) :: x !< Field projected onto H(Curl) + Grad(H^1) basis
integer(i4) :: i,jc,m
integer(i4), allocatable, dimension(:) :: j_hcurl,j_hgrad
real(r8) :: det,vol,bcc(3),goptmp(3,4)
real(r8), pointer, dimension(:) :: xcurl,xgrad
real(r8), allocatable, dimension(:,:) :: rop_curl,rop_grad
logical :: curved
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
CLASS(oft_h1_fem), POINTER :: hgrad_rep
DEBUG_STACK_PUSH
!---
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,hcurl_grad_rep%fields(1)%fe))CALL oft_abort("Incorrect Curl FE type","hcurl_grad_getmop",__FILE__)
IF(.NOT.oft_3D_h1_cast(hgrad_rep,hcurl_grad_rep%fields(2)%fe))CALL oft_abort("Incorrect Grad FE type","hcurl_grad_getmop",__FILE__)
!---Initialize vectors to zero
NULLIFY(xcurl,xgrad)
call x%set(0.d0)
call x%get_local(xcurl,1)
call x%get_local(xgrad,2)
!---Integerate over the volume
!$omp parallel default(firstprivate) shared(xcurl,xgrad) private(curved,det)
allocate(j_hcurl(hcurl_rep%nce),rop_curl(3,hcurl_rep%nce))
allocate(j_hgrad(hgrad_rep%nce),rop_grad(3,hgrad_rep%nce))
!$omp do schedule(guided)
do i=1,hgrad_rep%mesh%nc ! Loop over cells
  call hcurl_rep%ncdofs(i,j_hcurl) ! Get DOFs
  call hgrad_rep%ncdofs(i,j_hgrad) ! Get DOFs
  curved=cell_is_curved(hgrad_rep%mesh,i) ! Straight cell test
  do m=1,hcurl_rep%quad%np
    if(curved.OR.m==1)call hgrad_rep%mesh%jacobian(i,hcurl_rep%quad%pts(:,m),goptmp,vol)
    det=vol*hcurl_rep%quad%wts(m)
    call field%interp(i,hcurl_rep%quad%pts(:,m),goptmp,bcc)
    call oft_hcurl_eval_all(hcurl_rep,i,hcurl_rep%quad%pts(:,m),rop_curl,goptmp)
    call oft_h1_geval_all(hgrad_rep,i,hcurl_rep%quad%pts(:,m),rop_grad,goptmp)
    do jc=1,hcurl_rep%nce
      !$omp atomic
      xcurl(j_hcurl(jc))=xcurl(j_hcurl(jc))+DOT_PRODUCT(rop_curl(:,jc),bcc)*det
    end do
    do jc=1,hgrad_rep%nce
      !$omp atomic
      xgrad(j_hgrad(jc))=xgrad(j_hgrad(jc))+DOT_PRODUCT(rop_grad(:,jc),bcc)*det
    end do
  end do
end do
deallocate(j_hcurl,j_hgrad,rop_curl,rop_grad)
!$omp end parallel
call x%restore_local(xcurl,1,add=.TRUE.,wait=.TRUE.)
call x%restore_local(xgrad,2,add=.TRUE.)
deallocate(xcurl,xgrad)
DEBUG_STACK_POP
end subroutine oft_hcurl_grad_project
!---------------------------------------------------------------------------------
!> Boundary projection of a vector field onto a H(Curl) + Grad(H^1) basis
!!
!! @note This subroutine only performs the integration of the field with
!! boundary test functions for a H(Curl) + Grad(H^1) basis
!---------------------------------------------------------------------------------
SUBROUTINE oft_hcurl_grad_bproject(hcurl_grad_rep,field,x)
class(oft_fem_comp_type), intent(inout) :: hcurl_grad_rep
CLASS(fem_interp), INTENT(inout) :: field !< Vector field for projection
CLASS(oft_vector), INTENT(inout) :: x !< Field projected onto H(Curl) + Grad(H^1) basis
INTEGER(i4) :: i,m,jc,cf,face,cell,ptmap(3)
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: j_hcurl,j_hgrad
REAL(r8) :: vol,det,f(4),norm(3),etmp(3),goptmp(3,3),gop(3,4)
REAL(r8), POINTER, DIMENSION(:) :: xcurl,xgrad
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: rop_curl,rop_grad
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
CLASS(oft_h1_fem), POINTER :: hgrad_rep
CLASS(oft_mesh), POINTER :: mesh
CLASS(oft_bmesh), POINTER :: smesh
TYPE(oft_quad_type) :: quad
DEBUG_STACK_PUSH
!---
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,hcurl_grad_rep%fields(1)%fe))CALL oft_abort("Incorrect Curl FE type","oft_hcurl_grad_bproject",__FILE__)
IF(.NOT.oft_3D_h1_cast(hgrad_rep,hcurl_grad_rep%fields(2)%fe))CALL oft_abort("Incorrect Grad FE type","oft_hcurl_grad_bproject",__FILE__)
mesh=>hcurl_rep%mesh
smesh=>hcurl_rep%mesh%bmesh
CALL smesh%quad_rule(hcurl_rep%order*2+1,quad)
!---Initialize vectors to zero
NULLIFY(xcurl,xgrad)
call x%set(0.d0)
call x%get_local(xcurl,1)
call x%get_local(xgrad,2)
!---Operator integration loop
!$omp parallel default(firstprivate) shared(xcurl,xgrad) private(det)
allocate(j_hcurl(hcurl_rep%nce),rop_curl(3,hcurl_rep%nce))
allocate(j_hgrad(hgrad_rep%nce),rop_grad(3,hgrad_rep%nce))
!$omp do schedule(guided)
do i=1,smesh%nc
  CALL mesh%get_surf_map(i,cell,ptmap) ! Find parent cell and logical coordinate mapping
  !---Get local to global DOF mapping
  call hcurl_rep%ncdofs(cell,j_hcurl)
  call hgrad_rep%ncdofs(cell,j_hgrad)
  !---Get local reconstructed operators
  do m=1,quad%np ! Loop over quadrature points
    call smesh%jacobian(i,quad%pts(:,m),goptmp,vol)
    call smesh%norm(i,quad%pts(:,m),norm)
    det=vol*quad%wts(m)
    !---
    CALL mesh%surf_to_vol(quad%pts(:,m),ptmap,f)
    call mesh%jacobian(cell,f,gop,vol)
    call field%interp(cell,f,gop,etmp)
    call oft_hcurl_eval_all(hcurl_rep,cell,f,rop_curl,gop)
    call oft_h1_geval_all(hgrad_rep,cell,f,rop_grad,gop)
    do jc=1,hcurl_rep%nce
      !$omp atomic
      xcurl(j_hcurl(jc))=xcurl(j_hcurl(jc)) &
      + DOT_PRODUCT(cross_product(rop_curl(:,jc),etmp),norm)*det
    end do
    do jc=1,hgrad_rep%nce
      !$omp atomic
      xgrad(j_hgrad(jc))=xgrad(j_hgrad(jc)) &
      + DOT_PRODUCT(cross_product(rop_grad(:,jc),etmp),norm)*det
    end do
  end do
end do
deallocate(j_hcurl,j_hgrad,rop_curl,rop_grad)
!$omp end parallel
call x%restore_local(xcurl,1,add=.TRUE.,wait=.TRUE.)
call x%restore_local(xgrad,2,add=.TRUE.)
deallocate(xcurl,xgrad)
DEBUG_STACK_POP
END SUBROUTINE oft_hcurl_grad_bproject
!------------------------------------------------------------------------------
!> Setup matrix and solver with default
!------------------------------------------------------------------------------
subroutine divout_setup(self,ML_hcurl_grad_rep,bc,solver)
class(oft_hcurl_grad_divout), intent(inout) :: self
CLASS(oft_ml_fem_comp_type), target, intent(inout) :: ML_hcurl_grad_rep
character(LEN=*), intent(in) :: bc !< Boundary condition
class(oft_solver), target, optional, intent(in) :: solver
!---
CLASS(oft_matrix), POINTER :: lop
CLASS(oft_solver), POINTER :: linv
TYPE(oft_h1_zerob), POINTER :: bc_zerob
TYPE(oft_h1_zerogrnd), POINTER :: bc_zerogrnd
DEBUG_STACK_PUSH
self%ML_hcurl_full=>ML_hcurl_grad_rep
self%ML_curl=>ML_hcurl_grad_rep%ml_fields(1)%ml
self%ML_grad=>ML_hcurl_grad_rep%ml_fields(2)%ml
IF(PRESENT(solver))THEN
  self%solver=>solver
  self%internal_solver=.FALSE.
ELSE
  NULLIFY(lop)
  CALL oft_h1_getlop(self%ML_grad%current_level,lop,bc)
  CALL create_cg_solver(linv)
  linv%A=>lop
  linv%its=-3
  CALL create_diag_pre(linv%pre)
  self%solver=>linv
  self%internal_solver=.TRUE.
END IF
!
SELECT CASE(TRIM(bc)) 
  CASE('grnd')
    ALLOCATE(bc_zerogrnd)
    bc_zerogrnd%ML_H1_rep=>self%ML_grad
    self%bc=>bc_zerogrnd
  CASE('zero')
    ALLOCATE(bc_zerob)
    bc_zerob%ML_H1_rep=>self%ML_grad
    self%bc=>bc_zerob
  CASE('none')
    NULLIFY(self%bc)
  CASE DEFAULT
    CALL oft_abort("Invalid BC","divout_setup",__FILE__)
END SELECT
DEBUG_STACK_POP
end subroutine divout_setup
!------------------------------------------------------------------------------
!> Remove divergence from a H(Curl) + Grad(H^1) vector field by adding a gradient correction
!------------------------------------------------------------------------------
subroutine divout_apply(self,a)
class(oft_hcurl_grad_divout), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field for divergence cleaning
class(oft_vector), pointer :: u,g,tmp,tmp2
integer(i4) :: i,order_tmp
real(r8) :: uu
logical :: pm_save
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%solver))CALL oft_abort('No solver specified.','divout_apply',__FILE__)
IF(mod(self%count,self%app_freq)/=0)THEN
  DEBUG_STACK_POP
  RETURN
END IF
!---
call self%ML_grad%vec_create(g)
call self%ML_grad%vec_create(u)
!---
IF(ASSOCIATED(self%mop))THEN
  CALL hcurl_grad_gradtp(self%ML_grad%current_level,a,g)
ELSE
  CALL hcurl_grad_div(self%ML_hcurl_full%current_level,a,g)
END IF
uu=a%dot(a)
self%solver%atol=MAX(self%solver%atol,SQRT(uu*1.d-20))
call u%set(0.d0)
IF(ASSOCIATED(self%bnorm))CALL g%add(1.d0,-1.d0,self%bnorm)
call self%bc%apply(g)
!---
pm_save=oft_env%pm; oft_env%pm=self%pm
call self%solver%apply(u,g)
oft_env%pm=pm_save
!---
CALL u%scale(-1.d0)
CALL a%new(tmp)
CALL hcurl_grad_grad(self%ML_hcurl_full%current_level,u,tmp,self%keep_boundary)
IF(ASSOCIATED(self%mop))THEN
  CALL a%new(tmp2)
  CALL self%mop%apply(tmp,tmp2)
  CALL a%add(1.d0,1.d0,tmp2)
  CALL tmp2%delete()
  DEALLOCATE(tmp2)
ELSE
  CALL a%add(1.d0,1.d0,tmp)
END IF
CALL tmp%delete
call u%delete
call g%delete
DEALLOCATE(tmp,u,g)
DEBUG_STACK_POP
end subroutine divout_apply
!------------------------------------------------------------------------------
!> Clean-up internal storage for a oft_hcurl_grad_divout object
!------------------------------------------------------------------------------
subroutine divout_delete(self)
class(oft_hcurl_grad_divout), intent(inout) :: self
IF(ASSOCIATED(self%solver).AND.self%internal_solver)THEN
  call self%solver%A%delete
  deallocate(self%solver%A)
  call self%solver%pre%delete
  deallocate(self%solver%pre)
  call self%solver%delete
  deallocate(self%solver)
END IF
IF(ASSOCIATED(self%bc))THEN
  CALL self%bc%delete()
  DEALLOCATE(self%bc)
END IF
NULLIFY(self%ML_hcurl_full,self%ML_curl,self%ML_grad)
end subroutine divout_delete
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine zerograd_apply(self,a)
class(oft_hcurl_grad_gzero), intent(inout) :: self
class(oft_vector), intent(inout) :: a
real(r8), pointer, dimension(:) :: ugrad
DEBUG_STACK_PUSH
!---Get local field
NULLIFY(ugrad)
CALL a%get_local(ugrad,2)
ugrad=0.d0
CALL a%restore_local(ugrad,2)
DEALLOCATE(ugrad)
DEBUG_STACK_POP
end subroutine zerograd_apply
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine zerograd_delete(self)
class(oft_hcurl_grad_gzero), intent(inout) :: self
! Do nothing
end subroutine zerograd_delete
!------------------------------------------------------------------------------
!> Construct interpolation matrices on each MG level
!------------------------------------------------------------------------------
SUBROUTINE hcurl_grad_setup_interp(ML_hcurl_aug_rep,ML_h1_rep,create_full)
CLASS(oft_ml_fem_comp_type), intent(inout) :: ML_hcurl_aug_rep
CLASS(oft_ml_fem_type), intent(inout) :: ML_h1_rep
LOGICAL, OPTIONAL, INTENT(in) :: create_full
INTEGER(i4) :: i
LOGICAL :: full_interp
CLASS(oft_ml_fem_type), POINTER :: ML_grad,ML_curl
DEBUG_STACK_PUSH
full_interp=.FALSE.
IF(PRESENT(create_full))full_interp=create_full
!---
ML_curl=>ML_hcurl_aug_rep%ml_fields(1)%ml
ML_grad=>ML_hcurl_aug_rep%ml_fields(2)%ml
DO i=ML_hcurl_aug_rep%minlev+1,ML_hcurl_aug_rep%nlevels
  CALL ML_hcurl_aug_rep%set_level(i,propogate=.TRUE.)
  !---
  IF(ML_hcurl_aug_rep%level==ML_hcurl_aug_rep%blevel+1)CYCLE
  !---Setup interpolation operators for Grad(H^1) space
  IF(ML_curl%current_level%order==1)THEN
    CALL hgrad_ginterpmatrix(ML_grad%interp_matrices(ML_grad%level)%m)
    CALL ML_grad%interp_matrices(ML_grad%level)%m%assemble
  ELSE
    CALL ML_h1_rep%set_level(ML_hcurl_aug_rep%level+1)
    ML_grad%interp_graphs(ML_grad%level)%g=>ML_h1_rep%interp_graphs(ML_hcurl_aug_rep%level+1)%g
    ML_grad%interp_matrices(ML_grad%level)%m=>ML_h1_rep%interp_matrices(ML_hcurl_aug_rep%level+1)%m
  END IF
END DO
!---Create full H(Curl) + Grad(H^1) interpolation operator
IF(full_interp)CALL ML_hcurl_aug_rep%build_interp
DEBUG_STACK_POP
CONTAINS
!------------------------------------------------------------------------------
!> Construct interpolation matrix for polynomial levels
!------------------------------------------------------------------------------
SUBROUTINE hgrad_ginterpmatrix(mat)
class(oft_matrix), pointer, intent(inout) :: mat !< Interpolation matrix
INTEGER(i4) :: i,j,k,m,icors,ifine,jb,i_ind(1),j_ind(1)
INTEGER(i4) :: etmp(2),ftmp(3),fetmp(3),ctmp(4),fc,ed
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pmap,emap,fmap
CLASS(oft_afem_type), POINTER :: hgrad_cors => NULL()
CLASS(oft_afem_type), POINTER :: hgrad_fine => NULL()
class(oft_mesh), pointer :: cmesh,mesh
CLASS(oft_vector), POINTER :: hgrad_vec,hgrad_vec_cors
integer(i4) :: jcp(4),jfe(3),jce(6)
integer(i4), pointer :: lcdg(:),lfde(:,:),lede(:,:),lcde(:,:)
real(r8) :: f(4),incr,val,d(3),h_rop(3),goptmp(3,4),v,mop(1)
type(oft_graph_ptr), pointer :: graphs(:,:)
type(oft_graph), pointer :: interp_graph
!---
if(ML_grad%ml_mesh%level<1)call oft_abort('Invalid mesh level','hgrad_ginterpmatrix',__FILE__)
mesh=>ML_grad%ml_mesh%mesh
cmesh=>ML_grad%ml_mesh%meshes(ML_grad%ml_mesh%level-1)
if(cmesh%type/=1)CALL oft_abort("Only supported with tet meshes", &
  "hgrad_ginterpmatrix", __FILE__)
hgrad_fine=>ML_grad%levels(ML_grad%level)%fe
if(hgrad_fine%order/=2)then
  call oft_abort('Attempted geometric interpolation for pd /= 2','hgrad_ginterpmatrix',__FILE__)
end if
!---
hgrad_cors=>ML_grad%levels(ML_grad%level-1)%fe
lede=>ML_grad%ml_mesh%inter(ML_grad%ml_mesh%level-1)%lede
lfde=>ML_grad%ml_mesh%inter(ML_grad%ml_mesh%level-1)%lfde
lcdg=>ML_grad%ml_mesh%inter(ML_grad%ml_mesh%level-1)%lcdg
lcde=>ML_grad%ml_mesh%inter(ML_grad%ml_mesh%level-1)%lcde
ALLOCATE(ML_grad%interp_graphs(ML_grad%level)%g)
interp_graph=>ML_grad%interp_graphs(ML_grad%level)%g
!---Setup matrix sizes
interp_graph%nr=hgrad_fine%ne
interp_graph%nrg=hgrad_fine%global%ne
interp_graph%nc=hgrad_cors%ne
interp_graph%ncg=hgrad_cors%global%ne
!---Setup Matrix graph
ALLOCATE(interp_graph%kr(interp_graph%nr+1))
interp_graph%kr=0
interp_graph%nnz=cmesh%np+5*cmesh%ne+3*cmesh%nf+6*cmesh%nc
ALLOCATE(interp_graph%lc(interp_graph%nnz))
interp_graph%lc=0_i4
!---Construct linkage
!$omp parallel do
DO i=1,cmesh%np
  interp_graph%kr(i)=1
END DO
!$omp parallel do
DO i=1,cmesh%ne
  interp_graph%kr(cmesh%np+i)=3
  interp_graph%kr(mesh%np+lede(1,i))=1
  interp_graph%kr(mesh%np+lede(2,i))=1
END DO
!$omp parallel do
DO i=1,cmesh%nf
  interp_graph%kr(mesh%np+lfde(1,i))=1
  interp_graph%kr(mesh%np+lfde(2,i))=1
  interp_graph%kr(mesh%np+lfde(3,i))=1
END DO
!$omp parallel do
DO i=1,cmesh%nc
  interp_graph%kr(mesh%np+ABS(lcde(1,i)))=6
END DO
interp_graph%kr(interp_graph%nr+1)=interp_graph%nnz+1
do i=interp_graph%nr,1,-1 ! cumulative point to point count
  interp_graph%kr(i)=interp_graph%kr(i+1)-interp_graph%kr(i)
end do
if(interp_graph%kr(1)/=1)call oft_abort('Bad element to element count','hgrad_ginterpmatrix',__FILE__)
!$omp parallel do private(k)
DO i=1,cmesh%np
  k=interp_graph%kr(i)
  interp_graph%lc(k)=i
END DO
!$omp parallel do private(k)
DO i=1,cmesh%ne
  !---Daughter point
  k=interp_graph%kr(cmesh%np+i)
  interp_graph%lc(k)=cmesh%le(1,i)
  interp_graph%lc(k+1)=cmesh%le(2,i)
  interp_graph%lc(k+2)=cmesh%np+i
  !---Daughter edge 1
  k=interp_graph%kr(mesh%np+lede(1,i))
  interp_graph%lc(k)=cmesh%np+i
  !---Daughter edge 2
  k=interp_graph%kr(mesh%np+lede(2,i))
  interp_graph%lc(k)=cmesh%np+i
END DO
!$omp parallel do private(k,j,jfe)
DO i=1,cmesh%nf
  jfe=ABS(cmesh%lfe(:,i)) ! face edges
  !---
  DO j=1,3
    k=interp_graph%kr(mesh%np+lfde(j,i))
    interp_graph%lc(k)=cmesh%np+jfe(j)
  END DO
END DO
!$omp parallel do private(k,j,jce)
DO i=1,cmesh%nc
  jce=ABS(cmesh%lce(:,i)) ! face edges
  CALL sort_array(jce,6)
  !---
  k=interp_graph%kr(mesh%np+ABS(lcde(1,i)))
  DO j=0,5
    interp_graph%lc(k+j)=cmesh%np+jce(1+j)
  END DO
END DO
!------------------------------------------------------------------------------
! Construct matrix
!------------------------------------------------------------------------------
NULLIFY(hgrad_vec,hgrad_vec_cors)
CALL ML_grad%vec_create(hgrad_vec)
CALL ML_grad%vec_create(hgrad_vec_cors,ML_hcurl_aug_rep%level-1)
!---
ALLOCATE(graphs(1,1))
graphs(1,1)%g=>interp_graph
!---
CALL create_matrix(mat,graphs,hgrad_vec,hgrad_vec_cors)
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Construct matrix
allocate(pmap(cmesh%np))
CALL get_inverse_map(cmesh%lbp,cmesh%nbp,pmap,cmesh%np)
DO i=1,cmesh%np
  IF(cmesh%bp(i))THEN
    ! IF(.NOT.cmesh%linkage%lpo(pmap(i)))CYCLE
    IF(.NOT.cmesh%pstitch%leo(pmap(i)))CYCLE
  END IF
  !---
  i_ind=i
  j_ind=i
  mop=1.d0
  CALL mat%add_values(i_ind,j_ind,mop,1,1)
END DO
deallocate(pmap)
!---
allocate(emap(cmesh%ne))
CALL get_inverse_map(cmesh%lbe,cmesh%nbe,emap,cmesh%ne)
DO i=1,cmesh%ne
  IF(cmesh%be(i))THEN
    ! IF(.NOT.cmesh%linkage%leo(emap(i)))CYCLE
    IF(.NOT.cmesh%estitch%leo(emap(i)))CYCLE
  END IF
  !---Daughter point
  i_ind=cmesh%np+i
  mop=.5d0
  DO j=1,2
    j_ind=cmesh%le(j,i)
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
  END DO
  mop=-.5d0
  j_ind=cmesh%np+i
  CALL mat%add_values(i_ind,j_ind,mop,1,1)
  !---Daughter edges
  DO j=1,2
    i_ind=mesh%np+lede(j,i)
    j_ind=cmesh%np+i
    mop=.25d0
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
  END DO
END DO
deallocate(emap)
!---
allocate(fmap(cmesh%nf))
CALL get_inverse_map(cmesh%lbf,cmesh%nbf,fmap,cmesh%nf)
DO i=1,cmesh%nf
  IF(cmesh%bf(i))THEN
    ! IF(.NOT.cmesh%linkage%lfo(fmap(i)))CYCLE
    IF(.NOT.cmesh%fstitch%leo(fmap(i)))CYCLE
  END IF
  !---
  jfe=ABS(cmesh%lfe(:,i)) ! face edges
  !---
  mop=.25d0
  DO j=1,3
    i_ind=mesh%np+lfde(j,i)
    j_ind=cmesh%np+jfe(j)
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
  END DO
END DO
deallocate(fmap)
!---
DO i=1,cmesh%nc ! loop over coarse cells
  jce=ABS(cmesh%lce(:,i)) ! cell edge indices
  !---
  i_ind=mesh%np+ABS(lcde(1,i))
  DO j=1,6
    j_ind=cmesh%np+jce(j)
    IF(j==lcdg(i).OR.j==lcdg(i)+3)THEN
      mop=-.25d0
    ELSE
      mop=.25d0
    END IF
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
  END DO
END DO
END SUBROUTINE hgrad_ginterpmatrix
END SUBROUTINE hcurl_grad_setup_interp
!------------------------------------------------------------------------------
!> Interpolate a coarse level Lagrange scalar field to the next finest level
!!
!! @note The global Lagrange level in incremented by one in this subroutine
!------------------------------------------------------------------------------
subroutine ml_vecspace_interp(self,acors,afine)
class(oft_ml_hcurl_grad_vecspace), intent(inout) :: self
class(oft_vector), intent(inout) :: acors !< Vector to interpolate
class(oft_vector), intent(inout) :: afine !< Fine vector from interpolation
integer(i4) :: i
real(r8), pointer, dimension(:) :: agrad,acurl,tmp
class(oft_mesh), pointer :: mesh
CLASS(oft_ml_fem_type), POINTER :: ML_curl
DEBUG_STACK_PUSH
!---Step one level up
call self%ML_FE_rep%set_level(self%ML_FE_rep%level+1,propogate=.TRUE.)
call afine%set(0.d0)
if(self%ML_FE_rep%level==self%ML_FE_rep%blevel+1)then
  IF(.NOT.ASSOCIATED(self%base_pop))CALL oft_abort("Base transfer not defined","ml_vecspace_interp",__FILE__)
  call self%base_pop(acors,afine)
  DEBUG_STACK_POP
  return
end if
CALL self%ML_FE_rep%interp_matrices(self%ML_FE_rep%level)%m%apply(acors,afine)
!---Correct gradient subspace following geometric interpolation
ML_curl=>self%ML_FE_rep%ml_fields(1)%ml
IF(ML_curl%current_level%order==1)THEN
  SELECT TYPE(this=>ML_curl%current_level)
    CLASS IS(oft_fem_type)
      mesh=>this%mesh
    CLASS DEFAULT
      CALL oft_abort("Invalid FE type","ml_vecspace_interp",__FILE__)
  END SELECT
  NULLIFY(agrad,acurl)
  CALL afine%get_local(acurl,1)
  CALL afine%get_local(agrad,2)
  ALLOCATE(tmp(mesh%np))
  !---
  !$omp parallel if(mesh%np>OFT_OMP_VTHRESH)
  !$omp do
  DO i=1,mesh%np
    tmp(i)=agrad(i)
    agrad(i)=0.d0
  END DO
  !$omp do
  DO i=1,mesh%ne
   acurl(i) = acurl(i) + &
    (tmp(mesh%le(2,i))-tmp(mesh%le(1,i)))*SIGN(1_i8,mesh%global%le(i))
  END DO
  !$omp end parallel
  !---
  CALL afine%restore_local(acurl,1,wait=.TRUE.)
  CALL afine%restore_local(agrad,2)
  DEALLOCATE(acurl,agrad,tmp)
END IF
DEBUG_STACK_POP
end subroutine ml_vecspace_interp
!------------------------------------------------------------------------------
!> Transfer a base level Lagrange scalar field to the next MPI level
!------------------------------------------------------------------------------
subroutine base_pop(self,acors,afine)
class(oft_ml_fe_comp_vecspace), intent(inout) :: self
class(oft_vector), intent(inout) :: acors !< Vector to transfer
class(oft_vector), intent(inout) :: afine !< Fine vector from transfer
integer(i4), pointer, dimension(:) :: lptmp
integer(i4), pointer, dimension(:) :: lbege
integer(i4) :: i
real(r8), pointer, dimension(:) :: array_c,array_f
CLASS(oft_ml_fem_type), POINTER :: ML_grad
DEBUG_STACK_PUSH
!---
ML_grad=>self%ML_FE_rep%ml_fields(2)%ml
NULLIFY(array_c,array_f)
lbege=>ML_grad%ml_mesh%inter(ML_grad%ml_mesh%nbase)%lbege
CALL acors%get_local(array_c,1)
CALL afine%get_local(array_f,1)
!$omp parallel do
do i=1,afine%n
  array_f(i)=array_c(ABS(lbege(i)))
end do
CALL afine%restore_local(array_f,1)
DEALLOCATE(array_c,array_f)
!---
lptmp=>ML_grad%ml_mesh%meshes(ML_grad%ml_mesh%nbase+1)%base%lp
CALL acors%get_local(array_c,2)
CALL afine%get_local(array_f,2)
!$omp parallel do
do i=1,afine%n
  array_f(i)=array_c(lptmp(i))
end do
CALL afine%restore_local(array_f,2)
DEALLOCATE(array_c,array_f)
DEBUG_STACK_POP
end subroutine base_pop
!------------------------------------------------------------------------------
!> Interpolate a coarse level Lagrange scalar field to the next finest level
!!
!! @note The global Lagrange level in incremented by one in this subroutine
!------------------------------------------------------------------------------
subroutine ml_vecspace_inject(self,afine,acors)
class(oft_ml_hcurl_grad_vecspace), intent(inout) :: self
class(oft_vector), intent(inout) :: afine !< Fine vector from interpolation
class(oft_vector), intent(inout) :: acors !< Vector to interpolate
DEBUG_STACK_PUSH
! Step down level down
call self%ML_FE_rep%set_level(self%ML_FE_rep%level-1,propogate=.TRUE.)
call acors%set(0.d0)
if(self%ML_FE_rep%level==self%ML_FE_rep%blevel)then
  IF(.NOT.ASSOCIATED(self%base_push))CALL oft_abort("Base transfer not defined","ml_vecspace_inject",__FILE__)
  call self%base_push(afine,acors)
  DEBUG_STACK_POP
  return
end if
CALL self%ML_FE_rep%interp_matrices(self%ML_FE_rep%level+1)%m%applyT(afine,acors)
DEBUG_STACK_POP
end subroutine ml_vecspace_inject
!------------------------------------------------------------------------------
!> Transfer a MPI level Lagrange scalar field to the base level
!------------------------------------------------------------------------------
subroutine base_push(self,afine,acors)
class(oft_ml_fe_comp_vecspace), intent(inout) :: self
class(oft_vector), intent(inout) :: afine !< Vector to transfer
class(oft_vector), intent(inout) :: acors !< Fine vector from transfer
integer(i4), pointer, dimension(:) :: lptmp
integer(i4), pointer, dimension(:) :: lbege
integer(i4) :: i,j,ierr
real(r8), pointer, dimension(:) :: alias,array_c,array_f
CLASS(oft_ml_fem_type), POINTER :: Ml_curl,ML_grad
CLASS(oft_hcurl_fem), POINTER :: curl_rep
CLASS(oft_h1_fem), POINTER :: grad_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_hcurl_cast(curl_rep,self%ML_FE_rep%ml_fields(1)%ml%current_level))CALL oft_abort("Incorrect Curl FE type","base_push",__FILE__)
IF(.NOT.oft_3D_h1_cast(grad_rep,self%ML_FE_rep%ml_fields(2)%ml%current_level))CALL oft_abort("Incorrect Grad FE type","base_push",__FILE__)
!---
NULLIFY(array_c,array_f)
Ml_curl=>self%ML_FE_rep%ml_fields(1)%ml
ML_grad=>self%ML_FE_rep%ml_fields(2)%ml
lbege=>ML_grad%ml_mesh%inter(ML_grad%ml_mesh%nbase)%lbege
lptmp=>ML_grad%ml_mesh%meshes(ML_grad%ml_mesh%nbase+1)%base%lp
CALL acors%get_local(array_c)
CALL afine%get_local(array_f)
!---
allocate(alias(acors%n))
alias=0.d0
!$omp parallel do
do i=1,curl_rep%ne
  if(curl_rep%linkage%be(i))cycle
  alias(ABS(lbege(i)))=array_f(i)
end do
!$omp parallel do private(j)
do i=1,curl_rep%linkage%nbe
  j=curl_rep%linkage%lbe(i)
  if(.NOT.curl_rep%linkage%leo(i))cycle
  alias(ABS(lbege(j)))=array_f(j)
end do
!$omp parallel do
do i=1,grad_rep%ne
  if(grad_rep%linkage%be(i))cycle
  alias(Ml_curl%levels(self%ML_FE_rep%level-1)%fe%ne+lptmp(i))=array_f(curl_rep%ne+i)
end do
!$omp parallel do private(j)
do i=1,grad_rep%linkage%nbe
  j=grad_rep%linkage%lbe(i)
  if(.NOT.grad_rep%linkage%leo(i))cycle
  alias(Ml_curl%levels(self%ML_FE_rep%level-1)%fe%ne+lptmp(j))=array_f(curl_rep%ne+j)
end do
!---Global reduction over all processors
array_c=oft_mpi_sum(alias,acors%n)
call acors%restore_local(array_c)
deallocate(alias,array_c,array_f)
DEBUG_STACK_POP
end subroutine base_push
!------------------------------------------------------------------------------
!> Compute eigenvalues and smoothing coefficients for the mass matrix
!------------------------------------------------------------------------------
SUBROUTINE hcurl_grad_mop_eigs(ML_hcurl_aug_obj,minlev)
type(oft_ml_fem_comp_type), target, intent(inout) :: ML_hcurl_aug_obj
INTEGER(i4), INTENT(in) :: minlev
#ifdef HAVE_ARPACK
INTEGER(i4) :: i
REAL(r8) :: lam0
REAL(r8), ALLOCATABLE :: df(:)
CLASS(oft_vector), POINTER :: u
TYPE(oft_irlm_eigsolver) :: arsolver
CLASS(oft_matrix), POINTER :: md => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
TYPE(oft_hcurl_grad_zerob), TARGET :: bc_tmp
DEBUG_STACK_PUSH
!------------------------------------------------------------------------------
! Compute optimal smoother coefficients
!------------------------------------------------------------------------------
IF(oft_env%head_proc)WRITE(*,*)'Optimizing Jacobi damping for H(Curl) + Grad(H^1) mass matrix'
bc_tmp%ML_hcurl_grad_rep=>ML_hcurl_aug_obj
ALLOCATE(df(ML_hcurl_aug_obj%nlevels))
df=0.d0
DO i=minlev,ML_hcurl_aug_obj%nlevels
  CALL ML_hcurl_aug_obj%set_level(i,propogate=.TRUE.)
  !---Create fields
  CALL ML_hcurl_aug_obj%vec_create(u)
  !---Get Ev range
  NULLIFY(mop)
  CALL hcurl_grad_getmop(ML_hcurl_aug_obj%current_level,mop,'lop')
  CALL create_diagmatrix(md,mop%D)
  !---
  arsolver%A=>mop
  arsolver%M=>md
  arsolver%mode=2
  arsolver%tol=1.E-5_r8
  arsolver%bc=>bc_tmp
  CALL create_native_pre(arsolver%Minv, "jacobi")
  arsolver%Minv%A=>mop
  !---
  IF(oft_debug_print(1))WRITE(*,*)'  optimizing level = ',i
  CALL arsolver%max(u,lam0)
  df(i) = 1.8d0/lam0
  !---
  CALL u%delete
  CALL md%delete
  CALL mop%delete
END DO
!---Output
IF(oft_env%head_proc)THEN
  WRITE(*,'(A)',ADVANCE='NO')' df_mop = '
  DO i=1,ML_hcurl_aug_obj%nlevels-1
    WRITE(*,'(F5.3,A)',ADVANCE='NO')df(i),', '
  END DO
  WRITE(*,'(F5.3,A)')df(ML_hcurl_aug_obj%nlevels)
END IF
DEALLOCATE(df)
DEBUG_STACK_POP
#else
CALL oft_abort("Subroutine requires ARPACK", "lag_lop_eigs", __FILE__)
#endif
END SUBROUTINE hcurl_grad_mop_eigs
!------------------------------------------------------------------------------
!> Compute eigenvalues and smoothing coefficients for the operator 
!! H(Curl) + Grad(H^1) mass matrix
!------------------------------------------------------------------------------
SUBROUTINE hcurl_grad_getmop_pre(ML_hcurl_aug_obj,pre,mats,level,nlevels)
type(oft_ml_fem_comp_type), target, intent(inout) :: ML_hcurl_aug_obj
CLASS(oft_solver), POINTER, INTENT(out) :: pre
TYPE(oft_matrix_ptr), POINTER, INTENT(inout) :: mats(:)
INTEGER(i4), OPTIONAL, INTENT(in) :: level
INTEGER(i4), OPTIONAL, INTENT(in) :: nlevels
!--- ML structures for MG-preconditioner
TYPE(oft_matrix_ptr), POINTER :: ml_int(:)
INTEGER(i4), ALLOCATABLE :: levels(:)
REAL(r8), ALLOCATABLE :: df(:)
INTEGER(i4), ALLOCATABLE :: nu(:)
INTEGER(i4) :: minlev,toplev,nl
INTEGER(i4) :: i,j,levin,ierr
LOGICAL :: create_mats
CHARACTER(LEN=2) :: lev_char
TYPE(xml_node), POINTER :: pre_node
TYPE(oft_ml_hcurl_grad_vecspace), POINTER :: tmp_vecspace
#ifdef HAVE_XML
integer(i4) :: nnodes
TYPE(xml_node), POINTER :: hcurl_grad_node
#endif
DEBUG_STACK_PUSH
!---
minlev=1
toplev=ML_hcurl_aug_obj%level
levin=ML_hcurl_aug_obj%level
IF(PRESENT(level))toplev=level
IF(PRESENT(nlevels))minlev=toplev-nlevels+1
nl=toplev-minlev+1
!---
IF(minlev<ML_hcurl_aug_obj%minlev)CALL oft_abort('Minimum level is < minlev','hcurl_grad_getmop_pre',__FILE__)
IF(toplev>ML_hcurl_aug_obj%nlevels)CALL oft_abort('Maximum level is > nlevels','hcurl_grad_getmop_pre',__FILE__)
!------------------------------------------------------------------------------
! Create ML Matrices
!------------------------------------------------------------------------------
create_mats=.FALSE.
IF(.NOT.ASSOCIATED(mats))THEN
  create_mats=.TRUE.
  ALLOCATE(mats(nl))
END IF
ALLOCATE(ml_int(nl-1),levels(nl))
ALLOCATE(df(nl),nu(nl))
DO i=1,nl
  CALL ML_hcurl_aug_obj%set_level(minlev+(i-1))
  levels(i)=minlev+(i-1)
  df(i)=df_mop(levels(i))
  nu(i)=nu_mop(levels(i))
  IF(df(i)<-1.d90)THEN
    WRITE(lev_char,'(I2.2)')levels(i)
    CALL oft_abort('Smoother values not set for level: '//lev_char,'hcurl_grad_getmop_pre',__FILE__)
  END IF
  !---
  IF(create_mats)THEN
    NULLIFY(mats(i)%M)
    CALL hcurl_grad_getmop(ML_hcurl_aug_obj%current_level,mats(i)%M,'none')
  END IF
  IF(i>1)ml_int(i-1)%M=>ML_hcurl_aug_obj%interp_matrices(ML_hcurl_aug_obj%level)%m
END DO
CALL ML_hcurl_aug_obj%set_level(levin,propogate=.TRUE.)
!------------------------------------------------------------------------------
! Search for XML-spec
!------------------------------------------------------------------------------
NULLIFY(pre_node)
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml))THEN
  CALL xml_get_element(oft_env%xml,"hcurl_grad",hcurl_grad_node,ierr)
  IF(ierr==0)CALL xml_get_element(hcurl_grad_node,"mop",pre_node,ierr)
END IF
#endif
!------------------------------------------------------------------------------
! Setup preconditioner
!------------------------------------------------------------------------------
NULLIFY(pre)
ALLOCATE(tmp_vecspace)
tmp_vecspace%ML_FE_rep=>ML_hcurl_aug_obj
tmp_vecspace%base_pop=>base_pop
tmp_vecspace%base_push=>base_push
CALL create_mlpre(pre,mats(1:nl),levels,nlevels=nl,ml_vecspace=tmp_vecspace, &
  stype=1,df=df,nu=nu,xml_root=pre_node)
!------------------------------------------------------------------------------
! Cleanup
!------------------------------------------------------------------------------
DEALLOCATE(ml_int,levels,df,nu)
DEBUG_STACK_POP
END SUBROUTINE hcurl_grad_getmop_pre
!------------------------------------------------------------------------------
!> Evaluate the jump error in a field over internal faces
!!
!! @note Currently faces on domain boundaries are skipped, this is due to the
!! fact that evaluting the error would require costly communication.
!!
!! @return Jump error metric
!------------------------------------------------------------------------------
FUNCTION hcurl_grad_jump_error(hcurl_grad_fe,u,quad_order) RESULT(error)
CLASS(oft_fem_comp_type), INTENT(inout) :: hcurl_grad_fe
CLASS(oft_vector), INTENT(inout) :: u !< H(Curl) + Grad(H^1) vector field to evaluate
INTEGER(i4), INTENT(in) :: quad_order !< Desired quadrature order for integration
REAL(r8) :: error,reg_jump,reg_energy,goptmp(3,3)
REAL(r8) :: gop(3,4),vol,area,val(3),f(4),norm(3)
REAL(r8), ALLOCATABLE :: Bn(:,:),rop_curl(:,:),rop_grad(:,:)
REAL(r8), POINTER :: curl_vals(:),grad_vals(:)
INTEGER(i4) :: i,j,m,cf,jc,ptmap(3)
INTEGER(i4) :: cell(2),face(2)
INTEGER(i4), ALLOCATABLE :: j_hcurl(:),j_hgrad(:)
TYPE(oft_quad_type) :: quad
CLASS(oft_bmesh), POINTER :: mesh_tmp
CLASS(oft_mesh), POINTER :: mesh
CLASS(oft_hcurl_fem), POINTER :: curl_rep
CLASS(oft_h1_fem), POINTER :: grad_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_hcurl_cast(curl_rep,hcurl_grad_fe%fields(1)%fe))CALL oft_abort("Incorrect Curl FE type","hcurl_grad_jump_error",__FILE__)
IF(.NOT.oft_3D_h1_cast(grad_rep,hcurl_grad_fe%fields(2)%fe))CALL oft_abort("Incorrect Grad FE type","hcurl_grad_jump_error",__FILE__)
!---Set quadrature order
ALLOCATE(oft_trimesh::mesh_tmp)
CALL mesh_tmp%quad_rule(quad_order, quad)
!---
NULLIFY(curl_vals,grad_vals)
CALL u%get_local(curl_vals,1)
CALL u%get_local(grad_vals,2)
!---Evalute integral
reg_jump=0.d0; reg_energy=0.d0
!!$omp parallel do private(curved,goptmp,v,pt,det,bcc) reduction(+:reg_jump) reduction(+:reg_energy)
ALLOCATE(Bn(2,quad%np))
ALLOCATE(j_hcurl(curl_rep%nce),j_hgrad(grad_rep%nce))
ALLOCATE(rop_curl(3,curl_rep%nce),rop_grad(3,grad_rep%nce))
!---
mesh=>curl_rep%mesh
mesh_tmp%np=3
mesh_tmp%nc=1
ALLOCATE(mesh_tmp%lc(3,1),mesh_tmp%r(3,3))
mesh_tmp%lc(:,1)=(/1,2,3/)
DO i=1,mesh%nf
  IF(mesh%lfc(2,i)==0)CYCLE
  !---Get border cells
  cell=mesh%lfc(:,i)
  DO j=1,2
    !---
    IF(j==1)THEN
      DO jc=1,3
        mesh_tmp%r(:,jc)=mesh%r(:,mesh%lf(jc,i))
      END DO
      CALL mesh_tmp%jacobian(1,quad%pts(:,1),goptmp,area)
      CALL mesh_tmp%norm(1,quad%pts(:,1),norm)
    ELSE
      norm=-norm
    END IF
    !---
    DO cf=1,4
      IF(ABS(mesh%lcf(cf,cell(j)))==i)EXIT
    END DO
    ptmap=(/1,2,3/)
    CALL orient_listn_inv(mesh%lcfo(cf,cell(j)),ptmap,3_i4)
    !---Get curl dofs
    call curl_rep%ncdofs(cell(j),j_hcurl) ! get curl DOFs
    call grad_rep%ncdofs(cell(j),j_hgrad) ! get grad DOFs
    !---Get local reconstructed operators
    DO m=1,quad%np ! Loop over quadrature points
      !---
      f=0.d0; f(mesh%cell_fc(:,cf))=quad%pts(ptmap,m)
      CALL mesh%jacobian(cell(j),f,gop,vol)
      !---
      val=0.d0
      CALL oft_hcurl_eval_all(curl_rep,cell(j),f,rop_curl,gop)
      CALL oft_h1_geval_all(grad_rep,cell(j),f,rop_grad,gop)
      DO jc=1,curl_rep%nce
        val=val+curl_vals(j_hcurl(jc))*rop_curl(:,jc)
      END DO
      DO jc=1,grad_rep%nce
        val=val+grad_vals(j_hgrad(jc))*rop_grad(:,jc)
      END DO
      !---
      Bn(j,m)=DOT_PRODUCT(val,norm)
      IF(j==2)THEN
        reg_jump = reg_jump + (SUM(Bn(:,m))**2)*area*quad%wts(m)
        reg_energy = reg_energy + SUM(val**2)*area*quad%wts(m)
      END IF
    END DO
  END DO
END DO
DEALLOCATE(Bn,j_hcurl,j_hgrad)
DEALLOCATE(rop_curl,rop_grad)
DEALLOCATE(mesh_tmp%lc,mesh_tmp%r)
DEALLOCATE(curl_vals,grad_vals)
!---Reduction over processors
reg_jump=oft_mpi_sum(reg_jump)
reg_energy=oft_mpi_sum(reg_energy)
error=SQRT(reg_jump/reg_energy)
!---Delete quadrature object
CALL quad%delete
DEBUG_STACK_POP
END FUNCTION hcurl_grad_jump_error
end module oft_hcurl_grad_operators
