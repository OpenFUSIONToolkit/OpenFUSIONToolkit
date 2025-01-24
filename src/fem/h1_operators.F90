!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_h1_operators.F90
!
!> Nedelec H1 FE operator definitions
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
!! @ingroup doxy_oft_h1
!---------------------------------------------------------------------------
MODULE oft_h1_operators
USE oft_base
USE oft_sort, ONLY: sort_array
USE oft_quadrature
USE oft_mesh_type, ONLY: oft_mesh, mesh, oft_bmesh, smesh, cell_is_curved
USE oft_trimesh_type, ONLY: oft_trimesh
USE multigrid, ONLY: mg_mesh
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, oft_matrix, oft_matrix_ptr, &
  oft_graph_ptr
USE oft_solver_base, ONLY: oft_solver, oft_orthog, oft_bc_proto
USE oft_deriv_matrices, ONLY: create_diagmatrix
#ifdef HAVE_ARPACK
USE oft_arpack, ONLY: oft_irlm_eigsolver
#endif
USE oft_la_utils, ONLY: create_matrix, combine_matrices
USE oft_solver_utils, ONLY: create_mlpre, create_cg_solver, create_diag_pre, &
  create_native_pre
USE fem_base, ONLY: oft_afem_type, oft_fem_type, fem_max_levels
USE fem_utils, ONLY: fem_interp
USE oft_lag_basis, ONLY: oft_lagrange, oft_lag_geval_all
USE oft_lag_fields, ONLY: oft_lag_create
USE oft_h0_basis, ONLY: oft_h0, oft_h0_ops, ML_oft_h0_ops, oft_h0_geval_all, &
  oft_h0_d2eval, oft_h0_fem
USE oft_h0_fields, ONLY: oft_h0_create
USE oft_h0_operators, ONLY: h0_zerob, h0_zerogrnd, oft_h0_getlop, oft_h0_gop
USE oft_hcurl_basis, ONLY: oft_hcurl, ML_oft_hcurl, oft_bhcurl, oft_hcurl_eval_all, &
oft_hcurl_ceval_all, oft_hcurl_level, oft_hcurl_blevel, oft_hcurl_get_cgops, &
oft_hcurl_fem
USE oft_hcurl_operators, ONLY: oft_hcurl_ops, oft_hcurl_rop, oft_hcurl_cop
USE oft_h1_basis, ONLY: oft_hgrad, ML_oft_hgrad, h1_ops, oft_h1_ops, ML_oft_h1_ops, &
oft_h1_level, oft_h1_blevel, oft_h1_nlevels, oft_h1_set_level, ML_oft_h1, oft_h1, &
oft_h1_minlev
USE oft_h1_fields, ONLY: oft_h1_create, oft_hgrad_create
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------
! CLASS oft_h1_rinterp
!---------------------------------------------------------------------------
!> Interpolate a H1 vector field
!---------------------------------------------------------------------------
type, extends(fem_interp) :: oft_h1_rinterp
  class(oft_vector), pointer :: u => NULL() !< Field to interpolate
  integer(i4), pointer, dimension(:) :: cache_cell => NULL()
  real(r8), pointer, dimension(:) :: grad_vals => NULL() !< Local gradient values
  real(r8), pointer, dimension(:) :: curl_vals => NULL() !< Local curl values
  real(r8), pointer, dimension(:,:) :: cache_grad => NULL()
  real(r8), pointer, dimension(:,:) :: cache_curl => NULL()
  class(oft_h0_fem), pointer :: hgrad_rep => NULL() !< H1(Grad) FE representation
  class(oft_hcurl_fem), pointer :: hcurl_rep => NULL() !< H1(Curl) FE representation
contains
  !> Retrieve local values for interpolation
  procedure :: setup => h1_rinterp_setup
  !> Reconstruct field
  procedure :: interp => h1_rinterp
  !> Destroy temporary internal storage
  procedure :: delete => h1_rinterp_delete
end type oft_h1_rinterp
!---------------------------------------------------------------------------
! CLASS oft_h1_cinterp
!---------------------------------------------------------------------------
!> Interpolate \f$ \nabla \times \f$ of a H1 vector field
!---------------------------------------------------------------------------
type, extends(oft_h1_rinterp) :: oft_h1_cinterp
contains
  !> Reconstruct field
  procedure :: interp => h1_cinterp
end type oft_h1_cinterp
!---------------------------------------------------------------------------
! CLASS oft_h1_dinterp
!---------------------------------------------------------------------------
!> Interpolate \f$ \nabla \times \f$ of a H1 vector field
!---------------------------------------------------------------------------
type, extends(oft_h1_rinterp) :: oft_h1_dinterp
contains
  !> Reconstruct field
  procedure :: interp => h1_dinterp
end type oft_h1_dinterp
!---------------------------------------------------------------------------
! CLASS oft_h1_divout
!---------------------------------------------------------------------------
!> Clean the divergence from a H1 vector field
!!
!! Divergence is removed by adding a gradient field, such that \f$ f = f +
!! \nabla \phi \f$, where \f$ \nabla^2 \phi = - \nabla \cdot f \f$. Cleaning
!! may also be applied to a field which is pre-multiplied by the H1::MOP
!! in which case \f$ \nabla^2 \phi = - \nabla^T f \f$ and \f$ f = f + M \nabla
!! \phi \f$. The mass matrix version is applied if \c mop is associated with the
!! corresponding H1::MOP.
!---------------------------------------------------------------------------
type, extends(oft_orthog) :: oft_h1_divout
  integer(i4) :: count = 0 !< Number of times apply has been called
  integer(i4) :: app_freq = 1 !< Frequency to apply solver
  logical :: keep_boundary = .FALSE. !< Flag for keeping boundary gradients
  logical :: pm = .FALSE. !< Flag for solver convergence monitor
  class(oft_solver), pointer :: solver => NULL() !< Solver object for H0::LOP operator
  procedure(oft_bc_proto), pointer, nopass :: bc => NULL() !< Boundary condition
  class(oft_matrix), pointer :: mop => NULL() !< Mass matrix, applies divoutm if associated
  class(oft_vector), pointer :: bnorm => NULL() !< Normal field source on boundary
contains
  !> Setup matrix and solver with default
  procedure :: setup => h1_divout_setup
  !> Clean divergence from field
  procedure :: apply => h1_divout_apply
  !> Clean-up internal storage
  procedure :: delete => h1_divout_delete
end type oft_h1_divout
!---------------------------------------------------------------------------
! CLASS oft_h1_zerograd
!---------------------------------------------------------------------------
!> Orthogonalize a H1 vector field by zeroing the gradient subspace
!---------------------------------------------------------------------------
type, extends(oft_orthog) :: oft_h1_zerograd
contains
  !> Perform orthoganlization
  procedure :: apply => h1_zerograd_apply
  !> Clean-up internal variables
  procedure :: delete => h1_zerograd_delete
end type oft_h1_zerograd
!---Pre options
integer(i4), private :: nu_mop(fem_max_levels)=0
real(r8), private :: df_mop(fem_max_levels)=-1.d99
contains
!---------------------------------------------------------------------------
! SUBROUTINE: h1_mloptions
!---------------------------------------------------------------------------
!> Read-in options for the basic Nedelec H1(Curl) ML preconditioners
!---------------------------------------------------------------------------
subroutine h1_mloptions
integer(i4) :: ierr,io_unit
namelist/h1_op_options/df_mop,nu_mop
DEBUG_STACK_PUSH
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,h1_op_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr>0)THEN
  CALL oft_abort('Error reading "h1_op_options" group in input file','h1_mloptions',__FILE__)
END IF
IF(df_mop(1)<-1.d90)THEN
  IF(oft_env%head_proc)THEN
    WRITE(*,*)'No H1 MG smoother settings found:'
    WRITE(*,*)'  Using default values, which may result in convergence failure.'
  END IF
  nu_mop=2
  df_mop=.2d0
END IF
DEBUG_STACK_POP
end subroutine h1_mloptions
!---------------------------------------------------------------------------
! SUBROUTINE: h1_rinterp_setup
!---------------------------------------------------------------------------
!> Setup interpolator for H1 vector fields
!!
!! Fetches local representation used for interpolation from vector object
!!
!! @note Should only be used via class \ref oft_h1_rinterp or children
!---------------------------------------------------------------------------
subroutine h1_rinterp_setup(self)
class(oft_h1_rinterp), intent(inout) :: self
IF(ASSOCIATED(self%parent))THEN
  SELECT TYPE(this=>self%parent)
    CLASS IS(oft_h1_rinterp)
      self%curl_vals=>this%curl_vals
      self%grad_vals=>this%grad_vals
      self%hgrad_rep=>this%hgrad_rep
      self%hcurl_rep=>this%hcurl_rep
      self%cache_cell=>this%cache_cell
      self%cache_grad=>this%cache_grad
      self%cache_curl=>this%cache_curl
    CLASS DEFAULT
      CALL oft_abort('Parent interpolator must be of same class','h1_rinterp_setup',__FILE__)
  END SELECT
ELSE
  !---Get local slice
  CALL self%u%get_local(self%curl_vals,1)
  CALL self%u%get_local(self%grad_vals,2)
  self%hgrad_rep=>oft_hgrad
  self%hcurl_rep=>oft_hcurl
  IF(.NOT.ASSOCIATED(self%cache_cell))THEN
    ALLOCATE(self%cache_cell(0:oft_env%nthreads-1))
    ALLOCATE(self%cache_grad(self%hgrad_rep%nce,0:oft_env%nthreads-1))
    ALLOCATE(self%cache_curl(self%hcurl_rep%nce,0:oft_env%nthreads-1))
  END IF
  self%cache_cell=-1
  self%cache_grad=0.d0
  self%cache_curl=0.d0
END IF
end subroutine h1_rinterp_setup
!---------------------------------------------------------------------------
! SUBROUTINE: h1_rinterp_delete
!---------------------------------------------------------------------------
!> Destroy temporary internal storage
!!
!! @note Should only be used via class \ref oft_h1_rinterp or children
!---------------------------------------------------------------------------
subroutine h1_rinterp_delete(self)
class(oft_h1_rinterp), intent(inout) :: self
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
end subroutine h1_rinterp_delete
!---------------------------------------------------------------------------
! SUBROUTINE: h1_rinterp
!---------------------------------------------------------------------------
!> Reconstruct a Nedelec H1 vector field
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [3]
!---------------------------------------------------------------------------
subroutine h1_rinterp(self,cell,f,gop,val)
class(oft_h1_rinterp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: rop(3)
real(r8), allocatable :: hgrad_rop(:,:)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%grad_vals))CALL oft_abort('Setup has not been called!','h1_rinterp',__FILE__)
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
IF(ASSOCIATED(oft_h0_gop))THEN
  do jc=1,self%hgrad_rep%nce
    val=val+self%cache_grad(jc,oft_tid)*oft_h0_gop(:,jc)
  end do
ELSE
  ALLOCATE(hgrad_rop(3,self%hgrad_rep%nce))
  CALL oft_h0_geval_all(self%hgrad_rep,cell,f,hgrad_rop,gop)
  do jc=1,self%hgrad_rep%nce
    val=val+self%cache_grad(jc,oft_tid)*hgrad_rop(:,jc)
  end do
  DEALLOCATE(hgrad_rop)
END IF
DEBUG_STACK_POP
end subroutine h1_rinterp
!---------------------------------------------------------------------------
! SUBROUTINE: h1_cinterp
!---------------------------------------------------------------------------
!> Reconstruct \f$ \nabla \times \f$ of a Nedelec H1 vector field
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [3]
!---------------------------------------------------------------------------
subroutine h1_cinterp(self,cell,f,gop,val)
class(oft_h1_cinterp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: cop(3),cgop(3,6)
real(r8), allocatable :: hcurl_cop(:,:)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%grad_vals))CALL oft_abort('Setup has not been called!','h1_cinterp',__FILE__)
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
end subroutine h1_cinterp
!---------------------------------------------------------------------------
! SUBROUTINE: h1_dinterp
!---------------------------------------------------------------------------
!> Reconstruct \f$ \nabla \cdot \f$ of a Nedelec H1 vector field
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [1]
!---------------------------------------------------------------------------
subroutine h1_dinterp(self,cell,f,gop,val)
class(oft_h1_dinterp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: dop(6),g2op(6,10),Kmat(10,3)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%grad_vals))CALL oft_abort('Setup has not been called!','h1_dinterp',__FILE__)
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
CALL mesh%hessian(cell,f,g2op,Kmat)
do jc=1,self%hgrad_rep%nce
  call oft_h0_d2eval(self%hgrad_rep,cell,jc,f,dop,g2op)
  val=val+self%cache_grad(jc,oft_tid)*(dop(1)+dop(4)+dop(6))
end do
DEBUG_STACK_POP
end subroutine h1_dinterp
!---------------------------------------------------------------------------
! SUBROUTINE: h1_zerob
!---------------------------------------------------------------------------
!> Zero a Nedelec H1 vector field at all boundary nodes
!!
!! @param[in,out] a Field to be zeroed
!---------------------------------------------------------------------------
subroutine h1_zerob(a)
class(oft_vector), intent(inout) :: a
real(r8), pointer, dimension(:) :: agrad,acurl
integer(i4) :: i,j
DEBUG_STACK_PUSH
!---Cast to vector type
NULLIFY(agrad,acurl)
CALL a%get_local(acurl,1)
CALL a%get_local(agrad,2)
! Apply operator
do i=1,oft_hcurl%nbe
  j=oft_hcurl%lbe(i)
  if(oft_hcurl%global%gbe(j))acurl(j)=0.d0
end do
! Apply operator
do i=1,oft_hgrad%nbe
  j=oft_hgrad%lbe(i)
  if(oft_hgrad%global%gbe(j))agrad(j)=0.d0
end do
CALL a%restore_local(acurl,1)
CALL a%restore_local(agrad,2)
DEALLOCATE(acurl,agrad)
DEBUG_STACK_POP
end subroutine h1_zerob
!---------------------------------------------------------------------------
! SUBROUTINE: h1curl_zerob
!---------------------------------------------------------------------------
!> Zero the curl components of a Nedelec H1 vector field at all boundary nodes
!!
!! @param[in,out] a Field to be zeroed
!---------------------------------------------------------------------------
subroutine h1curl_zerob(a)
class(oft_vector), intent(inout) :: a
real(r8), pointer, dimension(:) :: acurl
integer(i4) :: i,j
DEBUG_STACK_PUSH
!---Cast to vector type
NULLIFY(acurl)
CALL a%get_local(acurl,1)
! Apply operator
do i=1,oft_hcurl%nbe
  j=oft_hcurl%lbe(i)
  if(oft_hcurl%global%gbe(j))acurl(j)=0.d0
end do
CALL a%restore_local(acurl,1)
DEALLOCATE(acurl)
DEBUG_STACK_POP
end subroutine h1curl_zerob
!---------------------------------------------------------------------------
! SUBROUTINE: h1grad_zerop
!---------------------------------------------------------------------------
!> Zero the gradient components of a Nedelec H1 vector field at all boundary nodes
!!
!! @param[in,out] a Field to be zeroed
!---------------------------------------------------------------------------
subroutine h1grad_zerop(a)
class(oft_vector), intent(inout) :: a
real(r8), pointer, dimension(:) :: agrad
integer(i4) :: i,j
DEBUG_STACK_PUSH
!---Cast to vector type
NULLIFY(agrad)
CALL a%get_local(agrad,2)
! Apply operator
do i=1,mesh%np
  agrad(i)=0.d0
end do
CALL a%restore_local(agrad,2)
DEALLOCATE(agrad)
DEBUG_STACK_POP
end subroutine h1grad_zerop
!---------------------------------------------------------------------------
! SUBROUTINE: h1_zeroi
!---------------------------------------------------------------------------
!> Zero a Nedelec H1 vector field at all interior nodes
!!
!! @param[in,out] a Field to be zeroed
!---------------------------------------------------------------------------
subroutine h1_zeroi(a)
class(oft_vector), intent(inout) :: a
real(r8), pointer, dimension(:) :: agrad,acurl
integer(i4) :: i
DEBUG_STACK_PUSH
!---Cast to vector type
NULLIFY(agrad,acurl)
CALL a%get_local(acurl,1)
CALL a%get_local(agrad,2)
! Apply operator
DO i=1,oft_hcurl%ne
  IF(oft_hcurl%global%gbe(i))CYCLE
  acurl(i)=0.d0
END DO
! Apply operator
DO i=1,oft_hgrad%ne
  IF(oft_hgrad%global%gbe(i))CYCLE
  agrad(i)=0.d0
END DO
CALL a%restore_local(acurl,1)
CALL a%restore_local(agrad,2)
DEALLOCATE(acurl,agrad)
DEBUG_STACK_POP
end subroutine h1_zeroi
!---------------------------------------------------------------------------
! SUBROUTINE: h1_mc
!---------------------------------------------------------------------------
!> Compute the 0-th order gradient due to a jump plane
!!
!! The jump is represented as a circular surface defined by a center possition
!! and surface normal. This method requires that the mesh contain a
!! matching internal surface, such that no edge crosses the jump plane.
!!
!! @note The radius of the surface is represented by \f$ \left| hcpv \right| \f$
!!
!! @param[in,out] a Jump field
!! @param[in] hcpc Jump plane center possition [3]
!! @param[in] hcpv Jump plane normal vector [3]
!---------------------------------------------------------------------------
subroutine h1_mc(a,hcpc,hcpv,new_tol)
class(oft_vector), intent(inout) :: a
real(r8), intent(in) :: hcpc(3),hcpv(3)
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
end subroutine h1_mc
!---------------------------------------------------------------------------
! SUBROUTINE: h1_bmc
!---------------------------------------------------------------------------
!> Compute the 0-th order gradient due to a jump plane
!!
!! The jump is represented as a circular surface defined by a center possition
!! and surface normal. This method requires that the mesh contain a
!! matching internal surface, such that no edge crosses the jump plane.
!!
!! @note The radius of the surface is represented by \f$ \left| hcpv \right| \f$
!!
!! @param[in,out] a Jump field
!! @param[in] hcpc Jump plane center possition [3]
!! @param[in] hcpv Jump plane normal vector [3]
!---------------------------------------------------------------------------
subroutine h1_bmc(a,hcpc,hcpv,new_tol)
class(oft_vector), intent(inout) :: a
real(r8), intent(in) :: hcpc(3),hcpv(3)
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
  IF(ABS(r1dv)<=tol.AND.ABS(r2dv)<=tol)CALL oft_abort('Bad edge found','h1_bmc',__FILE__)
  IF((r1dv>=tol.AND.r2dv<=-tol).OR.(r1dv<=-tol.AND.r2dv>=tol))THEN
    r1 = (r1+r2)/2.d0
    r1cv=SUM(cross_product(r1,hcpv)**2)
    IF(r1cv<1.d0)acurl(i)=SIGN(1.d0,r2dv-r1dv)*SIGN(INT(1,8),mesh%global%le(i))
  END IF
END DO
CALL a%restore_local(acurl,1)
DEALLOCATE(acurl)
DEBUG_STACK_POP
end subroutine h1_bmc
!---------------------------------------------------------------------------
! SUBROUTINE: h1_grad
!---------------------------------------------------------------------------
!> Add the gradient of a H0 scalar field to a H1 vector field
!!
!! @note By default the 0-th order gradient subspace is represented on the
!! H1(Curl) DOF, use the \c keep_boundary flag otherwise.
!!
!! @param[in,out] a Scalar field
!! @param[in,out] b Vector field for gradient
!! @param[in] keep_boundary Flag to keep 0-th order boundary component (optional)
!---------------------------------------------------------------------------
subroutine h1_grad(a,b,keep_boundary)
class(oft_vector), intent(inout) :: a
class(oft_vector), intent(inout) :: b
logical, OPTIONAL, INTENT(in) :: keep_boundary
real(r8), pointer, dimension(:) :: aloc
real(r8), pointer, dimension(:) :: bgrad,bcurl
integer(i4), allocatable :: emap(:)
integer(i4) :: i,j,k,l
real(r8) :: reg
LOGICAL :: zero_boundary
DEBUG_STACK_PUSH
NULLIFY(aloc,bcurl,bgrad)
zero_boundary=.FALSE.
IF(PRESENT(keep_boundary))zero_boundary=keep_boundary
!---Cast to H1 type and get aliases
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
DO i=mesh%np+1,oft_hgrad%ne
  bgrad(i)=bgrad(i)+aloc(i)
END DO
!---Off-diagonal part
DO i=1,mesh%ne
  ! bcurl(i)=bcurl(i) + &
  ! (aloc(mesh%le(2,i))-aloc(mesh%le(1,i)))*SIGN(1_i8,mesh%global%le(i))
  bcurl((i-1)*oft_hcurl%gstruct(2)+1)=bcurl((i-1)*oft_hcurl%gstruct(2)+1) + &
    (aloc(mesh%le(2,i))-aloc(mesh%le(1,i)))*SIGN(1_i8,mesh%global%le(i))
END DO
!---
CALL b%restore_local(bcurl,1,wait=.TRUE.)
CALL b%restore_local(bgrad,2)
DEALLOCATE(aloc,bcurl,bgrad)
DEBUG_STACK_POP
end subroutine h1_grad
!---------------------------------------------------------------------------
! SUBROUTINE: h1_gradtp
!---------------------------------------------------------------------------
!> Apply the transposed gradient operator to a Nedelec H1 vector field.
!!
!! @param[in,out] a Input field
!! @param[in,out] b \f$ G^{T} a \f$
!---------------------------------------------------------------------------
subroutine h1_gradtp(a,b)
class(oft_vector), intent(inout) :: a
class(oft_vector), intent(inout) :: b
real(r8), pointer, dimension(:) :: agrad,acurl
real(r8), pointer, dimension(:) :: bloc
integer(i4), allocatable, dimension(:) :: emap
integer(i4) :: i,j,k,l
real(r8) :: reg
DEBUG_STACK_PUSH
NULLIFY(acurl,agrad,bloc)
!---Cast to H1 type and get aliases
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
do i=1,oft_h0%ne
  bloc(i)=bloc(i)+agrad(i)
end do
CALL b%restore_local(bloc,add=.TRUE.)
DEALLOCATE(bloc,acurl,agrad)
DEBUG_STACK_POP
end subroutine h1_gradtp
!---------------------------------------------------------------------------
! SUBROUTINE: h1_div
!---------------------------------------------------------------------------
!> Apply the divergence operator to a H1 field
!---------------------------------------------------------------------------
subroutine h1_div(a,b)
class(oft_vector), intent(inout) :: a
class(oft_vector), intent(inout) :: b
integer(i4) :: i,m,jr
integer(i4), allocatable :: j_curl(:),j_grad(:)
real(r8) :: vol,det,goptmp(3,4),aloc(3)
real(r8), pointer, contiguous, dimension(:) :: ac_loc,ag_loc,bloc
real(r8), allocatable, dimension(:) :: ac_tmp,ag_tmp,btmp
real(r8), allocatable, dimension(:,:) :: rop_curl,rop_grad
logical :: curved
DEBUG_STACK_PUSH
!---------------------------------------------------------------------------
! Allocate Laplacian Op
!---------------------------------------------------------------------------
NULLIFY(ac_loc,ag_loc,bloc)
!---Cast to H1 type and get aliases
CALL a%get_local(ac_loc,1)
CALL a%get_local(ag_loc,2)
!---Zero output and get aliases
CALL b%set(0.d0)
CALL b%get_local(bloc)
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j_curl,j_grad,rop_curl,rop_grad,det,btmp,ac_tmp,ag_tmp, &
!$omp aloc,curved,goptmp,vol,m,jr)
allocate(j_curl(oft_hcurl%nce),j_grad(oft_hgrad%nce))
allocate(rop_grad(3,oft_hgrad%nce),rop_curl(3,oft_hcurl%nce))
allocate(ac_tmp(oft_hcurl%nce),ag_tmp(oft_hgrad%nce))
allocate(btmp(oft_hgrad%nce))
!$omp do schedule(guided)
do i=1,mesh%nc
  !---Get local to global DOF mapping
  call oft_hcurl%ncdofs(i,j_curl)
  call oft_hgrad%ncdofs(i,j_grad)
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  do jr=1,oft_hcurl%nce
    ac_tmp(jr)=ac_loc(j_curl(jr))
  end do
  do jr=1,oft_hgrad%nce
    ag_tmp(jr)=ag_loc(j_grad(jr))
  end do
  btmp=0.d0
  do m=1,oft_hgrad%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call mesh%jacobian(i,oft_hgrad%quad%pts(:,m),goptmp,vol)
    det=vol*oft_hgrad%quad%wts(m)
    call oft_hcurl_eval_all(oft_hcurl,i,oft_hgrad%quad%pts(:,m),rop_curl,goptmp)
    call oft_h0_geval_all(oft_hgrad,i,oft_hgrad%quad%pts(:,m),rop_grad,goptmp)
    !---Compute local operator contribution
    aloc = 0.d0
    do jr=1,oft_hcurl%nce
      aloc = aloc + rop_curl(:,jr)*ac_tmp(jr)
    end do
    do jr=1,oft_hgrad%nce
      aloc = aloc + rop_grad(:,jr)*ag_tmp(jr)
    end do
    do jr=1,oft_hgrad%nce
      btmp(jr)=btmp(jr)+DOT_PRODUCT(rop_grad(:,jr),aloc)*det
    end do
  end do
  !---Add local values to global vector
  do jr=1,oft_hgrad%nce
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
end subroutine h1_div
!---------------------------------------------------------------------------
! SUBROUTINE: h1_curltp
!---------------------------------------------------------------------------
!> Apply the curl transpose operator to a H1 field
!---------------------------------------------------------------------------
subroutine h1_curltp(a,b)
class(oft_vector), intent(inout) :: a
class(oft_vector), intent(inout) :: b
real(r8), pointer, dimension(:) :: agrad,acurl
real(r8), pointer, dimension(:) :: bcurl
integer(i4) :: i,jr,jc,m
integer(i4), allocatable :: j_curl(:),j_grad(:)
real(r8) :: goptmp(3,4),cgop(3,6)
real(r8) :: v,f(4),det,vol
logical :: curved
real(r8), allocatable :: rop_curl(:,:),cop_curl(:,:),rop_grad(:,:),btmp(:)
DEBUG_STACK_PUSH
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
NULLIFY(acurl,agrad,bcurl)
!---Cast to H1 type and get aliases
CALL a%get_local(acurl,1)
CALL a%get_local(agrad,2)
!---Zero output and get aliases
CALL b%set(0.d0)
CALL b%get_local(bcurl,1)
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j_curl,j_grad,rop_curl,cop_curl,rop_grad,det,btmp,curved,goptmp,cgop,m,v,jc,jr,f)
allocate(j_curl(oft_hcurl%nce)) ! Local DOF and matrix indices
allocate(j_grad(oft_hgrad%nce)) ! Local DOF and matrix indices
allocate(rop_grad(3,oft_hgrad%nce)) ! Reconstructed gradient operator
allocate(rop_curl(3,oft_hcurl%nce)) ! Reconstructed gradient operator
allocate(cop_curl(3,oft_hcurl%nce)) ! Reconstructed gradient operator
allocate(btmp(oft_hcurl%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,mesh%nc
  !---Get local to global DOF mapping
  call oft_hcurl%ncdofs(i,j_curl)
  call oft_hgrad%ncdofs(i,j_grad)
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  btmp=0.d0
  do m=1,oft_hcurl%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)then
      call mesh%jacobian(i,oft_hcurl%quad%pts(:,m),goptmp,v)
      CALL oft_hcurl_get_cgops(goptmp,cgop)
    end if
    det=v*oft_hcurl%quad%wts(m)
    call oft_hcurl_eval_all(oft_hcurl,i,oft_hcurl%quad%pts(:,m),rop_curl,goptmp)
    call oft_hcurl_ceval_all(oft_hcurl,i,oft_hcurl%quad%pts(:,m),cop_curl,cgop)
    call oft_h0_geval_all(oft_hgrad,i,oft_hcurl%quad%pts(:,m),rop_grad,goptmp)
    !---Compute local operator contribution
    do jr=1,oft_hcurl%nce
      do jc=1,oft_hcurl%nce
        btmp(jr)=btmp(jr)+DOT_PRODUCT(cop_curl(:,jr),rop_curl(:,jc))*acurl(j_curl(jc))*det
      end do
      do jc=1,oft_hgrad%nce
        btmp(jr)=btmp(jr)+DOT_PRODUCT(cop_curl(:,jr),rop_grad(:,jc))*agrad(j_grad(jc))*det
      end do
    end do
  end do
  !---Add local values to global vector
  do jr=1,oft_hcurl%nce
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
end subroutine h1_curltp
!---------------------------------------------------------------------------
! SUBROUTINE: h1_getmop
!---------------------------------------------------------------------------
!> Construct mass matrix for a H1 representation
!!
!! Supported boundary conditions
!! - \c 'none' Full matrix
!! - \c 'zerob' Dirichlet for all boundary DOF
!!
!! @param[in,out] mat Matrix object
!! @param[in] bc Boundary condition
!---------------------------------------------------------------------------
subroutine h1_getmop(mat,bc)
class(oft_matrix), pointer, intent(inout) :: mat
character(LEN=*), intent(in) :: bc
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j_curl(:),j_grad(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time
real(r8), allocatable, dimension(:,:) :: rop_curl,rop_grad
real(r8), allocatable, dimension(:,:) :: mop11,mop12,mop21,mop22
logical :: curved
CLASS(oft_vector), POINTER :: oft_h1_vec
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing H1::MOP'
  CALL mytimer%tick()
END IF
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL oft_h1%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!---------------------------------------------------------------------------
! Operator integration
!---------------------------------------------------------------------------
!$omp parallel private(j_curl,j_grad,rop_curl,rop_grad,det,mop11,mop12,mop21,mop22, &
!$omp curved,goptmp,m,vol,jc,jr)
allocate(j_curl(oft_hcurl%nce)) ! Local DOF and matrix indices
allocate(j_grad(oft_hgrad%nce)) ! Local DOF and matrix indices
allocate(rop_grad(3,oft_hgrad%nce)) ! Reconstructed gradient operator
allocate(rop_curl(3,oft_hcurl%nce)) ! Reconstructed gradient operator
allocate(mop11(oft_hcurl%nce,oft_hcurl%nce)) ! Local laplacian matrix
allocate(mop12(oft_hcurl%nce,oft_hgrad%nce)) ! Local laplacian matrix
allocate(mop21(oft_hgrad%nce,oft_hcurl%nce)) ! Local laplacian matrix
allocate(mop22(oft_hgrad%nce,oft_hgrad%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  mop11=0.d0; mop12=0.d0
  mop21=0.d0; mop22=0.d0
  do m=1,oft_hcurl%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call mesh%jacobian(i,oft_hcurl%quad%pts(:,m),goptmp,vol)
    det=vol*oft_hcurl%quad%wts(m)
    call oft_hcurl_eval_all(oft_hcurl,i,oft_hcurl%quad%pts(:,m),rop_curl,goptmp)
    call oft_h0_geval_all(oft_hgrad,i,oft_hcurl%quad%pts(:,m),rop_grad,goptmp)
    !---Compute local matrix contributions
    do jr=1,oft_hcurl%nce
      do jc=1,oft_hcurl%nce
        mop11(jr,jc) = mop11(jr,jc) + DOT_PRODUCT(rop_curl(:,jr),rop_curl(:,jc))*det
      end do
      do jc=1,oft_hgrad%nce
        mop12(jr,jc) = mop12(jr,jc) + DOT_PRODUCT(rop_curl(:,jr),rop_grad(:,jc))*det
      end do
    end do
    do jr=1,oft_hgrad%nce
      do jc=1,oft_hgrad%nce
        mop22(jr,jc) = mop22(jr,jc) + DOT_PRODUCT(rop_grad(:,jr),rop_grad(:,jc))*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call oft_hcurl%ncdofs(i,j_curl)
  call oft_hgrad%ncdofs(i,j_grad)
  mop21=TRANSPOSE(mop12)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("none")
      mop21(1:mesh%cell_np,:)=0.d0
      mop22(1:mesh%cell_np,:)=0.d0
    CASE("zerob")
      DO jr=1,oft_hcurl%nce
        IF(oft_hcurl%global%gbe(j_curl(jr)))THEN
          mop11(jr,:)=0.d0
          mop12(jr,:)=0.d0
        END IF
      END DO
      mop21(1:mesh%cell_np,:)=0.d0
      mop22(1:mesh%cell_np,:)=0.d0
      DO jr=mesh%cell_np,oft_hgrad%nce
        IF(oft_hgrad%global%gbe(j_grad(jr)))THEN
          mop12(jr,:)=0.d0
          mop22(jr,:)=0.d0
        END IF
      END DO
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical (h1_getmop_1)
  call mat%atomic_add_values(j_curl,j_curl,mop11,oft_hcurl%nce,oft_hcurl%nce,1,1)
  call mat%atomic_add_values(j_curl,j_grad,mop12,oft_hcurl%nce,oft_hgrad%nce,1,2)
  ! !$omp end critical (h1_getmop_1)
  ! !$omp critical (h1_getmop_2)
  call mat%atomic_add_values(j_grad,j_curl,mop21,oft_hgrad%nce,oft_hcurl%nce,2,1)
  call mat%atomic_add_values(j_grad,j_grad,mop22,oft_hgrad%nce,oft_hgrad%nce,2,2)
  ! !$omp end critical (h1_getmop_2)
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
    DO i=1,mesh%np
      IF(mesh%bp(i))CYCLE
      j_curl=i
      call mat%add_values(j_curl,j_curl,mop11,1,1,2,2)
    END DO
    DO i=1,mesh%nbp
      jr=mesh%lbp(i)
      ! IF(.NOT.mesh%linkage%lpo(i))CYCLE
      IF(.NOT.mesh%pstitch%leo(i))CYCLE
      j_curl=jr
      call mat%add_values(j_curl,j_curl,mop11,1,1,2,2)
    END DO
  CASE("zerob")
    mop11(1,1)=1.d0
    DO i=1,oft_hcurl%nbe
      jr=oft_hcurl%lbe(i)
      IF(.NOT.oft_hcurl%global%gbe(jr))CYCLE
      IF(.NOT.oft_hcurl%linkage%leo(i))CYCLE
      j_curl=jr
      call mat%add_values(j_curl,j_curl,mop11,1,1,1,1)
    END DO
    DO i=1,mesh%np
      IF(.NOT.oft_hgrad%linkage%be(i))CYCLE
      j_curl=i
      call mat%add_values(j_curl,j_curl,mop11,1,1,1,1)
    END DO
    DO i=1,oft_hgrad%nbe
      jr=oft_hgrad%lbe(i)
      IF(.NOT.oft_hgrad%global%gbe(jr))CYCLE
      IF(.NOT.oft_hgrad%linkage%leo(i))CYCLE
      j_curl=jr
      call mat%add_values(j_curl,j_curl,mop11,1,1,1,1)
    END DO
END SELECT
DEALLOCATE(j_curl,mop11)
CALL oft_h1_create(oft_h1_vec)
CALL mat%assemble(oft_h1_vec)
CALL oft_h1_vec%delete
DEALLOCATE(oft_h1_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine h1_getmop
!---------------------------------------------------------------------------
! SUBROUTINE: oft_h1_project
!---------------------------------------------------------------------------
!> Project a vector field onto a H1 basis
!!
!! @note This subroutine only performs the integration of the field with
!! test functions for a H1 basis. To retrieve the correct projection the
!! result must be multiplied by the inverse of H1::MOP.
!!
!! @param[in,out] field Vector field for projection
!! @param[in,out] x Field projected onto H1 basis
!---------------------------------------------------------------------------
subroutine oft_h1_project(field,x)
class(fem_interp), intent(inout) :: field
class(oft_vector), intent(inout) :: x
integer(i4) :: i,jc,m
integer(i4), allocatable, dimension(:) :: j_hcurl,j_hgrad
real(r8) :: det,vol,bcc(3),goptmp(3,4)
real(r8), pointer, dimension(:) :: xcurl,xgrad
real(r8), allocatable, dimension(:,:) :: rop_curl,rop_grad
logical :: curved
DEBUG_STACK_PUSH
!---Initialize vectors to zero
NULLIFY(xcurl,xgrad)
call x%set(0.d0)
call x%get_local(xcurl,1)
call x%get_local(xgrad,2)
!---Integerate over the volume
!$omp parallel default(firstprivate) shared(xcurl,xgrad) private(curved,det)
allocate(j_hcurl(oft_hcurl%nce),rop_curl(3,oft_hcurl%nce))
allocate(j_hgrad(oft_hgrad%nce),rop_grad(3,oft_hgrad%nce))
!$omp do schedule(guided)
do i=1,mesh%nc ! Loop over cells
  call oft_hcurl%ncdofs(i,j_hcurl) ! Get DOFs
  call oft_hgrad%ncdofs(i,j_hgrad) ! Get DOFs
  curved=cell_is_curved(mesh,i) ! Straight cell test
  do m=1,oft_hcurl%quad%np
    if(curved.OR.m==1)call mesh%jacobian(i,oft_hcurl%quad%pts(:,m),goptmp,vol)
    det=vol*oft_hcurl%quad%wts(m)
    call field%interp(i,oft_hcurl%quad%pts(:,m),goptmp,bcc)
    call oft_hcurl_eval_all(oft_hcurl,i,oft_hcurl%quad%pts(:,m),rop_curl,goptmp)
    call oft_h0_geval_all(oft_hgrad,i,oft_hcurl%quad%pts(:,m),rop_grad,goptmp)
    do jc=1,oft_hcurl%nce
      !$omp atomic
      xcurl(j_hcurl(jc))=xcurl(j_hcurl(jc))+DOT_PRODUCT(rop_curl(:,jc),bcc)*det
    end do
    do jc=1,oft_hgrad%nce
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
end subroutine oft_h1_project
!------------------------------------------------------------------------------
! SUBROUTINE: oft_h1_bproject
!------------------------------------------------------------------------------
!> Boundary projection of a vector field onto a H1 basis
!!
!! @note This subroutine only performs the integration of the field with
!! boundary test functions for a H1 basis.
!!
!! @param[in,out] field Vector field for projection
!! @param[in,out] x Field projected onto H1 basis
!------------------------------------------------------------------------------
SUBROUTINE oft_h1_bproject(field,x)
CLASS(fem_interp), INTENT(inout) :: field
CLASS(oft_vector), INTENT(inout) :: x
INTEGER(i4) :: i,m,jc,cf,face,cell,ptmap(3)
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: j_hcurl,j_hgrad
REAL(r8) :: vol,det,f(4),norm(3),etmp(3),goptmp(3,3),gop(3,4)
REAL(r8), POINTER, DIMENSION(:) :: xcurl,xgrad
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: rop_curl,rop_grad
DEBUG_STACK_PUSH
!---Initialize vectors to zero
NULLIFY(xcurl,xgrad)
call x%set(0.d0)
call x%get_local(xcurl,1)
call x%get_local(xgrad,2)
!---Operator integration loop
!$omp parallel default(firstprivate) shared(xcurl,xgrad) private(det)
allocate(j_hcurl(oft_hcurl%nce),rop_curl(3,oft_hcurl%nce))
allocate(j_hgrad(oft_hgrad%nce),rop_grad(3,oft_hgrad%nce))
!$omp do schedule(guided)
do i=1,smesh%nc
  CALL mesh%get_surf_map(i,cell,ptmap) ! Find parent cell and logical coordinate mapping
  !---Get local to global DOF mapping
  call oft_hcurl%ncdofs(cell,j_hcurl)
  call oft_hgrad%ncdofs(cell,j_hgrad)
  !---Get local reconstructed operators
  do m=1,oft_bhcurl%quad%np ! Loop over quadrature points
    call smesh%jacobian(i,oft_bhcurl%quad%pts(:,m),goptmp,vol)
    call smesh%norm(i,oft_bhcurl%quad%pts(:,m),norm)
    det=vol*oft_bhcurl%quad%wts(m)
    !---
    CALL mesh%surf_to_vol(oft_bhcurl%quad%pts(:,m),ptmap,f)
    call mesh%jacobian(cell,f,gop,vol)
    call field%interp(cell,f,gop,etmp)
    call oft_hcurl_eval_all(oft_hcurl,cell,f,rop_curl,gop)
    call oft_h0_geval_all(oft_hgrad,cell,f,rop_grad,gop)
    do jc=1,oft_hcurl%nce
      !$omp atomic
      xcurl(j_hcurl(jc))=xcurl(j_hcurl(jc)) &
      + DOT_PRODUCT(cross_product(rop_curl(:,jc),etmp),norm)*det
    end do
    do jc=1,oft_hgrad%nce
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
END SUBROUTINE oft_h1_bproject
!---------------------------------------------------------------------------
! SUBROUTINE: h1_divout_setup
!---------------------------------------------------------------------------
!> Setup matrix and solver with default
!!
!! @note Should only be used via class \ref oft_h1_divout
!!
!! @param[in] bc Boundary condition
!---------------------------------------------------------------------------
subroutine h1_divout_setup(self,bc)
class(oft_h1_divout), intent(inout) :: self
character(LEN=*), intent(in) :: bc
!---
CLASS(oft_matrix), POINTER :: lop
CLASS(oft_solver), POINTER :: linv
DEBUG_STACK_PUSH
NULLIFY(lop)
CALL oft_h0_getlop(lop,bc)
CALL create_cg_solver(linv)
linv%A=>lop
linv%its=-3
CALL create_diag_pre(linv%pre)
self%solver=>linv
IF(TRIM(bc)=='grnd')THEN
  self%bc=>h0_zerogrnd
ELSE
  self%bc=>h0_zerob
END IF
DEBUG_STACK_POP
end subroutine h1_divout_setup
!---------------------------------------------------------------------------
! SUBROUTINE: h1_divout_apply
!---------------------------------------------------------------------------
!> Remove divergence from a H1 vector field by adding a gradient correction
!!
!! @note Should only be used via class \ref oft_h1_divout
!!
!! @param[in,out] u Field for divergence cleaning
!---------------------------------------------------------------------------
subroutine h1_divout_apply(self,u)
class(oft_h1_divout), intent(inout) :: self
class(oft_vector), intent(inout) :: u
class(oft_vector), pointer :: a,g,tmp,tmp2
integer(i4) :: i,order_tmp
real(r8) :: uu
logical :: pm_save
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%solver))CALL oft_abort('No solver specified.','h1_divout_apply',__FILE__)
IF(mod(self%count,self%app_freq)/=0)THEN
  DEBUG_STACK_POP
  RETURN
END IF
!---
call oft_h0_create(g)
call oft_h0_create(a)
!---
IF(ASSOCIATED(self%mop))THEN
  CALL h1_gradtp(u,g)
ELSE
  CALL h1_div(u,g)
END IF
uu=u%dot(u)
self%solver%atol=MAX(self%solver%atol,SQRT(uu*1.d-20))
call a%set(0.d0)
IF(ASSOCIATED(self%bnorm))CALL g%add(1.d0,-1.d0,self%bnorm)
call self%bc(g)
!---
pm_save=oft_env%pm; oft_env%pm=self%pm
call self%solver%apply(a,g)
oft_env%pm=pm_save
!---
CALL a%scale(-1.d0)
CALL u%new(tmp)
CALL h1_grad(a,tmp,self%keep_boundary)
IF(ASSOCIATED(self%mop))THEN
  CALL u%new(tmp2)
  CALL self%mop%apply(tmp,tmp2)
  CALL u%add(1.d0,1.d0,tmp2)
  CALL tmp2%delete()
  DEALLOCATE(tmp2)
ELSE
  CALL u%add(1.d0,1.d0,tmp)
END IF
CALL tmp%delete
call a%delete
call g%delete
DEALLOCATE(tmp,a,g)
DEBUG_STACK_POP
end subroutine h1_divout_apply
!---------------------------------------------------------------------------
! SUBROUTINE: h1_divout_delete
!---------------------------------------------------------------------------
!> Clean-up internal storage for a oft_h1_divout object
!!
!! @note Should only be used via class \ref oft_h1_divout
!---------------------------------------------------------------------------
subroutine h1_divout_delete(self)
class(oft_h1_divout), intent(inout) :: self
IF(ASSOCIATED(self%solver))THEN
  call self%solver%A%delete
  deallocate(self%solver%A)
  call self%solver%pre%delete
  deallocate(self%solver%pre)
  call self%solver%delete
  deallocate(self%solver)
END IF
end subroutine h1_divout_delete
!---------------------------------------------------------------------------
! SUBROUTINE: h1_zerograd_apply
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
subroutine h1_zerograd_apply(self,u)
class(oft_h1_zerograd), intent(inout) :: self
class(oft_vector), intent(inout) :: u
real(r8), pointer, dimension(:) :: ugrad
DEBUG_STACK_PUSH
!---Get local field
NULLIFY(ugrad)
CALL u%get_local(ugrad,2)
ugrad=0.d0
CALL u%restore_local(ugrad,2)
DEALLOCATE(ugrad)
DEBUG_STACK_POP
end subroutine h1_zerograd_apply
!---------------------------------------------------------------------------
! SUBROUTINE: h1_zerograd_delete
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
subroutine h1_zerograd_delete(self)
class(oft_h1_zerograd), intent(inout) :: self

end subroutine h1_zerograd_delete
!---------------------------------------------------------------------------
! SUBROUTINE: h1_setup_interp
!---------------------------------------------------------------------------
!> Construct interpolation matrices on each MG level
!---------------------------------------------------------------------------
SUBROUTINE h1_setup_interp(create_full)
LOGICAL, OPTIONAL, INTENT(in) :: create_full
TYPE(oft_graph_ptr) :: graphs(2,2)
TYPE(oft_matrix_ptr) :: mats(2,2)
CLASS(oft_vector), POINTER :: fvec,cvec
INTEGER(i4) :: i
LOGICAL :: full_interp
DEBUG_STACK_PUSH
full_interp=.FALSE.
IF(PRESENT(create_full))full_interp=create_full
!---
DO i=oft_h1_minlev+1,oft_h1_nlevels
  CALL oft_h1_set_level(i)
  !---
  IF(oft_h1_level==oft_h1_blevel+1)THEN
    CYCLE
  END IF
  !---Setup interpolation operators for H1(Grad) space
  IF(oft_hcurl%order==1)THEN
    CALL hgrad_ginterpmatrix(ML_oft_hgrad%interp_matrices(ML_oft_hgrad%level)%m)
    oft_h1_ops%hgrad_interp=>ML_oft_hgrad%interp_matrices(ML_oft_hgrad%level)%m
    CALL oft_h1_ops%hgrad_interp%assemble
  ELSE
    ML_oft_hgrad%interp_graphs(ML_oft_hgrad%level)%g=>ML_oft_h0_ops(oft_h1_level+1)%interp_graph
    ML_oft_hgrad%interp_matrices(ML_oft_hgrad%level)%m=>ML_oft_h0_ops(oft_h1_level+1)%interp
    oft_h1_ops%hgrad_interp_graph=>ML_oft_h0_ops(oft_h1_level+1)%interp_graph
    oft_h1_ops%hgrad_interp=>ML_oft_h0_ops(oft_h1_level+1)%interp
  END IF
END DO
!---Create full H1 interpolation operator
IF(full_interp)THEN
  CALL ML_oft_h1%build_interp
  DO i=oft_h1_minlev+1,oft_h1_nlevels
    CALL oft_h1_set_level(i)
    !---
    if(oft_h1_level==oft_h1_blevel+1)then
      CYCLE
    end if
    oft_h1_ops%interp=>ML_oft_h1%interp_matrices(i)%m
  END DO
END IF
DEBUG_STACK_POP
END SUBROUTINE h1_setup_interp
!---------------------------------------------------------------------------
! SUBROUTINE: hgrad_ginterpmatrix
!---------------------------------------------------------------------------
!> Construct interpolation matrix for polynomial levels
!---------------------------------------------------------------------------
SUBROUTINE hgrad_ginterpmatrix(mat)
class(oft_matrix), pointer, intent(inout) :: mat
INTEGER(i4) :: i,j,k,m,icors,ifine,jb,i_ind(1),j_ind(1)
INTEGER(i4) :: etmp(2),ftmp(3),fetmp(3),ctmp(4),fc,ed
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pmap,emap,fmap
CLASS(oft_afem_type), POINTER :: hgrad_cors => NULL()
TYPE(h1_ops), POINTER :: ops
class(oft_mesh), pointer :: cmesh
CLASS(oft_vector), POINTER :: hgrad_vec,hgrad_vec_cors
integer(i4) :: jcp(4),jfe(3),jce(6)
integer(i4), pointer :: lcdg(:),lfde(:,:),lede(:,:),lcde(:,:)
real(r8) :: f(4),incr,val,d(3),h_rop(3),goptmp(3,4),v,mop(1)
type(oft_graph_ptr), pointer :: graphs(:,:)
DEBUG_STACK_PUSH
!---
if(mg_mesh%level<1)then
  call oft_abort('Invalid mesh level','hgrad_ginterpmatrix',__FILE__)
end if
cmesh=>mg_mesh%meshes(mg_mesh%level-1)
if(cmesh%type/=1)CALL oft_abort("Only supported with tet meshes", &
  "hgrad_ginterpmatrix", __FILE__)
if(oft_hgrad%order/=2)then
  call oft_abort('Attempted geometric interpolation for pd /= 2','hgrad_ginterpmatrix',__FILE__)
end if
!---
ops=>oft_h1_ops
hgrad_cors=>ML_oft_hgrad%levels(oft_h1_level-1)%fe
lede=>mg_mesh%inter(mg_mesh%level-1)%lede
lfde=>mg_mesh%inter(mg_mesh%level-1)%lfde
lcdg=>mg_mesh%inter(mg_mesh%level-1)%lcdg
lcde=>mg_mesh%inter(mg_mesh%level-1)%lcde
ALLOCATE(ML_oft_hgrad%interp_graphs(ML_oft_hgrad%level)%g)
ops%hgrad_interp_graph=>ML_oft_hgrad%interp_graphs(ML_oft_hgrad%level)%g
!---Setup matrix sizes
ops%hgrad_interp_graph%nr=oft_hgrad%ne
ops%hgrad_interp_graph%nrg=oft_hgrad%global%ne
ops%hgrad_interp_graph%nc=hgrad_cors%ne
ops%hgrad_interp_graph%ncg=hgrad_cors%global%ne
!---Setup Matrix graph
ALLOCATE(ops%hgrad_interp_graph%kr(ops%hgrad_interp_graph%nr+1))
ops%hgrad_interp_graph%kr=0
ops%hgrad_interp_graph%nnz=cmesh%np+5*cmesh%ne+3*cmesh%nf+6*cmesh%nc
ALLOCATE(ops%hgrad_interp_graph%lc(ops%hgrad_interp_graph%nnz))
ops%hgrad_interp_graph%lc=0_i4
!---Construct linkage
!$omp parallel do
DO i=1,cmesh%np
  ops%hgrad_interp_graph%kr(i)=1
END DO
!$omp parallel do
DO i=1,cmesh%ne
  ops%hgrad_interp_graph%kr(cmesh%np+i)=3
  ops%hgrad_interp_graph%kr(mesh%np+lede(1,i))=1
  ops%hgrad_interp_graph%kr(mesh%np+lede(2,i))=1
END DO
!$omp parallel do
DO i=1,cmesh%nf
  ops%hgrad_interp_graph%kr(mesh%np+lfde(1,i))=1
  ops%hgrad_interp_graph%kr(mesh%np+lfde(2,i))=1
  ops%hgrad_interp_graph%kr(mesh%np+lfde(3,i))=1
END DO
!$omp parallel do
DO i=1,cmesh%nc
  ops%hgrad_interp_graph%kr(mesh%np+ABS(lcde(1,i)))=6
END DO
ops%hgrad_interp_graph%kr(ops%hgrad_interp_graph%nr+1)=ops%hgrad_interp_graph%nnz+1
do i=ops%hgrad_interp_graph%nr,1,-1 ! cumulative point to point count
  ops%hgrad_interp_graph%kr(i)=ops%hgrad_interp_graph%kr(i+1)-ops%hgrad_interp_graph%kr(i)
end do
if(ops%hgrad_interp_graph%kr(1)/=1)call oft_abort('Bad element to element count','hgrad_ginterpmatrix',__FILE__)
!$omp parallel do private(k)
DO i=1,cmesh%np
  k=ops%hgrad_interp_graph%kr(i)
  ops%hgrad_interp_graph%lc(k)=i
END DO
!$omp parallel do private(k)
DO i=1,cmesh%ne
  !---Daughter point
  k=ops%hgrad_interp_graph%kr(cmesh%np+i)
  ops%hgrad_interp_graph%lc(k)=cmesh%le(1,i)
  ops%hgrad_interp_graph%lc(k+1)=cmesh%le(2,i)
  ops%hgrad_interp_graph%lc(k+2)=cmesh%np+i
  !---Daughter edge 1
  k=ops%hgrad_interp_graph%kr(mesh%np+lede(1,i))
  ops%hgrad_interp_graph%lc(k)=cmesh%np+i
  !---Daughter edge 2
  k=ops%hgrad_interp_graph%kr(mesh%np+lede(2,i))
  ops%hgrad_interp_graph%lc(k)=cmesh%np+i
END DO
!$omp parallel do private(k,j,jfe)
DO i=1,cmesh%nf
  jfe=ABS(cmesh%lfe(:,i)) ! face edges
  !---
  DO j=1,3
    k=ops%hgrad_interp_graph%kr(mesh%np+lfde(j,i))
    ops%hgrad_interp_graph%lc(k)=cmesh%np+jfe(j)
  END DO
END DO
!$omp parallel do private(k,j,jce)
DO i=1,cmesh%nc
  jce=ABS(cmesh%lce(:,i)) ! face edges
  CALL sort_array(jce,6)
  !---
  k=ops%hgrad_interp_graph%kr(mesh%np+ABS(lcde(1,i)))
  DO j=0,5
    ops%hgrad_interp_graph%lc(k+j)=cmesh%np+jce(1+j)
  END DO
END DO
!---------------------------------------------------------------------------
! Construct matrix
!---------------------------------------------------------------------------
CALL oft_hgrad_create(hgrad_vec)
CALL oft_hgrad_create(hgrad_vec_cors,oft_h1_level-1)
!---
ALLOCATE(graphs(1,1))
graphs(1,1)%g=>ops%hgrad_interp_graph
!---
CALL create_matrix(mat,graphs,hgrad_vec,hgrad_vec_cors)
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
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
DEBUG_STACK_POP
END SUBROUTINE hgrad_ginterpmatrix
!---------------------------------------------------------------------------
! SUBROUTINE: h1_interp
!---------------------------------------------------------------------------
!> Interpolate a coarse level Lagrange scalar field to the next finest level
!!
!! @note The global Lagrange level in incremented by one in this subroutine
!!
!! @param[in] acors Vector to interpolate
!! @param[in,out] afine Fine vector from interpolation
!---------------------------------------------------------------------------
subroutine h1_interp(acors,afine)
class(oft_vector), intent(inout) :: acors
class(oft_vector), intent(inout) :: afine
integer(i4) :: i
real(r8), pointer, dimension(:) :: agrad,acurl,tmp
class(oft_mesh), pointer :: cmesh
DEBUG_STACK_PUSH
!---Step one level up
call oft_h1_set_level(oft_h1_level+1)
call afine%set(0.d0)
!---
if(oft_hcurl_level==oft_hcurl_blevel+1)then
  call h1_base_pop(acors,afine)
  DEBUG_STACK_POP
  return
end if
CALL oft_h1_ops%interp%apply(acors,afine)
!---Correct gradient subspace following geometric interpolation
IF(oft_hcurl%order==1)THEN
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
end subroutine h1_interp
!---------------------------------------------------------------------------
! SUBROUTINE: h1_base_pop
!---------------------------------------------------------------------------
!> Transfer a base level Lagrange scalar field to the next MPI level
!!
!! @param[in] acors Vector to transfer
!! @param[in,out] afine Fine vector from transfer
!---------------------------------------------------------------------------
subroutine h1_base_pop(acors,afine)
class(oft_vector), intent(inout) :: acors
class(oft_vector), intent(inout) :: afine
integer(i4), pointer, dimension(:) :: lptmp
integer(i4), pointer, dimension(:) :: lbege
integer(i4) :: i
real(r8), pointer, dimension(:) :: array_c,array_f
DEBUG_STACK_PUSH
!---
NULLIFY(array_c,array_f)
lbege=>mg_mesh%inter(mg_mesh%nbase)%lbege
CALL acors%get_local(array_c,1)
CALL afine%get_local(array_f,1)
!$omp parallel do
do i=1,afine%n
  array_f(i)=array_c(ABS(lbege(i)))
end do
CALL afine%restore_local(array_f,1)
DEALLOCATE(array_c,array_f)
!---
lptmp=>mg_mesh%meshes(mg_mesh%nbase+1)%base%lp
CALL acors%get_local(array_c,2)
CALL afine%get_local(array_f,2)
!$omp parallel do
do i=1,afine%n
  array_f(i)=array_c(lptmp(i))
end do
CALL afine%restore_local(array_f,2)
DEALLOCATE(array_c,array_f)
DEBUG_STACK_POP
end subroutine h1_base_pop
!---------------------------------------------------------------------------
! SUBROUTINE: h1_inject
!---------------------------------------------------------------------------
!> Inject a fine level Lagrange scalar field to the next coarsest level
!!
!! @note The global Lagrange level in decremented by one in this subroutine
!!
!! @param[in] afine Vector to inject
!! @param[in,out] acors Coarse vector from injection
!---------------------------------------------------------------------------
subroutine h1_inject(afine,acors)
class(oft_vector), intent(inout) :: afine
class(oft_vector), intent(inout) :: acors
real(r8), pointer, dimension(:) :: agrad,acurl
class(oft_vector), pointer :: tmp
integer(i4) :: i,j,k
logical :: gcheck
DEBUG_STACK_PUSH
gcheck=(oft_hcurl%order==1)
! Step down level up
call oft_h1_set_level(oft_h1_level-1)
! Cast fine field
call acors%set(0.d0)
if(oft_hcurl_level==oft_hcurl_blevel)then
  call h1_base_push(afine,acors)
  DEBUG_STACK_POP
  return
end if
CALL ML_oft_h1_ops(oft_h1_level+1)%interp%applyT(afine,acors)
DEBUG_STACK_POP
end subroutine h1_inject
!---------------------------------------------------------------------------
! SUBROUTINE: h1_base_push
!---------------------------------------------------------------------------
!> Transfer a MPI level Lagrange scalar field to the base level
!!
!! @param[in] afine Vector to transfer
!! @param[in,out] acors Fine vector from transfer
!---------------------------------------------------------------------------
subroutine h1_base_push(afine,acors)
class(oft_vector), intent(inout) :: afine
class(oft_vector), intent(inout) :: acors
integer(i4), pointer, dimension(:) :: lptmp
integer(i4), pointer, dimension(:) :: lbege
integer(i4) :: i,j,ierr
real(r8), pointer, dimension(:) :: alias,array_c,array_f
DEBUG_STACK_PUSH
!---
NULLIFY(array_c,array_f)
lbege=>mg_mesh%inter(mg_mesh%nbase)%lbege
lptmp=>mg_mesh%meshes(mg_mesh%nbase+1)%base%lp
CALL acors%get_local(array_c)
CALL afine%get_local(array_f)
!---
allocate(alias(acors%n))
alias=0.d0
!$omp parallel do
do i=1,oft_hcurl%ne
  if(oft_hcurl%linkage%be(i))cycle
  alias(ABS(lbege(i)))=array_f(i)
end do
!$omp parallel do private(j)
do i=1,oft_hcurl%linkage%nbe
  j=oft_hcurl%linkage%lbe(i)
  if(.NOT.oft_hcurl%linkage%leo(i))cycle
  alias(ABS(lbege(j)))=array_f(j)
end do
!$omp parallel do
do i=1,oft_hgrad%ne
  if(oft_hgrad%linkage%be(i))cycle
  alias(ML_oft_hcurl%levels(oft_h1_level-1)%fe%ne+lptmp(i))=array_f(oft_hcurl%ne+i)
end do
!$omp parallel do private(j)
do i=1,oft_hgrad%linkage%nbe
  j=oft_hgrad%linkage%lbe(i)
  if(.NOT.oft_hgrad%linkage%leo(i))cycle
  alias(ML_oft_hcurl%levels(oft_h1_level-1)%fe%ne+lptmp(j))=array_f(oft_hcurl%ne+j)
end do
!---Global reduction over all processors
array_c=oft_mpi_sum(alias,acors%n)
call acors%restore_local(array_c)
deallocate(alias,array_c,array_f)
DEBUG_STACK_POP
end subroutine h1_base_push
!---------------------------------------------------------------------------
! SUBROUTINE: h1_mop_eigs
!---------------------------------------------------------------------------
!> Compute eigenvalues and smoothing coefficients for the operator H1::MOP
!---------------------------------------------------------------------------
SUBROUTINE h1_mop_eigs(minlev)
INTEGER(i4), INTENT(in) :: minlev
#ifdef HAVE_ARPACK
INTEGER(i4) :: i
REAL(r8) :: lam0
REAL(r8), ALLOCATABLE :: df(:)
CLASS(oft_vector), POINTER :: u
TYPE(oft_irlm_eigsolver) :: arsolver
CLASS(oft_matrix), POINTER :: md => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
DEBUG_STACK_PUSH
!---------------------------------------------------------------------------
! Compute optimal smoother coefficients
!---------------------------------------------------------------------------
IF(oft_env%head_proc)WRITE(*,*)'Optimizing Jacobi damping for H1::MOP'
ALLOCATE(df(oft_h1_nlevels))
df=0.d0
DO i=minlev,oft_h1_nlevels
  CALL oft_h1_set_level(i)
  !---Create fields
  CALL oft_h1_create(u)
  !---Get Ev range
  NULLIFY(mop)
  CALL h1_getmop(mop,'lop')
  CALL create_diagmatrix(md,mop%D)
  !---
  arsolver%A=>mop
  arsolver%M=>md
  arsolver%mode=2
  arsolver%tol=1.E-5_r8
  arsolver%bc=>h1_zerob
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
  DO i=1,oft_h1_nlevels-1
    WRITE(*,'(F5.3,A)',ADVANCE='NO')df(i),', '
  END DO
  WRITE(*,'(F5.3,A)')df(oft_h1_nlevels)
END IF
DEALLOCATE(df)
DEBUG_STACK_POP
#else
CALL oft_abort("Subroutine requires ARPACK", "lag_lop_eigs", __FILE__)
#endif
END SUBROUTINE h1_mop_eigs
!---------------------------------------------------------------------------
! SUBROUTINE: h1_getmop_pre
!---------------------------------------------------------------------------
!> Compute eigenvalues and smoothing coefficients for the operator H1::MOP
!---------------------------------------------------------------------------
SUBROUTINE h1_getmop_pre(pre,mats,level,nlevels)
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
INTEGER(i4) :: i,j,levin
LOGICAL :: create_mats
CHARACTER(LEN=2) :: lev_char
TYPE(fox_node), POINTER :: pre_node
#ifdef HAVE_XML
integer(i4) :: nnodes
TYPE(fox_node), POINTER :: h1_node
TYPE(fox_nodelist), POINTER :: current_nodes
#endif
DEBUG_STACK_PUSH
!---
minlev=1
toplev=oft_h1_level
levin=oft_h1_level
IF(PRESENT(level))toplev=level
IF(PRESENT(nlevels))minlev=toplev-nlevels+1
nl=toplev-minlev+1
!---
IF(minlev<1)CALL oft_abort('Minimum level is < 0','h1_getmop_pre',__FILE__)
IF(toplev>oft_h1_nlevels)CALL oft_abort('Maximum level is > h1_nlevels','h1_getmop_pre',__FILE__)
!---------------------------------------------------------------------------
! Create ML Matrices
!---------------------------------------------------------------------------
create_mats=.FALSE.
IF(.NOT.ASSOCIATED(mats))THEN
  create_mats=.TRUE.
  ALLOCATE(mats(nl))
END IF
ALLOCATE(ml_int(nl-1),levels(nl))
ALLOCATE(df(nl),nu(nl))
DO i=1,nl
  CALL oft_h1_set_level(minlev+(i-1))
  levels(i)=minlev+(i-1)
  df(i)=df_mop(levels(i))
  nu(i)=nu_mop(levels(i))
  IF(df(i)<-1.d90)THEN
    WRITE(lev_char,'(I2.2)')levels(i)
    CALL oft_abort('Smoother values not set for level: '//lev_char,'h1_getmop_pre',__FILE__)
  END IF
  !---
  IF(create_mats)THEN
    NULLIFY(mats(i)%M)
    CALL h1_getmop(mats(i)%M,'none')
  END IF
  IF(i>1)ml_int(i-1)%M=>oft_h1_ops%interp
END DO
CALL oft_h1_set_level(levin)
!---------------------------------------------------------------------------
! Search for XML-spec
!---------------------------------------------------------------------------
NULLIFY(pre_node)
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml))THEN
  !---Look for Lagrange node
  current_nodes=>fox_getElementsByTagName(oft_env%xml,"nedelec_h1")
  nnodes=fox_getLength(current_nodes)
  IF(nnodes>0)THEN
    h1_node=>fox_item(current_nodes,0)
    !---Look for lop node
    current_nodes=>fox_getElementsByTagName(h1_node,"mop")
    nnodes=fox_getLength(current_nodes)
    IF(nnodes>0)pre_node=>fox_item(current_nodes,0)
  END IF
END IF
#endif
!---------------------------------------------------------------------------
! Setup preconditioner
!---------------------------------------------------------------------------
NULLIFY(pre)
CALL create_mlpre(pre,mats(1:nl),levels,nlevels=nl,create_vec=oft_h1_create,interp=h1_interp, &
     inject=h1_inject,stype=1,df=df,nu=nu,xml_root=pre_node)
!---------------------------------------------------------------------------
! Cleanup
!---------------------------------------------------------------------------
DEALLOCATE(ml_int,levels,df,nu)
DEBUG_STACK_POP
END SUBROUTINE h1_getmop_pre
!---------------------------------------------------------------------------
! FUNCTION: h1_jump_error
!---------------------------------------------------------------------------
!> Evaluate the jump error in a field over internal faces
!!
!! @note Currently faces on domain boundaries are skipped, this is due to the
!! fact that evaluting the error would require costly communication.
!!
!! @param[in,out] u H1 vector field to evaluate
!! @param[in] quad_order Desired quadrature order for integration
!! @return Jump error metric
!---------------------------------------------------------------------------
FUNCTION h1_jump_error(u,quad_order) RESULT(error)
CLASS(oft_vector), INTENT(inout) :: u
INTEGER(i4), INTENT(in) :: quad_order
REAL(r8) :: error,reg_jump,reg_energy,goptmp(3,3)
REAL(r8) :: gop(3,4),vol,area,val(3),f(4),norm(3)
REAL(r8), ALLOCATABLE :: Bn(:,:),rop_curl(:,:),rop_grad(:,:)
REAL(r8), POINTER :: curl_vals(:),grad_vals(:)
INTEGER(i4) :: i,j,m,cf,jc,ptmap(3)
INTEGER(i4) :: cell(2),face(2)
INTEGER(i4), ALLOCATABLE :: j_hcurl(:),j_hgrad(:)
TYPE(oft_quad_type) :: quad
CLASS(oft_bmesh), POINTER :: mesh_tmp
DEBUG_STACK_PUSH
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
ALLOCATE(j_hcurl(oft_hcurl%nce),j_hgrad(oft_hgrad%nce))
ALLOCATE(rop_curl(3,oft_hcurl%nce),rop_grad(3,oft_hgrad%nce))
!---
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
    call oft_hcurl%ncdofs(cell(j),j_hcurl) ! get curl DOFs
    call oft_hgrad%ncdofs(cell(j),j_hgrad) ! get grad DOFs
    !---Get local reconstructed operators
    DO m=1,quad%np ! Loop over quadrature points
      !---
      f=0.d0; f(mesh%cell_fc(:,cf))=quad%pts(ptmap,m)
      CALL mesh%jacobian(cell(j),f,gop,vol)
      !---
      val=0.d0
      CALL oft_hcurl_eval_all(oft_hcurl,cell(j),f,rop_curl,gop)
      CALL oft_h0_geval_all(oft_hgrad,cell(j),f,rop_grad,gop)
      DO jc=1,oft_hcurl%nce
        val=val+curl_vals(j_hcurl(jc))*rop_curl(:,jc)
      END DO
      DO jc=1,oft_hgrad%nce
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
END FUNCTION h1_jump_error
end module oft_h1_operators
