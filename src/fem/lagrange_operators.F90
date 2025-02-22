!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_lag_operators.F90
!
!> Lagrange FE operator definitions
!! - Operator construction
!!   - MOP: mass matrix        \f$ \int \left( u^T v \right) dV \f$
!!   - LOP: laplacian matrix   \f$ \int \left( \nabla u^T \cdot \nabla v \right) dV \f$
!!   - PDOP: parallel diffusion matrix   \f$ \int \left( \nabla u^T \cdot \hat{b} \hat{b} \cdot \nabla v \right) dV \f$
!!   - VMOP: vector mass matrix    \f$ \int \left( u^T \cdot v \right) dV \f$
!! - Interpolation classes
!! - Field projection (to plotting mesh)
!! - Boundary conditions
!! - Multi-Grid setup and operators
!! - Default preconditioner setup
!!   - Optimal smoother calculation
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_lag
!---------------------------------------------------------------------------
MODULE oft_lag_operators
USE oft_base
USE oft_sort, ONLY: sort_array
USE oft_mesh_type, ONLY: oft_mesh, mesh, cell_is_curved
USE multigrid, ONLY: mg_mesh
!---
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr, &
  oft_graph, oft_graph_ptr
USE oft_deriv_matrices, ONLY: oft_diagmatrix, create_diagmatrix
USE oft_solver_base, ONLY: oft_solver
USE oft_la_utils, ONLY: create_vector, create_matrix, combine_matrices
USE oft_solver_utils, ONLY: create_mlpre, create_native_pre
#ifdef HAVE_ARPACK
USE oft_arpack, ONLY: oft_irlm_eigsolver
#endif
!---
USE fem_base, ONLY: oft_afem_type, oft_fem_type, oft_bfem_type, fem_max_levels
USE fem_utils, ONLY: fem_interp
USE oft_lag_basis, ONLY: oft_lagrange, oft_lagrange_level, oft_lagrange_nlevels, oft_lag_set_level, &
oft_lagrange_blevel, ML_oft_lagrange, oft_lagrange_ops, oft_lag_ops, ML_oft_lagrange_ops, &
oft_lag_eval_all, oft_lag_geval_all, oft_lag_eval, &
oft_lag_nodes, oft_lagrange_lev, ML_oft_vlagrange, oft_vlagrange, oft_blagrange, oft_blag_eval, &
oft_blag_geval, oft_lagrange_minlev, oft_lag_npos, oft_scalar_fem, oft_scalar_bfem
USE oft_lag_fields, ONLY: oft_lag_create, oft_blag_create, oft_lag_vcreate
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------
! CLASS oft_lag_rinterp
!---------------------------------------------------------------------------
!> Interpolate a Lagrange field
!---------------------------------------------------------------------------
type, extends(fem_interp) :: oft_lag_rinterp
  class(oft_vector), pointer :: u => NULL() !< Field for interpolation
  integer(i4), pointer, dimension(:) :: cache_cell => NULL()
  real(r8), pointer, dimension(:) :: vals => NULL() !< Local values
  real(r8), pointer, dimension(:,:) :: cache_vals => NULL()
  class(oft_scalar_fem), pointer :: lag_rep => NULL() !< Lagrange FE representation
contains
  !> Retrieve local values for interpolation
  procedure :: setup => lag_rinterp_setup
  !> Reconstruct field
  procedure :: interp => lag_rinterp
  !> Delete reconstruction object
  procedure :: delete => lag_rinterp_delete
end type oft_lag_rinterp
!---------------------------------------------------------------------------
! CLASS oft_lag_ginterp
!---------------------------------------------------------------------------
!> Interpolate \f$ \nabla \f$ of a Lagrange field
!---------------------------------------------------------------------------
type, extends(oft_lag_rinterp) :: oft_lag_ginterp
contains
  !> Reconstruct field
  procedure :: interp => lag_ginterp_apply
end type oft_lag_ginterp
! !---------------------------------------------------------------------------
! ! CLASS oft_lag_brinterp
! !---------------------------------------------------------------------------
! !> Interpolate a boundary Lagrange field
! !---------------------------------------------------------------------------
! type, extends(fem_interp) :: oft_lag_brinterp
!   class(oft_vector), pointer :: u => NULL() !< Field for interpolation
!   real(r8), pointer, dimension(:) :: vals => NULL() !< Local values
!   class(oft_scalar_bfem), pointer :: lag_rep => NULL() !< Lagrange FE representation
! contains
!   !> Retrieve local values for interpolation
!   procedure :: setup => lag_brinterp_setup
!   !> Reconstruct field
!   procedure :: interp => lag_brinterp
!   !> Delete reconstruction object
!   procedure :: delete => lag_brinterp_delete
! end type oft_lag_brinterp
!---------------------------------------------------------------------------
! CLASS oft_lag_vrinterp
!---------------------------------------------------------------------------
!> Interpolate a Lagrange vector field
!---------------------------------------------------------------------------
type, extends(fem_interp) :: oft_lag_vrinterp
  class(oft_vector), pointer :: u => NULL() !< Field for interpolation
  integer(i4), pointer, dimension(:) :: cache_cell => NULL()
  real(r8), pointer, dimension(:,:) :: vals => NULL() !< Local values
  real(r8), pointer, dimension(:,:,:) :: cache_vals => NULL()
  class(oft_scalar_fem), pointer :: lag_rep => NULL() !< Lagrange FE representation
contains
  !> Retrieve local values for interpolation
  procedure :: setup => lag_vrinterp_setup
  !> Reconstruct field
  procedure :: interp => lag_vrinterp
  !> Delete reconstruction object
  procedure :: delete => lag_vrinterp_delete
end type oft_lag_vrinterp
!---------------------------------------------------------------------------
! CLASS oft_lag_vcinterp
!---------------------------------------------------------------------------
!> Interpolate \f$ \nabla \times \f$ of a Lagrange vector field
!---------------------------------------------------------------------------
type, extends(oft_lag_vrinterp) :: oft_lag_vcinterp
contains
  !> Reconstruct field
  procedure :: interp => lag_vcinterp
end type oft_lag_vcinterp
!---------------------------------------------------------------------------
! CLASS oft_lag_vdinterp
!---------------------------------------------------------------------------
!> Interpolate \f$ \nabla \f$ of a Lagrange vector field
!---------------------------------------------------------------------------
type, extends(oft_lag_vrinterp) :: oft_lag_vdinterp
contains
  !> Reconstruct field
  procedure :: interp => lag_vdinterp
end type oft_lag_vdinterp
! !---------------------------------------------------------------------------
! ! CLASS oft_lag_bvrinterp
! !---------------------------------------------------------------------------
! !> Interpolate a boundary Lagrange vector field
! !---------------------------------------------------------------------------
! type, extends(fem_interp) :: oft_lag_bvrinterp
!   class(oft_vector), pointer :: u => NULL() !< Field for interpolation
!   real(r8), pointer, dimension(:,:) :: vals => NULL() !< Local values
!   class(oft_scalar_bfem), pointer :: lag_rep => NULL() !< Lagrange FE representation
! contains
!   !> Retrieve local values for interpolation
!   procedure :: setup => lag_bvrinterp_setup
!   !> Reconstruct field
!   procedure :: interp => lag_bvrinterp
!   !> Delete reconstruction object
!   procedure :: delete => lag_bvrinterp_delete
! end type oft_lag_bvrinterp
!---Pre options
integer(i4) :: nu_lop(fem_max_levels)=0
real(r8) :: df_lop(fem_max_levels)=-1.d99
integer(i4), private :: nu_pdop(fem_max_levels)=0
real(r8), private :: df_pdop(fem_max_levels)=-1.d99
!---Cache variables
REAL(r8), POINTER, DIMENSION(:) :: oft_lag_rop => NULL()
REAL(r8), POINTER, DIMENSION(:,:) :: oft_lag_gop => NULL()
!$omp threadprivate(oft_lag_rop,oft_lag_gop)
contains
!---------------------------------------------------------------------------
! SUBROUTINE: lag_mloptions
!---------------------------------------------------------------------------
!> Read-in options for the basic Lagrange ML preconditioners
!---------------------------------------------------------------------------
subroutine lag_mloptions()
integer(i4) :: ierr,io_unit
namelist/lag_op_options/df_lop,nu_lop,df_pdop,nu_pdop
DEBUG_STACK_PUSH
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,lag_op_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr>0)THEN
  CALL oft_abort('Error reading "lag_op_options" group in input file','lag_mloptions',__FILE__)
END IF
IF(df_lop(1)<-1.d90)THEN
  IF(oft_env%head_proc)THEN
    WRITE(*,*)'No Lagrange MG smoother settings found:'
    WRITE(*,*)'  Using default values, which may result in convergence failure.'
  END IF
  nu_lop=2
  df_lop=.2d0
  nu_pdop=2
  df_pdop=.2d0
END IF
DEBUG_STACK_POP
end subroutine lag_mloptions
!---------------------------------------------------------------------------
! SUBROUTINE: lag_rinterp_setup
!---------------------------------------------------------------------------
!> Setup interpolator for Lagrange scalar fields
!!
!! Fetches local representation used for interpolation from vector object
!!
!! @note Should only be used via class \ref oft_lag_rinterp or children
!---------------------------------------------------------------------------
subroutine lag_rinterp_setup(self)
class(oft_lag_rinterp), intent(inout) :: self
IF(ASSOCIATED(self%parent))THEN
  SELECT TYPE(this=>self%parent)
    CLASS IS(oft_lag_rinterp)
      self%vals=>this%vals
      self%lag_rep=>this%lag_rep
      self%cache_cell=>this%cache_cell
      self%cache_vals=>this%cache_vals
    CLASS DEFAULT
      CALL oft_abort('Parent interpolator must be of same class','oft_lag_rinterp',__FILE__)
  END SELECT
ELSE
  !---Get local slice
  CALL self%u%get_local(self%vals)
  self%lag_rep=>oft_lagrange
  IF(.NOT.ASSOCIATED(self%cache_cell))THEN
    ALLOCATE(self%cache_cell(0:oft_env%nthreads-1))
    ALLOCATE(self%cache_vals(self%lag_rep%nce,0:oft_env%nthreads-1))
  END IF
  self%cache_cell=-1
  self%cache_vals=0.d0
END IF
end subroutine lag_rinterp_setup
!---------------------------------------------------------------------------
! SUBROUTINE: lag_rinterp_delete
!---------------------------------------------------------------------------
!> Destroy temporary internal storage
!!
!! @note Should only be used via class \ref oft_lag_rinterp or children
!---------------------------------------------------------------------------
subroutine lag_rinterp_delete(self)
class(oft_lag_rinterp), intent(inout) :: self
IF(ASSOCIATED(self%parent))THEN
  NULLIFY(self%vals)
  NULLIFY(self%cache_cell,self%cache_vals)
  NULLIFY(self%lag_rep,self%u)
  NULLIFY(self%parent)
ELSE
  !---Destroy locals
  IF(ASSOCIATED(self%vals))DEALLOCATE(self%vals,self%cache_cell,self%cache_vals)
  NULLIFY(self%lag_rep,self%u)
END IF
end subroutine lag_rinterp_delete
!---------------------------------------------------------------------------
! SUBROUTINE: lag_rinterp
!---------------------------------------------------------------------------
!> Reconstruct a Lagrange scalar field
!!
!! @note Should only be used via class \ref oft_lag_rinterp
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [1]
!---------------------------------------------------------------------------
subroutine lag_rinterp(self,cell,f,gop,val)
class(oft_lag_rinterp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: rop(1)
real(r8), allocatable :: lag_rop(:)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!','lag_rinterp',__FILE__)
!---Pull cache
IF(self%cache_cell(oft_tid)/=cell)THEN
  allocate(j(self%lag_rep%nce))
  call self%lag_rep%ncdofs(cell,j)
  do jc=1,self%lag_rep%nce
    self%cache_vals(jc,oft_tid)=self%vals(j(jc))
  end do
  deallocate(j)
  self%cache_cell(oft_tid)=cell
END IF
!---Reconstruct field
val=0.d0
IF(ASSOCIATED(oft_lag_rop))THEN
  do jc=1,self%lag_rep%nce
    val=val+self%cache_vals(jc,oft_tid)*oft_lag_rop(jc)
  end do
ELSE
  ALLOCATE(lag_rop(self%lag_rep%nce))
  CALL oft_lag_eval_all(self%lag_rep,cell,f,lag_rop)
  do jc=1,self%lag_rep%nce
    val=val+self%cache_vals(jc,oft_tid)*lag_rop(jc)
  end do
  DEALLOCATE(lag_rop)
END IF
DEBUG_STACK_POP
end subroutine lag_rinterp
! !---------------------------------------------------------------------------
! ! SUBROUTINE: lag_brinterp_setup
! !---------------------------------------------------------------------------
! !> Setup interpolator for boundary Lagrange scalar fields
! !!
! !! Fetches local representation used for interpolation from vector object
! !!
! !! @note Should only be used via class \ref oft_lag_brinterp or children
! !---------------------------------------------------------------------------
! subroutine lag_brinterp_setup(self)
! class(oft_lag_brinterp), intent(inout) :: self
! !---Get local slice
! CALL self%u%get_local(self%vals)
! IF(.NOT.ASSOCIATED(self%lag_rep))self%lag_rep=>oft_blagrange
! end subroutine lag_brinterp_setup
! !---------------------------------------------------------------------------
! ! SUBROUTINE: lag_brinterp_delete
! !---------------------------------------------------------------------------
! !> Destroy temporary internal storage
! !!
! !! @note Should only be used via class \ref oft_lag_brinterp or children
! !---------------------------------------------------------------------------
! subroutine lag_brinterp_delete(self)
! class(oft_lag_brinterp), intent(inout) :: self
! !---Destroy locals
! IF(ASSOCIATED(self%vals))DEALLOCATE(self%vals)
! NULLIFY(self%lag_rep,self%u)
! end subroutine lag_brinterp_delete
! !---------------------------------------------------------------------------
! ! SUBROUTINE: lag_brinterp
! !---------------------------------------------------------------------------
! !> Reconstruct a boundary Lagrange scalar field
! !!
! !! @note Should only be used via class \ref oft_lag_brinterp
! !!
! !! @param[in] cell Cell for interpolation
! !! @param[in] f Possition in cell in logical coord [4]
! !! @param[in] gop Logical gradient vectors at f [3,4]
! !! @param[out] val Reconstructed field at f [1]
! !---------------------------------------------------------------------------
! subroutine lag_brinterp(self,cell,f,gop,val)
! class(oft_lag_brinterp), intent(inout) :: self
! integer(i4), intent(in) :: cell
! real(r8), intent(in) :: f(:)
! real(r8), intent(in) :: gop(3,4)
! real(r8), intent(out) :: val(:)
! integer(i4), allocatable :: j(:)
! integer(i4) :: jc
! real(r8) :: rop(1)
! DEBUG_STACK_PUSH
! !---
! IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!','lag_brinterp',__FILE__)
! !---Get dofs
! allocate(j(self%lag_rep%nfe))
! call self%lag_rep%nfdofs(cell,j) ! get DOFs
! !---Reconstruct field
! val=0.d0
! do jc=1,self%lag_rep%nfe
!   call oft_blag_eval(self%lag_rep,cell,jc,f,rop(1))
!   val=val+self%vals(j(jc))*rop
! end do
! deallocate(j)
! DEBUG_STACK_POP
! end subroutine lag_brinterp
!---------------------------------------------------------------------------
! SUBROUTINE: lag_ginterp_apply
!---------------------------------------------------------------------------
!> Reconstruct the gradient of a Lagrange scalar field
!!
!! @note Should only be used via class \ref oft_lag_ginterp
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed gradient at f [3]
!---------------------------------------------------------------------------
subroutine lag_ginterp_apply(self,cell,f,gop,val)
class(oft_lag_ginterp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: rop(3)
real(r8), allocatable :: lag_rop(:,:)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!','lag_ginterp_apply',__FILE__)
!---Pull cache
IF(self%cache_cell(oft_tid)/=cell)THEN
  allocate(j(self%lag_rep%nce))
  call self%lag_rep%ncdofs(cell,j)
  do jc=1,self%lag_rep%nce
    self%cache_vals(jc,oft_tid)=self%vals(j(jc))
  end do
  deallocate(j)
  self%cache_cell(oft_tid)=cell
END IF
!---Reconstruct field
val=0.d0
IF(ASSOCIATED(oft_lag_gop))THEN
  do jc=1,self%lag_rep%nce
    val=val+self%cache_vals(jc,oft_tid)*oft_lag_gop(:,jc)
  end do
ELSE
  ALLOCATE(lag_rop(3,self%lag_rep%nce))
  CALL oft_lag_geval_all(self%lag_rep,cell,f,lag_rop,gop)
  do jc=1,self%lag_rep%nce
    val=val+self%cache_vals(jc,oft_tid)*lag_rop(:,jc)
  end do
  DEALLOCATE(lag_rop)
END IF
DEBUG_STACK_POP
end subroutine lag_ginterp_apply
!---------------------------------------------------------------------------
! SUBROUTINE: lag_vrinterp_setup
!---------------------------------------------------------------------------
!> Setup interpolator for Lagrange vector fields
!!
!! Fetches local representation used for interpolation from vector object
!!
!! @note Should only be used via class \ref oft_lag_vrinterp or children
!---------------------------------------------------------------------------
subroutine lag_vrinterp_setup(self)
class(oft_lag_vrinterp), intent(inout) :: self
real(r8), pointer, dimension(:) :: vtmp
IF(ASSOCIATED(self%parent))THEN
  SELECT TYPE(this=>self%parent)
    CLASS IS(oft_lag_vrinterp)
      self%vals=>this%vals
      self%lag_rep=>this%lag_rep
      self%cache_cell=>this%cache_cell
      self%cache_vals=>this%cache_vals
    CLASS DEFAULT
      CALL oft_abort('Parent interpolator must be of same class','lag_vrinterp_setup',__FILE__)
  END SELECT
ELSE
  !---Get local slice
  IF(.NOT.ASSOCIATED(self%vals))ALLOCATE(self%vals(3,oft_lagrange%ne))
  vtmp=>self%vals(1,:)
  CALL self%u%get_local(vtmp,1)
  vtmp=>self%vals(2,:)
  CALL self%u%get_local(vtmp,2)
  vtmp=>self%vals(3,:)
  CALL self%u%get_local(vtmp,3)
  self%lag_rep=>oft_lagrange
  IF(.NOT.ASSOCIATED(self%cache_cell))THEN
    ALLOCATE(self%cache_cell(0:oft_env%nthreads-1))
    ALLOCATE(self%cache_vals(3,self%lag_rep%nce,0:oft_env%nthreads-1))
  END IF
  self%cache_cell=-1
  self%cache_vals=0.d0
END IF
end subroutine lag_vrinterp_setup
!---------------------------------------------------------------------------
! SUBROUTINE: lag_vrinterp_delete
!---------------------------------------------------------------------------
!> Setup interpolator for Lagrange vector fields
!!
!! Fetches local representation used for interpolation from vector object
!!
!! @note Should only be used via class \ref oft_lag_vrinterp or children
!---------------------------------------------------------------------------
subroutine lag_vrinterp_delete(self)
class(oft_lag_vrinterp), intent(inout) :: self
IF(ASSOCIATED(self%parent))THEN
  NULLIFY(self%vals)
  NULLIFY(self%cache_cell,self%cache_vals)
  NULLIFY(self%lag_rep,self%u)
  NULLIFY(self%parent)
ELSE
  !---Destroy locals
  IF(ASSOCIATED(self%vals))DEALLOCATE(self%vals,self%cache_cell,self%cache_vals)
  NULLIFY(self%lag_rep,self%u)
END IF
end subroutine lag_vrinterp_delete
!---------------------------------------------------------------------------
! SUBROUTINE: lag_vrinterp
!---------------------------------------------------------------------------
!> Reconstruct a Lagrange vector field
!!
!! @note Should only be used via class \ref oft_lag_vrinterp
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [3]
!---------------------------------------------------------------------------
subroutine lag_vrinterp(self,cell,f,gop,val)
class(oft_lag_vrinterp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: rop(1)
real(r8), allocatable  :: lag_rop(:)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!','lag_vrinterp',__FILE__)
!---Pull cache
IF(self%cache_cell(oft_tid)/=cell)THEN
  allocate(j(self%lag_rep%nce))
  call self%lag_rep%ncdofs(cell,j)
  do jc=1,self%lag_rep%nce
    self%cache_vals(:,jc,oft_tid)=self%vals(:,j(jc))
  end do
  deallocate(j)
  self%cache_cell(oft_tid)=cell
END IF
!---Reconstruct field
val=0.d0
IF(ASSOCIATED(oft_lag_rop))THEN
  do jc=1,self%lag_rep%nce
    val=val+self%cache_vals(:,jc,oft_tid)*oft_lag_rop(jc)
  end do
ELSE
  ALLOCATE(lag_rop(self%lag_rep%nce))
  CALL oft_lag_eval_all(self%lag_rep,cell,f,lag_rop)
  do jc=1,self%lag_rep%nce
    val=val+self%cache_vals(:,jc,oft_tid)*lag_rop(jc)
  end do
  DEALLOCATE(lag_rop)
END IF
DEBUG_STACK_POP
end subroutine lag_vrinterp
!---------------------------------------------------------------------------
! SUBROUTINE: lag_vcinterp
!---------------------------------------------------------------------------
!> Reconstruct the curl of a Lagrange vector field
!!
!! @note Should only be used via class \ref oft_lag_vcinterp
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed curl at f [3]
!---------------------------------------------------------------------------
subroutine lag_vcinterp(self,cell,f,gop,val)
class(oft_lag_vcinterp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: x,y,z
real(r8), allocatable  :: rop(:,:)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!','lag_vcinterp',__FILE__)
!---Get dofs
allocate(j(self%lag_rep%nce),rop(3,self%lag_rep%nce))
call self%lag_rep%ncdofs(cell,j) ! get DOFs
!---Reconstruct field
val=0.d0
call oft_lag_geval_all(self%lag_rep,cell,f,rop,gop)
do jc=1,self%lag_rep%nce
  x=self%vals(1,j(jc))
  y=self%vals(2,j(jc))
  z=self%vals(3,j(jc))
  val(1)=val(1) + z*rop(2,jc)-y*rop(3,jc)
  val(2)=val(2) + x*rop(3,jc)-z*rop(1,jc)
  val(3)=val(3) + y*rop(1,jc)-x*rop(2,jc)
end do
deallocate(j,rop)
DEBUG_STACK_POP
end subroutine lag_vcinterp
!---------------------------------------------------------------------------
! SUBROUTINE: lag_vdinterp
!---------------------------------------------------------------------------
!> Reconstruct \f$ \nabla v \f$ for a Lagrange vector field.
!!
!! The tensor is packed using reshape, to retrieve use
!!\code
!! dv = RESHAPE(val,(/3,3/))
!!\endcode
!!
!! @note Should only be used via class \ref oft_lag_vdinterp
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed tensor at f [9]
!---------------------------------------------------------------------------
subroutine lag_vdinterp(self,cell,f,gop,val)
class(oft_lag_vdinterp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: rop(3),dv(3,3),vtmp(9,1)
real(r8), allocatable :: lag_rop(:,:)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!','lag_vdinterp',__FILE__)
!---Pull cache
IF(self%cache_cell(oft_tid)/=cell)THEN
  allocate(j(self%lag_rep%nce))
  call self%lag_rep%ncdofs(cell,j)
  do jc=1,self%lag_rep%nce
    self%cache_vals(:,jc,oft_tid)=self%vals(:,j(jc))
  end do
  deallocate(j)
  self%cache_cell(oft_tid)=cell
END IF
!---Reconstruct field
dv=0.d0
IF(ASSOCIATED(oft_lag_gop))THEN
  do jc=1,self%lag_rep%nce
    dv(:,1)=dv(:,1)+oft_lag_gop(:,jc)*self%cache_vals(1,jc,oft_tid)
    dv(:,2)=dv(:,2)+oft_lag_gop(:,jc)*self%cache_vals(2,jc,oft_tid)
    dv(:,3)=dv(:,3)+oft_lag_gop(:,jc)*self%cache_vals(3,jc,oft_tid)
  end do
ELSE
  ALLOCATE(lag_rop(3,self%lag_rep%nce))
  CALL oft_lag_geval_all(self%lag_rep,cell,f,lag_rop,gop)
  do jc=1,self%lag_rep%nce
    dv(:,1)=dv(:,1)+self%cache_vals(1,jc,oft_tid)*lag_rop(:,jc)
    dv(:,2)=dv(:,2)+self%cache_vals(2,jc,oft_tid)*lag_rop(:,jc)
    dv(:,3)=dv(:,3)+self%cache_vals(3,jc,oft_tid)*lag_rop(:,jc)
  end do
  DEALLOCATE(lag_rop)
END IF
vtmp=RESHAPE(dv,(/9,1/))
val=vtmp(:,1)
DEBUG_STACK_POP
end subroutine lag_vdinterp
! !---------------------------------------------------------------------------
! ! SUBROUTINE: lag_bvrinterp_setup
! !---------------------------------------------------------------------------
! !> Setup interpolator for boundary Lagrange vector fields
! !!
! !! Fetches local representation used for interpolation from vector object
! !!
! !! @note Should only be used via class \ref oft_lag_bvrinterp or children
! !---------------------------------------------------------------------------
! subroutine lag_bvrinterp_setup(self)
! class(oft_lag_bvrinterp), intent(inout) :: self
! real(r8), pointer, dimension(:) :: vtmp
! !---Get local slice
! IF(.NOT.ASSOCIATED(self%vals))ALLOCATE(self%vals(3,oft_lagrange%ne))
! vtmp=>self%vals(1,:)
! CALL self%u%get_local(vtmp,1)
! vtmp=>self%vals(2,:)
! CALL self%u%get_local(vtmp,2)
! vtmp=>self%vals(3,:)
! CALL self%u%get_local(vtmp,3)
! IF(.NOT.ASSOCIATED(self%lag_rep))self%lag_rep=>oft_blagrange
! end subroutine lag_bvrinterp_setup
! !---------------------------------------------------------------------------
! ! SUBROUTINE: lag_bvrinterp_delete
! !---------------------------------------------------------------------------
! !> Destroy temporary internal storage
! !!
! !! @note Should only be used via class \ref oft_lag_bvrinterp or children
! !---------------------------------------------------------------------------
! subroutine lag_bvrinterp_delete(self)
! class(oft_lag_bvrinterp), intent(inout) :: self
! !---Destroy locals
! IF(ASSOCIATED(self%vals))DEALLOCATE(self%vals)
! NULLIFY(self%lag_rep,self%u)
! end subroutine lag_bvrinterp_delete
! !---------------------------------------------------------------------------
! ! SUBROUTINE: lag_bvrinterp
! !---------------------------------------------------------------------------
! !> Reconstruct a boundary Lagrange vector field
! !!
! !! @note Should only be used via class \ref oft_lag_bvrinterp
! !!
! !! @param[in] cell Cell for interpolation
! !! @param[in] f Possition in cell in logical coord [4]
! !! @param[in] gop Logical gradient vectors at f [3,4]
! !! @param[out] val Reconstructed field at f [3]
! !---------------------------------------------------------------------------
! subroutine lag_bvrinterp(self,cell,f,gop,val)
! class(oft_lag_bvrinterp), intent(inout) :: self
! integer(i4), intent(in) :: cell
! real(r8), intent(in) :: f(:)
! real(r8), intent(in) :: gop(3,4)
! real(r8), intent(out) :: val(:)
! integer(i4), allocatable :: j(:)
! integer(i4) :: jc
! real(r8) :: rop(1)
! DEBUG_STACK_PUSH
! !---
! IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!', &
! 'lag_bvrinterp',__FILE__)
! !---Get dofs
! allocate(j(self%lag_rep%nfe))
! call self%lag_rep%nfdofs(cell,j) ! get DOFs
! !---Reconstruct field
! val=0.d0
! do jc=1,self%lag_rep%nfe
!   call oft_blag_eval(self%lag_rep,cell,jc,f,rop(1))
!   val=val+self%vals(:,j(jc))*rop(1)
! end do
! deallocate(j)
! DEBUG_STACK_POP
! end subroutine lag_bvrinterp
!---------------------------------------------------------------------------
! SUBROUTINE: lag_zerob
!---------------------------------------------------------------------------
!> Zero a Lagrange scalar field at all boundary nodes
!!
!! @param[in,out] a Field to be zeroed
!---------------------------------------------------------------------------
subroutine lag_zerob(a)
class(oft_vector), intent(inout) :: a
integer(i4) :: i,j
real(r8), pointer, dimension(:) :: vloc
DEBUG_STACK_PUSH
!---Cast to vector type
NULLIFY(vloc)
call a%get_local(vloc)
!---Zero boundary values
!$omp parallel do private(j)
do i=1,oft_lagrange%nbe
  j=oft_lagrange%lbe(i)
  if(oft_lagrange%global%gbe(j))vloc(j)=0.d0
end do
!---
call a%restore_local(vloc)
deallocate(vloc)
DEBUG_STACK_POP
end subroutine lag_zerob
!---------------------------------------------------------------------------
! SUBROUTINE: lag_zerogrnd
!---------------------------------------------------------------------------
!> Zero a Lagrange scalar field at the global grounding node
!!
!! @note The possition of this node is defined by the mesh pointer igrnd in
!! mesh.
!!
!! @param[in,out] a Field to be zeroed
!---------------------------------------------------------------------------
subroutine lag_zerogrnd(a)
class(oft_vector), intent(inout) :: a
real(r8), pointer, dimension(:) :: aloc
integer(i4) :: i,j
DEBUG_STACK_PUSH
!---Zero boundary values
NULLIFY(aloc)
CALL a%get_local(aloc)
!---
if(mesh%igrnd(1)>0)aloc(mesh%igrnd(1))=0.d0
if(mesh%igrnd(2)>0)aloc(mesh%igrnd(2))=0.d0
CALL a%restore_local(aloc)
DEALLOCATE(aloc)
DEBUG_STACK_POP
end subroutine lag_zerogrnd
! !---------------------------------------------------------------------------
! ! SUBROUTINE: blag_zerob
! !---------------------------------------------------------------------------
! !> Zero a surface Lagrange scalar field at all edge nodes
! !!
! !! @param[in,out] a Field to be zeroed
! !---------------------------------------------------------------------------
! subroutine blag_zerob(a)
! class(oft_vector), intent(inout) :: a
! integer(i4) :: i,j
! real(r8), pointer, dimension(:) :: vloc
! DEBUG_STACK_PUSH
! NULLIFY(vloc)
! !---Cast to vector type
! NULLIFY(vloc)
! call a%get_local(vloc)
! !---Zero boundary values
! !$omp parallel do
! do i=1,oft_blagrange%nbe
!   vloc(oft_blagrange%lbe(i))=0.d0
! end do
! !---
! call a%restore_local(vloc)
! deallocate(vloc)
! DEBUG_STACK_POP
! end subroutine blag_zerob
! !---------------------------------------------------------------------------
! ! SUBROUTINE: blag_zeroe
! !---------------------------------------------------------------------------
! !> Zero a surface Lagrange scalar field at all edge nodes
! !!
! !! @param[in,out] a Field to be zeroed
! !---------------------------------------------------------------------------
! subroutine blag_zeroe(a)
! class(oft_vector), intent(inout) :: a
! integer(i4) :: i,j
! real(r8), pointer, dimension(:) :: vloc
! DEBUG_STACK_PUSH
! NULLIFY(vloc)
! !---Cast to vector type
! NULLIFY(vloc)
! call a%get_local(vloc)
! !---Zero boundary values
! !$omp parallel do private(j)
! do i=1,oft_blagrange%ne
!   j=oft_blagrange%parent%le(i)
!   if(oft_lagrange%bc(j)/=3)vloc(i)=0.d0
! end do
! !---
! call a%restore_local(vloc)
! deallocate(vloc)
! DEBUG_STACK_POP
! end subroutine blag_zeroe
!---------------------------------------------------------------------------
! SUBROUTINE: lag_vzerob
!---------------------------------------------------------------------------
!> Zero a surface Lagrange vector field at all edge nodes
!!
!! @param[in,out] a Field to be zeroed
!---------------------------------------------------------------------------
subroutine lag_vzerob(a)
class(oft_vector), intent(inout) :: a
integer(i4) :: i,j
real(r8), pointer :: vloc(:,:),vtmp(:)
DEBUG_STACK_PUSH
!---Cast to vector type
ALLOCATE(vloc(3,a%n))
DO i=1,3
  vtmp=>vloc(i,:)
  call a%get_local(vtmp,i)
END DO
!---Zero boundary values
!$omp parallel do private(j)
do i=1,oft_lagrange%nbe
  j=oft_lagrange%lbe(i)
  if(oft_lagrange%global%gbe(j))vloc(:,j)=0.d0
end do
!---
DO i=1,3
  vtmp=>vloc(i,:)
  call a%restore_local(vtmp,i,wait=(i/=3))
END DO
DEALLOCATE(vloc)
DEBUG_STACK_POP
end subroutine lag_vzerob
!---------------------------------------------------------------------------
! SUBROUTINE: lag_vzeron
!---------------------------------------------------------------------------
!> Zero normal component of a Lagrange vector field at boundary nodes
!!
!! @param[in,out] a Field to be zeroed
!---------------------------------------------------------------------------
subroutine lag_vzeron(a)
class(oft_vector), intent(inout) :: a
real(r8) :: vec(3),nn(3,3)
real(r8), pointer, dimension(:) :: x,y,z
integer(i4) :: i,j
DEBUG_STACK_PUSH
!---Cast to vector type
NULLIFY(x,y,z)
CALL a%get_local(x,1)
CALL a%get_local(y,2)
CALL a%get_local(z,3)
!---Zero boundary values
!$omp parallel do private(j,vec,nn)
do i=1,oft_lagrange%nbe
  j=oft_lagrange%lbe(i)
  if(oft_lagrange%global%gbe(j))then
    CALL lag_vbc_tensor(j,1,nn)
    vec=(/x(j),y(j),z(j)/)
    vec=MATMUL(nn,vec)
    x(j)=vec(1)
    y(j)=vec(2)
    z(j)=vec(3)
  end if
end do
!---
CALL a%restore_local(x,1,wait=.TRUE.)
CALL a%restore_local(y,2,wait=.TRUE.)
CALL a%restore_local(z,3)
DEALLOCATE(x,y,z)
DEBUG_STACK_POP
end subroutine lag_vzeron
!---------------------------------------------------------------------------
! SUBROUTINE: lag_vzerot
!---------------------------------------------------------------------------
!> Zero tangential component of a Lagrange vector field at boundary nodes
!!
!! @param[in,out] a Field to be zeroed
!---------------------------------------------------------------------------
subroutine lag_vzerot(a)
class(oft_vector), intent(inout) :: a
real(r8) :: vec(3),nn(3,3)
real(r8), pointer, dimension(:) :: x,y,z
integer(i4) :: i,j
DEBUG_STACK_PUSH
!---Cast to vector type
NULLIFY(x,y,z)
CALL a%get_local(x,1)
CALL a%get_local(y,2)
CALL a%get_local(z,3)
!---Zero boundary values
!$omp parallel do private(j,vec,nn)
do i=1,oft_lagrange%nbe
  j=oft_lagrange%lbe(i)
  if(oft_lagrange%global%gbe(j))then
    CALL lag_vbc_tensor(j,2,nn)
    vec=(/x(j),y(j),z(j)/)
    vec=MATMUL(nn,vec)
    x(j)=vec(1)
    y(j)=vec(2)
    z(j)=vec(3)
  end if
end do
!---
CALL a%restore_local(x,1,wait=.TRUE.)
CALL a%restore_local(y,2,wait=.TRUE.)
CALL a%restore_local(z,3)
DEALLOCATE(x,y,z)
DEBUG_STACK_POP
end subroutine lag_vzerot
!---------------------------------------------------------------------------
! SUBROUTINE: lag_vbc_tensor
!---------------------------------------------------------------------------
!> Get boundary condition tensor for desired node in a Lagrange vector field
!!
!! @param[in] j_lag Local index of Lagrange node
!! @param[in] bc_type Desired BC (1 -> norm, 2-> tang)
!! @param[out] nn BC tensor (3,3)
!---------------------------------------------------------------------------
subroutine lag_vbc_tensor(j_lag,bc_type,nn)
integer(i4), intent(in) :: j_lag
integer(i4), intent(in) :: bc_type
real(r8), intent(out) :: nn(3,3)
integer(i4) :: i,j
DEBUG_STACK_PUSH
IF(oft_lagrange%global%gbe(j_lag))THEN
  SELECT CASE(bc_type)
    CASE(1) ! Zero normal component
      IF(oft_lagrange%bc(j_lag)==3)THEN ! On face
        nn(:,1) = (/1.d0,0.d0,0.d0/) &
        - oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(1,j_lag)
        nn(:,2) = (/0.d0,1.d0,0.d0/) &
        - oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(2,j_lag)
        nn(:,3) = (/0.d0,0.d0,1.d0/) &
        - oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(3,j_lag)
      ELSE IF(oft_lagrange%bc(j_lag)==2)THEN ! On edge
        nn(:,1) = oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(1,j_lag)
        nn(:,2) = oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(2,j_lag)
        nn(:,3) = oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(3,j_lag)
      ELSE ! At vertex
        nn=0.d0
      END IF
    CASE(2) ! Zero tangential components
      IF(oft_lagrange%bc(j_lag)==3)THEN ! On face
        nn(:,1) = oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(1,j_lag)
        nn(:,2) = oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(2,j_lag)
        nn(:,3) = oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(3,j_lag)
      ELSE ! At edge or vertex
        nn=0.d0
      END IF
  END SELECT
ELSE
  nn(:,1)=(/1.d0,0.d0,0.d0/)
  nn(:,2)=(/0.d0,1.d0,0.d0/)
  nn(:,3)=(/0.d0,0.d0,1.d0/)
END IF
DEBUG_STACK_POP
end subroutine lag_vbc_tensor
!---------------------------------------------------------------------------
! SUBROUTINE: lag_vbc_diag
!---------------------------------------------------------------------------
!> Get diagonal entries for a given boundary condition and desired node in a
!! Lagrange vector field
!!
!! @param[in] j_lag Local index of Lagrange node
!! @param[in] bc_type Desired BC (1 -> norm, 2-> tang)
!! @param[out] nn Diagonal entries (3,3)
!---------------------------------------------------------------------------
subroutine lag_vbc_diag(j_lag,bc_type,nn)
integer(i4), intent(in) :: j_lag
integer(i4), intent(in) :: bc_type
real(r8), intent(out) :: nn(3,3)
integer(i4) :: i,j
DEBUG_STACK_PUSH
SELECT CASE(bc_type)
  CASE(1) ! Get normal projections
    !---Set diagonal entries for boundary terms
    IF(oft_lagrange%bc(j_lag)==3)THEN
      nn(:,1) = oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(1,j_lag)
      nn(:,2) = oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(2,j_lag)
      nn(:,3) = oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(3,j_lag)
    ELSE IF(oft_lagrange%bc(j_lag)==2)THEN
      nn(:,1) = (/1.d0,0.d0,0.d0/) &
      - oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(1,j_lag)
      nn(:,2) = (/0.d0,1.d0,0.d0/) &
      - oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(2,j_lag)
      nn(:,3) = (/0.d0,0.d0,1.d0/) &
      - oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(3,j_lag)
    ELSE
      nn(:,1)=(/1.d0,0.d0,0.d0/)
      nn(:,2)=(/0.d0,1.d0,0.d0/)
      nn(:,3)=(/0.d0,0.d0,1.d0/)
    END IF
  CASE(2) ! Get tangential projections
    IF(oft_lagrange%bc(j_lag)==3)THEN ! On face
      nn(:,1) = (/1.d0,0.d0,0.d0/) &
      - oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(1,j_lag)
      nn(:,2) = (/0.d0,1.d0,0.d0/) &
      - oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(2,j_lag)
      nn(:,3) = (/0.d0,0.d0,1.d0/) &
      - oft_lagrange%sn(:,j_lag)*oft_lagrange%sn(3,j_lag)
    ELSE
      nn(:,1)=(/1.d0,0.d0,0.d0/)
      nn(:,2)=(/0.d0,1.d0,0.d0/)
      nn(:,3)=(/0.d0,0.d0,1.d0/)
    END IF
END SELECT
DEBUG_STACK_POP
end subroutine lag_vbc_diag
!---------------------------------------------------------------------------
! SUBROUTINE: oft_lag_getmop
!---------------------------------------------------------------------------
!> Construct mass matrix for Lagrange scalar representation
!!
!! Supported boundary conditions
!! - \c 'none' Full matrix
!! - \c 'zerob' Dirichlet for all boundary DOF
!!
!! @param[in,out] mat Matrix object
!! @param[in] bc Boundary condition
!---------------------------------------------------------------------------
subroutine oft_lag_getmop(mat,bc)
class(oft_matrix), pointer, intent(inout) :: mat
character(LEN=*), intent(in) :: bc
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time
real(r8), allocatable :: rop(:),mop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_lag_vec
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing LAG::MOP'
  CALL mytimer%tick()
END IF
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL oft_lagrange%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,rop,det,mop,curved,goptmp,m,vol,jc,jr)
allocate(j(oft_lagrange%nce)) ! Local DOF and matrix indices
allocate(rop(oft_lagrange%nce)) ! Reconstructed gradient operator
allocate(mop(oft_lagrange%nce,oft_lagrange%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  mop=0.d0
  do m=1,oft_lagrange%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call mesh%jacobian(i,oft_lagrange%quad%pts(:,m),goptmp,vol)
    det=vol*oft_lagrange%quad%wts(m)
    call oft_lag_eval_all(oft_lagrange,i,oft_lagrange%quad%pts(:,m),rop)
    !---Compute local matrix contributions
    do jr=1,oft_lagrange%nce
      do jc=1,oft_lagrange%nce
        mop(jr,jc) = mop(jr,jc) + rop(jr)*rop(jc)*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call oft_lagrange%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,oft_lagrange%nce
        IF(oft_lagrange%global%gbe(j(jr)))mop(jr,:)=0.d0
      END DO
    CASE("grnd")
      IF(ANY(mesh%igrnd>0))THEN
        DO jr=1,oft_lagrange%nce
          IF(ANY(j(jr)==mesh%igrnd))mop(jr,:)=0.d0
        END DO
      END IF
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,mop,oft_lagrange%nce,oft_lagrange%nce)
  ! !$omp end critical
end do
deallocate(j,rop,mop)
!$omp end parallel
ALLOCATE(mop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
mop(1,1)=1.d0
SELECT CASE(TRIM(bc))
  CASE("zerob")
    DO i=1,oft_lagrange%nbe
      jr=oft_lagrange%lbe(i)
      IF(oft_lagrange%global%gbe(jr).AND.oft_lagrange%linkage%leo(i))THEN
        j=jr
        call mat%add_values(j,j,mop,1,1)
      END IF
    END DO
  CASE("grnd")
    IF(ANY(mesh%igrnd>0))THEN
      DO i=1,oft_lagrange%nbe
        jr=oft_lagrange%lbe(i)
        IF(oft_lagrange%linkage%leo(i).AND.ANY(jr==mesh%igrnd))THEN
          j=jr
          call mat%add_values(j,j,mop,1,1)
        END IF
      END DO
    END IF
END SELECT
DEALLOCATE(j,mop)
CALL oft_lag_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_lag_getmop
!---------------------------------------------------------------------------
! SUBROUTINE: oft_lag_getlop
!---------------------------------------------------------------------------
!> Construct laplacian matrix for Lagrange scalar representation
!!
!! Supported boundary conditions
!! - \c 'none' Full matrix
!! - \c 'zerob' Dirichlet for all boundary DOF
!! - \c 'grnd'  Dirichlet for only groundin point
!!
!! @param[in,out] mat Matrix object
!! @param[in] bc Boundary condition
!---------------------------------------------------------------------------
subroutine oft_lag_getlop(mat,bc)
class(oft_matrix), pointer, intent(inout) :: mat
character(LEN=*), intent(in) :: bc
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time
real(r8), allocatable :: gop(:,:),lop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_lag_vec
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing LAG::LOP'
  CALL mytimer%tick()
END IF
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL oft_lagrange%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!---------------------------------------------------------------------------
! Operator integration
!---------------------------------------------------------------------------
!$omp parallel private(j,gop,det,lop,curved,goptmp,m,vol,jc,jr)
allocate(j(oft_lagrange%nce)) ! Local DOF and matrix indices
allocate(gop(3,oft_lagrange%nce)) ! Reconstructed gradient operator
allocate(lop(oft_lagrange%nce,oft_lagrange%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  lop=0.d0
  do m=1,oft_lagrange%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call mesh%jacobian(i,oft_lagrange%quad%pts(:,m),goptmp,vol)
    det=vol*oft_lagrange%quad%wts(m)
    call oft_lag_geval_all(oft_lagrange,i,oft_lagrange%quad%pts(:,m),gop,goptmp)
    !---Compute local matrix contributions
    do jr=1,oft_lagrange%nce
      do jc=1,oft_lagrange%nce
        lop(jr,jc) = lop(jr,jc) + dot_product(gop(:,jr),gop(:,jc))*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call oft_lagrange%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,oft_lagrange%nce
        IF(oft_lagrange%global%gbe(j(jr)))lop(jr,:)=0.d0
      END DO
    CASE("grnd")
      IF(ANY(mesh%igrnd>0))THEN
        DO jr=1,oft_lagrange%nce
          IF(ANY(j(jr)==mesh%igrnd))lop(jr,:)=0.d0
        END DO
      END IF
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,lop,oft_lagrange%nce,oft_lagrange%nce)
  ! !$omp end critical
end do
deallocate(j,gop,lop)
!$omp end parallel
ALLOCATE(lop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
lop(1,1)=1.d0
SELECT CASE(TRIM(bc))
  CASE("zerob")
    DO i=1,oft_lagrange%nbe
      jr=oft_lagrange%lbe(i)
      IF(oft_lagrange%global%gbe(jr).AND.oft_lagrange%linkage%leo(i))THEN
        j=jr
        call mat%add_values(j,j,lop,1,1)
      END IF
    END DO
  CASE("grnd")
    IF(ANY(mesh%igrnd>0))THEN
      DO i=1,oft_lagrange%nbe
        jr=oft_lagrange%lbe(i)
        IF(oft_lagrange%linkage%leo(i).AND.ANY(jr==mesh%igrnd))THEN
          j=jr
          call mat%add_values(j,j,lop,1,1)
        END IF
      END DO
    END IF
END SELECT
DEALLOCATE(j,lop)
CALL oft_lag_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_lag_getlop
!---------------------------------------------------------------------------
! SUBROUTINE: oft_lag_getpdop
!---------------------------------------------------------------------------
!> Construct parallel diffusion matrix for Lagrange scalar representation
!!
!! Supported boundary conditions
!! - \c 'none' Full matrix
!! - \c 'zerob' Dirichlet for all boundary DOF
!! - \c 'grnd'  Dirichlet for only groundin point
!!
!! @param[in,out] mat Matrix object
!! @param[in,out] field Vector field defining \f$ \hat{b} \f$
!! @param[in] bc Boundary condition
!! @param[in] perp Value of perpendicular conductivity (optional)
!! @param[in] be_flag Flag for dirichlet nodes if different from boundary [ne] (optional)
!---------------------------------------------------------------------------
subroutine oft_lag_getpdop(mat,field,bc,perp,be_flag)
class(oft_matrix), pointer, intent(inout) :: mat
class(fem_interp), intent(inout) :: field
character(LEN=*), intent(in) :: bc
real(r8), optional, intent(in) :: perp
logical, optional, intent(in) :: be_flag(:)
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),B_nodal(3),par_diff,perp_diff,elapsed_time
real(r8), allocatable :: gop(:,:),pdop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_lag_vec
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing LAG::PDOP'
  CALL mytimer%tick()
END IF
!---
perp_diff=0.d0
if(present(perp))perp_diff=perp
par_diff=1.d0-perp_diff
!---
IF(TRIM(bc)=="list".AND.(.NOT.PRESENT(be_flag)))CALL &
oft_abort('Boundary flag array not provided.','oft_lag_getpdop',__FILE__)
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL oft_lagrange%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,gop,det,pdop,B_nodal,curved,goptmp,m,vol,jc,jr)
allocate(j(oft_lagrange%nce)) ! Local DOF and matrix indices
allocate(gop(3,oft_lagrange%nce)) ! Reconstructed gradient operator
allocate(pdop(oft_lagrange%nce,oft_lagrange%nce)) ! Local laplacian matrix
!$omp do
do i=1,mesh%nc
  !---Get local to global DOF mapping
  call oft_lagrange%ncdofs(i,j)
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  pdop=0.d0
  do m=1,oft_lagrange%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call mesh%jacobian(i,oft_lagrange%quad%pts(:,m),goptmp,vol)
    det=vol*oft_lagrange%quad%wts(m)
    call oft_lag_geval_all(oft_lagrange,i,oft_lagrange%quad%pts(:,m),gop,goptmp)
    !---Compute bhat
    call field%interp(i,oft_lagrange%quad%pts(:,m),goptmp,B_nodal)
    B_nodal=B_nodal/SQRT(SUM(B_nodal**2))
    !---Compute local matrix contributions
    do jr=1,oft_lagrange%nce
      do jc=1,oft_lagrange%nce
        pdop(jr,jc) = pdop(jr,jc) + &
        par_diff*DOT_PRODUCT(gop(:,jr),B_nodal)*DOT_PRODUCT(gop(:,jc),B_nodal)*det + &
        perp_diff*DOT_PRODUCT(gop(:,jr),gop(:,jc))*det
      end do
    end do
  end do
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,oft_lagrange%nce
        IF(oft_lagrange%global%gbe(j(jr)))pdop(jr,:)=0.d0
      END DO
    CASE("list")
      DO jr=1,oft_lagrange%nce
        IF(be_flag(j(jr)))pdop(jr,:)=0.d0
      END DO
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,pdop,oft_lagrange%nce,oft_lagrange%nce)
  ! !$omp end critical
end do
deallocate(j,gop,pdop)
!$omp end parallel
ALLOCATE(pdop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
SELECT CASE(TRIM(bc))
  CASE("zerob")
    pdop(1,1)=1.d0
    DO i=1,oft_lagrange%nbe
      jr=oft_lagrange%lbe(i)
      IF(.NOT.oft_lagrange%global%gbe(jr))CYCLE
      IF(.NOT.oft_lagrange%linkage%leo(i))CYCLE
      j=jr
      call mat%add_values(j,j,pdop,1,1)
    END DO
  CASE("list")
    pdop(1,1)=1.d0
    DO i=1,oft_lagrange%ne
      IF(oft_lagrange%be(i))CYCLE
      IF(.NOT.be_flag(i))CYCLE
      j=i
      call mat%add_values(j,j,pdop,1,1)
    END DO
    DO i=1,oft_lagrange%nbe
      IF(.NOT.oft_lagrange%linkage%leo(i))CYCLE
      IF(.NOT.be_flag(oft_lagrange%lbe(i)))CYCLE
      j=oft_lagrange%lbe(i)
      call mat%add_values(j,j,pdop,1,1)
    END DO
END SELECT
DEALLOCATE(j,pdop)
CALL oft_lag_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_lag_getpdop
!---------------------------------------------------------------------------
! SUBROUTINE: oft_lag_project
!---------------------------------------------------------------------------
!> Project a scalar field onto a lagrange basis
!!
!! @note This subroutine only performs the integration of the field with
!! test functions for a Lagrange basis. To retrieve the correct projection the
!! result must be multiplied by the inverse of LAG::MOP.
!!
!! @param[in,out] field Scalar field for projection
!! @param[in,out] x Field projected onto Lagrange basis
!---------------------------------------------------------------------------
subroutine oft_lag_project(field,x)
class(fem_interp), intent(inout) :: field
class(oft_vector), intent(inout) :: x
!---
integer(i4) :: i,jc,m
integer(i4), allocatable :: j(:)
real(r8) :: bcc(1),vol,det,goptmp(3,4)
real(r8), pointer :: xloc(:)
real(r8), allocatable :: rop(:)
logical :: curved
DEBUG_STACK_PUSH
!---Initialize vectors to zero
NULLIFY(xloc)
call x%set(0.d0)
call x%get_local(xloc)
!---Integerate over the volume
!$omp parallel private(j,rop,curved,m,goptmp,vol,det,bcc,jc)
!---Allocate cell DOF arrays
allocate(j(oft_lagrange%nce),rop(oft_lagrange%nce))
!$omp do schedule(guided)
do i=1,mesh%nc ! Loop over cells
  call oft_lagrange%ncdofs(i,j) ! Get DOFs
  curved=cell_is_curved(mesh,i) ! Straight cell test
  do m=1,oft_lagrange%quad%np
    if(curved.OR.m==1)call mesh%jacobian(i,oft_lagrange%quad%pts(:,m),goptmp,vol)
    det=vol*oft_lagrange%quad%wts(m)
    call field%interp(i,oft_lagrange%quad%pts(:,m),goptmp,bcc)
    call oft_lag_eval_all(oft_lagrange,i,oft_lagrange%quad%pts(:,m),rop)
    do jc=1,oft_lagrange%nce
      !$omp atomic
      xloc(j(jc))=xloc(j(jc))+rop(jc)*bcc(1)*det
    end do
  end do
end do
deallocate(j,rop)
!$omp end parallel
call x%restore_local(xloc,add=.TRUE.)
deallocate(xloc)
DEBUG_STACK_POP
end subroutine oft_lag_project
!---------------------------------------------------------------------------
! SUBROUTINE: oft_lag_project_div
!---------------------------------------------------------------------------
!> Project the divergence of a scalar field onto a lagrange basis
!!
!! @note This subroutine only performs the integration of the field with
!! test functions for a Lagrange basis. To retrieve the correct projection the
!! result must be multiplied by the inverse of LAG::MOP.
!!
!! @param[in,out] field Scalar field for projection
!! @param[in,out] x Field projected onto Lagrange basis
!---------------------------------------------------------------------------
subroutine oft_lag_project_div(field,x)
class(fem_interp), intent(inout) :: field
class(oft_vector), intent(inout) :: x
!---
real(r8) :: bcc(3),det,vol,goptmp(3,4)
real(r8), pointer :: xloc(:)
real(r8), allocatable :: gop(:,:)
integer(i4) :: i,jc,m
integer(i4), allocatable :: j(:)
logical :: curved
DEBUG_STACK_PUSH
!---Initialize vectors to zero
NULLIFY(xloc)
call x%set(0.d0)
call x%get_local(xloc)
!---Integerate over the volume
!$omp parallel private(j,gop,curved,m,goptmp,vol,det,bcc,jc)
!---Allocate cell DOF arrays
allocate(j(oft_lagrange%nce),gop(3,oft_lagrange%nce))
!$omp do schedule(guided)
do i=1,mesh%nc ! Loop over cells
  call oft_lagrange%ncdofs(i,j) ! Get DOFs
  curved=cell_is_curved(mesh,i) ! Straight cell test
  do m=1,oft_lagrange%quad%np
    if(curved.OR.m==1)call mesh%jacobian(i,oft_lagrange%quad%pts(:,m),goptmp,vol)
    det=vol*oft_lagrange%quad%wts(m)
    call field%interp(i,oft_lagrange%quad%pts(:,m),goptmp,bcc)
    call oft_lag_geval_all(oft_lagrange,i,oft_lagrange%quad%pts(:,m),gop,goptmp)
    do jc=1,oft_lagrange%nce
      !$omp atomic
      xloc(j(jc))=xloc(j(jc))+DOT_PRODUCT(gop(:,jc),bcc)*det
    end do
  end do
end do
deallocate(j,gop)
!$omp end parallel
call x%restore_local(xloc,add=.TRUE.)
deallocate(xloc)
DEBUG_STACK_POP
end subroutine oft_lag_project_div
! !---------------------------------------------------------------------------
! ! SUBROUTINE: oft_blag_getmop
! !---------------------------------------------------------------------------
! !> Construct mass matrix for a boundary Lagrange scalar representation
! !!
! !! Supported boundary conditions
! !! - \c 'none' Full matrix
! !! - \c 'zerob' Dirichlet for all boundary DOF
! !!
! !! @param[in,out] mat Matrix object
! !! @param[in] bc Boundary condition
! !---------------------------------------------------------------------------
! subroutine oft_blag_getmop(mat,bc)
! class(oft_matrix), pointer, intent(inout) :: mat
! character(LEN=*), intent(in) :: bc
! integer(i4) :: i,m,jr,jc
! integer(i4), allocatable :: j(:)
! real(r8) :: vol,det,goptmp(3,4),elapsed_time
! real(r8), allocatable :: rop(:),mop(:,:)
! logical :: curved
! CLASS(oft_vector), POINTER :: oft_lag_vec
! type(oft_timer) :: mytimer
! DEBUG_STACK_PUSH
! IF(oft_debug_print(1))THEN
!   WRITE(*,'(2X,A)')'Constructing Boundary LAG::MOP'
!   CALL mytimer%tick()
! END IF
! !---------------------------------------------------------------------------
! ! Allocate matrix
! !---------------------------------------------------------------------------
! IF(.NOT.ASSOCIATED(mat))THEN
!   CALL oft_blagrange%mat_create(mat)
! ELSE
!   CALL mat%zero
! END IF
! !---------------------------------------------------------------------------
! !
! !---------------------------------------------------------------------------
! !---Operator integration loop
! !$omp parallel private(j,rop,det,mop,curved,goptmp,m,vol,jc,jr)
! allocate(j(oft_blagrange%nfe)) ! Local DOF and matrix indices
! allocate(rop(oft_blagrange%nfe)) ! Reconstructed gradient operator
! allocate(mop(oft_blagrange%nfe,oft_blagrange%nfe)) ! Local laplacian matrix
! !$omp do
! do i=1,oft_blagrange%mesh%nf
!   !---Get local reconstructed operators
!   mop=0.d0
!   do m=1,oft_blagrange%quad%np ! Loop over quadrature points
!     call oft_blagrange%mesh%jacobian(i,oft_blagrange%quad%pts(:,m),goptmp,vol)
!     det=vol*oft_blagrange%quad%wts(m)
!     do jc=1,oft_blagrange%nfe ! Loop over degrees of freedom
!       call oft_blag_eval(oft_blagrange,i,jc,oft_blagrange%quad%pts(:,m),rop(jc))
!     end do
!     !---Compute local matrix contributions
!     do jr=1,oft_blagrange%nfe
!       do jc=1,oft_blagrange%nfe
!         mop(jr,jc) = mop(jr,jc) + rop(jr)*rop(jc)*det
!       end do
!     end do
!   end do
!   !---Get local to global DOF mapping
!   call oft_blagrange%nfdofs(i,j)
!   !---Add local values to global matrix
!   ! !$omp critical
!   call mat%atomic_add_values(j,j,mop,oft_blagrange%nfe,oft_blagrange%nfe)
!   ! !$omp end critical
! end do
! deallocate(j,rop,mop)
! !$omp end parallel
! CALL oft_blagrange%vec_create(oft_lag_vec)
! CALL mat%assemble(oft_lag_vec)
! CALL oft_lag_vec%delete
! DEALLOCATE(oft_lag_vec)
! IF(oft_debug_print(1))THEN
!   elapsed_time=mytimer%tock()
!   WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
! END IF
! DEBUG_STACK_POP
! end subroutine oft_blag_getmop
! !---------------------------------------------------------------------------
! ! SUBROUTINE: oft_blag_getlop
! !---------------------------------------------------------------------------
! !> Construct laplacian matrix for Lagrange scalar representation
! !!
! !! Supported boundary conditions
! !! - \c 'none' Full matrix
! !! - \c 'zerob' Dirichlet for all boundary DOF
! !! - \c 'grnd'  Dirichlet for only groundin point
! !!
! !! @param[in,out] mat Matrix object
! !! @param[in] bc Boundary condition
! !---------------------------------------------------------------------------
! subroutine oft_blag_getlop(mat,bc)
! class(oft_matrix), pointer, intent(inout) :: mat
! character(LEN=*), intent(in) :: bc
! integer(i4) :: i,m,jr,jc
! integer(i4), allocatable :: j(:)
! real(r8) :: vol,det,goptmp(3,4),elapsed_time
! real(r8), allocatable :: gop(:,:),lop(:,:)
! logical :: curved
! CLASS(oft_vector), POINTER :: oft_lag_vec
! type(oft_timer) :: mytimer
! DEBUG_STACK_PUSH
! IF(oft_debug_print(1))THEN
!   WRITE(*,'(2X,A)')'Constructing Boundary LAG::LOP'
!   CALL mytimer%tick()
! END IF
! !---------------------------------------------------------------------------
! ! Allocate matrix
! !---------------------------------------------------------------------------
! IF(.NOT.ASSOCIATED(mat))THEN
!   CALL oft_blagrange%mat_create(mat)
! ELSE
!   CALL mat%zero
! END IF
! !---------------------------------------------------------------------------
! ! Operator integration
! !---------------------------------------------------------------------------
! !$omp parallel private(j,gop,det,lop,curved,goptmp,m,vol,jc,jr)
! allocate(j(oft_blagrange%nfe)) ! Local DOF and matrix indices
! allocate(gop(3,oft_blagrange%nfe)) ! Reconstructed gradient operator
! allocate(lop(oft_blagrange%nfe,oft_blagrange%nfe)) ! Local laplacian matrix
! !$omp do
! do i=1,oft_blagrange%mesh%nf
!   !---Get local reconstructed operators
!   lop=0.d0
!   do m=1,oft_blagrange%quad%np ! Loop over quadrature points
!     call oft_blagrange%mesh%jacobian(i,oft_blagrange%quad%pts(:,m),goptmp,vol)
!     det=vol*oft_blagrange%quad%wts(m)
!     do jc=1,oft_blagrange%nfe ! Loop over degrees of freedom
!       call oft_blag_geval(oft_blagrange,i,jc,oft_blagrange%quad%pts(:,m),gop(:,jc),goptmp)
!     end do
!     !---Compute local matrix contributions
!     do jr=1,oft_blagrange%nfe
!       do jc=1,oft_blagrange%nfe
!         lop(jr,jc) = lop(jr,jc) + DOT_PRODUCT(gop(:,jr),gop(:,jc))*det
!       end do
!     end do
!   end do
!   !---Get local to global DOF mapping
!   call oft_blagrange%nfdofs(i,j)
!   !---Apply bc to local matrix
!   SELECT CASE(TRIM(bc))
!     CASE("zerob")
!       DO jr=1,oft_blagrange%nfe
!         IF(oft_blagrange%be(j(jr)))lop(jr,:)=0.d0
!       END DO
!     CASE("edges")
!       DO jr=1,oft_blagrange%nfe
!         jc=oft_blagrange%parent%le(j(jr))
!         IF(oft_lagrange%bc(jc)/=3)lop(jr,:)=0.d0
!       END DO
!   END SELECT
!   !---Add local values to global matrix
!   ! !$omp critical
!   call mat%atomic_add_values(j,j,lop,oft_blagrange%nfe,oft_blagrange%nfe)
!   ! !$omp end critical
! end do
! deallocate(j,gop,lop)
! !$omp end parallel
! ALLOCATE(lop(1,1),j(1))
! !---Set diagonal entries for dirichlet rows
! SELECT CASE(TRIM(bc))
!   CASE("zerob")
!     lop(1,1)=1.d0
!     DO i=1,oft_blagrange%nbe
!       IF(.NOT.oft_blagrange%linkage%leo(i))CYCLE
!       j=oft_blagrange%lbe(i)
!       call mat%add_values(j,j,lop,1,1)
!     END DO
!   CASE("edges")
!     lop(1,1)=1.d0
!     DO i=1,oft_blagrange%ne
!       IF(oft_blagrange%be(i))CYCLE
!       jr=oft_blagrange%parent%le(i)
!       IF(oft_lagrange%bc(jr)/=3)THEN
!         j=i
!         call mat%add_values(j,j,lop,1,1)
!       END IF
!     END DO
!     DO i=1,oft_blagrange%nbe
!       IF(.NOT.oft_blagrange%linkage%leo(i))CYCLE
!       jc=oft_blagrange%lbe(i)
!       jr=oft_blagrange%parent%le(jc)
!       IF(oft_lagrange%bc(jr)/=3)THEN
!         j=jc
!         call mat%add_values(j,j,lop,1,1)
!       END IF
!     END DO
! END SELECT
! DEALLOCATE(j,lop)
! CALL oft_blagrange%vec_create(oft_lag_vec)
! CALL mat%assemble(oft_lag_vec)
! CALL oft_lag_vec%delete
! DEALLOCATE(oft_lag_vec)
! IF(oft_debug_print(1))THEN
!   elapsed_time=mytimer%tock()
!   WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
! END IF
! DEBUG_STACK_POP
! end subroutine oft_blag_getlop
! !------------------------------------------------------------------------------
! ! SUBROUTINE: oft_blag_nproject
! !------------------------------------------------------------------------------
! !> Project the normal component of a vector field onto a boundary Lagrange basis
! !!
! !! @note This subroutine only performs the integration of the field with
! !! boundary test functions for a Lagrange basis.
! !!
! !! @param[in,out] field Vector field for projection
! !! @param[in,out] x Field projected onto boundary Lagrange basis
! !------------------------------------------------------------------------------
! SUBROUTINE oft_blag_nproject(field,x)
! CLASS(fem_interp), INTENT(inout) :: field
! CLASS(oft_vector), INTENT(inout) :: x
! INTEGER(i4) :: i,m,k,jc,cell,ptmap(3)
! INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: j
! REAL(r8) :: vol,det,flog(4),norm(3),etmp(3),sgop(3,3),vgop(3,4),rop
! REAL(r8), POINTER, DIMENSION(:) :: xloc
! DEBUG_STACK_PUSH
! !---Initialize vectors to zero
! NULLIFY(xloc)
! call x%set(0.d0)
! call x%get_local(xloc)
! !---Operator integration loop
! !$omp parallel default(firstprivate) shared(field,xloc)
! allocate(j(oft_blagrange%nfe))
! !$omp do schedule(guided)
! do i=1,smesh%nf
!   CALL mesh%get_surf_map(i,cell,ptmap) ! Find parent cell and logical coordinate mapping
!   call oft_blagrange%nfdofs(i,j) ! Get local to global DOF mapping
!   !---Loop over quadrature points
!   do m=1,oft_blagrange%quad%np
!     call smesh%jacobian(i,oft_blagrange%quad%pts(:,m),sgop,vol)
!     call smesh%norm(i,oft_blagrange%quad%pts(:,m),norm)
!     det=vol*oft_blagrange%quad%wts(m)
!     !---Evaluate in cell coordinates
!     CALL mesh%surf_to_vol(oft_blagrange%quad%pts(:,m),ptmap,flog)
!     call mesh%jacobian(cell,flog,vgop,vol)
!     call field%interp(cell,flog,vgop,etmp)
!     !---Project on to Lagrange basis
!     do jc=1,oft_blagrange%nfe
!       call oft_blag_eval(oft_blagrange,i,jc,oft_blagrange%quad%pts(:,m),rop)
!       !$omp atomic
!       xloc(j(jc))=xloc(j(jc)) + rop*DOT_PRODUCT(etmp,norm)*det
!     end do
!   end do
! end do
! deallocate(j)
! !$omp end parallel
! call x%restore_local(xloc,add=.TRUE.)
! deallocate(xloc)
! DEBUG_STACK_POP
! END SUBROUTINE oft_blag_nproject
!---------------------------------------------------------------------------
! SUBROUTINE: oft_lag_vgetmop
!---------------------------------------------------------------------------
!> Construct mass matrix for Lagrange vector representation
!!
!! Supported boundary conditions
!! - \c 'none' Full matrix
!! - \c 'all' Dirichlet for all components at boundary
!! - \c 'norm' Dirichlet for normal component at boundary
!! - \c 'tang' Dirichlet for tangential component at boundary
!!
!! @param[in,out] mat Matrix object
!! @param[in] bc Boundary condition
!---------------------------------------------------------------------------
subroutine oft_lag_vgetmop(mat,bc)
class(oft_matrix), pointer, intent(inout) :: mat
character(LEN=*), intent(in) :: bc
type :: block_mat
  real(r8), pointer, dimension(:,:) :: m => NULL()
end type block_mat
type(block_mat) :: mtmp(3,3)
integer(i4) :: i,m,jr,jc,jp,jn,vbc_type
integer(i4), allocatable :: j(:)
real(r8) :: u,vol,det,goptmp(3,4),mloc(3,3),nn(3,3),elapsed_time
real(r8), allocatable :: rop(:)
integer(i4), allocatable, dimension(:,:) :: lcache
logical :: curved
CLASS(oft_vector), POINTER :: oft_lag_vec
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing LAG_V::MOP'
  CALL mytimer%tick()
END IF
!---Set BC flag
SELECT CASE(TRIM(bc))
  CASE("none")
    vbc_type=-1
  CASE("all")
    vbc_type=0
  CASE("norm")
    vbc_type=1
  CASE("tang")
    vbc_type=2
  CASE DEFAULT
    CALL oft_abort("Unknown BC requested","oft_lag_vgetmop",__FILE__)
END SELECT
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL oft_vlagrange%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,rop,det,mtmp,nn,mloc,curved,goptmp,m,u,vol,jc,jr,lcache)
allocate(j(oft_lagrange%nce)) ! Local DOF and matrix indices
allocate(rop(oft_lagrange%nce)) ! Reconstructed gradient operator
allocate(lcache(oft_lagrange%nce,oft_lagrange%nce))
DO jr=1,3
  DO jc=1,3
    allocate(mtmp(jr,jc)%m(oft_lagrange%nce,oft_lagrange%nce)) ! Local laplacian matrix
  END DO
END DO
!$omp do schedule(guided)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local to global DOF mapping
  CALL oft_lagrange%ncdofs(i,j)
  !---Compute local matrix contributions
  DO jp=1,3
    DO jn=1,3
      IF(vbc_type<=0.AND.jp/=jn)CYCLE
      mtmp(jp,jn)%m = 0.d0
    END DO
  END DO
  !---Get local reconstructed operators
  DO m=1,oft_lagrange%quad%np ! Loop over quadrature points
    IF(curved.OR.m==1)CALL mesh%jacobian(i,oft_lagrange%quad%pts(:,m),goptmp,vol)
    det=vol*oft_lagrange%quad%wts(m)
    call oft_lag_eval_all(oft_lagrange,i,oft_lagrange%quad%pts(:,m),rop)
    !---Compute local matrix contributions
    DO jc=1,oft_lagrange%nce
      u=rop(jc)*det
      DO jr=1,oft_lagrange%nce
        DO jp=1,3
          mtmp(jp,jp)%m(jr,jc) = mtmp(jp,jp)%m(jr,jc) + rop(jr)*u
        END DO
      END DO
    END DO
  END DO
  !---Compute BC projection matrices
  DO jr=1,oft_lagrange%nce
    IF(oft_lagrange%global%gbe(j(jr)))THEN
      IF(vbc_type==0)THEN
        DO jp=1,3
          mtmp(jp,jp)%m(jr,:)=0.d0
        END DO
      ELSE IF(vbc_type>0)THEN
        DO jc=1,oft_lagrange%nce
          mloc=0.d0
          u=mtmp(1,1)%m(jr,jc)
          DO jp=1,3
            mloc(jp,jp)=u
          END DO
          CALL lag_vbc_tensor(j(jr),vbc_type,nn)
          mloc=MATMUL(nn,mloc)
          DO jp=1,3
            DO jn=1,3
              mtmp(jp,jn)%m(jr,jc) = mloc(jp,jn)
            END DO
          END DO
        END DO
      END IF
    END IF
  END DO
  ! !$omp critical
  !---Add local values to global matrix
  lcache(1,1)=0
  DO jp=1,3
    DO jn=1,3
      IF(vbc_type<=0.AND.jp/=jn)CYCLE
      CALL mat%atomic_add_values(j,j,mtmp(jp,jn)%m,oft_lagrange%nce,oft_lagrange%nce,jp,jn,lcache)
    END DO
  END DO
  ! !$omp end critical
END DO
!---Deallocate thread local variables
DO jr=1,3
  DO jc=1,3
    deallocate(mtmp(jr,jc)%m)
  END DO
END DO
deallocate(j,rop,lcache)
!$omp end parallel
ALLOCATE(j(1),rop(1))
!---Set diagonal entries for dirichlet rows
IF(vbc_type==0)THEN
  rop=1.d0
  DO i=1,oft_lagrange%nbe
    IF(.NOT.oft_lagrange%linkage%leo(i))CYCLE
    jr=oft_lagrange%lbe(i)
    IF(.NOT.oft_lagrange%global%gbe(jr))CYCLE
    j=jr
    CALL mat%add_values(j,j,rop,1,1,1,1)
    CALL mat%add_values(j,j,rop,1,1,2,2)
    CALL mat%add_values(j,j,rop,1,1,3,3)
  END DO
ELSE IF(vbc_type>0)THEN
  DO i=1,oft_lagrange%nbe
    IF(.NOT.oft_lagrange%linkage%leo(i))CYCLE
    jr=oft_lagrange%lbe(i)
    IF(.NOT.oft_lagrange%global%gbe(jr))CYCLE
    j=jr
    CALL lag_vbc_diag(jr,vbc_type,mloc)
    DO jr=1,3
      DO jc=1,3
        call mat%add_values(j,j,mloc(jr,jc),1,1,jr,jc)
      END DO
    END DO
  END DO
END IF
DEALLOCATE(j,rop)
CALL oft_lag_vcreate(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_lag_vgetmop
!---------------------------------------------------------------------------
! SUBROUTINE: oft_lag_vproject
!---------------------------------------------------------------------------
!> Project a vector field onto a lagrange basis
!!
!! @note This subroutine only performs the integration of the field with
!! test functions for a Lagrange basis. To retrieve the correct projection the
!! result must be multiplied by the inverse of LAG::VMOP.
!!
!! @param[in,out] field Vector field for projection
!! @param[in,out] x Field projected onto Lagrange basis
!---------------------------------------------------------------------------
subroutine oft_lag_vproject(field,x)
class(fem_interp), intent(inout) :: field
class(oft_vector), intent(inout) :: x
!---
real(r8) :: bcc(3),det,goptmp(3,4),vol
real(r8), pointer, dimension(:) :: xloc,yloc,zloc
real(r8), allocatable :: rop(:)
integer(i4) :: i,jc,m
integer(i4), allocatable :: j(:)
logical :: curved
DEBUG_STACK_PUSH
!---Initialize vectors to zero
NULLIFY(xloc,yloc,zloc)
call x%set(0.d0)
call x%get_local(xloc,1)
call x%get_local(yloc,2)
call x%get_local(zloc,3)
!---Integerate over the volume
!$omp parallel private(j,rop,curved,m,goptmp,vol,det,bcc,jc)
!---Allocate cell DOF arrays
allocate(j(oft_lagrange%nce),rop(oft_lagrange%nce))
!$omp do schedule(guided)
do i=1,mesh%nc ! Loop over cells
  call oft_lagrange%ncdofs(i,j) ! Get DOFs
  curved=cell_is_curved(mesh,i) ! Straight cell test
  do m=1,oft_lagrange%quad%np
    if(curved.OR.m==1)call mesh%jacobian(i,oft_lagrange%quad%pts(:,m),goptmp,vol)
    det=vol*oft_lagrange%quad%wts(m)
    call field%interp(i,oft_lagrange%quad%pts(:,m),goptmp,bcc)
    call oft_lag_eval_all(oft_lagrange,i,oft_lagrange%quad%pts(:,m),rop)
    do jc=1,oft_lagrange%nce
      !$omp atomic
      xloc(j(jc))=xloc(j(jc))+rop(jc)*bcc(1)*det
      !$omp atomic
      yloc(j(jc))=yloc(j(jc))+rop(jc)*bcc(2)*det
      !$omp atomic
      zloc(j(jc))=zloc(j(jc))+rop(jc)*bcc(3)*det
    end do
  end do
end do
deallocate(j,rop)
!$omp end parallel
call x%restore_local(xloc,1,add=.TRUE.,wait=.TRUE.)
call x%restore_local(yloc,2,add=.TRUE.,wait=.TRUE.)
call x%restore_local(zloc,3,add=.TRUE.)
deallocate(xloc,yloc,zloc)
DEBUG_STACK_POP
end subroutine oft_lag_vproject
!---------------------------------------------------------------------------
! SUBROUTINE: lag_div
!---------------------------------------------------------------------------
!> Compute the divergence of a Lagrange vector field.
!!
!! @param[in,out] a Input field
!! @param[out] reg \f$ \int_v \nabla \cdot a \; dV \f$
!---------------------------------------------------------------------------
subroutine lag_div(a,reg)
class(oft_vector), intent(inout) :: a
real(r8), intent(out) :: reg
real(r8), pointer, dimension(:) :: x,y,z
integer :: i,jc,m
integer(i4), allocatable :: j(:)
real(r8) :: goptmp(3,4),vol,det,div
real(r8), allocatable :: gop(:,:)
logical :: curved
DEBUG_STACK_PUSH
!---Cast to vector type
NULLIFY(x,y,z)
CALL a%get_local(x,1)
CALL a%get_local(y,2)
CALL a%get_local(z,3)
!---
reg=0.d0
allocate(j(oft_lagrange%nce),gop(3,oft_lagrange%nce))
do i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local to global DOF mapping
  call oft_lagrange%ncdofs(i,j)
  !---Get local reconstructed operators
  do m=1,oft_lagrange%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call mesh%jacobian(i,oft_lagrange%quad%pts(:,m),goptmp,vol)
    det=vol*oft_lagrange%quad%wts(m)
    call oft_lag_geval_all(oft_lagrange,i,oft_lagrange%quad%pts(:,m),gop,goptmp)
    div=0.d0
    do jc=1,oft_lagrange%nce ! Loop over degrees of freedom
      div=div+x(j(jc))*gop(1,jc)
      div=div+y(j(jc))*gop(2,jc)
      div=div+z(j(jc))*gop(3,jc)
    end do
    reg=reg+(div**2)*det
  end do
end do
deallocate(j,gop)
reg=oft_mpi_sum(reg)
DEALLOCATE(x,y,z)
DEBUG_STACK_POP
end subroutine lag_div
!---------------------------------------------------------------------------
! SUBROUTINE: lag_setup_interp
!---------------------------------------------------------------------------
!> Construct interpolation matrices on each MG level
!---------------------------------------------------------------------------
SUBROUTINE lag_setup_interp(create_vec)
LOGICAL, OPTIONAL, INTENT(in) :: create_vec
!---
class(oft_vector), pointer :: fvec,cvec
type(oft_graph_ptr) :: graphs(3,3)
type(oft_matrix_ptr) :: mats(3,3)
INTEGER(i4) :: i
LOGICAL :: vec_interp
DEBUG_STACK_PUSH
vec_interp=.FALSE.
IF(PRESENT(create_vec))vec_interp=create_vec
!---
DO i=oft_lagrange_minlev+1,oft_lagrange_nlevels
  CALL oft_lag_set_level(i)
  !---
  if(oft_lagrange_level==oft_lagrange_blevel+1)then
    CYCLE
  end if
  !---Setup interpolation
  if(oft_lagrange%order==1)then
    CALL lag_ginterpmatrix(ML_oft_lagrange%interp_matrices(ML_oft_lagrange%level)%m)
    oft_lagrange_ops%interp=>ML_oft_lagrange%interp_matrices(ML_oft_lagrange%level)%m
    CALL oft_lagrange_ops%interp%assemble
  else
    CALL lag_pinterpmatrix(ML_oft_lagrange%interp_matrices(ML_oft_lagrange%level)%m)
    oft_lagrange_ops%interp=>ML_oft_lagrange%interp_matrices(ML_oft_lagrange%level)%m
    CALL oft_lagrange_ops%interp%assemble
  end if
END DO
!---Create vector interpolation operator
IF(vec_interp)THEN
  CALL ML_oft_vlagrange%build_interp
  DO i=oft_lagrange_minlev+1,oft_lagrange_nlevels
    CALL oft_lag_set_level(i)
    !---
    if(oft_lagrange_level==oft_lagrange_blevel+1)then
      CYCLE
    end if
    oft_lagrange_ops%vinterp=>ML_oft_vlagrange%interp_matrices(i)%m
  END DO
END IF
DEBUG_STACK_POP
END SUBROUTINE lag_setup_interp
!---------------------------------------------------------------------------
! SUBROUTINE: lag_ginterpmatrix
!---------------------------------------------------------------------------
!> Construct interpolation matrix for polynomial levels
!---------------------------------------------------------------------------
SUBROUTINE lag_ginterpmatrix(mat)
class(oft_matrix), pointer, intent(inout) :: mat
INTEGER(i4) :: i,j,k,m,icors,ifine,jb,i_ind(1),j_ind(1)
INTEGER(i4) :: etmp(2),ftmp(3),fetmp(3),ctmp(4),fc,ed
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pmap,emap,fmap
REAL(r8) :: f(4),val,mop(1)
REAL(r8), POINTER, DIMENSION(:,:) :: ed_nodes,fc_nodes,c_nodes
CLASS(oft_afem_type), POINTER :: lag_cors => NULL()
TYPE(oft_lag_ops), POINTER :: ops
class(oft_mesh), pointer :: cmesh
CLASS(oft_vector), POINTER :: lag_vec,lag_vec_cors
type(oft_graph_ptr), pointer :: graphs(:,:)
DEBUG_STACK_PUSH
!---
if(mg_mesh%level<1)then
  call oft_abort('Invalid mesh level','lag_ginterpmatrix',__FILE__)
end if
cmesh=>mg_mesh%meshes(mg_mesh%level-1)
if(cmesh%type/=1)CALL oft_abort("Only supported with tet meshes", &
  "lag_ginterpmatrix", __FILE__)
if(oft_lagrange%order/=1)then
  call oft_abort('Attempted geometric interpolation for pd > 1','lag_ginterpmatrix',__FILE__)
end if
ops=>oft_lagrange_ops
lag_cors=>ML_oft_lagrange%levels(oft_lagrange_level-1)%fe
ALLOCATE(ML_oft_lagrange%interp_graphs(ML_oft_lagrange%level)%g)
ops%interp_graph=>ML_oft_lagrange%interp_graphs(ML_oft_lagrange%level)%g
!---Setup matrix sizes
ops%interp_graph%nr=oft_lagrange%ne
ops%interp_graph%nrg=oft_lagrange%global%ne
ops%interp_graph%nc=lag_cors%ne
ops%interp_graph%ncg=lag_cors%global%ne
!---Setup Matrix graph
ALLOCATE(ops%interp_graph%kr(ops%interp_graph%nr+1))
ops%interp_graph%nnz=cmesh%np+2*cmesh%ne
ALLOCATE(ops%interp_graph%lc(ops%interp_graph%nnz))
ops%interp_graph%lc=0_i4
!---Construct linkage
DO i=1,cmesh%np
  ops%interp_graph%kr(i)=1
END DO
DO i=1,cmesh%ne
  ops%interp_graph%kr(i+cmesh%np)=2
END DO
ops%interp_graph%kr(ops%interp_graph%nr+1)=ops%interp_graph%nnz+1
do i=ops%interp_graph%nr,1,-1 ! cumulative point to point count
  ops%interp_graph%kr(i)=ops%interp_graph%kr(i+1)-ops%interp_graph%kr(i)
end do
if(ops%interp_graph%kr(1)/=1)call oft_abort('Bad element to element count','oft_lag_ginterpmatrix',__FILE__)
DO i=1,cmesh%np
  ops%interp_graph%lc(ops%interp_graph%kr(i))=i
END DO
!---
DO i=1,cmesh%ne
  etmp=cmesh%le(:,i)
  !---
  ifine = cmesh%np+i
  jb=ops%interp_graph%kr(ifine)-1
  DO k=1,2
    ops%interp_graph%lc(jb+k)=etmp(k)
  END DO
END DO
!---------------------------------------------------------------------------
! Construct matrix
!---------------------------------------------------------------------------
CALL oft_lag_create(lag_vec)
CALL oft_lag_create(lag_vec_cors,oft_lagrange_level-1)
!---
ALLOCATE(graphs(1,1))
graphs(1,1)%g=>ops%interp_graph
!---
CALL create_matrix(mat,graphs,lag_vec,lag_vec_cors)
CALL lag_vec%delete
CALL lag_vec_cors%delete
DeALLOCATE(graphs,lag_vec,lag_vec_cors)
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
  etmp=cmesh%le(:,i)
  !---
  ifine = cmesh%np+i
  jb=ops%interp_graph%kr(ifine)-1
  DO k=1,2
    i_ind=ifine
    j_ind=etmp(k)
    mop=1.d0/2.d0
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
  END DO
END DO
deallocate(emap)
DEBUG_STACK_POP
END SUBROUTINE lag_ginterpmatrix
!---------------------------------------------------------------------------
! SUBROUTINE: lag_pinterpmatrix
!---------------------------------------------------------------------------
!> Construct interpolation matrix for polynomial levels
!---------------------------------------------------------------------------
SUBROUTINE lag_pinterpmatrix(mat)
class(oft_matrix), pointer, intent(inout) :: mat
INTEGER(i4) :: i,j,k,m,icors,ifine,jb,js,jn,i_ind(1),j_ind(1)
INTEGER(i4) :: etmp(2),fc,ed,cell,dof,offset
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pmap,emap,fmap,jcors,ftmp,fetmp,ctmp
REAL(r8) :: f(4),val,mop(1)
REAL(r8), POINTER, DIMENSION(:,:) :: ed_nodes,fc_nodes,c_nodes
CLASS(oft_scalar_fem), POINTER :: lag_cors => NULL()
TYPE(oft_lag_ops), POINTER :: ops
CLASS(oft_vector), POINTER :: lag_vec,lag_vec_cors
type(oft_graph_ptr), pointer :: graphs(:,:)
DEBUG_STACK_PUSH
allocate(ftmp(mesh%face_np),fetmp(mesh%face_np),ctmp(mesh%cell_np))
!---
ops=>oft_lagrange_ops
SELECT TYPE(this=>ML_oft_lagrange%levels(oft_lagrange_level-1)%fe)
  CLASS IS(oft_scalar_fem)
    lag_cors=>this
  CLASS DEFAULT
    CALL oft_abort("Error casting coarse level", "lag_pinterpmatrix", __FILE__)
END SELECT
! lag_cors=>ML_oft_lagrange%levels(oft_lagrange_level-1)%fe
!WRITE(*,*)lag_cors%gstruct
!WRITE(*,*)oft_lagrange%gstruct
CALL oft_lag_nodes(oft_lagrange%order,ed_nodes,fc_nodes,c_nodes)
ALLOCATE(ML_oft_lagrange%interp_graphs(ML_oft_lagrange%level)%g)
ops%interp_graph=>ML_oft_lagrange%interp_graphs(ML_oft_lagrange%level)%g
!---Setup matrix sizes
ops%interp_graph%nr=oft_lagrange%ne
ops%interp_graph%nrg=oft_lagrange%global%ne
ops%interp_graph%nc=lag_cors%ne
ops%interp_graph%ncg=lag_cors%global%ne
!---Setup Matrix graph
ALLOCATE(ops%interp_graph%kr(ops%interp_graph%nr+1))
ops%interp_graph%nnz=mesh%np
ops%interp_graph%nnz=ops%interp_graph%nnz + &
  (2+lag_cors%gstruct(2))*mesh%ne*oft_lagrange%gstruct(2)
ops%interp_graph%nnz=ops%interp_graph%nnz + &
  (mesh%face_np+lag_cors%gstruct(2)*mesh%face_np &
  + lag_cors%gstruct(3))*mesh%nf*oft_lagrange%gstruct(3)
ops%interp_graph%nnz=ops%interp_graph%nnz + &
  (mesh%cell_np+lag_cors%gstruct(2)*mesh%cell_ne &
  + lag_cors%gstruct(3)*mesh%cell_nf &
  + lag_cors%gstruct(4))*mesh%nc*oft_lagrange%gstruct(4)
ALLOCATE(ops%interp_graph%lc(ops%interp_graph%nnz))
ops%interp_graph%lc=0_i4
!---Construct matrix
DO i=1,mesh%np
  ops%interp_graph%kr(i)=1
END DO
!---
DO i=1,mesh%ne
  DO j=1,oft_lagrange%gstruct(2)
    ifine = mesh%np + (i-1)*oft_lagrange%gstruct(2) + j
    ops%interp_graph%kr(ifine)=2+lag_cors%gstruct(2)
  END DO
END DO
!---
DO i=1,mesh%nf
  DO j=1,oft_lagrange%gstruct(3)
    ifine = mesh%np+oft_lagrange%gstruct(2)*mesh%ne &
      + (i-1)*oft_lagrange%gstruct(3) + j
    ops%interp_graph%kr(ifine)= mesh%face_np &
      + mesh%face_np*lag_cors%gstruct(2)+lag_cors%gstruct(3)
  END DO
END DO
!---
DO i=1,mesh%nc
  DO j=1,oft_lagrange%gstruct(4)
    ifine = mesh%np + oft_lagrange%gstruct(2)*mesh%ne &
      + oft_lagrange%gstruct(3)*mesh%nf + (i-1)*oft_lagrange%gstruct(4) + j
    ops%interp_graph%kr(ifine) = mesh%cell_np &
      + mesh%cell_ne*lag_cors%gstruct(2) &
      + mesh%cell_nf*lag_cors%gstruct(3) + lag_cors%gstruct(4)
  END DO
END DO
ops%interp_graph%kr(ops%interp_graph%nr+1)=ops%interp_graph%nnz+1
do i=ops%interp_graph%nr,1,-1 ! cumulative point to point count
  ops%interp_graph%kr(i)=ops%interp_graph%kr(i+1)-ops%interp_graph%kr(i)
end do
if(ops%interp_graph%kr(1)/=1)call oft_abort('Bad element to element count','oft_lag_interpmatrix',__FILE__)
!---Construct matrix
DO i=1,mesh%np
  ops%interp_graph%lc(ops%interp_graph%kr(i))=i
END DO
!---
DO i=1,mesh%ne
  etmp=mesh%le(:,i)
  DO j=1,oft_lagrange%gstruct(2)
    ifine = mesh%np + (i-1)*oft_lagrange%gstruct(2) + j
    jb=ops%interp_graph%kr(ifine)-1
    DO k=1,2
      ops%interp_graph%lc(jb+k)=etmp(k)
    END DO
    DO k=1,lag_cors%gstruct(2)
      ops%interp_graph%lc(jb+2+k)=mesh%np + (i-1)*lag_cors%gstruct(2) + k
    END DO
    !---
    js=ops%interp_graph%kr(ifine)
    jn=ops%interp_graph%kr(ifine+1)-1
    CALL sort_array(ops%interp_graph%lc(js:jn),jn-js+1)
  END DO
END DO
!---
DO i=1,mesh%nf
  ftmp=mesh%lf(:,i)
  fetmp=mesh%lfe(:,i)
  DO j=1,oft_lagrange%gstruct(3)
    ifine = mesh%np+oft_lagrange%gstruct(2)*mesh%ne &
      + (i-1)*oft_lagrange%gstruct(3) + j
    jb=ops%interp_graph%kr(ifine)-1
    DO k=1,mesh%face_np
      ops%interp_graph%lc(jb+k)=ftmp(k)
      ! etmp=mesh%bmesh%face_ed(:,k)
      ed=fetmp(k)
      DO m=1,lag_cors%gstruct(2)
        ops%interp_graph%lc(jb+mesh%face_np+(m-1)*mesh%face_np+k)= &
          mesh%np + (ABS(ed)-1)*lag_cors%gstruct(2) + m
      END DO
    END DO
    offset=jb+mesh%face_np+mesh%face_np*lag_cors%gstruct(2)
    DO m=1,lag_cors%gstruct(3)
      ops%interp_graph%lc(offset+m)= &
        mesh%np + lag_cors%gstruct(2)*mesh%ne + (i-1)*lag_cors%gstruct(3) + m
    END DO
    !---
    js=ops%interp_graph%kr(ifine)
    jn=ops%interp_graph%kr(ifine+1)-1
    CALL sort_array(ops%interp_graph%lc(js:jn),jn-js+1)
  END DO
END DO
!---
DO i=1,mesh%nc
  ctmp=mesh%lc(:,i)
  DO j=1,oft_lagrange%gstruct(4)
    ifine = mesh%np + oft_lagrange%gstruct(2)*mesh%ne &
      + oft_lagrange%gstruct(3)*mesh%nf + (i-1)*oft_lagrange%gstruct(4) + j
    jb=ops%interp_graph%kr(ifine)-1
    DO k=1,mesh%cell_np
      ops%interp_graph%lc(jb+k)=ctmp(k)
    END DO
    offset=jb+mesh%cell_np
    DO k=1,mesh%cell_ne
      etmp=mesh%cell_ed(:,k)
      ed=mesh%lce(k,i)
      DO m=1,lag_cors%gstruct(2)
        ops%interp_graph%lc(offset+(m-1)*mesh%cell_ne+k)= &
          mesh%np + (ABS(ed)-1)*lag_cors%gstruct(2) + m
      END DO
    END DO
    offset=jb+mesh%cell_np+mesh%cell_ne*lag_cors%gstruct(2)
    DO k=1,mesh%cell_nf
      ftmp=mesh%cell_fc(:,k)
      DO m=1,lag_cors%gstruct(3)
        ops%interp_graph%lc(offset+(m-1)*mesh%cell_nf+k)=mesh%np &
        + lag_cors%gstruct(2)*mesh%ne + (ABS(mesh%lcf(k,i))-1)*lag_cors%gstruct(3) + m
      END DO
    END DO
    offset=jb+mesh%cell_np+mesh%cell_ne*lag_cors%gstruct(2)+mesh%cell_nf*lag_cors%gstruct(3)
    DO m=1,lag_cors%gstruct(4)
      ops%interp_graph%lc(offset+m)=mesh%np &
      + lag_cors%gstruct(2)*mesh%ne+lag_cors%gstruct(3)*mesh%nf &
      + (i-1)*lag_cors%gstruct(4) + m
    END DO
    !---
    js=ops%interp_graph%kr(ifine)
    jn=ops%interp_graph%kr(ifine+1)-1
    CALL sort_array(ops%interp_graph%lc(js:jn),jn-js+1)
  END DO
END DO
!---------------------------------------------------------------------------
! Construct matrix
!---------------------------------------------------------------------------
CALL oft_lag_create(lag_vec)
CALL oft_lag_create(lag_vec_cors,oft_lagrange_level-1)
!---
ALLOCATE(graphs(1,1))
graphs(1,1)%g=>ops%interp_graph
!---
CALL create_matrix(mat,graphs,lag_vec,lag_vec_cors)
CALL lag_vec%delete
CALL lag_vec_cors%delete
DeALLOCATE(graphs,lag_vec,lag_vec_cors)
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
allocate(pmap(mesh%np))
CALL get_inverse_map(mesh%lbp,mesh%nbp,pmap,mesh%np)
DO i=1,mesh%np
  IF(mesh%bp(i))THEN
    ! IF(.NOT.mesh%linkage%lpo(pmap(i)))CYCLE
    IF(.NOT.mesh%pstitch%leo(pmap(i)))CYCLE
  END IF
  i_ind=i
  j_ind=i
  mop=1.d0
  CALL mat%add_values(i_ind,j_ind,mop,1,1)
END DO
deallocate(pmap)
!---
allocate(emap(mesh%ne))
CALL get_inverse_map(mesh%lbe,mesh%nbe,emap,mesh%ne)
DO i=1,mesh%ne
  IF(mesh%be(i))THEN
    ! IF(.NOT.mesh%linkage%leo(emap(i)))CYCLE
    IF(.NOT.mesh%estitch%leo(emap(i)))CYCLE
  END IF
  cell=mesh%lec(mesh%kec(i))
  DO ed=1,mesh%cell_ne
    IF(ABS(mesh%lce(ed,cell))==i)EXIT
  END DO
  IF(ed>mesh%cell_ne)CALL oft_abort("Could not find edge", "", __FILE__)
  DO j=1,oft_lagrange%gstruct(2)
    ifine = mesh%np + (i-1)*oft_lagrange%gstruct(2) + j
    dof = mesh%cell_np + (ed-1)*oft_lagrange%gstruct(2) + j
    CALL oft_lag_npos(oft_lagrange,cell,dof,f)
    DO k=1,2
      dof=mesh%cell_ed(k,ed)
      CALL oft_lag_eval(lag_cors,cell,dof,f,val)
      i_ind=ifine
      j_ind=mesh%lc(dof,cell)
      mop=val
      CALL mat%add_values(i_ind,j_ind,mop,1,1)
    END DO
    DO k=1,lag_cors%gstruct(2)
      dof = mesh%cell_np + (ed-1)*lag_cors%gstruct(2) + k
      CALL oft_lag_eval(lag_cors,cell,dof,f,val)
      i_ind=ifine
      j_ind=mesh%np + (i-1)*lag_cors%gstruct(2) + k
      mop=val
      CALL mat%add_values(i_ind,j_ind,mop,1,1)
    END DO
  END DO
END DO
deallocate(emap)
!---
allocate(fmap(mesh%nf))
CALL get_inverse_map(mesh%lbf,mesh%nbf,fmap,mesh%nf)
DO i=1,mesh%nf
  IF(mesh%bf(i))THEN
    ! IF(.NOT.mesh%linkage%lfo(fmap(i)))CYCLE
    IF(.NOT.mesh%fstitch%leo(fmap(i)))CYCLE
  END IF
  cell=mesh%lfc(1,i)
  DO fc=1,mesh%cell_nf
    IF(ABS(mesh%lcf(fc,cell))==i)EXIT
  END DO
  IF(fc>mesh%cell_nf)CALL oft_abort("Could not find face", "", __FILE__)
  DO j=1,oft_lagrange%gstruct(3)
    ifine = mesh%np + oft_lagrange%gstruct(2)*mesh%ne &
          + (i-1)*oft_lagrange%gstruct(3) + j
    dof = mesh%cell_np + oft_lagrange%gstruct(2)*mesh%cell_ne &
          + (fc-1)*oft_lagrange%gstruct(3) + j
    CALL oft_lag_npos(oft_lagrange,cell,dof,f)
    DO k=1,mesh%face_np
      dof=mesh%cell_fc(k,fc)
      CALL oft_lag_eval(lag_cors,cell,dof,f,val)
      i_ind=ifine
      j_ind=mesh%lc(dof,cell)
      mop=val
      CALL mat%add_values(i_ind,j_ind,mop,1,1)
      DO m=1,lag_cors%gstruct(2)
        ed=ABS(mesh%cell_fe(k,fc))
        dof = mesh%cell_np + (ed-1)*lag_cors%gstruct(2) + m
        CALL oft_lag_eval(lag_cors,cell,dof,f,val)
        i_ind=ifine
        j_ind=mesh%np + (ABS(mesh%lce(ed,cell))-1)*lag_cors%gstruct(2)+m
        mop=val
        CALL mat%add_values(i_ind,j_ind,mop,1,1)
      END DO
    END DO
    DO k=1,lag_cors%gstruct(3)
      dof = mesh%cell_np + lag_cors%gstruct(2)*mesh%cell_ne &
          + (fc-1)*lag_cors%gstruct(3) + k
      CALL oft_lag_eval(lag_cors,cell,dof,f,val)
      i_ind=ifine
      j_ind=mesh%np+lag_cors%gstruct(2)*mesh%ne+(i-1)*lag_cors%gstruct(3)+k
      mop=val
      CALL mat%add_values(i_ind,j_ind,mop,1,1)
    END DO
  END DO
END DO
deallocate(fmap)
!---
allocate(jcors(lag_cors%nce))
DO i=1,mesh%nc
  CALL lag_cors%ncdofs(i, jcors)
  DO j=1,oft_lagrange%gstruct(4)
    ifine = mesh%np + oft_lagrange%gstruct(2)*mesh%ne &
          + oft_lagrange%gstruct(3)*mesh%nf + (i-1)*oft_lagrange%gstruct(4) + j
    dof = mesh%cell_np + oft_lagrange%gstruct(2)*mesh%cell_ne &
          + oft_lagrange%gstruct(3)*mesh%cell_nf &
          + j
    CALL oft_lag_npos(oft_lagrange,cell,dof,f)
    DO k=1,lag_cors%nce
      CALL oft_lag_eval(lag_cors,i,k,f,val)
      i_ind=ifine
      j_ind=jcors(k)
      mop=val
      CALL mat%add_values(i_ind,j_ind,mop,1,1)
    END DO
  END DO
END DO
!---
DEALLOCATE(ed_nodes,fc_nodes,c_nodes)
DEALLOCATE(ftmp,fetmp,ctmp)
DEBUG_STACK_POP
END SUBROUTINE lag_pinterpmatrix
!---------------------------------------------------------------------------
! SUBROUTINE: lag_interp
!---------------------------------------------------------------------------
!> Interpolate a coarse level Lagrange scalar field to the next finest level
!!
!! @note The global Lagrange level in incremented by one in this subroutine
!!
!! @param[in] acors Vector to interpolate
!! @param[in,out] afine Fine vector from interpolation
!---------------------------------------------------------------------------
subroutine lag_interp(acors,afine)
class(oft_vector), intent(inout) :: acors
class(oft_vector), intent(inout) :: afine
DEBUG_STACK_PUSH
!---Step one level up
call oft_lag_set_level(oft_lagrange_level+1)
call afine%set(0.d0)
!---
if(oft_lagrange_level==oft_lagrange_blevel+1)then
  call lag_base_pop(acors,afine)
  DEBUG_STACK_POP
  return
end if
CALL oft_lagrange_ops%interp%apply(acors,afine)
DEBUG_STACK_POP
end subroutine lag_interp
!---------------------------------------------------------------------------
! SUBROUTINE: lag_vinterp
!---------------------------------------------------------------------------
!> Interpolate a coarse level Lagrange scalar field to the next finest level
!!
!! @note The global Lagrange level in incremented by one in this subroutine
!!
!! @param[in] acors Vector to interpolate
!! @param[in,out] afine Fine vector from interpolation
!---------------------------------------------------------------------------
subroutine lag_vinterp(acors,afine)
class(oft_vector), intent(inout) :: acors
class(oft_vector), intent(inout) :: afine
DEBUG_STACK_PUSH
!---Step one level up
call oft_lag_set_level(oft_lagrange_level+1)
!---
if(oft_lagrange_level==oft_lagrange_blevel+1)then
  CALL oft_abort('MPI-Base transfer not supported','lag_vinterp',__FILE__)
!  call lag_base_pop(acors,afine)
  DEBUG_STACK_POP
  return
end if
CALL oft_lagrange_ops%vinterp%apply(acors,afine)
DEBUG_STACK_POP
end subroutine lag_vinterp
!---------------------------------------------------------------------------
! SUBROUTINE: lag_base_pop
!---------------------------------------------------------------------------
!> Transfer a base level Lagrange scalar field to the next MPI level
!!
!! @param[in] acors Vector to transfer
!! @param[in,out] afine Fine vector from transfer
!---------------------------------------------------------------------------
subroutine lag_base_pop(acors,afine)
class(oft_vector), intent(inout) :: acors
class(oft_vector), intent(inout) :: afine
integer(i4), pointer, dimension(:) :: lptmp
integer(i4) :: i
real(r8), pointer, dimension(:) :: array_c,array_f
DEBUG_STACK_PUSH
lptmp=>mg_mesh%meshes(mg_mesh%nbase+1)%base%lp
CALL acors%get_local(array_c)
CALL afine%get_local(array_f)
!$omp parallel do
do i=1,afine%n
  array_f(i)=array_c(lptmp(i))
end do
CALL afine%restore_local(array_f)
DEALLOCATE(array_c,array_f)
DEBUG_STACK_POP
end subroutine lag_base_pop
!---------------------------------------------------------------------------
! SUBROUTINE: lag_inject
!---------------------------------------------------------------------------
!> Inject a fine level Lagrange scalar field to the next coarsest level
!!
!! @note The global Lagrange level in decremented by one in this subroutine
!!
!! @param[in] afine Vector to inject
!! @param[in,out] acors Coarse vector from injection
!---------------------------------------------------------------------------
subroutine lag_inject(afine,acors)
class(oft_vector), intent(inout) :: afine
class(oft_vector), intent(inout) :: acors
integer(i4) :: i,j,k
logical :: gcheck
DEBUG_STACK_PUSH
gcheck=(oft_lagrange%order==1)
! Step down level up
call oft_lag_set_level(oft_lagrange_level-1)
! Cast fine field
call acors%set(0.d0)
if(oft_lagrange_level==oft_lagrange_blevel)then
  call lag_base_push(afine,acors)
  DEBUG_STACK_POP
  return
end if
CALL ML_oft_lagrange_ops(oft_lagrange_level+1)%interp%applyT(afine,acors)
DEBUG_STACK_POP
end subroutine lag_inject
!---------------------------------------------------------------------------
! SUBROUTINE: lag_vinject
!---------------------------------------------------------------------------
!> Inject a fine level Lagrange scalar field to the next coarsest level
!!
!! @note The global Lagrange level in decremented by one in this subroutine
!!
!! @param[in] afine Vector to inject
!! @param[in,out] acors Coarse vector from injection
!---------------------------------------------------------------------------
subroutine lag_vinject(afine,acors)
class(oft_vector), intent(inout) :: afine
class(oft_vector), intent(inout) :: acors
integer(i4) :: i,j,k
logical :: gcheck
DEBUG_STACK_PUSH
gcheck=(oft_lagrange%order==1)
! Step down level up
call oft_lag_set_level(oft_lagrange_level-1)
! Cast fine field
if(oft_lagrange_level==oft_lagrange_blevel)then
!  call lag_base_push(afine,acors)
  CALL oft_abort('MPI-Base transfer not supported','lag_vinterp',__FILE__)
  DEBUG_STACK_POP
  return
end if
CALL ML_oft_lagrange_ops(oft_lagrange_level+1)%vinterp%applyT(afine,acors)
DEBUG_STACK_POP
end subroutine lag_vinject
!---------------------------------------------------------------------------
! SUBROUTINE: lag_base_push
!---------------------------------------------------------------------------
!> Transfer a MPI level Lagrange scalar field to the base level
!!
!! @param[in] afine Vector to transfer
!! @param[in,out] acors Fine vector from transfer
!---------------------------------------------------------------------------
subroutine lag_base_push(afine,acors)
class(oft_vector), intent(inout) :: afine
class(oft_vector), intent(inout) :: acors
integer(i4), pointer :: lptmp(:)
integer(i4) :: i,j,ierr
real(r8), pointer, dimension(:) :: alias,array_c,array_f
CLASS(oft_afem_type), POINTER :: lag_fine => NULL()
DEBUG_STACK_PUSH
!---
lptmp=>mg_mesh%meshes(mg_mesh%nbase+1)%base%lp
CALL acors%get_local(array_c)
CALL afine%get_local(array_f)
lag_fine=>ML_oft_lagrange%levels(oft_lagrange_level+1)%fe
!---
allocate(alias(acors%n))
alias=0.d0
!$omp parallel do
do i=1,afine%n
  if(lag_fine%linkage%be(i))cycle
  alias(lptmp(i))=array_f(i)
end do
!$omp parallel do private(j)
do i=1,lag_fine%linkage%nbe
  j=lag_fine%linkage%lbe(i)
  if(.NOT.lag_fine%linkage%leo(i))cycle
  alias(lptmp(j))=array_f(j)
end do
!---Global reduction over all processors
array_c=oft_mpi_sum(alias,acors%n)
call acors%restore_local(array_c)
deallocate(alias,array_c,array_f)
DEBUG_STACK_POP
end subroutine lag_base_push
!---------------------------------------------------------------------------
! SUBROUTINE: lag_lop_eigs
!---------------------------------------------------------------------------
!> Compute eigenvalues and smoothing coefficients for the operator LAG::LOP
!---------------------------------------------------------------------------
SUBROUTINE lag_lop_eigs(minlev)
INTEGER(i4), INTENT(in) :: minlev
#ifdef HAVE_ARPACK
INTEGER(i4) :: i
REAL(r8) :: lam0
REAL(r8), ALLOCATABLE :: df(:)
CLASS(oft_vector), POINTER :: u
TYPE(oft_irlm_eigsolver) :: arsolver
CLASS(oft_matrix), POINTER :: md => NULL()
CLASS(oft_matrix), POINTER :: lop => NULL()
DEBUG_STACK_PUSH
!---------------------------------------------------------------------------
! Compute optimal smoother coefficients
!---------------------------------------------------------------------------
IF(oft_env%head_proc)WRITE(*,*)'Optimizing Jacobi damping for LAG::LOP'
ALLOCATE(df(oft_lagrange_nlevels))
df=0.d0
DO i=minlev,oft_lagrange_nlevels
  CALL oft_lag_set_level(i)
  !---Create fields
  CALL oft_lag_create(u)
  !---Get Ev range
  NULLIFY(lop)
  CALL oft_lag_getlop(lop,'zerob')
  CALL create_diagmatrix(md,lop%D)
  !---
  arsolver%A=>lop
  arsolver%M=>md
  arsolver%mode=2
  arsolver%tol=1.E-5_r8
  arsolver%bc=>lag_zerob
  CALL create_native_pre(arsolver%Minv, "jacobi")
  arsolver%Minv%A=>lop
  !---
  IF(oft_debug_print(1))WRITE(*,*)'  optimizing level = ',i
  CALL arsolver%max(u,lam0)
  df(i) = 1.8d0/lam0
  !---
  CALL u%delete
  CALL md%delete
  CALL lop%delete
END DO
!---Output
IF(oft_env%head_proc)THEN
  WRITE(*,'(A)',ADVANCE='NO')' df_lop = '
  DO i=1,oft_lagrange_nlevels-1
    WRITE(*,'(F5.3,A)',ADVANCE='NO')df(i),', '
  END DO
  WRITE(*,'(F5.3,A)')df(oft_lagrange_nlevels)
END IF
DEALLOCATE(df)
DEBUG_STACK_POP
#else
CALL oft_abort("Subroutine requires ARPACK", "lag_lop_eigs", __FILE__)
#endif
END SUBROUTINE lag_lop_eigs
!---------------------------------------------------------------------------
! SUBROUTINE: lag_getlop_pre
!---------------------------------------------------------------------------
!> Construct default MG preconditioner for LAG::LOP
!---------------------------------------------------------------------------
SUBROUTINE lag_getlop_pre(pre,mats,level,nlevels)
CLASS(oft_solver), POINTER, INTENT(out) :: pre
TYPE(oft_matrix_ptr), POINTER, INTENT(inout) :: mats(:)
INTEGER(i4), OPTIONAL, INTENT(in) :: level
INTEGER(i4), OPTIONAL, INTENT(in) :: nlevels
!--- ML structures for MG-preconditioner
TYPE(oft_matrix_ptr), POINTER :: ml_int(:)
INTEGER(i4), ALLOCATABLE :: levels(:)
REAL(r8), ALLOCATABLE :: df(:)
INTEGER(i4), ALLOCATABLE :: nu(:)
!---
INTEGER(i4) :: minlev,toplev,nl
INTEGER(i4) :: i,j,levin,ierr
LOGICAL :: create_mats
CHARACTER(LEN=2) :: lev_char
!---
TYPE(xml_node), POINTER :: pre_node
#ifdef HAVE_XML
integer(i4) :: nnodes
TYPE(xml_node), POINTER :: lag_node
#endif
DEBUG_STACK_PUSH
!---
minlev=1
toplev=oft_lagrange_level
levin=oft_lagrange_level
IF(PRESENT(level))toplev=level
IF(PRESENT(nlevels))minlev=toplev-nlevels+1
nl=toplev-minlev+1
!---
IF(minlev<1)CALL oft_abort('Minimum level is < 0','lag_getlop_pre',__FILE__)
IF(toplev>oft_lagrange_nlevels)CALL oft_abort('Maximum level is > lag_nlevels','lag_getlop_pre',__FILE__)
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
  CALL oft_lag_set_level(minlev+(i-1))
  levels(i)=minlev+(i-1)
  df(i)=df_lop(levels(i))
  nu(i)=nu_lop(levels(i))
  IF(df(i)<-1.d90)THEN
    WRITE(lev_char,'(I2.2)')levels(i)
    CALL oft_abort('Smoother values not set for level: '//lev_char,'lag_getlop_pre',__FILE__)
  END IF
  !---
  IF(create_mats)THEN
    NULLIFY(mats(i)%M)
    CALL oft_lag_getlop(mats(i)%M,'zerob')
  END IF
  IF(i>1)ml_int(i-1)%M=>oft_lagrange_ops%interp
END DO
CALL oft_lag_set_level(levin)
!---------------------------------------------------------------------------
! Search for XML-spec
!---------------------------------------------------------------------------
NULLIFY(pre_node)
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml))THEN
  CALL xml_get_element(oft_env%xml,"lagrange",lag_node,ierr)
  IF(ierr==0)CALL xml_get_element(lag_node,"jmlb",pre_node,ierr)
END IF
#endif
!---------------------------------------------------------------------------
! Setup preconditioner
!---------------------------------------------------------------------------
NULLIFY(pre)
CALL create_mlpre(pre,mats(1:nl),levels,nlevels=nl,create_vec=oft_lag_create,interp=lag_interp, &
     inject=lag_inject,bc=lag_zerob,stype=1,df=df,nu=nu,xml_root=pre_node)
!---------------------------------------------------------------------------
! Cleanup
!---------------------------------------------------------------------------
DEALLOCATE(ml_int,levels,df,nu)
DEBUG_STACK_POP
END SUBROUTINE lag_getlop_pre
end module oft_lag_operators
