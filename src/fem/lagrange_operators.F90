!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
MODULE oft_lag_operators
USE oft_base
USE oft_sort, ONLY: sort_array
USE oft_mesh_type, ONLY: oft_mesh, cell_is_curved
!---
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr, &
  oft_graph, oft_graph_ptr
USE oft_deriv_matrices, ONLY: oft_diagmatrix, create_diagmatrix
USE oft_solver_base, ONLY: oft_solver, oft_solver_bc
USE oft_la_utils, ONLY: create_matrix, combine_matrices
USE oft_solver_utils, ONLY: create_mlpre, create_native_pre
#ifdef HAVE_ARPACK
USE oft_arpack, ONLY: oft_irlm_eigsolver
#endif
!---
USE fem_base, ONLY: oft_afem_type, oft_fem_type, oft_bfem_type, fem_max_levels, oft_ml_fem_type, &
  oft_ml_fe_vecspace
USE fem_utils, ONLY: fem_interp
USE fem_composite, ONLY: oft_fem_comp_type, oft_ml_fem_comp_type
USE oft_lag_basis, ONLY: oft_lag_eval_all, oft_lag_geval_all, oft_lag_eval, &
  oft_lag_nodes, oft_blag_eval, &
  oft_blag_geval, oft_lag_npos, oft_scalar_fem, oft_scalar_bfem, &
  oft_3D_lagrange_cast
IMPLICIT NONE
#include "local.h"
!------------------------------------------------------------------------------
!> Interpolate a Lagrange field
!------------------------------------------------------------------------------
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
  !> Destroy interpolation object
  procedure :: delete => lag_rinterp_delete
end type oft_lag_rinterp
!------------------------------------------------------------------------------
!> Interpolate \f$ \nabla \f$ of a Lagrange field
!------------------------------------------------------------------------------
type, extends(oft_lag_rinterp) :: oft_lag_ginterp
contains
  !> Reconstruct field
  procedure :: interp => lag_ginterp_apply
end type oft_lag_ginterp
!------------------------------------------------------------------------------
!> Interpolate a Lagrange vector field
!------------------------------------------------------------------------------
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
  !> Destroy interpolation object
  procedure :: delete => lag_vrinterp_delete
end type oft_lag_vrinterp
!------------------------------------------------------------------------------
!> Interpolate \f$ \nabla \times \f$ of a Lagrange vector field
!------------------------------------------------------------------------------
type, extends(oft_lag_vrinterp) :: oft_lag_vcinterp
contains
  !> Reconstruct field
  procedure :: interp => lag_vcinterp
end type oft_lag_vcinterp
!------------------------------------------------------------------------------
!> Interpolate \f$ \nabla \f$ of a Lagrange vector field
!------------------------------------------------------------------------------
type, extends(oft_lag_vrinterp) :: oft_lag_vdinterp
contains
  !> Reconstruct field
  procedure :: interp => lag_vdinterp
end type oft_lag_vdinterp
!------------------------------------------------------------------------------
!> Zero Lagrange FE vector at all global boundary nodes
!------------------------------------------------------------------------------
type, extends(oft_solver_bc) :: oft_lag_zerob
  class(oft_ml_fem_type), pointer :: ML_lag_rep => NULL() !< FE representation
contains
  !> Zero field on boundary nodes
  procedure :: apply => zerob_apply
  !> Destroy BC object
  procedure :: delete => zerob_delete
end type oft_lag_zerob
!------------------------------------------------------------------------------
!> Zero Lagrange FE vector at all global grounding node(s)
!------------------------------------------------------------------------------
type, extends(oft_lag_zerob) :: oft_lag_zerogrnd
contains
  !> Zero field at grounding node(s)
  procedure :: apply => zerogrnd_apply
end type oft_lag_zerogrnd
!------------------------------------------------------------------------------
!> Zero all components of vector Lagrange FE vector at all global boundary nodes
!------------------------------------------------------------------------------
type, extends(oft_solver_bc) :: oft_vlag_zerob
  class(oft_ml_fem_comp_type), pointer :: ML_vlag_rep => NULL() !< FE representation
contains
  !> Zero all components of vector on boundary nodes
  procedure :: apply => vzerob_apply
  !> Destroy BC object
  procedure :: delete => vzerob_delete
end type oft_vlag_zerob
!------------------------------------------------------------------------------
!> Zero normal component of vector Lagrange FE vector at all global boundary nodes
!------------------------------------------------------------------------------
type, extends(oft_vlag_zerob) :: oft_vlag_zeron
contains
  !> Zero normal component of vector on boundary nodes
  procedure :: apply => vzeron_apply
end type oft_vlag_zeron
!------------------------------------------------------------------------------
!> Zero tangential component of vector Lagrange FE vector at all global boundary nodes
!------------------------------------------------------------------------------
type, extends(oft_vlag_zerob) :: oft_vlag_zerot
contains
  !> Zero tangential component of vector on boundary nodes
  procedure :: apply => vzerot_apply
end type oft_vlag_zerot
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
!------------------------------------------------------------------------------
!> Read-in options for the basic Lagrange ML preconditioners
!------------------------------------------------------------------------------
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
    CALL oft_warn("No Lagrange MG smoother settings found:")
    CALL oft_warn("  Using default values, which may result in convergence failure.")
  END IF
  nu_lop=2
  df_lop=.2d0
  nu_pdop=2
  df_pdop=.2d0
END IF
DEBUG_STACK_POP
end subroutine lag_mloptions
!------------------------------------------------------------------------------
!> Setup interpolator for Lagrange scalar fields
!!
!! Fetches local representation used for interpolation from vector object
!------------------------------------------------------------------------------
subroutine lag_rinterp_setup(self,lag_rep)
class(oft_lag_rinterp), intent(inout) :: self
class(oft_afem_type), target, intent(inout) :: lag_rep
IF(ASSOCIATED(self%parent))THEN
  SELECT TYPE(this=>self%parent)
    CLASS IS(oft_lag_rinterp)
      self%vals=>this%vals
      self%lag_rep=>this%lag_rep
      self%mesh=>this%mesh
      self%cache_cell=>this%cache_cell
      self%cache_vals=>this%cache_vals
    CLASS DEFAULT
      CALL oft_abort('Parent interpolator must be of same class','oft_lag_rinterp',__FILE__)
  END SELECT
ELSE
  !---Get local slice
  CALL self%u%get_local(self%vals)
  !
  IF(.NOT.oft_3D_lagrange_cast(self%lag_rep,lag_rep))CALL oft_abort("Incorrect FE type","lag_rinterp_setup",__FILE__)
  self%mesh=>self%lag_rep%mesh
  IF(.NOT.ASSOCIATED(self%cache_cell))THEN
    ALLOCATE(self%cache_cell(0:oft_env%nthreads-1))
    ALLOCATE(self%cache_vals(self%lag_rep%nce,0:oft_env%nthreads-1))
  END IF
  self%cache_cell=-1
  self%cache_vals=0.d0
END IF
end subroutine lag_rinterp_setup
!------------------------------------------------------------------------------
!> Destroy temporary internal storage and nullify references
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Reconstruct a Lagrange scalar field
!------------------------------------------------------------------------------
subroutine lag_rinterp(self,cell,f,gop,val)
class(oft_lag_rinterp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [4]
real(r8), intent(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
real(r8), intent(out) :: val(:) !< Reconstructed field at f [1]
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
!------------------------------------------------------------------------------
!> Reconstruct the gradient of a Lagrange scalar field
!------------------------------------------------------------------------------
subroutine lag_ginterp_apply(self,cell,f,gop,val)
class(oft_lag_ginterp), intent(inout) :: self !< 
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [4]
real(r8), intent(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
real(r8), intent(out) :: val(:) !< Reconstructed gradient at f [3]
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
!------------------------------------------------------------------------------
!> Setup interpolator for Lagrange vector fields
!!
!! Fetches local representation used for interpolation from vector object
!------------------------------------------------------------------------------
subroutine lag_vrinterp_setup(self,lag_rep)
class(oft_lag_vrinterp), intent(inout) :: self
class(oft_afem_type), target, intent(inout) :: lag_rep
real(r8), pointer, dimension(:) :: vtmp
IF(ASSOCIATED(self%parent))THEN
  SELECT TYPE(this=>self%parent)
    CLASS IS(oft_lag_vrinterp)
      self%vals=>this%vals
      self%lag_rep=>this%lag_rep
      self%mesh=>self%mesh
      self%cache_cell=>this%cache_cell
      self%cache_vals=>this%cache_vals
    CLASS DEFAULT
      CALL oft_abort('Parent interpolator must be of same class','lag_vrinterp_setup',__FILE__)
  END SELECT
ELSE
  !---Get local slice
  IF(.NOT.ASSOCIATED(self%vals))ALLOCATE(self%vals(3,lag_rep%ne))
  vtmp=>self%vals(1,:)
  CALL self%u%get_local(vtmp,1)
  vtmp=>self%vals(2,:)
  CALL self%u%get_local(vtmp,2)
  vtmp=>self%vals(3,:)
  CALL self%u%get_local(vtmp,3)
  IF(.NOT.oft_3D_lagrange_cast(self%lag_rep,lag_rep))CALL oft_abort("Incorrect FE type","lag_vrinterp_setup",__FILE__)
  self%mesh=>self%lag_rep%mesh
  IF(.NOT.ASSOCIATED(self%cache_cell))THEN
    ALLOCATE(self%cache_cell(0:oft_env%nthreads-1))
    ALLOCATE(self%cache_vals(3,self%lag_rep%nce,0:oft_env%nthreads-1))
  END IF
  self%cache_cell=-1
  self%cache_vals=0.d0
END IF
end subroutine lag_vrinterp_setup
!------------------------------------------------------------------------------
!> Destroy temporary internal storage and nullify references
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Reconstruct a Lagrange vector field
!------------------------------------------------------------------------------
subroutine lag_vrinterp(self,cell,f,gop,val)
class(oft_lag_vrinterp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [4]
real(r8), intent(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
real(r8), intent(out) :: val(:) !< Reconstructed field at f [3]
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
!------------------------------------------------------------------------------
!> Reconstruct the curl of a Lagrange vector field
!------------------------------------------------------------------------------
subroutine lag_vcinterp(self,cell,f,gop,val)
class(oft_lag_vcinterp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [4]
real(r8), intent(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
real(r8), intent(out) :: val(:) !< Reconstructed field at f [3]
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
!------------------------------------------------------------------------------
!> Reconstruct \f$ \nabla v \f$ for a Lagrange vector field.
!!
!! The tensor is packed using reshape, to retrieve use
!!\code
!! dv = RESHAPE(val,(/3,3/))
!!\endcode
!------------------------------------------------------------------------------
subroutine lag_vdinterp(self,cell,f,gop,val)
class(oft_lag_vdinterp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [4]
real(r8), intent(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
real(r8), intent(out) :: val(:) !< Reconstructed field at f [3]
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
!------------------------------------------------------------------------------
!> Zero a Lagrange scalar field at all boundary nodes
!!
!! @param[in,out] a Field to be zeroed
!------------------------------------------------------------------------------
subroutine zerob_apply(self,a)
class(oft_lag_zerob), intent(inout) :: self
class(oft_vector), intent(inout) :: a
integer(i4) :: i,j
real(r8), pointer, dimension(:) :: vloc
DEBUG_STACK_PUSH
!---Cast to vector type
NULLIFY(vloc)
call a%get_local(vloc)
!---Zero boundary values
!$omp parallel do private(j)
do i=1,self%ML_lag_rep%current_level%nbe
  j=self%ML_lag_rep%current_level%lbe(i)
  if(self%ML_lag_rep%current_level%global%gbe(j))vloc(j)=0.d0
end do
!---
call a%restore_local(vloc)
deallocate(vloc)
DEBUG_STACK_POP
end subroutine zerob_apply
!------------------------------------------------------------------------------
!> Destroy temporary internal storage and nullify references
!------------------------------------------------------------------------------
subroutine zerob_delete(self)
class(oft_lag_zerob), intent(inout) :: self
NULLIFY(self%ML_lag_rep)
end subroutine zerob_delete
!------------------------------------------------------------------------------
!> Zero a Lagrange scalar field at the global grounding node
!!
!! @note The possition of this node is defined by the mesh pointer igrnd in
!! mesh
!------------------------------------------------------------------------------
subroutine zerogrnd_apply(self,a)
class(oft_lag_zerogrnd), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to be zeroed
real(r8), pointer, dimension(:) :: aloc
integer(i4) :: i,j
CLASS(oft_scalar_fem), POINTER ::lag_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_lagrange_cast(lag_rep,self%ML_lag_rep%current_level))CALL oft_abort("Incorrect FE type","zerogrnd_apply",__FILE__)
!---Zero boundary values
NULLIFY(aloc)
CALL a%get_local(aloc)
!---
if(lag_rep%mesh%igrnd(1)>0)aloc(lag_rep%mesh%igrnd(1))=0.d0
if(lag_rep%mesh%igrnd(2)>0)aloc(lag_rep%mesh%igrnd(2))=0.d0
CALL a%restore_local(aloc)
DEALLOCATE(aloc)
DEBUG_STACK_POP
end subroutine zerogrnd_apply
!------------------------------------------------------------------------------
!> Zero a surface Lagrange vector field at all edge nodes
!------------------------------------------------------------------------------
subroutine vzerob_apply(self,a)
class(oft_vlag_zerob), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to be zeroed
integer(i4) :: i,j
real(r8), pointer :: vloc(:,:),vtmp(:)
CLASS(oft_afem_type), POINTER :: lag_rep
DEBUG_STACK_PUSH
!---Cast to vector type
ALLOCATE(vloc(3,a%n))
DO i=1,3
  vtmp=>vloc(i,:)
  call a%get_local(vtmp,i)
END DO
!---Zero boundary values
lag_rep=>self%ML_vlag_rep%current_level%fields(1)%fe
!$omp parallel do private(j)
do i=1,lag_rep%nbe
  j=lag_rep%lbe(i)
  if(lag_rep%global%gbe(j))vloc(:,j)=0.d0
end do
!---
DO i=1,3
  vtmp=>vloc(i,:)
  call a%restore_local(vtmp,i,wait=(i/=3))
END DO
DEALLOCATE(vloc)
DEBUG_STACK_POP
end subroutine vzerob_apply
!------------------------------------------------------------------------------
!> Destroy temporary internal storage and nullify references
!------------------------------------------------------------------------------
subroutine vzerob_delete(self)
class(oft_vlag_zerob), intent(inout) :: self
NULLIFY(self%ML_vlag_rep)
end subroutine vzerob_delete
!------------------------------------------------------------------------------
!> Zero normal component of a Lagrange vector field at boundary nodes
!------------------------------------------------------------------------------
subroutine vzeron_apply(self,a)
class(oft_vlag_zeron), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to be zeroed
real(r8) :: vec(3),nn(3,3)
real(r8), pointer, dimension(:) :: x,y,z
integer(i4) :: i,j
CLASS(oft_scalar_fem), POINTER :: lag_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_lagrange_cast(lag_rep,self%ML_vlag_rep%current_level%fields(1)%fe))CALL oft_abort("Incorrect FE type","vzeron_apply",__FILE__)
!---Cast to vector type
NULLIFY(x,y,z)
CALL a%get_local(x,1)
CALL a%get_local(y,2)
CALL a%get_local(z,3)
!---Zero boundary values
!$omp parallel do private(j,vec,nn)
do i=1,lag_rep%nbe
  j=lag_rep%lbe(i)
  if(lag_rep%global%gbe(j))then
    CALL lag_vbc_tensor(lag_rep,j,1,nn)
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
end subroutine vzeron_apply
!------------------------------------------------------------------------------
!> Zero tangential component of a Lagrange vector field at boundary nodes
!------------------------------------------------------------------------------
subroutine vzerot_apply(self,a)
class(oft_vlag_zerot), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to be zeroed
real(r8) :: vec(3),nn(3,3)
real(r8), pointer, dimension(:) :: x,y,z
integer(i4) :: i,j
CLASS(oft_scalar_fem), POINTER :: lag_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_lagrange_cast(lag_rep,self%ML_vlag_rep%current_level%fields(1)%fe))CALL oft_abort("Incorrect FE type","vzerot_apply",__FILE__)
!---Cast to vector type
NULLIFY(x,y,z)
CALL a%get_local(x,1)
CALL a%get_local(y,2)
CALL a%get_local(z,3)
!---Zero boundary values
!$omp parallel do private(j,vec,nn)
do i=1,lag_rep%nbe
  j=lag_rep%lbe(i)
  if(lag_rep%global%gbe(j))then
    CALL lag_vbc_tensor(lag_rep,j,2,nn)
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
end subroutine vzerot_apply
!------------------------------------------------------------------------------
!> Get boundary condition tensor for desired node in a Lagrange vector field
!------------------------------------------------------------------------------
subroutine lag_vbc_tensor(lag_rep,j_lag,bc_type,nn)
class(oft_scalar_fem), target, intent(inout) :: lag_rep
integer(i4), intent(in) :: j_lag !< Local index of Lagrange node
integer(i4), intent(in) :: bc_type !< Desired BC (1 -> norm, 2-> tang)
real(r8), intent(out) :: nn(3,3) !< BC tensor (3,3)
integer(i4) :: i,j
DEBUG_STACK_PUSH
IF(lag_rep%global%gbe(j_lag))THEN
  SELECT CASE(bc_type)
    CASE(1) ! Zero normal component
      IF(lag_rep%bc(j_lag)==3)THEN ! On face
        nn(:,1) = (/1.d0,0.d0,0.d0/) &
          - lag_rep%sn(:,j_lag)*lag_rep%sn(1,j_lag)
        nn(:,2) = (/0.d0,1.d0,0.d0/) &
          - lag_rep%sn(:,j_lag)*lag_rep%sn(2,j_lag)
        nn(:,3) = (/0.d0,0.d0,1.d0/) &
          - lag_rep%sn(:,j_lag)*lag_rep%sn(3,j_lag)
      ELSE IF(lag_rep%bc(j_lag)==2)THEN ! On edge
        nn(:,1) = lag_rep%sn(:,j_lag)*lag_rep%sn(1,j_lag)
        nn(:,2) = lag_rep%sn(:,j_lag)*lag_rep%sn(2,j_lag)
        nn(:,3) = lag_rep%sn(:,j_lag)*lag_rep%sn(3,j_lag)
      ELSE ! At vertex
        nn=0.d0
      END IF
    CASE(2) ! Zero tangential components
      IF(lag_rep%bc(j_lag)==3)THEN ! On face
        nn(:,1) = lag_rep%sn(:,j_lag)*lag_rep%sn(1,j_lag)
        nn(:,2) = lag_rep%sn(:,j_lag)*lag_rep%sn(2,j_lag)
        nn(:,3) = lag_rep%sn(:,j_lag)*lag_rep%sn(3,j_lag)
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
!------------------------------------------------------------------------------
!> Get diagonal entries for a given boundary condition and desired node in a
!! Lagrange vector field
!------------------------------------------------------------------------------
subroutine lag_vbc_diag(lag_rep,j_lag,bc_type,nn)
class(oft_scalar_fem), target, intent(inout) :: lag_rep
integer(i4), intent(in) :: j_lag !< Local index of Lagrange node
integer(i4), intent(in) :: bc_type !< Desired BC (1 -> norm, 2-> tang)
real(r8), intent(out) :: nn(3,3) !< Diagonal entries (3,3)
integer(i4) :: i,j
DEBUG_STACK_PUSH
SELECT CASE(bc_type)
  CASE(1) ! Get normal projections
    !---Set diagonal entries for boundary terms
    IF(lag_rep%bc(j_lag)==3)THEN
      nn(:,1) = lag_rep%sn(:,j_lag)*lag_rep%sn(1,j_lag)
      nn(:,2) = lag_rep%sn(:,j_lag)*lag_rep%sn(2,j_lag)
      nn(:,3) = lag_rep%sn(:,j_lag)*lag_rep%sn(3,j_lag)
    ELSE IF(lag_rep%bc(j_lag)==2)THEN
      nn(:,1) = (/1.d0,0.d0,0.d0/) &
        - lag_rep%sn(:,j_lag)*lag_rep%sn(1,j_lag)
      nn(:,2) = (/0.d0,1.d0,0.d0/) &
        - lag_rep%sn(:,j_lag)*lag_rep%sn(2,j_lag)
      nn(:,3) = (/0.d0,0.d0,1.d0/) &
        - lag_rep%sn(:,j_lag)*lag_rep%sn(3,j_lag)
    ELSE
      nn(:,1)=(/1.d0,0.d0,0.d0/)
      nn(:,2)=(/0.d0,1.d0,0.d0/)
      nn(:,3)=(/0.d0,0.d0,1.d0/)
    END IF
  CASE(2) ! Get tangential projections
    IF(lag_rep%bc(j_lag)==3)THEN ! On face
      nn(:,1) = (/1.d0,0.d0,0.d0/) &
        - lag_rep%sn(:,j_lag)*lag_rep%sn(1,j_lag)
      nn(:,2) = (/0.d0,1.d0,0.d0/) &
        - lag_rep%sn(:,j_lag)*lag_rep%sn(2,j_lag)
      nn(:,3) = (/0.d0,0.d0,1.d0/) &
        - lag_rep%sn(:,j_lag)*lag_rep%sn(3,j_lag)
    ELSE
      nn(:,1)=(/1.d0,0.d0,0.d0/)
      nn(:,2)=(/0.d0,1.d0,0.d0/)
      nn(:,3)=(/0.d0,0.d0,1.d0/)
    END IF
END SELECT
DEBUG_STACK_POP
end subroutine lag_vbc_diag
!------------------------------------------------------------------------------
!> Construct mass matrix for Lagrange scalar representation
!!
!! Supported boundary conditions
!! - `'none'` Full matrix
!! - `'zerob'` Dirichlet for all boundary DOF
!------------------------------------------------------------------------------
subroutine oft_lag_getmop(fe_rep,mat,bc)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix object
character(LEN=*), intent(in) :: bc !< Boundary condition
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time
real(r8), allocatable :: rop(:),mop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_lag_vec
type(oft_timer) :: mytimer
CLASS(oft_scalar_fem), POINTER :: lag_rep
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing LAG::MOP'
  CALL mytimer%tick()
END IF
IF(.NOT.oft_3D_lagrange_cast(lag_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_lag_getmop",__FILE__)
!------------------------------------------------------------------------------
! Allocate matrix
!------------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL lag_rep%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,rop,det,mop,curved,goptmp,m,vol,jc,jr)
allocate(j(lag_rep%nce)) ! Local DOF and matrix indices
allocate(rop(lag_rep%nce)) ! Reconstructed gradient operator
allocate(mop(lag_rep%nce,lag_rep%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,lag_rep%mesh%nc
  curved=cell_is_curved(lag_rep%mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  mop=0.d0
  do m=1,lag_rep%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call lag_rep%mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
    det=vol*lag_rep%quad%wts(m)
    call oft_lag_eval_all(lag_rep,i,lag_rep%quad%pts(:,m),rop)
    !---Compute local matrix contributions
    do jr=1,lag_rep%nce
      do jc=1,lag_rep%nce
        mop(jr,jc) = mop(jr,jc) + rop(jr)*rop(jc)*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call lag_rep%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,lag_rep%nce
        IF(lag_rep%global%gbe(j(jr)))mop(jr,:)=0.d0
      END DO
    CASE("grnd")
      IF(ANY(lag_rep%mesh%igrnd>0))THEN
        DO jr=1,lag_rep%nce
          IF(ANY(j(jr)==lag_rep%mesh%igrnd))mop(jr,:)=0.d0
        END DO
      END IF
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,mop,lag_rep%nce,lag_rep%nce)
  ! !$omp end critical
end do
deallocate(j,rop,mop)
!$omp end parallel
ALLOCATE(mop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
mop(1,1)=1.d0
SELECT CASE(TRIM(bc))
  CASE("zerob")
    DO i=1,lag_rep%nbe
      jr=lag_rep%lbe(i)
      IF(lag_rep%global%gbe(jr).AND.lag_rep%linkage%leo(i))THEN
        j=jr
        call mat%add_values(j,j,mop,1,1)
      END IF
    END DO
  CASE("grnd")
    IF(ANY(lag_rep%mesh%igrnd>0))THEN
      DO i=1,lag_rep%nbe
        jr=lag_rep%lbe(i)
        IF(lag_rep%linkage%leo(i).AND.ANY(jr==lag_rep%mesh%igrnd))THEN
          j=jr
          call mat%add_values(j,j,mop,1,1)
        END IF
      END DO
    END IF
END SELECT
DEALLOCATE(j,mop)
CALL lag_rep%vec_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_lag_getmop
!------------------------------------------------------------------------------
!> Construct laplacian matrix for Lagrange scalar representation
!!
!! Supported boundary conditions
!! - `'none'` Full matrix
!! - `'zerob'` Dirichlet for all boundary DOF
!! - `'grnd'`  Dirichlet for only groundin point
!------------------------------------------------------------------------------
subroutine oft_lag_getlop(fe_rep,mat,bc)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix object
character(LEN=*), intent(in) :: bc !< Boundary condition
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time
real(r8), allocatable :: gop(:,:),lop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_lag_vec
type(oft_timer) :: mytimer
CLASS(oft_scalar_fem), POINTER :: lag_rep
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing LAG::LOP'
  CALL mytimer%tick()
END IF
IF(.NOT.oft_3D_lagrange_cast(lag_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_lag_getlop",__FILE__)
!------------------------------------------------------------------------------
! Allocate matrix
!------------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL lag_rep%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!------------------------------------------------------------------------------
! Operator integration
!------------------------------------------------------------------------------
!$omp parallel private(j,gop,det,lop,curved,goptmp,m,vol,jc,jr)
allocate(j(lag_rep%nce)) ! Local DOF and matrix indices
allocate(gop(3,lag_rep%nce)) ! Reconstructed gradient operator
allocate(lop(lag_rep%nce,lag_rep%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,lag_rep%mesh%nc
  curved=cell_is_curved(lag_rep%mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  lop=0.d0
  do m=1,lag_rep%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call lag_rep%mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
    det=vol*lag_rep%quad%wts(m)
    call oft_lag_geval_all(lag_rep,i,lag_rep%quad%pts(:,m),gop,goptmp)
    !---Compute local matrix contributions
    do jr=1,lag_rep%nce
      do jc=1,lag_rep%nce
        lop(jr,jc) = lop(jr,jc) + dot_product(gop(:,jr),gop(:,jc))*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call lag_rep%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,lag_rep%nce
        IF(lag_rep%global%gbe(j(jr)))lop(jr,:)=0.d0
      END DO
    CASE("grnd")
      IF(ANY(lag_rep%mesh%igrnd>0))THEN
        DO jr=1,lag_rep%nce
          IF(ANY(j(jr)==lag_rep%mesh%igrnd))lop(jr,:)=0.d0
        END DO
      END IF
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,lop,lag_rep%nce,lag_rep%nce)
  ! !$omp end critical
end do
deallocate(j,gop,lop)
!$omp end parallel
ALLOCATE(lop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
lop(1,1)=1.d0
SELECT CASE(TRIM(bc))
  CASE("zerob")
    DO i=1,lag_rep%nbe
      jr=lag_rep%lbe(i)
      IF(lag_rep%global%gbe(jr).AND.lag_rep%linkage%leo(i))THEN
        j=jr
        call mat%add_values(j,j,lop,1,1)
      END IF
    END DO
  CASE("grnd")
    IF(ANY(lag_rep%mesh%igrnd>0))THEN
      DO i=1,lag_rep%nbe
        jr=lag_rep%lbe(i)
        IF(lag_rep%linkage%leo(i).AND.ANY(jr==lag_rep%mesh%igrnd))THEN
          j=jr
          call mat%add_values(j,j,lop,1,1)
        END IF
      END DO
    END IF
END SELECT
DEALLOCATE(j,lop)
CALL lag_rep%vec_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_lag_getlop
!------------------------------------------------------------------------------
!> Construct parallel diffusion matrix for Lagrange scalar representation
!!
!! Supported boundary conditions
!! - `'none'` Full matrix
!! - `'zerob'` Dirichlet for all boundary DOF
!! - `'grnd'`  Dirichlet for only groundin point
!------------------------------------------------------------------------------
subroutine oft_lag_getpdop(fe_rep,mat,field,bc,perp,be_flag)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix object
class(fem_interp), intent(inout) :: field !< Vector field defining \f$ \hat{b} \f$
character(LEN=*), intent(in) :: bc !< Boundary condition
real(r8), optional, intent(in) :: perp !< Value of perpendicular conductivity (optional)
logical, optional, intent(in) :: be_flag(:) !< Flag for dirichlet nodes if different from boundary [ne] (optional)
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),B_nodal(3),par_diff,perp_diff,elapsed_time
real(r8), allocatable :: gop(:,:),pdop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_lag_vec
CLASS(oft_scalar_fem), POINTER :: lag_rep
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing LAG::PDOP'
  CALL mytimer%tick()
END IF
IF(.NOT.oft_3D_lagrange_cast(lag_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_lag_getpdop",__FILE__)
!---
perp_diff=0.d0
if(present(perp))perp_diff=perp
par_diff=1.d0-perp_diff
!---
IF(TRIM(bc)=="list".AND.(.NOT.PRESENT(be_flag)))CALL &
oft_abort('Boundary flag array not provided.','oft_lag_getpdop',__FILE__)
!------------------------------------------------------------------------------
! Allocate matrix
!------------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL lag_rep%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,gop,det,pdop,B_nodal,curved,goptmp,m,vol,jc,jr)
allocate(j(lag_rep%nce)) ! Local DOF and matrix indices
allocate(gop(3,lag_rep%nce)) ! Reconstructed gradient operator
allocate(pdop(lag_rep%nce,lag_rep%nce)) ! Local laplacian matrix
!$omp do
do i=1,lag_rep%mesh%nc
  !---Get local to global DOF mapping
  call lag_rep%ncdofs(i,j)
  curved=cell_is_curved(lag_rep%mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  pdop=0.d0
  do m=1,lag_rep%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call lag_rep%mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
    det=vol*lag_rep%quad%wts(m)
    call oft_lag_geval_all(lag_rep,i,lag_rep%quad%pts(:,m),gop,goptmp)
    !---Compute bhat
    call field%interp(i,lag_rep%quad%pts(:,m),goptmp,B_nodal)
    B_nodal=B_nodal/SQRT(SUM(B_nodal**2))
    !---Compute local matrix contributions
    do jr=1,lag_rep%nce
      do jc=1,lag_rep%nce
        pdop(jr,jc) = pdop(jr,jc) + &
        par_diff*DOT_PRODUCT(gop(:,jr),B_nodal)*DOT_PRODUCT(gop(:,jc),B_nodal)*det + &
        perp_diff*DOT_PRODUCT(gop(:,jr),gop(:,jc))*det
      end do
    end do
  end do
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,lag_rep%nce
        IF(lag_rep%global%gbe(j(jr)))pdop(jr,:)=0.d0
      END DO
    CASE("list")
      DO jr=1,lag_rep%nce
        IF(be_flag(j(jr)))pdop(jr,:)=0.d0
      END DO
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,pdop,lag_rep%nce,lag_rep%nce)
  ! !$omp end critical
end do
deallocate(j,gop,pdop)
!$omp end parallel
ALLOCATE(pdop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
SELECT CASE(TRIM(bc))
  CASE("zerob")
    pdop(1,1)=1.d0
    DO i=1,lag_rep%nbe
      jr=lag_rep%lbe(i)
      IF(.NOT.lag_rep%global%gbe(jr))CYCLE
      IF(.NOT.lag_rep%linkage%leo(i))CYCLE
      j=jr
      call mat%add_values(j,j,pdop,1,1)
    END DO
  CASE("list")
    pdop(1,1)=1.d0
    DO i=1,lag_rep%ne
      IF(lag_rep%be(i))CYCLE
      IF(.NOT.be_flag(i))CYCLE
      j=i
      call mat%add_values(j,j,pdop,1,1)
    END DO
    DO i=1,lag_rep%nbe
      IF(.NOT.lag_rep%linkage%leo(i))CYCLE
      IF(.NOT.be_flag(lag_rep%lbe(i)))CYCLE
      j=lag_rep%lbe(i)
      call mat%add_values(j,j,pdop,1,1)
    END DO
END SELECT
DEALLOCATE(j,pdop)
CALL lag_rep%vec_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_lag_getpdop
!------------------------------------------------------------------------------
!> Project a scalar field onto a lagrange basis
!!
!! @note This subroutine only performs the integration of the field with
!! test functions for a Lagrange basis. To retrieve the correct projection the
!! result must be multiplied by the inverse of LAG::MOP
!------------------------------------------------------------------------------
subroutine oft_lag_project(fe_rep,field,x)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(fem_interp), intent(inout) :: field !< Scalar field for projection
class(oft_vector), intent(inout) :: x !< Field projected onto Lagrange basis
!---
integer(i4) :: i,jc,m
integer(i4), allocatable :: j(:)
real(r8) :: bcc(1),vol,det,goptmp(3,4)
real(r8), pointer :: xloc(:)
real(r8), allocatable :: rop(:)
logical :: curved
CLASS(oft_scalar_fem), POINTER :: lag_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_lagrange_cast(lag_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_lag_project",__FILE__)
!---Initialize vectors to zero
NULLIFY(xloc)
call x%set(0.d0)
call x%get_local(xloc)
!---Integerate over the volume
!$omp parallel private(j,rop,curved,m,goptmp,vol,det,bcc,jc)
!---Allocate cell DOF arrays
allocate(j(lag_rep%nce),rop(lag_rep%nce))
!$omp do schedule(guided)
do i=1,lag_rep%mesh%nc ! Loop over cells
  call lag_rep%ncdofs(i,j) ! Get DOFs
  curved=cell_is_curved(lag_rep%mesh,i) ! Straight cell test
  do m=1,lag_rep%quad%np
    if(curved.OR.m==1)call lag_rep%mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
    det=vol*lag_rep%quad%wts(m)
    call field%interp(i,lag_rep%quad%pts(:,m),goptmp,bcc)
    call oft_lag_eval_all(lag_rep,i,lag_rep%quad%pts(:,m),rop)
    do jc=1,lag_rep%nce
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
!------------------------------------------------------------------------------
!> Project the divergence of a scalar field onto a lagrange basis
!!
!! @note This subroutine only performs the integration of the field with
!! test functions for a Lagrange basis. To retrieve the correct projection the
!! result must be multiplied by the inverse of LAG::MOP
!------------------------------------------------------------------------------
subroutine oft_lag_project_div(fe_rep,field,x)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(fem_interp), intent(inout) :: field !< Scalar field for projection
class(oft_vector), intent(inout) :: x !< Field projected onto Lagrange basis
!---
real(r8) :: bcc(3),det,vol,goptmp(3,4)
real(r8), pointer :: xloc(:)
real(r8), allocatable :: gop(:,:)
integer(i4) :: i,jc,m
integer(i4), allocatable :: j(:)
logical :: curved
CLASS(oft_scalar_fem), POINTER :: lag_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_lagrange_cast(lag_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_lag_project_div",__FILE__)
!---Initialize vectors to zero
NULLIFY(xloc)
call x%set(0.d0)
call x%get_local(xloc)
!---Integerate over the volume
!$omp parallel private(j,gop,curved,m,goptmp,vol,det,bcc,jc)
!---Allocate cell DOF arrays
allocate(j(lag_rep%nce),gop(3,lag_rep%nce))
!$omp do schedule(guided)
do i=1,lag_rep%mesh%nc ! Loop over cells
  call lag_rep%ncdofs(i,j) ! Get DOFs
  curved=cell_is_curved(lag_rep%mesh,i) ! Straight cell test
  do m=1,lag_rep%quad%np
    if(curved.OR.m==1)call lag_rep%mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
    det=vol*lag_rep%quad%wts(m)
    call field%interp(i,lag_rep%quad%pts(:,m),goptmp,bcc)
    call oft_lag_geval_all(lag_rep,i,lag_rep%quad%pts(:,m),gop,goptmp)
    do jc=1,lag_rep%nce
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
!------------------------------------------------------------------------------
!> Construct mass matrix for Lagrange vector representation
!!
!! Supported boundary conditions
!! - `'none'` Full matrix
!! - `'all'` Dirichlet for all components at boundary
!! - `'norm'` Dirichlet for normal component at boundary
!! - `'tang'` Dirichlet for tangential component at boundary
!------------------------------------------------------------------------------
subroutine oft_lag_vgetmop(vlag_rep,mat,bc)
class(oft_fem_comp_type), target, intent(inout) :: vlag_rep
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix object
character(LEN=*), intent(in) :: bc !< Boundary condition
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
class(oft_scalar_fem), pointer :: lag_rep
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing LAG_V::MOP'
  CALL mytimer%tick()
END IF
IF(.NOT.oft_3D_lagrange_cast(lag_rep,vlag_rep%fields(1)%fe))CALL oft_abort("Incorrect FE type","oft_lag_vgetmop",__FILE__)
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
!------------------------------------------------------------------------------
! Allocate matrix
!------------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL vlag_rep%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,rop,det,mtmp,nn,mloc,curved,goptmp,m,u,vol,jc,jr,lcache)
allocate(j(lag_rep%nce)) ! Local DOF and matrix indices
allocate(rop(lag_rep%nce)) ! Reconstructed gradient operator
allocate(lcache(lag_rep%nce,lag_rep%nce))
DO jr=1,3
  DO jc=1,3
    allocate(mtmp(jr,jc)%m(lag_rep%nce,lag_rep%nce)) ! Local laplacian matrix
  END DO
END DO
!$omp do schedule(guided)
DO i=1,lag_rep%mesh%nc
  curved=cell_is_curved(lag_rep%mesh,i) ! Straight cell test
  !---Get local to global DOF mapping
  CALL lag_rep%ncdofs(i,j)
  !---Compute local matrix contributions
  DO jp=1,3
    DO jn=1,3
      IF(vbc_type<=0.AND.jp/=jn)CYCLE
      mtmp(jp,jn)%m = 0.d0
    END DO
  END DO
  !---Get local reconstructed operators
  DO m=1,lag_rep%quad%np ! Loop over quadrature points
    IF(curved.OR.m==1)CALL lag_rep%mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
    det=vol*lag_rep%quad%wts(m)
    call oft_lag_eval_all(lag_rep,i,lag_rep%quad%pts(:,m),rop)
    !---Compute local matrix contributions
    DO jc=1,lag_rep%nce
      u=rop(jc)*det
      DO jr=1,lag_rep%nce
        DO jp=1,3
          mtmp(jp,jp)%m(jr,jc) = mtmp(jp,jp)%m(jr,jc) + rop(jr)*u
        END DO
      END DO
    END DO
  END DO
  !---Compute BC projection matrices
  DO jr=1,lag_rep%nce
    IF(lag_rep%global%gbe(j(jr)))THEN
      IF(vbc_type==0)THEN
        DO jp=1,3
          mtmp(jp,jp)%m(jr,:)=0.d0
        END DO
      ELSE IF(vbc_type>0)THEN
        DO jc=1,lag_rep%nce
          mloc=0.d0
          u=mtmp(1,1)%m(jr,jc)
          DO jp=1,3
            mloc(jp,jp)=u
          END DO
          CALL lag_vbc_tensor(lag_rep,j(jr),vbc_type,nn)
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
      CALL mat%atomic_add_values(j,j,mtmp(jp,jn)%m,lag_rep%nce,lag_rep%nce,jp,jn,lcache)
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
  DO i=1,lag_rep%nbe
    IF(.NOT.lag_rep%linkage%leo(i))CYCLE
    jr=lag_rep%lbe(i)
    IF(.NOT.lag_rep%global%gbe(jr))CYCLE
    j=jr
    CALL mat%add_values(j,j,rop,1,1,1,1)
    CALL mat%add_values(j,j,rop,1,1,2,2)
    CALL mat%add_values(j,j,rop,1,1,3,3)
  END DO
ELSE IF(vbc_type>0)THEN
  DO i=1,lag_rep%nbe
    IF(.NOT.lag_rep%linkage%leo(i))CYCLE
    jr=lag_rep%lbe(i)
    IF(.NOT.lag_rep%global%gbe(jr))CYCLE
    j=jr
    CALL lag_vbc_diag(lag_rep,jr,vbc_type,mloc)
    DO jr=1,3
      DO jc=1,3
        call mat%add_values(j,j,mloc(jr,jc),1,1,jr,jc)
      END DO
    END DO
  END DO
END IF
DEALLOCATE(j,rop)
CALL vlag_rep%vec_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_lag_vgetmop
!------------------------------------------------------------------------------
!> Project a vector field onto a lagrange basis
!!
!! @note This subroutine only performs the integration of the field with
!! test functions for a Lagrange basis. To retrieve the correct projection the
!! result must be multiplied by the inverse of LAG::VMOP
!------------------------------------------------------------------------------
subroutine oft_lag_vproject(fe_rep,field,x)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(fem_interp), intent(inout) :: field !< Vector field for projection
class(oft_vector), intent(inout) :: x !< Field projected onto Lagrange basis
!---
real(r8) :: bcc(3),det,goptmp(3,4),vol
real(r8), pointer, dimension(:) :: xloc,yloc,zloc
real(r8), allocatable :: rop(:)
integer(i4) :: i,jc,m
integer(i4), allocatable :: j(:)
logical :: curved
CLASS(oft_scalar_fem), POINTER :: lag_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_lagrange_cast(lag_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_lag_vproject",__FILE__)
!---Initialize vectors to zero
NULLIFY(xloc,yloc,zloc)
call x%set(0.d0)
call x%get_local(xloc,1)
call x%get_local(yloc,2)
call x%get_local(zloc,3)
!---Integerate over the volume
!$omp parallel private(j,rop,curved,m,goptmp,vol,det,bcc,jc)
!---Allocate cell DOF arrays
allocate(j(lag_rep%nce),rop(lag_rep%nce))
!$omp do schedule(guided)
do i=1,lag_rep%mesh%nc ! Loop over cells
  call lag_rep%ncdofs(i,j) ! Get DOFs
  curved=cell_is_curved(lag_rep%mesh,i) ! Straight cell test
  do m=1,lag_rep%quad%np
    if(curved.OR.m==1)call lag_rep%mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
    det=vol*lag_rep%quad%wts(m)
    call field%interp(i,lag_rep%quad%pts(:,m),goptmp,bcc)
    call oft_lag_eval_all(lag_rep,i,lag_rep%quad%pts(:,m),rop)
    do jc=1,lag_rep%nce
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
!------------------------------------------------------------------------------
!> Compute the divergence of a Lagrange vector field
!------------------------------------------------------------------------------
subroutine lag_div(fe_rep,a,reg)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(oft_vector), intent(inout) :: a !< Input field
real(r8), intent(out) :: reg !< \f$ \int_v \nabla \cdot a \; dV \f$
real(r8), pointer, dimension(:) :: x,y,z
integer :: i,jc,m
integer(i4), allocatable :: j(:)
real(r8) :: goptmp(3,4),vol,det,div
real(r8), allocatable :: gop(:,:)
logical :: curved
CLASS(oft_scalar_fem), POINTER :: lag_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_lagrange_cast(lag_rep,fe_rep))CALL oft_abort("Incorrect FE type","lag_div",__FILE__)
!---Cast to vector type
NULLIFY(x,y,z)
CALL a%get_local(x,1)
CALL a%get_local(y,2)
CALL a%get_local(z,3)
!---
reg=0.d0
allocate(j(lag_rep%nce),gop(3,lag_rep%nce))
do i=1,lag_rep%mesh%nc
  curved=cell_is_curved(lag_rep%mesh,i) ! Straight cell test
  !---Get local to global DOF mapping
  call lag_rep%ncdofs(i,j)
  !---Get local reconstructed operators
  do m=1,lag_rep%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call lag_rep%mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
    det=vol*lag_rep%quad%wts(m)
    call oft_lag_geval_all(lag_rep,i,lag_rep%quad%pts(:,m),gop,goptmp)
    div=0.d0
    do jc=1,lag_rep%nce ! Loop over degrees of freedom
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
!------------------------------------------------------------------------------
!> Construct interpolation matrices on each MG level
!------------------------------------------------------------------------------
SUBROUTINE lag_setup_interp(ML_lag_rep,ML_vlag_rep)
CLASS(oft_ml_fem_type), intent(inout) :: ML_lag_rep
CLASS(oft_ml_fem_comp_type), optional, intent(inout) :: ML_vlag_rep
INTEGER(i4) :: i
DEBUG_STACK_PUSH
!---
DO i=ML_lag_rep%minlev+1,ML_lag_rep%nlevels
  CALL ML_lag_rep%set_level(i)
  !---
  if(ML_lag_rep%level==ML_lag_rep%blevel+1)then
    CYCLE
  end if
  !---Setup interpolation
  if(ML_lag_rep%current_level%order==1)then
    CALL lag_ginterpmatrix(ML_lag_rep%interp_matrices(ML_lag_rep%level)%m)
  else
    CALL lag_pinterpmatrix(ML_lag_rep%interp_matrices(ML_lag_rep%level)%m)
  end if
  CALL ML_lag_rep%interp_matrices(ML_lag_rep%level)%m%assemble
END DO
!---Create vector interpolation operator
IF(PRESENT(ML_vlag_rep))THEN
  CALL ML_vlag_rep%build_interp
END IF
DEBUG_STACK_POP
CONTAINS
!------------------------------------------------------------------------------
!> Construct interpolation matrix for polynomial levels
!------------------------------------------------------------------------------
SUBROUTINE lag_ginterpmatrix(mat)
class(oft_matrix), pointer, intent(inout) :: mat !< Needs docs
INTEGER(i4) :: i,j,k,m,icors,ifine,jb,i_ind(1),j_ind(1)
INTEGER(i4) :: etmp(2),ftmp(3),fetmp(3),ctmp(4),fc,ed
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pmap,emap,fmap
REAL(r8) :: f(4),val,mop(1)
REAL(r8), POINTER, DIMENSION(:,:) :: ed_nodes,fc_nodes,c_nodes
CLASS(oft_afem_type), POINTER :: lag_cors => NULL()
CLASS(oft_scalar_fem), POINTER :: lag_fine => NULL()
class(oft_mesh), pointer :: cmesh
CLASS(oft_vector), POINTER :: lag_vec_fine,lag_vec_cors
type(oft_graph_ptr), pointer :: graphs(:,:)
type(oft_graph), POINTER :: interp_graph
!---
if(ML_lag_rep%ml_mesh%level<1)call oft_abort('Invalid mesh level','lag_ginterpmatrix',__FILE__)
cmesh=>ML_lag_rep%ml_mesh%meshes(ML_lag_rep%ml_mesh%level-1)
if(cmesh%type/=1)CALL oft_abort("Only supported with tet meshes", &
  "lag_ginterpmatrix", __FILE__)
IF(.NOT.oft_3D_lagrange_cast(lag_fine,ML_lag_rep%current_level))CALL oft_abort("Incorrect FE type","lag_ginterpmatrix",__FILE__)
if(lag_fine%order/=1)then
  call oft_abort('Attempted geometric interpolation for pd > 1','lag_ginterpmatrix',__FILE__)
end if
lag_cors=>ML_lag_rep%levels(ML_lag_rep%level-1)%fe
ALLOCATE(ML_lag_rep%interp_graphs(ML_lag_rep%level)%g)
interp_graph=>ML_lag_rep%interp_graphs(ML_lag_rep%level)%g
!---Setup matrix sizes
interp_graph%nr=lag_fine%ne
interp_graph%nrg=lag_fine%global%ne
interp_graph%nc=lag_cors%ne
interp_graph%ncg=lag_cors%global%ne
!---Setup Matrix graph
ALLOCATE(interp_graph%kr(interp_graph%nr+1))
interp_graph%nnz=cmesh%np+2*cmesh%ne
ALLOCATE(interp_graph%lc(interp_graph%nnz))
interp_graph%lc=0_i4
!---Construct linkage
DO i=1,cmesh%np
  interp_graph%kr(i)=1
END DO
DO i=1,cmesh%ne
  interp_graph%kr(i+cmesh%np)=2
END DO
interp_graph%kr(interp_graph%nr+1)=interp_graph%nnz+1
do i=interp_graph%nr,1,-1 ! cumulative point to point count
  interp_graph%kr(i)=interp_graph%kr(i+1)-interp_graph%kr(i)
end do
if(interp_graph%kr(1)/=1)call oft_abort('Bad element to element count','oft_lag_ginterpmatrix',__FILE__)
DO i=1,cmesh%np
  interp_graph%lc(interp_graph%kr(i))=i
END DO
!---
DO i=1,cmesh%ne
  etmp=cmesh%le(:,i)
  !---
  ifine = cmesh%np+i
  jb=interp_graph%kr(ifine)-1
  DO k=1,2
    interp_graph%lc(jb+k)=etmp(k)
  END DO
END DO
!------------------------------------------------------------------------------
! Construct matrix
!------------------------------------------------------------------------------
NULLIFY(lag_vec_fine,lag_vec_cors)
CALL ML_lag_rep%vec_create(lag_vec_fine)
CALL ML_lag_rep%vec_create(lag_vec_cors,ML_lag_rep%level-1)
!---
ALLOCATE(graphs(1,1))
graphs(1,1)%g=>interp_graph
!---
CALL create_matrix(mat,graphs,lag_vec_fine,lag_vec_cors)
CALL lag_vec_fine%delete
CALL lag_vec_cors%delete
DeALLOCATE(graphs,lag_vec_fine,lag_vec_cors)
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
  jb=interp_graph%kr(ifine)-1
  DO k=1,2
    i_ind=ifine
    j_ind=etmp(k)
    mop=1.d0/2.d0
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
  END DO
END DO
deallocate(emap)
END SUBROUTINE lag_ginterpmatrix
!------------------------------------------------------------------------------
!> Construct interpolation matrix for polynomial levels
!------------------------------------------------------------------------------
SUBROUTINE lag_pinterpmatrix(mat)
class(oft_matrix), pointer, intent(inout) :: mat !< Needs docs
INTEGER(i4) :: i,j,k,m,icors,ifine,jb,js,jn,i_ind(1),j_ind(1)
INTEGER(i4) :: etmp(2),fc,ed,cell,dof,offset
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pmap,emap,fmap,jcors,ftmp,fetmp,ctmp
REAL(r8) :: f(4),val,mop(1)
REAL(r8), POINTER, DIMENSION(:,:) :: ed_nodes,fc_nodes,c_nodes
CLASS(oft_scalar_fem), POINTER :: lag_cors => NULL()
CLASS(oft_scalar_fem), POINTER :: lag_fine => NULL()
CLASS(oft_vector), POINTER :: lag_vec_fine,lag_vec_cors
type(oft_graph_ptr), pointer :: graphs(:,:)
type(oft_graph), POINTER :: interp_graph
class(oft_mesh), pointer :: mesh
IF(.NOT.oft_3D_lagrange_cast(lag_fine,ML_lag_rep%current_level))CALL oft_abort("Incorrect fine FE type","lag_pinterpmatrix",__FILE__)
mesh=>lag_fine%mesh
allocate(ftmp(mesh%face_np),fetmp(mesh%face_np),ctmp(mesh%cell_np))
!---
IF(.NOT.oft_3D_lagrange_cast(lag_cors,ML_lag_rep%levels(ML_lag_rep%level-1)%fe))CALL oft_abort("Incorrect coarse FE type","lag_pinterpmatrix",__FILE__)
CALL oft_lag_nodes(lag_fine%order,ed_nodes,fc_nodes,c_nodes)
ALLOCATE(ML_lag_rep%interp_graphs(ML_lag_rep%level)%g)
interp_graph=>ML_lag_rep%interp_graphs(ML_lag_rep%level)%g
!---Setup matrix sizes
interp_graph%nr=lag_fine%ne
interp_graph%nrg=lag_fine%global%ne
interp_graph%nc=lag_cors%ne
interp_graph%ncg=lag_cors%global%ne
!---Setup Matrix graph
ALLOCATE(interp_graph%kr(interp_graph%nr+1))
interp_graph%nnz=mesh%np
interp_graph%nnz=interp_graph%nnz + &
  (2+lag_cors%gstruct(2))*mesh%ne*lag_fine%gstruct(2)
interp_graph%nnz=interp_graph%nnz + &
  (mesh%face_np+lag_cors%gstruct(2)*mesh%face_np &
  + lag_cors%gstruct(3))*mesh%nf*lag_fine%gstruct(3)
interp_graph%nnz=interp_graph%nnz + &
  (mesh%cell_np+lag_cors%gstruct(2)*mesh%cell_ne &
  + lag_cors%gstruct(3)*mesh%cell_nf &
  + lag_cors%gstruct(4))*mesh%nc*lag_fine%gstruct(4)
ALLOCATE(interp_graph%lc(interp_graph%nnz))
interp_graph%lc=0_i4
!---Construct matrix
DO i=1,mesh%np
  interp_graph%kr(i)=1
END DO
!---
DO i=1,mesh%ne
  DO j=1,lag_fine%gstruct(2)
    ifine = mesh%np + (i-1)*lag_fine%gstruct(2) + j
    interp_graph%kr(ifine)=2+lag_cors%gstruct(2)
  END DO
END DO
!---
DO i=1,mesh%nf
  DO j=1,lag_fine%gstruct(3)
    ifine = mesh%np+lag_fine%gstruct(2)*mesh%ne &
      + (i-1)*lag_fine%gstruct(3) + j
    interp_graph%kr(ifine)= mesh%face_np &
      + mesh%face_np*lag_cors%gstruct(2)+lag_cors%gstruct(3)
  END DO
END DO
!---
DO i=1,mesh%nc
  DO j=1,lag_fine%gstruct(4)
    ifine = mesh%np + lag_fine%gstruct(2)*mesh%ne &
      + lag_fine%gstruct(3)*mesh%nf + (i-1)*lag_fine%gstruct(4) + j
    interp_graph%kr(ifine) = mesh%cell_np &
      + mesh%cell_ne*lag_cors%gstruct(2) &
      + mesh%cell_nf*lag_cors%gstruct(3) + lag_cors%gstruct(4)
  END DO
END DO
interp_graph%kr(interp_graph%nr+1)=interp_graph%nnz+1
do i=interp_graph%nr,1,-1 ! cumulative point to point count
  interp_graph%kr(i)=interp_graph%kr(i+1)-interp_graph%kr(i)
end do
if(interp_graph%kr(1)/=1)call oft_abort('Bad element to element count','oft_lag_interpmatrix',__FILE__)
!---Construct matrix
DO i=1,mesh%np
  interp_graph%lc(interp_graph%kr(i))=i
END DO
!---
DO i=1,mesh%ne
  etmp=mesh%le(:,i)
  DO j=1,lag_fine%gstruct(2)
    ifine = mesh%np + (i-1)*lag_fine%gstruct(2) + j
    jb=interp_graph%kr(ifine)-1
    DO k=1,2
      interp_graph%lc(jb+k)=etmp(k)
    END DO
    DO k=1,lag_cors%gstruct(2)
      interp_graph%lc(jb+2+k)=mesh%np + (i-1)*lag_cors%gstruct(2) + k
    END DO
    !---
    js=interp_graph%kr(ifine)
    jn=interp_graph%kr(ifine+1)-1
    CALL sort_array(interp_graph%lc(js:jn),jn-js+1)
  END DO
END DO
!---
DO i=1,mesh%nf
  ftmp=mesh%lf(:,i)
  fetmp=mesh%lfe(:,i)
  DO j=1,lag_fine%gstruct(3)
    ifine = mesh%np+lag_fine%gstruct(2)*mesh%ne &
      + (i-1)*lag_fine%gstruct(3) + j
    jb=interp_graph%kr(ifine)-1
    DO k=1,mesh%face_np
      interp_graph%lc(jb+k)=ftmp(k)
      ! etmp=mesh%bmesh%face_ed(:,k)
      ed=fetmp(k)
      DO m=1,lag_cors%gstruct(2)
        interp_graph%lc(jb+mesh%face_np+(m-1)*mesh%face_np+k)= &
          mesh%np + (ABS(ed)-1)*lag_cors%gstruct(2) + m
      END DO
    END DO
    offset=jb+mesh%face_np+mesh%face_np*lag_cors%gstruct(2)
    DO m=1,lag_cors%gstruct(3)
      interp_graph%lc(offset+m)= &
        mesh%np + lag_cors%gstruct(2)*mesh%ne + (i-1)*lag_cors%gstruct(3) + m
    END DO
    !---
    js=interp_graph%kr(ifine)
    jn=interp_graph%kr(ifine+1)-1
    CALL sort_array(interp_graph%lc(js:jn),jn-js+1)
  END DO
END DO
!---
DO i=1,mesh%nc
  ctmp=mesh%lc(:,i)
  DO j=1,lag_fine%gstruct(4)
    ifine = mesh%np + lag_fine%gstruct(2)*mesh%ne &
      + lag_fine%gstruct(3)*mesh%nf + (i-1)*lag_fine%gstruct(4) + j
    jb=interp_graph%kr(ifine)-1
    DO k=1,mesh%cell_np
      interp_graph%lc(jb+k)=ctmp(k)
    END DO
    offset=jb+mesh%cell_np
    DO k=1,mesh%cell_ne
      etmp=mesh%cell_ed(:,k)
      ed=mesh%lce(k,i)
      DO m=1,lag_cors%gstruct(2)
        interp_graph%lc(offset+(m-1)*mesh%cell_ne+k)= &
          mesh%np + (ABS(ed)-1)*lag_cors%gstruct(2) + m
      END DO
    END DO
    offset=jb+mesh%cell_np+mesh%cell_ne*lag_cors%gstruct(2)
    DO k=1,mesh%cell_nf
      ftmp=mesh%cell_fc(:,k)
      DO m=1,lag_cors%gstruct(3)
        interp_graph%lc(offset+(m-1)*mesh%cell_nf+k)=mesh%np &
        + lag_cors%gstruct(2)*mesh%ne + (ABS(mesh%lcf(k,i))-1)*lag_cors%gstruct(3) + m
      END DO
    END DO
    offset=jb+mesh%cell_np+mesh%cell_ne*lag_cors%gstruct(2)+mesh%cell_nf*lag_cors%gstruct(3)
    DO m=1,lag_cors%gstruct(4)
      interp_graph%lc(offset+m)=mesh%np &
      + lag_cors%gstruct(2)*mesh%ne+lag_cors%gstruct(3)*mesh%nf &
      + (i-1)*lag_cors%gstruct(4) + m
    END DO
    !---
    js=interp_graph%kr(ifine)
    jn=interp_graph%kr(ifine+1)-1
    CALL sort_array(interp_graph%lc(js:jn),jn-js+1)
  END DO
END DO
!------------------------------------------------------------------------------
! Construct matrix
!------------------------------------------------------------------------------
NULLIFY(lag_vec_fine,lag_vec_cors)
CALL ML_lag_rep%vec_create(lag_vec_fine)
CALL ML_lag_rep%vec_create(lag_vec_cors,ML_lag_rep%level-1)
!---
ALLOCATE(graphs(1,1))
graphs(1,1)%g=>interp_graph
!---
CALL create_matrix(mat,graphs,lag_vec_fine,lag_vec_cors)
CALL lag_vec_fine%delete
CALL lag_vec_cors%delete
DeALLOCATE(graphs,lag_vec_fine,lag_vec_cors)
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
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
  DO j=1,lag_fine%gstruct(2)
    ifine = mesh%np + (i-1)*lag_fine%gstruct(2) + j
    dof = mesh%cell_np + (ed-1)*lag_fine%gstruct(2) + j
    CALL oft_lag_npos(lag_fine,cell,dof,f)
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
  DO j=1,lag_fine%gstruct(3)
    ifine = mesh%np + lag_fine%gstruct(2)*mesh%ne &
          + (i-1)*lag_fine%gstruct(3) + j
    dof = mesh%cell_np + lag_fine%gstruct(2)*mesh%cell_ne &
          + (fc-1)*lag_fine%gstruct(3) + j
    CALL oft_lag_npos(lag_fine,cell,dof,f)
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
  DO j=1,lag_fine%gstruct(4)
    ifine = mesh%np + lag_fine%gstruct(2)*mesh%ne &
          + lag_fine%gstruct(3)*mesh%nf + (i-1)*lag_fine%gstruct(4) + j
    dof = mesh%cell_np + lag_fine%gstruct(2)*mesh%cell_ne &
          + lag_fine%gstruct(3)*mesh%cell_nf &
          + j
    CALL oft_lag_npos(lag_fine,cell,dof,f)
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
END SUBROUTINE lag_pinterpmatrix
END SUBROUTINE lag_setup_interp
!------------------------------------------------------------------------------
!> Transfer a base level Lagrange scalar field to the next MPI level
!------------------------------------------------------------------------------
subroutine lag_base_pop(self,acors,afine)
class(oft_ml_fe_vecspace), intent(inout) :: self
class(oft_vector), intent(inout) :: acors !< Vector to transfer
class(oft_vector), intent(inout) :: afine !< Fine vector from transfer
integer(i4), pointer, dimension(:) :: lptmp
integer(i4) :: i
real(r8), pointer, dimension(:) :: array_c,array_f
DEBUG_STACK_PUSH
lptmp=>self%ML_FE_rep%ml_mesh%meshes(self%ML_FE_rep%ml_mesh%nbase+1)%base%lp
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
!------------------------------------------------------------------------------
!> Transfer a MPI level Lagrange scalar field to the base level
!------------------------------------------------------------------------------
subroutine lag_base_push(self,afine,acors)
class(oft_ml_fe_vecspace), intent(inout) :: self
class(oft_vector), intent(inout) :: afine !< Needs docs
class(oft_vector), intent(inout) :: acors !< Needs docs
integer(i4), pointer :: lptmp(:)
integer(i4) :: i,j,ierr
real(r8), pointer, dimension(:) :: alias,array_c,array_f
CLASS(oft_afem_type), POINTER :: lag_fine => NULL()
DEBUG_STACK_PUSH
!---
lptmp=>self%ML_FE_rep%ml_mesh%meshes(self%ML_FE_rep%ml_mesh%nbase+1)%base%lp
CALL acors%get_local(array_c)
CALL afine%get_local(array_f)
lag_fine=>self%ML_FE_rep%levels(self%ML_FE_rep%level+1)%fe
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
!------------------------------------------------------------------------------
!> Compute eigenvalues and smoothing coefficients for the operator LAG::LOP
!------------------------------------------------------------------------------
SUBROUTINE lag_lop_eigs(ML_lag_rep,minlev)
type(oft_ml_fem_type), target, intent(inout) :: ML_lag_rep
INTEGER(i4), INTENT(in) :: minlev !< Needs docs
#ifdef HAVE_ARPACK
INTEGER(i4) :: i
REAL(r8) :: lam0
REAL(r8), ALLOCATABLE :: df(:)
CLASS(oft_vector), POINTER :: u
TYPE(oft_irlm_eigsolver) :: arsolver
CLASS(oft_matrix), POINTER :: md => NULL()
CLASS(oft_matrix), POINTER :: lop => NULL()
TYPE(oft_lag_zerob), TARGET :: bc_tmp
DEBUG_STACK_PUSH
bc_tmp%ML_lag_rep=>ML_lag_rep
!------------------------------------------------------------------------------
! Compute optimal smoother coefficients
!------------------------------------------------------------------------------
IF(oft_env%head_proc)WRITE(*,*)'Optimizing Jacobi damping for LAG::LOP'
ALLOCATE(df(ML_lag_rep%nlevels))
df=0.d0
DO i=minlev,ML_lag_rep%nlevels
  CALL ML_lag_rep%set_level(i)
  !---Create fields
  CALL ML_lag_rep%vec_create(u)
  !---Get Ev range
  NULLIFY(lop)
  SELECT TYPE(this=>ML_lag_rep%current_level)
  CLASS IS(oft_scalar_fem)
    CALL oft_lag_getlop(this,lop,'zerob')
  CLASS DEFAULT
    CALL oft_abort("Error getting current FE rep","lag_lop_eigs",__FILE__)
  END SELECT
  ! CALL oft_lag_getlop(lag_rep,lop,'zerob')
  CALL create_diagmatrix(md,lop%D)
  !---
  arsolver%A=>lop
  arsolver%M=>md
  arsolver%mode=2
  arsolver%tol=1.E-5_r8
  arsolver%bc=>bc_tmp
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
  DO i=1,ML_lag_rep%nlevels-1
    WRITE(*,'(F5.3,A)',ADVANCE='NO')df(i),', '
  END DO
  WRITE(*,'(F5.3,A)')df(ML_lag_rep%nlevels)
END IF
DEALLOCATE(df)
DEBUG_STACK_POP
#else
CALL oft_abort("Subroutine requires ARPACK", "lag_lop_eigs", __FILE__)
#endif
END SUBROUTINE lag_lop_eigs
!------------------------------------------------------------------------------
!> Construct default MG preconditioner for LAG::LOP
!------------------------------------------------------------------------------
SUBROUTINE lag_getlop_pre(ML_lag_rep,pre,mats,level,nlevels)
type(oft_ml_fem_type), target, intent(inout) :: ML_lag_rep
CLASS(oft_solver), POINTER, INTENT(out) :: pre !< Needs docs
TYPE(oft_matrix_ptr), POINTER, INTENT(inout) :: mats(:) !< Needs docs
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Needs docs
INTEGER(i4), OPTIONAL, INTENT(in) :: nlevels !< Needs docs
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
TYPE(oft_lag_zerob), POINTER :: bc_tmp
TYPE(oft_ml_fe_vecspace), POINTER :: tmp_vecspace
!---
TYPE(xml_node), POINTER :: pre_node
#ifdef HAVE_XML
integer(i4) :: nnodes
TYPE(xml_node), POINTER :: lag_node
#endif
DEBUG_STACK_PUSH
!---
minlev=1
toplev=ML_lag_rep%level
levin=ML_lag_rep%level
IF(PRESENT(level))toplev=level
IF(PRESENT(nlevels))minlev=toplev-nlevels+1
nl=toplev-minlev+1
!---
IF(minlev<ML_lag_rep%minlev)CALL oft_abort('Minimum level is < minlev','lag_getlop_pre',__FILE__)
IF(toplev>ML_lag_rep%nlevels)CALL oft_abort('Maximum level is > lag_nlevels','lag_getlop_pre',__FILE__)
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
  CALL ML_lag_rep%set_level(minlev+(i-1))
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
    CALL oft_lag_getlop(ML_lag_rep%current_level,mats(i)%M,'zerob')
  END IF
  IF(i>1)ml_int(i-1)%M=>ML_lag_rep%interp_matrices(ML_lag_rep%level)%m !oft_lagrange_ops%interp
END DO
CALL ML_lag_rep%set_level(levin)
!------------------------------------------------------------------------------
! Search for XML-spec
!------------------------------------------------------------------------------
NULLIFY(pre_node)
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml))THEN
  CALL xml_get_element(oft_env%xml,"lagrange",lag_node,ierr)
  IF(ierr==0)CALL xml_get_element(lag_node,"lop",pre_node,ierr)
END IF
#endif
!------------------------------------------------------------------------------
! Setup preconditioner
!------------------------------------------------------------------------------
ALLOCATE(bc_tmp)
bc_tmp%ML_lag_rep=>ML_lag_rep
NULLIFY(pre)
ALLOCATE(tmp_vecspace)
tmp_vecspace%ML_FE_rep=>ML_lag_rep
tmp_vecspace%base_pop=>lag_base_pop
tmp_vecspace%base_push=>lag_base_push
CALL create_mlpre(pre,mats(1:nl),levels,nlevels=nl,ml_vecspace=tmp_vecspace, &
       bc=bc_tmp,stype=1,df=df,nu=nu,xml_root=pre_node)
!------------------------------------------------------------------------------
! Cleanup
!------------------------------------------------------------------------------
DEALLOCATE(ml_int,levels,df,nu)
DEBUG_STACK_POP
END SUBROUTINE lag_getlop_pre
end module oft_lag_operators
