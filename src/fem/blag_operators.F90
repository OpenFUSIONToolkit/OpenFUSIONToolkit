!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_blag_operators.F90
!
!> Surface lagrange FE operator definitions
!! - Operator construction
!!   - MOP: mass matrix        \f$ \int \left( u^T v \right) dV \f$
!!   - LOP: laplacian matrix   \f$ \int \left( \nabla u^T \cdot \nabla v \right) dV \f$
!! - Interpolation classes
!! - Field projection (to plotting mesh)
!! - Boundary conditions
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_lag
!---------------------------------------------------------------------------
MODULE oft_blag_operators
USE oft_base
USE oft_sort, ONLY: sort_array
USE oft_mesh_type, ONLY: oft_mesh, oft_bmesh, cell_is_curved
!---
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr, &
  oft_graph, oft_graph_ptr
USE oft_deriv_matrices, ONLY: oft_diagmatrix, create_diagmatrix
USE oft_solver_base, ONLY: oft_solver, oft_solver_bc
#ifdef HAVE_ARPACK
USE oft_arpack, ONLY: oft_irlm_eigsolver
#endif
!---
USE fem_base, ONLY: oft_fem_type, oft_bfem_type, fem_max_levels, oft_ml_fem_type, &
  oft_afem_type
USE fem_utils, ONLY: fem_interp, bfem_interp
USE oft_lag_basis, ONLY: oft_lag_eval_all, oft_lag_geval_all, oft_lag_eval, oft_blag_d2eval, &
  oft_lag_nodes, oft_blag_eval, oft_2D_lagrange_cast, &
  oft_blag_geval, oft_lag_npos, oft_scalar_fem, oft_scalar_bfem
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------
!> Interpolate a surface Lagrange field
!---------------------------------------------------------------------------
type, extends(bfem_interp) :: oft_lag_brinterp
  logical :: own_vals = .TRUE. !< 
  class(oft_vector), pointer :: u => NULL() !< Field for interpolation
  real(r8), pointer, dimension(:) :: vals => NULL() !< Local values
  class(oft_scalar_bfem), pointer :: lag_rep => NULL() !< Lagrange FE representation
contains
  !> Retrieve local values for interpolation
  procedure :: setup => lag_brinterp_setup
  procedure :: shared_setup => lag_brinterp_share
  !> Reconstruct field
  procedure :: interp => lag_brinterp
  !> Delete reconstruction object
  procedure :: delete => lag_brinterp_delete
end type oft_lag_brinterp
!---------------------------------------------------------------------------
!> Interpolate \f$ \nabla \f$ of a Lagrange field
!---------------------------------------------------------------------------
type, extends(oft_lag_brinterp) :: oft_lag_bginterp
contains
  !> Reconstruct field
  procedure :: interp => lag_bginterp
end type oft_lag_bginterp
!---------------------------------------------------------------------------
!> Interpolate \f$ \frac{\partial }{\partial x_i \partial x_j} \f$ of a Lagrange field
!---------------------------------------------------------------------------
type, extends(oft_lag_brinterp) :: oft_lag_bg2interp
contains
  !> Reconstruct field
  procedure :: interp => lag_bg2interp
end type oft_lag_bg2interp
!---------------------------------------------------------------------------
!> Interpolate a boundary Lagrange vector field
!---------------------------------------------------------------------------
type, extends(bfem_interp) :: oft_lag_bvrinterp
  class(oft_vector), pointer :: u => NULL() !< Field for interpolation
  real(r8), pointer, dimension(:,:) :: vals => NULL() !< Local values
  class(oft_scalar_bfem), pointer :: lag_rep => NULL() !< Lagrange FE representation
contains
  !> Retrieve local values for interpolation
  procedure :: setup => lag_bvrinterp_setup
  !> Reconstruct field
  procedure :: interp => lag_bvrinterp
  !> Delete reconstruction object
  procedure :: delete => lag_bvrinterp_delete
end type oft_lag_bvrinterp
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
type, extends(oft_solver_bc) :: oft_blag_zerob
  class(oft_ml_fem_type), pointer :: ML_lag_rep => NULL() !< FE representation
contains
  procedure :: apply => zerob_apply
  procedure :: delete => zerob_delete
end type oft_blag_zerob
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
type, extends(oft_blag_zerob) :: oft_blag_zerogrnd
contains
  procedure :: apply => zerogrnd_apply
end type oft_blag_zerogrnd
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
type, extends(oft_blag_zerob) :: oft_blag_zeroe
  INTEGER(i4), POINTER, DIMENSION(:) :: parent_geom_flag
contains
  procedure :: apply => zeroe_apply
end type oft_blag_zeroe
!---Pre options
integer(i4) :: nu_lop_surf(fem_max_levels)=0
real(r8) :: df_lop_surf(fem_max_levels)=-1.d99
!---Cache variables
REAL(r8), POINTER, DIMENSION(:) :: oft_blag_rop => NULL()
REAL(r8), POINTER, DIMENSION(:,:) :: oft_blag_gop => NULL()
!$omp threadprivate(oft_blag_rop,oft_blag_gop)
contains
!---------------------------------------------------------------------------
!> Setup interpolator for boundary Lagrange scalar fields
!!
!! Fetches local representation used for interpolation from solution vector
!!
!! @note Should only be used via class \ref oft_lag_brinterp or children
!---------------------------------------------------------------------------
subroutine lag_brinterp_setup(self,lag_rep)
class(oft_lag_brinterp), intent(inout) :: self
class(oft_afem_type), target, intent(inout) :: lag_rep
IF(.NOT.oft_2D_lagrange_cast(self%lag_rep,lag_rep))CALL oft_abort("Incorrect FE type","lag_brinterp_setup",__FILE__)
self%mesh=>self%lag_rep%mesh
!---Get local slice
CALL self%u%get_local(self%vals)
end subroutine lag_brinterp_setup
!---------------------------------------------------------------------------
!> Setup interpolator by linking to another interpolator of the same class
!!
!! @note Should only be used via class \ref oft_lag_brinterp or children
!---------------------------------------------------------------------------
subroutine lag_brinterp_share(self,source_obj)
class(oft_lag_brinterp), intent(inout) :: self
class(oft_lag_brinterp), intent(in) :: source_obj
IF(.NOT.ASSOCIATED(source_obj%lag_rep).OR..NOT.ASSOCIATED(source_obj%vals))THEN
  CALL oft_abort("Source object not setup","lag_brinterp_share",__FILE__)
END IF
self%vals=>source_obj%vals
self%lag_rep=>source_obj%lag_rep
self%mesh=>source_obj%mesh
self%u=>source_obj%u
self%own_vals=.FALSE.
end subroutine lag_brinterp_share
!---------------------------------------------------------------------------
!> Destroy temporary internal storage
!!
!! @note Should only be used via class \ref oft_lag_brinterp or children
!---------------------------------------------------------------------------
subroutine lag_brinterp_delete(self)
class(oft_lag_brinterp), intent(inout) :: self
!---Destroy locals
IF(self%own_vals)THEN
  IF(ASSOCIATED(self%vals))DEALLOCATE(self%vals)
ELSE
  NULLIFY(self%vals)
END IF
NULLIFY(self%lag_rep,self%mesh,self%u)
end subroutine lag_brinterp_delete
!---------------------------------------------------------------------------
!> Reconstruct a surface Lagrange scalar field
!!
!! @note Should only be used via class \ref oft_lag_brinterp
!---------------------------------------------------------------------------
subroutine lag_brinterp(self,cell,f,gop,val)
class(oft_lag_brinterp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [3]
real(r8), intent(in) :: gop(3,3) !< Logical gradient vectors at `f` [3,3]
real(r8), intent(out) :: val(:) !< Reconstructed field at `f` [1]
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: rop(1)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!','lag_brinterp',__FILE__)
!---Get dofs
allocate(j(self%lag_rep%nce))
call self%lag_rep%ncdofs(cell,j) ! get DOFs
!---Reconstruct field
val=0.d0
do jc=1,self%lag_rep%nce
  call oft_blag_eval(self%lag_rep,cell,jc,f,rop(1))
  val=val+self%vals(j(jc))*rop
end do
deallocate(j)
DEBUG_STACK_POP
end subroutine lag_brinterp
!---------------------------------------------------------------------------
!> Reconstruct the gradient of a surface Lagrange scalar field
!!
!! @note Should only be used via class \ref lag_bginterp
!---------------------------------------------------------------------------
subroutine lag_bginterp(self,cell,f,gop,val)
class(oft_lag_bginterp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [3]
real(r8), intent(in) :: gop(3,3) !< Logical gradient vectors at `f` [3,3]
real(r8), intent(out) :: val(:) !< Reconstructed field at `f` [3]
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: rop(3)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!','lag_bginterp',__FILE__)
!---Get dofs
allocate(j(self%lag_rep%nce))
call self%lag_rep%ncdofs(cell,j) ! get DOFs
!---Reconstruct field
val=0.d0
do jc=1,self%lag_rep%nce
  call oft_blag_geval(self%lag_rep,cell,jc,f,rop,gop)
  val=val+self%vals(j(jc))*rop
end do
deallocate(j)
DEBUG_STACK_POP
end subroutine lag_bginterp
!---------------------------------------------------------------------------
!> Reconstruct the Hessian of a surface Lagrange scalar field
!!
!! @note Should only be used via class \ref oft_lag_bg2interp
!---------------------------------------------------------------------------
subroutine lag_bg2interp(self,cell,f,gop,val)
class(oft_lag_bg2interp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [3]
real(r8), intent(in) :: gop(3,3) !< Logical gradient vectors at `f` [3,3]
real(r8), intent(out) :: val(:) !< Reconstructed field at `f` [6]
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: rop(3),r2op(6)
real(8) :: g2op(6,6),Kmat(6,3)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!','lag_bg2interp',__FILE__)
!---Get dofs
allocate(j(self%lag_rep%nce))
call self%lag_rep%ncdofs(cell,j) ! get DOFs
!---Reconstruct field
call self%lag_rep%mesh%hessian(cell,f,g2op,Kmat)
val=0.d0
do jc=1,self%lag_rep%nce
  CALL oft_blag_geval(self%lag_rep,cell,jc,f,rop,gop)
  call oft_blag_d2eval(self%lag_rep,cell,jc,f,r2op,g2op)
  val=val+self%vals(j(jc))*(r2op-MATMUL(g2op,MATMUL(Kmat,rop)))
end do
deallocate(j)
DEBUG_STACK_POP
end subroutine lag_bg2interp
!---------------------------------------------------------------------------
!> Setup interpolator for boundary Lagrange vector fields
!!
!! Fetches local representation used for interpolation from vector object
!!
!! @note Should only be used via class \ref oft_lag_bvrinterp or children
!---------------------------------------------------------------------------
subroutine lag_bvrinterp_setup(self,lag_rep)
class(oft_lag_bvrinterp), intent(inout) :: self
class(oft_afem_type), target, intent(inout) :: lag_rep
real(r8), pointer, dimension(:) :: vtmp
IF(.NOT.oft_2D_lagrange_cast(self%lag_rep,lag_rep))CALL oft_abort("Incorrect FE type","lag_bvrinterp_setup",__FILE__)
self%mesh=>self%lag_rep%mesh
!---Get local slice
IF(.NOT.ASSOCIATED(self%vals))ALLOCATE(self%vals(3,self%lag_rep%ne))
vtmp=>self%vals(1,:)
CALL self%u%get_local(vtmp,1)
vtmp=>self%vals(2,:)
CALL self%u%get_local(vtmp,2)
vtmp=>self%vals(3,:)
CALL self%u%get_local(vtmp,3)
end subroutine lag_bvrinterp_setup
!---------------------------------------------------------------------------
!> Destroy temporary internal storage
!!
!! @note Should only be used via class \ref oft_lag_bvrinterp or children
!---------------------------------------------------------------------------
subroutine lag_bvrinterp_delete(self)
class(oft_lag_bvrinterp), intent(inout) :: self
!---Destroy locals
IF(ASSOCIATED(self%vals))DEALLOCATE(self%vals)
NULLIFY(self%lag_rep,self%u)
end subroutine lag_bvrinterp_delete
!---------------------------------------------------------------------------
!> Reconstruct a boundary Lagrange vector field
!!
!! @note Should only be used via class \ref oft_lag_bvrinterp
!---------------------------------------------------------------------------
subroutine lag_bvrinterp(self,cell,f,gop,val)
class(oft_lag_bvrinterp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [3]
real(r8), intent(in) :: gop(3,3) !< Logical gradient vectors at f [3,3]
real(r8), intent(out) :: val(:) !< Reconstructed field at f [3]
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: rop(1)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!', &
'lag_bvrinterp',__FILE__)
!---Get dofs
allocate(j(self%lag_rep%nce))
call self%lag_rep%ncdofs(cell,j) ! get DOFs
!---Reconstruct field
val=0.d0
do jc=1,self%lag_rep%nce
  call oft_blag_eval(self%lag_rep,cell,jc,f,rop(1))
  val=val+self%vals(:,j(jc))*rop(1)
end do
deallocate(j)
DEBUG_STACK_POP
end subroutine lag_bvrinterp
!---------------------------------------------------------------------------
!> Zero a surface Lagrange scalar field at all boundary nodes
!---------------------------------------------------------------------------
subroutine zerob_apply(self,a)
class(oft_blag_zerob), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to be zeroed
integer(i4) :: i,j
real(r8), pointer, dimension(:) :: vloc
DEBUG_STACK_PUSH
NULLIFY(vloc)
!---Cast to vector type
NULLIFY(vloc)
call a%get_local(vloc)
!---Zero boundary values
!$omp parallel do private(j)
do i=1,self%ML_lag_rep%current_level%nbe
  j=self%ML_lag_rep%current_level%lbe(i)
  ! vloc(j)=0.d0
  IF(self%ML_lag_rep%current_level%global%gbe(j))vloc(j)=0.d0
end do
!---
call a%restore_local(vloc)
deallocate(vloc)
DEBUG_STACK_POP
end subroutine zerob_apply
!---------------------------------------------------------------------------
!> Zero a surface Lagrange scalar field at all boundary nodes
!---------------------------------------------------------------------------
subroutine zerob_delete(self)
class(oft_blag_zerob), intent(inout) :: self
NULLIFY(self%ML_lag_rep)
end subroutine zerob_delete
!---------------------------------------------------------------------------
!> Zero a surface Lagrange scalar field at the mesh "grounding" node
!!
!! @note Presently the first boundary node is used as the "grounding" node
!---------------------------------------------------------------------------
subroutine zerogrnd_apply(self,a)
class(oft_blag_zerogrnd), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to be zeroed
integer(i4) :: i,j
real(r8), pointer, dimension(:) :: vloc
DEBUG_STACK_PUSH
NULLIFY(vloc)
!---Cast to vector type
NULLIFY(vloc)
call a%get_local(vloc)
vloc(self%ML_lag_rep%current_level%lbe(1))=0.d0
! !---Zero boundary values
! !$omp parallel do
! do i=1,oft_blagrange%nbe
!   vloc(oft_blagrange%lbe(i))=0.d0
! end do
call a%restore_local(vloc)
deallocate(vloc)
DEBUG_STACK_POP
end subroutine zerogrnd_apply
!---------------------------------------------------------------------------
!> Zero a surface Lagrange scalar field at all edge nodes
!---------------------------------------------------------------------------
subroutine zeroe_apply(self,a)
class(oft_blag_zeroe), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to be zeroed
integer(i4) :: i,j
real(r8), pointer, dimension(:) :: vloc
class(oft_scalar_bfem), pointer :: this
DEBUG_STACK_PUSH
IF(.NOT.oft_2D_lagrange_cast(this,self%ML_lag_rep%current_level))CALL oft_abort("Incorrect FE type","zeroe_apply",__FILE__)
!---Cast to vector type
NULLIFY(vloc)
call a%get_local(vloc)
!---Zero boundary values
!$omp parallel do private(j)
do i=1,this%ne
  j=this%parent%le(i)
  if(self%parent_geom_flag(j)/=3)vloc(i)=0.d0
end do
!---
call a%restore_local(vloc)
deallocate(vloc)
DEBUG_STACK_POP
end subroutine zeroe_apply
!---------------------------------------------------------------------------
!> Construct mass matrix for a boundary Lagrange scalar representation
!!
!! Supported boundary conditions
!! - `'none'` Full matrix
!! - `'zerob'` Dirichlet for all boundary DOF
!---------------------------------------------------------------------------
subroutine oft_blag_getmop(fe_rep,mat,bc)
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
CLASS(oft_scalar_bfem), POINTER :: lag_rep
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing Boundary LAG::MOP'
  CALL mytimer%tick()
END IF
IF(.NOT.oft_2D_lagrange_cast(lag_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_blag_getmop",__FILE__)
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL lag_rep%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,rop,det,mop,curved,goptmp,m,vol,jc,jr)
allocate(j(lag_rep%nce)) ! Local DOF and matrix indices
allocate(rop(lag_rep%nce)) ! Reconstructed gradient operator
allocate(mop(lag_rep%nce,lag_rep%nce)) ! Local laplacian matrix
!$omp do
do i=1,lag_rep%mesh%nc
  !---Get local reconstructed operators
  mop=0.d0
  do m=1,lag_rep%quad%np ! Loop over quadrature points
    call lag_rep%mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
    det=vol*lag_rep%quad%wts(m)
    do jc=1,lag_rep%nce ! Loop over degrees of freedom
      call oft_blag_eval(lag_rep,i,jc,lag_rep%quad%pts(:,m),rop(jc))
    end do
    !---Compute local matrix contributions
    do jr=1,lag_rep%nce
      do jc=1,lag_rep%nce
        mop(jr,jc) = mop(jr,jc) + rop(jr)*rop(jc)*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call lag_rep%ncdofs(i,j)
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,mop,lag_rep%nce,lag_rep%nce)
  ! !$omp end critical
end do
deallocate(j,rop,mop)
!$omp end parallel
CALL lag_rep%vec_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_blag_getmop
!---------------------------------------------------------------------------
!> Construct laplacian matrix for Lagrange scalar representation
!!
!! Supported boundary conditions
!! - `'none'` Full matrix
!! - `'zerob'` Dirichlet for all boundary DOF
!! - `'grnd'`  Dirichlet for only groundin point
!---------------------------------------------------------------------------
subroutine oft_blag_getlop(fe_rep,mat,bc,parent_geom_flag)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix object
character(LEN=*), intent(in) :: bc !< Boundary condition
INTEGER(i4), optional, intent(in) :: parent_geom_flag(:) !< Parent FE geometry type flag
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time
real(r8), allocatable :: gop(:,:),lop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_lag_vec
type(oft_timer) :: mytimer
CLASS(oft_scalar_bfem), POINTER :: lag_rep
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing Boundary LAG::LOP'
  CALL mytimer%tick()
END IF
IF(.NOT.oft_2D_lagrange_cast(lag_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_blag_getlop",__FILE__)
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL lag_rep%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!---------------------------------------------------------------------------
! Operator integration
!---------------------------------------------------------------------------
!$omp parallel private(j,gop,det,lop,curved,goptmp,m,vol,jc,jr)
allocate(j(lag_rep%nce)) ! Local DOF and matrix indices
allocate(gop(3,lag_rep%nce)) ! Reconstructed gradient operator
allocate(lop(lag_rep%nce,lag_rep%nce)) ! Local laplacian matrix
!$omp do
do i=1,lag_rep%mesh%nc
  !---Get local reconstructed operators
  lop=0.d0
  do m=1,lag_rep%quad%np ! Loop over quadrature points
    call lag_rep%mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
    det=vol*lag_rep%quad%wts(m)
    do jc=1,lag_rep%nce ! Loop over degrees of freedom
      call oft_blag_geval(lag_rep,i,jc,lag_rep%quad%pts(:,m),gop(:,jc),goptmp)
    end do
    !---Compute local matrix contributions
    do jr=1,lag_rep%nce
      do jc=1,lag_rep%nce
        lop(jr,jc) = lop(jr,jc) + DOT_PRODUCT(gop(:,jr),gop(:,jc))*det
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
    CASE("edges")
      DO jr=1,lag_rep%nce
        jc=lag_rep%parent%le(j(jr))
        IF(parent_geom_flag(jc)/=3)lop(jr,:)=0.d0
      END DO
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
      IF(.NOT.lag_rep%linkage%leo(i))CYCLE
      jr=lag_rep%lbe(i)
      IF(lag_rep%global%gbe(jr))THEN
        j=jr
        call mat%add_values(j,j,lop,1,1)
      END IF
    END DO
  CASE("grnd")
    IF(ANY(lag_rep%mesh%igrnd>0))THEN
      DO i=1,lag_rep%nbe
        IF(.NOT.lag_rep%linkage%leo(i))CYCLE
        jr=lag_rep%lbe(i)
        IF(ANY(jr==lag_rep%mesh%igrnd))THEN
          j=jr
          call mat%add_values(j,j,lop,1,1)
        END IF
      END DO
    END IF
  CASE("edges")
    DO i=1,lag_rep%ne
      IF(lag_rep%be(i))CYCLE
      jr=lag_rep%parent%le(i)
      IF(parent_geom_flag(jr)/=3)THEN
        j=i
        call mat%add_values(j,j,lop,1,1)
      END IF
    END DO
    DO i=1,lag_rep%nbe
      IF(.NOT.lag_rep%linkage%leo(i))CYCLE
      jc=lag_rep%lbe(i)
      jr=lag_rep%parent%le(jc)
      IF(parent_geom_flag(jr)/=3)THEN
        j=jc
        call mat%add_values(j,j,lop,1,1)
      END IF
    END DO
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
end subroutine oft_blag_getlop
!------------------------------------------------------------------------------
!> Project a scalar field onto a boundary Lagrange basis
!!
!! @note This subroutine only performs the integration of the field with
!! boundary test functions for a Lagrange basis.
!------------------------------------------------------------------------------
SUBROUTINE oft_blag_project(fe_rep,field,x)
class(oft_afem_type), target, intent(inout) :: fe_rep
CLASS(bfem_interp), INTENT(inout) :: field !< Scalar field for projection
CLASS(oft_vector), INTENT(inout) :: x !< Field projected onto boundary Lagrange basis
INTEGER(i4) :: i,m,k,jc,cell,ptmap(3)
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: j
REAL(r8) :: vol,det,etmp(1),sgop(3,3),rop
REAL(r8), POINTER, DIMENSION(:) :: xloc
CLASS(oft_scalar_bfem), POINTER :: lag_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_2D_lagrange_cast(lag_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_blag_project",__FILE__)
!---Initialize vectors to zero
NULLIFY(xloc)
call x%set(0.d0)
call x%get_local(xloc)
!---Operator integration loop
!$omp parallel default(firstprivate) shared(field,xloc) private(det)
allocate(j(lag_rep%nce))
!$omp do schedule(guided)
do i=1,lag_rep%mesh%nc
  call lag_rep%ncdofs(i,j) ! Get local to global DOF mapping
  !---Loop over quadrature points
  do m=1,lag_rep%quad%np
    call lag_rep%mesh%jacobian(i,lag_rep%quad%pts(:,m),sgop,vol)
    det=vol*lag_rep%quad%wts(m)
    call field%interp(i,lag_rep%quad%pts(:,m),sgop,etmp)
    !---Project on to Lagrange basis
    do jc=1,lag_rep%nce
      call oft_blag_eval(lag_rep,i,jc,lag_rep%quad%pts(:,m),rop)
      !$omp atomic
      xloc(j(jc))=xloc(j(jc)) + rop*etmp(1)*det
    end do
  end do
end do
deallocate(j)
!$omp end parallel
call x%restore_local(xloc,add=.TRUE.)
deallocate(xloc)
DEBUG_STACK_POP
END SUBROUTINE oft_blag_project
!------------------------------------------------------------------------------
!> Project a vector field onto a boundary Lagrange basis
!!
!! @note This subroutine only performs the integration of the field with
!! boundary test functions for a Lagrange basis.
!------------------------------------------------------------------------------
SUBROUTINE oft_blag_vproject(fe_rep,field,x,y,z)
class(oft_afem_type), target, intent(inout) :: fe_rep
CLASS(bfem_interp), INTENT(inout) :: field !< Vector field for projection
CLASS(oft_vector), INTENT(inout) :: x !< Field projected onto boundary Lagrange basis
CLASS(oft_vector), INTENT(inout) :: y !< Field projected onto boundary Lagrange basis
CLASS(oft_vector), INTENT(inout) :: z !< Field projected onto boundary Lagrange basis
INTEGER(i4) :: i,m,k,jc,cell,ptmap(3)
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: j
REAL(r8) :: vol,det,etmp(3),sgop(3,3),rop
REAL(r8), POINTER, DIMENSION(:) :: xloc,yloc,zloc
CLASS(oft_scalar_bfem), POINTER :: lag_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_2D_lagrange_cast(lag_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_blag_vproject",__FILE__)
!---Initialize vectors to zero
NULLIFY(xloc,yloc,zloc)
call x%set(0.d0)
call x%get_local(xloc)
call y%set(0.d0)
call y%get_local(yloc)
call z%set(0.d0)
call z%get_local(zloc)
!---Operator integration loop
!$omp parallel default(firstprivate) shared(field,xloc,yloc,zloc) private(det)
allocate(j(lag_rep%nce))
!$omp do schedule(guided)
do i=1,lag_rep%mesh%nc
  call lag_rep%ncdofs(i,j) ! Get local to global DOF mapping
  !---Loop over quadrature points
  do m=1,lag_rep%quad%np
    call lag_rep%mesh%jacobian(i,lag_rep%quad%pts(:,m),sgop,vol)
    det=vol*lag_rep%quad%wts(m)
    call field%interp(i,lag_rep%quad%pts(:,m),sgop,etmp)
    !---Project on to Lagrange basis
    do jc=1,lag_rep%nce
      call oft_blag_eval(lag_rep,i,jc,lag_rep%quad%pts(:,m),rop)
      !$omp atomic
      xloc(j(jc))=xloc(j(jc))+rop*etmp(1)*det
      !$omp atomic
      yloc(j(jc))=yloc(j(jc))+rop*etmp(2)*det
      !$omp atomic
      zloc(j(jc))=zloc(j(jc))+rop*etmp(3)*det
    end do
  end do
end do
deallocate(j)
!$omp end parallel
call x%restore_local(xloc,add=.TRUE.)
call y%restore_local(yloc,add=.TRUE.)
call z%restore_local(zloc,add=.TRUE.)
deallocate(xloc,yloc,zloc)
DEBUG_STACK_POP
END SUBROUTINE oft_blag_vproject
!------------------------------------------------------------------------------
!> Project the normal component of a vector field onto a boundary Lagrange basis
!!
!! @note This subroutine only performs the integration of the field with
!! boundary test functions for a Lagrange basis.
!------------------------------------------------------------------------------
SUBROUTINE oft_blag_nproject(fe_rep,vmesh,field,x)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(oft_mesh), target, intent(inout) :: vmesh
CLASS(fem_interp), INTENT(inout) :: field !< Vector field for projection
CLASS(oft_vector), INTENT(inout) :: x !< Field projected onto boundary Lagrange basis
INTEGER(i4) :: i,m,k,jc,cell,ptmap(3)
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: j
REAL(r8) :: vol,det,flog(4),norm(3),etmp(3),sgop(3,3),vgop(3,4),rop
REAL(r8), POINTER, DIMENSION(:) :: xloc
CLASS(oft_scalar_bfem), POINTER :: lag_rep
CLASS(oft_scalar_fem), POINTER :: vlag_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_2D_lagrange_cast(lag_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_blag_nproject",__FILE__)
!---Initialize vectors to zero
NULLIFY(xloc)
call x%set(0.d0)
call x%get_local(xloc)
!---Operator integration loop
!$omp parallel default(firstprivate) shared(field,xloc) private(det)
allocate(j(lag_rep%nce))
!$omp do schedule(guided)
do i=1,lag_rep%mesh%nc
  CALL vmesh%get_surf_map(i,cell,ptmap) ! Find parent cell and logical coordinate mapping
  call lag_rep%ncdofs(i,j) ! Get local to global DOF mapping
  !---Loop over quadrature points
  do m=1,lag_rep%quad%np
    call lag_rep%mesh%jacobian(i,lag_rep%quad%pts(:,m),sgop,vol)
    call lag_rep%mesh%norm(i,lag_rep%quad%pts(:,m),norm)
    det=vol*lag_rep%quad%wts(m)
    !---Evaluate in cell coordinates
    CALL vmesh%surf_to_vol(lag_rep%quad%pts(:,m),ptmap,flog)
    call vmesh%jacobian(cell,flog,vgop,vol)
    call field%interp(cell,flog,vgop,etmp)
    !---Project on to Lagrange basis
    do jc=1,lag_rep%nce
      call oft_blag_eval(lag_rep,i,jc,lag_rep%quad%pts(:,m),rop)
      !$omp atomic
      xloc(j(jc))=xloc(j(jc)) + rop*DOT_PRODUCT(etmp,norm)*det
    end do
  end do
end do
deallocate(j)
!$omp end parallel
call x%restore_local(xloc,add=.TRUE.)
deallocate(xloc)
DEBUG_STACK_POP
END SUBROUTINE oft_blag_nproject
end module oft_blag_operators
