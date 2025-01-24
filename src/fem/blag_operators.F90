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
USE oft_mesh_type, ONLY: oft_mesh, mesh, oft_bmesh, smesh, cell_is_curved
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
USE fem_base, ONLY: oft_fem_type, oft_bfem_type, fem_max_levels
USE fem_utils, ONLY: fem_interp, bfem_interp
USE oft_lag_basis, ONLY: oft_lagrange, oft_lagrange_level, oft_lagrange_nlevels, oft_lag_set_level, &
oft_lagrange_blevel, ML_oft_lagrange, oft_lagrange_ops, oft_lag_ops, ML_oft_lagrange_ops, &
oft_lag_eval_all, oft_lag_geval_all, oft_lag_eval, oft_blag_d2eval, &
oft_lag_nodes, oft_lagrange_lev, ML_oft_vlagrange, oft_vlagrange, oft_blagrange, oft_blag_eval, &
oft_blag_geval, oft_lagrange_minlev, oft_lag_npos, oft_scalar_fem, oft_scalar_bfem
USE oft_lag_fields, ONLY: oft_lag_create, oft_blag_create, oft_lag_vcreate
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
subroutine lag_brinterp_setup(self)
class(oft_lag_brinterp), intent(inout) :: self
!---Get local slice
CALL self%u%get_local(self%vals)
IF(.NOT.ASSOCIATED(self%lag_rep))self%lag_rep=>oft_blagrange
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
NULLIFY(self%lag_rep,self%u)
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
call smesh%hessian(cell,f,g2op,Kmat)
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
subroutine lag_bvrinterp_setup(self)
class(oft_lag_bvrinterp), intent(inout) :: self
real(r8), pointer, dimension(:) :: vtmp
!---Get local slice
IF(.NOT.ASSOCIATED(self%vals))ALLOCATE(self%vals(3,oft_lagrange%ne))
vtmp=>self%vals(1,:)
CALL self%u%get_local(vtmp,1)
vtmp=>self%vals(2,:)
CALL self%u%get_local(vtmp,2)
vtmp=>self%vals(3,:)
CALL self%u%get_local(vtmp,3)
IF(.NOT.ASSOCIATED(self%lag_rep))self%lag_rep=>oft_blagrange
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
subroutine blag_zerob(a)
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
do i=1,oft_blagrange%nbe
  j=oft_blagrange%lbe(i)
  ! vloc(j)=0.d0
  IF(oft_blagrange%global%gbe(j))vloc(j)=0.d0
end do
!---
call a%restore_local(vloc)
deallocate(vloc)
DEBUG_STACK_POP
end subroutine blag_zerob
!---------------------------------------------------------------------------
!> Zero a surface Lagrange scalar field at the mesh "grounding" node
!!
!! @note Presently the first boundary node is used as the "grounding" node
!---------------------------------------------------------------------------
subroutine blag_zerogrnd(a)
class(oft_vector), intent(inout) :: a !< Field to be zeroed
integer(i4) :: i,j
real(r8), pointer, dimension(:) :: vloc
DEBUG_STACK_PUSH
NULLIFY(vloc)
!---Cast to vector type
NULLIFY(vloc)
call a%get_local(vloc)
vloc(oft_blagrange%lbe(1))=0.d0
! !---Zero boundary values
! !$omp parallel do
! do i=1,oft_blagrange%nbe
!   vloc(oft_blagrange%lbe(i))=0.d0
! end do
call a%restore_local(vloc)
deallocate(vloc)
DEBUG_STACK_POP
end subroutine blag_zerogrnd
!---------------------------------------------------------------------------
!> Zero a surface Lagrange scalar field at all edge nodes
!---------------------------------------------------------------------------
subroutine blag_zeroe(a)
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
do i=1,oft_blagrange%ne
  j=oft_blagrange%parent%le(i)
  if(oft_lagrange%bc(j)/=3)vloc(i)=0.d0
end do
!---
call a%restore_local(vloc)
deallocate(vloc)
DEBUG_STACK_POP
end subroutine blag_zeroe
!---------------------------------------------------------------------------
!> Construct mass matrix for a boundary Lagrange scalar representation
!!
!! Supported boundary conditions
!! - \c 'none' Full matrix
!! - \c 'zerob' Dirichlet for all boundary DOF
!---------------------------------------------------------------------------
subroutine oft_blag_getmop(mat,bc)
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix object
character(LEN=*), intent(in) :: bc !< Boundary condition
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time
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
  CALL oft_blagrange%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,rop,det,mop,curved,goptmp,m,vol,jc,jr)
allocate(j(oft_blagrange%nce)) ! Local DOF and matrix indices
allocate(rop(oft_blagrange%nce)) ! Reconstructed gradient operator
allocate(mop(oft_blagrange%nce,oft_blagrange%nce)) ! Local laplacian matrix
!$omp do
do i=1,oft_blagrange%mesh%nc
  !---Get local reconstructed operators
  mop=0.d0
  do m=1,oft_blagrange%quad%np ! Loop over quadrature points
    call oft_blagrange%mesh%jacobian(i,oft_blagrange%quad%pts(:,m),goptmp,vol)
    det=vol*oft_blagrange%quad%wts(m)
    do jc=1,oft_blagrange%nce ! Loop over degrees of freedom
      call oft_blag_eval(oft_blagrange,i,jc,oft_blagrange%quad%pts(:,m),rop(jc))
    end do
    !---Compute local matrix contributions
    do jr=1,oft_blagrange%nce
      do jc=1,oft_blagrange%nce
        mop(jr,jc) = mop(jr,jc) + rop(jr)*rop(jc)*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call oft_blagrange%ncdofs(i,j)
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,mop,oft_blagrange%nce,oft_blagrange%nce)
  ! !$omp end critical
end do
deallocate(j,rop,mop)
!$omp end parallel
CALL oft_blagrange%vec_create(oft_lag_vec)
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
!! - \c 'none' Full matrix
!! - \c 'zerob' Dirichlet for all boundary DOF
!! - \c 'grnd'  Dirichlet for only groundin point
!---------------------------------------------------------------------------
subroutine oft_blag_getlop(mat,bc)
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix object
character(LEN=*), intent(in) :: bc !< Boundary condition
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time
real(r8), allocatable :: gop(:,:),lop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_lag_vec
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing Boundary LAG::LOP'
  CALL mytimer%tick()
END IF
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL oft_blagrange%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!---------------------------------------------------------------------------
! Operator integration
!---------------------------------------------------------------------------
!$omp parallel private(j,gop,det,lop,curved,goptmp,m,vol,jc,jr)
allocate(j(oft_blagrange%nce)) ! Local DOF and matrix indices
allocate(gop(3,oft_blagrange%nce)) ! Reconstructed gradient operator
allocate(lop(oft_blagrange%nce,oft_blagrange%nce)) ! Local laplacian matrix
!$omp do
do i=1,oft_blagrange%mesh%nc
  !---Get local reconstructed operators
  lop=0.d0
  do m=1,oft_blagrange%quad%np ! Loop over quadrature points
    call oft_blagrange%mesh%jacobian(i,oft_blagrange%quad%pts(:,m),goptmp,vol)
    det=vol*oft_blagrange%quad%wts(m)
    do jc=1,oft_blagrange%nce ! Loop over degrees of freedom
      call oft_blag_geval(oft_blagrange,i,jc,oft_blagrange%quad%pts(:,m),gop(:,jc),goptmp)
    end do
    !---Compute local matrix contributions
    do jr=1,oft_blagrange%nce
      do jc=1,oft_blagrange%nce
        lop(jr,jc) = lop(jr,jc) + DOT_PRODUCT(gop(:,jr),gop(:,jc))*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call oft_blagrange%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,oft_blagrange%nce
        IF(oft_blagrange%global%gbe(j(jr)))lop(jr,:)=0.d0
      END DO
    CASE("grnd")
      IF(ANY(smesh%igrnd>0))THEN
        DO jr=1,oft_blagrange%nce
          IF(ANY(j(jr)==smesh%igrnd))lop(jr,:)=0.d0
        END DO
      END IF
    CASE("edges")
      DO jr=1,oft_blagrange%nce
        jc=oft_blagrange%parent%le(j(jr))
        IF(oft_lagrange%bc(jc)/=3)lop(jr,:)=0.d0
      END DO
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,lop,oft_blagrange%nce,oft_blagrange%nce)
  ! !$omp end critical
end do
deallocate(j,gop,lop)
!$omp end parallel
ALLOCATE(lop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
lop(1,1)=1.d0
SELECT CASE(TRIM(bc))
  CASE("zerob")
    DO i=1,oft_blagrange%nbe
      IF(.NOT.oft_blagrange%linkage%leo(i))CYCLE
      jr=oft_blagrange%lbe(i)
      IF(oft_blagrange%global%gbe(jr))THEN
        j=jr
        call mat%add_values(j,j,lop,1,1)
      END IF
    END DO
  CASE("grnd")
    IF(ANY(smesh%igrnd>0))THEN
      DO i=1,oft_blagrange%nbe
        IF(.NOT.oft_blagrange%linkage%leo(i))CYCLE
        jr=oft_blagrange%lbe(i)
        IF(ANY(jr==smesh%igrnd))THEN
          j=jr
          call mat%add_values(j,j,lop,1,1)
        END IF
      END DO
    END IF
  CASE("edges")
    DO i=1,oft_blagrange%ne
      IF(oft_blagrange%be(i))CYCLE
      jr=oft_blagrange%parent%le(i)
      IF(oft_lagrange%bc(jr)/=3)THEN
        j=i
        call mat%add_values(j,j,lop,1,1)
      END IF
    END DO
    DO i=1,oft_blagrange%nbe
      IF(.NOT.oft_blagrange%linkage%leo(i))CYCLE
      jc=oft_blagrange%lbe(i)
      jr=oft_blagrange%parent%le(jc)
      IF(oft_lagrange%bc(jr)/=3)THEN
        j=jc
        call mat%add_values(j,j,lop,1,1)
      END IF
    END DO
END SELECT
DEALLOCATE(j,lop)
CALL oft_blagrange%vec_create(oft_lag_vec)
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
SUBROUTINE oft_blag_project(field,x)
CLASS(bfem_interp), INTENT(inout) :: field !< Scalar field for projection
CLASS(oft_vector), INTENT(inout) :: x !< Field projected onto boundary Lagrange basis
INTEGER(i4) :: i,m,k,jc,cell,ptmap(3)
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: j
REAL(r8) :: vol,det,etmp(1),sgop(3,3),rop
REAL(r8), POINTER, DIMENSION(:) :: xloc
DEBUG_STACK_PUSH
!---Initialize vectors to zero
NULLIFY(xloc)
call x%set(0.d0)
call x%get_local(xloc)
!---Operator integration loop
!$omp parallel default(firstprivate) shared(field,xloc) private(det)
allocate(j(oft_blagrange%nce))
!$omp do schedule(guided)
do i=1,smesh%nc
  call oft_blagrange%ncdofs(i,j) ! Get local to global DOF mapping
  !---Loop over quadrature points
  do m=1,oft_blagrange%quad%np
    call smesh%jacobian(i,oft_blagrange%quad%pts(:,m),sgop,vol)
    det=vol*oft_blagrange%quad%wts(m)
    call field%interp(i,oft_blagrange%quad%pts(:,m),sgop,etmp)
    !---Project on to Lagrange basis
    do jc=1,oft_blagrange%nce
      call oft_blag_eval(oft_blagrange,i,jc,oft_blagrange%quad%pts(:,m),rop)
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
SUBROUTINE oft_blag_vproject(field,x,y,z)
CLASS(bfem_interp), INTENT(inout) :: field !< Vector field for projection
CLASS(oft_vector), INTENT(inout) :: x !< Field projected onto boundary Lagrange basis
CLASS(oft_vector), INTENT(inout) :: y !< Field projected onto boundary Lagrange basis
CLASS(oft_vector), INTENT(inout) :: z !< Field projected onto boundary Lagrange basis
INTEGER(i4) :: i,m,k,jc,cell,ptmap(3)
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: j
REAL(r8) :: vol,det,etmp(3),sgop(3,3),rop
REAL(r8), POINTER, DIMENSION(:) :: xloc,yloc,zloc
DEBUG_STACK_PUSH
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
allocate(j(oft_blagrange%nce))
!$omp do schedule(guided)
do i=1,smesh%nc
  call oft_blagrange%ncdofs(i,j) ! Get local to global DOF mapping
  !---Loop over quadrature points
  do m=1,oft_blagrange%quad%np
    call smesh%jacobian(i,oft_blagrange%quad%pts(:,m),sgop,vol)
    det=vol*oft_blagrange%quad%wts(m)
    call field%interp(i,oft_blagrange%quad%pts(:,m),sgop,etmp)
    !---Project on to Lagrange basis
    do jc=1,oft_blagrange%nce
      call oft_blag_eval(oft_blagrange,i,jc,oft_blagrange%quad%pts(:,m),rop)
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
SUBROUTINE oft_blag_nproject(field,x)
CLASS(fem_interp), INTENT(inout) :: field !< Vector field for projection
CLASS(oft_vector), INTENT(inout) :: x !< Field projected onto boundary Lagrange basis
INTEGER(i4) :: i,m,k,jc,cell,ptmap(3)
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: j
REAL(r8) :: vol,det,flog(4),norm(3),etmp(3),sgop(3,3),vgop(3,4),rop
REAL(r8), POINTER, DIMENSION(:) :: xloc
DEBUG_STACK_PUSH
!---Initialize vectors to zero
NULLIFY(xloc)
call x%set(0.d0)
call x%get_local(xloc)
!---Operator integration loop
!$omp parallel default(firstprivate) shared(field,xloc) private(det)
allocate(j(oft_blagrange%nce))
!$omp do schedule(guided)
do i=1,smesh%nc
  CALL mesh%get_surf_map(i,cell,ptmap) ! Find parent cell and logical coordinate mapping
  call oft_blagrange%ncdofs(i,j) ! Get local to global DOF mapping
  !---Loop over quadrature points
  do m=1,oft_blagrange%quad%np
    call smesh%jacobian(i,oft_blagrange%quad%pts(:,m),sgop,vol)
    call smesh%norm(i,oft_blagrange%quad%pts(:,m),norm)
    det=vol*oft_blagrange%quad%wts(m)
    !---Evaluate in cell coordinates
    CALL mesh%surf_to_vol(oft_blagrange%quad%pts(:,m),ptmap,flog)
    call mesh%jacobian(cell,flog,vgop,vol)
    call field%interp(cell,flog,vgop,etmp)
    !---Project on to Lagrange basis
    do jc=1,oft_blagrange%nce
      call oft_blag_eval(oft_blagrange,i,jc,oft_blagrange%quad%pts(:,m),rop)
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
