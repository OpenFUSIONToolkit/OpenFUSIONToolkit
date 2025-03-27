!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!> @file oft_hcurl_operators.F90
!
!> H(Curl) FE operator definitions
!! - Operator construction
!!   - MOP: mass matrix     \f$ \int \left( u^T \cdot v \right) dV \f$
!!   - KOP: helicity matrix \f$ \int \left( \nabla \times u^T \cdot v \right) dV \f$
!!   - WOP: energy matrix   \f$ \int \left( \nabla \times u^T \cdot \nabla \times v \right) dV \f$
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
module oft_hcurl_operators
USE oft_base
USE oft_sort, ONLY: sort_array
USE oft_mesh_type, ONLY: oft_mesh, oft_bmesh, cell_is_curved
!---
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, oft_matrix, oft_matrix_ptr, &
  oft_graph_ptr, oft_graph
USE oft_deriv_matrices, ONLY: oft_diagmatrix, create_diagmatrix
USE oft_solver_base, ONLY: oft_solver, oft_solver_bc
#ifdef HAVE_ARPACK
USE oft_arpack, ONLY: oft_irlm_eigsolver
#endif
USE oft_la_utils, ONLY: create_matrix
USE oft_solver_utils, ONLY: create_mlpre, create_cg_solver, create_diag_pre, &
  create_native_pre
!---
USE fem_base, ONLY: oft_afem_type, oft_fem_type, fem_max_levels, oft_ml_fem_type, &
  oft_ml_fe_vecspace
USE fem_utils, ONLY: fem_interp
USE oft_lag_basis, ONLY: oft_lag_geval_all, oft_scalar_fem, &
  oft_3D_lagrange_cast
USE oft_lag_operators, ONLY: oft_lag_getlop, oft_lag_zerob, oft_lag_zerogrnd
USE oft_hcurl_basis, ONLY: oft_hcurl_eval_all, &
  oft_hcurl_ceval_all, oft_hcurl_get_cgops, &
  oft_hcurl_fem, oft_hcurl_bfem, oft_3D_hcurl_cast, oft_2D_hcurl_cast
IMPLICIT NONE
#include "local.h"
!------------------------------------------------------------------------------
!> Orthogonalize a H(Curl) vector against a library of given modes
!------------------------------------------------------------------------------
type, extends(oft_solver_bc) :: oft_hcurl_orthog
  integer(i4) :: nm = 0 !< Number of modes to orthogonalize against
  type(oft_vector_ptr), pointer, dimension(:,:) :: orthog => NULL() !< Library of modes
  class(oft_matrix), pointer :: wop => NULL() !< H(Curl) \f$ \nabla \times \nabla \times \f$, used as metric
  class(oft_ml_fem_type), pointer :: ML_hcurl_rep => NULL()
contains
  !> Perform orthoganlization
  procedure :: apply => hcurl_orthog_apply
  !> Clean-up internal variables
  procedure :: delete => hcurl_orthog_delete
end type oft_hcurl_orthog
!------------------------------------------------------------------------------
!> Interpolate a H(Curl) field
!------------------------------------------------------------------------------
type, extends(fem_interp) :: oft_hcurl_rinterp
  class(oft_vector), pointer :: u => NULL() !< Field to interpolate
  real(r8), pointer, dimension(:) :: vals => NULL() !< Local values
  class(oft_hcurl_fem), pointer :: hcurl_rep => NULL() !< H(Curl) FE representation
contains
  !> Retrieve local values for interpolation
  procedure :: setup => hcurl_rinterp_setup
  !> Reconstruct field
  procedure :: interp => hcurl_rinterp
  !> Destroy temporary internal storage
  procedure :: delete => hcurl_rinterp_delete
end type oft_hcurl_rinterp
!------------------------------------------------------------------------------
!> Interpolate \f$ \nabla \times \f$ of a H(Curl) field
!------------------------------------------------------------------------------
type, extends(oft_hcurl_rinterp) :: oft_hcurl_cinterp
contains
  !> Reconstruct field
  procedure :: interp => hcurl_cinterp
end type oft_hcurl_cinterp
!------------------------------------------------------------------------------
!> Clean the divergence from a H(Curl) vector field
!!
!! Divergence is removed by adding a gradient field, such that \f$ f = f +
!! \nabla \phi \f$, where \f$ \nabla^2 \phi = - \nabla \cdot f \f$. Cleaning
!! may also be applied to a field which is pre-multiplied by the H(Curl)::MOP
!! in which case \f$ \nabla^2 \phi = - \nabla^T f \f$ and \f$ f = f + M \nabla
!! \phi \f$. The mass matrix version is applied if `mop` is associated with the
!! corresponding H(Curl)::MOP.
!!
!! @note This only removes the 0-th order component, which is sufficient for
!! orthogonalization against the H(Curl)::WOP null space. Higher order cleaning
!! requires the full H(Curl) + Grad(H^1) vector space.
!------------------------------------------------------------------------------
type, extends(oft_solver_bc) :: oft_hcurl_divout
  integer(i4) :: count=0 !< Number of times apply has been called
  integer(i4) :: app_freq=1 !< Frequency to apply solver
  logical :: pm = .FALSE. !< Flag for solver convergence monitor
  logical :: internal_solver = .FALSE. !< Solver was constructed internally?
  class(oft_solver), pointer :: solver => NULL() !< Solver object for LAG::LOP operator
  class(oft_solver_bc), pointer :: bc => NULL() !< Boundary condition
  class(oft_matrix), pointer :: mop => NULL() !< Mass matrix, applies divoutm if associated
  class(oft_ml_fem_type), pointer :: ML_hcurl_rep => NULL() !< Lagrange FE space
  class(oft_ml_fem_type), pointer :: ML_lag_rep => NULL() !< Lagrange FE space
contains
  !> Setup matrix and solver with default
  procedure :: setup => hcurl_divout_setup
  !> Clean divergence from field
  procedure :: apply => hcurl_divout_apply
  !> Clean-up internal storage
  procedure :: delete => hcurl_divout_delete
end type oft_hcurl_divout
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(oft_solver_bc) :: oft_hcurl_zerob
  class(oft_ml_fem_type), pointer :: ML_hcurl_rep => NULL() !< FE representation
contains
  procedure :: apply => zerob_apply
  procedure :: delete => zerob_delete
end type oft_hcurl_zerob
!---Pre options
integer(i4) :: nu_wop(fem_max_levels)=0
real(r8) :: df_wop(fem_max_levels)=-1.d99
integer(i4) :: nu_jmlb(fem_max_levels)=20
!---Cache variables
REAL(r8), POINTER, DIMENSION(:,:) :: oft_hcurl_rop => NULL()
REAL(r8), POINTER, DIMENSION(:,:) :: oft_hcurl_cop => NULL()
!$omp threadprivate(oft_hcurl_rop,oft_hcurl_cop)
contains
!------------------------------------------------------------------------------
!> Read-in options for the basic H(Curl) ML preconditioners
!------------------------------------------------------------------------------
subroutine hcurl_mloptions(ML_hcurl_obj)
class(oft_ml_fem_type), intent(inout) :: ML_hcurl_obj
integer(i4) :: i,ierr,io_unit
namelist/hcurl_op_options/df_wop,nu_wop,nu_jmlb
DEBUG_STACK_PUSH
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,hcurl_op_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr>0)THEN
  CALL oft_abort('Error reading "hcurl_op_options" group in input file','hcurl_mloptions',__FILE__)
END IF
IF(df_wop(1)<-1.d90)THEN
  IF(oft_env%head_proc)THEN
    CALL oft_warn("No H(Curl) MG smoother settings found:")
    CALL oft_warn("  Using default values, which may result in convergence failure.")
  END IF
  DO i=ML_hcurl_obj%minlev,ML_hcurl_obj%nlevels
    CALL ML_hcurl_obj%set_level(i)
    SELECT CASE(ML_hcurl_obj%current_level%order)
    CASE(1)
      df_wop(i)=0.6d0
    CASE(2)
      df_wop(i)=0.3d0
    CASE(3)
      df_wop(i)=0.25d0
    CASE DEFAULT
      df_wop(i)=0.2d0
    END SELECT
    nu_wop(i)=MIN(64,2**(ML_hcurl_obj%nlevels-i))
  END DO
  nu_wop(ML_hcurl_obj%minlev)=64
  nu_wop(ML_hcurl_obj%nlevels)=1
END IF
DEBUG_STACK_POP
end subroutine hcurl_mloptions
!------------------------------------------------------------------------------
!> Setup interpolator for H(Curl) fields
!!
!! Fetches local representation used for interpolation from vector object
!------------------------------------------------------------------------------
subroutine hcurl_rinterp_setup(self,hcurl_rep)
class(oft_hcurl_rinterp), intent(inout) :: self
class(oft_afem_type), target, intent(inout) :: hcurl_rep
!---
SELECT TYPE(hcurl_rep)
CLASS IS(oft_hcurl_fem)
  self%hcurl_rep=>hcurl_rep
CLASS DEFAULT
  CALL oft_abort("Incorrect FE type","hcurl_rinterp_setup",__FILE__)
END SELECT
self%mesh=>self%hcurl_rep%mesh
!---Get local slice
CALL self%u%get_local(self%vals)
end subroutine hcurl_rinterp_setup
!------------------------------------------------------------------------------
!> Destroy temporary internal storage
!------------------------------------------------------------------------------
subroutine hcurl_rinterp_delete(self)
class(oft_hcurl_rinterp), intent(inout) :: self
!---Deallocate local storage
IF(ASSOCIATED(self%vals))DEALLOCATE(self%vals)
NULLIFY(self%hcurl_rep,self%u)
end subroutine hcurl_rinterp_delete
!------------------------------------------------------------------------------
!> Reconstruct a H(Curl) vector field
!------------------------------------------------------------------------------
subroutine hcurl_rinterp(self,cell,f,gop,val)
class(oft_hcurl_rinterp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [4]
real(r8), intent(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
real(r8), intent(out) :: val(:) !< Reconstructed field at f [3]
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8), allocatable  :: rop(:,:)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!','hcurl_rinterp',__FILE__)
!---Get dofs
allocate(j(self%hcurl_rep%nce),rop(3,self%hcurl_rep%nce))
call self%hcurl_rep%ncdofs(cell,j) ! get DOFs
!---Reconstruct field
call oft_hcurl_eval_all(self%hcurl_rep,cell,f,rop,gop)
val=0.d0
do jc=1,self%hcurl_rep%nce
  val=val+self%vals(j(jc))*rop(:,jc)
end do
deallocate(j,rop)
DEBUG_STACK_POP
end subroutine hcurl_rinterp
!------------------------------------------------------------------------------
!> Reconstruct the curl of a H(Curl) vector field
!------------------------------------------------------------------------------
subroutine hcurl_cinterp(self,cell,f,gop,val)
class(oft_hcurl_cinterp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [4]
real(r8), intent(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
real(r8), intent(out) :: val(:) !< Reconstructed field at f [3]
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: cgop(3,6)
real(r8), allocatable  :: cop(:,:)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!','hcurl_cinterp',__FILE__)
!---Get dofs
allocate(j(self%hcurl_rep%nce),cop(3,self%hcurl_rep%nce))
call self%hcurl_rep%ncdofs(cell,j) ! get DOFs
!---Reconstruct field
call oft_hcurl_get_cgops(gop,cgop)
call oft_hcurl_ceval_all(self%hcurl_rep,cell,f,cop,cgop)
val=0.d0
do jc=1,self%hcurl_rep%nce
  val=val+self%vals(j(jc))*cop(:,jc)
end do
deallocate(j,cop)
DEBUG_STACK_POP
end subroutine hcurl_cinterp
!------------------------------------------------------------------------------
!> Apply the divergence operator to a H(Curl) field
!!
!! @note This subroutine computes the divergence of a H(Curl) field as
!! projected on to a linear Lagrange scalar basis.
!------------------------------------------------------------------------------
subroutine hcurl_div(hcurl_fe,lag_fe,a,b)
class(oft_afem_type), intent(inout) :: hcurl_fe
class(oft_afem_type), intent(inout) :: lag_fe
class(oft_vector), intent(inout) :: a !< Needs docs
class(oft_vector), intent(inout) :: b !< Needs docs
real(r8), pointer, dimension(:) :: aloc
real(r8), pointer, dimension(:) :: bloc
integer(i4) :: i,jr,jc,m
integer(i4), allocatable :: j_curl(:),j_grad(:)
real(r8) :: goptmp(3,4)
real(r8) :: v,f(4),det,vol
logical :: curved
real(r8), allocatable :: rop_curl(:,:),rop_grad(:,:),btmp(:)
CLASS(oft_scalar_fem), POINTER :: lag_rep
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_lagrange_cast(lag_rep,lag_fe))CALL oft_abort("Incorrect Lagrange FE type","hcurl_div",__FILE__)
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,hcurl_fe))CALL oft_abort("Incorrect HCurl FE type","hcurl_div",__FILE__)
!------------------------------------------------------------------------------
! Allocate Laplacian Op
!------------------------------------------------------------------------------
NULLIFY(aloc,bloc)
!---Get local values
CALL a%get_local(aloc)
!---Zero output and get aliases
CALL b%set(0.d0)
CALL b%get_local(bloc)
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j_curl,j_grad,rop_curl,rop_grad,det,btmp,curved,goptmp,m,v,jc,jr,f)
allocate(j_curl(hcurl_rep%nce)) ! Local DOF and matrix indices
allocate(j_grad(lag_rep%nce)) ! Local DOF and matrix indices
allocate(rop_grad(3,lag_rep%nce)) ! Reconstructed gradient operator
allocate(rop_curl(3,hcurl_rep%nce)) ! Reconstructed gradient operator
allocate(btmp(lag_rep%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,hcurl_rep%mesh%nc
  !---Get local to global DOF mapping
  call hcurl_rep%ncdofs(i,j_curl)
  call lag_rep%ncdofs(i,j_grad)
  curved=cell_is_curved(hcurl_rep%mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  btmp=0.d0
  do m=1,hcurl_rep%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call hcurl_rep%mesh%jacobian(i,hcurl_rep%quad%pts(:,m),goptmp,v)
    det=v*hcurl_rep%quad%wts(m)
    call oft_hcurl_eval_all(hcurl_rep,i,hcurl_rep%quad%pts(:,m),rop_curl,goptmp)
    call oft_lag_geval_all(lag_rep,i,hcurl_rep%quad%pts(:,m),rop_grad,goptmp)
    !---Compute local operator contribution
    do jr=1,lag_rep%nce
      do jc=1,hcurl_rep%nce
        btmp(jr)=btmp(jr)+DOT_PRODUCT(rop_grad(:,jr),rop_curl(:,jc))*aloc(j_curl(jc))*det
      end do
    end do
  end do
  !---Add local values to global vector
  do jr=1,lag_rep%nce
    !$omp atomic
    bloc(j_grad(jr))=bloc(j_grad(jr))+btmp(jr)
  end do
end do
deallocate(j_curl,j_grad)
deallocate(rop_curl,rop_grad,btmp)
!$omp end parallel
CALL b%restore_local(bloc,add=.TRUE.)
deallocate(aloc,bloc)
DEBUG_STACK_POP
end subroutine hcurl_div
!------------------------------------------------------------------------------
!> Add the gradient of a linear Lagrange scalar field to a H(Curl) field
!------------------------------------------------------------------------------
subroutine hcurl_grad(hcurl_fe,a,b)
class(oft_afem_type), intent(inout) :: hcurl_fe
class(oft_vector), intent(inout) :: a !< Needs docs
class(oft_vector), intent(inout) :: b !< Needs docs
real(r8), pointer, dimension(:) :: aloc,bloc
integer(i4), allocatable :: emap(:)
integer :: i,j,k,l
real(r8) :: reg
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,hcurl_fe))CALL oft_abort("Incorrect HCurl FE type","hcurl_grad",__FILE__)
!---Get local values
NULLIFY(aloc,bloc)
CALL a%get_local(aloc)
CALL b%get_local(bloc)
!---Get boundary map
do i=1,hcurl_rep%mesh%ne
  bloc((i-1)*hcurl_rep%gstruct(2)+1)=bloc((i-1)*hcurl_rep%gstruct(2)+1) + &
    (aloc(hcurl_rep%mesh%le(2,i))-aloc(hcurl_rep%mesh%le(1,i)))*SIGN(1_i8,hcurl_rep%mesh%global%le(i))
enddo
!---
CALL b%restore_local(bloc)
DEALLOCATE(aloc,bloc)
DEBUG_STACK_POP
end subroutine hcurl_grad
!------------------------------------------------------------------------------
!> Apply the transposed gradient operator to a H(Curl) vector field.
!!
!! @note Only the 0th order component is computed from the discrete gradient
!! operator
!------------------------------------------------------------------------------
subroutine hcurl_gradtp(hcurl_fe,a,b)
class(oft_afem_type), intent(inout) :: hcurl_fe
class(oft_vector), intent(inout) :: a !< Input field
class(oft_vector), intent(inout) :: b !< \f$ G^{T} a \f$
real(r8), pointer, dimension(:) :: aloc,bloc
integer(i4), allocatable, dimension(:) :: emap
integer :: i,j,k
class(oft_mesh), pointer :: mesh
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,hcurl_fe))CALL oft_abort("Incorrect HCurl FE type","hcurl_gradtp",__FILE__)
!---Get local values
NULLIFY(aloc,bloc)
CALL a%get_local(aloc)
CALL b%get_local(bloc)
!---
mesh=>hcurl_rep%mesh
allocate(emap(mesh%ne))
CALL get_inverse_map(mesh%lbe,mesh%nbe,emap,mesh%ne)
!$omp parallel do private(j,k)
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
      bloc(i)=bloc(i)+aloc(k)*(1-2*min(1,abs(mesh%le(1,k)-i)))
    else
      bloc(i)=bloc(i)-aloc(k)*(1-2*min(1,abs(mesh%le(1,k)-i)))
    end if
  enddo
enddo
deallocate(emap)
CALL b%restore_local(bloc,add=.TRUE.)
DEALLOCATE(aloc,bloc)
DEBUG_STACK_POP
end subroutine hcurl_gradtp
!------------------------------------------------------------------------------
!> Zero the tangential component of a H(Curl) vector field on the boundary
!------------------------------------------------------------------------------
subroutine zerob_apply(self,a)
class(oft_hcurl_zerob), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to be zeroed
real(r8), pointer, dimension(:) :: aloc
integer(i4) :: i,j
DEBUG_STACK_PUSH
! Cast to vector type
NULLIFY(aloc)
CALL a%get_local(aloc)
! Apply operator
do i=1,self%ML_hcurl_rep%current_level%nbe
  j=self%ML_hcurl_rep%current_level%lbe(i)
  if(self%ML_hcurl_rep%current_level%global%gbe(j))aloc(j)=0.d0
end do
CALL a%restore_local(aloc)
DEALLOCATE(aloc)
DEBUG_STACK_POP
end subroutine zerob_apply
!------------------------------------------------------------------------------
!> Zero the tangential component of a H(Curl) vector field on the boundary
!------------------------------------------------------------------------------
subroutine zerob_delete(self)
class(oft_hcurl_zerob), intent(inout) :: self
NULLIFY(self%ML_hcurl_rep)
end subroutine zerob_delete
!------------------------------------------------------------------------------
!> Construct mass matrix for a H(Curl) representation
!!
!! Supported boundary conditions
!! - `'none'` Full matrix
!! - `'zerob'` Dirichlet for all boundary DOF
!------------------------------------------------------------------------------
subroutine oft_hcurl_getmop(fe_rep,mat,bc)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix object
character(LEN=*), intent(in) :: bc !< Boundary condition
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time
real(r8), allocatable :: rop(:,:),mop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_hcurl_vec
type(oft_timer) :: mytimer
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing H(Curl)::MOP'
  CALL mytimer%tick()
END IF
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_hcurl_getmop",__FILE__)
!------------------------------------------------------------------------------
! Allocate matrix
!------------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL hcurl_rep%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,rop,det,mop,curved,goptmp,m,vol,jc,jr)
allocate(j(hcurl_rep%nce)) ! Local DOF and matrix indices
allocate(rop(3,hcurl_rep%nce)) ! Reconstructed gradient operator
allocate(mop(hcurl_rep%nce,hcurl_rep%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,hcurl_rep%mesh%nc
  curved=cell_is_curved(hcurl_rep%mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  mop=0.d0
  do m=1,hcurl_rep%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call hcurl_rep%mesh%jacobian(i,hcurl_rep%quad%pts(:,m),goptmp,vol)
    det=vol*hcurl_rep%quad%wts(m)
    call oft_hcurl_eval_all(hcurl_rep,i,hcurl_rep%quad%pts(:,m),rop,goptmp)
    !---Compute local matrix contributions
    do jr=1,hcurl_rep%nce
      do jc=1,hcurl_rep%nce
        mop(jr,jc) = mop(jr,jc) + DOT_PRODUCT(rop(:,jr),rop(:,jc))*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call hcurl_rep%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,hcurl_rep%nce
        IF(hcurl_rep%global%gbe(j(jr)))THEN
          mop(jr,:)=0.d0
        END IF
      END DO
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,mop,hcurl_rep%nce,hcurl_rep%nce)
  ! !$omp end critical
end do
deallocate(j,rop,mop)
!$omp end parallel
ALLOCATE(mop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
SELECT CASE(TRIM(bc))
  CASE("zerob")
    mop(1,1)=1.d0
    DO i=1,hcurl_rep%nbe
      jr=hcurl_rep%lbe(i)
      IF(.NOT.hcurl_rep%global%gbe(jr))CYCLE
      IF(.NOT.hcurl_rep%linkage%leo(i))CYCLE
      j=jr
      call mat%add_values(j,j,mop,1,1)
    END DO
END SELECT
DEALLOCATE(j,mop)
CALL hcurl_rep%vec_create(oft_hcurl_vec)
CALL mat%assemble(oft_hcurl_vec)
CALL oft_hcurl_vec%delete
DEALLOCATE(oft_hcurl_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_hcurl_getmop
!------------------------------------------------------------------------------
!> Construct helicity matrix for a H(Curl) representation
!!
!! Supported boundary conditions
!! - `'none'` Full matrix
!! - `'zerob'` Dirichlet for all boundary DOF
!------------------------------------------------------------------------------
subroutine oft_hcurl_getkop(fe_rep,mat,bc)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix object
character(LEN=*), intent(in) :: bc !< Boundary condition
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),cgop(3,6),elapsed_time
real(r8), allocatable :: rop(:,:),cop(:,:),kop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_hcurl_vec
type(oft_timer) :: mytimer
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing H(Curl)::KOP'
  CALL mytimer%tick()
END IF
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_hcurl_getkop",__FILE__)
!------------------------------------------------------------------------------
! Allocate matrix
!------------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL hcurl_rep%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,rop,cop,det,kop,curved,goptmp,cgop,m,vol,jc,jr)
allocate(j(hcurl_rep%nce)) ! Local DOF and matrix indices
allocate(rop(3,hcurl_rep%nce)) ! Reconstructed gradient operator
allocate(cop(3,hcurl_rep%nce)) ! Reconstructed gradient operator
allocate(kop(hcurl_rep%nce,hcurl_rep%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,hcurl_rep%mesh%nc
  curved=cell_is_curved(hcurl_rep%mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  kop=0.d0
  do m=1,hcurl_rep%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)then
      call hcurl_rep%mesh%jacobian(i,hcurl_rep%quad%pts(:,m),goptmp,vol)
      call oft_hcurl_get_cgops(goptmp,cgop)
    end if
    det=vol*hcurl_rep%quad%wts(m)
    call oft_hcurl_eval_all(hcurl_rep,i,hcurl_rep%quad%pts(:,m),rop,goptmp)
    call oft_hcurl_ceval_all(hcurl_rep,i,hcurl_rep%quad%pts(:,m),cop,cgop)
    !---Compute local matrix contributions
    do jr=1,hcurl_rep%nce
      do jc=1,hcurl_rep%nce
        kop(jr,jc) = kop(jr,jc) + DOT_PRODUCT(rop(:,jr),cop(:,jc))*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call hcurl_rep%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,hcurl_rep%nce
        IF(hcurl_rep%global%gbe(j(jr)))kop(jr,:)=0.d0
      END DO
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,kop,hcurl_rep%nce,hcurl_rep%nce)
  ! !$omp end critical
end do
deallocate(j,rop,cop,kop)
!$omp end parallel
ALLOCATE(kop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
SELECT CASE(TRIM(bc))
  CASE("zerob")
    kop(1,1)=1.d0
    DO i=1,hcurl_rep%nbe
      jr=hcurl_rep%lbe(i)
      IF(.NOT.hcurl_rep%global%gbe(jr))CYCLE
      IF(.NOT.hcurl_rep%linkage%leo(i))CYCLE
      j=jr
      call mat%add_values(j,j,kop,1,1)
    END DO
END SELECT
DEALLOCATE(j,kop)
CALL hcurl_rep%vec_create(oft_hcurl_vec)
CALL mat%assemble(oft_hcurl_vec)
CALL oft_hcurl_vec%delete
DEALLOCATE(oft_hcurl_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_hcurl_getkop
!------------------------------------------------------------------------------
!> Construct energy matrix for a H(Curl) representation
!!
!! Supported boundary conditions
!! - `'none'` Full matrix
!! - `'zerob'` Dirichlet for all boundary DOF
!------------------------------------------------------------------------------
subroutine oft_hcurl_getwop(fe_rep,mat,bc)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix object
character(LEN=*), intent(in) :: bc !< Boundary condition
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),cgop(3,6),elapsed_time
real(r8), allocatable :: cop(:,:),wop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_hcurl_vec
type(oft_timer) :: mytimer
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing H(Curl)::WOP'
  CALL mytimer%tick()
END IF
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_hcurl_getwop",__FILE__)
!------------------------------------------------------------------------------
! Allocate Laplacian Op
!------------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL hcurl_rep%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,cop,det,wop,curved,goptmp,cgop,m,vol,jc,jr)
allocate(j(hcurl_rep%nce)) ! Local DOF and matrix indices
allocate(cop(3,hcurl_rep%nce)) ! Reconstructed gradient operator
allocate(wop(hcurl_rep%nce,hcurl_rep%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,hcurl_rep%mesh%nc
  curved=cell_is_curved(hcurl_rep%mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  wop=0.d0
  do m=1,hcurl_rep%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)then
      call hcurl_rep%mesh%jacobian(i,hcurl_rep%quad%pts(:,m),goptmp,vol)
      call oft_hcurl_get_cgops(goptmp,cgop)
    end if
    det=vol*hcurl_rep%quad%wts(m)
    call oft_hcurl_ceval_all(hcurl_rep,i,hcurl_rep%quad%pts(:,m),cop,cgop)
    !---Compute local matrix contributions
    do jr=1,hcurl_rep%nce
      do jc=1,hcurl_rep%nce
        wop(jr,jc) = wop(jr,jc) + DOT_PRODUCT(cop(:,jr),cop(:,jc))*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call hcurl_rep%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,hcurl_rep%nce
        IF(hcurl_rep%global%gbe(j(jr)))wop(jr,:)=0.d0
      END DO
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,wop,hcurl_rep%nce,hcurl_rep%nce)
  ! !$omp end critical
end do
deallocate(j,cop,wop)
!$omp end parallel
ALLOCATE(wop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
SELECT CASE(TRIM(bc))
  CASE("zerob")
    wop(1,1)=1.d0
    DO i=1,hcurl_rep%nbe
      jr=hcurl_rep%lbe(i)
      IF(.NOT.hcurl_rep%global%gbe(jr))CYCLE
      IF(.NOT.hcurl_rep%linkage%leo(i))CYCLE
      j=jr
      call mat%add_values(j,j,wop,1,1)
    END DO
END SELECT
DEALLOCATE(j,wop)
CALL hcurl_rep%vec_create(oft_hcurl_vec)
CALL mat%assemble(oft_hcurl_vec)
CALL oft_hcurl_vec%delete
DEALLOCATE(oft_hcurl_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_hcurl_getwop
!------------------------------------------------------------------------------
!> Construct force-free response matrix for a H(Curl) representation
!!
!! Supported boundary conditions
!! - `'none'` Full matrix
!! - `'zerob'` Dirichlet for all boundary DOF
!------------------------------------------------------------------------------
subroutine oft_hcurl_getjmlb(fe_rep,mat,alam,bc)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix object
real(r8), intent(in) :: alam !< Lambda for response
character(LEN=*), intent(in) :: bc !< Boundary condition
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),cgop(3,6),elapsed_time
real(r8), allocatable :: rop(:,:),cop(:,:),wop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_hcurl_vec
type(oft_timer) :: mytimer
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing H(Curl)::JMLB'
  CALL mytimer%tick()
END IF
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_hcurl_getjmlb",__FILE__)
!------------------------------------------------------------------------------
! Allocate Laplacian Op
!------------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL hcurl_rep%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,rop,cop,det,wop,curved,goptmp,cgop,m,vol,jc,jr)
allocate(j(hcurl_rep%nce)) ! Local DOF and matrix indices
allocate(rop(3,hcurl_rep%nce)) ! Reconstructed gradient operator
allocate(cop(3,hcurl_rep%nce)) ! Reconstructed gradient operator
allocate(wop(hcurl_rep%nce,hcurl_rep%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,hcurl_rep%mesh%nc
  curved=cell_is_curved(hcurl_rep%mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  wop=0.d0
  do m=1,hcurl_rep%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)then
      call hcurl_rep%mesh%jacobian(i,hcurl_rep%quad%pts(:,m),goptmp,vol)
      call oft_hcurl_get_cgops(goptmp,cgop)
    end if
    det=vol*hcurl_rep%quad%wts(m)
    call oft_hcurl_eval_all(hcurl_rep,i,hcurl_rep%quad%pts(:,m),rop,goptmp)
    call oft_hcurl_ceval_all(hcurl_rep,i,hcurl_rep%quad%pts(:,m),cop,cgop)
    !---Compute local matrix contributions
    do jr=1,hcurl_rep%nce
      do jc=1,hcurl_rep%nce
        wop(jr,jc) = wop(jr,jc) + (DOT_PRODUCT(cop(:,jr),cop(:,jc))-alam*DOT_PRODUCT(rop(:,jr),cop(:,jc)))*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call hcurl_rep%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,hcurl_rep%nce
        IF(hcurl_rep%global%gbe(j(jr)))wop(jr,:)=0.d0
      END DO
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,wop,hcurl_rep%nce,hcurl_rep%nce)
  ! !$omp end critical
end do
deallocate(j,rop,cop,wop)
!$omp end parallel
ALLOCATE(wop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
SELECT CASE(TRIM(bc))
  CASE("zerob")
    wop(1,1)=1.d0
    DO i=1,hcurl_rep%nbe
      jr=hcurl_rep%lbe(i)
      IF(.NOT.hcurl_rep%global%gbe(jr))CYCLE
      IF(.NOT.hcurl_rep%linkage%leo(i))CYCLE
      j=jr
      call mat%add_values(j,j,wop,1,1)
    END DO
END SELECT
DEALLOCATE(j,wop)
CALL hcurl_rep%vec_create(oft_hcurl_vec)
CALL mat%assemble(oft_hcurl_vec)
CALL oft_hcurl_vec%delete
DEALLOCATE(oft_hcurl_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_hcurl_getjmlb
!------------------------------------------------------------------------------
!> Project a vector field onto a H(Curl) basis
!!
!! @note This subroutine only performs the integration of the field with
!! test functions for a H(Curl) basis. To retrieve the correct projection the
!! result must be multiplied by the inverse of H(Curl)::MOP
!------------------------------------------------------------------------------
subroutine oft_hcurl_project(fe_rep,field,x)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(fem_interp), intent(inout) :: field !< Vector field for projection
class(oft_vector), intent(inout) :: x !< Field projected onto H(Curl) basis
!---
integer(i4) :: i,jc,m
integer(i4), allocatable, dimension(:) :: j
real(r8) :: det,vol,bcc(3),goptmp(3,4)
real(r8), pointer, dimension(:) :: xloc
real(r8), allocatable, dimension(:,:) :: rop
logical :: curved
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_hcurl_project",__FILE__)
!---Initialize vectors to zero
NULLIFY(xloc)
call x%set(0.d0)
call x%get_local(xloc)
!---Integerate over the volume
!$omp parallel default(firstprivate) shared(xloc) private(curved,det)
allocate(j(hcurl_rep%nce),rop(3,hcurl_rep%nce))
!$omp do schedule(guided)
do i=1,hcurl_rep%mesh%nc ! Loop over cells
  call hcurl_rep%ncdofs(i,j) ! Get DOFs
  curved=cell_is_curved(hcurl_rep%mesh,i) ! Straight cell test
  do m=1,hcurl_rep%quad%np
    if(curved.OR.m==1)call hcurl_rep%mesh%jacobian(i,hcurl_rep%quad%pts(:,m),goptmp,vol)
    det=vol*hcurl_rep%quad%wts(m)
    call field%interp(i,hcurl_rep%quad%pts(:,m),goptmp,bcc)
    call oft_hcurl_eval_all(hcurl_rep,i,hcurl_rep%quad%pts(:,m),rop,goptmp)
    do jc=1,hcurl_rep%nce
      !$omp atomic
      xloc(j(jc))=xloc(j(jc))+DOT_PRODUCT(rop(:,jc),bcc)*det
    end do
  end do
end do
deallocate(j,rop)
!$omp end parallel
call x%restore_local(xloc,add=.TRUE.)
deallocate(xloc)
DEBUG_STACK_POP
end subroutine oft_hcurl_project
!---------------------------------------------------------------------------------
!> Compute the boundary term for integration by parts of the curl operator using
!! a HCurl basis
!---------------------------------------------------------------------------------
SUBROUTINE oft_hcurl_bcurl(fe_rep,bfe_rep,field,x)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(oft_afem_type), target, intent(inout) :: bfe_rep
CLASS(fem_interp), INTENT(inout) :: field !< Vector field for projection
CLASS(oft_vector), INTENT(inout) :: x
!< \f$ \int \left( \textbf{u}^T \times \textbf{F} \right) \cdot \textbf{dS} \f$
INTEGER(i4) :: i,m,k,jc,cell,ptmap(3)
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: j
REAL(r8) :: vol,det,flog(4),norm(3),etmp(3),sgop(3,3),vgop(3,4)
REAL(r8), POINTER, DIMENSION(:) :: xloc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: rop
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
CLASS(oft_hcurl_bfem), POINTER :: bhcurl_rep
CLASS(oft_mesh), POINTER :: mesh
CLASS(oft_bmesh), POINTER :: smesh
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_hcurl_cast(hcurl_rep,fe_rep))CALL oft_abort("Incorrect 3D FE type","oft_hcurl_project",__FILE__)
IF(.NOT.oft_2D_hcurl_cast(bhcurl_rep,bfe_rep))CALL oft_abort("Incorrect 2D FE type","oft_hcurl_project",__FILE__)
mesh=>hcurl_rep%mesh
smesh=>bhcurl_rep%mesh
!---Initialize vectors to zero
NULLIFY(xloc)
call x%set(0.d0)
call x%get_local(xloc)
!---Operator integration loop
det=0.d0
!$omp parallel default(firstprivate) shared(field,xloc) private(det)
allocate(j(hcurl_rep%nce),rop(3,hcurl_rep%nce))
!$omp do schedule(guided)
do i=1,smesh%nc
  CALL mesh%get_surf_map(i,cell,ptmap) ! Find parent cell and logical coordinate mapping
  call hcurl_rep%ncdofs(cell,j) ! Get local to global DOF mapping
  !---Loop over quadrature points
  do m=1,bhcurl_rep%quad%np
    call smesh%jacobian(i,bhcurl_rep%quad%pts(:,m),sgop,vol)
    call smesh%norm(i,bhcurl_rep%quad%pts(:,m),norm)
    det=vol*bhcurl_rep%quad%wts(m)
    !---Evaluate in cell coordinates
    CALL mesh%surf_to_vol(bhcurl_rep%quad%pts(:,m),ptmap,flog)
    call mesh%jacobian(cell,flog,vgop,vol)
    call field%interp(cell,flog,vgop,etmp)
    !---Project on to HCurl basis
    call oft_hcurl_eval_all(hcurl_rep,cell,flog,rop,vgop)
    do jc=1,hcurl_rep%nce
      !$omp atomic
      xloc(j(jc))=xloc(j(jc)) + DOT_PRODUCT(cross_product(rop(:,jc),etmp),norm)*det
    end do
  end do
end do
deallocate(j,rop)
!$omp end parallel
call x%restore_local(xloc,add=.TRUE.)
deallocate(xloc)
DEBUG_STACK_POP
END SUBROUTINE oft_hcurl_bcurl
!------------------------------------------------------------------------------
!> Setup matrix and solver with default
!------------------------------------------------------------------------------
subroutine hcurl_divout_setup(self,ML_hcurl_rep,ML_lag_rep,bc,solver)
class(oft_hcurl_divout), intent(inout) :: self
class(oft_ml_fem_type), target, intent(inout) :: ML_hcurl_rep
CLASS(oft_ml_fem_type), target, intent(inout) :: ML_lag_rep
character(LEN=*), intent(in) :: bc !< Boundary condition
class(oft_solver), target, optional, intent(in) :: solver
TYPE(oft_lag_zerob), POINTER :: bc_zerob
TYPE(oft_lag_zerogrnd), POINTER :: bc_zerogrnd
DEBUG_STACK_PUSH
self%ML_hcurl_rep=>ML_hcurl_rep
self%ML_lag_rep=>ML_lag_rep
IF(PRESENT(solver))THEN
  self%solver=>solver
  self%internal_solver=.FALSE.
ELSE
  CALL create_cg_solver(self%solver)
  self%solver%its=-3
  CALL oft_lag_getlop(ML_lag_rep%current_level,self%solver%A,bc)
  CALL create_diag_pre(self%solver%pre)
  self%internal_solver=.TRUE.
END IF
!
SELECT CASE(TRIM(bc)) 
CASE('grnd')
  ALLOCATE(bc_zerogrnd)
  bc_zerogrnd%ML_lag_rep=>ML_lag_rep
  self%bc=>bc_zerogrnd
CASE('zero')
  ALLOCATE(bc_zerob)
  bc_zerob%ML_lag_rep=>ML_lag_rep
  self%bc=>bc_zerob
CASE('none')
  NULLIFY(self%bc)
CASE DEFAULT
  CALL oft_abort("Invalid BC","hcurl_divout_setup",__FILE__)
END SELECT
DEBUG_STACK_POP
end subroutine hcurl_divout_setup
!------------------------------------------------------------------------------
!> Remove divergence from a H(Curl) vector field by adding a gradient correction
!------------------------------------------------------------------------------
subroutine hcurl_divout_apply(self,a)
class(oft_hcurl_divout), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field for divergence cleaning
class(oft_vector), pointer :: u,g,tmp,tmp2
integer(i4) :: i,order_tmp
real(r8) :: uu
logical :: pm_save
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%solver))CALL oft_abort('No solver specified.','hcurl_divout_apply',__FILE__)
self%count=self%count+1
IF(mod(self%count,self%app_freq)/=0)THEN
  DEBUG_STACK_POP
  RETURN
END IF
!---
call self%ML_lag_rep%vec_create(g)
call self%ML_lag_rep%vec_create(u)
!---
IF(ASSOCIATED(self%mop))THEN
  CALL hcurl_gradtp(self%ML_hcurl_rep%current_level,a,g)
ELSE
  CALL hcurl_div(self%ML_hcurl_rep%current_level,self%ML_lag_rep%current_level,a,g)
END IF
uu=a%dot(a)
self%solver%atol=MAX(self%solver%atol,SQRT(uu*1.d-20))
call u%set(0.d0)
call self%bc%apply(g)
!---
pm_save=oft_env%pm; oft_env%pm=self%pm
call self%solver%apply(u,g)
oft_env%pm=pm_save
!---
CALL u%scale(-1.d0)
CALL a%new(tmp)
CALL hcurl_grad(self%ML_hcurl_rep%current_level,u,tmp)
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
deallocate(tmp,u,g)
DEBUG_STACK_POP
end subroutine hcurl_divout_apply
!------------------------------------------------------------------------------
!> Clean-up internal storage for a oft_hcurl_divout object
!------------------------------------------------------------------------------
subroutine hcurl_divout_delete(self)
class(oft_hcurl_divout), intent(inout) :: self
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
NULLIFY(self%ML_hcurl_rep,self%ML_lag_rep)
end subroutine hcurl_divout_delete
!------------------------------------------------------------------------------
!> Orthogonalize a H(Curl) vector against a library of given modes.
!!
!! @note Used as a member function of oft_hcurl_orthog only
!------------------------------------------------------------------------------
subroutine hcurl_orthog_apply(self,a)
class(oft_hcurl_orthog), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to orthogonalize
class(oft_vector), pointer :: b
real(r8) :: c
integer(i4) :: i
DEBUG_STACK_PUSH
!---Get temporary variable
call a%new(b)
do i=1,self%nm
  !---Compute coupling
  call self%wop%apply(self%orthog(i,self%ML_hcurl_rep%level)%f,b)
  c=a%dot(b)
  !---Remove coupling
  call a%add(1.d0,-c,self%orthog(i,self%ML_hcurl_rep%level)%f)
end do
!---Delete temporary variable
call b%delete
deallocate(b)
DEBUG_STACK_POP
end subroutine hcurl_orthog_apply
!------------------------------------------------------------------------------
!> Clean-up internal storage for a oft_hcurl_orthog object
!!
!! @note Used as a member function of oft_hcurl_orthog only
!------------------------------------------------------------------------------
subroutine hcurl_orthog_delete(self)
class(oft_hcurl_orthog), intent(inout) :: self
end subroutine hcurl_orthog_delete
!------------------------------------------------------------------------------
!> Construct interpolation matrices on each MG level
!------------------------------------------------------------------------------
SUBROUTINE hcurl_setup_interp(ML_hcurl_rep)
CLASS(oft_ml_fem_type), intent(inout) :: ML_hcurl_rep
INTEGER(i4) :: i
DEBUG_STACK_PUSH
!---
DO i=ML_hcurl_rep%minlev+1,ML_hcurl_rep%nlevels
  CALL ML_hcurl_rep%set_level(i)
  !---
  if(ML_hcurl_rep%level==ML_hcurl_rep%blevel+1)then
    CYCLE
  end if
  !---Setup interpolation
  if(ML_hcurl_rep%current_level%order==1)then
    CALL hcurl_ginterpmatrix(ML_hcurl_rep%interp_matrices(ML_hcurl_rep%level)%m)
    CALL ML_hcurl_rep%interp_matrices(ML_hcurl_rep%level)%m%assemble
  else
    CALL hcurl_pinterpmatrix(ML_hcurl_rep%interp_matrices(ML_hcurl_rep%level)%m)
    CALL ML_hcurl_rep%interp_matrices(ML_hcurl_rep%level)%m%assemble
  end if
END DO
DEBUG_STACK_POP
CONTAINS
!------------------------------------------------------------------------------
!> Construct interpolation matrix for polynomial levels
!------------------------------------------------------------------------------
SUBROUTINE hcurl_ginterpmatrix(mat)
class(oft_matrix), pointer, intent(inout) :: mat !< Interpolation matrix
INTEGER(i4) :: i,j,k,m,icors,ifine,jb,i_ind(1),j_ind(1)
INTEGER(i4) :: etmp(2),ftmp(3),fetmp(3),ctmp(4),fc,ed
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pmap,emap,fmap
CLASS(oft_hcurl_fem), POINTER :: hcurl_cors => NULL()
CLASS(oft_hcurl_fem), POINTER :: hcurl_fine => NULL()
class(oft_mesh), pointer :: cmesh
CLASS(oft_vector), POINTER :: hcurl_vec_fine,hcurl_vec_cors
integer(i4) :: jfe(3),jce(6)
integer(i4), pointer :: lcdg(:),lfde(:,:),lede(:,:),lcde(:,:)
real(r8) :: f(4),incr,val,d(3),goptmp(3,4),v,mop(1),h_rop(3,6)
type(oft_graph_ptr), pointer :: graphs(:,:)
type(oft_graph), POINTER :: interp_graph
CLASS(oft_mesh), POINTER :: mesh
DEBUG_STACK_PUSH
!---
if(ML_hcurl_rep%ml_mesh%level<1)then
  call oft_abort('Invalid mesh level','hcurl_ginterpmatrix',__FILE__)
end if
mesh=>ML_hcurl_rep%ml_mesh%mesh
cmesh=>ML_hcurl_rep%ml_mesh%meshes(ML_hcurl_rep%ml_mesh%level-1)
if(cmesh%type/=1)CALL oft_abort("Only supported with tet meshes", &
  "hcurl_ginterpmatrix", __FILE__)
! ops=>oft_hcurl_ops
SELECT TYPE(this=>ML_hcurl_rep%levels(ML_hcurl_rep%level)%fe)
  CLASS IS(oft_hcurl_fem)
    hcurl_fine=>this
  CLASS DEFAULT
    CALL oft_abort("Error casting fine level", "hcurl_ginterpmatrix", __FILE__)
END SELECT
if(hcurl_fine%order/=1)then
  call oft_abort('Attempted geometric interpolation for pd > 1', &
    'hcurl_ginterpmatrix',__FILE__)
end if
SELECT TYPE(this=>ML_hcurl_rep%levels(ML_hcurl_rep%level-1)%fe)
  CLASS IS(oft_hcurl_fem)
    hcurl_cors=>this
  CLASS DEFAULT
    CALL oft_abort("Error casting coarse level", "hcurl_ginterpmatrix", __FILE__)
END SELECT
! hcurl_fine=>ML_hcurl_rep%levels(ML_hcurl_rep%level-1)%fe
lede=>ML_hcurl_rep%ml_mesh%inter(ML_hcurl_rep%ml_mesh%level-1)%lede
lfde=>ML_hcurl_rep%ml_mesh%inter(ML_hcurl_rep%ml_mesh%level-1)%lfde
lcdg=>ML_hcurl_rep%ml_mesh%inter(ML_hcurl_rep%ml_mesh%level-1)%lcdg
lcde=>ML_hcurl_rep%ml_mesh%inter(ML_hcurl_rep%ml_mesh%level-1)%lcde
ALLOCATE(ML_hcurl_rep%interp_graphs(ML_hcurl_rep%level)%g)
interp_graph=>ML_hcurl_rep%interp_graphs(ML_hcurl_rep%level)%g
!---Setup matrix sizes
interp_graph%nr=hcurl_fine%ne
interp_graph%nrg=hcurl_fine%global%ne
interp_graph%nc=hcurl_cors%ne
interp_graph%ncg=hcurl_cors%global%ne
!---Setup Matrix graph
ALLOCATE(interp_graph%kr(interp_graph%nr+1))
interp_graph%nnz=2*cmesh%ne+9*cmesh%nf+6*cmesh%nc
ALLOCATE(interp_graph%lc(interp_graph%nnz))
interp_graph%lc=0_i4
!---Construct linkage
DO i=1,cmesh%ne
  interp_graph%kr(lede(1,i))=1
  interp_graph%kr(lede(2,i))=1
END DO
DO i=1,cmesh%nf
  jfe=ABS(cmesh%lfe(:,i)) ! face edges
  interp_graph%kr(ABS(lfde(1,i)))=3
  interp_graph%kr(ABS(lfde(2,i)))=3
  interp_graph%kr(ABS(lfde(3,i)))=3
END DO
DO i=1,cmesh%nc
  interp_graph%kr(ABS(lcde(1,i)))=6
END DO
interp_graph%kr(interp_graph%nr+1)=interp_graph%nnz+1
do i=interp_graph%nr,1,-1 ! cumulative point to point count
  interp_graph%kr(i)=interp_graph%kr(i+1)-interp_graph%kr(i)
end do
if(interp_graph%kr(1)/=1)call oft_abort('Bad element to element count','oft_hcurl_ginterpmatrix',__FILE__)
DO i=1,cmesh%ne
  k=interp_graph%kr(lede(1,i))
  interp_graph%lc(k)=i
  k=interp_graph%kr(lede(2,i))
  interp_graph%lc(k)=i
END DO
DO i=1,cmesh%nf
  jfe=ABS(cmesh%lfe(:,i)) ! face edges
  CALL sort_array(jfe,3)
  !---
  DO j=1,3
    k=interp_graph%kr(ABS(lfde(j,i)))
    interp_graph%lc(k)=jfe(1)
    interp_graph%lc(k+1)=jfe(2)
    interp_graph%lc(k+2)=jfe(3)
  END DO
END DO
DO i=1,cmesh%nc
  jce=ABS(cmesh%lce(:,i)) ! face edges
  CALL sort_array(jce,6)
  !---
  k=interp_graph%kr(ABS(lcde(1,i)))
  DO j=0,5
    interp_graph%lc(k+j)=jce(1+j)
  END DO
END DO
!------------------------------------------------------------------------------
! Construct matrix
!------------------------------------------------------------------------------
NULLIFY(hcurl_vec_fine,hcurl_vec_cors)
CALL ML_hcurl_rep%vec_create(hcurl_vec_fine)
CALL ML_hcurl_rep%vec_create(hcurl_vec_cors,ML_hcurl_rep%level-1)
!---
ALLOCATE(graphs(1,1))
graphs(1,1)%g=>interp_graph
!---
CALL create_matrix(mat,graphs,hcurl_vec_fine,hcurl_vec_cors)
CALL hcurl_vec_fine%delete
CALL hcurl_vec_cors%delete
DEALLOCATE(graphs,hcurl_vec_fine,hcurl_vec_cors)
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Construct matrix
allocate(emap(cmesh%ne))
CALL get_inverse_map(cmesh%lbe,cmesh%nbe,emap,cmesh%ne)
DO i=1,cmesh%ne
  IF(cmesh%be(i))THEN
    ! IF(.NOT.cmesh%linkage%leo(emap(i)))CYCLE
    IF(.NOT.cmesh%estitch%leo(emap(i)))CYCLE
  END IF
  !---
  i_ind=lede(1,i)
  j_ind=i
  mop=1.d0*SIGN(1_i8,cmesh%global%le(j_ind(1)))*SIGN(1_i8,mesh%global%le(i_ind(1)))
  CALL mat%add_values(i_ind,j_ind,mop/2.d0,1,1)
  !---
  i_ind=lede(2,i)
  j_ind=i
  mop=1.d0*SIGN(1_i8,cmesh%global%le(j_ind(1)))*SIGN(1_i8,mesh%global%le(i_ind(1)))
  CALL mat%add_values(i_ind,j_ind,-mop/2.d0,1,1)
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
  i_ind=lfde(1,i)
  j_ind=jfe(1)
  mop=1.d0/4.d0*SIGN(1_i8,cmesh%global%le(j_ind(1)))*SIGN(1_i8,mesh%global%le(i_ind(1)))
  CALL mat%add_values(i_ind,j_ind,mop,1,1)
  j_ind=jfe(2)
  mop=1.d0/4.d0*SIGN(1_i8,cmesh%global%le(j_ind(1)))*SIGN(1_i8,mesh%global%le(i_ind(1)))
  CALL mat%add_values(i_ind,j_ind,mop,1,1)
  j_ind=jfe(3)
  mop=-1.d0/4.d0*SIGN(1_i8,cmesh%global%le(j_ind(1)))*SIGN(1_i8,mesh%global%le(i_ind(1)))
  CALL mat%add_values(i_ind,j_ind,mop,1,1)
  !---
  i_ind=lfde(2,i)
  j_ind=jfe(1)
  mop=1.d0/4.d0*SIGN(1_i8,cmesh%global%le(j_ind(1)))*SIGN(1_i8,mesh%global%le(i_ind(1)))
  CALL mat%add_values(i_ind,j_ind,mop,1,1)
  j_ind=jfe(2)
  mop=1.d0/4.d0*SIGN(1_i8,cmesh%global%le(j_ind(1)))*SIGN(1_i8,mesh%global%le(i_ind(1)))
  CALL mat%add_values(i_ind,j_ind,mop,1,1)
  j_ind=jfe(3)
  mop=1.d0/4.d0*SIGN(1_i8,cmesh%global%le(j_ind(1)))*SIGN(1_i8,mesh%global%le(i_ind(1)))
  CALL mat%add_values(i_ind,j_ind,mop,1,1)
  !---
  i_ind=lfde(3,i)
  j_ind=jfe(2)
  mop=1.d0/4.d0*SIGN(1_i8,cmesh%global%le(j_ind(1)))*SIGN(1_i8,mesh%global%le(i_ind(1)))
  CALL mat%add_values(i_ind,j_ind,mop,1,1)
  j_ind=jfe(1)
  mop=-1.d0/4.d0*SIGN(1_i8,cmesh%global%le(j_ind(1)))*SIGN(1_i8,mesh%global%le(i_ind(1)))
  CALL mat%add_values(i_ind,j_ind,mop,1,1)
  j_ind=jfe(3)
  mop=1.d0/4.d0*SIGN(1_i8,cmesh%global%le(j_ind(1)))*SIGN(1_i8,mesh%global%le(i_ind(1)))
  CALL mat%add_values(i_ind,j_ind,mop,1,1)
END DO
deallocate(fmap)
f=1.d0/4.d0
DO i=1,cmesh%nc ! loop over coarse cells
  jce=ABS(cmesh%lce(:,i)) ! cell edge indices
  d=(cmesh%r(:,cmesh%le(1,jce(lcdg(i)+3)))-cmesh%r(:,cmesh%le(1,jce(lcdg(i)))) &
  +cmesh%r(:,cmesh%le(2,jce(lcdg(i)+3)))-cmesh%r(:,cmesh%le(2,jce(lcdg(i)))))/2
  CALL cmesh%jacobian(i,f,goptmp,v)
  !---
  CALL oft_hcurl_eval_all(hcurl_cors,i,f,h_rop,goptmp)
  DO j=1,6
    i_ind=ABS(lcde(1,i))
    j_ind=jce(j)
    mop=SIGN(1_i4,lcde(1,i))*DOT_PRODUCT(h_rop(:,j),d)*SIGN(1_i8,mesh%global%le(i_ind(1)))
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
  END DO
END DO
DEBUG_STACK_POP
END SUBROUTINE hcurl_ginterpmatrix
!------------------------------------------------------------------------------
!> Construct interpolation matrix for polynomial levels
!------------------------------------------------------------------------------
SUBROUTINE hcurl_pinterpmatrix(mat)
class(oft_matrix), pointer, intent(inout) :: mat !< Interolation matrix
INTEGER(i4) :: i,j,k,m,icors,ifine,jb,js,jn,i_ind(1),j_ind(1)
INTEGER(i4) :: etmp(2),ftmp(3),fetmp(3),ctmp(4),fc,ed
INTEGER(i4) :: offsetc,offsetf
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pmap,emap,fmap
REAL(r8) :: f(4),val,mop(1)
TYPE(oft_fem_type), POINTER :: hcurl_cors => NULL()
CLASS(oft_hcurl_fem), POINTER :: hcurl_fine => NULL()
CLASS(oft_vector), POINTER :: hcurl_vec_fine,hcurl_vec_cors
type(oft_graph_ptr), pointer :: graphs(:,:)
type(oft_graph), POINTER :: interp_graph
CLASS(oft_mesh), POINTER :: mesh
DEBUG_STACK_PUSH
!---
SELECT TYPE(this=>ML_hcurl_rep%levels(ML_hcurl_rep%level)%fe)
CLASS IS(oft_hcurl_fem)
  hcurl_fine=>this
  mesh=>this%mesh
CLASS DEFAULT
  CALL oft_abort("Error getting coarse FE object","hcurl_pinterpmatrix",__FILE__)
END SELECT
SELECT TYPE(this=>ML_hcurl_rep%levels(ML_hcurl_rep%level-1)%fe)
CLASS IS(oft_fem_type)
  hcurl_cors=>this
CLASS DEFAULT
  CALL oft_abort("Error getting coarse FE object","hcurl_pinterpmatrix",__FILE__)
END SELECT
ALLOCATE(ML_hcurl_rep%interp_graphs(ML_hcurl_rep%level)%g)
interp_graph=>ML_hcurl_rep%interp_graphs(ML_hcurl_rep%level)%g
!---Setup matrix sizes
interp_graph%nr=hcurl_fine%ne
interp_graph%nrg=hcurl_fine%global%ne
interp_graph%nc=hcurl_cors%ne
interp_graph%ncg=hcurl_cors%global%ne
!---Setup Matrix graph
ALLOCATE(interp_graph%kr(interp_graph%nr+1))
interp_graph%kr=0
interp_graph%nnz=hcurl_cors%ne
ALLOCATE(interp_graph%lc(interp_graph%nnz))
interp_graph%lc=0_i4
!---Construct matrix
do i=1,mesh%ne
  do j=1,hcurl_cors%gstruct(2)
    offsetf=(i-1)*hcurl_fine%gstruct(2)
    interp_graph%kr(j+offsetf)=1
  end do
end do
!---
do i=1,mesh%nf
  do j=1,hcurl_cors%gstruct(3)
    offsetf=mesh%ne*hcurl_fine%gstruct(2)+(i-1)*hcurl_fine%gstruct(3)
    interp_graph%kr(j+offsetf)=1
  end do
end do
!---
do i=1,mesh%nc
  do j=1,hcurl_cors%gstruct(4)
    offsetf=mesh%ne*hcurl_fine%gstruct(2)+hcurl_fine%gstruct(3)*mesh%nf+(i-1)*hcurl_fine%gstruct(4)
    interp_graph%kr(j+offsetf)=1
  end do
end do
interp_graph%kr(interp_graph%nr+1)=interp_graph%nnz+1
do i=interp_graph%nr,1,-1 ! cumulative point to point count
  interp_graph%kr(i)=interp_graph%kr(i+1)-interp_graph%kr(i)
end do
if(interp_graph%kr(1)/=1)call oft_abort('Bad element to element count','hcurl_fine_pinterpmatrix',__FILE__)
!---Construct matrix
do i=1,mesh%ne
  do j=1,hcurl_cors%gstruct(2)
    offsetf=(i-1)*hcurl_fine%gstruct(2)
    offsetc=(i-1)*hcurl_cors%gstruct(2)
    interp_graph%lc(interp_graph%kr(j+offsetf))=j+offsetc
  end do
end do
!---
do i=1,mesh%nf
  do j=1,hcurl_cors%gstruct(3)
    offsetf=mesh%ne*hcurl_fine%gstruct(2)+(i-1)*hcurl_fine%gstruct(3)
    offsetc=mesh%ne*hcurl_cors%gstruct(2)+(i-1)*hcurl_cors%gstruct(3)
    interp_graph%lc(interp_graph%kr(j+offsetf))=j+offsetc
  end do
end do
!---
do i=1,mesh%nc
  do j=1,hcurl_cors%gstruct(4)
    offsetf=mesh%ne*hcurl_fine%gstruct(2)+hcurl_fine%gstruct(3)*mesh%nf+(i-1)*hcurl_fine%gstruct(4)
    offsetc=mesh%ne*hcurl_cors%gstruct(2)+hcurl_cors%gstruct(3)*mesh%nf+(i-1)*hcurl_cors%gstruct(4)
    interp_graph%lc(interp_graph%kr(j+offsetf))=j+offsetc
  end do
end do
!------------------------------------------------------------------------------
! Construct matrix
!------------------------------------------------------------------------------
NULLIFY(hcurl_vec_fine,hcurl_vec_cors)
CALL ML_hcurl_rep%vec_create(hcurl_vec_fine)
CALL ML_hcurl_rep%vec_create(hcurl_vec_cors,ML_hcurl_rep%level-1)
!---
ALLOCATE(graphs(1,1))
graphs(1,1)%g=>interp_graph
!---
CALL create_matrix(mat,graphs,hcurl_vec_fine,hcurl_vec_cors)
CALL hcurl_vec_fine%delete
CALL hcurl_vec_cors%delete
DEALLOCATE(graphs,hcurl_vec_fine,hcurl_vec_cors)
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---
allocate(emap(mesh%ne))
CALL get_inverse_map(mesh%lbe,mesh%nbe,emap,mesh%ne)
DO i=1,mesh%ne
  IF(mesh%be(i))THEN
    ! IF(.NOT.mesh%linkage%leo(emap(i)))CYCLE
    IF(.NOT.mesh%estitch%leo(emap(i)))CYCLE
  END IF
  DO j=1,hcurl_cors%gstruct(2)
    i_ind=j+(i-1)*hcurl_fine%gstruct(2)
    j_ind=j+(i-1)*hcurl_cors%gstruct(2)
    mop=1.d0
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
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
  DO j=1,hcurl_cors%gstruct(3)
    i_ind=j+mesh%ne*hcurl_fine%gstruct(2)+(i-1)*hcurl_fine%gstruct(3)
    j_ind=j+mesh%ne*hcurl_cors%gstruct(2)+(i-1)*hcurl_cors%gstruct(3)
    mop=1.d0
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
  END DO
END DO
deallocate(fmap)
!---
DO i=1,mesh%nc
  DO j=1,hcurl_cors%gstruct(4)
    i_ind=j+mesh%ne*hcurl_fine%gstruct(2)+hcurl_fine%gstruct(3)*mesh%nf+(i-1)*hcurl_fine%gstruct(4)
    j_ind=j+mesh%ne*hcurl_cors%gstruct(2)+hcurl_cors%gstruct(3)*mesh%nf+(i-1)*hcurl_cors%gstruct(4)
    mop=1.d0
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
  END DO
END DO
DEBUG_STACK_POP
END SUBROUTINE hcurl_pinterpmatrix
END SUBROUTINE hcurl_setup_interp
!------------------------------------------------------------------------------
!> Transfer a base level H(Curl) vector field to the next MPI level
!------------------------------------------------------------------------------
subroutine hcurl_base_pop(self,acors,afine)
class(oft_ml_fe_vecspace), intent(inout) :: self
class(oft_vector), intent(inout) :: acors !< Vector to transfer
class(oft_vector), intent(inout) :: afine !< Fine vector from transfer
integer(i4), pointer :: lbege(:)
integer(i4) :: i
real(r8), pointer, dimension(:) :: array_c,array_f
DEBUG_STACK_PUSH
lbege=>self%ML_FE_rep%ml_mesh%inter(self%ML_FE_rep%ml_mesh%nbase)%lbege
CALL acors%get_local(array_c)
CALL afine%get_local(array_f)
!$omp parallel do
do i=1,afine%n
  array_f(i)=array_c(ABS(lbege(i)))
end do
CALL afine%restore_local(array_f)
DEALLOCATE(array_c,array_f)
DEBUG_STACK_POP
end subroutine hcurl_base_pop
!------------------------------------------------------------------------------
!> Transfer a MPI level H(Curl) vector field to the base level
!------------------------------------------------------------------------------
subroutine hcurl_base_push(self,afine,acors)
class(oft_ml_fe_vecspace), intent(inout) :: self
class(oft_vector), intent(inout) :: afine !< Vector to transfer
class(oft_vector), intent(inout) :: acors !< Fine vector from transfer
integer(i4), pointer :: lbege(:)
integer(i4) :: i,j,ierr
real(r8), pointer, dimension(:) :: alias,array_c,array_f
CLASS(oft_afem_type), POINTER :: hcurl_fine => NULL()
DEBUG_STACK_PUSH
!---
lbege=>self%ML_FE_rep%ml_mesh%inter(self%ML_FE_rep%ml_mesh%nbase)%lbege
CALL acors%get_local(array_c)
CALL afine%get_local(array_f)
hcurl_fine=>self%ML_FE_rep%levels(self%ML_FE_rep%level+1)%fe
!---
allocate(alias(acors%n))
alias=0.d0
!$omp parallel do
do i=1,afine%n
  if(hcurl_fine%linkage%be(i))cycle
  alias(ABS(lbege(i)))=array_f(i)
end do
!$omp parallel do private(j)
do i=1,hcurl_fine%linkage%nbe
  j=hcurl_fine%linkage%lbe(i)
  if(.NOT.hcurl_fine%linkage%leo(i))cycle
  alias(ABS(lbege(j)))=array_f(j)
end do
!---Global reduction over all processors
array_c=oft_mpi_sum(alias,acors%n)
call acors%restore_local(array_c)
deallocate(alias,array_c,array_f)
DEBUG_STACK_POP
end subroutine hcurl_base_push
!------------------------------------------------------------------------------
!> Compute eigenvalues and smoothing coefficients for the operator H(Curl)::WOP
!------------------------------------------------------------------------------
SUBROUTINE hcurl_wop_eigs(ML_hcurl_rep,minlev)
type(oft_ml_fem_type), target, intent(inout) :: ML_hcurl_rep
INTEGER(i4), INTENT(in) :: minlev !< Needs docs
#ifdef HAVE_ARPACK
INTEGER(i4) :: i
REAL(r8) :: lam0
REAL(r8), ALLOCATABLE :: df(:)
CLASS(oft_vector), POINTER :: u
TYPE(oft_irlm_eigsolver) :: arsolver
CLASS(oft_matrix), POINTER :: md => NULL()
CLASS(oft_matrix), POINTER :: wop => NULL()
TYPE(oft_hcurl_zerob), TARGET :: bc_tmp
DEBUG_STACK_PUSH
!------------------------------------------------------------------------------
! Compute optimal smoother coefficients
!------------------------------------------------------------------------------
IF(oft_env%head_proc)WRITE(*,*)'Optimizing Jacobi damping for H(Curl)::WOP'
bc_tmp%ML_hcurl_rep=>ML_hcurl_rep
ALLOCATE(df(ML_hcurl_rep%nlevels))
df=0.d0
DO i=minlev,ML_hcurl_rep%nlevels
  CALL ML_hcurl_rep%set_level(i)
  !---Create fields
  CALL ML_hcurl_rep%vec_create(u)
  !---Get Ev range
  NULLIFY(wop)
  SELECT TYPE(this=>ML_hcurl_rep%current_level)
  CLASS IS(oft_hcurl_fem)
    CALL oft_hcurl_getwop(this,wop,'zerob')
  CLASS DEFAULT
    CALL oft_abort("Error getting current FE rep","hcurl_wop_eigs",__FILE__)
  END SELECT
  ! CALL oft_hcurl_getwop(oft_hcurl,wop,'zerob')
  CALL create_diagmatrix(md,wop%D)
  !---
  arsolver%A=>wop
  arsolver%M=>md
  arsolver%mode=2
  arsolver%tol=1.E-5_r8
  arsolver%bc=>bc_tmp
  CALL create_native_pre(arsolver%Minv, "jacobi")
  arsolver%Minv%A=>wop
  !---
  IF(oft_debug_print(1))WRITE(*,*)'  optimizing level = ',i
  CALL arsolver%max(u,lam0)
  df(i) = 1.8d0/lam0
  !---
  CALL u%delete
  CALL md%delete
  CALL wop%delete
END DO
!---Output
IF(oft_env%head_proc)THEN
  WRITE(*,'(A)',ADVANCE='NO')' df_wop = '
  DO i=1,ML_hcurl_rep%nlevels-1
    WRITE(*,'(F5.3,A)',ADVANCE='NO')df(i),', '
  END DO
  WRITE(*,'(F5.3,A)')df(ML_hcurl_rep%nlevels)
END IF
DEALLOCATE(df)
DEBUG_STACK_POP
#else
CALL oft_abort("Subroutine requires ARPACK", "lag_lop_eigs", __FILE__)
#endif
END SUBROUTINE hcurl_wop_eigs
!------------------------------------------------------------------------------
!> Construct default MG preconditioner for H(Curl)::WOP
!------------------------------------------------------------------------------
SUBROUTINE hcurl_getwop_pre(ML_hcurl_rep,pre,mats,level,nlevels)
type(oft_ml_fem_type), target, intent(inout) :: ML_hcurl_rep
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
TYPE(oft_hcurl_zerob), POINTER :: bc_tmp
TYPE(oft_ml_fe_vecspace), POINTER :: tmp_vecspace
!---
TYPE(xml_node), POINTER :: pre_node
#ifdef HAVE_XML
integer(i4) :: nnodes
TYPE(xml_node), POINTER :: hcurl_node
#endif
DEBUG_STACK_PUSH
!---
minlev=1
toplev=ML_hcurl_rep%level
levin=ML_hcurl_rep%level
IF(PRESENT(level))toplev=level
IF(PRESENT(nlevels))minlev=toplev-nlevels+1
nl=toplev-minlev+1
!---
IF(minlev<ML_hcurl_rep%minlev)CALL oft_abort('Minimum level is < minlev','hcurl_getwop_pre',__FILE__)
IF(toplev>ML_hcurl_rep%nlevels)CALL oft_abort('Maximum level is > hcurl_nlevels','hcurl_getwop_pre',__FILE__)
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
  CALL ML_hcurl_rep%set_level(minlev+(i-1))
  levels(i)=minlev+(i-1)
  df(i)=df_wop(levels(i))
  nu(i)=nu_wop(levels(i))
  IF(df(i)<-1.d90)THEN
    WRITE(lev_char,'(I2.2)')levels(i)
    CALL oft_abort('Smoother values not set for level: '//lev_char,'hcurl_getwop_pre',__FILE__)
  END IF
  !---
  IF(create_mats)THEN
    NULLIFY(mats(i)%M)
    CALL oft_hcurl_getwop(ML_hcurl_rep%current_level,mats(i)%M,'zerob')
  END IF
  IF(i>1)ml_int(i-1)%M=>ML_hcurl_rep%interp_matrices(ML_hcurl_rep%level)%m !oft_hcurl_ops%interp
END DO
CALL ML_hcurl_rep%set_level(levin)
!------------------------------------------------------------------------------
! Search for XML-spec
!------------------------------------------------------------------------------
NULLIFY(pre_node)
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml))THEN
  CALL xml_get_element(oft_env%xml,"hcurl",hcurl_node,ierr)
  IF(ierr==0)CALL xml_get_element(hcurl_node,"wop",pre_node,ierr)
END IF
#endif
!------------------------------------------------------------------------------
! Setup preconditioner
!------------------------------------------------------------------------------
ALLOCATE(bc_tmp)
bc_tmp%ML_hcurl_rep=>ML_hcurl_rep
NULLIFY(pre)
ALLOCATE(tmp_vecspace)
tmp_vecspace%ML_FE_rep=>ML_hcurl_rep
tmp_vecspace%base_pop=>hcurl_base_pop
tmp_vecspace%base_push=>hcurl_base_push
CALL create_mlpre(pre,mats(1:nl),levels,nlevels=nl,ml_vecspace=tmp_vecspace, &
      bc=bc_tmp,stype=1,df=df,nu=nu,xml_root=pre_node)
!------------------------------------------------------------------------------
! Cleanup
!------------------------------------------------------------------------------
DEALLOCATE(ml_int,levels,df,nu)
DEBUG_STACK_POP
END SUBROUTINE hcurl_getwop_pre
!------------------------------------------------------------------------------
!> Construct default MG preconditioner for H(Curl)::JMLB
!------------------------------------------------------------------------------
SUBROUTINE hcurl_getjmlb_pre(ML_hcurl_rep,pre,mats,alam,level,nlevels)
type(oft_ml_fem_type), target, intent(inout) :: ML_hcurl_rep
CLASS(oft_solver), POINTER, INTENT(out) :: pre !< Needs docs
TYPE(oft_matrix_ptr), POINTER, INTENT(inout) :: mats(:) !< Needs docs
REAL(r8), INTENT(in) :: alam !< Needs docs
INTEGER(i4), OPTIONAL, INTENT(in) :: level !< Needs docs
INTEGER(i4), OPTIONAL, INTENT(in) :: nlevels !< Needs docs
!--- ML structures for MG-preconditioner
TYPE(oft_matrix_ptr), POINTER :: ml_int(:)
INTEGER(i4), ALLOCATABLE :: levels(:)
INTEGER(i4), ALLOCATABLE :: nu(:)
!---
INTEGER(i4) :: minlev,toplev,nl
INTEGER(i4) :: i,j,levin,ierr
LOGICAL :: create_mats
CHARACTER(LEN=2) :: lev_char
CLASS(oft_vector), POINTER :: oft_hcurl_vec
TYPE(oft_hcurl_zerob), POINTER :: bc_tmp
TYPE(oft_ml_fe_vecspace), POINTER :: tmp_vecspace
!---
TYPE(xml_node), POINTER :: pre_node
#ifdef HAVE_XML
integer(i4) :: nnodes
TYPE(xml_node), POINTER :: hcurl_node
#endif
DEBUG_STACK_PUSH
!---
minlev=1
toplev=ML_hcurl_rep%level
levin=ML_hcurl_rep%level
IF(PRESENT(level))toplev=level
IF(PRESENT(nlevels))minlev=toplev-nlevels+1
nl=toplev-minlev+1
!---
IF(minlev<ML_hcurl_rep%minlev)CALL oft_abort('Minimum level is < minlev','hcurl_getjmlb_pre',__FILE__)
IF(toplev>ML_hcurl_rep%nlevels)CALL oft_abort('Maximum level is > hcurl_nlevels','hcurl_getjmlb_pre',__FILE__)
!------------------------------------------------------------------------------
! Create ML Matrices
!------------------------------------------------------------------------------
create_mats=.FALSE.
IF(.NOT.ASSOCIATED(mats))THEN
  create_mats=.TRUE.
  ALLOCATE(mats(nl))
END IF
ALLOCATE(ml_int(nl-1),levels(nl),nu(nl))
DO i=1,nl
  CALL ML_hcurl_rep%set_level(minlev+(i-1))
  levels(i)=minlev+(i-1)
  nu(i)=nu_jmlb(levels(i))
  !---
  IF(create_mats)THEN
    NULLIFY(mats(i)%M)
    CALL oft_hcurl_getjmlb(ML_hcurl_rep%current_level,mats(i)%M,alam,'zerob')
  END IF
  IF(i>1)ml_int(i-1)%M=>ML_hcurl_rep%interp_matrices(ML_hcurl_rep%level)%m !oft_hcurl_ops%interp
END DO
CALL ML_hcurl_rep%set_level(levin)
!------------------------------------------------------------------------------
! Search for XML-spec
!------------------------------------------------------------------------------
NULLIFY(pre_node)
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml))THEN
  CALL xml_get_element(oft_env%xml,"hcurl",hcurl_node,ierr)
  IF(ierr==0)CALL xml_get_element(hcurl_node,"jmlb",pre_node,ierr)
END IF
#endif
!------------------------------------------------------------------------------
! Setup preconditioner
!------------------------------------------------------------------------------
ALLOCATE(bc_tmp)
bc_tmp%ML_hcurl_rep=>ML_hcurl_rep
NULLIFY(pre)
ALLOCATE(tmp_vecspace)
tmp_vecspace%ML_FE_rep=>ML_hcurl_rep
tmp_vecspace%base_pop=>hcurl_base_pop
tmp_vecspace%base_push=>hcurl_base_push
CALL create_mlpre(pre,mats(1:nl),levels,nlevels=nl,ml_vecspace=tmp_vecspace, &
      bc=bc_tmp,stype=2,nu=nu,xml_root=pre_node)
!------------------------------------------------------------------------------
! Cleanup
!------------------------------------------------------------------------------
DEALLOCATE(ml_int,levels,nu)
DEBUG_STACK_POP
END SUBROUTINE hcurl_getjmlb_pre
end module oft_hcurl_operators
