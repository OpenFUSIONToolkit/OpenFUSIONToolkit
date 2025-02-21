!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_hcurl_operators.F90
!
!> Nedelec H1(Curl) FE operator definitions
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
!---------------------------------------------------------------------------
module oft_hcurl_operators
USE oft_base
USE oft_sort, ONLY: sort_array
USE oft_mesh_type, ONLY: oft_mesh, mesh, smesh, cell_is_curved
USE multigrid, ONLY: mg_mesh
!---
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, oft_matrix, oft_matrix_ptr, &
  oft_graph_ptr
USE oft_deriv_matrices, ONLY: oft_diagmatrix, create_diagmatrix
USE oft_solver_base, ONLY: oft_solver, oft_orthog, oft_bc_proto
#ifdef HAVE_ARPACK
USE oft_arpack, ONLY: oft_irlm_eigsolver
#endif
USE oft_la_utils, ONLY: create_matrix
USE oft_solver_utils, ONLY: create_mlpre, create_cg_solver, create_diag_pre, &
  create_native_pre
!---
USE fem_base, ONLY: oft_afem_type, oft_fem_type, fem_max_levels
USE fem_utils, ONLY: fem_interp
USE oft_lag_basis, ONLY: oft_lagrange, oft_lag_geval_all
USE oft_lag_fields, ONLY: oft_lag_create
USE oft_lag_operators, ONLY: oft_lag_getlop, lag_zerob, lag_zerogrnd
USE oft_hcurl_basis, ONLY: oft_hcurl, ML_oft_hcurl, oft_bhcurl, oft_hcurl_level, oft_hcurl_blevel, &
oft_hcurl_nlevels, oft_nedelec_ops, oft_hcurl_ops, ML_oft_hcurl_ops, oft_hcurl_eval_all, &
oft_hcurl_ceval_all, oft_hcurl_set_level, oft_hcurl_lev, oft_hcurl_minlev, oft_hcurl_get_cgops, &
oft_hcurl_fem
USE oft_hcurl_fields, ONLY: oft_hcurl_create
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------
! CLASS oft_hcurl_orthog
!---------------------------------------------------------------------------
!> Orthogonalize a H1(Curl) vector against a library of given modes
!---------------------------------------------------------------------------
type, extends(oft_orthog) :: oft_hcurl_orthog
  integer(i4) :: nm = 0 !< Number of modes to orthogonalize against
  type(oft_vector_ptr), pointer, dimension(:,:) :: orthog => NULL() !< Library of modes
  class(oft_matrix), pointer :: wop => NULL() !< H1(Curl)::WOP, used as metric
contains
  !> Perform orthoganlization
  procedure :: apply => hcurl_orthog_apply
  !> Clean-up internal variables
  procedure :: delete => hcurl_orthog_delete
end type oft_hcurl_orthog
!---------------------------------------------------------------------------
! CLASS oft_hcurl_rinterp
!---------------------------------------------------------------------------
!> Interpolate a H1(Curl) field
!---------------------------------------------------------------------------
type, extends(fem_interp) :: oft_hcurl_rinterp
  class(oft_vector), pointer :: u => NULL() !< Field to interpolate
  real(r8), pointer, dimension(:) :: vals => NULL() !< Local values
  class(oft_hcurl_fem), pointer :: hcurl_rep => NULL() !< H1(Curl) FE representation
contains
  !> Retrieve local values for interpolation
  procedure :: setup => hcurl_rinterp_setup
  !> Reconstruct field
  procedure :: interp => hcurl_rinterp
  !> Destroy temporary internal storage
  procedure :: delete => hcurl_rinterp_delete
end type oft_hcurl_rinterp
!---------------------------------------------------------------------------
! CLASS oft_hcurl_cinterp
!---------------------------------------------------------------------------
!> Interpolate \f$ \nabla \times \f$ of a H1(Curl) field
!---------------------------------------------------------------------------
type, extends(oft_hcurl_rinterp) :: oft_hcurl_cinterp
contains
  !> Reconstruct field
  procedure :: interp => hcurl_cinterp
end type oft_hcurl_cinterp
!---------------------------------------------------------------------------
! CLASS oft_hcurl_divout
!---------------------------------------------------------------------------
!> Clean the divergence from a H1(Curl) vector field
!!
!! Divergence is removed by adding a gradient field, such that \f$ f = f +
!! \nabla \phi \f$, where \f$ \nabla^2 \phi = - \nabla \cdot f \f$. Cleaning
!! may also be applied to a field which is pre-multiplied by the H1(Curl)::MOP
!! in which case \f$ \nabla^2 \phi = - \nabla^T f \f$ and \f$ f = f + M \nabla
!! \phi \f$. The mass matrix version is applied if \c mop is associated with the
!! corresponding H1(Curl)::MOP.
!!
!! @note This only removes the 0-th order component, which is sufficient for
!! orthogonalization against the H1(Curl)::WOP null space. Higher order cleaning
!! requires the full H1 vector space.
!---------------------------------------------------------------------------
type, extends(oft_orthog) :: oft_hcurl_divout
  integer(i4) :: count=0 !< Number of times apply has been called
  integer(i4) :: app_freq=1 !< Frequency to apply solver
  logical :: pm = .FALSE. !< Flag for solver convergence monitor
  class(oft_solver), pointer :: solver => NULL() !< Solver object for LAG::LOP operator
  procedure(oft_bc_proto), pointer, nopass :: bc => NULL() !< Boundary condition
  class(oft_matrix), pointer :: mop => NULL() !< Mass matrix, applies divoutm if associated
contains
  !> Setup matrix and solver with default
  procedure :: setup => hcurl_divout_setup
  !> Clean divergence from field
  procedure :: apply => hcurl_divout_apply
  !> Clean-up internal storage
  procedure :: delete => hcurl_divout_delete
end type oft_hcurl_divout
!---Pre options
integer(i4) :: nu_wop(fem_max_levels)=0
real(r8) :: df_wop(fem_max_levels)=-1.d99
integer(i4) :: nu_jmlb(fem_max_levels)=20
!---Cache variables
REAL(r8), POINTER, DIMENSION(:,:) :: oft_hcurl_rop => NULL()
REAL(r8), POINTER, DIMENSION(:,:) :: oft_hcurl_cop => NULL()
!$omp threadprivate(oft_hcurl_rop,oft_hcurl_cop)
contains
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_mloptions
!---------------------------------------------------------------------------
!> Read-in options for the basic Nedelec H1(Curl) ML preconditioners
!---------------------------------------------------------------------------
subroutine hcurl_mloptions
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
    WRITE(*,*)'No H(Curl) MG smoother settings found:'
    WRITE(*,*)'  Using default values, which may result in convergence failure.'
  END IF
  DO i=oft_hcurl_minlev,oft_hcurl_nlevels
    CALL oft_hcurl_set_level(i)
    SELECT CASE(oft_hcurl%order)
    CASE(1)
      df_wop(i)=0.6d0
    CASE(2)
      df_wop(i)=0.3d0
    CASE(3)
      df_wop(i)=0.25d0
    CASE DEFAULT
      df_wop(i)=0.2d0
    END SELECT
    nu_wop(i)=MIN(64,2**(oft_hcurl_nlevels-i))
  END DO
  nu_wop(oft_hcurl_minlev)=64
  nu_wop(oft_hcurl_nlevels)=1
END IF
DEBUG_STACK_POP
end subroutine hcurl_mloptions
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_rinterp_setup
!---------------------------------------------------------------------------
!> Setup interpolator for H1(Curl) fields
!!
!! Fetches local representation used for interpolation from vector object
!!
!! @note Should only be used via class \ref oft_hcurl_rinterp or children
!---------------------------------------------------------------------------
subroutine hcurl_rinterp_setup(self)
class(oft_hcurl_rinterp), intent(inout) :: self
!---Get local slice
CALL self%u%get_local(self%vals)
self%hcurl_rep=>oft_hcurl
end subroutine hcurl_rinterp_setup
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_rinterp_delete
!---------------------------------------------------------------------------
!> Destroy temporary internal storage
!!
!! @note Should only be used via class \ref oft_hcurl_rinterp or children
!---------------------------------------------------------------------------
subroutine hcurl_rinterp_delete(self)
class(oft_hcurl_rinterp), intent(inout) :: self
!---Deallocate local storage
IF(ASSOCIATED(self%vals))DEALLOCATE(self%vals)
NULLIFY(self%hcurl_rep,self%u)
end subroutine hcurl_rinterp_delete
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_rinterp
!---------------------------------------------------------------------------
!> Reconstruct a Nedelec H1(Curl) vector field
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [3]
!---------------------------------------------------------------------------
subroutine hcurl_rinterp(self,cell,f,gop,val)
class(oft_hcurl_rinterp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
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
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_cinterp
!---------------------------------------------------------------------------
!> Reconstruct the curl of a Nedelec H1(Curl) vector field
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed gradient at f [3]
!---------------------------------------------------------------------------
subroutine hcurl_cinterp(self,cell,f,gop,val)
class(oft_hcurl_cinterp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
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
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_div
!---------------------------------------------------------------------------
!> Apply the divergence operator to a H1(Curl) field
!!
!! @note This subroutine computes the divergence of a H1(Curl) field as
!! projected on to a linear Lagrange scalar basis.
!---------------------------------------------------------------------------
subroutine hcurl_div(a,b)
class(oft_vector), intent(inout) :: a
class(oft_vector), intent(inout) :: b
real(r8), pointer, dimension(:) :: aloc
real(r8), pointer, dimension(:) :: bloc
integer(i4) :: i,jr,jc,m
integer(i4), allocatable :: j_curl(:),j_grad(:)
real(r8) :: goptmp(3,4)
real(r8) :: v,f(4),det,vol
logical :: curved
real(r8), allocatable :: rop_curl(:,:),rop_grad(:,:),btmp(:)
DEBUG_STACK_PUSH
!---------------------------------------------------------------------------
! Allocate Laplacian Op
!---------------------------------------------------------------------------
NULLIFY(aloc,bloc)
!---Cast to H1 type and get aliases
CALL a%get_local(aloc)
!---Zero output and get aliases
CALL b%set(0.d0)
CALL b%get_local(bloc)
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j_curl,j_grad,rop_curl,rop_grad,det,btmp,curved,goptmp,m,v,jc,jr,f)
allocate(j_curl(oft_hcurl%nce)) ! Local DOF and matrix indices
allocate(j_grad(oft_lagrange%nce)) ! Local DOF and matrix indices
allocate(rop_grad(3,oft_lagrange%nce)) ! Reconstructed gradient operator
allocate(rop_curl(3,oft_hcurl%nce)) ! Reconstructed gradient operator
allocate(btmp(oft_lagrange%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,mesh%nc
  !---Get local to global DOF mapping
  call oft_hcurl%ncdofs(i,j_curl)
  call oft_lagrange%ncdofs(i,j_grad)
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  btmp=0.d0
  do m=1,oft_hcurl%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call mesh%jacobian(i,oft_hcurl%quad%pts(:,m),goptmp,v)
    det=v*oft_hcurl%quad%wts(m)
    call oft_hcurl_eval_all(oft_hcurl,i,oft_hcurl%quad%pts(:,m),rop_curl,goptmp)
    call oft_lag_geval_all(oft_lagrange,i,oft_hcurl%quad%pts(:,m),rop_grad,goptmp)
    !---Compute local operator contribution
    do jr=1,oft_lagrange%nce
      do jc=1,oft_hcurl%nce
        btmp(jr)=btmp(jr)+DOT_PRODUCT(rop_grad(:,jr),rop_curl(:,jc))*aloc(j_curl(jc))*det
      end do
    end do
  end do
  !---Add local values to global vector
  do jr=1,oft_lagrange%nce
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
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_grad
!---------------------------------------------------------------------------
!> Add the gradient of a linear Lagrange scalar field to a H1(Curl) field
!---------------------------------------------------------------------------
subroutine hcurl_grad(a,b)
class(oft_vector), intent(inout) :: a
class(oft_vector), intent(inout) :: b
real(r8), pointer, dimension(:) :: aloc,bloc
integer(i4), allocatable :: emap(:)
integer :: i,j,k,l
real(r8) :: reg
DEBUG_STACK_PUSH
!---Get local values
NULLIFY(aloc,bloc)
CALL a%get_local(aloc)
CALL b%get_local(bloc)
!---Get boundary map
do i=1,mesh%ne
  bloc((i-1)*oft_hcurl%gstruct(2)+1)=bloc((i-1)*oft_hcurl%gstruct(2)+1) + &
  (aloc(mesh%le(2,i))-aloc(mesh%le(1,i)))*SIGN(1_i8,mesh%global%le(i))
enddo
!---
CALL b%restore_local(bloc)
DEALLOCATE(aloc,bloc)
DEBUG_STACK_POP
end subroutine hcurl_grad
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_gradtp
!---------------------------------------------------------------------------
!> Apply the transposed gradient operator to a Nedelec H1(Curl) vector field.
!!
!! @note Only the 0th order component is computed from the discrete gradient
!! operator
!!
!! @param[in,out] a Input field
!! @param[in,out] b \f$ G^{T} a \f$
!---------------------------------------------------------------------------
subroutine hcurl_gradtp(a,b)
class(oft_vector), intent(inout) :: a
class(oft_vector), intent(inout) :: b
real(r8), pointer, dimension(:) :: aloc,bloc
integer(i4), allocatable, dimension(:) :: emap
integer :: i,j,k
DEBUG_STACK_PUSH
!---Get local values
NULLIFY(aloc,bloc)
CALL a%get_local(aloc)
CALL b%get_local(bloc)
!---
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
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_zerob
!---------------------------------------------------------------------------
!> Zero the tangential component of a Nedelec H1(Curl) vector field on the boundary
!!
!! @param[in,out] a Field to be zeroed
!---------------------------------------------------------------------------
subroutine hcurl_zerob(a)
class(oft_vector), intent(inout) :: a
real(r8), pointer, dimension(:) :: aloc
integer(i4) :: i,j
DEBUG_STACK_PUSH
! Cast to vector type
NULLIFY(aloc)
CALL a%get_local(aloc)
! Apply operator
do i=1,oft_hcurl%nbe
  j=oft_hcurl%lbe(i)
  if(oft_hcurl%global%gbe(j))aloc(j)=0.d0
end do
CALL a%restore_local(aloc)
DEALLOCATE(aloc)
DEBUG_STACK_POP
end subroutine hcurl_zerob
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_getmop
!---------------------------------------------------------------------------
!> Construct mass matrix for a H1(Curl) representation
!!
!! Supported boundary conditions
!! - \c 'none' Full matrix
!! - \c 'zerob' Dirichlet for all boundary DOF
!!
!! @param[in,out] mat Matrix object
!! @param[in] bc Boundary condition
!---------------------------------------------------------------------------
subroutine oft_hcurl_getmop(mat,bc)
class(oft_matrix), pointer, intent(inout) :: mat
character(LEN=*), intent(in) :: bc
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time
real(r8), allocatable :: rop(:,:),mop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_hcurl_vec
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing H1(Curl)::MOP'
  CALL mytimer%tick()
END IF
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL oft_hcurl%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,rop,det,mop,curved,goptmp,m,vol,jc,jr)
allocate(j(oft_hcurl%nce)) ! Local DOF and matrix indices
allocate(rop(3,oft_hcurl%nce)) ! Reconstructed gradient operator
allocate(mop(oft_hcurl%nce,oft_hcurl%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  mop=0.d0
  do m=1,oft_hcurl%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call mesh%jacobian(i,oft_hcurl%quad%pts(:,m),goptmp,vol)
    det=vol*oft_hcurl%quad%wts(m)
    call oft_hcurl_eval_all(oft_hcurl,i,oft_hcurl%quad%pts(:,m),rop,goptmp)
    !---Compute local matrix contributions
    do jr=1,oft_hcurl%nce
      do jc=1,oft_hcurl%nce
        mop(jr,jc) = mop(jr,jc) + DOT_PRODUCT(rop(:,jr),rop(:,jc))*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call oft_hcurl%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,oft_hcurl%nce
        IF(oft_hcurl%global%gbe(j(jr)))THEN
          mop(jr,:)=0.d0
        END IF
      END DO
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,mop,oft_hcurl%nce,oft_hcurl%nce)
  ! !$omp end critical
end do
deallocate(j,rop,mop)
!$omp end parallel
ALLOCATE(mop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
SELECT CASE(TRIM(bc))
  CASE("zerob")
    mop(1,1)=1.d0
    DO i=1,oft_hcurl%nbe
      jr=oft_hcurl%lbe(i)
      IF(.NOT.oft_hcurl%global%gbe(jr))CYCLE
      IF(.NOT.oft_hcurl%linkage%leo(i))CYCLE
      j=jr
      call mat%add_values(j,j,mop,1,1)
    END DO
END SELECT
DEALLOCATE(j,mop)
CALL oft_hcurl_create(oft_hcurl_vec)
CALL mat%assemble(oft_hcurl_vec)
CALL oft_hcurl_vec%delete
DEALLOCATE(oft_hcurl_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_hcurl_getmop
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_getkop
!---------------------------------------------------------------------------
!> Construct helicity matrix for a H1(Curl) representation
!!
!! Supported boundary conditions
!! - \c 'none' Full matrix
!! - \c 'zerob' Dirichlet for all boundary DOF
!!
!! @param[in,out] mat Matrix object
!! @param[in] bc Boundary condition
!---------------------------------------------------------------------------
subroutine oft_hcurl_getkop(mat,bc)
class(oft_matrix), pointer, intent(inout) :: mat
character(LEN=*), intent(in) :: bc
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),cgop(3,6),elapsed_time
real(r8), allocatable :: rop(:,:),cop(:,:),kop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_hcurl_vec
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing H1(Curl)::KOP'
  CALL mytimer%tick()
END IF
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL oft_hcurl%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,rop,cop,det,kop,curved,goptmp,cgop,m,vol,jc,jr)
allocate(j(oft_hcurl%nce)) ! Local DOF and matrix indices
allocate(rop(3,oft_hcurl%nce)) ! Reconstructed gradient operator
allocate(cop(3,oft_hcurl%nce)) ! Reconstructed gradient operator
allocate(kop(oft_hcurl%nce,oft_hcurl%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  kop=0.d0
  do m=1,oft_hcurl%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)then
      call mesh%jacobian(i,oft_hcurl%quad%pts(:,m),goptmp,vol)
      call oft_hcurl_get_cgops(goptmp,cgop)
    end if
    det=vol*oft_hcurl%quad%wts(m)
    call oft_hcurl_eval_all(oft_hcurl,i,oft_hcurl%quad%pts(:,m),rop,goptmp)
    call oft_hcurl_ceval_all(oft_hcurl,i,oft_hcurl%quad%pts(:,m),cop,cgop)
    !---Compute local matrix contributions
    do jr=1,oft_hcurl%nce
      do jc=1,oft_hcurl%nce
        kop(jr,jc) = kop(jr,jc) + DOT_PRODUCT(rop(:,jr),cop(:,jc))*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call oft_hcurl%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,oft_hcurl%nce
        IF(oft_hcurl%global%gbe(j(jr)))kop(jr,:)=0.d0
      END DO
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,kop,oft_hcurl%nce,oft_hcurl%nce)
  ! !$omp end critical
end do
deallocate(j,rop,cop,kop)
!$omp end parallel
ALLOCATE(kop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
SELECT CASE(TRIM(bc))
  CASE("zerob")
    kop(1,1)=1.d0
    DO i=1,oft_hcurl%nbe
      jr=oft_hcurl%lbe(i)
      IF(.NOT.oft_hcurl%global%gbe(jr))CYCLE
      IF(.NOT.oft_hcurl%linkage%leo(i))CYCLE
      j=jr
      call mat%add_values(j,j,kop,1,1)
    END DO
END SELECT
DEALLOCATE(j,kop)
CALL oft_hcurl_create(oft_hcurl_vec)
CALL mat%assemble(oft_hcurl_vec)
CALL oft_hcurl_vec%delete
DEALLOCATE(oft_hcurl_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_hcurl_getkop
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_getwop
!---------------------------------------------------------------------------
!> Construct energy matrix for a H1(Curl) representation
!!
!! Supported boundary conditions
!! - \c 'none' Full matrix
!! - \c 'zerob' Dirichlet for all boundary DOF
!!
!! @param[in,out] mat Matrix object
!! @param[in] bc Boundary condition
!---------------------------------------------------------------------------
subroutine oft_hcurl_getwop(mat,bc)
class(oft_matrix), pointer, intent(inout) :: mat
character(LEN=*), intent(in) :: bc
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),cgop(3,6),elapsed_time
real(r8), allocatable :: cop(:,:),wop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_hcurl_vec
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing H1(Curl)::WOP'
  CALL mytimer%tick()
END IF
!---------------------------------------------------------------------------
! Allocate Laplacian Op
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL oft_hcurl%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,cop,det,wop,curved,goptmp,cgop,m,vol,jc,jr)
allocate(j(oft_hcurl%nce)) ! Local DOF and matrix indices
allocate(cop(3,oft_hcurl%nce)) ! Reconstructed gradient operator
allocate(wop(oft_hcurl%nce,oft_hcurl%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  wop=0.d0
  do m=1,oft_hcurl%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)then
      call mesh%jacobian(i,oft_hcurl%quad%pts(:,m),goptmp,vol)
      call oft_hcurl_get_cgops(goptmp,cgop)
    end if
    det=vol*oft_hcurl%quad%wts(m)
    call oft_hcurl_ceval_all(oft_hcurl,i,oft_hcurl%quad%pts(:,m),cop,cgop)
    !---Compute local matrix contributions
    do jr=1,oft_hcurl%nce
      do jc=1,oft_hcurl%nce
        wop(jr,jc) = wop(jr,jc) + DOT_PRODUCT(cop(:,jr),cop(:,jc))*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call oft_hcurl%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,oft_hcurl%nce
        IF(oft_hcurl%global%gbe(j(jr)))wop(jr,:)=0.d0
      END DO
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,wop,oft_hcurl%nce,oft_hcurl%nce)
  ! !$omp end critical
end do
deallocate(j,cop,wop)
!$omp end parallel
ALLOCATE(wop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
SELECT CASE(TRIM(bc))
  CASE("zerob")
    wop(1,1)=1.d0
    DO i=1,oft_hcurl%nbe
      jr=oft_hcurl%lbe(i)
      IF(.NOT.oft_hcurl%global%gbe(jr))CYCLE
      IF(.NOT.oft_hcurl%linkage%leo(i))CYCLE
      j=jr
      call mat%add_values(j,j,wop,1,1)
    END DO
END SELECT
DEALLOCATE(j,wop)
CALL oft_hcurl_create(oft_hcurl_vec)
CALL mat%assemble(oft_hcurl_vec)
CALL oft_hcurl_vec%delete
DEALLOCATE(oft_hcurl_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_hcurl_getwop
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_getjmlb
!---------------------------------------------------------------------------
!> Construct force-free response matrix for a H1(Curl) representation
!!
!! Supported boundary conditions
!! - \c 'none' Full matrix
!! - \c 'zerob' Dirichlet for all boundary DOF
!!
!! @param[in,out] mat Matrix object
!! @param[in] alam Lambda for response
!! @param[in] bc Boundary condition
!---------------------------------------------------------------------------
subroutine oft_hcurl_getjmlb(mat,alam,bc)
class(oft_matrix), pointer, intent(inout) :: mat
real(r8), intent(in) :: alam
character(LEN=*), intent(in) :: bc
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),cgop(3,6),elapsed_time
real(r8), allocatable :: rop(:,:),cop(:,:),wop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_hcurl_vec
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing H1(Curl)::JMLB'
  CALL mytimer%tick()
END IF
!---------------------------------------------------------------------------
! Allocate Laplacian Op
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL oft_hcurl%mat_create(mat)
ELSE
  CALL mat%zero
END IF
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,rop,cop,det,wop,curved,goptmp,cgop,m,vol,jc,jr)
allocate(j(oft_hcurl%nce)) ! Local DOF and matrix indices
allocate(rop(3,oft_hcurl%nce)) ! Reconstructed gradient operator
allocate(cop(3,oft_hcurl%nce)) ! Reconstructed gradient operator
allocate(wop(oft_hcurl%nce,oft_hcurl%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  wop=0.d0
  do m=1,oft_hcurl%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)then
      call mesh%jacobian(i,oft_hcurl%quad%pts(:,m),goptmp,vol)
      call oft_hcurl_get_cgops(goptmp,cgop)
    end if
    det=vol*oft_hcurl%quad%wts(m)
    call oft_hcurl_eval_all(oft_hcurl,i,oft_hcurl%quad%pts(:,m),rop,goptmp)
    call oft_hcurl_ceval_all(oft_hcurl,i,oft_hcurl%quad%pts(:,m),cop,cgop)
    !---Compute local matrix contributions
    do jr=1,oft_hcurl%nce
      do jc=1,oft_hcurl%nce
        wop(jr,jc) = wop(jr,jc) + (DOT_PRODUCT(cop(:,jr),cop(:,jc))-alam*DOT_PRODUCT(rop(:,jr),cop(:,jc)))*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call oft_hcurl%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,oft_hcurl%nce
        IF(oft_hcurl%global%gbe(j(jr)))wop(jr,:)=0.d0
      END DO
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,wop,oft_hcurl%nce,oft_hcurl%nce)
  ! !$omp end critical
end do
deallocate(j,rop,cop,wop)
!$omp end parallel
ALLOCATE(wop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
SELECT CASE(TRIM(bc))
  CASE("zerob")
    wop(1,1)=1.d0
    DO i=1,oft_hcurl%nbe
      jr=oft_hcurl%lbe(i)
      IF(.NOT.oft_hcurl%global%gbe(jr))CYCLE
      IF(.NOT.oft_hcurl%linkage%leo(i))CYCLE
      j=jr
      call mat%add_values(j,j,wop,1,1)
    END DO
END SELECT
DEALLOCATE(j,wop)
CALL oft_hcurl_create(oft_hcurl_vec)
CALL mat%assemble(oft_hcurl_vec)
CALL oft_hcurl_vec%delete
DEALLOCATE(oft_hcurl_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_hcurl_getjmlb
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_project
!---------------------------------------------------------------------------
!> Project a vector field onto a H1(Curl) basis
!!
!! @note This subroutine only performs the integration of the field with
!! test functions for a H1(Curl) basis. To retrieve the correct projection the
!! result must be multiplied by the inverse of H1(Curl)::MOP.
!!
!! @param[in,out] field Vector field for projection
!! @param[in,out] x Field projected onto H1(Curl) basis
!---------------------------------------------------------------------------
subroutine oft_hcurl_project(field,x)
class(fem_interp), intent(inout) :: field
class(oft_vector), intent(inout) :: x
!---
integer(i4) :: i,jc,m
integer(i4), allocatable, dimension(:) :: j
real(r8) :: det,vol,bcc(3),goptmp(3,4)
real(r8), pointer, dimension(:) :: xloc
real(r8), allocatable, dimension(:,:) :: rop
logical :: curved
DEBUG_STACK_PUSH
!---Initialize vectors to zero
NULLIFY(xloc)
call x%set(0.d0)
call x%get_local(xloc)
!---Integerate over the volume
!$omp parallel default(firstprivate) shared(xloc) private(curved,det)
allocate(j(oft_hcurl%nce),rop(3,oft_hcurl%nce))
!$omp do schedule(guided)
do i=1,mesh%nc ! Loop over cells
  call oft_hcurl%ncdofs(i,j) ! Get DOFs
  curved=cell_is_curved(mesh,i) ! Straight cell test
  do m=1,oft_hcurl%quad%np
    if(curved.OR.m==1)call mesh%jacobian(i,oft_hcurl%quad%pts(:,m),goptmp,vol)
    det=vol*oft_hcurl%quad%wts(m)
    call field%interp(i,oft_hcurl%quad%pts(:,m),goptmp,bcc)
    call oft_hcurl_eval_all(oft_hcurl,i,oft_hcurl%quad%pts(:,m),rop,goptmp)
    do jc=1,oft_hcurl%nce
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
!------------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_bcurl
!
!> Compute the boundary term for integration by parts of the curl operator using
!! a HCurl basis
!------------------------------------------------------------------------------
SUBROUTINE oft_hcurl_bcurl(field,x)
CLASS(fem_interp), INTENT(inout) :: field !< Vector field for projection
CLASS(oft_vector), INTENT(inout) :: x
!< \f$ \int \left( \textbf{u}^T \times \textbf{F} \right) \cdot \textbf{dS} \f$
INTEGER(i4) :: i,m,k,jc,cell,ptmap(3)
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: j
REAL(r8) :: vol,det,flog(4),norm(3),etmp(3),sgop(3,3),vgop(3,4)
REAL(r8), POINTER, DIMENSION(:) :: xloc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: rop
DEBUG_STACK_PUSH
!---Initialize vectors to zero
NULLIFY(xloc)
call x%set(0.d0)
call x%get_local(xloc)
!---Operator integration loop
det=0.d0
!$omp parallel default(firstprivate) shared(field,xloc) private(det)
allocate(j(oft_hcurl%nce),rop(3,oft_hcurl%nce))
!$omp do schedule(guided)
do i=1,smesh%nc
  CALL mesh%get_surf_map(i,cell,ptmap) ! Find parent cell and logical coordinate mapping
  call oft_hcurl%ncdofs(cell,j) ! Get local to global DOF mapping
  !---Loop over quadrature points
  do m=1,oft_bhcurl%quad%np
    call smesh%jacobian(i,oft_bhcurl%quad%pts(:,m),sgop,vol)
    call smesh%norm(i,oft_bhcurl%quad%pts(:,m),norm)
    det=vol*oft_bhcurl%quad%wts(m)
    !---Evaluate in cell coordinates
    CALL mesh%surf_to_vol(oft_bhcurl%quad%pts(:,m),ptmap,flog)
    call mesh%jacobian(cell,flog,vgop,vol)
    call field%interp(cell,flog,vgop,etmp)
    !---Project on to HCurl basis
    call oft_hcurl_eval_all(oft_hcurl,cell,flog,rop,vgop)
    do jc=1,oft_hcurl%nce
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
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_divout_setup
!---------------------------------------------------------------------------
!> Setup matrix and solver with default
!!
!! @note Should only be used via class \ref oft_hcurl_divout
!!
!! @param[in] bc Boundary condition
!---------------------------------------------------------------------------
subroutine hcurl_divout_setup(self,bc)
class(oft_hcurl_divout), intent(inout) :: self
character(LEN=*), intent(in) :: bc
DEBUG_STACK_PUSH
CALL create_cg_solver(self%solver)
self%solver%its=-3
CALL oft_lag_getlop(self%solver%A,bc)
CALL create_diag_pre(self%solver%pre)
IF(TRIM(bc)=='grnd')THEN
  self%bc=>lag_zerogrnd
ELSE
  self%bc=>lag_zerob
END IF
DEBUG_STACK_POP
end subroutine hcurl_divout_setup
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_divout_apply
!---------------------------------------------------------------------------
!> Remove divergence from a H1(Curl) vector field by adding a gradient correction
!!
!! @note Should only be used via class \ref oft_hcurl_divout
!!
!! @param[in,out] u Field for divergence cleaning
!---------------------------------------------------------------------------
subroutine hcurl_divout_apply(self,u)
class(oft_hcurl_divout), intent(inout) :: self
class(oft_vector), intent(inout) :: u
class(oft_vector), pointer :: a,g,tmp,tmp2
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
call oft_lag_create(g)
call oft_lag_create(a)
!---
IF(ASSOCIATED(self%mop))THEN
  CALL hcurl_gradtp(u,g)
ELSE
  CALL hcurl_div(u,g)
END IF
uu=u%dot(u)
self%solver%atol=MAX(self%solver%atol,SQRT(uu*1.d-20))
call a%set(0.d0)
call self%bc(g)
!---
pm_save=oft_env%pm; oft_env%pm=self%pm
call self%solver%apply(a,g)
oft_env%pm=pm_save
!---
CALL a%scale(-1.d0)
CALL u%new(tmp)
CALL hcurl_grad(a,tmp)
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
deallocate(tmp,a,g)
DEBUG_STACK_POP
end subroutine hcurl_divout_apply
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_divout_delete
!---------------------------------------------------------------------------
!> Clean-up internal storage for a oft_hcurl_divout object
!!
!! @note Should only be used via class \ref oft_hcurl_divout
!---------------------------------------------------------------------------
subroutine hcurl_divout_delete(self)
class(oft_hcurl_divout), intent(inout) :: self
end subroutine hcurl_divout_delete
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_orthog_apply
!---------------------------------------------------------------------------
!> Orthogonalize a H1(Curl) vector against a library of given modes.
!!
!! @note Used as a member function of oft_hcurl_orthog only
!!
!! @param[in,out] u Field to orthogonalize
!---------------------------------------------------------------------------
subroutine hcurl_orthog_apply(self,u)
class(oft_hcurl_orthog), intent(inout) :: self
class(oft_vector), intent(inout) :: u
class(oft_vector), pointer :: b
real(r8) :: c
integer(i4) :: i
DEBUG_STACK_PUSH
!---Get temporary variable
call u%new(b)
do i=1,self%nm
  !---Compute coupling
  call self%wop%apply(self%orthog(i,oft_hcurl_level)%f,b)
  c=u%dot(b)
  !---Remove coupling
  call u%add(1.d0,-c,self%orthog(i,oft_hcurl_level)%f)
end do
!---Delete temporary variable
call b%delete
deallocate(b)
DEBUG_STACK_POP
end subroutine hcurl_orthog_apply
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_orthog_delete
!---------------------------------------------------------------------------
!> Clean-up internal storage for a oft_hcurl_orthog object
!!
!! @note Used as a member function of oft_hcurl_orthog only
!---------------------------------------------------------------------------
subroutine hcurl_orthog_delete(self)
class(oft_hcurl_orthog), intent(inout) :: self
end subroutine hcurl_orthog_delete
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_setup_interp
!---------------------------------------------------------------------------
!> Construct interpolation matrices on each MG level
!---------------------------------------------------------------------------
SUBROUTINE hcurl_setup_interp
INTEGER(i4) :: i
DEBUG_STACK_PUSH
!---
DO i=oft_hcurl_minlev+1,oft_hcurl_nlevels
  CALL oft_hcurl_set_level(i)
  !---
  if(oft_hcurl_level==oft_hcurl_blevel+1)then
    CYCLE
  end if
  !---Setup interpolation
  if(oft_hcurl%order==1)then
    CALL hcurl_ginterpmatrix(ML_oft_hcurl%interp_matrices(ML_oft_hcurl%level)%m)
    oft_hcurl_ops%interp=>ML_oft_hcurl%interp_matrices(ML_oft_hcurl%level)%m
    CALL oft_hcurl_ops%interp%assemble
  else
    CALL hcurl_pinterpmatrix(ML_oft_hcurl%interp_matrices(ML_oft_hcurl%level)%m)
    oft_hcurl_ops%interp=>ML_oft_hcurl%interp_matrices(ML_oft_hcurl%level)%m
    CALL oft_hcurl_ops%interp%assemble
  end if
END DO
DEBUG_STACK_POP
END SUBROUTINE hcurl_setup_interp
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_ginterpmatrix
!---------------------------------------------------------------------------
!> Construct interpolation matrix for polynomial levels
!---------------------------------------------------------------------------
SUBROUTINE hcurl_ginterpmatrix(mat)
class(oft_matrix), pointer, intent(inout) :: mat
INTEGER(i4) :: i,j,k,m,icors,ifine,jb,i_ind(1),j_ind(1)
INTEGER(i4) :: etmp(2),ftmp(3),fetmp(3),ctmp(4),fc,ed
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pmap,emap,fmap
CLASS(oft_hcurl_fem), POINTER :: hcurl_cors => NULL()
! TYPE(oft_fem_type), POINTER :: hcurl_cors => NULL()
TYPE(oft_nedelec_ops), POINTER :: ops
class(oft_mesh), pointer :: cmesh
CLASS(oft_vector), POINTER :: hcurl_vec,hcurl_vec_cors
integer(i4) :: jfe(3),jce(6)
integer(i4), pointer :: lcdg(:),lfde(:,:),lede(:,:),lcde(:,:)
real(r8) :: f(4),incr,val,d(3),goptmp(3,4),v,mop(1),h_rop(3,6)
type(oft_graph_ptr), pointer :: graphs(:,:)
DEBUG_STACK_PUSH
!---
if(mg_mesh%level<1)then
  call oft_abort('Invalid mesh level','hcurl_ginterpmatrix',__FILE__)
end if
cmesh=>mg_mesh%meshes(mg_mesh%level-1)
if(cmesh%type/=1)CALL oft_abort("Only supported with tet meshes", &
  "hcurl_ginterpmatrix", __FILE__)
if(oft_hcurl%order/=1)then
  call oft_abort('Attempted geometric interpolation for pd > 1', &
    'hcurl_ginterpmatrix',__FILE__)
end if
ops=>oft_hcurl_ops
SELECT TYPE(this=>ML_oft_hcurl%levels(oft_hcurl_level-1)%fe)
  CLASS IS(oft_hcurl_fem)
    hcurl_cors=>this
  CLASS DEFAULT
    CALL oft_abort("Error casting coarse level", "hcurl_ginterpmatrix", __FILE__)
END SELECT
! hcurl_cors=>ML_oft_hcurl%levels(oft_hcurl_level-1)%fe
lede=>mg_mesh%inter(mg_mesh%level-1)%lede
lfde=>mg_mesh%inter(mg_mesh%level-1)%lfde
lcdg=>mg_mesh%inter(mg_mesh%level-1)%lcdg
lcde=>mg_mesh%inter(mg_mesh%level-1)%lcde
ALLOCATE(ML_oft_hcurl%interp_graphs(ML_oft_hcurl%level)%g)
ops%interp_graph=>ML_oft_hcurl%interp_graphs(ML_oft_hcurl%level)%g
!---Setup matrix sizes
ops%interp_graph%nr=oft_hcurl%ne
ops%interp_graph%nrg=oft_hcurl%global%ne
ops%interp_graph%nc=hcurl_cors%ne
ops%interp_graph%ncg=hcurl_cors%global%ne
!---Setup Matrix graph
ALLOCATE(ops%interp_graph%kr(ops%interp_graph%nr+1))
ops%interp_graph%nnz=2*cmesh%ne+9*cmesh%nf+6*cmesh%nc
ALLOCATE(ops%interp_graph%lc(ops%interp_graph%nnz))
ops%interp_graph%lc=0_i4
!---Construct linkage
DO i=1,cmesh%ne
  ops%interp_graph%kr(lede(1,i))=1
  ops%interp_graph%kr(lede(2,i))=1
END DO
DO i=1,cmesh%nf
  jfe=ABS(cmesh%lfe(:,i)) ! face edges
  ops%interp_graph%kr(ABS(lfde(1,i)))=3
  ops%interp_graph%kr(ABS(lfde(2,i)))=3
  ops%interp_graph%kr(ABS(lfde(3,i)))=3
END DO
DO i=1,cmesh%nc
  ops%interp_graph%kr(ABS(lcde(1,i)))=6
END DO
ops%interp_graph%kr(ops%interp_graph%nr+1)=ops%interp_graph%nnz+1
do i=ops%interp_graph%nr,1,-1 ! cumulative point to point count
  ops%interp_graph%kr(i)=ops%interp_graph%kr(i+1)-ops%interp_graph%kr(i)
end do
if(ops%interp_graph%kr(1)/=1)call oft_abort('Bad element to element count','oft_hcurl_ginterpmatrix',__FILE__)
DO i=1,cmesh%ne
  k=ops%interp_graph%kr(lede(1,i))
  ops%interp_graph%lc(k)=i
  k=ops%interp_graph%kr(lede(2,i))
  ops%interp_graph%lc(k)=i
END DO
DO i=1,cmesh%nf
  jfe=ABS(cmesh%lfe(:,i)) ! face edges
  CALL sort_array(jfe,3)
  !---
  DO j=1,3
    k=ops%interp_graph%kr(ABS(lfde(j,i)))
    ops%interp_graph%lc(k)=jfe(1)
    ops%interp_graph%lc(k+1)=jfe(2)
    ops%interp_graph%lc(k+2)=jfe(3)
  END DO
END DO
DO i=1,cmesh%nc
  jce=ABS(cmesh%lce(:,i)) ! face edges
  CALL sort_array(jce,6)
  !---
  k=ops%interp_graph%kr(ABS(lcde(1,i)))
  DO j=0,5
    ops%interp_graph%lc(k+j)=jce(1+j)
  END DO
END DO
!---------------------------------------------------------------------------
! Construct matrix
!---------------------------------------------------------------------------
CALL oft_hcurl_create(hcurl_vec)
CALL oft_hcurl_create(hcurl_vec_cors,oft_hcurl_level-1)
!---
ALLOCATE(graphs(1,1))
graphs(1,1)%g=>ops%interp_graph
!---
CALL create_matrix(mat,graphs,hcurl_vec,hcurl_vec_cors)
CALL hcurl_vec%delete
CALL hcurl_vec_cors%delete
DEALLOCATE(graphs,hcurl_vec,hcurl_vec_cors)
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
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
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_pinterpmatrix
!---------------------------------------------------------------------------
!> Construct interpolation matrix for polynomial levels
!---------------------------------------------------------------------------
SUBROUTINE hcurl_pinterpmatrix(mat)
class(oft_matrix), pointer, intent(inout) :: mat
INTEGER(i4) :: i,j,k,m,icors,ifine,jb,js,jn,i_ind(1),j_ind(1)
INTEGER(i4) :: etmp(2),ftmp(3),fetmp(3),ctmp(4),fc,ed
INTEGER(i4) :: offsetc,offsetf
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pmap,emap,fmap
REAL(r8) :: f(4),val,mop(1)
TYPE(oft_fem_type), POINTER :: hcurl_cors => NULL()
TYPE(oft_nedelec_ops), POINTER :: ops
CLASS(oft_vector), POINTER :: hcurl_vec,hcurl_vec_cors
type(oft_graph_ptr), pointer :: graphs(:,:)
DEBUG_STACK_PUSH
!---
ops=>oft_hcurl_ops
SELECT TYPE(this=>ML_oft_hcurl%levels(oft_hcurl_level-1)%fe)
CLASS IS(oft_fem_type)
  hcurl_cors=>this
CLASS DEFAULT
  CALL oft_abort("Error getting coarse FE object","hcurl_pinterpmatrix",__FILE__)
END SELECT
ALLOCATE(ML_oft_hcurl%interp_graphs(ML_oft_hcurl%level)%g)
ops%interp_graph=>ML_oft_hcurl%interp_graphs(ML_oft_hcurl%level)%g
!---Setup matrix sizes
ops%interp_graph%nr=oft_hcurl%ne
ops%interp_graph%nrg=oft_hcurl%global%ne
ops%interp_graph%nc=hcurl_cors%ne
ops%interp_graph%ncg=hcurl_cors%global%ne
!---Setup Matrix graph
ALLOCATE(ops%interp_graph%kr(ops%interp_graph%nr+1))
ops%interp_graph%kr=0
ops%interp_graph%nnz=hcurl_cors%ne
ALLOCATE(ops%interp_graph%lc(ops%interp_graph%nnz))
ops%interp_graph%lc=0_i4
!---Construct matrix
do i=1,mesh%ne
  do j=1,hcurl_cors%gstruct(2)
    offsetf=(i-1)*oft_hcurl%gstruct(2)
    ops%interp_graph%kr(j+offsetf)=1
  end do
end do
!---
do i=1,mesh%nf
  do j=1,hcurl_cors%gstruct(3)
    offsetf=mesh%ne*oft_hcurl%gstruct(2)+(i-1)*oft_hcurl%gstruct(3)
    ops%interp_graph%kr(j+offsetf)=1
  end do
end do
!---
do i=1,mesh%nc
  do j=1,hcurl_cors%gstruct(4)
    offsetf=mesh%ne*oft_hcurl%gstruct(2)+oft_hcurl%gstruct(3)*mesh%nf+(i-1)*oft_hcurl%gstruct(4)
    ops%interp_graph%kr(j+offsetf)=1
  end do
end do
ops%interp_graph%kr(ops%interp_graph%nr+1)=ops%interp_graph%nnz+1
do i=ops%interp_graph%nr,1,-1 ! cumulative point to point count
  ops%interp_graph%kr(i)=ops%interp_graph%kr(i+1)-ops%interp_graph%kr(i)
end do
if(ops%interp_graph%kr(1)/=1)call oft_abort('Bad element to element count','oft_hcurl_pinterpmatrix',__FILE__)
!---Construct matrix
do i=1,mesh%ne
  do j=1,hcurl_cors%gstruct(2)
    offsetf=(i-1)*oft_hcurl%gstruct(2)
    offsetc=(i-1)*hcurl_cors%gstruct(2)
    ops%interp_graph%lc(ops%interp_graph%kr(j+offsetf))=j+offsetc
  end do
end do
!---
do i=1,mesh%nf
  do j=1,hcurl_cors%gstruct(3)
    offsetf=mesh%ne*oft_hcurl%gstruct(2)+(i-1)*oft_hcurl%gstruct(3)
    offsetc=mesh%ne*hcurl_cors%gstruct(2)+(i-1)*hcurl_cors%gstruct(3)
    ops%interp_graph%lc(ops%interp_graph%kr(j+offsetf))=j+offsetc
  end do
end do
!---
do i=1,mesh%nc
  do j=1,hcurl_cors%gstruct(4)
    offsetf=mesh%ne*oft_hcurl%gstruct(2)+oft_hcurl%gstruct(3)*mesh%nf+(i-1)*oft_hcurl%gstruct(4)
    offsetc=mesh%ne*hcurl_cors%gstruct(2)+hcurl_cors%gstruct(3)*mesh%nf+(i-1)*hcurl_cors%gstruct(4)
    ops%interp_graph%lc(ops%interp_graph%kr(j+offsetf))=j+offsetc
  end do
end do
!---------------------------------------------------------------------------
! Construct matrix
!---------------------------------------------------------------------------
CALL oft_hcurl_create(hcurl_vec)
CALL oft_hcurl_create(hcurl_vec_cors,oft_hcurl_level-1)
!---
ALLOCATE(graphs(1,1))
graphs(1,1)%g=>ops%interp_graph
!---
CALL create_matrix(mat,graphs,hcurl_vec,hcurl_vec_cors)
CALL hcurl_vec%delete
CALL hcurl_vec_cors%delete
DEALLOCATE(graphs,hcurl_vec,hcurl_vec_cors)
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---
allocate(emap(mesh%ne))
CALL get_inverse_map(mesh%lbe,mesh%nbe,emap,mesh%ne)
DO i=1,mesh%ne
  IF(mesh%be(i))THEN
    ! IF(.NOT.mesh%linkage%leo(emap(i)))CYCLE
    IF(.NOT.mesh%estitch%leo(emap(i)))CYCLE
  END IF
  DO j=1,hcurl_cors%gstruct(2)
    i_ind=j+(i-1)*oft_hcurl%gstruct(2)
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
    i_ind=j+mesh%ne*oft_hcurl%gstruct(2)+(i-1)*oft_hcurl%gstruct(3)
    j_ind=j+mesh%ne*hcurl_cors%gstruct(2)+(i-1)*hcurl_cors%gstruct(3)
    mop=1.d0
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
  END DO
END DO
deallocate(fmap)
!---
DO i=1,mesh%nc
  DO j=1,hcurl_cors%gstruct(4)
    i_ind=j+mesh%ne*oft_hcurl%gstruct(2)+oft_hcurl%gstruct(3)*mesh%nf+(i-1)*oft_hcurl%gstruct(4)
    j_ind=j+mesh%ne*hcurl_cors%gstruct(2)+hcurl_cors%gstruct(3)*mesh%nf+(i-1)*hcurl_cors%gstruct(4)
    mop=1.d0
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
  END DO
END DO
DEBUG_STACK_POP
END SUBROUTINE hcurl_pinterpmatrix
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_interp
!---------------------------------------------------------------------------
!> Interpolate a coarse level H1(Curl) vector field to the next finest level
!!
!! @note The global H1(Curl) level in incremented by one in this subroutine
!!
!! @param[in] acors Vector to interpolate
!! @param[in,out] afine Fine vector from interpolation
!---------------------------------------------------------------------------
subroutine hcurl_interp(acors,afine)
class(oft_vector), intent(inout) :: acors
class(oft_vector), intent(inout) :: afine
DEBUG_STACK_PUSH
!---Step one level up
call oft_hcurl_set_level(oft_hcurl_level+1)
call afine%set(0.d0)
!---
if(oft_hcurl_level==oft_hcurl_blevel+1)then
  call hcurl_base_pop(acors,afine)
  DEBUG_STACK_POP
  return
end if
CALL oft_hcurl_ops%interp%apply(acors,afine)
DEBUG_STACK_POP
end subroutine hcurl_interp
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_base_pop
!---------------------------------------------------------------------------
!> Transfer a base level H1(Curl) vector field to the next MPI level
!!
!! @param[in] acors Vector to transfer
!! @param[in,out] afine Fine vector from transfer
!---------------------------------------------------------------------------
subroutine hcurl_base_pop(acors,afine)
class(oft_vector), intent(inout) :: acors
class(oft_vector), intent(inout) :: afine
integer(i4), pointer :: lbege(:)
integer(i4) :: i
real(r8), pointer, dimension(:) :: array_c,array_f
DEBUG_STACK_PUSH
lbege=>mg_mesh%inter(mg_mesh%nbase)%lbege
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
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_inject
!---------------------------------------------------------------------------
!> Inject a fine level H1(Curl) vector field to the next coarsest level
!!
!! @note The global H1(Curl) level in decremented by one in this subroutine
!!
!! @param[in] afine Vector to inject
!! @param[in,out] acors Coarse vector from injection
!---------------------------------------------------------------------------
subroutine hcurl_inject(afine,acors)
class(oft_vector), intent(inout) :: afine
class(oft_vector), intent(inout) :: acors
integer(i4) :: i,j,k
logical :: gcheck
DEBUG_STACK_PUSH
gcheck=(oft_hcurl%order==1)
! Step down level up
call oft_hcurl_set_level(oft_hcurl_level-1)
! Cast fine field
call acors%set(0.d0)
if(oft_hcurl_level==oft_hcurl_blevel)then
  call hcurl_base_push(afine,acors)
  DEBUG_STACK_POP
  return
end if
CALL ML_oft_hcurl_ops(oft_hcurl_level+1)%interp%applyT(afine,acors)
DEBUG_STACK_POP
end subroutine hcurl_inject
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_base_push
!---------------------------------------------------------------------------
!> Transfer a MPI level H1(Curl) vector field to the base level
!!
!! @param[in] afine Vector to transfer
!! @param[in,out] acors Fine vector from transfer
!---------------------------------------------------------------------------
subroutine hcurl_base_push(afine,acors)
class(oft_vector), intent(inout) :: afine
class(oft_vector), intent(inout) :: acors
integer(i4), pointer :: lbege(:)
integer(i4) :: i,j,ierr
real(r8), pointer, dimension(:) :: alias,array_c,array_f
CLASS(oft_afem_type), POINTER :: hcurl_fine => NULL()
DEBUG_STACK_PUSH
!---
lbege=>mg_mesh%inter(mg_mesh%nbase)%lbege
CALL acors%get_local(array_c)
CALL afine%get_local(array_f)
hcurl_fine=>ML_oft_hcurl%levels(oft_hcurl_level+1)%fe
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
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_wop_eigs
!---------------------------------------------------------------------------
!> Compute eigenvalues and smoothing coefficients for the operator H1(Curl)::WOP
!---------------------------------------------------------------------------
SUBROUTINE hcurl_wop_eigs(minlev)
INTEGER(i4), INTENT(in) :: minlev
#ifdef HAVE_ARPACK
INTEGER(i4) :: i
REAL(r8) :: lam0
REAL(r8), ALLOCATABLE :: df(:)
CLASS(oft_vector), POINTER :: u
TYPE(oft_irlm_eigsolver) :: arsolver
CLASS(oft_matrix), POINTER :: md => NULL()
CLASS(oft_matrix), POINTER :: wop => NULL()
DEBUG_STACK_PUSH
!---------------------------------------------------------------------------
! Compute optimal smoother coefficients
!---------------------------------------------------------------------------
IF(oft_env%head_proc)WRITE(*,*)'Optimizing Jacobi damping for H1(Curl)::WOP'
ALLOCATE(df(oft_hcurl_nlevels))
df=0.d0
DO i=minlev,oft_hcurl_nlevels
  CALL oft_hcurl_set_level(i)
  !---Create fields
  CALL oft_hcurl_create(u)
  !---Get Ev range
  NULLIFY(wop)
  CALL oft_hcurl_getwop(wop,'zerob')
  CALL create_diagmatrix(md,wop%D)
  !---
  arsolver%A=>wop
  arsolver%M=>md
  arsolver%mode=2
  arsolver%tol=1.E-5_r8
  arsolver%bc=>hcurl_zerob
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
  DO i=1,oft_hcurl_nlevels-1
    WRITE(*,'(F5.3,A)',ADVANCE='NO')df(i),', '
  END DO
  WRITE(*,'(F5.3,A)')df(oft_hcurl_nlevels)
END IF
DEALLOCATE(df)
DEBUG_STACK_POP
#else
CALL oft_abort("Subroutine requires ARPACK", "lag_lop_eigs", __FILE__)
#endif
END SUBROUTINE hcurl_wop_eigs
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_getwop_pre
!---------------------------------------------------------------------------
!> Construct default MG preconditioner for H1(Curl)::WOP
!---------------------------------------------------------------------------
SUBROUTINE hcurl_getwop_pre(pre,mats,level,nlevels)
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
INTEGER(i4) :: i,j,levin
LOGICAL :: create_mats
CHARACTER(LEN=2) :: lev_char
!---
TYPE(fox_node), POINTER :: pre_node
#ifdef HAVE_XML
integer(i4) :: nnodes
TYPE(fox_node), POINTER :: hcurl_node
TYPE(fox_nodelist), POINTER :: current_nodes
#endif
DEBUG_STACK_PUSH
!---
minlev=1
toplev=oft_hcurl_level
levin=oft_hcurl_level
IF(PRESENT(level))toplev=level
IF(PRESENT(nlevels))minlev=toplev-nlevels+1
nl=toplev-minlev+1
!---
IF(minlev<1)CALL oft_abort('Minimum level is < 0','hcurl_getwop_pre',__FILE__)
IF(toplev>oft_hcurl_nlevels)CALL oft_abort('Maximum level is > hcurl_nlevels','hcurl_getwop_pre',__FILE__)
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
  CALL oft_hcurl_set_level(minlev+(i-1))
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
    CALL oft_hcurl_getwop(mats(i)%M,'zerob')
  END IF
  IF(i>1)ml_int(i-1)%M=>oft_hcurl_ops%interp
END DO
CALL oft_hcurl_set_level(levin)
!---------------------------------------------------------------------------
! Search for XML-spec
!---------------------------------------------------------------------------
NULLIFY(pre_node)
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml))THEN
  !---Look for Lagrange node
  current_nodes=>fox_getElementsByTagName(oft_env%xml,"nedelec_h1curl")
  nnodes=fox_getLength(current_nodes)
  IF(nnodes>0)THEN
    hcurl_node=>fox_item(current_nodes,0)
    !---Look for lop node
    current_nodes=>fox_getElementsByTagName(hcurl_node,"wop")
    nnodes=fox_getLength(current_nodes)
    IF(nnodes>0)pre_node=>fox_item(current_nodes,0)
  END IF
END IF
#endif
!---------------------------------------------------------------------------
! Setup preconditioner
!---------------------------------------------------------------------------
NULLIFY(pre)
CALL create_mlpre(pre,mats(1:nl),levels,nlevels=nl,create_vec=oft_hcurl_create, &
     interp=hcurl_interp,inject=hcurl_inject,bc=hcurl_zerob,stype=1,df=df,nu=nu,xml_root=pre_node)
!---------------------------------------------------------------------------
! Cleanup
!---------------------------------------------------------------------------
DEALLOCATE(ml_int,levels,df,nu)
DEBUG_STACK_POP
END SUBROUTINE hcurl_getwop_pre
!---------------------------------------------------------------------------
! SUBROUTINE: hcurl_getjmlb_pre
!---------------------------------------------------------------------------
!> Construct default MG preconditioner for H1(Curl)::JMLB
!---------------------------------------------------------------------------
SUBROUTINE hcurl_getjmlb_pre(pre,mats,alam,level,nlevels)
CLASS(oft_solver), POINTER, INTENT(out) :: pre
TYPE(oft_matrix_ptr), POINTER, INTENT(inout) :: mats(:)
REAL(r8), INTENT(in) :: alam
INTEGER(i4), OPTIONAL, INTENT(in) :: level
INTEGER(i4), OPTIONAL, INTENT(in) :: nlevels
!--- ML structures for MG-preconditioner
TYPE(oft_matrix_ptr), POINTER :: ml_int(:)
INTEGER(i4), ALLOCATABLE :: levels(:)
INTEGER(i4), ALLOCATABLE :: nu(:)
!---
INTEGER(i4) :: minlev,toplev,nl
INTEGER(i4) :: i,j,levin
LOGICAL :: create_mats
CHARACTER(LEN=2) :: lev_char
CLASS(oft_vector), POINTER :: oft_hcurl_vec
!---
TYPE(fox_node), POINTER :: pre_node
#ifdef HAVE_XML
integer(i4) :: nnodes
TYPE(fox_node), POINTER :: hcurl_node
TYPE(fox_nodelist), POINTER :: current_nodes
#endif
DEBUG_STACK_PUSH
!---
minlev=1
toplev=oft_hcurl_level
levin=oft_hcurl_level
IF(PRESENT(level))toplev=level
IF(PRESENT(nlevels))minlev=toplev-nlevels+1
nl=toplev-minlev+1
!---
IF(minlev<1)CALL oft_abort('Minimum level is < 0','hcurl_getjmlb_pre',__FILE__)
IF(toplev>oft_hcurl_nlevels)CALL oft_abort('Maximum level is > hcurl_nlevels','hcurl_getjmlb_pre',__FILE__)
!---------------------------------------------------------------------------
! Create ML Matrices
!---------------------------------------------------------------------------
create_mats=.FALSE.
IF(.NOT.ASSOCIATED(mats))THEN
  create_mats=.TRUE.
  ALLOCATE(mats(nl))
END IF
ALLOCATE(ml_int(nl-1),levels(nl),nu(nl))
DO i=1,nl
  CALL oft_hcurl_set_level(minlev+(i-1))
  levels(i)=minlev+(i-1)
  nu(i)=nu_jmlb(levels(i))
  !---
  IF(create_mats)THEN
    NULLIFY(mats(i)%M)
    CALL oft_hcurl_getjmlb(mats(i)%M,alam,'zerob')
  END IF
  IF(i>1)ml_int(i-1)%M=>oft_hcurl_ops%interp
END DO
CALL oft_hcurl_set_level(levin)
!---------------------------------------------------------------------------
! Search for XML-spec
!---------------------------------------------------------------------------
NULLIFY(pre_node)
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml))THEN
  !---Look for Lagrange node
  current_nodes=>fox_getElementsByTagName(oft_env%xml,"nedelec_h1curl")
  nnodes=fox_getLength(current_nodes)
  IF(nnodes>0)THEN
    hcurl_node=>fox_item(current_nodes,0)
    !---Look for lop node
    current_nodes=>fox_getElementsByTagName(hcurl_node,"jmlb")
    nnodes=fox_getLength(current_nodes)
    IF(nnodes>0)pre_node=>fox_item(current_nodes,0)
  END IF
END IF
#endif
!---------------------------------------------------------------------------
! Setup preconditioner
!---------------------------------------------------------------------------
NULLIFY(pre)
CALL create_mlpre(pre,mats(1:nl),levels,nlevels=nl,create_vec=oft_hcurl_create, &
     interp=hcurl_interp,inject=hcurl_inject,bc=hcurl_zerob,stype=2,nu=nu,xml_root=pre_node)
!---------------------------------------------------------------------------
! Cleanup
!---------------------------------------------------------------------------
DEALLOCATE(ml_int,levels,nu)
DEBUG_STACK_POP
END SUBROUTINE hcurl_getjmlb_pre
end module oft_hcurl_operators
