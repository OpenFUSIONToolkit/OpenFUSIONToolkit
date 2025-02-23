!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_h0_operators.F90
!
!> Nedelec H0 FE operator definitions
!! - Operator construction
!!   - MOP: mass matrix        \f$ \int \left( u^T v \right) dV \f$
!!   - LOP: laplacian matrix   \f$ \int \left( \nabla u^T \cdot \nabla v \right) dV \f$
!! - Interpolation classes
!! - Boundary conditions
!! - Multi-Grid setup and operators
!! - Default preconditioner setup
!!   - Optimal smoother calculation
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_h0
!---------------------------------------------------------------------------
MODULE oft_h0_operators
USE oft_base
USE oft_mesh_type, ONLY: oft_mesh, mesh, cell_is_curved
USE multigrid, ONLY: mg_mesh
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr, oft_graph_ptr
USE oft_deriv_matrices, ONLY: oft_diagmatrix, create_diagmatrix
USE oft_solver_base, ONLY: oft_solver
USE oft_la_utils, ONLY: create_matrix
USE oft_solver_utils, ONLY: create_mlpre, create_native_pre
#ifdef HAVE_ARPACK
USE oft_arpack, ONLY: oft_irlm_eigsolver
#endif
USE fem_base, ONLY: oft_afem_type, oft_fem_type, fem_max_levels
USE fem_utils, ONLY: fem_interp
USE oft_h0_basis, ONLY: oft_h0, ML_oft_h0, oft_h0_level, oft_h0_blevel, oft_h0_nlevels, &
h0_ops, oft_h0_ops, ML_oft_h0_ops, oft_h0_eval_all, oft_h0_geval_all, oft_h0_set_level, &
oft_h0_minlev, oft_h0_fem
USE oft_h0_fields, ONLY: oft_h0_create
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------
! CLASS fem_interp
!---------------------------------------------------------------------------
!> Interpolate a H0 field
!---------------------------------------------------------------------------
type, extends(fem_interp) :: oft_h0_rinterp
  class(oft_vector), pointer :: u => NULL() !< Field to interpolate
  real(r8), pointer, dimension(:) :: vals => NULL() !< Local values
  class(oft_h0_fem), pointer :: h0_rep => NULL() !< FE representation
contains
  !> Retrieve local values for interpolation
  procedure :: setup => h0_rinterp_setup
  !> Reconstruct field
  procedure :: interp => h0_rinterp
  !> Destroy temporary internal storage
  procedure :: delete => h0_rinterp_delete
end type oft_h0_rinterp
!---------------------------------------------------------------------------
! CLASS oft_h0_ginterp
!---------------------------------------------------------------------------
!> Interpolate \f$ \nabla \f$ of a H0 field
!---------------------------------------------------------------------------
type, extends(oft_h0_rinterp) :: oft_h0_ginterp
contains
  !> Reconstruct field
  procedure :: interp => h0_ginterp_apply
end type oft_h0_ginterp
!---Pre options
integer(i4) :: nu_lop(fem_max_levels)=0
real(r8) :: df_lop(fem_max_levels)=-1.d99
!---Cache variables
REAL(r8), POINTER, DIMENSION(:) :: oft_h0_rop => NULL()
REAL(r8), POINTER, DIMENSION(:,:) :: oft_h0_gop => NULL()
!$omp threadprivate(oft_h0_rop,oft_h0_gop)
contains
!---------------------------------------------------------------------------
! SUBROUTINE: h0_mloptions
!---------------------------------------------------------------------------
!> Read-in options for the basic Nedelec H0 ML preconditioners
!---------------------------------------------------------------------------
subroutine h0_mloptions
integer(i4) :: ierr,io_unit
namelist/h0_op_options/df_lop,nu_lop
DEBUG_STACK_PUSH
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,h0_op_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr>0)THEN
  CALL oft_abort('Error reading "h0_op_options" group in input file','h0_mloptions',__FILE__)
END IF
IF(df_lop(1)<-1.d90)THEN
  IF(oft_env%head_proc)THEN
    WRITE(*,*)'No H0/H(Grad) MG smoother settings found:'
    WRITE(*,*)'  Using default values, which may result in convergence failure.'
  END IF
  nu_lop=2
  df_lop=.2d0
END IF
DEBUG_STACK_POP
end subroutine h0_mloptions
!---------------------------------------------------------------------------
! SUBROUTINE: h0_rinterp_setup
!---------------------------------------------------------------------------
!> Setup interpolator for H0 scalar fields
!!
!! Fetches local representation used for interpolation from vector object
!!
!! @note Should only be used via class \ref oft_h0_rinterp or children
!---------------------------------------------------------------------------
subroutine h0_rinterp_setup(self)
class(oft_h0_rinterp), intent(inout) :: self
!---Get local slice
CALL self%u%get_local(self%vals)
self%h0_rep=>oft_h0
end subroutine h0_rinterp_setup
!---------------------------------------------------------------------------
! SUBROUTINE: h0_rinterp_delete
!---------------------------------------------------------------------------
!> Destroy temporary internal storage
!!
!! @note Should only be used via class \ref oft_h0_rinterp or children
!---------------------------------------------------------------------------
subroutine h0_rinterp_delete(self)
class(oft_h0_rinterp), intent(inout) :: self
!---Deallocate local storage
IF(ASSOCIATED(self%vals))DEALLOCATE(self%vals)
NULLIFY(self%h0_rep,self%u)
end subroutine h0_rinterp_delete
!---------------------------------------------------------------------------
! SUBROUTINE: h0_rinterp
!---------------------------------------------------------------------------
!> Reconstruct a Nedelec H0 scalar field
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [1]
!---------------------------------------------------------------------------
subroutine h0_rinterp(self,cell,f,gop,val)
class(oft_h0_rinterp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8), allocatable :: rop(:)
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!','h0_rinterp',__FILE__)
!---Get dofs
allocate(j(self%h0_rep%nce),rop(self%h0_rep%nce))
call self%h0_rep%ncdofs(cell,j) ! get DOFs
!---Reconstruct field
CALL oft_h0_eval_all(self%h0_rep,cell,f,rop)
val=0.d0
do jc=1,self%h0_rep%nce
  val=val+self%vals(j(jc))*rop(jc)
end do
deallocate(rop,j)
DEBUG_STACK_POP
end subroutine h0_rinterp
!---------------------------------------------------------------------------
! SUBROUTINE: h0_ginterp_apply
!---------------------------------------------------------------------------
!> Reconstruct the gradient of a Nedelec H0 scalar field
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed gradient at f [3]
!---------------------------------------------------------------------------
subroutine h0_ginterp_apply(self,cell,f,gop,val)
class(oft_h0_ginterp), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8), allocatable :: rop(:,:)
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!','h0_ginterp_apply',__FILE__)
!---Get dofs
allocate(j(self%h0_rep%nce),rop(3,self%h0_rep%nce))
call self%h0_rep%ncdofs(cell,j) ! get DOFs
!---Reconstruct field
CALL oft_h0_geval_all(self%h0_rep,cell,f,rop,gop)
val=0.d0
do jc=1,self%h0_rep%nce
  val=val+self%vals(j(jc))*rop(:,jc)
end do
DEALLOCATE(rop,j)
DEBUG_STACK_POP
end subroutine h0_ginterp_apply
!---------------------------------------------------------------------------
! SUBROUTINE: h0_zerob
!---------------------------------------------------------------------------
!> Zero a Nedelec H0 scalar field at all boundary nodes
!!
!! @param[in,out] a Field to be zeroed
!---------------------------------------------------------------------------
subroutine h0_zerob(a)
class(oft_vector), intent(inout) :: a
real(r8), pointer, dimension(:) :: aloc
integer(i4) :: i,j
DEBUG_STACK_PUSH
! Cast to vector type
NULLIFY(aloc)
CALL a%get_local(aloc)
! Apply operator
do i=1,oft_h0%nbe
  j=oft_h0%lbe(i)
  if(oft_h0%global%gbe(j))aloc(j)=0.d0
end do
CALL a%restore_local(aloc)
DEALLOCATE(aloc)
DEBUG_STACK_POP
end subroutine h0_zerob
!---------------------------------------------------------------------------
! SUBROUTINE: h0_zerogrnd
!---------------------------------------------------------------------------
!> Zero a Nedelec H0 scalar field at the global grounding node
!!
!! @note The possition of this node is defined by the mesh pointer igrnd in
!! mesh.
!!
!! @param[in,out] a Field to be zeroed
!---------------------------------------------------------------------------
subroutine h0_zerogrnd(a)
class(oft_vector), intent(inout) :: a
real(r8), pointer, dimension(:) :: aloc
integer(i4) :: i,j
DEBUG_STACK_PUSH
! Cast to vector type
NULLIFY(aloc)
CALL a%get_local(aloc)
!---
if(mesh%igrnd(1)>0)aloc(mesh%igrnd(1))=0.d0
if(mesh%igrnd(2)>0)aloc(mesh%igrnd(2))=0.d0
CALL a%restore_local(aloc)
DEALLOCATE(aloc)
DEBUG_STACK_POP
end subroutine h0_zerogrnd
!---------------------------------------------------------------------------
! SUBROUTINE: h0_zeroi
!---------------------------------------------------------------------------
!> Zero a Nedelec H0 scalar field at all interior nodes
!!
!! @param[in,out] a Field to be zeroed
!---------------------------------------------------------------------------
subroutine h0_zeroi(a)
class(oft_vector), intent(inout) :: a
real(r8), pointer, dimension(:) :: aloc
integer(i4) :: i
DEBUG_STACK_PUSH
!---Cast to vector type
NULLIFY(aloc)
CALL a%get_local(aloc)
! Apply operator
DO i=1,oft_h0%ne
  IF(oft_h0%global%gbe(i))CYCLE
  aloc(i)=0.d0
END DO
CALL a%restore_local(aloc)
DEALLOCATE(aloc)
DEBUG_STACK_POP
end subroutine h0_zeroi
!---------------------------------------------------------------------------
! SUBROUTINE: oft_h0_getmop
!---------------------------------------------------------------------------
!> Construct mass matrix for H0 scalar representation
!!
!! Supported boundary conditions
!! - \c 'none' Full matrix
!! - \c 'zerob' Dirichlet for all boundary DOF
!!
!! @param[in,out] mat Matrix object
!! @param[in] bc Boundary condition
!---------------------------------------------------------------------------
subroutine oft_h0_getmop(mat,bc)
class(oft_matrix), pointer, intent(inout) :: mat
character(LEN=*), intent(in) :: bc
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time
real(r8), allocatable :: rop(:),mop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_h0_vec
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing H0::MOP'
  CALL mytimer%tick()
END IF
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL oft_h0%mat_create(mat)
ELSE
  CALL mat%zero()
END IF
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,rop,det,mop,curved,goptmp,m,vol,jc,jr)
allocate(j(oft_h0%nce)) ! Local DOF and matrix indices
allocate(rop(oft_h0%nce)) ! Reconstructed gradient operator
allocate(mop(oft_h0%nce,oft_h0%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  mop=0.d0
  do m=1,oft_h0%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call mesh%jacobian(i,oft_h0%quad%pts(:,m),goptmp,vol)
    det=vol*oft_h0%quad%wts(m)
    CALL oft_h0_eval_all(oft_h0,i,oft_h0%quad%pts(:,m),rop)
    !---Compute local matrix contributions
    do jr=1,oft_h0%nce
      do jc=1,oft_h0%nce
        mop(jr,jc) = mop(jr,jc) + rop(jr)*rop(jc)*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call oft_h0%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,oft_h0%nce
        IF(oft_h0%global%gbe(j(jr)))mop(jr,:)=0.d0
      END DO
    CASE("grnd")
      IF(ANY(mesh%igrnd>0))THEN
        DO jr=1,oft_h0%nce
          IF(ANY(mesh%igrnd==j(jr)))mop(jr,:)=0.d0
        END DO
      END IF
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,mop,oft_h0%nce,oft_h0%nce)
  ! !$omp end critical
end do
deallocate(j,rop,mop)
!$omp end parallel
ALLOCATE(mop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
mop(1,1)=1.d0
SELECT CASE(TRIM(bc))
  CASE("zerob")
    DO i=1,oft_h0%nbe
      jr=oft_h0%lbe(i)
      IF(oft_h0%linkage%leo(i).AND.oft_h0%global%gbe(jr))THEN
        j=jr
        call mat%add_values(j,j,mop,1,1)
      END IF
    END DO
  CASE("grnd")
    IF(ANY(mesh%igrnd>0))THEN
      DO i=1,oft_h0%nbe
        jr=oft_h0%lbe(i)
        IF(oft_h0%linkage%leo(i).AND.ANY(jr==mesh%igrnd))THEN
          j=jr
          call mat%add_values(j,j,mop,1,1)
        END IF
      END DO
    END IF
END SELECT
DEALLOCATE(j,mop)
CALL oft_h0_create(oft_h0_vec)
CALL mat%assemble(oft_h0_vec)
CALL oft_h0_vec%delete
DEALLOCATE(oft_h0_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_h0_getmop
!---------------------------------------------------------------------------
! SUBROUTINE: oft_h0_getlop
!---------------------------------------------------------------------------
!> Construct laplacian matrix for H0 scalar representation
!!
!! Supported boundary conditions
!! - \c 'none' Full matrix
!! - \c 'zerob' Dirichlet for all boundary DOF
!! - \c 'grnd'  Dirichlet for only groundin point
!!
!! @param[in,out] mat Matrix object
!! @param[in] bc Boundary condition
!---------------------------------------------------------------------------
subroutine oft_h0_getlop(mat,bc)
class(oft_matrix), pointer, intent(inout) :: mat
character(LEN=*), intent(in) :: bc
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time
real(r8), allocatable :: gop(:,:),lop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: oft_h0_vec
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing H0::LOP'
  CALL mytimer%tick()
END IF
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL oft_h0%mat_create(mat)
ELSE
  CALL mat%zero()
END IF
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,gop,det,lop,curved,goptmp,m,vol,jc,jr)
allocate(j(oft_h0%nce)) ! Local DOF and matrix indices
allocate(gop(3,oft_h0%nce)) ! Reconstructed gradient operator
allocate(lop(oft_h0%nce,oft_h0%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  lop=0.d0
  do m=1,oft_h0%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call mesh%jacobian(i,oft_h0%quad%pts(:,m),goptmp,vol)
    det=vol*oft_h0%quad%wts(m)
    CALL oft_h0_geval_all(oft_h0,i,oft_h0%quad%pts(:,m),gop,goptmp)
    !---Compute local matrix contributions
    do jr=1,oft_h0%nce
      do jc=1,oft_h0%nce
        lop(jr,jc) = lop(jr,jc) + DOT_PRODUCT(gop(:,jr),gop(:,jc))*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call oft_h0%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,oft_h0%nce
        IF(oft_h0%global%gbe(j(jr)))lop(jr,:)=0.d0
      END DO
    CASE("grnd")
      IF(ANY(mesh%igrnd>0))THEN
        DO jr=1,oft_h0%nce
          IF(ANY(mesh%igrnd==j(jr)))lop(jr,:)=0.d0
        END DO
      END IF
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,lop,oft_h0%nce,oft_h0%nce)
  ! !$omp end critical
end do
deallocate(j,gop,lop)
!$omp end parallel
ALLOCATE(lop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
lop(1,1)=1.d0
SELECT CASE(TRIM(bc))
  CASE("zerob")
    DO i=1,oft_h0%nbe
      jr=oft_h0%lbe(i)
      IF(oft_h0%linkage%leo(i).AND.oft_h0%global%gbe(jr))THEN
        j=jr
        call mat%add_values(j,j,lop,1,1)
      END IF
    END DO
  CASE("grnd")
    IF(ANY(mesh%igrnd>0))THEN
      DO i=1,oft_h0%nbe
        jr=oft_h0%lbe(i)
        IF(oft_h0%linkage%leo(i).AND.ANY(jr==mesh%igrnd))THEN
          j=jr
          call mat%add_values(j,j,lop,1,1)
        END IF
      END DO
    END IF
END SELECT
DEALLOCATE(j,lop)
CALL oft_h0_create(oft_h0_vec)
CALL mat%assemble(oft_h0_vec)
CALL oft_h0_vec%delete
DEALLOCATE(oft_h0_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_h0_getlop
!---------------------------------------------------------------------------
! SUBROUTINE: oft_h0_project
!---------------------------------------------------------------------------
!> Project a scalar field onto a H0 basis
!!
!! @note This subroutine only performs the integration of the field with
!! test functions for a H0 basis. To retrieve the correct projection the
!! result must be multiplied by the inverse of H0::MOP.
!!
!! @param[in,out] field Vector field for projection
!! @param[in,out] x Field projected onto H0 basis
!---------------------------------------------------------------------------
subroutine oft_h0_project(field,x)
class(fem_interp), intent(inout) :: field
class(oft_vector), intent(inout) :: x
!---
real(r8) :: bcc(1),det,goptmp(3,4),vol
real(r8), pointer :: xloc(:)
real(r8), allocatable :: rop(:)
integer(i4) :: i,jc,m
integer(i4), allocatable :: j(:)
logical :: curved
DEBUG_STACK_PUSH
!---Initialize vectors to zero
NULLIFY(xloc)
call x%set(0.d0)
call x%get_local(xloc)
!---Integerate over the volume
!$omp parallel private(j,rop,curved,m,goptmp,vol,det,bcc,jc)
allocate(j(oft_h0%nce),rop(oft_h0%nce))
!$omp do schedule(guided)
do i=1,mesh%nc ! Loop over cells
  call oft_h0%ncdofs(i,j) ! Get DOFs
  curved=cell_is_curved(mesh,i) ! Straight cell test
  do m=1,oft_h0%quad%np
    if(curved.OR.m==1)call mesh%jacobian(i,oft_h0%quad%pts(:,m),goptmp,vol)
    det=vol*oft_h0%quad%wts(m)
    call field%interp(i,oft_h0%quad%pts(:,m),goptmp,bcc)
    CALL oft_h0_eval_all(oft_h0,i,oft_h0%quad%pts(:,m),rop)
    do jc=1,oft_h0%nce
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
end subroutine oft_h0_project
!---------------------------------------------------------------------------
! SUBROUTINE: h0_setup_interp
!---------------------------------------------------------------------------
!> Construct interpolation matrices for transfer between H0 finite element
!! spaces
!---------------------------------------------------------------------------
SUBROUTINE h0_setup_interp
INTEGER(i4) :: i
DEBUG_STACK_PUSH
!---
DO i=oft_h0_minlev+1,oft_h0_nlevels
  CALL oft_h0_set_level(i)
  !---
  if(oft_h0_level==oft_h0_blevel+1)CYCLE
  !---Setup interpolation
  if(oft_h0%order==1)then
    CALL h0_ginterpmatrix(ML_oft_h0%interp_matrices(ML_oft_h0%level)%m)
    oft_h0_ops%interp=>ML_oft_h0%interp_matrices(ML_oft_h0%level)%m
    CALL oft_h0_ops%interp%assemble
  else
    CALL h0_pinterpmatrix(ML_oft_h0%interp_matrices(ML_oft_h0%level)%m)
    oft_h0_ops%interp=>ML_oft_h0%interp_matrices(ML_oft_h0%level)%m
    CALL oft_h0_ops%interp%assemble
  end if
END DO
DEBUG_STACK_POP
END SUBROUTINE h0_setup_interp
!---------------------------------------------------------------------------
! SUBROUTINE: h0_ginterpmatrix
!---------------------------------------------------------------------------
!> Construct interpolation matrix for transfer between geometric levels
!! of H0 finite element space
!---------------------------------------------------------------------------
SUBROUTINE h0_ginterpmatrix(mat)
class(oft_matrix), pointer, intent(inout) :: mat
INTEGER(i4) :: i,j,k,m,icors,ifine,jb,i_ind(1),j_ind(1)
INTEGER(i4) :: etmp(2),ftmp(3),fetmp(3),ctmp(4),fc,ed
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pmap,emap,fmap
CLASS(oft_afem_type), POINTER :: h0_cors => NULL()
TYPE(h0_ops), POINTER :: ops
class(oft_mesh), pointer :: cmesh
CLASS(oft_vector), POINTER :: h0_vec,h0_vec_cors
integer(i4) :: jfe(3),jce(6)
integer(i4), pointer :: lfde(:,:),lede(:,:)
real(r8) :: f(4),incr,val,d(3),h_rop(3),goptmp(3,4),v,mop(1)
type(oft_graph_ptr), pointer :: graphs(:,:)
DEBUG_STACK_PUSH
!---
if(mg_mesh%level<1)then
  call oft_abort('Invalid mesh level','h0_ginterpmatrix',__FILE__)
end if
cmesh=>mg_mesh%meshes(mg_mesh%level-1)
if(cmesh%type/=1)CALL oft_abort("Only supported with tet meshes", &
  "h0_ginterpmatrix", __FILE__)
if(oft_h0%order/=1)then
  call oft_abort('Attempted geometric interpolation for pd > 1','h0_ginterpmatrix',__FILE__)
end if
ops=>oft_h0_ops
h0_cors=>ML_oft_h0%levels(oft_h0_level-1)%fe
lede=>mg_mesh%inter(mg_mesh%level-1)%lede
lfde=>mg_mesh%inter(mg_mesh%level-1)%lfde
ALLOCATE(ML_oft_h0%interp_graphs(ML_oft_h0%level)%g)
ops%interp_graph=>ML_oft_h0%interp_graphs(ML_oft_h0%level)%g
!---Setup matrix sizes
ops%interp_graph%nr=oft_h0%ne
ops%interp_graph%nrg=oft_h0%global%ne
ops%interp_graph%nc=h0_cors%ne
ops%interp_graph%ncg=h0_cors%global%ne
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
if(ops%interp_graph%kr(1)/=1)call oft_abort('Bad element to element count','oft_h0_ginterpmatrix',__FILE__)
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
CALL oft_h0_create(h0_vec)
CALL oft_h0_create(h0_vec_cors,oft_h0_level-1)
!---
ALLOCATE(graphs(1,1))
graphs(1,1)%g=>ops%interp_graph
!---
CALL create_matrix(mat,graphs,h0_vec,h0_vec_cors)
CALL h0_vec%delete
CALL h0_vec_cors%delete
DEALLOCATE(graphs,h0_vec,h0_vec_cors)
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
END SUBROUTINE h0_ginterpmatrix
!---------------------------------------------------------------------------
! SUBROUTINE: h0_pinterpmatrix
!---------------------------------------------------------------------------
!> Construct interpolation matrix for transfer between polynomial levels
!! of H0 finite element space
!---------------------------------------------------------------------------
SUBROUTINE h0_pinterpmatrix(mat)
class(oft_matrix), pointer, intent(inout) :: mat
INTEGER(i4) :: i,j,k,m,icors,ifine,jb,js,jn,i_ind(1),j_ind(1)
INTEGER(i4) :: etmp(2),ftmp(3),fetmp(3),ctmp(4),fc,ed
INTEGER(i4) :: offsetc,offsetf
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pmap,emap,fmap
REAL(r8) :: f(4),val,mop(1)
CLASS(oft_fem_type), POINTER :: h0_cors => NULL()
TYPE(h0_ops), POINTER :: ops
CLASS(oft_vector), POINTER :: h0_vec,h0_vec_cors
type(oft_graph_ptr), pointer :: graphs(:,:)
DEBUG_STACK_PUSH
!---
ops=>oft_h0_ops
SELECT TYPE(this=>ML_oft_h0%levels(oft_h0_level-1)%fe)
CLASS IS(oft_fem_type)
  h0_cors=>this
CLASS DEFAULT
  CALL oft_abort("Error getting coarse FE object","h0_pinterpmatrix",__FILE__)
END SELECT
ALLOCATE(ML_oft_h0%interp_graphs(ML_oft_h0%level)%g)
ops%interp_graph=>ML_oft_h0%interp_graphs(ML_oft_h0%level)%g
!---Setup matrix sizes
ops%interp_graph%nr=oft_h0%ne
ops%interp_graph%nrg=oft_h0%global%ne
ops%interp_graph%nc=h0_cors%ne
ops%interp_graph%ncg=h0_cors%global%ne
!---Setup Matrix graph
ALLOCATE(ops%interp_graph%kr(ops%interp_graph%nr+1))
ops%interp_graph%kr=0
ops%interp_graph%nnz=h0_cors%ne
ALLOCATE(ops%interp_graph%lc(ops%interp_graph%nnz))
ops%interp_graph%lc=0_i4
!---Construct matrix
do i=1,mesh%np
  ops%interp_graph%kr(i)=1
end do
!---
do i=1,mesh%ne
  do j=1,h0_cors%gstruct(2)
    offsetf=mesh%np+(i-1)*oft_h0%gstruct(2)
    ops%interp_graph%kr(j+offsetf)=1
  end do
end do
!---
do i=1,mesh%nf
  do j=1,h0_cors%gstruct(3)
    offsetf=mesh%np+oft_h0%gstruct(2)*mesh%ne+(i-1)*oft_h0%gstruct(3)
    ops%interp_graph%kr(j+offsetf)=1
  end do
end do
!---
do i=1,mesh%nc
  do j=1,h0_cors%gstruct(4)
    offsetf=mesh%np+oft_h0%gstruct(2)*mesh%ne+oft_h0%gstruct(3)*mesh%nf+(i-1)*oft_h0%gstruct(4)
    ops%interp_graph%kr(j+offsetf)=1
  end do
end do
ops%interp_graph%kr(ops%interp_graph%nr+1)=ops%interp_graph%nnz+1
do i=ops%interp_graph%nr,1,-1 ! cumulative point to point count
  ops%interp_graph%kr(i)=ops%interp_graph%kr(i+1)-ops%interp_graph%kr(i)
end do
if(ops%interp_graph%kr(1)/=1)call oft_abort('Bad element to element count','oft_h0_pinterpmatrix',__FILE__)
!---Construct matrix
do i=1,mesh%np
  ops%interp_graph%lc(ops%interp_graph%kr(i))=i
end do
!---
do i=1,mesh%ne
  do j=1,h0_cors%gstruct(2)
    offsetf=mesh%np+(i-1)*oft_h0%gstruct(2)
    offsetc=mesh%np+(i-1)*h0_cors%gstruct(2)
    ops%interp_graph%lc(ops%interp_graph%kr(j+offsetf))=j+offsetc
  end do
end do
!---
do i=1,mesh%nf
  do j=1,h0_cors%gstruct(3)
    offsetf=mesh%np+oft_h0%gstruct(2)*mesh%ne+(i-1)*oft_h0%gstruct(3)
    offsetc=mesh%np+h0_cors%gstruct(2)*mesh%ne+(i-1)*h0_cors%gstruct(3)
    ops%interp_graph%lc(ops%interp_graph%kr(j+offsetf))=j+offsetc
  end do
end do
!---
do i=1,mesh%nc
  do j=1,h0_cors%gstruct(4)
    offsetf=mesh%np+oft_h0%gstruct(2)*mesh%ne+oft_h0%gstruct(3)*mesh%nf+(i-1)*oft_h0%gstruct(4)
    offsetc=mesh%np+h0_cors%gstruct(2)*mesh%ne+h0_cors%gstruct(3)*mesh%nf+(i-1)*h0_cors%gstruct(4)
    ops%interp_graph%lc(ops%interp_graph%kr(j+offsetf))=j+offsetc
  end do
end do
!---------------------------------------------------------------------------
! Construct matrix
!---------------------------------------------------------------------------
CALL oft_h0_create(h0_vec)
CALL oft_h0_create(h0_vec_cors,oft_h0_level-1)
!---
ALLOCATE(graphs(1,1))
graphs(1,1)%g=>ops%interp_graph
!---
CALL create_matrix(mat,graphs,h0_vec,h0_vec_cors)
CALL h0_vec%delete
CALL h0_vec_cors%delete
DEALLOCATE(graphs,h0_vec,h0_vec_cors)
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!---
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
  DO j=1,h0_cors%gstruct(2)
    i_ind=j+mesh%np+(i-1)*oft_h0%gstruct(2)
    j_ind=j+mesh%np+(i-1)*h0_cors%gstruct(2)
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
  DO j=1,h0_cors%gstruct(3)
    i_ind=j+mesh%np+oft_h0%gstruct(2)*mesh%ne+(i-1)*oft_h0%gstruct(3)
    j_ind=j+mesh%np+h0_cors%gstruct(2)*mesh%ne+(i-1)*h0_cors%gstruct(3)
    mop=1.d0
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
  END DO
END DO
deallocate(fmap)
!---
DO i=1,mesh%nc
  DO j=1,h0_cors%gstruct(4)
    i_ind=j+mesh%np+oft_h0%gstruct(2)*mesh%ne+oft_h0%gstruct(3)*mesh%nf+(i-1)*oft_h0%gstruct(4)
    j_ind=j+mesh%np+h0_cors%gstruct(2)*mesh%ne+h0_cors%gstruct(3)*mesh%nf+(i-1)*h0_cors%gstruct(4)
    mop=1.d0
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
  END DO
END DO
DEBUG_STACK_POP
END SUBROUTINE h0_pinterpmatrix
!---------------------------------------------------------------------------
! SUBROUTINE: h0_interp
!---------------------------------------------------------------------------
!> Interpolate a coarse level H0 scalar field to the next finest level
!!
!! @note The global H0 level in incremented by one in this subroutine
!!
!! @param[in] acors Vector to interpolate
!! @param[in,out] afine Fine vector from interpolation
!---------------------------------------------------------------------------
subroutine h0_interp(acors,afine)
class(oft_vector), intent(inout) :: acors
class(oft_vector), intent(inout) :: afine
DEBUG_STACK_PUSH
!---Step one level up
call oft_h0_set_level(oft_h0_level+1)
call afine%set(0.d0)
!---
if(oft_h0_level==oft_h0_blevel+1)then
  call h0_base_pop(acors,afine)
  DEBUG_STACK_POP
  return
end if
CALL oft_h0_ops%interp%apply(acors,afine)
DEBUG_STACK_POP
end subroutine h0_interp
!---------------------------------------------------------------------------
! SUBROUTINE: h0_base_pop
!---------------------------------------------------------------------------
!> Transfer a base level H0 scalar field to the next MPI level
!!
!! @param[in] acors Vector to transfer
!! @param[in,out] afine Fine vector from transfer
!---------------------------------------------------------------------------
subroutine h0_base_pop(acors,afine)
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
end subroutine h0_base_pop
!---------------------------------------------------------------------------
! SUBROUTINE: h0_inject
!---------------------------------------------------------------------------
!> Inject a fine level H0 scalar field to the next coarsest level
!!
!! @note The global H0 level in decremented by one in this subroutine
!!
!! @param[in] afine Vector to inject
!! @param[in,out] acors Coarse vector from injection
!---------------------------------------------------------------------------
subroutine h0_inject(afine,acors)
class(oft_vector), intent(inout) :: afine
class(oft_vector), intent(inout) :: acors
integer(i4) :: i,j,k
logical :: gcheck
DEBUG_STACK_PUSH
gcheck=(oft_h0%order==1)
! Step down level up
call oft_h0_set_level(oft_h0_level-1)
! Cast fine field
call acors%set(0.d0)
if(oft_h0_level==oft_h0_blevel)then
  call h0_base_push(afine,acors)
  DEBUG_STACK_POP
  return
end if
CALL ML_oft_h0_ops(oft_h0_level+1)%interp%applyT(afine,acors)
DEBUG_STACK_POP
end subroutine h0_inject
!---------------------------------------------------------------------------
! SUBROUTINE: h0_base_push
!---------------------------------------------------------------------------
!> Transfer a MPI level H0 scalar field to the base level
!!
!! @param[in] afine Vector to transfer
!! @param[in,out] acors Fine vector from transfer
!---------------------------------------------------------------------------
subroutine h0_base_push(afine,acors)
class(oft_vector), intent(inout) :: afine
class(oft_vector), intent(inout) :: acors
integer(i4), pointer :: lptmp(:)
integer(i4) :: i,j,ierr
real(r8), pointer, dimension(:) :: alias,array_c,array_f
CLASS(oft_afem_type), POINTER :: h0_fine => NULL()
DEBUG_STACK_PUSH
!---
lptmp=>mg_mesh%meshes(mg_mesh%nbase+1)%base%lp
CALL acors%get_local(array_c)
CALL afine%get_local(array_f)
h0_fine=>ML_oft_h0%levels(oft_h0_level+1)%fe
!---
allocate(alias(acors%n))
alias=0.d0
!$omp parallel do
do i=1,afine%n
  if(h0_fine%linkage%be(i))cycle
  alias(lptmp(i))=array_f(i)
end do
!$omp parallel do private(j)
do i=1,h0_fine%linkage%nbe
  j=h0_fine%linkage%lbe(i)
  if(.NOT.h0_fine%linkage%leo(i))cycle
  alias(lptmp(j))=array_f(j)
end do
!---Global reduction over all processors
array_c=oft_mpi_sum(alias,acors%n)
call acors%restore_local(array_c)
deallocate(alias,array_c,array_f)
DEBUG_STACK_POP
end subroutine h0_base_push
!---------------------------------------------------------------------------
! SUBROUTINE: h0_lop_eigs
!---------------------------------------------------------------------------
!> Compute eigenvalues and smoothing coefficients for the operator H0::LOP
!---------------------------------------------------------------------------
SUBROUTINE h0_lop_eigs(minlev)
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
IF(oft_env%head_proc)WRITE(*,*)'Optimizing Jacobi damping for H0::LOP'
ALLOCATE(df(oft_h0_nlevels))
df=0.d0
DO i=minlev,oft_h0_nlevels
  CALL oft_h0_set_level(i)
  !---Create fields
  CALL oft_h0_create(u)
  !---Get Ev range
  NULLIFY(lop)
  CALL oft_h0_getlop(lop,'grnd')
  CALL create_diagmatrix(md,lop%D)
  !---
  arsolver%A=>lop
  arsolver%M=>md
  arsolver%mode=2
  arsolver%tol=1.E-5_r8
  arsolver%bc=>h0_zerogrnd
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
  DO i=1,oft_h0_nlevels-1
    WRITE(*,'(F5.3,A)',ADVANCE='NO')df(i),', '
  END DO
  WRITE(*,'(F5.3,A)')df(oft_h0_nlevels)
END IF
DEALLOCATE(df)
DEBUG_STACK_POP
#else
CALL oft_abort("Subroutine requires ARPACK", "lag_lop_eigs", __FILE__)
#endif
END SUBROUTINE h0_lop_eigs
!---------------------------------------------------------------------------
! SUBROUTINE: h0_getlop_pre
!---------------------------------------------------------------------------
!> Compute eigenvalues and smoothing coefficients for the operator H0::LOP
!---------------------------------------------------------------------------
SUBROUTINE h0_getlop_pre(pre,mats,level,nlevels)
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
TYPE(xml_node), POINTER :: h0_node
#endif
DEBUG_STACK_PUSH
!---
minlev=1
toplev=oft_h0_level
levin=oft_h0_level
IF(PRESENT(level))toplev=level
IF(PRESENT(nlevels))minlev=toplev-nlevels+1
nl=toplev-minlev+1
!---
IF(minlev<1)CALL oft_abort('Minimum level is < 0','h0_getlop_pre',__FILE__)
IF(toplev>oft_h0_nlevels)CALL oft_abort('Maximum level is > h0_nlevels','h0_getlop_pre',__FILE__)
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
  CALL oft_h0_set_level(minlev+(i-1))
  levels(i)=minlev+(i-1)
  df(i)=df_lop(levels(i))
  nu(i)=nu_lop(levels(i))
  IF(df(i)<-1.d90)THEN
    WRITE(lev_char,'(I2.2)')levels(i)
    CALL oft_abort('Smoother values not set for level: '//lev_char,'h0_getlop_pre',__FILE__)
  END IF
  !---
  IF(create_mats)THEN
    NULLIFY(mats(i)%M)
    CALL oft_h0_getlop(mats(i)%M,'grnd')
  END IF
  IF(i>1)ml_int(i-1)%M=>oft_h0_ops%interp
END DO
CALL oft_h0_set_level(levin)
!---------------------------------------------------------------------------
! Search for XML-spec
!---------------------------------------------------------------------------
NULLIFY(pre_node)
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml))THEN
  !---Look for Lagrange node
  CALL xml_get_element(oft_env%xml,"nedelec_h0",h0_node,ierr)
  IF(ierr==0)CALL xml_get_element(h0_node,"lop",pre_node,ierr)
END IF
#endif
!---------------------------------------------------------------------------
! Setup preconditioner
!---------------------------------------------------------------------------
NULLIFY(pre)
CALL create_mlpre(pre,mats(1:nl),levels,nlevels=nl,create_vec=oft_h0_create,interp=h0_interp, &
     inject=h0_inject,bc=h0_zerogrnd,stype=1,df=df,nu=nu,xml_root=pre_node)
!---------------------------------------------------------------------------
! Cleanup
!---------------------------------------------------------------------------
DEALLOCATE(ml_int,levels,df,nu)
DEBUG_STACK_POP
END SUBROUTINE h0_getlop_pre
end module oft_h0_operators
