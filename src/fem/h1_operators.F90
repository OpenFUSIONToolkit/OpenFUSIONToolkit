!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!> @file oft_h1_operators.F90
!
!> H^1 FE operator definitions
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
!! @ingroup doxy_oft_h1
!------------------------------------------------------------------------------
MODULE oft_h1_operators
USE oft_base
USE oft_mesh_type, ONLY: oft_mesh, cell_is_curved
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr, oft_graph_ptr, &
  oft_graph
USE oft_deriv_matrices, ONLY: oft_diagmatrix, create_diagmatrix
USE oft_solver_base, ONLY: oft_solver, oft_solver_bc
USE oft_la_utils, ONLY: create_matrix
USE oft_solver_utils, ONLY: create_mlpre, create_native_pre
#ifdef HAVE_ARPACK
USE oft_arpack, ONLY: oft_irlm_eigsolver
#endif
USE fem_base, ONLY: oft_afem_type, oft_fem_type, fem_max_levels, oft_ml_fem_type, &
  oft_ml_fe_vecspace
USE fem_utils, ONLY: fem_interp
USE oft_h1_basis, ONLY: oft_h1_eval_all, oft_h1_geval_all, &
  oft_h1_fem, oft_3D_h1_cast
IMPLICIT NONE
#include "local.h"
!------------------------------------------------------------------------------
!> Interpolate a H^1 field
!------------------------------------------------------------------------------
type, extends(fem_interp) :: oft_h1_rinterp
  class(oft_vector), pointer :: u => NULL() !< Field to interpolate
  real(r8), pointer, dimension(:) :: vals => NULL() !< Local values
  class(oft_h1_fem), pointer :: h1_rep => NULL() !< FE representation
contains
  !> Retrieve local values for interpolation
  procedure :: setup => rinterp_setup
  !> Reconstruct field
  procedure :: interp => rinterp_apply
  !> Destroy temporary internal storage
  procedure :: delete => rinterp_delete
end type oft_h1_rinterp
!------------------------------------------------------------------------------
!> Interpolate \f$ \nabla \f$ of a H^1 field
!------------------------------------------------------------------------------
type, extends(oft_h1_rinterp) :: oft_h1_ginterp
contains
  !> Reconstruct field
  procedure :: interp => ginterp_apply
end type oft_h1_ginterp
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(oft_solver_bc) :: oft_h1_zerob
  class(oft_ml_fem_type), pointer :: ML_H1_rep => NULL() !< FE representation
contains
  procedure :: apply => zerob_apply
  procedure :: delete => zerob_delete
end type oft_h1_zerob
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(oft_h1_zerob) :: oft_h1_zerogrnd
contains
  procedure :: apply => zerogrnd_apply
end type oft_h1_zerogrnd
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(oft_h1_zerob) :: oft_h1_zeroi
contains
  procedure :: apply => zeroi_apply
end type oft_h1_zeroi
!---Pre options
integer(i4) :: nu_lop(fem_max_levels)=0
real(r8) :: df_lop(fem_max_levels)=-1.d99
!---Cache variables
REAL(r8), POINTER, DIMENSION(:) :: oft_h1_rop => NULL()
REAL(r8), POINTER, DIMENSION(:,:) :: oft_h1_gop => NULL()
!$omp threadprivate(oft_h1_rop,oft_h1_gop)
contains
!------------------------------------------------------------------------------
!> Read-in options for the basic H^1 ML preconditioners
!------------------------------------------------------------------------------
subroutine h1_mloptions
integer(i4) :: ierr,io_unit
namelist/h1_op_options/df_lop,nu_lop
DEBUG_STACK_PUSH
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,h1_op_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr>0)THEN
  CALL oft_abort('Error reading "h1_op_options" group in input file','h1_mloptions',__FILE__)
END IF
IF(df_lop(1)<-1.d90)THEN
  IF(oft_env%head_proc)THEN
    CALL oft_warn("No H^1 MG smoother settings found:")
    CALL oft_warn("  Using default values, which may result in convergence failure.")
  END IF
  nu_lop=2
  df_lop=.2d0
END IF
DEBUG_STACK_POP
end subroutine h1_mloptions
!------------------------------------------------------------------------------
!> Setup interpolator for H^1 scalar fields
!!
!! Fetches local representation used for interpolation from vector object
!------------------------------------------------------------------------------
subroutine rinterp_setup(self,h1_rep)
class(oft_h1_rinterp), intent(inout) :: self
class(oft_afem_type), target, intent(inout) :: h1_rep
IF(.NOT.oft_3D_h1_cast(self%h1_rep,h1_rep))CALL oft_abort("Incorrect FE type","h1_rinterp_setup",__FILE__)
self%mesh=>self%h1_rep%mesh
!---Get local slice
CALL self%u%get_local(self%vals)
end subroutine rinterp_setup
!------------------------------------------------------------------------------
!> Destroy temporary internal storage
!------------------------------------------------------------------------------
subroutine rinterp_delete(self)
class(oft_h1_rinterp), intent(inout) :: self
!---Deallocate local storage
IF(ASSOCIATED(self%vals))DEALLOCATE(self%vals)
NULLIFY(self%h1_rep,self%u)
end subroutine rinterp_delete
!------------------------------------------------------------------------------
!> Reconstruct a H^1 scalar field
!------------------------------------------------------------------------------
subroutine rinterp_apply(self,cell,f,gop,val)
class(oft_h1_rinterp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [4]
real(r8), intent(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
real(r8), intent(out) :: val(:) !< Reconstructed field at f [1]
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8), allocatable :: rop(:)
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!','rinterp_apply',__FILE__)
!---Get dofs
allocate(j(self%h1_rep%nce),rop(self%h1_rep%nce))
call self%h1_rep%ncdofs(cell,j) ! get DOFs
!---Reconstruct field
CALL oft_h1_eval_all(self%h1_rep,cell,f,rop)
val=0.d0
do jc=1,self%h1_rep%nce
  val=val+self%vals(j(jc))*rop(jc)
end do
deallocate(rop,j)
DEBUG_STACK_POP
end subroutine rinterp_apply
!------------------------------------------------------------------------------
!> Reconstruct the gradient of a H^1 scalar field
!------------------------------------------------------------------------------
subroutine ginterp_apply(self,cell,f,gop,val)
class(oft_h1_ginterp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [4]
real(r8), intent(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
real(r8), intent(out) :: val(:) !< Reconstructed field at f [3]
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8), allocatable :: rop(:,:)
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%vals))CALL oft_abort('Setup has not been called!','ginterp_apply',__FILE__)
!---Get dofs
allocate(j(self%h1_rep%nce),rop(3,self%h1_rep%nce))
call self%h1_rep%ncdofs(cell,j) ! get DOFs
!---Reconstruct field
CALL oft_h1_geval_all(self%h1_rep,cell,f,rop,gop)
val=0.d0
do jc=1,self%h1_rep%nce
  val=val+self%vals(j(jc))*rop(:,jc)
end do
DEALLOCATE(rop,j)
DEBUG_STACK_POP
end subroutine ginterp_apply
!------------------------------------------------------------------------------
!> Zero a H^1 scalar field at all boundary nodes
!------------------------------------------------------------------------------
subroutine zerob_apply(self,a)
class(oft_h1_zerob), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to be zeroed
real(r8), pointer, dimension(:) :: aloc
integer(i4) :: i,j
class(oft_afem_type), pointer :: fe_rep
DEBUG_STACK_PUSH
! Cast to vector type
NULLIFY(aloc)
CALL a%get_local(aloc)
! Apply operator
fe_rep=>self%ML_H1_rep%current_level
do i=1,fe_rep%nbe
  j=fe_rep%lbe(i)
  if(fe_rep%global%gbe(j))aloc(j)=0.d0
end do
CALL a%restore_local(aloc)
DEALLOCATE(aloc)
DEBUG_STACK_POP
end subroutine zerob_apply
!------------------------------------------------------------------------------
!> Zero a H^1 scalar field at all boundary nodes
!------------------------------------------------------------------------------
subroutine zerob_delete(self)
class(oft_h1_zerob), intent(inout) :: self
NULLIFY(self%ML_H1_rep)
end subroutine zerob_delete
!------------------------------------------------------------------------------
!> Zero a H^1 scalar field at the global grounding node
!!
!! @note The possition of this node is defined by the mesh pointer igrnd in
!! mesh
!------------------------------------------------------------------------------
subroutine zerogrnd_apply(self,a)
class(oft_h1_zerogrnd), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to be zeroed
real(r8), pointer, dimension(:) :: aloc
integer(i4) :: i,j
CLASS(oft_h1_fem), POINTER :: this
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_h1_cast(this,self%ML_H1_rep%current_level))CALL oft_abort("Invalid fe object","zerogrnd_apply",__FILE__)
! Cast to vector type
NULLIFY(aloc)
CALL a%get_local(aloc)
!---
if(this%mesh%igrnd(1)>0)aloc(this%mesh%igrnd(1))=0.d0
if(this%mesh%igrnd(2)>0)aloc(this%mesh%igrnd(2))=0.d0
CALL a%restore_local(aloc)
DEALLOCATE(aloc)
DEBUG_STACK_POP
end subroutine zerogrnd_apply
!------------------------------------------------------------------------------
!> Zero a H^1 scalar field at all interior nodes
!------------------------------------------------------------------------------
subroutine zeroi_apply(self,a)
class(oft_h1_zeroi), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Field to be zeroed
real(r8), pointer, dimension(:) :: aloc
integer(i4) :: i
DEBUG_STACK_PUSH
!---Cast to vector type
NULLIFY(aloc)
CALL a%get_local(aloc)
! Apply operator
DO i=1,self%ML_H1_rep%current_level%ne
  IF(self%ML_H1_rep%current_level%global%gbe(i))CYCLE
  aloc(i)=0.d0
END DO
CALL a%restore_local(aloc)
DEALLOCATE(aloc)
DEBUG_STACK_POP
end subroutine zeroi_apply
!------------------------------------------------------------------------------
!> Construct mass matrix for H^1 scalar representation
!!
!! Supported boundary conditions
!! - `'none'` Full matrix
!! - `'zerob'` Dirichlet for all boundary DOF
!------------------------------------------------------------------------------
subroutine oft_h1_getmop(fe_rep,mat,bc)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix object
character(LEN=*), intent(in) :: bc !< Boundary condition
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time
real(r8), allocatable :: rop(:),mop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: h1_vec
CLASS(oft_h1_fem), POINTER :: h1_rep
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing H^1::MOP'
  CALL mytimer%tick()
END IF
IF(.NOT.oft_3D_h1_cast(h1_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_h1_getmop",__FILE__)
!------------------------------------------------------------------------------
! Allocate matrix
!------------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL h1_rep%mat_create(mat)
ELSE
  CALL mat%zero()
END IF
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,rop,det,mop,curved,goptmp,m,vol,jc,jr)
allocate(j(h1_rep%nce)) ! Local DOF and matrix indices
allocate(rop(h1_rep%nce)) ! Reconstructed gradient operator
allocate(mop(h1_rep%nce,h1_rep%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,h1_rep%mesh%nc
  curved=cell_is_curved(h1_rep%mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  mop=0.d0
  do m=1,h1_rep%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call h1_rep%mesh%jacobian(i,h1_rep%quad%pts(:,m),goptmp,vol)
    det=vol*h1_rep%quad%wts(m)
    CALL oft_h1_eval_all(h1_rep,i,h1_rep%quad%pts(:,m),rop)
    !---Compute local matrix contributions
    do jr=1,h1_rep%nce
      do jc=1,h1_rep%nce
        mop(jr,jc) = mop(jr,jc) + rop(jr)*rop(jc)*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call h1_rep%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,h1_rep%nce
        IF(h1_rep%global%gbe(j(jr)))mop(jr,:)=0.d0
      END DO
    CASE("grnd")
      IF(ANY(h1_rep%mesh%igrnd>0))THEN
        DO jr=1,h1_rep%nce
          IF(ANY(h1_rep%mesh%igrnd==j(jr)))mop(jr,:)=0.d0
        END DO
      END IF
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,mop,h1_rep%nce,h1_rep%nce)
  ! !$omp end critical
end do
deallocate(j,rop,mop)
!$omp end parallel
ALLOCATE(mop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
mop(1,1)=1.d0
SELECT CASE(TRIM(bc))
  CASE("zerob")
    DO i=1,h1_rep%nbe
      jr=h1_rep%lbe(i)
      IF(h1_rep%linkage%leo(i).AND.h1_rep%global%gbe(jr))THEN
        j=jr
        call mat%add_values(j,j,mop,1,1)
      END IF
    END DO
  CASE("grnd")
    IF(ANY(h1_rep%mesh%igrnd>0))THEN
      DO i=1,h1_rep%nbe
        jr=h1_rep%lbe(i)
        IF(h1_rep%linkage%leo(i).AND.ANY(jr==h1_rep%mesh%igrnd))THEN
          j=jr
          call mat%add_values(j,j,mop,1,1)
        END IF
      END DO
    END IF
END SELECT
DEALLOCATE(j,mop)
CALL h1_rep%vec_create(h1_vec)
CALL mat%assemble(h1_vec)
CALL h1_vec%delete
DEALLOCATE(h1_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_h1_getmop
!------------------------------------------------------------------------------
!> Construct laplacian matrix for H^1 scalar representation
!!
!! Supported boundary conditions
!! - `'none'` Full matrix
!! - `'zerob'` Dirichlet for all boundary DOF
!! - `'grnd'`  Dirichlet for only groundin point
!------------------------------------------------------------------------------
subroutine oft_h1_getlop(fe_rep,mat,bc)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(oft_matrix), pointer, intent(inout) :: mat !< Matrix object
character(LEN=*), intent(in) :: bc !< Boundary condition
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,4),elapsed_time
real(r8), allocatable :: gop(:,:),lop(:,:)
logical :: curved
CLASS(oft_vector), POINTER :: h1_vec
CLASS(oft_h1_fem), POINTER :: h1_rep
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2X,A)')'Constructing H^1::LOP'
  CALL mytimer%tick()
END IF
IF(.NOT.oft_3D_h1_cast(h1_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_h1_getlop",__FILE__)
!------------------------------------------------------------------------------
! Allocate matrix
!------------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
  CALL h1_rep%mat_create(mat)
ELSE
  CALL mat%zero()
END IF
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!---Operator integration loop
!$omp parallel private(j,gop,det,lop,curved,goptmp,m,vol,jc,jr)
allocate(j(h1_rep%nce)) ! Local DOF and matrix indices
allocate(gop(3,h1_rep%nce)) ! Reconstructed gradient operator
allocate(lop(h1_rep%nce,h1_rep%nce)) ! Local laplacian matrix
!$omp do schedule(guided)
do i=1,h1_rep%mesh%nc
  curved=cell_is_curved(h1_rep%mesh,i) ! Straight cell test
  !---Get local reconstructed operators
  lop=0.d0
  do m=1,h1_rep%quad%np ! Loop over quadrature points
    if(curved.OR.m==1)call h1_rep%mesh%jacobian(i,h1_rep%quad%pts(:,m),goptmp,vol)
    det=vol*h1_rep%quad%wts(m)
    CALL oft_h1_geval_all(h1_rep,i,h1_rep%quad%pts(:,m),gop,goptmp)
    !---Compute local matrix contributions
    do jr=1,h1_rep%nce
      do jc=1,h1_rep%nce
        lop(jr,jc) = lop(jr,jc) + DOT_PRODUCT(gop(:,jr),gop(:,jc))*det
      end do
    end do
  end do
  !---Get local to global DOF mapping
  call h1_rep%ncdofs(i,j)
  !---Apply bc to local matrix
  SELECT CASE(TRIM(bc))
    CASE("zerob")
      DO jr=1,h1_rep%nce
        IF(h1_rep%global%gbe(j(jr)))lop(jr,:)=0.d0
      END DO
    CASE("grnd")
      IF(ANY(h1_rep%mesh%igrnd>0))THEN
        DO jr=1,h1_rep%nce
          IF(ANY(h1_rep%mesh%igrnd==j(jr)))lop(jr,:)=0.d0
        END DO
      END IF
  END SELECT
  !---Add local values to global matrix
  ! !$omp critical
  call mat%atomic_add_values(j,j,lop,h1_rep%nce,h1_rep%nce)
  ! !$omp end critical
end do
deallocate(j,gop,lop)
!$omp end parallel
ALLOCATE(lop(1,1),j(1))
!---Set diagonal entries for dirichlet rows
lop(1,1)=1.d0
SELECT CASE(TRIM(bc))
  CASE("zerob")
    DO i=1,h1_rep%nbe
      jr=h1_rep%lbe(i)
      IF(h1_rep%linkage%leo(i).AND.h1_rep%global%gbe(jr))THEN
        j=jr
        call mat%add_values(j,j,lop,1,1)
      END IF
    END DO
  CASE("grnd")
    IF(ANY(h1_rep%mesh%igrnd>0))THEN
      DO i=1,h1_rep%nbe
        jr=h1_rep%lbe(i)
        IF(h1_rep%linkage%leo(i).AND.ANY(jr==h1_rep%mesh%igrnd))THEN
          j=jr
          call mat%add_values(j,j,lop,1,1)
        END IF
      END DO
    END IF
END SELECT
DEALLOCATE(j,lop)
CALL h1_rep%vec_create(h1_vec)
CALL mat%assemble(h1_vec)
CALL h1_vec%delete
DEALLOCATE(h1_vec)
IF(oft_debug_print(1))THEN
  elapsed_time=mytimer%tock()
  WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine oft_h1_getlop
!------------------------------------------------------------------------------
!> Project a scalar field onto a H^1 basis
!!
!! @note This subroutine only performs the integration of the field with
!! test functions for a H^1 basis. To retrieve the correct projection the
!! result must be multiplied by the inverse of H^1::MOP.
!!
!! @param[in,out] field Vector field for projection
!! @param[in,out] x Field projected onto H^1 basis
!------------------------------------------------------------------------------
subroutine oft_h1_project(fe_rep,field,x)
class(oft_afem_type), target, intent(inout) :: fe_rep
class(fem_interp), intent(inout) :: field
class(oft_vector), intent(inout) :: x
!---
real(r8) :: bcc(1),det,goptmp(3,4),vol
real(r8), pointer :: xloc(:)
real(r8), allocatable :: rop(:)
integer(i4) :: i,jc,m
integer(i4), allocatable :: j(:)
logical :: curved
CLASS(oft_h1_fem), POINTER :: h1_rep
DEBUG_STACK_PUSH
IF(.NOT.oft_3D_h1_cast(h1_rep,fe_rep))CALL oft_abort("Incorrect FE type","oft_h1_project",__FILE__)
!---Initialize vectors to zero
NULLIFY(xloc)
call x%set(0.d0)
call x%get_local(xloc)
!---Integerate over the volume
!$omp parallel private(j,rop,curved,m,goptmp,vol,det,bcc,jc)
allocate(j(h1_rep%nce),rop(h1_rep%nce))
!$omp do schedule(guided)
do i=1,h1_rep%mesh%nc ! Loop over cells
  call h1_rep%ncdofs(i,j) ! Get DOFs
  curved=cell_is_curved(h1_rep%mesh,i) ! Straight cell test
  do m=1,h1_rep%quad%np
    if(curved.OR.m==1)call h1_rep%mesh%jacobian(i,h1_rep%quad%pts(:,m),goptmp,vol)
    det=vol*h1_rep%quad%wts(m)
    call field%interp(i,h1_rep%quad%pts(:,m),goptmp,bcc)
    CALL oft_h1_eval_all(h1_rep,i,h1_rep%quad%pts(:,m),rop)
    do jc=1,h1_rep%nce
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
end subroutine oft_h1_project
!------------------------------------------------------------------------------
!> Construct interpolation matrices for transfer between H^1 finite element
!! spaces
!------------------------------------------------------------------------------
SUBROUTINE h1_setup_interp(ML_h1_rep)
CLASS(oft_ml_fem_type), intent(inout) :: ML_h1_rep
INTEGER(i4) :: i
DEBUG_STACK_PUSH
!---
DO i=ML_h1_rep%minlev+1,ML_h1_rep%nlevels
  CALL ML_h1_rep%set_level(i)
  !---
  if(ML_h1_rep%level==ML_h1_rep%blevel+1)CYCLE
  !---Setup interpolation
  if(ML_h1_rep%current_level%order==1)then
    CALL build_ginterpmatrix(ML_h1_rep%interp_matrices(ML_h1_rep%level)%m)
    CALL ML_h1_rep%interp_matrices(ML_h1_rep%level)%m%assemble
  else
    CALL build_pinterpmatrix(ML_h1_rep%interp_matrices(ML_h1_rep%level)%m)
    CALL ML_h1_rep%interp_matrices(ML_h1_rep%level)%m%assemble
  end if
END DO
DEBUG_STACK_POP
CONTAINS
!------------------------------------------------------------------------------
!> Construct interpolation matrix for transfer between geometric levels
!! of H^1 finite element space
!------------------------------------------------------------------------------
SUBROUTINE build_ginterpmatrix(mat)
class(oft_matrix), pointer, intent(inout) :: mat
INTEGER(i4) :: i,j,k,m,icors,ifine,jb,i_ind(1),j_ind(1)
INTEGER(i4) :: etmp(2),ftmp(3),fetmp(3),ctmp(4),fc,ed
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pmap,emap,fmap
CLASS(oft_afem_type), POINTER :: h1_cors => NULL()
CLASS(oft_afem_type), POINTER :: h1_fine => NULL()
class(oft_mesh), pointer :: cmesh
CLASS(oft_vector), POINTER :: h1_vec_fine,h1_vec_cors
integer(i4) :: jfe(3),jce(6)
integer(i4), pointer :: lfde(:,:),lede(:,:)
real(r8) :: f(4),incr,val,d(3),h_rop(3),goptmp(3,4),v,mop(1)
type(oft_graph_ptr), pointer :: graphs(:,:)
type(oft_graph), POINTER :: interp_graph
!---
if(ML_h1_rep%ml_mesh%level<1)then
  call oft_abort('Invalid mesh level','h1_setup_interp::build_ginterpmatrix',__FILE__)
end if
cmesh=>ML_h1_rep%ml_mesh%meshes(ML_h1_rep%ml_mesh%level-1)
if(cmesh%type/=1)CALL oft_abort("Only supported with tet meshes", &
  "h1_setup_interp::build_ginterpmatrix", __FILE__)
h1_fine=>ML_h1_rep%levels(ML_h1_rep%level)%fe
if(h1_fine%order/=1)then
  call oft_abort('Attempted geometric interpolation for pd > 1','h1_setup_interp::build_ginterpmatrix',__FILE__)
end if
h1_cors=>ML_h1_rep%levels(ML_h1_rep%level-1)%fe
lede=>ML_h1_rep%ml_mesh%inter(ML_h1_rep%ml_mesh%level-1)%lede
lfde=>ML_h1_rep%ml_mesh%inter(ML_h1_rep%ml_mesh%level-1)%lfde
ALLOCATE(ML_h1_rep%interp_graphs(ML_h1_rep%level)%g)
interp_graph=>ML_h1_rep%interp_graphs(ML_h1_rep%level)%g
!---Setup matrix sizes
interp_graph%nr=h1_fine%ne
interp_graph%nrg=h1_fine%global%ne
interp_graph%nc=h1_cors%ne
interp_graph%ncg=h1_cors%global%ne
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
if(interp_graph%kr(1)/=1)call oft_abort('Bad element to element count','h1_setup_interp::build_ginterpmatrix',__FILE__)
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
NULLIFY(h1_vec_fine,h1_vec_cors)
CALL ML_h1_rep%vec_create(h1_vec_fine)
CALL ML_h1_rep%vec_create(h1_vec_cors,ML_h1_rep%level-1)
!---
ALLOCATE(graphs(1,1))
graphs(1,1)%g=>interp_graph
!---
CALL create_matrix(mat,graphs,h1_vec_fine,h1_vec_cors)
CALL h1_vec_fine%delete
CALL h1_vec_cors%delete
DEALLOCATE(graphs,h1_vec_fine,h1_vec_cors)
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
END SUBROUTINE build_ginterpmatrix
!------------------------------------------------------------------------------
!> Construct interpolation matrix for transfer between polynomial levels
!! of H^1 finite element space
!------------------------------------------------------------------------------
SUBROUTINE build_pinterpmatrix(mat)
class(oft_matrix), pointer, intent(inout) :: mat
INTEGER(i4) :: i,j,k,m,icors,ifine,jb,js,jn,i_ind(1),j_ind(1)
INTEGER(i4) :: etmp(2),ftmp(3),fetmp(3),ctmp(4),fc,ed
INTEGER(i4) :: offsetc,offsetf
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: pmap,emap,fmap
REAL(r8) :: f(4),val,mop(1)
CLASS(oft_h1_fem), POINTER :: h1_cors => NULL()
CLASS(oft_h1_fem), POINTER :: h1_fine => NULL()
CLASS(oft_vector), POINTER :: h1_vec_fine,h1_vec_cors
type(oft_graph_ptr), pointer :: graphs(:,:)
class(oft_mesh), pointer :: mesh
type(oft_graph), POINTER :: interp_graph
!---
IF(.NOT.oft_3D_h1_cast(h1_fine,ML_h1_rep%levels(ML_h1_rep%level)%fe))CALL oft_abort("Error casting fine level", "h1_setup_interp::build_pinterpmatrix",__FILE__)
IF(.NOT.oft_3D_h1_cast(h1_cors,ML_h1_rep%levels(ML_h1_rep%level-1)%fe))CALL oft_abort("Error casting coarse level", "h1_setup_interp::build_pinterpmatrix",__FILE__)
mesh=>h1_fine%mesh
ALLOCATE(ML_h1_rep%interp_graphs(ML_h1_rep%level)%g)
interp_graph=>ML_h1_rep%interp_graphs(ML_h1_rep%level)%g
!---Setup matrix sizes
interp_graph%nr=h1_fine%ne
interp_graph%nrg=h1_fine%global%ne
interp_graph%nc=h1_cors%ne
interp_graph%ncg=h1_cors%global%ne
!---Setup Matrix graph
ALLOCATE(interp_graph%kr(interp_graph%nr+1))
interp_graph%kr=0
interp_graph%nnz=h1_cors%ne
ALLOCATE(interp_graph%lc(interp_graph%nnz))
interp_graph%lc=0_i4
!---Construct matrix
do i=1,mesh%np
  interp_graph%kr(i)=1
end do
!---
do i=1,mesh%ne
  do j=1,h1_cors%gstruct(2)
    offsetf=mesh%np+(i-1)*h1_fine%gstruct(2)
    interp_graph%kr(j+offsetf)=1
  end do
end do
!---
do i=1,mesh%nf
  do j=1,h1_cors%gstruct(3)
    offsetf=mesh%np+h1_fine%gstruct(2)*mesh%ne+(i-1)*h1_fine%gstruct(3)
    interp_graph%kr(j+offsetf)=1
  end do
end do
!---
do i=1,mesh%nc
  do j=1,h1_cors%gstruct(4)
    offsetf=mesh%np+h1_fine%gstruct(2)*mesh%ne+h1_fine%gstruct(3)*mesh%nf+(i-1)*h1_fine%gstruct(4)
    interp_graph%kr(j+offsetf)=1
  end do
end do
interp_graph%kr(interp_graph%nr+1)=interp_graph%nnz+1
do i=interp_graph%nr,1,-1 ! cumulative point to point count
  interp_graph%kr(i)=interp_graph%kr(i+1)-interp_graph%kr(i)
end do
if(interp_graph%kr(1)/=1)call oft_abort('Bad element to element count','h1_setup_interp::build_pinterpmatrix',__FILE__)
!---Construct matrix
do i=1,mesh%np
  interp_graph%lc(interp_graph%kr(i))=i
end do
!---
do i=1,mesh%ne
  do j=1,h1_cors%gstruct(2)
    offsetf=mesh%np+(i-1)*h1_fine%gstruct(2)
    offsetc=mesh%np+(i-1)*h1_cors%gstruct(2)
    interp_graph%lc(interp_graph%kr(j+offsetf))=j+offsetc
  end do
end do
!---
do i=1,mesh%nf
  do j=1,h1_cors%gstruct(3)
    offsetf=mesh%np+h1_fine%gstruct(2)*mesh%ne+(i-1)*h1_fine%gstruct(3)
    offsetc=mesh%np+h1_cors%gstruct(2)*mesh%ne+(i-1)*h1_cors%gstruct(3)
    interp_graph%lc(interp_graph%kr(j+offsetf))=j+offsetc
  end do
end do
!---
do i=1,mesh%nc
  do j=1,h1_cors%gstruct(4)
    offsetf=mesh%np+h1_fine%gstruct(2)*mesh%ne+h1_fine%gstruct(3)*mesh%nf+(i-1)*h1_fine%gstruct(4)
    offsetc=mesh%np+h1_cors%gstruct(2)*mesh%ne+h1_cors%gstruct(3)*mesh%nf+(i-1)*h1_cors%gstruct(4)
    interp_graph%lc(interp_graph%kr(j+offsetf))=j+offsetc
  end do
end do
!------------------------------------------------------------------------------
! Construct matrix
!------------------------------------------------------------------------------
NULLIFY(h1_vec_fine,h1_vec_cors)
CALL ML_h1_rep%vec_create(h1_vec_fine)
CALL ML_h1_rep%vec_create(h1_vec_cors,ML_h1_rep%level-1)
!---
ALLOCATE(graphs(1,1))
graphs(1,1)%g=>interp_graph
!---
CALL create_matrix(mat,graphs,h1_vec_fine,h1_vec_cors)
CALL h1_vec_fine%delete
CALL h1_vec_cors%delete
DEALLOCATE(graphs,h1_vec_fine,h1_vec_cors)
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
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
  DO j=1,h1_cors%gstruct(2)
    i_ind=j+mesh%np+(i-1)*h1_fine%gstruct(2)
    j_ind=j+mesh%np+(i-1)*h1_cors%gstruct(2)
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
  DO j=1,h1_cors%gstruct(3)
    i_ind=j+mesh%np+h1_fine%gstruct(2)*mesh%ne+(i-1)*h1_fine%gstruct(3)
    j_ind=j+mesh%np+h1_cors%gstruct(2)*mesh%ne+(i-1)*h1_cors%gstruct(3)
    mop=1.d0
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
  END DO
END DO
deallocate(fmap)
!---
DO i=1,mesh%nc
  DO j=1,h1_cors%gstruct(4)
    i_ind=j+mesh%np+h1_fine%gstruct(2)*mesh%ne+h1_fine%gstruct(3)*mesh%nf+(i-1)*h1_fine%gstruct(4)
    j_ind=j+mesh%np+h1_cors%gstruct(2)*mesh%ne+h1_cors%gstruct(3)*mesh%nf+(i-1)*h1_cors%gstruct(4)
    mop=1.d0
    CALL mat%add_values(i_ind,j_ind,mop,1,1)
  END DO
END DO
END SUBROUTINE build_pinterpmatrix
END SUBROUTINE h1_setup_interp
!------------------------------------------------------------------------------
!> Transfer a base level H^1 scalar field to the next MPI level
!------------------------------------------------------------------------------
subroutine h1_base_pop(self,acors,afine)
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
end subroutine h1_base_pop
!------------------------------------------------------------------------------
!> Transfer a MPI level H^1 scalar field to the base level
!------------------------------------------------------------------------------
subroutine h1_base_push(self,afine,acors)
class(oft_ml_fe_vecspace), intent(inout) :: self
class(oft_vector), intent(inout) :: afine !< Vector to transfer
class(oft_vector), intent(inout) :: acors !< Fine vector from transfer
integer(i4), pointer :: lptmp(:)
integer(i4) :: i,j,ierr
real(r8), pointer, dimension(:) :: alias,array_c,array_f
CLASS(oft_afem_type), POINTER :: h1_fine => NULL()
DEBUG_STACK_PUSH
!---
lptmp=>self%ML_FE_rep%ml_mesh%meshes(self%ML_FE_rep%ml_mesh%nbase+1)%base%lp
CALL acors%get_local(array_c)
CALL afine%get_local(array_f)
h1_fine=>self%ML_FE_rep%levels(self%ML_FE_rep%level+1)%fe
!---
allocate(alias(acors%n))
alias=0.d0
!$omp parallel do
do i=1,afine%n
  if(h1_fine%linkage%be(i))cycle
  alias(lptmp(i))=array_f(i)
end do
!$omp parallel do private(j)
do i=1,h1_fine%linkage%nbe
  j=h1_fine%linkage%lbe(i)
  if(.NOT.h1_fine%linkage%leo(i))cycle
  alias(lptmp(j))=array_f(j)
end do
!---Global reduction over all processors
array_c=oft_mpi_sum(alias,acors%n)
call acors%restore_local(array_c)
deallocate(alias,array_c,array_f)
DEBUG_STACK_POP
end subroutine h1_base_push
!------------------------------------------------------------------------------
!> Compute eigenvalues and smoothing coefficients for the operator H^1::LOP
!------------------------------------------------------------------------------
SUBROUTINE h1_lop_eigs(ML_h1_rep,minlev)
type(oft_ml_fem_type), target, intent(inout) :: ML_h1_rep
INTEGER(i4), INTENT(in) :: minlev
#ifdef HAVE_ARPACK
INTEGER(i4) :: i
REAL(r8) :: lam0
REAL(r8), ALLOCATABLE :: df(:)
CLASS(oft_vector), POINTER :: u
TYPE(oft_irlm_eigsolver) :: arsolver
CLASS(oft_matrix), POINTER :: md => NULL()
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_h1_fem), POINTER :: h1_obj
TYPE(oft_h1_zerob), TARGET :: bc_tmp
DEBUG_STACK_PUSH
bc_tmp%ML_H1_rep=>ML_h1_rep
!------------------------------------------------------------------------------
! Compute optimal smoother coefficients
!------------------------------------------------------------------------------
IF(oft_env%head_proc)WRITE(*,*)'Optimizing Jacobi damping for H^1::LOP'
ALLOCATE(df(ML_h1_rep%nlevels))
df=0.d0
DO i=minlev,ML_h1_rep%nlevels
  CALL ML_h1_rep%set_level(i)
  !---Create fields
  CALL ML_h1_rep%vec_create(u)
  !---Get Ev range
  NULLIFY(lop)
  IF(.NOT.oft_3D_h1_cast(h1_obj,ML_h1_rep%current_level))CALL oft_abort("Error getting current FE rep","h1_lop_eigs",__FILE__)
  CALL oft_h1_getlop(h1_obj,lop,'grnd')
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
  DO i=1,ML_h1_rep%nlevels-1
    WRITE(*,'(F5.3,A)',ADVANCE='NO')df(i),', '
  END DO
  WRITE(*,'(F5.3,A)')df(ML_h1_rep%nlevels)
END IF
DEALLOCATE(df)
DEBUG_STACK_POP
#else
CALL oft_abort("Subroutine requires ARPACK", "lag_lop_eigs", __FILE__)
#endif
END SUBROUTINE h1_lop_eigs
!------------------------------------------------------------------------------
!> Compute eigenvalues and smoothing coefficients for the operator H^1::LOP
!------------------------------------------------------------------------------
SUBROUTINE h1_getlop_pre(ML_h1_rep,pre,mats,bc_type,level,nlevels)
type(oft_ml_fem_type), target, intent(inout) :: ML_h1_rep
CLASS(oft_solver), POINTER, INTENT(out) :: pre
TYPE(oft_matrix_ptr), POINTER, INTENT(inout) :: mats(:)
CHARACTER(LEN=*), INTENT(in) :: bc_type !< Boundary condition
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
TYPE(oft_h1_zerob), POINTER :: zerob_bc
TYPE(oft_h1_zerogrnd), POINTER :: zerogrnd_bc
TYPE(oft_ml_fe_vecspace), POINTER :: tmp_vecspace
!---
TYPE(xml_node), POINTER :: pre_node
#ifdef HAVE_XML
integer(i4) :: nnodes
TYPE(xml_node), POINTER :: h1_node
#endif
DEBUG_STACK_PUSH
!---
minlev=1
toplev=ML_h1_rep%level
levin=ML_h1_rep%level
IF(PRESENT(level))toplev=level
IF(PRESENT(nlevels))minlev=toplev-nlevels+1
nl=toplev-minlev+1
!---
IF(minlev<ML_h1_rep%minlev)CALL oft_abort('Minimum level is < minlev','h1_getlop_pre',__FILE__)
IF(toplev>ML_h1_rep%nlevels)CALL oft_abort('Maximum level is > nlevels','h1_getlop_pre',__FILE__)
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
  CALL ML_h1_rep%set_level(minlev+(i-1))
  levels(i)=minlev+(i-1)
  df(i)=df_lop(levels(i))
  nu(i)=nu_lop(levels(i))
  IF(df(i)<-1.d90)THEN
    WRITE(lev_char,'(I2.2)')levels(i)
    CALL oft_abort('Smoother values not set for level: '//lev_char,'h1_getlop_pre',__FILE__)
  END IF
  !---
  IF(create_mats)THEN
    NULLIFY(mats(i)%M)
    CALL oft_h1_getlop(ML_h1_rep%current_level,mats(i)%M,'grnd')
  END IF
  IF(i>1)ml_int(i-1)%M=>ML_h1_rep%interp_matrices(ML_h1_rep%level)%m
END DO
CALL ML_h1_rep%set_level(levin)
!------------------------------------------------------------------------------
! Search for XML-spec
!------------------------------------------------------------------------------
NULLIFY(pre_node)
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml))THEN
  !---Look for H^1 FE node
  CALL xml_get_element(oft_env%xml,"h1",h1_node,ierr)
  IF(ierr==0)CALL xml_get_element(h1_node,"lop",pre_node,ierr)
END IF
#endif
!------------------------------------------------------------------------------
! Setup preconditioner
!------------------------------------------------------------------------------
NULLIFY(pre)
ALLOCATE(tmp_vecspace)
tmp_vecspace%ML_FE_rep=>ML_h1_rep
tmp_vecspace%base_pop=>h1_base_pop
tmp_vecspace%base_push=>h1_base_push
SELECT CASE(TRIM(bc_type))
  CASE("zerob")
    ALLOCATE(zerob_bc)
    zerob_bc%ML_H1_rep=>ML_h1_rep
    CALL create_mlpre(pre,mats(1:nl),levels,nlevels=nl,ml_vecspace=tmp_vecspace, &
      bc=zerob_bc,stype=1,df=df,nu=nu,xml_root=pre_node)
  CASE("grnd")
    ALLOCATE(zerogrnd_bc)
    zerogrnd_bc%ML_H1_rep=>ML_h1_rep
    CALL create_mlpre(pre,mats(1:nl),levels,nlevels=nl,ml_vecspace=tmp_vecspace, &
      bc=zerogrnd_bc,stype=1,df=df,nu=nu,xml_root=pre_node)
  CASE DEFAULT
    CALL oft_abort("Unknown BC type","h1_getlop_pre",__FILE__)
END SELECT
!------------------------------------------------------------------------------
! Cleanup
!------------------------------------------------------------------------------
DEALLOCATE(ml_int,levels,df,nu)
DEBUG_STACK_POP
END SUBROUTINE h1_getlop_pre
end module oft_h1_operators
