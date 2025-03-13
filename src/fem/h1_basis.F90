!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_h1_basis.F90
!
!> @defgroup doxy_oft_h1 H^1 FE for conforming FE sequence
!! Scalar H^1 finite element implementation for the Open FUSION Toolkit
!! @ingroup doxy_oft_fem
!
!> Base scalar H^1 FE class and basis evaluation
!! - FE Construction
!! - Basis evaluation
!!   - Interpolation
!!   - Gradient
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_h1
!---------------------------------------------------------------------------
MODULE oft_h1_basis
USE oft_base
USE oft_lag_poly
USE oft_mesh_type, ONLY: oft_mesh, oft_bmesh
USE oft_mesh_local_util, ONLY: mesh_local_orient, oriented_cell, &
  oriented_edges, oriented_faces
USE oft_hexmesh_type, ONLY: hex_get_bary, hex_get_bary_gop, &
  hex_bary_pfcoords, hex_bary_efcoords, hex_bary_ecoords, hex_bary_fcoords
USE multigrid, ONLY: multigrid_mesh, multigrid_level
USE oft_la_utils, ONLY: oft_matrix, oft_graph
USE fem_base, ONLY: oft_fem_type, oft_ml_fem_type, oft_bfem_type, oft_afem_type
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
type, extends(oft_fem_type) :: oft_h1_fem
  INTEGER(i4), POINTER, DIMENSION(:,:) :: indsf => NULL() !< Needs docs
  INTEGER(i4), POINTER, DIMENSION(:,:) :: indsc => NULL() !< Needs docs
end type oft_h1_fem
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
type, extends(oft_bfem_type) :: oft_h1_bfem
  INTEGER(i4), POINTER, DIMENSION(:,:) :: indsf => NULL() !< Needs docs
end type oft_h1_bfem
!---Global Variables
integer(i4), parameter :: oft_h1_id = 2 !< FE type ID
contains
!------------------------------------------------------------------------------
!> Cast abstract FE type to 3D H^1 finite element type
!!
!! The source matrix must be @ref oft_h1_fem or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION oft_3D_h1_cast(self,source) RESULT(success)
CLASS(oft_h1_fem), POINTER, INTENT(out) :: self !< Reference to source object with desired class
CLASS(oft_afem_type), TARGET, INTENT(in) :: source !< Source object to reference
LOGICAL :: success !< Cast success flag
DEBUG_STACK_PUSH
SELECT TYPE(source)
  CLASS IS(oft_h1_fem)
    self=>source
    success=.TRUE.
  CLASS DEFAULT
    NULLIFY(self)
    success=.FALSE.
END SELECT
DEBUG_STACK_POP
END FUNCTION oft_3D_h1_cast
!------------------------------------------------------------------------------
!> Cast abstract FE type to 2D H^1 finite element type
!!
!! The source matrix must be @ref oft_h1_bfem or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION oft_2D_h1_cast(self,source) RESULT(success)
CLASS(oft_h1_bfem), POINTER, INTENT(out) :: self !< Reference to source object with desired class
CLASS(oft_afem_type), TARGET, INTENT(in) :: source !< Source object to reference
LOGICAL :: success !< Cast success flag
DEBUG_STACK_PUSH
SELECT TYPE(source)
  CLASS IS(oft_h1_bfem)
    self=>source
    success=.TRUE.
  CLASS DEFAULT
    NULLIFY(self)
    success=.FALSE.
END SELECT
DEBUG_STACK_POP
END FUNCTION oft_2D_h1_cast
!---------------------------------------------------------------------------
!> Construct H^1 scalar FE on each mesh level
!!
!! @note Highest supported representation is Quartic
!---------------------------------------------------------------------------
subroutine oft_h1_setup(mg_mesh,order,ML_h1_obj,ML_bh1_obj,minlev)
type(multigrid_mesh), target, intent(inout) :: mg_mesh
integer(i4), intent(in) :: order !< Order of representation desired
TYPE(oft_ml_fem_type), optional, INTENT(inout) :: ML_h1_obj
TYPE(oft_ml_fem_type), optional, INTENT(inout) :: ML_bh1_obj
integer(i4), optional, intent(in) :: minlev !< Lowest level to construct
integer(i4) :: i,j,k,nlevels,minlev_out
REAL(r8), POINTER, DIMENSION(:) :: xnodes
DEBUG_STACK_PUSH
minlev_out=1
IF(PRESENT(minlev))minlev_out=minlev
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(A)')'**** Creating H^1 FE space'
  WRITE(*,'(2X,A,I4)')'Order  = ',order
  WRITE(*,'(2X,A,I4)')'Minlev = ',minlev_out
END IF
!---Allocate multigrid operators
nlevels=mg_mesh%mgdim+(order-1)
IF(minlev_out<0)minlev_out=nlevels
IF(PRESENT(ML_h1_obj))THEN
  IF(.NOT.ASSOCIATED(mg_mesh%meshes))CALL oft_abort("No volume mesh available for 3D elements","oft_h1_setup",__FILE__)
  ML_h1_obj%nlevels=nlevels
  ML_h1_obj%minlev=minlev_out
  ML_h1_obj%ml_mesh=>mg_mesh
ELSE
  IF(.NOT.PRESENT(ML_bh1_obj))THEN
    WRITE(*,*)'No H1 FE objects requested, returning'
    DEBUG_STACK_POP
    RETURN
  END IF
END IF
IF(PRESENT(ML_bh1_obj))THEN
  IF(.NOT.ASSOCIATED(mg_mesh%smeshes))CALL oft_abort("No surface mesh available for 2D elements","oft_h1_setup",__FILE__)
  ML_bh1_obj%nlevels=nlevels
  ML_bh1_obj%minlev=minlev_out
  ML_bh1_obj%ml_mesh=>mg_mesh
END IF
! Set linear elements
do i=1,mg_mesh%mgdim-1
  IF(i<minlev_out)CYCLE
  CALL multigrid_level(mg_mesh,i)
  IF(PRESENT(ML_h1_obj))THEN
    CALL oft_h1_setup_vol(ML_h1_obj%levels(i)%fe,mg_mesh%mesh,1)
    IF(mg_mesh%level==mg_mesh%nbase)ML_h1_obj%blevel=i
    CALL ML_h1_obj%set_level(i)
    IF(mg_mesh%level==mg_mesh%nbase)ML_h1_obj%blevel=i
  END IF
  IF(PRESENT(ML_bh1_obj))THEN
    CALL oft_h1_setup_surf(ML_bh1_obj%levels(i)%fe,mg_mesh%smesh,1)
    IF(mg_mesh%level==mg_mesh%nbase)ML_bh1_obj%blevel=i
  END IF
end do
call multigrid_level(mg_mesh,mg_mesh%mgdim)
! Set high order elements
do i=1,order
  IF(mg_mesh%mgdim+i-1<minlev_out)CYCLE
  IF(PRESENT(ML_h1_obj))THEN
    CALL oft_h1_setup_vol(ML_h1_obj%levels(mg_mesh%mgdim+i-1)%fe,mg_mesh%mesh,i)
    CALL ML_h1_obj%set_level(mg_mesh%mgdim+i-1)
  END IF
  IF(PRESENT(ML_bh1_obj))THEN
    CALL oft_h1_setup_surf(ML_bh1_obj%levels(mg_mesh%mgdim+i-1)%fe,mg_mesh%smesh,i)
  END IF
end do
IF(PRESENT(ML_h1_obj))CALL ML_h1_obj%set_level(ML_h1_obj%nlevels)
IF(PRESENT(ML_bh1_obj))CALL ML_bh1_obj%set_level(ML_bh1_obj%nlevels)
IF(oft_env%head_proc)WRITE(*,*)
DEBUG_STACK_POP
end subroutine oft_h1_setup
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine oft_h1_setup_vol(self,tmesh,order)
class(oft_afem_type), pointer, intent(out) :: self !< Needs docs
class(oft_mesh), target, intent(in) :: tmesh !< Needs docs
integer(i4), intent(in) :: order !< Order of representation desired
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2A)')oft_indent,'Creating 3D H^1 FE space'
  WRITE(*,'(A,2X,A,I4)')oft_indent,'Order  = ',order
END IF
CALL oft_increase_indent
!---
ALLOCATE(oft_h1_fem::self)
SELECT TYPE(self)
CLASS IS(oft_h1_fem)
  self%mesh=>tmesh
  self%order=order
  self%dim=1
  self%type=oft_h1_id
  IF(self%mesh%type==3)THEN
    CALL hpoly_2d_grid(self%order-1, self%indsf)
    CALL hpoly_3d_grid(self%order-1, self%indsc)
    select case(self%order)
      case(1)
        self%gstruct=(/1,0,0,0/)
      case(2)
        self%gstruct=(/1,1,1,1/)
      case(3)
        self%gstruct=(/1,2,4,8/)
      case(4)
        self%gstruct=(/1,3,9,27/)
      case(5)
        self%gstruct=(/1,4,16,64/)
      case default
        call oft_abort('Invalid polynomial degree (npmax=4 for hex grids)','oft_h1_setup_vol',__FILE__)
    end select
  ELSE
    select case(self%order)
      case(1)
        self%gstruct=(/1,0,0,0/)
      case(2)
        self%gstruct=(/1,1,0,0/)
      case(3)
        self%gstruct=(/1,2,1,0/)
      case(4)
        self%gstruct=(/1,3,3,1/)
      case(5)
        self%gstruct=(/1,4,6,4/)
      case default
        call oft_abort('Invalid polynomial degree (npmax=5)','oft_h1_setup_vol',__FILE__)
    end select
  END IF
CLASS DEFAULT
  CALL oft_abort("Error allocating H^1 FE object","oft_h1_setup_vol",__FILE__)
END SELECT
call self%setup(self%order*2+1)
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine oft_h1_setup_vol
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine oft_h1_setup_surf(self,tmesh,order)
class(oft_afem_type), pointer, intent(out) :: self !< Needs docs
class(oft_bmesh), target, intent(in) :: tmesh !< Needs docs
integer(i4), intent(in) :: order !< Order of representation desired
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2A)')oft_indent,'Creating 2D H^1 FE space'
  WRITE(*,'(A,2X,A,I4)')oft_indent,'Order  = ',order
END IF
CALL oft_increase_indent
!---
ALLOCATE(oft_h1_bfem::self)
SELECT TYPE(self)
CLASS IS(oft_h1_bfem)
  self%mesh=>tmesh
  self%order=order
  self%dim=1
  self%type=oft_h1_id
  IF(self%mesh%type==3)THEN
    select case(self%order)
      case(1)
        self%gstruct=(/1,0,0/)
      case(2)
        self%gstruct=(/1,1,1/)
      case(3)
        self%gstruct=(/1,2,4/)
      case(4)
        self%gstruct=(/1,3,9/)
      case(5)
        self%gstruct=(/1,4,16/)
      case default
        call oft_abort('Invalid polynomial degree (npmax=5)','oft_h1_setup',__FILE__)
    end select
  ELSE
    select case(self%order)
      case(1)
        self%gstruct=(/1,0,0/)
      case(2)
        self%gstruct=(/1,1,0/)
      case(3)
        self%gstruct=(/1,2,1/)
      case(4)
        self%gstruct=(/1,3,3/)
      case(5)
        self%gstruct=(/1,4,6/)
      case default
        call oft_abort('Invalid polynomial degree (npmax=5)','oft_h1_setup',__FILE__)
    end select
  END IF
CLASS DEFAULT
  CALL oft_abort("Error allocating H^1 FE object","oft_h1_setup_surf",__FILE__)
END SELECT
call self%setup(self%order*2+1)
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine oft_h1_setup_surf
!---------------------------------------------------------------------------
!> Evaluate H^1 interpolation function
!!
!! @note Evaluation is performed in logical coordinates with the resulting
!! vector in physical coordinates
!---------------------------------------------------------------------------
subroutine oft_h1_eval(self,cell,dof,f,val)
class(oft_h1_fem), intent(in) :: self !< H^1 type for evaluation
integer(i4), intent(in) :: cell !< Cell for evaluation
integer(i4), intent(in) :: dof !< Element to evaluate
real(r8), intent(in) :: f(:) !< Position in cell in logical space
real(r8), intent(out) :: val !< Value of interpolation function (dof) at point (f)
integer(i4) :: i,ed,etmp(2),fc,ftmp(3),finds(9),ind,fhtmp(4),inds(4)
real(r8) :: fhex(6)
DEBUG_STACK_PUSH
IF(self%mesh%type==3)THEN
  fhex=hex_get_bary(f)
  select case(self%cmap(dof)%type)
    case(1)
      val=1.d0
      DO i=1,3
        val=val*fhex(hex_bary_pfcoords(i,self%cmap(dof)%el))
      END DO
    case(2)
      etmp=hex_bary_ecoords(:,self%cmap(dof)%el)
      call orient_list2(self%mesh%lce(self%cmap(dof)%el,cell), etmp)
      !
      val=hpoly_bary(self%cmap(dof)%ind, fhex(etmp(1)))*fhex(etmp(2))
      ! Fixed is always linear
      DO i=1,2
        val=val*fhex(hex_bary_efcoords(i,self%cmap(dof)%el))
      END DO
    case(3)
      fhtmp=hex_bary_fcoords(:,self%cmap(dof)%el)
      IF(self%mesh%lcfo(self%cmap(dof)%el,cell)<0)fhtmp=fhtmp((/2,3,4,1/))
      call orient_listn(self%mesh%lcfo(self%cmap(dof)%el,cell), fhtmp, 4_i4)
      !
      val=hpoly_2d(self%indsf(:,self%cmap(dof)%ind), fhex(fhtmp(1:2)))
      ! Fixed is always linear
      val=val*fhex(self%cmap(dof)%el)*fhex(fhtmp(3))*fhex(fhtmp(4))
    case(4)
      val=hpoly_3d(self%indsc(:,self%cmap(dof)%ind), fhex(1:3))* &
        fhex(4)*fhex(5)*fhex(6)
  end select
ELSE
  IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
  select case(self%cmap(dof)%type)
    case(1)
      val=f(self%cmap(dof)%el)
    case(2)
      etmp=oriented_edges(:,self%cmap(dof)%el)
      call oft_h1_evale(self%order,etmp,self%cmap(dof)%ind,f,val)
    case(3)
      ftmp=oriented_faces(:,self%cmap(dof)%el)
      call oft_h1_evalf(self%order,ftmp,self%cmap(dof)%ind,f,val)
    case(4)
      call oft_h1_evalc(self%order,self%cmap(dof)%ind,f,val)
    case default
      call oft_abort('Invalid geometry type!','oft_h1_eval',__FILE__)
  end select
END IF
DEBUG_STACK_POP
end subroutine oft_h1_eval
!---------------------------------------------------------------------------
!> Evaluate H^1 interpolation function on the boundary
!!
!! @note Evaluation is performed in logical coordinates with the resulting
!! vector in physical coordinates
!---------------------------------------------------------------------------
subroutine oft_bh1_eval(self,face,dof,f,val)
class(oft_bfem_type), intent(in) :: self !< H^1 type for evaluation (bfem)
integer(i4), intent(in) :: face !< Face for evaluation
integer(i4), intent(in) :: dof !< Element to evaluate
real(r8), intent(in) :: f(:) !< Position on face in logical space [4]
real(r8), intent(out) :: val !< Value of interpolation function (dof) at point (f)
integer(i4) :: ed,etmp(2),fc,ftmp(3)
DEBUG_STACK_PUSH
select case(self%cmap(dof)%type)
  case(1)
    val=f(self%cmap(dof)%el)
  case(2)
    ed=self%mesh%lce(self%cmap(dof)%el,face)
    etmp=self%mesh%cell_ed(:,self%cmap(dof)%el)
    call orient_list2(ed,etmp)
    call oft_h1_evale(self%order,etmp,self%cmap(dof)%ind,f,val)
  case(3)
    fc=self%mesh%lco(face)
    ftmp=(/1,2,3/)
    call orient_listn(fc,ftmp,3_i4)
    call oft_h1_evalf(self%order,ftmp,self%cmap(dof)%ind,f,val)
end select
DEBUG_STACK_POP
end subroutine oft_bh1_eval
!---------------------------------------------------------------------------
!> Evaluate point based interpolation functions
!---------------------------------------------------------------------------
subroutine oft_h1_evalp(order,pt,f,val)
integer(i4), intent(in) :: order !< Needs docs
integer(i4), intent(in) :: pt !< Needs docs
real(r8), intent(in) :: f(:) !< Needs docs
real(r8), intent(out) :: val !< Needs docs
DEBUG_STACK_PUSH
val=f(pt)
DEBUG_STACK_POP
end subroutine oft_h1_evalp
!---------------------------------------------------------------------------
!> Evaluate edge based interpolation functions
!---------------------------------------------------------------------------
subroutine oft_h1_evale(order,ed,dof,f,val)
integer(i4), intent(in) :: order !< Needs docs
integer(i4), intent(in) :: ed(2) !< Needs docs
integer(i4), intent(in) :: dof !< Needs docs
real(r8), intent(in) :: f(:) !< Needs docs
real(r8), intent(out) :: val !< Needs docs
real(r8) :: x1,x2
DEBUG_STACK_PUSH
x1=f(ed(1)); x2=f(ed(2))
select case(dof)
  case(1)
    val=-2.d0*x1*x2
  case(2)
    val=2.d0*x1*x2*(-x1 + x2)
  case(3)
    val=2.d0*x1*x2*(-x1**2 + 3.d0*x1*x2 - x2**2)
  case(4)
    val=2.d0*x1*x2*(-x1**3 + 6.d0*(x1**2)*x2 - 6.d0*x1*(x2**2) + x2**3)
  case default
    call oft_abort('Invalid DOF!','oft_h1_evale',__FILE__)
end select
DEBUG_STACK_POP
end subroutine oft_h1_evale
!---------------------------------------------------------------------------
!> Evaluate face based interpolation functions
!---------------------------------------------------------------------------
subroutine oft_h1_evalf(order,fc,dof,f,val)
integer(i4), intent(in) :: order !< Needs docs
integer(i4), intent(in) :: fc(3) !< Needs docs
integer(i4), intent(in) :: dof !< Needs docs
real(r8), intent(in) :: f(:) !< Needs docs
real(r8), intent(out) :: val !< Needs docs
real(r8) :: x1,x2,x3
DEBUG_STACK_PUSH
x1=f(fc(1)); x2=f(fc(2)); x3=f(fc(3))
select case(dof)
  case(1)
    val=-2.d0*x1*x2*x3
  case(2)
    val=2.d0*x1*x2*x3*(x1 + x2 - x3)
  case(3)
    val=2.d0*x1*x2*x3*(-x1 + x2)
  case(4)
    val=2.d0*x1*x2*x3*(-x1**2 - 2.d0*x1*x2 + 4.d0*x1*x3 - x2**2 + 4.d0*x2*x3 - x3**2)
  case(5)
    val=2.d0*x1*x2*x3*(x1**2 - x1*x3 - x2**2 + x2*x3)
  case(6)
    val=2.d0*x1*x2*x3*(-x1**2 + 3.d0*x1*x2 - x2**2)
  case default
    call oft_abort('Invalid DOF!','oft_h1_evalf',__FILE__)
end select
DEBUG_STACK_POP
end subroutine oft_h1_evalf
!---------------------------------------------------------------------------
!> Evaluate cell based interpolation functions
!---------------------------------------------------------------------------
subroutine oft_h1_evalc(order,dof,f,val)
integer(i4), intent(in) :: order !< Needs docs
integer(i4), intent(in) :: dof !< Needs docs
real(r8), intent(in) :: f(:) !< Needs docs
real(r8), intent(out) :: val !< Needs docs
real(r8) :: x1,x2,x3,x4
DEBUG_STACK_PUSH
x1=f(1); x2=f(2); x3=f(3); x4=f(4)
select case(dof)
  case(1)
    val=-2.d0*x1*x2*x3*x4
  case(2)
    val=2.d0*x1*x2*x3*x4*(x1 + x2 + x3 - x4)
  case(3)
    val=2.d0*x1*x2*x3*x4*(x1 + x2 - x3)
  case(4)
    val=2.d0*x1*x2*x3*x4*(-x1 + x2)
  case default
    call oft_abort('Invalid DOF!','oft_h1_evalc',__FILE__)
end select
DEBUG_STACK_POP
end subroutine oft_h1_evalc
!---------------------------------------------------------------------------
!> Evaluate all H^1 interpolation functions
!!
!! @note Evaluation is performed in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_h1_eval_all(self,cell,f,rop)
class(oft_h1_fem), intent(in) :: self !< H^1 type for evaluation
integer(i4), intent(in) :: cell !< Cell for evaluation
real(r8), intent(in) :: f(:) !< Position in cell in logical space
real(r8), contiguous, intent(out) :: rop(:) !< Value of interpolation functions at point (f) [ncdofs]
integer(i4) :: i,j,etmp(2),fhtmp(4),offset
real(r8) :: fhex(6),vtmp
DEBUG_STACK_PUSH
IF(self%mesh%type==3)THEN
  ! DO i=1,self%nce
  !   call oft_h1_eval(self,cell,i,f,rop(i))
  ! END DO
  fhex=hex_get_bary(f)
  DO i=1,8
    rop(i)=PRODUCT(fhex(hex_bary_pfcoords(:,i)))
  END DO
  IF(self%gstruct(2)>0)THEN
    !---Edges
    DO i=1,12
      offset=(i-1)*self%gstruct(2)+8
      etmp=hex_bary_ecoords(:,i)
      call orient_list2(self%mesh%lce(i,cell), etmp)
      vtmp=PRODUCT(fhex(hex_bary_efcoords(:,i)))*fhex(etmp(2))
      DO j=1,self%gstruct(2)
        rop(offset+j)=hpoly_bary(j, fhex(etmp(1)))*vtmp
      END DO
    END DO
  END IF
  IF(self%gstruct(3)>0)THEN
    !---Faces
    DO i=1,6
      offset=(i-1)*self%gstruct(3)+12*self%gstruct(2)+8
      fhtmp=hex_bary_fcoords(:,i)
      IF(self%mesh%lcfo(i,cell)<0)fhtmp=fhtmp((/2,3,4,1/))
      call orient_listn(self%mesh%lcfo(i,cell), fhtmp, 4_i4)
      vtmp=fhex(i)*fhex(fhtmp(3))*fhex(fhtmp(4))
      DO j=1,self%gstruct(3)
        rop(offset+j)=hpoly_2d(self%indsf(1:2,j), fhex(fhtmp(1:2)))*vtmp
      END DO
    END DO
  END IF
  IF(self%gstruct(4)>0)THEN
    !---Cell
    offset=6*self%gstruct(3)+12*self%gstruct(2)+8
    vtmp=fhex(4)*fhex(5)*fhex(6)
    DO j=1,self%gstruct(4)
      rop(offset+j)=hpoly_3d(self%indsc(1:3,j), fhex(1:3))*vtmp
    END DO
  END IF
ELSE
  IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
  select case(self%order)
    case(1)
      rop=f
    case(2)
      call oft_h1_eval_all2()!self,cell,f,rop)
    case(3)
      call oft_h1_eval_all3()!self,cell,f,rop)
    case(4)
      call oft_h1_eval_all4()!self,cell,f,rop)
    case(5)
      call oft_h1_eval_all5()!self,cell,f,rop)
    case default ! Fall back to normal evaluation
      DO i=1,self%nce
        call oft_h1_eval(self,cell,i,f,rop(i))
      END DO
  end select
END IF
DEBUG_STACK_POP
contains
!---------------------------------------------------------------------------
!> Evaluate all H^1 interpolation functions (quadratic)
!!
!! @note Evaluation is performed in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_h1_eval_all2()!self,cell,f,rop)
! class(oft_fem_type), intent(in) :: self !< H^1 type for evaluation
! integer(i4), intent(in) :: cell !< Cell for evaluation
! real(r8), intent(in) :: f(4) !< Position in cell in logical space
! real(r8), intent(out) :: rop(10) !< Value of interpolation functions at point (f) [ncdofs]
integer(i4) :: i
real(r8) :: x1,x2
DEBUG_STACK_PUSH
DO i=1,4
  rop(i) = f(i)
END DO
DO i=1,6
  x1 = f(self%mesh%cell_ed(1,i)); x2 = f(self%mesh%cell_ed(2,i))
  rop(i+4) = -2.d0*x1*x2
END DO
DEBUG_STACK_POP
end subroutine oft_h1_eval_all2
!---------------------------------------------------------------------------
!> Evaluate all H^1 interpolation functions (cubic)
!!
!! @note Evaluation is performed in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_h1_eval_all3()!self,cell,f,rop)
! class(oft_fem_type), intent(in) :: self !< H^1 type for evaluation
! integer(i4), intent(in) :: cell !< Cell for evaluation
! real(r8), intent(in) :: f(4) !< Position in cell in logical space
! real(r8), intent(out) :: rop(20) !< Value of interpolation functions at point (f) [ncdofs]
integer(i4) :: i
real(r8) :: x1,x2,x3
DEBUG_STACK_PUSH
DO i=1,4
  rop(i) = f(i)
END DO
DO i=1,6
  x1 = f(oriented_edges(1,i)); x2 = f(oriented_edges(2,i))
  rop((i-1)*2+5)  = -2.d0*x1*x2
  rop((i-1)*2+6) = 2.d0*x1*x2*(-x1 + x2)
END DO
DO i=1,4
  x1 = f(self%mesh%cell_fc(1,i))
  x2 = f(self%mesh%cell_fc(2,i))
  x3 = f(self%mesh%cell_fc(3,i))
  rop(i+16) = -2.d0*x1*x2*x3
END DO
DEBUG_STACK_POP
end subroutine oft_h1_eval_all3
!---------------------------------------------------------------------------
!> Evaluate all H^1 interpolation functions (quartic)
!!
!! @note Evaluation is performed in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_h1_eval_all4()!self,cell,f,rop)
! class(oft_fem_type), intent(in) :: self !< H^1 type for evaluation
! integer(i4), intent(in) :: cell !< Cell for evaluation
! real(r8), intent(in) :: f(4) !< Position in cell in logical space
! real(r8), intent(out) :: rop(35) !< Value of interpolation functions at point (f) [ncdofs]
integer(i4) :: i
real(r8) :: x1,x2,x3,x4
DEBUG_STACK_PUSH
DO i=1,4
  rop(i) = f(i)
END DO
DO i=1,6
  x1 = f(oriented_edges(1,i)); x2 = f(oriented_edges(2,i))
  rop((i-1)*3+5)  = -2.d0*x1*x2
  rop((i-1)*3+6) = 2.d0*x1*x2*(-x1 + x2)
  rop((i-1)*3+7) = 2.d0*x1*x2*(-x1**2 + 3.d0*x1*x2 - x2**2)
END DO
DO i=1,4
  x1 = f(oriented_faces(1,i)); x2 = f(oriented_faces(2,i)); x3 = f(oriented_faces(3,i))
  rop((i-1)*3+23) = -2.d0*x1*x2*x3
  rop((i-1)*3+24) = 2.d0*x1*x2*x3*(x1 + x2 - x3)
  rop((i-1)*3+25) = 2.d0*x1*x2*x3*(-x1 + x2)
END DO
x1 = f(1); x2 = f(2); x3 = f(3); x4 = f(4)
rop(35) = -2.d0*x1*x2*x3*x4
DEBUG_STACK_POP
end subroutine oft_h1_eval_all4
!---------------------------------------------------------------------------
!> Evaluate all H^1 interpolation functions (quartic)
!!
!! @note Evaluation is performed in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_h1_eval_all5()!self,cell,f,rop)
! class(oft_fem_type), intent(in) :: self !< H^1 type for evaluation
! integer(i4), intent(in) :: cell !< Cell for evaluation
! real(r8), intent(in) :: f(4) !< Position in cell in logical space
! real(r8), intent(out) :: rop(56) !< Value of interpolation functions at point (f) [ncdofs]
integer(i4) :: i
real(r8) :: x1,x2,x3,x4
DEBUG_STACK_PUSH
DO i=1,4
  rop(i) = f(i)
END DO
DO i=1,6
  x1 = f(oriented_edges(1,i)); x2 = f(oriented_edges(2,i))
  rop((i-1)*4+5)  = -2.d0*x1*x2
  rop((i-1)*4+6) = 2.d0*x1*x2*(-x1 + x2)
  rop((i-1)*4+7) = 2.d0*x1*x2*(-x1**2 + 3.d0*x1*x2 - x2**2)
  rop((i-1)*4+8) = 2.d0*x1*x2*(-x1**3 + 6.d0*(x1**2)*x2 - 6.d0*x1*(x2**2) + x2**3)
END DO
DO i=1,4
  x1 = f(oriented_faces(1,i)); x2 = f(oriented_faces(2,i)); x3 = f(oriented_faces(3,i))
  rop((i-1)*6+29) = -2.d0*x1*x2*x3
  rop((i-1)*6+30) = 2.d0*x1*x2*x3*(x1 + x2 - x3)
  rop((i-1)*6+31) = 2.d0*x1*x2*x3*(-x1 + x2)
  rop((i-1)*6+32) = 2.d0*x1*x2*x3*(-x1**2 - 2.d0*x1*x2 + 4.d0*x1*x3 - x2**2 + 4.d0*x2*x3 - x3**2)
  rop((i-1)*6+33) = 2.d0*x1*x2*x3*(x1**2 - x1*x3 - x2**2 + x2*x3)
  rop((i-1)*6+34) = 2.d0*x1*x2*x3*(-x1**2 + 3.d0*x1*x2 - x2**2)
END DO
x1 = f(1); x2 = f(2); x3 = f(3); x4 = f(4)
rop(53) = -2.d0*x1*x2*x3*x4
rop(54) = 2.d0*x1*x2*x3*x4*(x1 + x2 + x3 - x4)
rop(55) = 2.d0*x1*x2*x3*x4*(x1 + x2 - x3)
rop(56) = 2.d0*x1*x2*x3*x4*(-x1 + x2)
DEBUG_STACK_POP
end subroutine oft_h1_eval_all5
end subroutine oft_h1_eval_all
!---------------------------------------------------------------------------
!> Evaluate H^1 gradient function
!!
!! @note Evaluation is performed in logical coordinates with the resulting
!! vector in, and gradient with respect to, physical coordinates
!---------------------------------------------------------------------------
subroutine oft_h1_geval(self,cell,dof,f,val,gop)
class(oft_h1_fem), intent(in) :: self !< H^1 type for evaluation
integer(i4), intent(in) :: cell !< Cell for evaluation
integer(i4), intent(in) :: dof !< Element to evaluate
real(r8), intent(in) :: f(:) !< Position in cell in logical space
real(r8), intent(out) :: val(:) !< Gradient of H^1 element (dof) at point (f) [3]
real(r8), intent(in) :: gop(:,:) !< Cell Jacobian matrix at point (f) [3,4]
integer(i4) :: ed,etmp(2),fc,ftmp(3),i,j,ind,finds(9),fhtmp(4)
real(r8) :: cofs(4),fhex(6),gbary(3,6),dtmp,vtmp(4)
DEBUG_STACK_PUSH
IF(self%mesh%type==3)THEN
  val=0.d0
  fhex=hex_get_bary(f)
  gbary=hex_get_bary_gop(gop)
  select case(self%cmap(dof)%type)
    case(1)
      DO i=1,3
        dtmp=1.d0
        DO j=1,3
          IF(i==j)CYCLE
          dtmp=dtmp*fhex(hex_bary_pfcoords(j,self%cmap(dof)%el))
        END DO
        val = val + gbary(:,hex_bary_pfcoords(i,self%cmap(dof)%el))*dtmp
      END DO
    case(2)
      etmp=hex_bary_ecoords(:,self%cmap(dof)%el)
      call orient_list2(self%mesh%lce(self%cmap(dof)%el,cell), etmp)
      !
      dtmp=dhpoly_bary(self%cmap(dof)%ind, fhex(etmp(1)))*fhex(etmp(2))
      ! Fixed is always linear
      DO j=1,2
        dtmp=dtmp*fhex(hex_bary_efcoords(j,self%cmap(dof)%el))
      END DO
      val = val + gbary(:,etmp(1))*dtmp
      !
      dtmp=hpoly_bary(self%cmap(dof)%ind, fhex(etmp(1)))
      ! Fixed is always linear
      DO j=1,2
        dtmp=dtmp*fhex(hex_bary_efcoords(j,self%cmap(dof)%el))
      END DO
      val = val + gbary(:,etmp(2))*dtmp
      !
      dtmp=dtmp*fhex(etmp(2))
      val = val + gbary(:,hex_bary_efcoords(1,self%cmap(dof)%el))*dtmp* &
        fhex(hex_bary_efcoords(2,self%cmap(dof)%el))
      val = val + gbary(:,hex_bary_efcoords(2,self%cmap(dof)%el))*dtmp* &
        fhex(hex_bary_efcoords(1,self%cmap(dof)%el))
    case(3)
      fhtmp=hex_bary_fcoords(:,self%cmap(dof)%el)
      IF(self%mesh%lcfo(self%cmap(dof)%el,cell)<0)fhtmp=fhtmp((/2,3,4,1/))
      call orient_listn(self%mesh%lcfo(self%cmap(dof)%el,cell), fhtmp, 4_i4)
      !
      vtmp(1:3)=dhpoly_2d(self%indsf(:,self%cmap(dof)%ind), fhex(fhtmp(1:2)))
      DO i=1,2
        dtmp=vtmp(i+1)*fhex(self%cmap(dof)%el)*fhex(fhtmp(3))*fhex(fhtmp(4))
        val = val + gbary(:,fhtmp(i))*dtmp
      END DO
      !
      val = val + gbary(:,fhtmp(3))*vtmp(1)*fhex(self%cmap(dof)%el)*fhex(fhtmp(4))
      val = val + gbary(:,fhtmp(4))*vtmp(1)*fhex(self%cmap(dof)%el)*fhex(fhtmp(3))
      ! Fixed is always linear
      val = val + gbary(:,self%cmap(dof)%el)*vtmp(1)*fhex(fhtmp(3))*fhex(fhtmp(4))
    case(4)
      vtmp=dhpoly_3d(self%indsc(:,self%cmap(dof)%ind), fhex(1:3))
      DO i=1,3
        val=val+gbary(:,i)*vtmp(i+1)*fhex(4)*fhex(5)*fhex(6)
      END DO
      val=val+gbary(:,4)*vtmp(1)*fhex(5)*fhex(6)
      val=val+gbary(:,5)*vtmp(1)*fhex(4)*fhex(6)
      val=val+gbary(:,6)*vtmp(1)*fhex(4)*fhex(5)
  end select
ELSE
  IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
  cofs=0.d0
  select case(self%cmap(dof)%type)
    case(1)
      cofs(self%cmap(dof)%el)=1.d0
    case(2)
      etmp=oriented_edges(:,self%cmap(dof)%el)
      call oft_h1_gevale(self%order,etmp,self%cmap(dof)%ind,f,cofs)
    case(3)
      ftmp=oriented_faces(:,self%cmap(dof)%el)
      call oft_h1_gevalf(self%order,ftmp,self%cmap(dof)%ind,f,cofs)
    case(4)
      call oft_h1_gevalc(self%order,self%cmap(dof)%ind,f,cofs)
  end select
  !---Sum contributions
  val=0.d0
  do i=1,4
    val=val+gop(:,i)*cofs(i)
  end do
END IF
DEBUG_STACK_POP
end subroutine oft_h1_geval
!---------------------------------------------------------------------------
!> Evaluate H^1 gradient function on the boundary
!!
!! @note Evaluation is performed in logical coordinates with the resulting
!! vector in, and gradient with respect to, physical coordinates
!---------------------------------------------------------------------------
subroutine oft_bh1_geval(self,face,dof,f,val,gop)
class(oft_bfem_type), intent(in) :: self !< H^1 type for evaluation (bfem)
integer(i4), intent(in) :: face !< Face for evaluation
integer(i4), intent(in) :: dof !< Element to evaluate
real(r8), intent(in) :: f(:) !< Position on face in logical space [4]
real(r8), intent(out) :: val(3) !< Gradient of H^1 element (dof) at point (f) [3]
real(r8), optional, intent(in) :: gop(:,:) !< Face Jacobian matrix at point (f) [3,3]
real(r8) :: grads(3,4),cofs(4)
integer(i4) :: ed,etmp(2),fc,ftmp(3),i
DEBUG_STACK_PUSH
grads=1.d0
if(present(gop))grads(:,1:3)=gop
val=0.d0; cofs=0.d0
select case(self%cmap(dof)%type)
  case(1)
    cofs(self%cmap(dof)%el)=1.d0
  case(2)
    ed=self%mesh%lce(self%cmap(dof)%el,face)
    etmp=self%mesh%cell_ed(:,self%cmap(dof)%el)
    call orient_list2(ed,etmp)
    call oft_h1_gevale(self%order,etmp,self%cmap(dof)%ind,f,cofs)
  case(3)
    fc=self%mesh%lco(face)
    ftmp=(/1,2,3/)
    call orient_listn(fc,ftmp,3_i4)
    call oft_h1_gevalf(self%order,ftmp,self%cmap(dof)%ind,f,cofs)
end select
!---Sum contributions
do i=1,4
  val=val+grads(:,i)*cofs(i)
end do
DEBUG_STACK_POP
end subroutine oft_bh1_geval
!---------------------------------------------------------------------------
!> Evaluate point based gradient functions
!---------------------------------------------------------------------------
subroutine oft_h1_gevalp(order,pt,f,val)
integer(i4), intent(in) :: order !< Needs docs
integer(i4), intent(in) :: pt !< Needs docs
real(r8), intent(in) :: f(:) !< Needs docs
real(r8), intent(out) :: val(4) !< Needs docs
DEBUG_STACK_PUSH
val(pt)=1.d0
DEBUG_STACK_POP
end subroutine oft_h1_gevalp
!---------------------------------------------------------------------------
!> Evaluate edge based curl functions
!---------------------------------------------------------------------------
subroutine oft_h1_gevale(order,ed,dof,f,val)
integer(i4), intent(in) :: order !< Needs docs
integer(i4), intent(in) :: ed(2) !< Needs docs
integer(i4), intent(in) :: dof !< Needs docs
real(r8), intent(in) :: f(4) !< Needs docs
real(r8), intent(out) :: val(4) !< Needs docs
real(r8) :: x1,x2,y1,y2
DEBUG_STACK_PUSH
x1=f(ed(1)); x2=f(ed(2))
select case(dof)
  case(1)
    y1=-2.d0*x2
    y2=-2.d0*x1
  case(2)
    y1=2.d0*x2*(-2.d0*x1 + x2)
    y2=2.d0*x1*(-x1 + 2.d0*x2)
  case(3)
    y1=2.d0*x2*(-3.d0*x1**2 + 6.d0*x1*x2 - x2**2)
    y2=2.d0*x1*(-x1**2 + 6.d0*x1*x2 - 3.d0*x2**2)
  case(4)
    y1=2.d0*x2*(-4.d0*x1**3 + 18.d0*(x1**2)*x2 - 12.d0*x1*(x2**2) + x2**3)
    y2=2.d0*x1*(-x1**3 + 12.d0*(x1**2)*x2 - 18.d0*x1*(x2**2) + 4.d0*x2**3)
  case default
    call oft_abort('Invalid DOF!','oft_h1_gevale',__FILE__)
end select
val(ed(1))=y1; val(ed(2))=y2
DEBUG_STACK_POP
end subroutine oft_h1_gevale
!---------------------------------------------------------------------------
!> Evaluate face based curl functions
!---------------------------------------------------------------------------
subroutine oft_h1_gevalf(order,fc,dof,f,val)
integer(i4), intent(in) :: order !< Needs docs
integer(i4), intent(in) :: fc(3) !< Needs docs
integer(i4), intent(in) :: dof !< Needs docs
real(r8), intent(in) :: f(:) !< Needs docs
real(r8), intent(out) :: val(4) !< Needs docs
real(r8) :: x1,x2,x3,y1,y2,y3
DEBUG_STACK_PUSH
x1=f(fc(1)); x2=f(fc(2)); x3=f(fc(3))
select case(dof)
  case(1)
    y1=-2.d0*x2*x3
    y2=-2.d0*x1*x3
    y3=-2.d0*x1*x2
  case(2)
    y1=2.d0*x2*x3*(2.d0*x1 + x2 - x3)
    y2=2.d0*x1*x3*(x1 + 2.d0*x2 - x3)
    y3=2.d0*x1*x2*(x1 + x2 - 2.d0*x3)
  case(3)
    y1=2.d0*x2*x3*(-2.d0*x1 + x2)
    y2=2.d0*x1*x3*(-x1 + 2.d0*x2)
    y3=2.d0*x1*x2*(-x1 + x2)
  case(4)
    y1=2.d0*x2*x3*(-3.d0*x1**2 - 4.d0*x1*x2 + 8.d0*x1*x3 - x2**2 + 4.d0*x2*x3 - x3**2)
    y2=2.d0*x1*x3*(-x1**2 - 4.d0*x1*x2 + 4.d0*x1*x3 - 3.d0*x2**2 + 8.d0*x2*x3 - x3**2)
    y3=2.d0*x1*x2*(-x1**2 - 2.d0*x1*x2 + 8.d0*x1*x3 - x2**2 + 8.d0*x2*x3 - 3.d0*x3**2)
  case(5)
    y1=2.d0*x2*x3*(3.d0*x1**2 - 2.d0*x1*x3 - x2**2 + x2*x3)
    y2=2.d0*x1*x3*(x1**2 - x1*x3 - 3.d0*x2**2 + 2.d0*x2*x3)
    y3=2.d0*x1*x2*(x1**2 - 2.d0*x1*x3 - x2**2 + 2.d0*x2*x3)
  case(6)
    y1=2.d0*x2*x3*(-3.d0*x1**2 + 6.d0*x1*x2 - x2**2)
    y2=2.d0*x1*x3*(-x1**2 + 6.d0*x1*x2 - 3.d0*x2**2)
    y3=2.d0*x1*x2*(-x1**2 + 3.d0*x1*x2 - x2**2)
  case default
    call oft_abort('Invalid DOF!','oft_h1_gevalf',__FILE__)
end select
val(fc(1))=y1; val(fc(2))=y2; val(fc(3))=y3
DEBUG_STACK_POP
end subroutine oft_h1_gevalf
!---------------------------------------------------------------------------
!> Evaluate cell based curl functions
!---------------------------------------------------------------------------
subroutine oft_h1_gevalc(order,dof,f,val)
integer(i4), intent(in) :: order !< Needs docs
integer(i4), intent(in) :: dof !< Needs docs
real(r8), intent(in) :: f(:) !< Needs docs
real(r8), intent(out) :: val(4) !< Needs docs
real(r8) :: x1,x2,x3,x4
DEBUG_STACK_PUSH
x1=f(1); x2=f(2); x3=f(3); x4=f(4)
select case(dof)
  case(1)
    val(1)=-2.d0*x2*x3*x4
    val(2)=-2.d0*x1*x3*x4
    val(3)=-2.d0*x1*x2*x4
    val(4)=-2.d0*x1*x2*x3
  case(2)
    val(1)=2.d0*x2*x3*x4*(2.d0*x1 + x2 + x3 - x4)
    val(2)=2.d0*x1*x3*x4*(x1 + 2.d0*x2 + x3 - x4)
    val(3)=2.d0*x1*x2*x4*(x1 + x2 + 2.d0*x3 - x4)
    val(4)=2.d0*x1*x2*x3*(x1 + x2 + x3 - 2.d0*x4)
  case(3)
    val(1)=2.d0*x2*x3*x4*(2.d0*x1 + x2 - x3)
    val(2)=2.d0*x1*x3*x4*(x1 + 2.d0*x2 - x3)
    val(3)=2.d0*x1*x2*x4*(x1 + x2 - 2.d0*x3)
    val(4)=2.d0*x1*x2*x3*(x1 + x2 - x3)
  case(4)
    val(1)=2.d0*x2*x3*x4*(-2.d0*x1 + x2)
    val(2)=2.d0*x1*x3*x4*(-x1 + 2.d0*x2)
    val(3)=2.d0*x1*x2*x4*(-x1 + x2)
    val(4)=2.d0*x1*x2*x3*(-x1 + x2)
  case default
    call oft_abort('Invalid DOF!','oft_h1_gevalc',__FILE__)
end select
DEBUG_STACK_POP
end subroutine oft_h1_gevalc
!---------------------------------------------------------------------------
!> Evaluate all H^1 interpolation functions
!!
!! @note Evaluation is performed in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_h1_geval_all(self,cell,f,rop,gop)
class(oft_h1_fem), intent(in) :: self !< H^1 type for evaluation
integer(i4), intent(in) :: cell !< Cell for evaluation
real(r8), intent(in) :: f(:) !< Position in cell in logical space
real(r8), contiguous, intent(out) :: rop(:,:) !< Value of interpolation functions at point (f) [3,ncdofs]
real(r8), intent(in) :: gop(:,:) !< Cell Jacobian matrix at point (f) [3,4]
integer(i4) :: i,j,k,etmp(2),fhtmp(4),offset
real(r8) :: fhex(6),gbary(3,6),val(3),cords(4),dtmp,vtmp(4)
DEBUG_STACK_PUSH
IF(self%mesh%type==3)THEN
  fhex=hex_get_bary(f)
  gbary=hex_get_bary_gop(gop)
  DO i=1,8
    cords(1:3)=fhex(hex_bary_pfcoords(:,i))
    rop(:,i)=gbary(:,hex_bary_pfcoords(1,i))*cords(2)*cords(3) &
            +gbary(:,hex_bary_pfcoords(2,i))*cords(1)*cords(3) &
            +gbary(:,hex_bary_pfcoords(3,i))*cords(1)*cords(2)
  END DO
  IF(self%gstruct(2)>0)THEN
    !---Edges
    DO i=1,12
      offset=(i-1)*self%gstruct(2)+8
      etmp=hex_bary_ecoords(:,i)
      call orient_list2(self%mesh%lce(i,cell), etmp)
      cords(1:2)=fhex(etmp)
      cords(3:4)=fhex(hex_bary_efcoords(:,i))
      DO j=1,self%gstruct(2)
        val = gbary(:,etmp(1))*dhpoly_bary(j, cords(1))*cords(2)*cords(3)*cords(4)
        !
        dtmp=hpoly_bary(j, cords(1))
        val = val + gbary(:,etmp(2))*dtmp*cords(3)*cords(4)
        !
        dtmp=dtmp*cords(2)
        val = val + gbary(:,hex_bary_efcoords(1,i))*dtmp*cords(4)
        val = val + gbary(:,hex_bary_efcoords(2,i))*dtmp*cords(3)
        !
        rop(:,offset+j)=val
      END DO
    END DO
  END IF
  IF(self%gstruct(3)>0)THEN
    !---Faces
    DO i=1,6
      offset=(i-1)*self%gstruct(3)+12*self%gstruct(2)+8
      fhtmp=hex_bary_fcoords(:,i)
      IF(self%mesh%lcfo(i,cell)<0)fhtmp=fhtmp((/2,3,4,1/))
      call orient_listn(self%mesh%lcfo(i,cell), fhtmp, 4_i4)
      cords=fhex(fhtmp)
      DO j=1,self%gstruct(3)
        val=0.d0
        vtmp(1:3)=dhpoly_2d(self%indsf(1:2,j), cords(1:2))
        DO k=1,2
          dtmp=vtmp(k+1)*fhex(i)*cords(3)*cords(4)
          val = val + gbary(:,fhtmp(k))*dtmp
        END DO
        !
        val = val + gbary(:,fhtmp(3))*vtmp(1)*fhex(i)*cords(4)
        val = val + gbary(:,fhtmp(4))*vtmp(1)*fhex(i)*cords(3)
        ! Fixed is always linear
        val = val + gbary(:,i)*vtmp(1)*cords(3)*cords(4)
        !
        rop(:,offset+j)=val
      END DO
    END DO
  END IF
  IF(self%gstruct(4)>0)THEN
    !---Cell
    offset=6*self%gstruct(3)+12*self%gstruct(2)+8
    DO j=1,self%gstruct(4)
      val=0.d0
      vtmp=dhpoly_3d(self%indsc(1:3,j), fhex(1:3))
      DO k=1,3
        val=val+gbary(:,k)*vtmp(k+1)*fhex(4)*fhex(5)*fhex(6)
      END DO
      val=val+gbary(:,4)*vtmp(1)*fhex(5)*fhex(6)
      val=val+gbary(:,5)*vtmp(1)*fhex(4)*fhex(6)
      val=val+gbary(:,6)*vtmp(1)*fhex(4)*fhex(5)
      !
      rop(:,offset+j)=val
    END DO
  END IF
ELSE
  IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
  select case(self%order)
    case(1)
      DO i=1,4
        rop(:,i)=gop(:,i)
      END DO
    case(2)
      call oft_h1_geval_all2()!self,cell,f,rop,gop)
    case(3)
      call oft_h1_geval_all3()!self,cell,f,rop,gop)
    case(4)
      call oft_h1_geval_all4()!self,cell,f,rop,gop)
    case(5)
      call oft_h1_geval_all5()!self,cell,f,rop,gop)
    case default ! Fall back to normal evaluation
      DO i=1,self%nce
        call oft_h1_geval(self,cell,i,f,rop(:,i),gop)
      END DO
  end select
END IF
DEBUG_STACK_POP
contains
!---------------------------------------------------------------------------
!> Evaluate all H^1 interpolation functions (quadratic)
!!
!! @note Evaluation is performed in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_h1_geval_all2()!self,cell,f,rop,gop)
! class(oft_fem_type), intent(in) :: self !< H^1 type for evaluation
! integer(i4), intent(in) :: cell !< Cell for evaluation
! real(r8), intent(in) :: f(4) !< Position in cell in logical space
! real(r8), intent(out) :: rop(3,10) !< Value of interpolation functions at point (f) [ncdofs]
! real(r8), intent(in) :: gop(3,4) !< Cell Jacobian matrix at point (f) [3,4]
integer(i4) :: i,etmp(2)
real(r8) :: x1,x2,u1(3),u2(3)
DEBUG_STACK_PUSH
DO i=1,4
  rop(:,i)=gop(:,i)
END DO
!---
DO i=1,6
  etmp=self%mesh%cell_ed(:,i)
  x1 = f(etmp(1)); x2 = f(etmp(2))
  u1 = gop(:,etmp(1)); u2 = gop(:,etmp(2))
  !---
  rop(:,i+4) = -2.d0*x2*u1 &
               -2.d0*x1*u2
END DO
DEBUG_STACK_POP
end subroutine oft_h1_geval_all2
!---------------------------------------------------------------------------
!> Evaluate all H^1 interpolation functions (cubic)
!!
!! @note Evaluation is performed in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_h1_geval_all3()!self,cell,f,rop,gop)
! class(oft_fem_type), intent(in) :: self !< H^1 type for evaluation
! integer(i4), intent(in) :: cell !< Cell for evaluation
! real(r8), intent(in) :: f(4) !< Position in cell in logical space
! real(r8), intent(out) :: rop(3,20) !< Value of interpolation functions at point (f) [ncdofs]
! real(r8), intent(in) :: gop(3,4) !< Cell Jacobian matrix at point (f) [3,4]
integer(i4) :: i,etmp(2),ftmp(3)
real(r8) :: x1,x2,x3,u1(3),u2(3),u3(3)
DEBUG_STACK_PUSH
DO i=1,4
  rop(:,i)=gop(:,i)
END DO
DO i=1,6
  etmp=oriented_edges(:,i)
  x1 = f(etmp(1)); x2 = f(etmp(2))
  u1 = gop(:,etmp(1)); u2 = gop(:,etmp(2))
  !---
  rop(:,(i-1)*2+5) = -2.d0*x2*u1 &
               -2.d0*x1*u2
  !---
  rop(:,(i-1)*2+6) = 2.d0*x2*(-2.d0*x1 + x2)*u1 &
              + 2.d0*x1*(-x1 + 2.d0*x2)*u2
END DO
DO i=1,4
  ftmp=self%mesh%cell_fc(:,i)
  x1 = f(ftmp(1)); x2 = f(ftmp(2)); x3 = f(ftmp(3))
  u1 = gop(:,ftmp(1)); u2 = gop(:,ftmp(2)); u3 = gop(:,ftmp(3))
  !---
  rop(:,i+16) = -2.d0*x2*x3*u1 &
                -2.d0*x1*x3*u2 &
                -2.d0*x1*x2*u3
END DO
DEBUG_STACK_POP
end subroutine oft_h1_geval_all3
!---------------------------------------------------------------------------
!> Evaluate all H^1 interpolation functions (quartic)
!!
!! @note Evaluation is performed in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_h1_geval_all4()!self,cell,f,rop,gop)
! class(oft_fem_type), intent(in) :: self !< H^1 type for evaluation
! integer(i4), intent(in) :: cell !< Cell for evaluation
! real(r8), intent(in) :: f(4) !< Position in cell in logical space
! real(r8), intent(out) :: rop(3,35) !< Value of interpolation functions at point (f) [ncdofs]
! real(r8), intent(in) :: gop(3,4) !< Cell Jacobian matrix at point (f) [3,4]
integer(i4) :: i,etmp(2),ftmp(3)
real(r8) :: x1,x2,x3,x4,u1(3),u2(3),u3(3)
DEBUG_STACK_PUSH
DO i=1,4
  rop(:,i)=gop(:,i)
END DO
DO i=1,6
  etmp=oriented_edges(:,i)
  x1 = f(etmp(1)); x2 = f(etmp(2))
  u1 = gop(:,etmp(1)); u2 = gop(:,etmp(2))
  !---
  rop(:,(i-1)*3+5) = -2.d0*x2*u1 &
               -2.d0*x1*u2
  !---
  rop(:,(i-1)*3+6) = 2.d0*x2*(-2.d0*x1 + x2)*u1 &
              + 2.d0*x1*(-x1 + 2.d0*x2)*u2
  !---
  rop(:,(i-1)*3+7) = 2.d0*x2*(-3.d0*x1**2 + 6.d0*x1*x2 - x2**2)*u1 &
              + 2.d0*x1*(-x1**2 + 6.d0*x1*x2 - 3.d0*x2**2)*u2
END DO
DO i=1,4
  ftmp=oriented_faces(:,i)
  x1 = f(ftmp(1)); x2 = f(ftmp(2)); x3 = f(ftmp(3))
  u1 = gop(:,ftmp(1)); u2 = gop(:,ftmp(2)); u3 = gop(:,ftmp(3))
  !---
  rop(:,(i-1)*3+23) = -2.d0*x2*x3*u1 &
                -2.d0*x1*x3*u2 &
                -2.d0*x1*x2*u3
  !---
  rop(:,(i-1)*3+24) = 2.d0*x2*x3*(2.d0*x1 + x2 - x3)*u1 &
              + 2.d0*x1*x3*(x1 + 2.d0*x2 - x3)*u2 &
              + 2.d0*x1*x2*(x1 + x2 - 2.d0*x3)*u3
  !---
  rop(:,(i-1)*3+25) = 2.d0*x2*x3*(-2.d0*x1 + x2)*u1 &
              + 2.d0*x1*x3*(-x1 + 2.d0*x2)*u2 &
              + 2.d0*x1*x2*(-x1 + x2)*u3
END DO
x1 = f(1); x2 = f(2); x3 = f(3); x4 = f(4)
!
rop(:,35) = -2.d0*x2*x3*x4*gop(:,1) &
            -2.d0*x1*x3*x4*gop(:,2) &
            -2.d0*x1*x2*x4*gop(:,3) &
            -2.d0*x1*x2*x3*gop(:,4)
DEBUG_STACK_POP
end subroutine oft_h1_geval_all4
!---------------------------------------------------------------------------
!> Evaluate all H^1 interpolation functions (quartic)
!!
!! @note Evaluation is performed in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_h1_geval_all5()!self,cell,f,rop,gop)
! class(oft_fem_type), intent(in) :: self !< H^1 type for evaluation
! integer(i4), intent(in) :: cell !< Cell for evaluation
! real(r8), intent(in) :: f(4) !< Position in cell in logical space
! real(r8), intent(out) :: rop(3,56) !< Value of interpolation functions at point (f) [ncdofs]
! real(r8), intent(in) :: gop(3,4) !< Cell Jacobian matrix at point (f) [3,4]
integer(i4) :: i,etmp(2),ftmp(3)
real(r8) :: x1,x2,x3,x4,u1(3),u2(3),u3(3)
DEBUG_STACK_PUSH
DO i=1,4
  rop(:,i)=gop(:,i)
END DO
DO i=1,6
  etmp=oriented_edges(:,i)
  x1 = f(etmp(1)); x2 = f(etmp(2))
  u1 = gop(:,etmp(1)); u2 = gop(:,etmp(2))
  !---
  rop(:,(i-1)*4+5) = -2.d0*x2*u1 &
               -2.d0*x1*u2
  !---
  rop(:,(i-1)*4+6) = 2.d0*x2*(-2.d0*x1 + x2)*u1 &
              + 2.d0*x1*(-x1 + 2.d0*x2)*u2
  !---
  rop(:,(i-1)*4+7) = 2.d0*x2*(-3.d0*x1**2 + 6.d0*x1*x2 - x2**2)*u1 &
              + 2.d0*x1*(-x1**2 + 6.d0*x1*x2 - 3.d0*x2**2)*u2
  !---
  rop(:,(i-1)*4+8) = 2.d0*x2*(-4.d0*x1**3 + 18.d0*(x1**2)*x2 - 12.d0*x1*(x2**2) + x2**3)*u1 &
              + 2.d0*x1*(-x1**3 + 12.d0*(x1**2)*x2 - 18.d0*x1*(x2**2) + 4.d0*x2**3)*u2
END DO
DO i=1,4
  ftmp=oriented_faces(:,i)
  x1 = f(ftmp(1)); x2 = f(ftmp(2)); x3 = f(ftmp(3))
  u1 = gop(:,ftmp(1)); u2 = gop(:,ftmp(2)); u3 = gop(:,ftmp(3))
  !---
  rop(:,(i-1)*6+29) = -2.d0*x2*x3*u1 &
                -2.d0*x1*x3*u2 &
                -2.d0*x1*x2*u3
  !---
  rop(:,(i-1)*6+30) = 2.d0*x2*x3*(2.d0*x1 + x2 - x3)*u1 &
              + 2.d0*x1*x3*(x1 + 2.d0*x2 - x3)*u2 &
              + 2.d0*x1*x2*(x1 + x2 - 2.d0*x3)*u3
  !---
  rop(:,(i-1)*6+31) = 2.d0*x2*x3*(-2.d0*x1 + x2)*u1 &
              + 2.d0*x1*x3*(-x1 + 2.d0*x2)*u2 &
              + 2.d0*x1*x2*(-x1 + x2)*u3
  !---
  rop(:,(i-1)*6+32) = 2.d0*x2*x3*(-3.d0*x1**2 - 4.d0*x1*x2 + 8.d0*x1*x3 - x2**2 + 4.d0*x2*x3 - x3**2)*u1 &
              + 2.d0*x1*x3*(-x1**2 - 4.d0*x1*x2 + 4.d0*x1*x3 - 3.d0*x2**2 + 8.d0*x2*x3 - x3**2)*u2 &
              + 2.d0*x1*x2*(-x1**2 - 2.d0*x1*x2 + 8.d0*x1*x3 - x2**2 + 8.d0*x2*x3 - 3.d0*x3**2)*u3
  !---
  rop(:,(i-1)*6+33) = 2.d0*x2*x3*(3.d0*x1**2 - 2.d0*x1*x3 - x2**2 + x2*x3)*u1 &
              + 2.d0*x1*x3*(x1**2 - x1*x3 - 3.d0*x2**2 + 2.d0*x2*x3)*u2 &
              + 2.d0*x1*x2*(x1**2 - 2.d0*x1*x3 - x2**2 + 2.d0*x2*x3)*u3
  !---
  rop(:,(i-1)*6+34) = 2.d0*x2*x3*(-3.d0*x1**2 + 6.d0*x1*x2 - x2**2)*u1 &
              + 2.d0*x1*x3*(-x1**2 + 6.d0*x1*x2 - 3.d0*x2**2)*u2 &
              + 2.d0*x1*x2*(-x1**2 + 3.d0*x1*x2 - x2**2)*u3
END DO
x1 = f(1); x2 = f(2); x3 = f(3); x4 = f(4)
!
rop(:,53) = -2.d0*x2*x3*x4*gop(:,1) &
            -2.d0*x1*x3*x4*gop(:,2) &
            -2.d0*x1*x2*x4*gop(:,3) &
            -2.d0*x1*x2*x3*gop(:,4)
!
rop(:,54) = 2.d0*x2*x3*x4*(2.d0*x1 + x2 + x3 - x4)*gop(:,1) &
          + 2.d0*x1*x3*x4*(x1 + 2.d0*x2 + x3 - x4)*gop(:,2) &
          + 2.d0*x1*x2*x4*(x1 + x2 + 2.d0*x3 - x4)*gop(:,3) &
          + 2.d0*x1*x2*x3*(x1 + x2 + x3 - 2.d0*x4)*gop(:,4)
!
rop(:,55) = 2.d0*x2*x3*x4*(2.d0*x1 + x2 - x3)*gop(:,1) &
          + 2.d0*x1*x3*x4*(x1 + 2.d0*x2 - x3)*gop(:,2) &
          + 2.d0*x1*x2*x4*(x1 + x2 - 2.d0*x3)*gop(:,3) &
          + 2.d0*x1*x2*x3*(x1 + x2 - x3)*gop(:,4)
!
rop(:,56) = 2.d0*x2*x3*x4*(-2.d0*x1 + x2)*gop(:,1) &
          + 2.d0*x1*x3*x4*(-x1 + 2.d0*x2)*gop(:,2) &
          + 2.d0*x1*x2*x4*(-x1 + x2)*gop(:,3) &
          + 2.d0*x1*x2*x3*(-x1 + x2)*gop(:,4)
DEBUG_STACK_POP
end subroutine oft_h1_geval_all5
end subroutine oft_h1_geval_all
!---------------------------------------------------------------------------
!> Evaluate H^1 gradient function
!!
!! @note Evaluation is performed in logical coordinates with the resulting
!! gradient with respect to physical coordinates
!---------------------------------------------------------------------------
subroutine oft_h1_d2eval(self,cell,dof,f,val,g2op)
class(oft_fem_type), intent(in) :: self !< H^1 type for evaluation
integer(i4), intent(in) :: cell !< Cell for evaluation
integer(i4), intent(in) :: dof !< Element to evaluate
real(r8), intent(in) :: f(:) !< Position in cell in logical space
real(r8), intent(out) :: val(6) !< Gradient of H^1 element (dof) at point (f) [3]
real(r8), intent(in) :: g2op(6,10) !< Cell Jacobian matrix at point (f) [3,4]
real(r8) :: grad(2),cofs(10)
integer(i4) :: ed,etmp(2),fc,ftmp(3),i
IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
val=0.d0; cofs=0.d0
select case(self%cmap(dof)%type)
  case(1)
    cofs=0.d0
  case(2)
    etmp=oriented_edges(:,self%cmap(dof)%el)
    call oft_h1_d2evale(self%order,etmp,self%cmap(dof)%el,self%cmap(dof)%ind,f,cofs)
  case(3)
    ftmp=oriented_faces(:,self%cmap(dof)%el)
    call oft_h1_d2evalf(self%order,ftmp,self%cmap(dof)%el,self%cmap(dof)%ind,f,cofs)
  case(4)
    call oft_h1_d2evalc(self%order,self%cmap(dof)%ind,f,cofs)
  case default
    call oft_abort("Invalid dof type","oft_h1_d2eval",__FILE__)
end select
!---Sum contributions
do i=1,10
  val=val+g2op(:,i)*cofs(i)
end do
end subroutine oft_h1_d2eval
!---------------------------------------------------------------------------
!> @brief Evaluate edge based gradient functions
!---------------------------------------------------------------------------
subroutine oft_h1_d2evale(order,ed,el,dof,f,val)
integer(i4), intent(in) :: order !< Needs docs
integer(i4), intent(in) :: ed(2) !< Needs docs
integer(i4), intent(in) :: el !< Needs docs
integer(i4), intent(in) :: dof !< Needs docs
real(r8), intent(in) :: f(:)
real(r8), intent(out) :: val(10)
real(r8) :: x1,x2
integer(i4), parameter :: ptmp(4)=(/1,5,8,10/)
integer(i4), parameter :: etmp(6)=(/4,7,9,6,3,2/)
val=0.d0
x1=f(ed(1)); x2=f(ed(2))
select case(dof)
  case(1)
    val(etmp(el))=-2.d0
  case(2)
    val(ptmp(ed(1)))=-4.d0*x2
    val(etmp(el))   =4.d0*(x2 - x1)
    val(ptmp(ed(2)))=4.d0*x1
  case(3)
    val(ptmp(ed(1)))=12.d0*x2*(x2 - x1)
    val(etmp(el))   =-6.d0*(x2**2 + x1**2) + 24.d0*x1*x2
    val(ptmp(ed(2)))=12.d0*x1*(x1 - x2)
  case(4)
    val(ptmp(ed(1)))=-24.d0*x2*(x2**2 - 3.d0*x1*x2 + x1**2)
    val(etmp(el))   =8.d0*(x2**3 - x1**3) + 72.d0*x2*x1**2 - 72.d0*x1*x2**2
    val(ptmp(ed(2)))=24.d0*x2*(x2**2 - 3.d0*x1*x2 + x1**2)
  case default
    call oft_abort("Invalid dof","oft_h1_d2evale",__FILE__)
end select
end subroutine oft_h1_d2evale
!---------------------------------------------------------------------------
!> @brief Evaluate cell based gradient functions
!---------------------------------------------------------------------------
subroutine oft_h1_d2evalf(order,fc,el,dof,f,val)
integer(i4), intent(in) :: order !< Needs docs
integer(i4), intent(in) :: fc(3) !< Needs docs
integer(i4), intent(in) :: el !< Needs docs
integer(i4), intent(in) :: dof !< Needs docs
real(r8), intent(in) :: f(:)
real(r8), intent(out) :: val(10)
real(r8) :: x1,x2,x3
integer(i4) :: dmap(6)
integer(i4), parameter :: ftmp(4,4)=RESHAPE((/1,2,3,4,2,5,6,7,3,6,8,9,4,7,9,10/),(/4,4/))
val=0.d0
dmap=(/ftmp(fc(1),fc(1)),ftmp(fc(1),fc(2)),ftmp(fc(1),fc(3)),ftmp(fc(2),fc(2)),ftmp(fc(2),fc(3)), &
ftmp(fc(3),fc(3))/)
x1=f(fc(1)); x2=f(fc(2)); x3=f(fc(3))
select case(dof)
  case(1)
    val(dmap(2))=-2.d0*x3
    val(dmap(3))=-2.d0*x2
    val(dmap(5))=-2.d0*x1
  case(2)
    val(dmap(1))=4.d0*x2*x3
    val(dmap(2))=2.d0*x3*(2.d0*(x1 + x2) - x3)
    val(dmap(3))=2.d0*x2*(2.d0*(x1 - x3) + x2)
    val(dmap(4))=4.d0*x1*x3
    val(dmap(5))=2.d0*x1*(2.d0*(x2 - x3) + x1)
    val(dmap(6))=-4.d0*x1*x2
  case(3)
    val(dmap(1))=-4.d0*x2*x3
    val(dmap(2))=4.d0*x3*(-x1 + x2)
    val(dmap(3))=2.d0*x2*(-2.d0*x1 + x2)
    val(dmap(4))=4.d0*x1*x3
    val(dmap(5))=2.d0*x1*(-x1 + 2.d0*x2)
    val(dmap(6))=0.d0
  case(4)
    val(dmap(1))=4.d0*x2*x3*(-3.d0*x1 - 2.d0*x2 + 4.d0*x3)
    val(dmap(2))=2.d0*x3*(-3.d0*x1**2 - 8.d0*x1*x2 + 8.d0*x1*x3 - 3.d0*x2**2 + 8.d0*x2*x3 - x3**2)
    val(dmap(3))=2.d0*x2*(-3.d0*x1**2 - 4.d0*x1*x2 + 16.d0*x1*x3 - x2**2 + 8.d0*x2*x3 - 3.d0*x3**2)
    val(dmap(4))=4.d0*x1*x3*(-2.d0*x1 - 3.d0*x2 + 4.d0*x3)
    val(dmap(5))=2.d0*x1*(-x1**2 - 4.d0*x1*x2 + 8.d0*x1*x3 - 3.d0*x2**2 + 16.d0*x2*x3 - 3.d0*x3**2)
    val(dmap(6))=4.d0*x1*x2*(4.d0*x1 + 4.d0*x2 - 3.d0*x3)
  case(5)
    val(dmap(1))=4.d0*x2*x3*(3.d0*x1 - x3)
    val(dmap(2))=2.d0*x3*(3.d0*x1**2 - 2.d0*x1*x3 - 3.d0*x2**2 + 2.d0*x2*x3)
    val(dmap(3))=2.d0*x2*(3.d0*x1**2 - 4.d0*x1*x3 - x2**2 + 2.d0*x2*x3)
    val(dmap(4))=4.d0*x1*x3*(-3.d0*x2 + x3)
    val(dmap(5))=2.d0*x1*(x1**2 - 2.d0*x1*x3 - 3.d0*x2**2 + 4.d0*x2*x3)
    val(dmap(6))=4.d0*x1*x2*(-x1 + x2)
  case(6)
    val(dmap(1))=12.d0*x2*x3*(-x1 + x2)
    val(dmap(2))=6.d0*x3*(-x1**2 + 4.d0*x1*x2 - x2**2)
    val(dmap(3))=2.d0*x2*(-3.d0*x1**2 + 6.d0*x1*x2 - x2**2)
    val(dmap(4))=12.d0*x1*x3*(x1 - x2)
    val(dmap(5))=2.d0*x1*(-x1**2 + 6.d0*x1*x2 - 3.d0*x2**2)
    val(dmap(6))=0.d0
  case default
    call oft_abort("Invalid dof","oft_h1_d2evalf",__FILE__)
end select
end subroutine oft_h1_d2evalf
!---------------------------------------------------------------------------
!> @brief Evaluate cell based gradient functions
!---------------------------------------------------------------------------
subroutine oft_h1_d2evalc(order,dof,f,val)
integer(i4), intent(in) :: order !< Needs docs
integer(i4), intent(in) :: dof !< Needs docs
real(r8), intent(in) :: f(:) !< Needs docs
real(r8), intent(out) :: val(10) !< Needs docs
real(r8) :: x1,x2,x3,x4
val=0.d0
x1=f(1); x2=f(2); x3=f(3); x4=f(4)
select case(dof)
  case(1)
    val(1) =0.d0
    val(2) =-2.d0*f(3)*f(4)
    val(3) =-2.d0*f(2)*f(4)
    val(4) =-2.d0*f(2)*f(3)
    val(5) =0.d0
    val(6) =-2.d0*f(1)*f(4)
    val(7) =-2.d0*f(1)*f(3)
    val(8) =0.d0
    val(9) =-2.d0*f(1)*f(2)
    val(10)=0.d0
  case(2)
    val(1) =4.d0*x2*x3*x4
    val(2) =2.d0*x3*x4*(2.d0*x1 + 2.d0*x2 + x3 - x4)
    val(3) =2.d0*x2*x4*(2.d0*x1 + x2 + 2.d0*x3 - x4)
    val(4) =2.d0*x2*x3*(2.d0*x1 + x2 + x3 - 2.d0*x4)
    val(5) =4.d0*x1*x3*x4
    val(6) =2.d0*x1*x4*(x1 + 2.d0*x2 + 2.d0*x3 - x4)
    val(7) =2.d0*x1*x3*(x1 + 2.d0*x2 + x3 - 2.d0*x4)
    val(8) =4.d0*x1*x2*x4
    val(9) =2.d0*x1*x2*(x1 + x2 + 2.d0*x3 - 2.d0*x4)
    val(10)=-4.d0*x1*x2*x3
  case(3)
    val(1) =4.d0*x2*x3*x4
    val(2) =2.d0*x3*x4*(2.d0*x1 + 2.d0*x2 - x3)
    val(3) =2.d0*x2*x4*(2.d0*x1 + x2 - 2.d0*x3)
    val(4) =2.d0*x2*x3*(2.d0*x1 + x2 - x3)
    val(5) =4.d0*x1*x3*x4
    val(6) =2.d0*x1*x4*(x1 + 2.d0*x2 - 2.d0*x3)
    val(7) =2.d0*x1*x3*(x1 + 2.d0*x2 - x3)
    val(8) =-4.d0*x1*x2*x4
    val(9) =2.d0*x1*x2*(x1 + x2 - 2.d0*x3)
    val(10)=0.d0
  case(4)
    val(1) =-4.d0*x2*x3*x4
    val(2) =4.d0*x3*x4*(-x1 + x2)
    val(3) =2.d0*x2*x4*(-2.d0*x1 + x2)
    val(4) =2.d0*x2*x3*(-2.d0*x1 + x2)
    val(5) =4.d0*x1*x3*x4
    val(6) =2.d0*x1*x4*(-x1 + 2.d0*x2)
    val(7) =2.d0*x1*x3*(-x1 + 2.d0*x2)
    val(8) =0.d0
    val(9) =2.d0*x1*x2*(-x1 + x2)
    val(10)=0.d0
  case default
    call oft_abort("Invalid dof","oft_h1_d2evalc",__FILE__)
end select
end subroutine oft_h1_d2evalc
end module oft_h1_basis
