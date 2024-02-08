!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_hcurl_basis.F90
!
!> @defgroup doxy_oft_hcurl Nedelec H1(Curl)
!! Nedelec H1(Curl) finite element implementation for the Open FUSION Toolkit
!! @ingroup doxy_oft_fem
!
!> Base Nedelec H1(Curl) FE class and basis evaluation
!! - FE Construction
!! - Basis evaluation
!!   - Interpolation
!!   - Curl
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_hcurl
!---------------------------------------------------------------------------
MODULE oft_hcurl_basis
! USE timer
USE oft_base
USE oft_lag_poly
USE oft_mesh_type, ONLY: mesh, smesh
USE oft_mesh_local_util, ONLY: mesh_local_orient, oriented_cell, &
  oriented_edges, oriented_faces
USE oft_hexmesh_type, ONLY: hex_bary_ecoords, hex_bary_efcoords, hex_bary_fcoords, &
  hex_get_bary, hex_get_bary_gop, hex_get_bary_cgop
USE multigrid, ONLY: mg_mesh, multigrid_level
USE oft_la_base, ONLY: oft_matrix, oft_graph
USE fem_base, ONLY: oft_fem_type, oft_ml_fem_type, oft_bfem_type
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
type :: oft_nedelec_ops
  type(oft_graph), pointer :: interp_graph => NULL() !< Interpolation graph
  class(oft_matrix), pointer :: interp => NULL() !< Interpolation graph
end type oft_nedelec_ops
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
type, extends(oft_fem_type) :: oft_hcurl_fem
  INTEGER(i4), POINTER, DIMENSION(:,:) :: indsf => NULL()
  INTEGER(i4), POINTER, DIMENSION(:,:) :: indsc => NULL()
end type oft_hcurl_fem
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
type, extends(oft_bfem_type) :: oft_hcurl_bfem
  INTEGER(i4), POINTER, DIMENSION(:,:) :: indsf => NULL()
end type oft_hcurl_bfem
!---Global Variables
integer(i4), parameter :: cgop_map(4,4) = RESHAPE((/0,-1,-2,-3,1,0,-4,-5,2,4,0,-6,3,5,6,0/),(/4,4/))
integer(i4), parameter :: oft_hcurl_id = 3 !< FE type ID
integer(i4) :: oft_hcurl_blevel = 0 !< Highest level on base meshes
integer(i4) :: oft_hcurl_lev = 0 !< Active FE level
integer(i4) :: oft_hcurl_level = 0 !< Active FE level
integer(i4) :: oft_hcurl_lin_level = 0 !< Highest linear element level
integer(i4) :: oft_hcurl_minlev = 0 !<
integer(i4) :: oft_hcurl_nlevels = 0 !< Number of total levels
!
type(oft_hcurl_bfem), pointer :: oft_bhcurl !< Active FE representation
type(oft_hcurl_bfem), pointer :: oft_bhcurl_lin !< Highest linear element representation
type(oft_ml_fem_type), TARGET :: ML_oft_bhcurl !< ML container for all FE representations
!
class(oft_hcurl_fem), pointer :: oft_hcurl !< Active FE representation
class(oft_hcurl_fem), pointer :: oft_hcurl_lin !< Highest linear element representation
type(oft_ml_fem_type), TARGET :: ML_oft_hcurl !< ML container for all FE representations
!
type(oft_nedelec_ops), pointer :: oft_hcurl_ops !< Active operators
type(oft_nedelec_ops), pointer :: oft_hcurl_ops_lin !< Highest linear element operators
type(oft_nedelec_ops), pointer :: ML_oft_hcurl_ops(:) !< ML container for all operators
!
logical, private :: hex_mesh = .FALSE.
contains
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_set_level
!---------------------------------------------------------------------------
!> Set the current level for Nedelec H1(Curl) FE
!!
!! @param[in] level Desired level
!---------------------------------------------------------------------------
subroutine oft_hcurl_set_level(level)
integer(i4), intent(in) :: level
DEBUG_STACK_PUSH
if(level>oft_hcurl_nlevels.OR.level<=0)then
  write(*,*)level
  call oft_abort('Invalid FE level','oft_hcurl_set_level',__FILE__)
end if
if(level<mg_mesh%mgdim)then
  call multigrid_level(level)
else
  call multigrid_level(mg_mesh%mgdim)
end if
CALL ML_oft_hcurl%set_level(level)
SELECT TYPE(this=>ML_oft_hcurl%current_level)
  CLASS IS(oft_hcurl_fem)
    oft_hcurl=>this
  CLASS DEFAULT
    CALL oft_abort("Error setting HCurl level", "oft_hcurl_set_level", __FILE__)
END SELECT
! oft_hcurl=>ML_oft_hcurl%current_level
CALL ML_oft_bhcurl%set_level(level)
SELECT TYPE(this=>ML_oft_bhcurl%current_level)
  CLASS IS(oft_hcurl_bfem)
    oft_bhcurl=>this
  CLASS DEFAULT
    CALL oft_abort("Error setting boundary HCurl level", "oft_hcurl_set_level", __FILE__)
END SELECT
!---
oft_hcurl_level=level
oft_hcurl_lev=level
if(oft_hcurl_level>oft_hcurl_blevel.AND.oft_hcurl_blevel>0)oft_hcurl_lev=level-1
oft_hcurl_ops=>ML_oft_hcurl_ops(level)
DEBUG_STACK_POP
end subroutine oft_hcurl_set_level
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_setup
!---------------------------------------------------------------------------
!> Construct Nedelec H1(Curl) FE on each mesh level
!!
!! @note Highest supported representation is quadratic.
!!
!! @param[in] order Order of representation desired
!---------------------------------------------------------------------------
subroutine oft_hcurl_setup(order,minlev)
integer(i4), intent(in) :: order
integer(i4), optional, intent(in) :: minlev
integer(i4) :: i,j
DEBUG_STACK_PUSH
oft_hcurl_minlev=1
IF(PRESENT(minlev))oft_hcurl_minlev=minlev
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(A)')'**** Creating Nedelec H1(Curl) FE space'
  WRITE(*,'(2X,A,I4)')'Order  = ',order
  WRITE(*,'(2X,A,I4)')'Minlev = ',oft_hcurl_minlev
END IF
IF(mesh%type==3)hex_mesh=.TRUE.
!---Allocate multigrid operators
oft_hcurl_nlevels=mg_mesh%mgdim+(order-1)
IF(oft_hcurl_minlev<0)oft_hcurl_minlev=oft_hcurl_nlevels
ML_oft_hcurl%nlevels=oft_hcurl_nlevels
ML_oft_bhcurl%nlevels=oft_hcurl_nlevels
allocate(ML_oft_hcurl_ops(oft_hcurl_nlevels))
!---Set linear elements
do i=1,mg_mesh%mgdim-1
  IF(i<oft_hcurl_minlev)CYCLE
  ALLOCATE(oft_hcurl_fem::ML_oft_hcurl%levels(i)%fe)
  ALLOCATE(oft_hcurl_bfem::ML_oft_bhcurl%levels(i)%fe)
  call oft_hcurl_set_level(i)
  if(mg_mesh%level==mg_mesh%nbase)THEN
    ML_oft_hcurl%blevel=i
    ML_oft_bhcurl%blevel=i
    oft_hcurl_blevel=i
  END IF
  !---
  oft_hcurl%mesh=>mesh
  oft_hcurl%order=1
  oft_hcurl%dim=1
  oft_hcurl%type=oft_hcurl_id
  oft_hcurl%gstruct=(/0,1,0,0/)
  call oft_hcurl%setup(3)
  !---
  oft_bhcurl%mesh=>smesh
  oft_bhcurl%order=1
  oft_bhcurl%dim=1
  oft_bhcurl%type=oft_hcurl_id
  oft_bhcurl%gstruct=(/0,1,0/)
  call oft_bhcurl%setup(3)
end do
call multigrid_level(mg_mesh%mgdim)
!---Set high order elements
do i=1,order
  IF(mg_mesh%mgdim+i-1<oft_hcurl_minlev)CYCLE
  ALLOCATE(oft_hcurl_fem::ML_oft_hcurl%levels(mg_mesh%mgdim+i-1)%fe)
  ALLOCATE(oft_hcurl_bfem::ML_oft_bhcurl%levels(mg_mesh%mgdim+i-1)%fe)
  call oft_hcurl_set_level(mg_mesh%mgdim+i-1)
  !---
  oft_hcurl%mesh=>mesh
  oft_hcurl%order=i
  oft_hcurl%dim=1
  oft_hcurl%type=oft_hcurl_id
  IF(hex_mesh)THEN
    CALL hcurl_2d_grid(oft_hcurl%order-1, oft_hcurl%indsf)
    CALL hcurl_3d_grid(oft_hcurl%order-1, oft_hcurl%indsc)
    select case(oft_hcurl%order)
      case(1)
        oft_hcurl%gstruct=(/0,1,0,0/)
      case(2)
        oft_hcurl%gstruct=(/0,1,(oft_hcurl%order-1)**2 + 2*(oft_hcurl%order-1), &
                                2*(oft_hcurl%order-1)**3 + 3*(oft_hcurl%order-1)**2/)
      case(3)
        oft_hcurl%gstruct=(/0,1,(oft_hcurl%order-1)**2 + 2*(oft_hcurl%order-1), &
                                2*(oft_hcurl%order-1)**3 + 3*(oft_hcurl%order-1)**2/)
      case(4)
        oft_hcurl%gstruct=(/0,1,(oft_hcurl%order-1)**2 + 2*(oft_hcurl%order-1), &
                                2*(oft_hcurl%order-1)**3 + 3*(oft_hcurl%order-1)**2/)
      case default
        call oft_abort('Invalid polynomial degree (npmax=1 for hex grids)','oft_hcurl_setup',__FILE__)
    end select
  ELSE
    select case(oft_hcurl%order)
      case(1)
        oft_hcurl%gstruct=(/0,1,0,0/)
      case(2)
        oft_hcurl%gstruct=(/0,1,2,0/)
      case(3)
        oft_hcurl%gstruct=(/0,1,5,3/)
      case(4)
        oft_hcurl%gstruct=(/0,1,9,11/)
      case default
        call oft_abort('Invalid polynomial degree (npmax=4)','oft_hcurl_setup',__FILE__)
    end select
  END IF
  call oft_hcurl%setup(i*2+1)
  !---
  oft_bhcurl%mesh=>smesh
  oft_bhcurl%order=i
  oft_bhcurl%dim=1
  oft_bhcurl%type=oft_hcurl_id
  select case(oft_bhcurl%order)
    case(1)
      oft_bhcurl%gstruct=(/0,1,0/)
    case(2)
      oft_bhcurl%gstruct=(/0,1,2/)
    case(3)
      oft_bhcurl%gstruct=(/0,1,5/)
    case(4)
      oft_bhcurl%gstruct=(/0,1,9/)
    case default
      call oft_abort('Invalid polynomial degree (npmax=4)','oft_hcurl_setup',__FILE__)
  end select
  call oft_bhcurl%setup(i*2+1)
end do
IF(mg_mesh%mgdim>=oft_hcurl_minlev)THEN
  oft_hcurl_lin_level=mg_mesh%mgdim
  SELECT TYPE(this=>ML_oft_hcurl%levels(mg_mesh%mgdim)%fe)
    CLASS IS(oft_hcurl_fem)
      oft_hcurl_lin=>this
    CLASS DEFAULT
      CALL oft_abort("Error casting HCurl object", "oft_hcurl_setup", __FILE__)
  END SELECT
  ! oft_hcurl_lin=>ML_oft_hcurl%levels(mg_mesh%mgdim)%fe
  oft_hcurl_ops_lin=>ML_oft_hcurl_ops(mg_mesh%mgdim)
ELSE
  oft_hcurl_lin_level=-1
END IF
CALL oft_hcurl_set_level(oft_hcurl_nlevels)
IF(oft_env%head_proc)WRITE(*,*)
DEBUG_STACK_POP
end subroutine oft_hcurl_setup
!------------------------------------------------------------------------------
! SUBROUTINE: hcurl_2d_grid
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
subroutine hcurl_2d_grid(order,inds)
integer(i4), intent(in) :: order
integer(i4), pointer, intent(out) :: inds(:,:)
integer(i4) :: i,j,k,m
ALLOCATE(inds(3,order**2 + 2*order))
inds=-1
m=0
DO k=1,order
  DO i=1,k
    DO j=1,k
      IF(i<k.AND.j<k)CYCLE
      IF(j==1)THEN
        m=m+1
        inds(:,m)=(/i,i,1/)
        m=m+1
        inds(:,m)=(/i,i,2/)
      END IF
      m=m+1
      inds(:,m)=(/i,j,0/)
    END DO
  END DO
END DO
end subroutine hcurl_2d_grid
!------------------------------------------------------------------------------
! SUBROUTINE: hcurl_3d_grid
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
subroutine hcurl_3d_grid(order,inds)
integer(i4), intent(in) :: order
integer(i4), pointer, intent(out) :: inds(:,:)
integer(i4) :: i,j,k,l,m
ALLOCATE(inds(4,2*order**3 + 3*order**2))
inds=-1
m=0
DO k=1,order
  DO i=1,k
    DO j=1,k
      DO l=1,k
        IF(i<k.AND.j<k.AND.l<k)CYCLE
        IF(l==1)THEN
          m=m+1
          inds(:,m)=(/0,i,j,2/)
          m=m+1
          inds(:,m)=(/i,0,j,3/)
          m=m+1
          inds(:,m)=(/i,j,0,4/)
        END IF
        m=m+1
        inds(:,m)=(/i,j,l,0/)
        m=m+1
        inds(:,m)=(/i,j,l,1/)
      END DO
    END DO
  END DO
END DO
end subroutine hcurl_3d_grid
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_eval
!---------------------------------------------------------------------------
!> Evaluate Nedelec H1(Curl) interpolation function in the interior
!!
!! @note Evaluation is performed in logical coordinates with the resulting
!! vector in physical coordinates.
!!
!! @param[in] self Nedelec type for evaluation (fem)
!! @param[in] cell Cell for evaluation
!! @param[in] dof Element to evaluate
!! @param[in] f Position in cell in logical space [4]
!! @param[out] val Value of interpolation function (dof) at point (f) [3]
!! @param[in] gop Cell Jacobian matrix at point (f) [3,4]
!---------------------------------------------------------------------------
subroutine oft_hcurl_eval(self,cell,dof,f,val,gop)
class(oft_hcurl_fem), intent(in) :: self
integer(i4), intent(in) :: cell,dof
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(:,:)
real(r8), intent(out) :: val(:)
real(r8) :: cofs(4),fhex(6),gbary(3,6),dtmp,cords(4),f1(3),f2(3),f3(3),vtmp(4)
integer(i4) :: ed,etmp(2),fc,ftmp(3),i,j,fhtmp(4),ind,form
DEBUG_STACK_PUSH
IF(hex_mesh)THEN
  val=0.d0
  fhex=hex_get_bary(f)
  gbary=hex_get_bary_gop(gop)
  select case(self%cmap(dof)%type)
    case(2)
      etmp=hex_bary_ecoords(:,self%cmap(dof)%el)
      call orient_list2(self%mesh%lce(self%cmap(dof)%el,cell),etmp)
      cords(1:2)=fhex(hex_bary_efcoords(:,self%cmap(dof)%el))
      val=gbary(:,etmp(2))*cords(1)*cords(2)
    case(3)
      fhtmp=hex_bary_fcoords(:,self%cmap(dof)%el)
      IF(self%mesh%lcfo(self%cmap(dof)%el,cell)<0)fhtmp=fhtmp((/2,3,4,1/))
      call orient_listn(self%mesh%lcfo(self%cmap(dof)%el,cell), fhtmp, 4_i4)
      ind=self%cmap(dof)%ind
      form=self%indsf(3,ind)
      cords=fhex(fhtmp)
      IF(form==0)THEN
        vtmp(1:3)=dhpoly_2d(self%indsf(1:2,ind), cords(1:2))
        ! Part 1
        dtmp = (vtmp(2)*cords(3) - vtmp(1))*cords(4)*fhex(self%cmap(dof)%el)
        val = gbary(:,fhtmp(1))*dtmp
        ! Part 2
        dtmp = (vtmp(3)*cords(4) - vtmp(1))*cords(3)*fhex(self%cmap(dof)%el)
        val = val - gbary(:,fhtmp(2))*dtmp
      ELSE IF(form==1)THEN
        dtmp = hpoly_bary(self%indsf(2,ind), cords(2))*cords(4)*fhex(self%cmap(dof)%el)
        val=gbary(:,fhtmp(1))*dtmp
      ELSE IF(form==2)THEN
        dtmp = hpoly_bary(self%indsf(1,ind), cords(1))*cords(3)*fhex(self%cmap(dof)%el)
        val=gbary(:,fhtmp(2))*dtmp
      END IF
    case(4)
      ind=self%cmap(dof)%ind
      form=self%indsc(4,ind)
      IF(form<2)THEN
        vtmp=dhpoly_3d(self%indsc(1:3,ind), fhex(1:3))
        ! Part 1
        dtmp=(vtmp(2)*fhex(6) - vtmp(1))*fhex(4)*fhex(5)
        f1 = gbary(:,1)*dtmp
        ! Part 2
        dtmp=(vtmp(3)*fhex(4) - vtmp(1))*fhex(5)*fhex(6)
        f2 = gbary(:,2)*dtmp
        ! Part 3
        dtmp=(vtmp(4)*fhex(5) - vtmp(1))*fhex(4)*fhex(6)
        f3 = gbary(:,3)*dtmp
      END IF
      IF(form==0)THEN
        val = f1 - f2 + f3
      ELSE IF(form==1)THEN
        val = f1 - f2 - f3
      ELSE IF(form==2)THEN
        dtmp = hpoly_2d(self%indsc(2:3,ind), fhex(2:3))*fhex(4)*fhex(5)
        val=gbary(:,1)*dtmp
      ELSE IF(form==3)THEN
        dtmp = hpoly_2d(self%indsc((/1,3/),ind), fhex((/1,3/)))*fhex(5)*fhex(6)
        val=gbary(:,2)*dtmp
      ELSE IF(form==4)THEN
        dtmp = hpoly_2d(self%indsc(1:2,ind), fhex(1:2))*fhex(4)*fhex(6)
        val=gbary(:,3)*dtmp
      END IF
    case default
      call oft_abort('Invalid element geometry','oft_hcurl_ceval',__FILE__)
  end select
ELSE
  IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
  val=0.d0; cofs=0.d0
  select case(self%cmap(dof)%type)
    case(2)
      etmp=oriented_edges(:,self%cmap(dof)%el)
      cofs(etmp(1)) = -f(etmp(2))
      cofs(etmp(2)) = f(etmp(1))
    case(3)
      ftmp=oriented_faces(:,self%cmap(dof)%el)
      call oft_hcurl_evalf(self%order,ftmp,self%cmap(dof)%ind,f,cofs)
    case(4)
      call oft_hcurl_evalc(self%order,self%cmap(dof)%ind,f,cofs)
    case default
      call oft_abort('Invalid element geometry','oft_hcurl_eval',__FILE__)
  end select
  !---Sum contributions
  do i=1,4
    val=val+gop(:,i)*cofs(i)
  end do
END IF
DEBUG_STACK_POP
end subroutine oft_hcurl_eval
!---------------------------------------------------------------------------
! SUBROUTINE: oft_bhcurl_eval
!---------------------------------------------------------------------------
!> Evaluate Nedelec H1(Curl) interpolation function on the boundary
!!
!! @note Evaluation is performed in logical coordinates with the resulting
!! vector in physical coordinates.
!!
!! @param[in] self Nedelec type for evaluation (bfem)
!! @param[in] face Face for evaluation
!! @param[in] dof Element to evaluate
!! @param[in] f Position on face in logical space [4]
!! @param[out] val Value of interpolation function (dof) at point (f) [3]
!! @param[in] gop Face Jacobian matrix at point (f) [3,3]
!---------------------------------------------------------------------------
subroutine oft_bhcurl_eval(self,face,dof,f,val,gop)
class(oft_bfem_type), intent(in) :: self
integer(i4), intent(in) :: face,dof
real(r8), intent(in) :: f(:)
real(r8), optional, intent(in) :: gop(:,:)
real(r8), intent(out) :: val(3)
real(r8) :: grads(3,4),cofs(4)
integer(i4) :: ed,etmp(2),fc,ftmp(3),i
DEBUG_STACK_PUSH
grads=1.d0
if(present(gop))grads(:,1:3)=gop
val=0.d0; cofs=0.d0
select case(self%cmap(dof)%type)
  case(1)
    call oft_abort('Invalid element geometry (no point DOFs)','oft_bhcurl_eval',__FILE__)
  case(2)
    ed=self%mesh%lce(self%cmap(dof)%el,face)
    etmp=self%mesh%cell_ed(:,self%cmap(dof)%el)
    call orient_list2(ed,etmp)
    call oft_hcurl_evale(self%order,etmp,self%cmap(dof)%ind,f,cofs)
  case(3)
    fc=self%mesh%lco(face)
    ftmp=(/1,2,3/)
    call orient_listn(fc,ftmp,3_i4)
    call oft_hcurl_evalf(self%order,ftmp,self%cmap(dof)%ind,f,cofs)
end select
!---Sum contributions
do i=1,4
  val=val+grads(:,i)*cofs(i)
end do
DEBUG_STACK_POP
end subroutine oft_bhcurl_eval
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_evale
!---------------------------------------------------------------------------
!> Evaluate edge based interpolation functions
!---------------------------------------------------------------------------
subroutine oft_hcurl_evale(order,ed,dof,f,val)
integer(i4), intent(in) :: order
integer(i4), intent(in) :: ed(2),dof
real(r8), intent(in) :: f(:)
real(r8), intent(out) :: val(4)
real(r8) :: f1,f2
DEBUG_STACK_PUSH
f1=f(ed(1))
f2=f(ed(2))
!---Linear DOF
val(ed(1)) = -f2
val(ed(2)) = f1
DEBUG_STACK_POP
end subroutine oft_hcurl_evale
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_evalf
!---------------------------------------------------------------------------
!> Evaluate face based interpolation functions
!---------------------------------------------------------------------------
subroutine oft_hcurl_evalf(order,fc,dof,f,val)
integer(i4), intent(in) :: order
integer(i4), intent(in) :: fc(3),dof
real(r8), intent(in) :: f(:)
real(r8), intent(out) :: val(4)
real(r8) :: f1,f2,f3,y1,y2,y3
DEBUG_STACK_PUSH
f1=f(fc(1)); f2=f(fc(2)); f3=f(fc(3))
SELECT CASE(dof)
!---Quadratic DOF
CASE(1)
  y1 = -2*f2*f3
  y2 = -2*f1*f3
  y3 = 2*f1*f2
CASE(2)
  y1 = f2*f3
  y2 = -f1*f3
  y3 = 0.d0
!---Cubic DOF
CASE(3)
  y1 = 2*f2*f3*(f2-f3)
  y2 = 2*f1*f3*(f1-f3)
  y3 = 2*f1*f2*(-f1-f2+2*f3)
CASE(4)
  y1 =  f2*f3*(-f1 - f2 + f3)
  y2 =  -f1*f3*(-f1 - f2 + f3)
  y3 =  0.d0
CASE(5)
  y1 = 2*f2*f3*(-2*f1+f2)
  y2 = 2*f1*f3*(-f1+2*f2)
  y3 = 2*f1*f2*(f1-f2)
!---Quartic DOF
CASE(6)
  y1 =  f2*f3*(2.d0*f1**2 - 2.d0*f2**2 + 8.d0*f2*f3 - 2.d0*f3**2)
  y2 =  f1*f3*(-2.d0*f1**2 + 8.d0*f1*f3 + 2.d0*f2**2 - 2.d0*f3**2)
  y3 =  f1*f2*(2.d0*f1**2 + 4.d0*f1*f2 - 16.d0*f1*f3 + 2.d0*f2**2 - 16.d0*f2*f3 + 6.d0*f3**2)
CASE(7)
  y1 =  f2*f3*((1.5d0*(f1 + f2 - f3)**2)/(f1 + f2 + f3)**2 - 0.5d0)*(f1 + f2 + f3)**2
  y2 =  -f1*f3*(1.5d0*(f1 + f2 - f3)**2 - 0.5d0*(f1 + f2 + f3)**2)
  y3 =  0.d0
CASE(8)
  y1 =  f2*f3*(2.d0*f1**2 + 4.d0*f1*f2 - 4.d0*f1*f3 - 2.d0*f2**2 + 2.d0*f2*f3)
  y2 =  f1*f3*(2.d0*f1**2 - 4.d0*f1*f2 - 2.d0*f1*f3 - 2.d0*f2**2 + 4.d0*f2*f3)
  y3 =  f1*f2*(-2.d0*f1**2 + 4.d0*f1*f3 + 2.d0*f2**2 - 4.d0*f2*f3)
CASE(9)
  y1 =  f2*f3*(-6.d0*f1**2 + 12.d0*f1*f2 - 2.d0*f2**2)
  y2 =  f1*f3*(-2.d0*f1**2 + 12.d0*f1*f2 - 6.d0*f2**2)
  y3 =  f1*f2*(2.d0*f1**2 - 6.d0*f1*f2 + 2.d0*f2**2)
END SELECT
val(fc(1))=y1; val(fc(2))=y2; val(fc(3))=y3
DEBUG_STACK_POP
end subroutine oft_hcurl_evalf
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_evalc
!---------------------------------------------------------------------------
!> Evaluate cell based interpolation functions
!---------------------------------------------------------------------------
subroutine oft_hcurl_evalc(order,dof,f,val)
integer(i4), intent(in) :: order
integer(i4), intent(in) :: dof
real(r8), intent(in) :: f(:)
real(r8), intent(out) :: val(4)
real(r8) :: f1,f2,f3,f4
DEBUG_STACK_PUSH
!---
f1=f(1)
f2=f(2)
f3=f(3)
f4=f(4)
!---
SELECT CASE(dof)
!---Cubic DOF
CASE(1)
  val(1) =  -2*f2*f3*f4
  val(2) =  -2*f1*f3*f4
  val(3) =  2*f1*f2*f4
  val(4) =  -2*f1*f2*f3
CASE(2)
  val(1) =  -2*f2*f3*f4
  val(2) =  -2*f1*f3*f4
  val(3) =  2*f1*f2*f4
  val(4) =  2*f1*f2*f3
CASE(3)
  val(1) =  f2*f3*f4
  val(2) =  -f1*f3*f4
  val(3) =  0.d0
  val(4) =  0.d0
!---Quartic DOF
CASE(4)
  val(1) =  f2*f3*f4*(4.d0*f1 + 2.d0*f2 + 2.d0*f3 - 2.d0*f4)
  val(2) =  f1*f3*f4*(2.d0*f1 + 4.d0*f2 + 2.d0*f3 - 2.d0*f4)
  val(3) =  f1*f2*f4*(-2.d0*f1 - 2.d0*f2 + 2.d0*f4)
  val(4) =  f1*f2*f3*(2.d0*f1 + 2.d0*f2 + 2.d0*f3 - 4.d0*f4)
CASE(5)
  val(1) =  f2*f3*f4*(2.d0*f2 + 2.d0*f3 - 2.d0*f4)
  val(2) =  f1*f3*f4*(2.d0*f1 + 2.d0*f3 - 2.d0*f4)
  val(3) =  f1*f2*f4*(-2.d0*f1 - 2.d0*f2 - 4.d0*f3 + 2.d0*f4)
  val(4) =  f1*f2*f3*(-2.d0*f1 - 2.d0*f2 - 2.d0*f3 + 4.d0*f4)
CASE(6)
  val(1) =  f2*f3*f4*(-f1 - f2 - f3 + f4)
  val(2) =  f1*f3*f4*(f1 + f2 + f3 - f4)
  val(3) =  0.d0
  val(4) =  0.d0
CASE(7)
  val(1) =  f2*f3*f4*(2.d0*f2 - 2.d0*f3)
  val(2) =  f1*f3*f4*(2.d0*f1 - 2.d0*f3)
  val(3) =  f1*f2*f4*(-2.d0*f1 - 2.d0*f2 + 4.d0*f3)
  val(4) =  f1*f2*f3*(2.d0*f1 + 2.d0*f2 - 2.d0*f3)
CASE(8)
  val(1) =  f2*f3*f4*(2.d0*f2 - 2.0*f3)
  val(2) =  f1*f3*f4*(2.d0*f1 - 2.0*f3)
  val(3) =  f1*f2*f4*(-2.d0*f1 - 2.d0*f2 + 4.d0*f3)
  val(4) =  f1*f2*f3*(-2.d0*f1 - 2.d0*f2 + 2.d0*f3)
CASE(9)
  val(1) =  f2*f3*f4*(-f1 - f2 + f3)
  val(2) =  f1*f3*f4*(f1 + f2 - f3)
  val(3) =  0.d0
  val(4) =  0.d0
CASE(10)
  val(1) =  f2*f3*f4*(-4.d0*f1 + 2.d0*f2)
  val(2) =  f1*f3*f4*(-2.d0*f1 + 4.d0*f2)
  val(3) =  f1*f2*f4*(2.d0*f1 - 2.d0*f2)
  val(4) =  f1*f2*f3*(-2.d0*f1 + 2.d0*f2)
CASE(11)
  val(1) =  f2*f3*f4*(-4.d0*f1 + 2.d0*f2)
  val(2) =  f1*f3*f4*(-2.d0*f1 + 4.d0*f2)
  val(3) =  f1*f2*f4*(2.d0*f1 - 2.d0*f2)
  val(4) =  f1*f2*f3*(2.d0*f1 - 2.d0*f2)
END SELECT
DEBUG_STACK_POP
end subroutine oft_hcurl_evalc
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_eval_all
!---------------------------------------------------------------------------
!> Evaluate all lagrange interpolation functions
!!
!! @note Evaluation is performed in logical coordinates.
!!
!! @param[in] self Lagrange type for evaluation
!! @param[in] cell Cell for evaluation
!! @param[in] f Position in cell in logical space
!! @param[out] rop Value of interpolation functions at point (f) [3,ncdofs]
!! @param[in] gop Cell Jacobian matrix at point (f) [3,4]
!---------------------------------------------------------------------------
subroutine oft_hcurl_eval_all(self,cell,f,rop,gop)
class(oft_hcurl_fem), intent(in) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(4)
real(r8), intent(in) :: gop(3,4)
real(r8), contiguous, intent(out) :: rop(:,:)
integer(i4) :: i,j,etmp(2),fhtmp(4),offset
real(r8) :: fhex(6),gbary(3,6),dtmp,cords(4),f1(3),f2(3),f3(3),vtmp(4)
DEBUG_STACK_PUSH
IF(hex_mesh)THEN
  fhex=hex_get_bary(f)
  gbary=hex_get_bary_gop(gop)
  !---Edges
  DO i=1,12
    offset=(i-1)*self%gstruct(2)
    etmp=hex_bary_ecoords(:,i)
    call orient_list2(self%mesh%lce(i,cell), etmp)
    cords(1:2)=fhex(etmp)
    cords(3:4)=fhex(hex_bary_efcoords(:,i))
    rop(:,offset+1)=gbary(:,etmp(2))*cords(3)*cords(4)
  END DO
  IF(self%gstruct(3)>0)THEN
    !---Faces
    DO i=1,6
      offset=(i-1)*self%gstruct(3)+12*self%gstruct(2)
      fhtmp=hex_bary_fcoords(:,i)
      IF(self%mesh%lcfo(i,cell)<0)fhtmp=fhtmp((/2,3,4,1/))
      call orient_listn(self%mesh%lcfo(i,cell), fhtmp, 4_i4)
      cords=fhex(fhtmp)
      DO j=1,self%gstruct(3)
        IF(self%indsf(3,j)==0)THEN
          vtmp(1:3) = dhpoly_2d(self%indsf(1:2,j), cords(1:2))
          ! Part 1
          dtmp = (vtmp(2)*cords(3) - vtmp(1))*cords(4)*fhex(i)
          rop(:,offset+j) = gbary(:,fhtmp(1))*dtmp
          ! Part 2
          dtmp = (vtmp(3)*cords(4) - vtmp(1))*cords(3)*fhex(i)
          rop(:,offset+j) = rop(:,offset+j) - gbary(:,fhtmp(2))*dtmp
        ELSE IF(self%indsf(3,j)==1)THEN
          dtmp = hpoly_bary(self%indsf(2,j), cords(2))*cords(4)*fhex(i)
          rop(:,offset+j)=gbary(:,fhtmp(1))*dtmp
        ELSE IF(self%indsf(3,j)==2)THEN
          dtmp = hpoly_bary(self%indsf(1,j), cords(1))*cords(3)*fhex(i)
          rop(:,offset+j)=gbary(:,fhtmp(2))*dtmp
        END IF
      END DO
    END DO
  END IF
  IF(self%gstruct(4)>0)THEN
    !---Cell
    offset=6*self%gstruct(3)+12*self%gstruct(2)
    DO j=1,self%gstruct(4)
      IF(self%indsc(4,j)<2)THEN
        vtmp = dhpoly_3d(self%indsc(1:3,j), fhex(1:3))
        ! Part 1
        f1 = gbary(:,1)*(vtmp(2)*fhex(6) - vtmp(1))*fhex(4)*fhex(5)
        ! Part 2
        f2 = gbary(:,2)*(vtmp(3)*fhex(4) - vtmp(1))*fhex(5)*fhex(6)
        ! Part 3
        f3 = gbary(:,3)*(vtmp(4)*fhex(5) - vtmp(1))*fhex(4)*fhex(6)
      END IF
      IF(self%indsc(4,j)==0)THEN
        rop(:,offset+j) = f1 - f2 + f3
      ELSE IF(self%indsc(4,j)==1)THEN
        rop(:,offset+j) = f1 - f2 - f3
      ELSE IF(self%indsc(4,j)==2)THEN
        dtmp = hpoly_2d(self%indsc(2:3,j), fhex(2:3))*fhex(4)*fhex(5)
        rop(:,offset+j)=gbary(:,1)*dtmp
      ELSE IF(self%indsc(4,j)==3)THEN
        dtmp = hpoly_2d(self%indsc((/1,3/),j), fhex((/1,3/)))*fhex(5)*fhex(6)
        rop(:,offset+j)=gbary(:,2)*dtmp
      ELSE IF(self%indsc(4,j)==4)THEN
        dtmp = hpoly_2d(self%indsc(1:2,j), fhex(1:2))*fhex(4)*fhex(6)
        rop(:,offset+j)=gbary(:,3)*dtmp
      END IF
    END DO
  END IF
ELSE
  IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
  select case(self%order)
    case(1)
      DO i=1,6
        etmp=oriented_edges(:,i)
        rop(:,i) = -f(etmp(2))*gop(:,etmp(1)) &
                 +  f(etmp(1))*gop(:,etmp(2))
      END DO
    case(2)
      call oft_hcurl_eval_all2(self,cell,f,rop,gop)
    case(3)
      call oft_hcurl_eval_all3(self,cell,f,rop,gop)
    case(4)
      call oft_hcurl_eval_all4(self,cell,f,rop,gop)
    case default ! Fall back to normal evaluation
      DO i=1,self%nce
        call oft_hcurl_eval(self,cell,i,f,rop(:,i),gop)
      END DO
  end select
END IF
DEBUG_STACK_POP
end subroutine oft_hcurl_eval_all
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_eval_all2
!---------------------------------------------------------------------------
!> Evaluate all lagrange interpolation functions (quadratic)
!!
!! @note Evaluation is performed in logical coordinates.
!!
!! @param[in] self Lagrange type for evaluation
!! @param[in] cell Cell for evaluation
!! @param[in] f Position in cell in logical space
!! @param[out] rop Value of interpolation functions at point (f) [ncdofs]
!---------------------------------------------------------------------------
subroutine oft_hcurl_eval_all2(self,cell,f,rop,gop)
class(oft_fem_type), intent(in) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(4)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: rop(3,14)
integer(i4) :: i,etmp(2),ftmp(3)
real(r8) :: f1,f2,f3,u1(3),u2(3),u3(3)
DEBUG_STACK_PUSH
DO i=1,6
  etmp=oriented_edges(:,i)
  rop(:,i) = -f(etmp(2))*gop(:,etmp(1)) &
           +  f(etmp(1))*gop(:,etmp(2))
END DO
DO i=1,4
  ftmp=oriented_faces(:,i)
  f1 = f(ftmp(1)); f2 = f(ftmp(2)); f3 = f(ftmp(3))
  u1 = gop(:,ftmp(1)); u2 = gop(:,ftmp(2)); u3 = gop(:,ftmp(3))
  !---Quadratic DOF
  rop(:,(i-1)*2+7) = -2*f2*f3*u1 &
               -2*f1*f3*u2 &
             +  2*f1*f2*u3
  !
  rop(:,(i-1)*2+8) = f2*f3*u1 &
               -f1*f3*u2 &
              + 0*u3
END DO
DEBUG_STACK_POP
end subroutine oft_hcurl_eval_all2
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_eval_all3
!---------------------------------------------------------------------------
!> Evaluate all lagrange interpolation functions (cubic)
!!
!! @note Evaluation is performed in logical coordinates.
!!
!! @param[in] self Lagrange type for evaluation
!! @param[in] cell Cell for evaluation
!! @param[in] f Position in cell in logical space
!! @param[out] rop Value of interpolation functions at point (f) [ncdofs]
!---------------------------------------------------------------------------
subroutine oft_hcurl_eval_all3(self,cell,f,rop,gop)
class(oft_fem_type), intent(in) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(4)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: rop(3,29)
integer(i4) :: i,etmp(2),ftmp(3)
real(r8) :: f1,f2,f3,f4,u1(3),u2(3),u3(3)
DEBUG_STACK_PUSH
DO i=1,6
  etmp=oriented_edges(:,i)
  rop(:,i) = -f(etmp(2))*gop(:,etmp(1)) &
           +  f(etmp(1))*gop(:,etmp(2))
END DO
DO i=1,4
  ftmp=oriented_faces(:,i)
  f1 = f(ftmp(1)); f2 = f(ftmp(2)); f3 = f(ftmp(3))
  u1 = gop(:,ftmp(1)); u2 = gop(:,ftmp(2)); u3 = gop(:,ftmp(3))
  !---Quadratic DOF
  rop(:,(i-1)*5+7) = -2*f2*f3*u1 &
               -2*f1*f3*u2 &
             +  2*f1*f2*u3
  !
  rop(:,(i-1)*5+8) = f2*f3*u1 &
               -f1*f3*u2 &
              + 0*u3
  !---Cubic DOF
  rop(:,(i-1)*5+9) = 2*f2*f3*(f2-f3)*u1 &
              + 2*f1*f3*(f1-f3)*u2 &
              + 2*f1*f2*(-f1-f2+2*f3)*u3
  !
  rop(:,(i-1)*5+10) = f2*f3*(-f1 - f2 + f3)*u1 &
               -f1*f3*(-f1 - f2 + f3)*u2 &
              + 0*u3
  !
  rop(:,(i-1)*5+11) = 2*f2*f3*(-2*f1+f2)*u1 &
              + 2*f1*f3*(-f1+2*f2)*u2 &
              + 2*f1*f2*(f1-f2)*u3
END DO
f1 = f(1); f2 = f(2); f3 = f(3); f4 = f(4)
!---Cubic DOF
rop(:,27) = -2*f2*f3*f4*gop(:,1) &
            -2*f1*f3*f4*gop(:,2) &
          +  2*f1*f2*f4*gop(:,3) &
            -2*f1*f2*f3*gop(:,4)
!
rop(:,28) = -2*f2*f3*f4*gop(:,1) &
            -2*f1*f3*f4*gop(:,2) &
          +  2*f1*f2*f4*gop(:,3) &
          +  2*f1*f2*f3*gop(:,4)
!
rop(:,29) = f2*f3*f4*gop(:,1) &
            -f1*f3*f4*gop(:,2) &
          +  0*gop(:,3) &
          +  0*gop(:,4)
DEBUG_STACK_POP
end subroutine oft_hcurl_eval_all3
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_eval_all4
!---------------------------------------------------------------------------
!> Evaluate all lagrange interpolation functions (quartic)
!!
!! @note Evaluation is performed in logical coordinates.
!!
!! @param[in] self Lagrange type for evaluation
!! @param[in] cell Cell for evaluation
!! @param[in] f Position in cell in logical space
!! @param[out] rop Value of interpolation functions at point (f) [ncdofs]
!---------------------------------------------------------------------------
subroutine oft_hcurl_eval_all4(self,cell,f,rop,gop)
class(oft_fem_type), intent(in) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(4)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: rop(3,53)
integer(i4) :: i,etmp(2),ftmp(3)
real(r8) :: f1,f2,f3,f4,u1(3),u2(3),u3(3)
DEBUG_STACK_PUSH
DO i=1,6
  etmp=oriented_edges(:,i)
  rop(:,i) = -f(etmp(2))*gop(:,etmp(1)) &
           +  f(etmp(1))*gop(:,etmp(2))
END DO
DO i=1,4
  ftmp=oriented_faces(:,i)
  f1 = f(ftmp(1)); f2 = f(ftmp(2)); f3 = f(ftmp(3))
  u1 = gop(:,ftmp(1)); u2 = gop(:,ftmp(2)); u3 = gop(:,ftmp(3))
  !---Quadratic DOF
  rop(:,(i-1)*9+7) = -2*f2*f3*u1 &
               -2*f1*f3*u2 &
             +  2*f1*f2*u3
  !
  rop(:,(i-1)*9+8) = f2*f3*u1 &
               -f1*f3*u2 &
              + 0*u3
  !---Cubic DOF
  rop(:,(i-1)*9+9) = 2*f2*f3*(f2-f3)*u1 &
              + 2*f1*f3*(f1-f3)*u2 &
              + 2*f1*f2*(-f1-f2+2*f3)*u3
  !
  rop(:,(i-1)*9+10) = f2*f3*(-f1 - f2 + f3)*u1 &
               -f1*f3*(-f1 - f2 + f3)*u2 &
              + 0*u3
  !
  rop(:,(i-1)*9+11) = 2*f2*f3*(-2*f1+f2)*u1 &
              + 2*f1*f3*(-f1+2*f2)*u2 &
              + 2*f1*f2*(f1-f2)*u3
  !---Quartic DOF
  rop(:,(i-1)*9+12) = f2*f3*(2.d0*f1**2 - 2.d0*f2**2 + 8.d0*f2*f3 - 2.d0*f3**2)*u1 &
              + f1*f3*(-2.d0*f1**2 + 8.d0*f1*f3 + 2.d0*f2**2 - 2.d0*f3**2)*u2 &
              + f1*f2*(2.d0*f1**2 + 4.d0*f1*f2 - 16.d0*f1*f3 + 2.d0*f2**2 - 16.d0*f2*f3 + 6.d0*f3**2)*u3
  !
  rop(:,(i-1)*9+13) = (f2*f3*((1.5d0*(f1 + f2 - f3)**2)/(f1 + f2 + f3)**2 - 0.5d0)*(f1 + f2 + f3)**2)*u1 &
               -f1*f3*(1.5d0*(f1 + f2 - f3)**2 - 0.5d0*(f1 + f2 + f3)**2)*u2 &
              + 0*u3
  !
  rop(:,(i-1)*9+14) = f2*f3*(2.d0*f1**2 + 4.d0*f1*f2 - 4.d0*f1*f3 - 2.d0*f2**2 + 2.d0*f2*f3)*u1 &
              + f1*f3*(2.d0*f1**2 - 4.d0*f1*f2 - 2.d0*f1*f3 - 2.d0*f2**2 + 4.d0*f2*f3)*u2 &
              + f1*f2*(-2.d0*f1**2 + 4.d0*f1*f3 + 2.d0*f2**2 - 4.d0*f2*f3)*u3
  !
  rop(:,(i-1)*9+15) = f2*f3*(-6.d0*f1**2 + 12.d0*f1*f2 - 2.d0*f2**2)*u1 &
              + f1*f3*(-2.d0*f1**2 + 12.d0*f1*f2 - 6.d0*f2**2)*u2 &
              + f1*f2*(2.d0*f1**2 - 6.d0*f1*f2 + 2.d0*f2**2)*u3
END DO
f1 = f(1); f2 = f(2); f3 = f(3); f4 = f(4)
!---Cubic DOF
rop(:,43) = -2*f2*f3*f4*gop(:,1) &
            -2*f1*f3*f4*gop(:,2) &
          +  2*f1*f2*f4*gop(:,3) &
            -2*f1*f2*f3*gop(:,4)
!
rop(:,44) = -2*f2*f3*f4*gop(:,1) &
            -2*f1*f3*f4*gop(:,2) &
          +  2*f1*f2*f4*gop(:,3) &
          +  2*f1*f2*f3*gop(:,4)
!
rop(:,45) = f2*f3*f4*gop(:,1) &
            -f1*f3*f4*gop(:,2) &
          +  0*gop(:,3) &
          +  0*gop(:,4)
!---Quartic DOF
rop(:,46) = f2*f3*f4*(4.d0*f1 + 2.d0*f2 + 2.d0*f3 - 2.d0*f4)*gop(:,1) &
          + f1*f3*f4*(2.d0*f1 + 4.d0*f2 + 2.d0*f3 - 2.d0*f4)*gop(:,2) &
          + f1*f2*f4*(-2.d0*f1 - 2.d0*f2 + 2.d0*f4)*gop(:,3) &
          + f1*f2*f3*(2.d0*f1 + 2.d0*f2 + 2.d0*f3 - 4.d0*f4)*gop(:,4)
!
rop(:,47) = f2*f3*f4*(2.d0*f2 + 2.d0*f3 - 2.d0*f4)*gop(:,1) &
          + f1*f3*f4*(2.d0*f1 + 2.d0*f3 - 2.d0*f4)*gop(:,2) &
          + f1*f2*f4*(-2.d0*f1 - 2.d0*f2 - 4.d0*f3 + 2.d0*f4)*gop(:,3) &
          + f1*f2*f3*(-2.d0*f1 - 2.d0*f2 - 2.d0*f3 + 4.d0*f4)*gop(:,4)
!
rop(:,48) = f2*f3*f4*(-f1 - f2 - f3 + f4)*gop(:,1) &
          + f1*f3*f4*(f1 + f2 + f3 - f4)*gop(:,2) &
          + 0*gop(:,3) &
          + 0*gop(:,4)
!
rop(:,49) = f2*f3*f4*(2.d0*f2 - 2.d0*f3)*gop(:,1) &
          + f1*f3*f4*(2.d0*f1 - 2.d0*f3)*gop(:,2) &
          + f1*f2*f4*(-2.d0*f1 - 2.d0*f2 + 4.d0*f3)*gop(:,3) &
          + f1*f2*f3*(2.d0*f1 + 2.d0*f2 - 2.d0*f3)*gop(:,4)
!
rop(:,50) = f2*f3*f4*(2.d0*f2 - 2.0*f3)*gop(:,1) &
          + f1*f3*f4*(2.d0*f1 - 2.0*f3)*gop(:,2) &
          + f1*f2*f4*(-2.d0*f1 - 2.d0*f2 + 4.d0*f3)*gop(:,3) &
          + f1*f2*f3*(-2.d0*f1 - 2.d0*f2 + 2.d0*f3)*gop(:,4)
!
rop(:,51) = f2*f3*f4*(-f1 - f2 + f3)*gop(:,1) &
          + f1*f3*f4*(f1 + f2 - f3)*gop(:,2) &
          + 0*gop(:,3) &
          + 0*gop(:,4)
!
rop(:,52) = f2*f3*f4*(-4.d0*f1 + 2.d0*f2)*gop(:,1) &
          + f1*f3*f4*(-2.d0*f1 + 4.d0*f2)*gop(:,2) &
          + f1*f2*f4*(2.d0*f1 - 2.d0*f2)*gop(:,3) &
          + f1*f2*f3*(-2.d0*f1 + 2.d0*f2)*gop(:,4)
!
rop(:,53) = f2*f3*f4*(-4.d0*f1 + 2.d0*f2)*gop(:,1) &
          + f1*f3*f4*(-2.d0*f1 + 4.d0*f2)*gop(:,2) &
          + f1*f2*f4*(2.d0*f1 - 2.d0*f2)*gop(:,3) &
          + f1*f2*f3*(2.d0*f1 - 2.d0*f2)*gop(:,4)
DEBUG_STACK_POP
end subroutine oft_hcurl_eval_all4
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_ceval
!---------------------------------------------------------------------------
!> Evaluate Nedelec H1(Curl) curl function in the interior
!!
!! @note Evaluation is performed in logical coordinates with the resulting
!! vector in, and curl with respect to, physical coordinates.
!!
!! @param[in] self Nedelec type for evaluation
!! @param[in] cell Cell for evaluation
!! @param[in] dof Element to evaluate
!! @param[in] f Position in cell in logical space [4]
!! @param[out] val Curl of nedelec element (dof) at point (f) [3]
!! @param[in] gop Cell Jacobian matrix at point (f) [3,4]
!---------------------------------------------------------------------------
subroutine oft_hcurl_ceval(self,cell,dof,f,val,gop)
class(oft_hcurl_fem), intent(in) :: self
integer(i4), intent(in) :: cell,dof
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
integer(i4) :: i,j,ed,etmp(2),fc,ftmp(3),fhtmp(4),ind,form
real(r8) :: fhex(6),gbary(3,6),dtmp,hcgop(3,3)
real(r8) :: cords(4),f1(3),f2(3),f3(3),vec(3,3),vtmp(4)
DEBUG_STACK_PUSH
IF(hex_mesh)THEN
  val=0.d0
  fhex=hex_get_bary(f)
  gbary=hex_get_bary_gop(gop)
  hcgop(:,1)=cross_product(gop(:,1), gop(:,2)) ! (1 x 2)
  hcgop(:,2)=cross_product(gop(:,1), gop(:,3)) ! (1 x 3)
  hcgop(:,3)=cross_product(gop(:,2), gop(:,3)) ! (2 x 3)
  select case(self%cmap(dof)%type)
    case(2)
      etmp=hex_bary_ecoords(:,self%cmap(dof)%el)
      call orient_list2(self%mesh%lce(self%cmap(dof)%el,cell),etmp)
      cords(1:2)=fhex(etmp)
      cords(3:4)=fhex(hex_bary_efcoords(:,self%cmap(dof)%el))
      val=hex_get_bary_cgop(hcgop,hex_bary_efcoords(1,self%cmap(dof)%el),etmp(2))*cords(4) &
         +hex_get_bary_cgop(hcgop,hex_bary_efcoords(2,self%cmap(dof)%el),etmp(2))*cords(3)
    case(3)
      fhtmp=hex_bary_fcoords(:,self%cmap(dof)%el)
      IF(self%mesh%lcfo(self%cmap(dof)%el,cell)<0)fhtmp=fhtmp((/2,3,4,1/))
      call orient_listn(self%mesh%lcfo(self%cmap(dof)%el,cell), fhtmp, 4_i4)
      cords=fhex(fhtmp)
      vec(:,1)=hex_get_bary_cgop(hcgop,self%cmap(dof)%el,fhtmp(1))
      vec(:,2)=hex_get_bary_cgop(hcgop,self%cmap(dof)%el,fhtmp(2))
      vec(:,3)=hex_get_bary_cgop(hcgop,fhtmp(1),fhtmp(2))
      ind=self%cmap(dof)%ind
      form=self%indsf(3,ind)
      IF(form==0)THEN
        vtmp(1:3)=dhpoly_2d(self%indsf(1:2,ind), cords(1:2))
        ! Part 1
        dtmp=(vtmp(2)*cords(3) - vtmp(1))*cords(4)
        val = vec(:,1)*dtmp
        !
        dtmp= vtmp(2)*cords(3) - vtmp(1) &
            - d2hpoly_2d(self%indsf(1:2,ind), 1,2, cords(1:2))*cords(3)*cords(4) &
            + vtmp(3)*cords(4)
        val = val + vec(:,3)*dtmp*fhex(self%cmap(dof)%el)
        ! Part 2
        dtmp=(vtmp(3)*cords(4) - vtmp(1))*cords(3)
        val = val - vec(:,2)*dtmp
        !
        dtmp= vtmp(3)*cords(4) - vtmp(1) &
            - d2hpoly_2d(self%indsf(1:2,ind), 2,1, cords(1:2))*cords(3)*cords(4) &
            + vtmp(2)*cords(3)
        val = val + vec(:,3)*dtmp*fhex(self%cmap(dof)%el)
      ELSE IF(form==1)THEN
        !
        dtmp = hpoly_bary(self%indsf(2,ind), cords(2))*cords(4)
        val = vec(:,1)*dtmp
        !
        dtmp = dhpoly_bary(self%indsf(2,ind), cords(2))*cords(4) &
             - hpoly_bary(self%indsf(2,ind), cords(2))
        val = val - vec(:,3)*dtmp*fhex(self%cmap(dof)%el)
      ELSE IF(form==2)THEN
        !
        dtmp = hpoly_bary(self%indsf(1,ind), cords(1))*cords(3)
        val = vec(:,2)*dtmp
        !
        dtmp = dhpoly_bary(self%indsf(1,ind), cords(1))*cords(3) &
             - hpoly_bary(self%indsf(1,ind), cords(1))
        val = val + vec(:,3)*dtmp*fhex(self%cmap(dof)%el)
      END IF
    case(4)
      vec(:,1)=-hcgop(:,3) ! (-3 x -2)
      vec(:,2)=hcgop(:,2) ! (-3 x 1)
      vec(:,3)=hcgop(:,1) ! (-2 x 1)
      ind=self%cmap(dof)%ind
      form=self%indsc(4,ind)
      IF(form<2)THEN
        vtmp=dhpoly_3d(self%indsc(1:3,ind), fhex(1:3))
        ! Part 1
        ! 1,2 6,2 1,4 6,4
        dtmp=(d2hpoly_3d(self%indsc(1:3,ind), 1, 2, fhex(1:3))*fhex(4)*fhex(6) &
            - vtmp(3)*fhex(4) - vtmp(2)*fhex(6) + vtmp(1))*fhex(5)
        f1 = -vec(:,1)*dtmp
        ! 1,3 6,3 1,5 6,5
        dtmp=(d2hpoly_3d(self%indsc(1:3,ind), 1, 3, fhex(1:3))*fhex(5)*fhex(6) &
            - vtmp(4)*fhex(5) - vtmp(2)*fhex(6) + vtmp(1))*fhex(4)
        f1 = f1 - vec(:,2)*dtmp
        ! Part 2
        ! 2,1 4,1 2,6 4,6
        dtmp=(d2hpoly_3d(self%indsc(1:3,ind), 2, 1, fhex(1:3))*fhex(4)*fhex(6) &
            - vtmp(2)*fhex(6) - vtmp(3)*fhex(4) + vtmp(1))*fhex(5)
        f2 = vec(:,1)*dtmp
        ! 2,3 4,3 2,5 4,5
        dtmp=(d2hpoly_3d(self%indsc(1:3,ind), 2, 3, fhex(1:3))*fhex(4)*fhex(5) &
            - vtmp(4)*fhex(5) - vtmp(3)*fhex(4) + vtmp(1))*fhex(6)
        f2 = f2 - vec(:,3)*dtmp
        ! Part 3
        ! 3,1 5,1 3,6 5,6
        dtmp=(d2hpoly_3d(self%indsc(1:3,ind), 3, 1, fhex(1:3))*fhex(5)*fhex(6) &
            - vtmp(2)*fhex(6) - vtmp(4)*fhex(5) + vtmp(1))*fhex(4)
        f3 = vec(:,2)*dtmp
        ! 3,2 5,2 3,4 5,4
        dtmp=(d2hpoly_3d(self%indsc(1:3,ind), 3, 2, fhex(1:3))*fhex(4)*fhex(5) &
            - vtmp(3)*fhex(4) - vtmp(4)*fhex(5) + vtmp(1))*fhex(6)
        f3 = f3 + vec(:,3)*dtmp
      END IF
      IF(form==0)THEN
        val = f1 - f2 + f3
      ELSE IF(form==1)THEN
        val = f1 - f2 - f3
      ELSE IF(form==2)THEN
        vtmp(1:3)=dhpoly_2d(self%indsc(2:3,ind), fhex(2:3))
        !
        dtmp = (vtmp(2)*fhex(4) - vtmp(1))*fhex(5)
        val = -vec(:,1)*dtmp
        !
        dtmp = (vtmp(3)*fhex(5) - vtmp(1))*fhex(4)
        val = val - vec(:,2)*dtmp
      ELSE IF(form==3)THEN
        vtmp(1:3)=dhpoly_2d(self%indsc((/1,3/),ind), fhex((/1,3/)))
        !
        dtmp = (vtmp(2)*fhex(6) - vtmp(1))*fhex(5)
        val = vec(:,1)*dtmp
        !
        dtmp = (vtmp(3)*fhex(5) - vtmp(1))*fhex(6)
        val = val -vec(:,3)*dtmp
      ELSE IF(form==4)THEN
        vtmp(1:3)=dhpoly_2d(self%indsc(1:2,ind), fhex(1:2))
        !
        dtmp = (vtmp(2)*fhex(6) - vtmp(1))*fhex(4)
        val = vec(:,2)*dtmp
        !
        dtmp = (vtmp(3)*fhex(4) - vtmp(1))*fhex(6)
        val = val + vec(:,3)*dtmp
      END IF
  case default
    call oft_abort('Invalid element geometry','oft_hcurl_ceval',__FILE__)
  end select
ELSE
  IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
  select case(self%cmap(dof)%type)
    case(2)
      etmp=oriented_edges(:,self%cmap(dof)%el)
      val = 2.d0*cross_product(gop(:,etmp(1)),gop(:,etmp(2)))
    case(3)
      ftmp=oriented_faces(:,self%cmap(dof)%el)
      call oft_hcurl_cevalf(self%order,ftmp,self%cmap(dof)%ind,f,gop,val)
    case(4)
      call oft_hcurl_cevalc(self%order,self%cmap(dof)%ind,f,gop,val)
    case default
      call oft_abort('Invalid element geometry','oft_hcurl_ceval',__FILE__)
  end select
END IF
DEBUG_STACK_POP
end subroutine oft_hcurl_ceval
!---------------------------------------------------------------------------
! SUBROUTINE: oft_bhcurl_ceval
!---------------------------------------------------------------------------
!> Evaluate Nedelec H1(Curl) curl function on the boundary
!!
!! @note Evaluation is performed in logical coordinates with the resulting
!! vector in, and curl with respect to, physical coordinates.
!!
!! @param[in] self Nedelec type for evaluation (bfem)
!! @param[in] face Face for evaluation
!! @param[in] dof Element to evaluate
!! @param[in] f Position on face in logical space [4]
!! @param[out] val Curl of nedelec element (dof) at point (f) [3]
!! @param[in] gop Face Jacobian matrix at point (f) [3,3]
!---------------------------------------------------------------------------
subroutine oft_bhcurl_ceval(self,face,dof,f,val,gop)
class(oft_bfem_type), intent(in) :: self
integer(i4), intent(in) :: face,dof
real(r8), intent(in) :: f(:)
real(r8), optional, intent(in) :: gop(:,:)
real(r8), intent(out) :: val(3)
real(r8) :: grads(3,4)
integer(i4) :: ed,etmp(2),fc,ftmp(3)
DEBUG_STACK_PUSH
grads=1.d0
if(present(gop))grads(:,1:3)=gop
select case(self%cmap(dof)%type)
  case(1)
    call oft_abort('Invalid element geometry (no point DOFs)','oft_bhcurl_ceval',__FILE__)
  case(2)
    ed=self%mesh%lce(self%cmap(dof)%el,face)
    etmp=self%mesh%cell_ed(:,self%cmap(dof)%el)
    call orient_list2(ed,etmp)
    call oft_hcurl_cevale(self%order,etmp,self%cmap(dof)%ind,f,grads,val)
  case(3)
    fc=self%mesh%lco(face)
    ftmp=(/1,2,3/)
    call orient_listn(fc,ftmp,3_i4)
    call oft_hcurl_cevalf(self%order,ftmp,self%cmap(dof)%ind,f,grads,val)
end select
DEBUG_STACK_POP
end subroutine oft_bhcurl_ceval
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_cevale
!---------------------------------------------------------------------------
!> Evaluate edge based curl functions
!---------------------------------------------------------------------------
subroutine oft_hcurl_cevale(order,ed,dof,f,grads,val)
integer(i4), intent(in) :: order
integer(i4), intent(in) :: ed(2),dof
real(r8), intent(in) :: f(:),grads(3,4)
real(r8), intent(out) :: val(3)
real(r8) :: g1xg2
DEBUG_STACK_PUSH
!---Linear DOF
g1xg2 =  2.d0
!---Set contribution
val = g1xg2*cross_product(grads(:,ed(1)),grads(:,ed(2)))
DEBUG_STACK_POP
end subroutine oft_hcurl_cevale
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_cevalf
!---------------------------------------------------------------------------
!> Evaluate face based curl functions
!---------------------------------------------------------------------------
subroutine oft_hcurl_cevalf(order,fc,dof,f,grads,val)
integer(i4), intent(in) :: order
integer(i4), intent(in) :: fc(3),dof
real(r8), intent(in) :: f(:),grads(3,4)
real(r8), intent(out) :: val(3)
real(r8) :: f1,f2,f3
real(r8) :: g1xg2,g1xg3,g2xg3
DEBUG_STACK_PUSH
!---Get local variables
f1=f(fc(1)); g1xg2=0.d0
f2=f(fc(2)); g1xg3=0.d0
f3=f(fc(3)); g2xg3=0.d0
!---
SELECT CASE(dof)
!---Quadratic DOF
CASE(1)
  g1xg2 = 0.d0
  g1xg3 = 4*f2
  g2xg3 = 4*f1
CASE(2)
  g1xg2 = -2*f3
  g1xg3 = -f2
  g2xg3 = f1
!---Cubic DOF
CASE(3)
  g1xg2 = 4*f3*(f1-f2)
  g1xg3 = 4*f2*(-f1-f2+2*f3)
  g2xg3 = 4*f1*(-f1-f2+2*f3)
CASE(4)
  g1xg2 = f3*(3*f1+3*f2-2*f3)
  g1xg3 = f2*(f1+f2-2*f3)
  g2xg3 = f1*(-f1-f2+2*f3)
CASE(5)
  g1xg2 = 0.d0
  g1xg3 = 4*f2*(2*f1-f2)
  g2xg3 = 4*f1*(f1-2*f2)
!---Quartic DOF
CASE(6)
  g1xg2 =  f3*(-8.d0*f1**2 + 16.d0*f1*f3 + 8.d0*f2**2 - 16.d0*f2*f3)
  g1xg3 =  f2*(4.d0*f1**2 + 8.d0*f1*f2 - 32.d0*f1*f3 + 4.d0*f2**2 - 32.d0*f2*f3 + 12.d0*f3**2)
  g2xg3 =  f1*(4.d0*f1**2 + 8.d0*f1*f2 - 32.d0*f1*f3 + 4.d0*f2**2 - 32.d0*f2*f3 + 12.d0*f3**2)
CASE(7)
  g1xg2 =  f3*(-4.d0*f1**2 - 8.d0*f1*f2 + 12.d0*f1*f3 - 4.d0*f2**2 + 12.d0*f2*f3 - 2.d0*f3**2)
  g1xg3 =  f2*(-1.d0*f1**2 - 2.d0*f1*f2 + 8.d0*f1*f3 - 1.d0*f2**2 + 8.d0*f2*f3 - 3.d0*f3**2)
  g2xg3 =  f1*(1.d0*f1**2 + 2.d0*f1*f2 - 8.d0*f1*f3 + 1.d0*f2**2 - 8.d0*f2*f3 + 3.d0*f3**2)
CASE(8)
  g1xg2 =  f3*(4.d0*f1**2 - 16.d0*f1*f2 + 4.d0*f2**2)
  g1xg3 =  f2*(-8.d0*f1**2 - 4.d0*f1*f2 + 16.d0*f1*f3 + 4.d0*f2**2 - 8.d0*f2*f3)
  g2xg3 =  f1*(-4.d0*f1**2 + 4.d0*f1*f2 + 8.d0*f1*f3 + 8.d0*f2**2 - 16.d0*f2*f3)
CASE(9)
  g1xg2 =  0.d0
  g1xg3 =  f2*(12.d0*f1**2 - 24.d0*f1*f2 + 4.d0*f2**2)
  g2xg3 =  f1*(4.d0*f1**2 - 24.d0*f1*f2 + 12.d0*f2**2)
END SELECT
!---Sum contributions
val = g1xg2*cross_product(grads(:,fc(1)),grads(:,fc(2))) &
    + g1xg3*cross_product(grads(:,fc(1)),grads(:,fc(3))) &
    + g2xg3*cross_product(grads(:,fc(2)),grads(:,fc(3)))
DEBUG_STACK_POP
end subroutine oft_hcurl_cevalf
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_cevalc
!---------------------------------------------------------------------------
!> Evaluate cell based curl functions
!---------------------------------------------------------------------------
subroutine oft_hcurl_cevalc(order,dof,f,grads,val)
integer(i4), intent(in) :: order
integer(i4), intent(in) :: dof
real(r8), intent(in) :: f(:),grads(3,4)
real(r8), intent(out) :: val(3)
real(r8) :: f1,f2,f3,f4
real(r8) :: g1xg2,g1xg3,g1xg4,g2xg3,g2xg4,g3xg4
DEBUG_STACK_PUSH
!---Get local variables
f1=f(1); g1xg2=0.d0; g1xg3=0.d0; g1xg4=0.d0
f2=f(2); g2xg3=0.d0; g2xg4=0.d0
f3=f(3); g3xg4=0.d0
f4=f(4)
!---
SELECT CASE(dof)
!---Cubic DOF
CASE(1)
  g1xg2 = 0.d0
  g1xg3 = 4*f2*f4
  g1xg4 = 0.d0
  g2xg3 = 4*f1*f4
  g2xg4 = 0.d0
  g3xg4 = -4*f1*f2
CASE(2)
  g1xg2 = 0.d0
  g1xg3 = 4*f2*f4
  g1xg4 = 4*f2*f3
  g2xg3 = 4*f1*f4
  g2xg4 = 4*f1*f3
  g3xg4 = 0.d0
CASE(3)
  g1xg2 = -2*f3*f4
  g1xg3 = -f2*f4
  g1xg4 = -f2*f3
  g2xg3 = f1*f4
  g2xg4 = f1*f3
  g3xg4 = 0.d0
!---Quartic DOF
CASE(4)
  g1xg2 =  0.d0
  g1xg3 =  f2*f4*(-8.d0*f1 - 4.d0*f2 - 4.d0*f3 + 4.d0*f4)
  g1xg4 =  0.d0
  g2xg3 =  f1*f4*(-4.d0*f1 - 8.d0*f2 - 4.d0*f3 + 4.d0*f4)
  g2xg4 =  0.d0
  g3xg4 =  f1*f2*(4.d0*f1 + 4.d0*f2 + 4.d0*f3 - 8.d0*f4)
CASE(5)
  g1xg2 =  f3*f4*(4.d0*f1 - 4.d0*f2)
  g1xg3 =  f2*f4*(-4.d0*f1 - 4.d0*f2 - 8.d0*f3 + 4.d0*f4)
  g1xg4 =  f2*f3*(-4.d0*f1 - 4.d0*f2 - 4.d0*f3 + 8.d0*f4)
  g2xg3 =  f1*f4*(-4.d0*f1 - 4.d0*f2 - 8.d0*f3 + 4.d0*f4)
  g2xg4 =  f1*f3*(-4.d0*f1 - 4.d0*f2 - 4.d0*f3 + 8.d0*f4)
  g3xg4 =  0.d0
CASE(6)
  g1xg2 =  f3*f4*(3.d0*f1 + 3.d0*f2 + 2.d0*f3 - 2.d0*f4)
  g1xg3 =  f2*f4*(f1 + f2 + 2.d0*f3 - f4)
  g1xg4 =  f2*f3*(f1 + f2 + f3 - 2*f4)
  g2xg3 =  f1*f4*(-f1 - f2 - 2.d0*f3 + f4)
  g2xg4 =  f1*f3*(-f1 - f2 - f3 + 2.d0*f4)
  g3xg4 =  0.d0
CASE(7)
  g1xg2 =  f3*f4*(4.d0*f1 - 4.d0*f2)
  g1xg3 =  f2*f4*(-4.d0*f1 - 4.d0*f2 + 8.d0*f3)
  g1xg4 =  4.d0*f1*f2*f3
  g2xg3 =  f1*f4*(-4.d0*f1 - 4.d0*f2 + 8.d0*f3)
  g2xg4 =  4.d0*f1*f2*f3
  g3xg4 =  f1*f2*(4.d0*f1 + 4.d0*f2 - 8.d0*f3)
CASE(8)
  g1xg2 =  f3*f4*(4.d0*f1 - 4.d0*f2)
  g1xg3 =  f2*f4*(-4.d0*f1 - 4.d0*f2 + 8.d0*f3)
  g1xg4 =  f2*f3*(-4.d0*f1 - 4.d0*f2 + 4.d0*f3)
  g2xg3 =  f1*f4*(-4.d0*f1 - 4.d0*f2 + 8.d0*f3)
  g2xg4 =  f1*f3*(-4.d0*f1 - 4.d0*f2 + 4.d0*f3)
  g3xg4 =  0.d0
CASE(9)
  g1xg2 =  f3*f4*(3.d0*f1 + 3.d0*f2 - 2.d0*f3)
  g1xg3 =  f2*f4*(f1 + f2 - 2.d0*f3)
  g1xg4 =  f2*f3*(f1 + f2 - f3)
  g2xg3 =  f1*f4*(-f1 - f2 + 2.d0*f3)
  g2xg4 =  f1*f3*(-f1 - f2 + f3)
  g3xg4 =  0.d0
CASE(10)
  g1xg2 =  0.d0
  g1xg3 =  f2*f4*(8.d0*f1 - 4.d0*f2)
  g1xg4 =  0.d0
  g2xg3 =  f1*f4*(4.d0*f1 - 8.d0*f2)
  g2xg4 =  0.d0
  g3xg4 =  f1*f2*(-4.d0*f1 + 4.d0*f2)
CASE(11)
  g1xg2 =  0.d0
  g1xg3 =  f2*f4*(8.d0*f1 - 4.d0*f2)
  g1xg4 =  f2*f3*(8.d0*f1 - 4.d0*f2)
  g2xg3 =  f1*f4*(4.d0*f1 - 8.d0*f2)
  g2xg4 =  f1*f3*(4.d0*f1 - 8.d0*f2)
  g3xg4 =  0.d0
END SELECT
!---Sum contributions
val = g1xg2*cross_product(grads(:,1),grads(:,2)) &
    + g1xg3*cross_product(grads(:,1),grads(:,3)) &
    + g1xg4*cross_product(grads(:,1),grads(:,4)) &
    + g2xg3*cross_product(grads(:,2),grads(:,3)) &
    + g2xg4*cross_product(grads(:,2),grads(:,4)) &
    + g3xg4*cross_product(grads(:,3),grads(:,4))
DEBUG_STACK_POP
end subroutine oft_hcurl_cevalc
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_get_cgops
!---------------------------------------------------------------------------
!> Evaluate all lagrange interpolation functions
!!
!! @note Evaluation is performed in logical coordinates.
!!
!! @param[in] self Lagrange type for evaluation
!! @param[in] cell Cell for evaluation
!! @param[in] f Position in cell in logical space
!! @param[out] rop Value of interpolation functions at point (f) [3,ncdofs]
!! @param[in] gop Cell Jacobian matrix at point (f) [3,4]
!---------------------------------------------------------------------------
subroutine oft_hcurl_get_cgops(gop,cgop)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: cgop(3,6)
cgop(:,1) = cross_product(gop(:,1),gop(:,2))
cgop(:,2) = cross_product(gop(:,1),gop(:,3))
cgop(:,3) = cross_product(gop(:,1),gop(:,4))
cgop(:,4) = cross_product(gop(:,2),gop(:,3))
cgop(:,5) = cross_product(gop(:,2),gop(:,4))
cgop(:,6) = cross_product(gop(:,3),gop(:,4))
end subroutine oft_hcurl_get_cgops
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_ceval_all
!---------------------------------------------------------------------------
!> Evaluate all lagrange interpolation functions
!!
!! @note Evaluation is performed in logical coordinates.
!!
!! @param[in] self Lagrange type for evaluation
!! @param[in] cell Cell for evaluation
!! @param[in] f Position in cell in logical space
!! @param[out] rop Value of interpolation functions at point (f) [3,ncdofs]
!! @param[in] gop Cell Jacobian matrix at point (f) [3,4]
!---------------------------------------------------------------------------
subroutine oft_hcurl_ceval_all(self,cell,f,rop,cgop)
class(oft_hcurl_fem), intent(in) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(4)
real(r8), intent(in) :: cgop(3,6)
real(r8), contiguous, intent(out) :: rop(:,:)
integer(i4) :: i,j,etmp(2),fhtmp(4),offset
real(r8) :: gop(3,4),jac,vec(3,3),hcgop(3,3),vtmp(4)
real(r8) :: val(3),cords(4),fhex(6),gbary(3,6),dtmp,f1(3),f2(3),f3(3)
DEBUG_STACK_PUSH
IF(hex_mesh)THEN
  fhex=hex_get_bary(f)
  gbary=hex_get_bary_gop(gop)
  hcgop(:,1)=cgop(:,1) ! (1 x 2)
  hcgop(:,2)=cgop(:,2) ! (1 x 3)
  hcgop(:,3)=cgop(:,4) ! (2 x 3)
  !---Edges
  DO i=1,12
    offset=(i-1)*self%gstruct(2)
    etmp=hex_bary_ecoords(:,i)
    call orient_list2(self%mesh%lce(i,cell), etmp)
    !
    cords(1:2)=fhex(etmp)
    cords(3:4)=fhex(hex_bary_efcoords(:,i))
    vec(:,1)=hex_get_bary_cgop(hcgop,hex_bary_efcoords(1,i),etmp(2))*cords(4) &
            +hex_get_bary_cgop(hcgop,hex_bary_efcoords(2,i),etmp(2))*cords(3)
    rop(:,offset+1)=vec(:,1)
  END DO
  IF(self%gstruct(3)>0)THEN
    !---Faces
    DO i=1,6
      offset=(i-1)*self%gstruct(3)+12*self%gstruct(2)
      fhtmp=hex_bary_fcoords(:,i)
      IF(self%mesh%lcfo(i,cell)<0)fhtmp=fhtmp((/2,3,4,1/))
      call orient_listn(self%mesh%lcfo(i,cell), fhtmp, 4_i4)
      !
      cords=fhex(fhtmp)
      vec(:,1)=hex_get_bary_cgop(hcgop,i,fhtmp(1))
      vec(:,2)=hex_get_bary_cgop(hcgop,i,fhtmp(2))
      vec(:,3)=hex_get_bary_cgop(hcgop,fhtmp(1),fhtmp(2))
      DO j=1,self%gstruct(3)
        IF(self%indsf(3,j)==0)THEN
          vtmp(1:3)=dhpoly_2d(self%indsf(1:2,j), cords(1:2))
          ! Form 1
          dtmp=(vtmp(2)*cords(3) - vtmp(1))*cords(4)
          val = vec(:,1)*dtmp
          !
          dtmp= vtmp(2)*cords(3) - vtmp(1) &
              - d2hpoly_2d(self%indsf(1:2,j), 1,2, cords(1:2))*cords(3)*cords(4) &
              + vtmp(3)*cords(4)
          val = val + vec(:,3)*dtmp*fhex(i)
          ! Part 2
          dtmp=(vtmp(3)*cords(4) - vtmp(1))*cords(3)
          val = val - vec(:,2)*dtmp
          !
          dtmp= vtmp(3)*cords(4) - vtmp(1) &
              - d2hpoly_2d(self%indsf(1:2,j), 1,2, cords(1:2))*cords(3)*cords(4) &
              + vtmp(2)*cords(3)
          val = val + vec(:,3)*dtmp*fhex(i)
        ELSE IF(self%indsf(3,j)==1)THEN
          vtmp(1)=hpoly_bary(self%indsf(2,j), cords(2))
          !
          val = vec(:,1)*vtmp(1)*cords(4)
          !
          dtmp = dhpoly_bary(self%indsf(2,j), cords(2))*cords(4) &
               - vtmp(1)
          val = val - vec(:,3)*dtmp*fhex(i)
        ELSE IF(self%indsf(3,j)==2)THEN
          vtmp(1)=hpoly_bary(self%indsf(1,j), cords(1))
          !
          val = vec(:,2)*vtmp(1)*cords(3)
          !
          dtmp = dhpoly_bary(self%indsf(1,j), cords(1))*cords(3) &
               - vtmp(1)
          val = val + vec(:,3)*dtmp*fhex(i)
        END IF
        !
        rop(:,offset+j)=val
      END DO
    END DO
  END IF
  IF(self%gstruct(4)>0)THEN
    !---Cell
    offset=6*self%gstruct(3)+12*self%gstruct(2)
    vec(:,1)=-hcgop(:,3) ! (-3 x -2)
    vec(:,2)=hcgop(:,2) ! (-3 x 1)
    vec(:,3)=hcgop(:,1) ! (-2 x 1)
    DO j=1,self%gstruct(4)
      IF(self%indsc(4,j)<2)THEN
        ! Cached function evaluations
        vtmp=dhpoly_3d(self%indsc(1:3,j), fhex(1:3))
        ! Part 1
        ! 1,2 6,2 1,4 6,4
        dtmp=(d2hpoly_3d(self%indsc(1:3,j), 1, 2, fhex(1:3))*fhex(4)*fhex(6) &
            - vtmp(3)*fhex(4) - vtmp(2)*fhex(6) + vtmp(1))*fhex(5)
        f1 = -vec(:,1)*dtmp
        ! 1,3 6,3 1,5 6,5
        dtmp=(d2hpoly_3d(self%indsc(1:3,j), 1, 3, fhex(1:3))*fhex(5)*fhex(6) &
            - vtmp(4)*fhex(5) - vtmp(2)*fhex(6) + vtmp(1))*fhex(4)
        f1 = f1 - vec(:,2)*dtmp
        ! Part 2
        ! 2,1 4,1 2,6 4,6
        dtmp=(d2hpoly_3d(self%indsc(1:3,j), 1, 2, fhex(1:3))*fhex(4)*fhex(6) &
            - vtmp(2)*fhex(6) - vtmp(3)*fhex(4) + vtmp(1))*fhex(5)
        f2 = vec(:,1)*dtmp
        ! 2,3 4,3 2,5 4,5
        dtmp=(d2hpoly_3d(self%indsc(1:3,j), 2, 3, fhex(1:3))*fhex(4)*fhex(5) &
            - vtmp(4)*fhex(5) - vtmp(3)*fhex(4) + vtmp(1))*fhex(6)
        f2 = f2 - vec(:,3)*dtmp
        ! Part 3
        ! 3,1 5,1 3,6 5,6
        dtmp=(d2hpoly_3d(self%indsc(1:3,j), 1, 3, fhex(1:3))*fhex(5)*fhex(6) &
            - vtmp(2)*fhex(6) - vtmp(4)*fhex(5) + vtmp(1))*fhex(4)
        f3 = vec(:,2)*dtmp
        ! 3,2 5,2 3,4 5,4
        dtmp=(d2hpoly_3d(self%indsc(1:3,j), 2, 3, fhex(1:3))*fhex(4)*fhex(5) &
            - vtmp(3)*fhex(4) - vtmp(4)*fhex(5) + vtmp(1))*fhex(6)
        f3 = f3 + vec(:,3)*dtmp
      END IF
      IF(self%indsc(4,j)==0)THEN
        val = f1 - f2 + f3
      ELSE IF(self%indsc(4,j)==1)THEN
        val = f1 - f2 - f3
      ELSE IF(self%indsc(4,j)==2)THEN
        vtmp(1:3)=dhpoly_2d(self%indsc(2:3,j), fhex(2:3))
        !
        dtmp = (vtmp(2)*fhex(4) - vtmp(1))*fhex(5)
        val = -vec(:,1)*dtmp
        !
        dtmp = (vtmp(3)*fhex(5) - vtmp(1))*fhex(4)
        val = val - vec(:,2)*dtmp
      ELSE IF(self%indsc(4,j)==3)THEN
        vtmp(1:3)=dhpoly_2d(self%indsc((/1,3/),j), fhex((/1,3/)))
        !
        dtmp = (vtmp(2)*fhex(6) - vtmp(1))*fhex(5)
        val = vec(:,1)*dtmp
        !
        dtmp = (vtmp(3)*fhex(5) - vtmp(1))*fhex(6)
        val = val -vec(:,3)*dtmp
      ELSE IF(self%indsc(4,j)==4)THEN
        vtmp(1:3)=dhpoly_2d(self%indsc(1:2,j), fhex(1:2))
        !
        dtmp = (vtmp(2)*fhex(6) - vtmp(1))*fhex(4)
        val = vec(:,2)*dtmp
        !
        dtmp = (vtmp(3)*fhex(4) - vtmp(1))*fhex(6)
        val = val + vec(:,3)*dtmp
      END IF
      rop(:,offset+j)=val
    END DO
  END IF
ELSE
  IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
  select case(self%order)
    case(1)
      DO i=1,6
        etmp = oriented_edges(:,i)
        rop(:,i) = 2.d0*cgop(:,ABS(cgop_map(etmp(1),etmp(2))))*SIGN(1,cgop_map(etmp(1),etmp(2)))
      END DO
    case(2)
      call oft_hcurl_ceval_all2(self,cell,f,rop,cgop)
    case(3)
      call oft_hcurl_ceval_all3(self,cell,f,rop,cgop)
    case(4)
      call oft_hcurl_ceval_all4(self,cell,f,rop,cgop)
    case default ! Fall back to normal evaluation
      CALL oft_abort('BAD ORDER','oft_hcurl_ceval_all',__FILE__)
  end select
END IF
DEBUG_STACK_POP
end subroutine oft_hcurl_ceval_all
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_ceval_all2
!---------------------------------------------------------------------------
!> Evaluate all lagrange interpolation functions (quadratic)
!!
!! @note Evaluation is performed in logical coordinates.
!!
!! @param[in] self Lagrange type for evaluation
!! @param[in] cell Cell for evaluation
!! @param[in] f Position in cell in logical space
!! @param[out] rop Value of interpolation functions at point (f) [ncdofs]
!---------------------------------------------------------------------------
subroutine oft_hcurl_ceval_all2(self,cell,f,rop,cgop)
class(oft_fem_type), intent(in) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(4)
real(r8), intent(in) :: cgop(3,6)
real(r8), intent(out) :: rop(3,14)
integer(i4) :: i,etmp(2),ftmp(3)
real(r8) :: f1,f2,f3
real(r8) :: g1xg2(3),g1xg3(3),g2xg3(3)
DEBUG_STACK_PUSH
DO i=1,6
  etmp=oriented_edges(:,i)
  rop(:,i) = 2.d0*cgop(:,ABS(cgop_map(etmp(1),etmp(2))))*SIGN(1,cgop_map(etmp(1),etmp(2)))
END DO
DO i=1,4
  ftmp=oriented_faces(:,i)
  f1 = f(ftmp(1)); f2 = f(ftmp(2)); f3 = f(ftmp(3))
  g1xg2 = cgop(:,ABS(cgop_map(ftmp(1),ftmp(2))))*SIGN(1,cgop_map(ftmp(1),ftmp(2)))
  g1xg3 = cgop(:,ABS(cgop_map(ftmp(1),ftmp(3))))*SIGN(1,cgop_map(ftmp(1),ftmp(3)))
  g2xg3 = cgop(:,ABS(cgop_map(ftmp(2),ftmp(3))))*SIGN(1,cgop_map(ftmp(2),ftmp(3)))
  !---Quadratic DOF
  rop(:,(i-1)*2+7) = 0*g1xg2 &
             + 4*f2*g1xg3 &
             + 4*f1*g2xg3
  !
  rop(:,(i-1)*2+8) = -2*f3*g1xg2 &
              - f2*g1xg3 &
              + f1*g2xg3
END DO
DEBUG_STACK_POP
end subroutine oft_hcurl_ceval_all2
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_ceval_all3
!---------------------------------------------------------------------------
!> Evaluate all lagrange interpolation functions (cubic)
!!
!! @note Evaluation is performed in logical coordinates.
!!
!! @param[in] self Lagrange type for evaluation
!! @param[in] cell Cell for evaluation
!! @param[in] f Position in cell in logical space
!! @param[out] rop Value of interpolation functions at point (f) [ncdofs]
!---------------------------------------------------------------------------
subroutine oft_hcurl_ceval_all3(self,cell,f,rop,cgop)
class(oft_fem_type), intent(in) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(4)
real(r8), intent(in) :: cgop(3,6)
real(r8), intent(out) :: rop(3,29)
integer(i4) :: i,etmp(2),ftmp(3)
real(r8) :: f1,f2,f3,f4
real(r8) :: g1xg2(3),g1xg3(3),g2xg3(3)
DEBUG_STACK_PUSH
DO i=1,6
  etmp=oriented_edges(:,i)
  rop(:,i) = 2.d0*cgop(:,ABS(cgop_map(etmp(1),etmp(2))))*SIGN(1,cgop_map(etmp(1),etmp(2)))
END DO
DO i=1,4
  ftmp=oriented_faces(:,i)
  f1 = f(ftmp(1)); f2 = f(ftmp(2)); f3 = f(ftmp(3))
  g1xg2 = cgop(:,ABS(cgop_map(ftmp(1),ftmp(2))))*SIGN(1,cgop_map(ftmp(1),ftmp(2)))
  g1xg3 = cgop(:,ABS(cgop_map(ftmp(1),ftmp(3))))*SIGN(1,cgop_map(ftmp(1),ftmp(3)))
  g2xg3 = cgop(:,ABS(cgop_map(ftmp(2),ftmp(3))))*SIGN(1,cgop_map(ftmp(2),ftmp(3)))
  !---Quadratic DOF
  rop(:,(i-1)*5+7) = 0*g1xg2 &
              + 4*f2*g1xg3 &
              + 4*f1*g2xg3
  !
  rop(:,(i-1)*5+8) = -2*f3*g1xg2 &
              - f2*g1xg3 &
              + f1*g2xg3
  !---Cubic DOF
  rop(:,(i-1)*5+9) = 4*f3*(f1-f2)*g1xg2 &
              + 4*f2*(-f1-f2+2*f3)*g1xg3 &
              + 4*f1*(-f1-f2+2*f3)*g2xg3
  !
  rop(:,(i-1)*5+10) = f3*(3*f1+3*f2-2*f3)*g1xg2 &
              + f2*(f1+f2-2*f3)*g1xg3 &
              + f1*(-f1-f2+2*f3)*g2xg3
  !
  rop(:,(i-1)*5+11) = 0*g1xg2 &
              + 4*f2*(2*f1-f2)*g1xg3 &
              + 4*f1*(f1-2*f2)*g2xg3
END DO
f1 = f(1); f2 = f(2); f3 = f(3); f4 = f(4)
!---Cubic DOF
rop(:,27) = 0*cgop(:,1) &
          + 4*f2*f4*cgop(:,2) &
          + 0*cgop(:,3) &
          + 4*f1*f4*cgop(:,4) &
          + 0*cgop(:,5) &
          - 4*f1*f2*cgop(:,6)
!
rop(:,28) = 0*cgop(:,1) &
          + 4*f2*f4*cgop(:,2) &
          + 4*f2*f3*cgop(:,3) &
          + 4*f1*f4*cgop(:,4) &
          + 4*f1*f3*cgop(:,5) &
          + 0*cgop(:,6)
!
rop(:,29) = -2*f3*f4*cgop(:,1) &
          - f2*f4*cgop(:,2) &
          - f2*f3*cgop(:,3) &
          + f1*f4*cgop(:,4) &
          + f1*f3*cgop(:,5) &
          + 0*cgop(:,6)
DEBUG_STACK_POP
end subroutine oft_hcurl_ceval_all3
!---------------------------------------------------------------------------
! SUBROUTINE: oft_hcurl_ceval_all4
!---------------------------------------------------------------------------
!> Evaluate all lagrange interpolation functions (quartic)
!!
!! @note Evaluation is performed in logical coordinates.
!!
!! @param[in] self Lagrange type for evaluation
!! @param[in] cell Cell for evaluation
!! @param[in] f Position in cell in logical space
!! @param[out] rop Value of interpolation functions at point (f) [ncdofs]
!---------------------------------------------------------------------------
subroutine oft_hcurl_ceval_all4(self,cell,f,rop,cgop)
class(oft_fem_type), intent(in) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(4)
real(r8), intent(in) :: cgop(3,6)
real(r8), intent(out) :: rop(3,53)
integer(i4) :: i,etmp(2),ftmp(3)
real(r8) :: f1,f2,f3,f4
real(r8) :: g1xg2(3),g1xg3(3),g2xg3(3)
DEBUG_STACK_PUSH
DO i=1,6
  etmp=oriented_edges(:,i)
  rop(:,i) = 2.d0*cgop(:,ABS(cgop_map(etmp(1),etmp(2))))*SIGN(1,cgop_map(etmp(1),etmp(2)))
END DO
DO i=1,4
  ftmp=oriented_faces(:,i)
  f1 = f(ftmp(1)); f2 = f(ftmp(2)); f3 = f(ftmp(3))
  g1xg2 = cgop(:,ABS(cgop_map(ftmp(1),ftmp(2))))*SIGN(1,cgop_map(ftmp(1),ftmp(2)))
  g1xg3 = cgop(:,ABS(cgop_map(ftmp(1),ftmp(3))))*SIGN(1,cgop_map(ftmp(1),ftmp(3)))
  g2xg3 = cgop(:,ABS(cgop_map(ftmp(2),ftmp(3))))*SIGN(1,cgop_map(ftmp(2),ftmp(3)))
  !---Quadratic DOF
  rop(:,(i-1)*9+7) = 0*g1xg2 &
              + 4*f2*g1xg3 &
              + 4*f1*g2xg3
  !
  rop(:,(i-1)*9+8) = -2*f3*g1xg2 &
              - f2*g1xg3 &
              + f1*g2xg3
  !---Cubic DOF
  rop(:,(i-1)*9+9) = 4*f3*(f1-f2)*g1xg2 &
              + 4*f2*(-f1-f2+2*f3)*g1xg3 &
              + 4*f1*(-f1-f2+2*f3)*g2xg3
  !
  rop(:,(i-1)*9+10) = f3*(3*f1+3*f2-2*f3)*g1xg2 &
              + f2*(f1+f2-2*f3)*g1xg3 &
              + f1*(-f1-f2+2*f3)*g2xg3
  !
  rop(:,(i-1)*9+11) = 0*g1xg2 &
              + 4*f2*(2*f1-f2)*g1xg3 &
              + 4*f1*(f1-2*f2)*g2xg3
  !---Quartic DOF
  rop(:,(i-1)*9+12) = f3*(-8.d0*f1**2 + 16.d0*f1*f3 + 8.d0*f2**2 - 16.d0*f2*f3)*g1xg2 &
              + f2*(4.d0*f1**2 + 8.d0*f1*f2 - 32.d0*f1*f3 + 4.d0*f2**2 - 32.d0*f2*f3 + 12.d0*f3**2)*g1xg3 &
              + f1*(4.d0*f1**2 + 8.d0*f1*f2 - 32.d0*f1*f3 + 4.d0*f2**2 - 32.d0*f2*f3 + 12.d0*f3**2)*g2xg3
  !
  rop(:,(i-1)*9+13) = f3*(-4.d0*f1**2 - 8.d0*f1*f2 + 12.d0*f1*f3 - 4.d0*f2**2 + 12.d0*f2*f3 - 2.d0*f3**2)*g1xg2 &
              + f2*(-1.d0*f1**2 - 2.d0*f1*f2 + 8.d0*f1*f3 - 1.d0*f2**2 + 8.d0*f2*f3 - 3.d0*f3**2)*g1xg3 &
              + f1*(1.d0*f1**2 + 2.d0*f1*f2 - 8.d0*f1*f3 + 1.d0*f2**2 - 8.d0*f2*f3 + 3.d0*f3**2)*g2xg3
  !
  rop(:,(i-1)*9+14) = f3*(4.d0*f1**2 - 16.d0*f1*f2 + 4.d0*f2**2)*g1xg2 &
              + f2*(-8.d0*f1**2 - 4.d0*f1*f2 + 16.d0*f1*f3 + 4.d0*f2**2 - 8.d0*f2*f3)*g1xg3 &
              + f1*(-4.d0*f1**2 + 4.d0*f1*f2 + 8.d0*f1*f3 + 8.d0*f2**2 - 16.d0*f2*f3)*g2xg3
  !
  rop(:,(i-1)*9+15) = 0*g1xg2 &
              + f2*(12.d0*f1**2 - 24.d0*f1*f2 + 4.d0*f2**2)*g1xg3 &
              + f1*(4.d0*f1**2 - 24.d0*f1*f2 + 12.d0*f2**2)*g2xg3
END DO
f1 = f(1); f2 = f(2); f3 = f(3); f4 = f(4)
!---Cubic DOF
rop(:,43) = 0*cgop(:,1) &
          + 4*f2*f4*cgop(:,2) &
          + 0*cgop(:,3) &
          + 4*f1*f4*cgop(:,4) &
          + 0*cgop(:,5) &
          - 4*f1*f2*cgop(:,6)
!
rop(:,44) = 0*cgop(:,1) &
          + 4*f2*f4*cgop(:,2) &
          + 4*f2*f3*cgop(:,3) &
          + 4*f1*f4*cgop(:,4) &
          + 4*f1*f3*cgop(:,5) &
          + 0*cgop(:,6)
!
rop(:,45) = -2*f3*f4*cgop(:,1) &
          - f2*f4*cgop(:,2) &
          - f2*f3*cgop(:,3) &
          + f1*f4*cgop(:,4) &
          + f1*f3*cgop(:,5) &
          + 0*cgop(:,6)
!---Quartic DOF
rop(:,46) = 0*cgop(:,1) &
          + f2*f4*(-8.d0*f1 - 4.d0*f2 - 4.d0*f3 + 4.d0*f4)*cgop(:,2) &
          + 0*cgop(:,3) &
          + f1*f4*(-4.d0*f1 - 8.d0*f2 - 4.d0*f3 + 4.d0*f4)*cgop(:,4) &
          + 0*cgop(:,5) &
          + f1*f2*(4.d0*f1 + 4.d0*f2 + 4.d0*f3 - 8.d0*f4)*cgop(:,6)
!
rop(:,47) = f3*f4*(4.d0*f1 - 4.d0*f2)*cgop(:,1) &
          + f2*f4*(-4.d0*f1 - 4.d0*f2 - 8.d0*f3 + 4.d0*f4)*cgop(:,2) &
          + f2*f3*(-4.d0*f1 - 4.d0*f2 - 4.d0*f3 + 8.d0*f4)*cgop(:,3) &
          + f1*f4*(-4.d0*f1 - 4.d0*f2 - 8.d0*f3 + 4.d0*f4)*cgop(:,4) &
          + f1*f3*(-4.d0*f1 - 4.d0*f2 - 4.d0*f3 + 8.d0*f4)*cgop(:,5) &
          + 0*cgop(:,6)
!
rop(:,48) = f3*f4*(3.d0*f1 + 3.d0*f2 + 2.d0*f3 - 2.d0*f4)*cgop(:,1) &
          + f2*f4*(f1 + f2 + 2.d0*f3 - f4)*cgop(:,2) &
          + f2*f3*(f1 + f2 + f3 - 2*f4)*cgop(:,3) &
          + f1*f4*(-f1 - f2 - 2.d0*f3 + f4)*cgop(:,4) &
          + f1*f3*(-f1 - f2 - f3 + 2.d0*f4)*cgop(:,5) &
          + 0*cgop(:,6)
!
rop(:,49) = f3*f4*(4.d0*f1 - 4.d0*f2)*cgop(:,1) &
          + f2*f4*(-4.d0*f1 - 4.d0*f2 + 8.d0*f3)*cgop(:,2) &
          + 4.d0*f1*f2*f3*cgop(:,3) &
          + f1*f4*(-4.d0*f1 - 4.d0*f2 + 8.d0*f3)*cgop(:,4) &
          + 4.d0*f1*f2*f3*cgop(:,5) &
          + f1*f2*(4.d0*f1 + 4.d0*f2 - 8.d0*f3)*cgop(:,6)
!
rop(:,50) = f3*f4*(4.d0*f1 - 4.d0*f2)*cgop(:,1) &
          + f2*f4*(-4.d0*f1 - 4.d0*f2 + 8.d0*f3)*cgop(:,2) &
          + f2*f3*(-4.d0*f1 - 4.d0*f2 + 4.d0*f3)*cgop(:,3) &
          + f1*f4*(-4.d0*f1 - 4.d0*f2 + 8.d0*f3)*cgop(:,4) &
          + f1*f3*(-4.d0*f1 - 4.d0*f2 + 4.d0*f3)*cgop(:,5) &
          + 0*cgop(:,6)
!
rop(:,51) = f3*f4*(3.d0*f1 + 3.d0*f2 - 2.d0*f3)*cgop(:,1) &
          + f2*f4*(f1 + f2 - 2.d0*f3)*cgop(:,2) &
          + f2*f3*(f1 + f2 - f3)*cgop(:,3) &
          + f1*f4*(-f1 - f2 + 2.d0*f3)*cgop(:,4) &
          + f1*f3*(-f1 - f2 + f3)*cgop(:,5) &
          + 0*cgop(:,6)
!
rop(:,52) = 0*cgop(:,1) &
          + f2*f4*(8.d0*f1 - 4.d0*f2)*cgop(:,2) &
          + 0*cgop(:,3) &
          + f1*f4*(4.d0*f1 - 8.d0*f2)*cgop(:,4) &
          + 0*cgop(:,5) &
          + f1*f2*(-4.d0*f1 + 4.d0*f2)*cgop(:,6)
!
rop(:,53) = 0*cgop(:,1) &
          + f2*f4*(8.d0*f1 - 4.d0*f2)*cgop(:,2) &
          + f2*f3*(8.d0*f1 - 4.d0*f2)*cgop(:,3) &
          + f1*f4*(4.d0*f1 - 8.d0*f2)*cgop(:,4) &
          + f1*f3*(4.d0*f1 - 8.d0*f2)*cgop(:,5) &
          + 0*cgop(:,6)
DEBUG_STACK_POP
end subroutine oft_hcurl_ceval_all4
end module oft_hcurl_basis
