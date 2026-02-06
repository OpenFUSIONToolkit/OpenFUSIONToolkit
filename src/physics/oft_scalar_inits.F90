!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_scalar_inits.F90
!
!> Field initializations and evaluation for common scalar analytic field types
!!
!! @authors Chris Hansen
!! @date September 2012
!! @ingroup doxy_oft_physics
!---------------------------------------------------------------------------
module oft_scalar_inits
use oft_base
use oft_mesh_type, only: oft_mesh
use fem_utils, only: fem_interp, bfem_interp
implicit none
!---------------------------------------------------------------------------
!> Interpolation class for a uniform vector field
!---------------------------------------------------------------------------
type, extends(fem_interp) :: poss_scalar_field
  procedure(poss_scalar_eval), pointer, nopass :: func => NULL()
contains
  !> Reconstruct magnetic field
  procedure :: interp => poss_scalar_interp
end type poss_scalar_field
!---------------------------------------------------------------------------
!> Interpolation class for a uniform vector field
!---------------------------------------------------------------------------
type, extends(bfem_interp) :: poss_scalar_bfield
  procedure(poss_scalar_eval), pointer, nopass :: func => NULL()
contains
  !> Reconstruct magnetic field
  procedure :: interp => poss_scalar_binterp
end type poss_scalar_bfield
!
INTERFACE
  SUBROUTINE poss_scalar_eval(pt,val)
  IMPORT i4, r8
  REAL(r8), INTENT(in) :: pt(3)
  REAL(r8), INTENT(out) :: val
  END SUBROUTINE poss_scalar_eval
END INTERFACE
!---------------------------------------------------------------------------
!> Interpolation class for an axisymmetric gaussian source in toroidal geometry
!!
!! In toroidal coordinates defined by the class the scalar field is defined
!! as
!! \f[ S = e^{-\frac{r^2}{\lambda}} \f]
!---------------------------------------------------------------------------
type, extends(fem_interp) :: oft_scalar_torus
  real(r8) :: r0 = 1.d0 !< Major radius
  real(r8) :: z = 0.d0 !< Z-axis offset
  real(r8) :: lam = 1.d0 !< Source width parameter
contains
  !> Reconstruct field
  procedure :: interp => torus_interp
end type oft_scalar_torus
contains
!---------------------------------------------------------------------------
!> Return a uniform vector field
!---------------------------------------------------------------------------
subroutine poss_scalar_interp(self,cell,f,gop,val)
class(poss_scalar_field), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
real(r8) :: pt(3)
IF(.NOT.ASSOCIATED(self%func))CALL oft_abort("No eval function specified", &
  "poss_scalar_interp", __FILE__)
pt=self%mesh%log2phys(cell,f)
CALL self%func(pt,val(1))
end subroutine poss_scalar_interp
!---------------------------------------------------------------------------
!> Return a uniform vector field
!---------------------------------------------------------------------------
subroutine poss_scalar_binterp(self,cell,f,gop,val)
class(poss_scalar_bfield), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,3)
real(r8), intent(out) :: val(:)
real(r8) :: pt(3)
IF(.NOT.ASSOCIATED(self%func))CALL oft_abort("No eval function specified", &
  "poss_scalar_binterp", __FILE__)
pt=self%mesh%log2phys(cell,f)
CALL self%func(pt,val(1))
end subroutine poss_scalar_binterp
!---------------------------------------------------------------------------
! SUBROUTINE torus_interp
!---------------------------------------------------------------------------
!> Evaluate torus source
!---------------------------------------------------------------------------
subroutine torus_interp(self,cell,f,gop,val)
class(oft_scalar_torus), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [4]
real(r8), intent(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
real(r8), intent(out) :: val(:) !< Reconstructed field at f [1]
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: rop(1),pt(3),r,z,d
!---Get coordinates
pt=self%mesh%log2phys(cell,f)
r=sqrt(pt(1)**2+pt(2)**2)
z=pt(3)
!---Evaluate function
d=sqrt((r-self%r0)**2 + (z-self%z)**2)
val(1)=exp(-(d**2)/self%lam)
end subroutine torus_interp
end module oft_scalar_inits
