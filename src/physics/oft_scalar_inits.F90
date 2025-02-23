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
use multigrid, only: mg_mesh
use fem_utils, only: fem_interp
implicit none
!---------------------------------------------------------------------------
! CLASS oft_scalar_torus
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
! SUBROUTINE torus_interp
!---------------------------------------------------------------------------
!> Evaluate torus source
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [1]
!---------------------------------------------------------------------------
subroutine torus_interp(self,cell,f,gop,val)
class(oft_scalar_torus), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: rop(1),pt(3),r,z,d
!---Get coordinates
pt=mg_mesh%mesh%log2phys(cell,f)
r=sqrt(pt(1)**2+pt(2)**2)
z=pt(3)
!---Evaluate function
d=sqrt((r-self%r0)**2 + (z-self%z)**2)
val(1)=exp(-(d**2)/self%lam)
end subroutine torus_interp
end module oft_scalar_inits
