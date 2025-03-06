!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_alfven.F90
!
!> Module for initializing a traveling alfven wave
!---------------------------------------------------------------------------
MODULE test_phys_helpers
USE oft_base
USE oft_mesh_type, ONLY: oft_mesh
USE fem_utils, ONLY: fem_interp
USE mhd_utils, ONLY: mu0
IMPLICIT NONE
!---------------------------------------------------------------------------
!> Interpolation class for sound wave initialization
!---------------------------------------------------------------------------
type, extends(fem_interp) :: sound_eig
  character(LEN=1) :: field = 'n' !< Field component to initialize
  logical :: diff = .FALSE. !< Compute deviation from equilibrium
  real(r8) :: k_dir(3) = (/0.d0,0.d0,1.d0/) !< Direction of wave propogation
  real(r8) :: r0(3) = (/0.d0,0.d0,0.d0/)  !< Zero-phase position
  real(r8) :: lam = 2.d0 !< Wavelength
  real(r8) :: delta = 1.d-3 !< Relative size of perturbation (<<1)
contains
  !> Reconstruct sound wave fields
  procedure :: interp => sound_eig_interp
end type sound_eig
!---------------------------------------------------------------------------
!> Interpolation class for alfven wave initialization
!---------------------------------------------------------------------------
type, extends(fem_interp) :: alfven_eig
  real(r8) :: v_dir(3) = (/1.d0,0.d0,0.d0/) !< Direction of velocity perturbation
  real(r8) :: k_dir(3) = (/0.d0,0.d0,1.d0/) !< Direction of wave propogation
  real(r8) :: r0(3) = (/0.d0,0.d0,0.d0/)  !< Zero-phase position
  real(r8) :: lam = 2.d0 !< Wavelength
contains
  !> Reconstruct alfven wave fields
  procedure :: interp => alfven_eig_interp
end type alfven_eig
CONTAINS
!---------------------------------------------------------------------------
!> Interpolate the desired component of a traveling sound wave
!---------------------------------------------------------------------------
subroutine sound_eig_interp(self,cell,f,gop,val)
class(sound_eig), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
real(r8) :: pt(3),i,ar,az,s,c
!---
pt=self%mesh%log2phys(cell,f)
IF(self%field=='n')THEN
  val=(1.d0+self%delta*SIN(DOT_PRODUCT(pt-self%r0,self%k_dir)*2.d0*pi/self%lam))**(3.d0/5.d0)
  IF(self%diff)val=val-1.d0
ELSE IF(self%field=='t')THEN
  val=(1.d0+self%delta*SIN(DOT_PRODUCT(pt-self%r0,self%k_dir)*2.d0*pi/self%lam))**(2.d0/5.d0)
  IF(self%diff)val=val-1.d0
ELSE IF(self%field=='v')THEN
  val=self%delta*self%k_dir*SIN(DOT_PRODUCT(pt-self%r0,self%k_dir)*2.d0*pi/self%lam)
ELSE
  CALL oft_abort('Unknown field component','sound_eig_interp',__FILE__)
END IF
end subroutine sound_eig_interp
!---------------------------------------------------------------------------
!> Interpolate the vector perturbation of a traveling alfven wave
!---------------------------------------------------------------------------
subroutine alfven_eig_interp(self,cell,f,gop,val)
class(alfven_eig), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
real(r8) :: pt(3),i,ar,az,s,c
!---
pt=self%mesh%log2phys(cell,f)
val=self%v_dir*SIN(DOT_PRODUCT(pt-self%r0,self%k_dir)*2.d0*pi/self%lam)
end subroutine alfven_eig_interp
END MODULE test_phys_helpers