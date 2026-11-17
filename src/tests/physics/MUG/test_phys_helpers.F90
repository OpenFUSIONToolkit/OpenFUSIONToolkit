!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_alfven.F90
!
!> Module for initializing a traveling alfven wave
!------------------------------------------------------------------------------
MODULE test_phys_helpers
USE oft_base
USE oft_mesh_type, ONLY: oft_mesh
USE oft_la_base, ONLY: oft_vector
USE fem_base, ONLY: oft_afem_type
USE fem_utils, ONLY: fem_interp
USE oft_lag_operators, ONLY: oft_lag_vproject
USE diagnostic, ONLY: vec_energy
USE mhd_utils, ONLY: mu0
IMPLICIT NONE
!------------------------------------------------------------------------------
!> Interpolation class for sound wave initialization
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Interpolation class for alfven wave initialization
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Interpolate the desired component of a traveling sound wave
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Interpolate the vector perturbation of a traveling alfven wave
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Project a vector Lagrange field onto sine and cosine plane-wave modes
!!
!! For each Cartesian direction e_i, computes the L2 projection coefficient of
!! the field onto e_i sin(2 pi k.r/lam) and e_i cos(2 pi k.r/lam)
!------------------------------------------------------------------------------
subroutine lag_vec_mode_proj(fe_rep,mesh,field,k_dir,lam,quad_order,proj)
class(oft_afem_type), target, intent(inout) :: fe_rep !< Scalar Lagrange FE representation
class(oft_mesh), target, intent(inout) :: mesh !< Mesh
class(oft_vector), intent(inout) :: field !< Field to project
real(r8), intent(in) :: k_dir(3) !< Direction of wave propogation
real(r8), intent(in) :: lam !< Wavelength
integer(i4), intent(in) :: quad_order !< Quadrature order for mode normalization
real(r8), intent(out) :: proj(2,3) !< Projections [sin/cos,component]
integer(i4) :: i,j
class(oft_vector), pointer :: rhs
type(alfven_eig) :: mode
!---
CALL field%new(rhs)
mode%mesh=>mesh
mode%k_dir=k_dir
mode%lam=lam
DO i=1,3
  mode%v_dir=0.d0
  mode%v_dir(i)=1.d0
  DO j=1,2
    mode%r0=0.d0
    IF(j==2)mode%r0=-k_dir*lam/4.d0 ! Quarter wavelength shift gives cosine
    CALL oft_lag_vproject(fe_rep,mode,rhs)
    proj(j,i)=field%dot(rhs)/vec_energy(mesh,mode,quad_order)
  END DO
END DO
CALL rhs%delete
DEALLOCATE(rhs)
end subroutine lag_vec_mode_proj
END MODULE test_phys_helpers
