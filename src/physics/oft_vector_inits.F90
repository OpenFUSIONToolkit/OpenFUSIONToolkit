!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_vector_inits.F90
!
!> Field initializations and evaluation for common vector analytic field types
!!
!! @authors Chris Hansen
!! @date September 2012
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
module oft_vector_inits
use oft_base
use oft_mesh_type, only: oft_mesh
use fem_utils, only: fem_interp
use mhd_utils, only: mu0
implicit none
!------------------------------------------------------------------------------
!> Interpolation class for a uniform vector field
!------------------------------------------------------------------------------
type, extends(fem_interp) :: uniform_field
  integer(i4) :: n=3
  real(r8) :: val(3) = (/1.d0,0.d0,0.d0/) !< Field to initialize
contains
  !> Reconstruct magnetic field
  procedure :: interp => uniform_field_interp
end type uniform_field
!------------------------------------------------------------------------------
!> Interpolation class for a uniform vector field
!------------------------------------------------------------------------------
type, extends(fem_interp) :: poss_vec_field
  integer(i4) :: n=3
  procedure(poss_vec_eval), pointer, nopass :: func => NULL()
contains
  !> Reconstruct magnetic field
  procedure :: interp => poss_vec_interp
end type poss_vec_field
!------------------------------------------------------------------------------
!> Evaluate analytic fields for the tuna can spheromak
!------------------------------------------------------------------------------
type, extends(fem_interp) :: cyl_taylor
  real(r8) :: scale = 1.d0 !< Global scale of field
  real(r8) :: zmin = 1.d0 !< Lowest z-location on mesh
  real(r8) :: zmax = 0.d0 !< Highest z-location on mesh
  real(r8) :: rmax = 1.d0 !< Maximum r-location on mesh
  integer(i4) :: rplane(2) = (/1,2/) !< Poloidal plane components
  integer(i4) :: zaxis = 3 !< Toroidal axis component
  real(r8), private :: akr = 0.d0 !< Radial scale value
  real(r8), private :: akz = 0.d0 !< Axial scale value
  real(r8), private :: alm = 0.d0 !< Global scale value
contains
  !> Setup reconstruction for current mesh
  procedure :: setup => cyl_taylor_setup
  !> Reconstruct magnetic field
  procedure :: interp => cyl_taylor_interp
end type cyl_taylor
!------------------------------------------------------------------------------
!> Evaluate analytic fields for a set of straight infinite coils
!------------------------------------------------------------------------------
type, extends(fem_interp) :: inf_coils
  integer(i4) :: ncoils = 0 !< Number of coils in set
  real(r8), pointer, dimension(:,:) :: axis => NULL() !< Direction of extent
  real(r8), pointer, dimension(:,:) :: center => NULL() !< Possition on coil
  real(r8), pointer, dimension(:) :: current => NULL() !< Current in coil
contains
  !> Setup coil interpolation class
  procedure :: setup => inf_coils_setup
  !> Evalute cummulative field from all coils
  procedure :: interp => inf_coils_interp
end type inf_coils
!------------------------------------------------------------------------------
!> Field corresponding to a poloidal circulation in toroidal corrdinates
!------------------------------------------------------------------------------
type, extends(fem_interp) :: tor_radial
  integer(i4) :: n = 1 !< Toroidal mode number
  integer(i4) :: m = 0 !< Poloidal mode number
  real(r8) :: phase0 = 0.d0 !< Mode phase
  real(r8) :: R0 = 1.d0 !< Radial center for toroidal coordinates
  real(r8) :: Z0 = 0.d0 !< Vertical center for toroidal coordinates
  real(r8) :: rc = 1.d0 !< Center of circulation layer
  real(r8) :: rwidth = 1.d0 !< Width of circulation layer
contains
  !> Reconstruct field
  procedure :: interp => tor_radial_interp
end type tor_radial
!
INTERFACE
  SUBROUTINE poss_vec_eval(pt,val,n)
  IMPORT i4, r8
  REAL(r8), INTENT(in) :: pt(3)
  REAL(r8), INTENT(out) :: val(n)
  INTEGER(i4), INTENT(in) :: n
  END SUBROUTINE poss_vec_eval
END INTERFACE
contains
!------------------------------------------------------------------------------
!> Return a uniform vector field
!------------------------------------------------------------------------------
subroutine uniform_field_interp(self,cell,f,gop,val)
class(uniform_field), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
val(1:self%n)=self%val(1:self%n)
end subroutine uniform_field_interp
!------------------------------------------------------------------------------
!> Return a uniform vector field
!------------------------------------------------------------------------------
subroutine poss_vec_interp(self,cell,f,gop,val)
class(poss_vec_field), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
real(r8) :: pt(3)
IF(.NOT.ASSOCIATED(self%func))CALL oft_abort("No eval function specified", &
  "poss_vec_interp", __FILE__)
pt=self%mesh%log2phys(cell,f)
CALL self%func(pt,val,self%n)
end subroutine poss_vec_interp
!------------------------------------------------------------------------------
!> Setup analytic Taylor state interpolator for a cylindrical geometry
!------------------------------------------------------------------------------
subroutine cyl_taylor_setup(self,mesh)
class(cyl_taylor), intent(inout) :: self
class(oft_mesh), target, intent(inout) :: mesh
integer(i4) :: i,ierr
real(r8) :: zmin,zmax,rmax,pt(3),val(3)
!---
IF(self%zaxis>3.OR.self%zaxis<0)CALL oft_abort('Invalid z-axis','cyl_taylor_setup',__FILE__)
IF(ANY(self%rplane>3).OR.ANY(self%rplane<0))CALL oft_abort('Invalid xy-plane','cyl_taylor_setup',__FILE__)
self%mesh=>mesh
!---
zmin=1.d99
zmax=-1.d99; rmax=-1.d99
!$omp parallel do reduction(min:zmin) reduction(max:zmax) reduction(max:rmax)
DO i=1,self%mesh%np
  zmin=MIN(zmin,self%mesh%r(self%zaxis,i))
  zmax=MAX(zmax,self%mesh%r(self%zaxis,i))
  rmax=MAX(rmax,SQRT(SUM(self%mesh%r(self%rplane,i)**2)))
END DO
#ifdef HAVE_MPI
CALL MPI_ALLREDUCE(zmin,self%zmin,1,OFT_MPI_R8,MPI_MIN,oft_env%COMM,ierr)
#else
self%zmin=zmin
#endif
self%zmax=oft_mpi_max(zmax)
self%rmax=oft_mpi_max(rmax)
!---
self%akr=3.83170597020751d0/self%rmax
self%akz=pi/(self%zmax-self%zmin)
self%alm=SQRT(self%akr**2+self%akz**2)
!---
pt=0.d0
pt(self%zaxis)=(self%zmax+self%zmin)/2.d0
CALL cyl_taylor_eval(self,pt,val)
self%scale=1.d0/SUM(val**2)
end subroutine cyl_taylor_setup
!------------------------------------------------------------------------------
!> Evalute analytic Taylor state fields for a cylindrical geometry
!------------------------------------------------------------------------------
subroutine cyl_taylor_interp(self,cell,f,gop,val)
class(cyl_taylor), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
real(r8) :: pt(3),i,ar,az,s,c
!---
pt=self%mesh%log2phys(cell,f)
CALL cyl_taylor_eval(self,pt,val)
val=val*self%scale
end subroutine cyl_taylor_interp
!------------------------------------------------------------------------------
!> Evalute analytic Taylor state fields for a cylindrical geometry
!------------------------------------------------------------------------------
subroutine cyl_taylor_eval(self,pt,val)
class(cyl_taylor), intent(inout) :: self
real(r8), intent(in) :: pt(:)
real(r8), intent(out) :: val(:)
real(r8) :: i,ar,az,s,c
!---
ar=self%akr*SQRT(SUM(pt(self%rplane)**2))
az=self%akz*(pt(self%zaxis)-self%zmin)
if(ar==0)then
  i=.5d0
else
  i=dbesj1(ar)/ar
end if
s=SIN(az); c=COS(az)
!---
val(self%rplane(2))=(self%alm*pt(self%rplane(1))*s-self%akz*pt(self%rplane(2))*c)*i
val(self%rplane(1))=-(self%alm*pt(self%rplane(2))*s+self%akz*pt(self%rplane(1))*c)*i
val(self%zaxis)=dbesj0(ar)*s
end subroutine cyl_taylor_eval
!------------------------------------------------------------------------------
!> Setup infinite coil interpolation class
!------------------------------------------------------------------------------
subroutine inf_coils_setup(self,mesh)
class(inf_coils), intent(inout) :: self
class(oft_mesh), target, intent(inout) :: mesh
integer(i4) :: i
!---
IF(.NOT.ASSOCIATED(self%axis))CALL oft_abort('No axes set.','inf_coils_setup',__FILE__)
IF(.NOT.ASSOCIATED(self%center))CALL oft_abort('No centers set.','inf_coils_setup',__FILE__)
IF(.NOT.ASSOCIATED(self%current))CALL oft_abort('No currents set.','inf_coils_setup',__FILE__)
self%mesh=>mesh
!---
DO i=1,self%ncoils
  self%axis(:,i) = self%axis(:,i)/magnitude(self%axis(:,i))
END DO
end subroutine inf_coils_setup
!------------------------------------------------------------------------------
!> Evalute cummulative field from all coils
!------------------------------------------------------------------------------
subroutine inf_coils_interp(self,cell,f,gop,val)
class(inf_coils), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
integer(i4) :: i
real(r8) :: pt(3),r(3)
!---
pt=self%mesh%log2phys(cell,f)
!---
val=0.d0
DO i=1,self%ncoils
  r=pt-self%center(:,i)
  val = val + self%current(i)*cross_product(self%axis(:,i),r)/magnitude(r)**2
END DO
val = val*mu0/(4.d0*pi)
end subroutine inf_coils_interp
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine tor_radial_interp(self,cell,f,gop,val)
class(tor_radial), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
integer(i4) :: i
real(r8) :: pt(3),r,rhat(3),that(3),phi,theta
!---Get position in physical coordinates
pt=self%mesh%log2phys(cell,f)
!---Convert to toroidal coordinates
r=magnitude(pt(1:2))
phi=ATAN2(pt(2),pt(1))
theta=ATAN2(pt(3),r-self%R0)
!---Get poloidal slice unit vectors
rhat=pt-self%R0*(/COS(phi),SIN(phi),0.d0/)
r=magnitude(rhat)
rhat=rhat/r
that=cross_product(rhat,(/-SIN(phi),COS(phi),0.d0/))
!---Compute field
val = rhat*COS(self%m*theta + self%n*phi + self%phase0)*COS((r-self%rc)*pi/self%rwidth)*self%rwidth/pi &
    + that*SIN(self%m*theta + self%n*phi + self%phase0)*SIN((r-self%rc)*pi/self%rwidth)/REAL(self%m,8)
IF(ABS(r-self%rc)>self%rwidth/2.d0)val=0.d0
end subroutine tor_radial_interp
end module oft_vector_inits
