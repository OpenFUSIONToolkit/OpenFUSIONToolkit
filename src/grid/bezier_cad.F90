!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file bezier_cad.F90
!
!> CAD utility functions and class definition for reconstruction of
!! rational bezier CAD objects.
!!
!! This module supports internal refinement of T3D generated meshes while maintaining
!! the true CAD boundary. Refinement uses MINPACK to minimize the weighted sum of
!! distances between constraint points with the point constrained to the CAD object.
!!
!! @author Chris Hansen
!! @date June 2010
!! @ingroup doxy_oft_grid
!--------------------------------------------------------------------------------
module bezier_cad
use oft_base
implicit none
#include "local.h"
INTEGER(i4), PARAMETER :: cad_ngrid = 20 !< Number of grid points to use for object meshes
!------------------------------------------------------------------------------
!> CAD entity class
!------------------------------------------------------------------------------
type, abstract :: cad_entity
  integer(i4) :: id !< Input ID of CAD object
contains
  !> Physical point location
  procedure(cad_dummy_find), deferred :: locate
  !> Parametric point evaluation
  procedure(cad_dummy_eval), deferred :: eval
end type cad_entity
!------------------------------------------------------------------------------
!> CAD vertex class
!------------------------------------------------------------------------------
type, extends(cad_entity) :: cad_vertex
  real(r8), pointer :: pt(:) !< Vertex location [3]
contains
  !> Physical point location
  procedure :: locate => cad_vertex_find
  !> Parametric point evaluation
  procedure :: eval => cad_vertex_eval
  !> Reflect vertex
  procedure :: reflect => cad_vertex_reflect
end type cad_vertex
!------------------------------------------------------------------------------
!> CAD curve class
!!
!! A curve is defined in terms of the parametric coordinate \f$ u \f$ as
!! \f[ r (u) = \frac{\sum_0^n x_i \omega_i f^n_i (u) }{\sum_0^n \omega_i
!! f^n_i (u)} \f]
!! where \f$ n+1 \f$ is the order of the curve, \f$ \omega_i, x_i \f$ are
!! the weights and positions of the control points, and \f$ f^k_l \f$ is
!! the Bernstein polynomial of degree \f$ k \f$ and kind \f$ l \f$.
!------------------------------------------------------------------------------
type, extends(cad_entity) :: cad_curve
  integer(i4) :: order !< Curve order
  real(r8), pointer :: pt(:,:) !< Vertex locations [3,order]
  real(r8), pointer :: wt(:) !< Weight values [order]
  REAL(r8) :: rgrid(3,cad_ngrid) !< Object mesh
contains
  !> Generate grid over curve domain
  procedure :: grid => cad_curve_grid
  !> Physical point location
  procedure :: locate => cad_curve_find
  !> Parametric point evaluation
  procedure :: eval => cad_curve_eval
  !> Reflect curve
  procedure :: reflect => cad_curve_reflect
end type cad_curve
!------------------------------------------------------------------------------
!> CAD surface class
!!
!! A surface is defined in terms of the parametric coordinates \f$ u,v \f$ as
!! \f[ r (u,v) = \frac{ \sum_{i=0}^n \sum_{j=0}^m x_{ij} \omega_{ij} f^n_i (u)
!! f^m_j (v) }{\sum_{i=0}^n \sum_{j=0}^m \omega_{ij} f^n_i (u) f^m_j (v)} \f]
!! where \f$ n+1, m+1 \f$ are the orders of each bounding curve, \f$ \omega_{ij},
!! x_{ij} \f$ are the weights and positions of the control points, and \f$ f^k_l \f$
!! is the Bernstein polynomial of degree \f$ k \f$ and kind \f$ l \f$.
!------------------------------------------------------------------------------
type, extends(cad_entity) :: cad_surf
  integer(i4) :: order(2) !< Surface order [2]
  real(r8), pointer :: pt(:,:,:) !< Vertex locations [3,order(1),order(2)]
  real(r8), pointer :: wt(:,:) !< Weight values [order(1),order(2)]
  REAL(r8) :: rgrid(3,cad_ngrid,cad_ngrid) !< Object mesh
contains
  !> Generate grid over surface domain
  procedure :: grid => cad_surf_grid
  !> Physical point location
  procedure :: locate => cad_surf_find
  !> Parametric point evaluation
  procedure :: eval => cad_surf_eval
  !> Reflect surface
  procedure :: reflect => cad_surf_reflect
end type cad_surf
!------------------------------------------------------------------------------
!> List of CAD entities
!------------------------------------------------------------------------------
type :: cad_entity_ptr
  class(cad_entity), pointer :: wo => NULL() !< CAD object
end type cad_entity_ptr
!---CAD function prototypes
abstract interface
!------------------------------------------------------------------------------
!> Map the parametric possition on an entity to physical coordinates
!! - (u,v) -> (x,y,z)
!------------------------------------------------------------------------------
  subroutine cad_dummy_eval(self,pt,u,v)
    import :: cad_entity, r8
    class(cad_entity), intent(in) :: self
    real(r8), intent(out) :: pt(3) !< Position vector [3]
    real(r8), intent(in) :: u !< Parametric coordinate 1
    real(r8), intent(in) :: v !< Parametric coordinate 2
  end subroutine cad_dummy_eval
!------------------------------------------------------------------------------
!> Find the parametric representation of a boundary point in CAD representation
!! - (x,y,z) -> (u,v)
!------------------------------------------------------------------------------
  subroutine cad_dummy_find(self,pt,u,v)
    import :: cad_entity, r8
    class(cad_entity), intent(in) :: self
    real(r8), intent(in) :: pt(3) !< Position vector [3]
    real(r8), intent(out) :: u !< Parametric coordinate 1
    real(r8), intent(out) :: v !< Parametric coordinate 2
  end subroutine cad_dummy_find
end interface
private
!---Fitting variables
CLASS(cad_curve), POINTER :: active_curve => NULL() !< Active curve for MINPACK fitting
CLASS(cad_surf), POINTER :: active_surf => NULL() !< Active surface for MINPACK fitting
REAL(r8) :: active_endpts(3,3) !< Active constraint points for MINPACK fitting
REAL(r8) :: active_wts(3) !< Active constraint weights for MINPACK fitting
!$omp threadprivate(active_endpts,active_wts,active_curve,active_surf)
public cad_entity, cad_vertex, cad_curve, cad_surf, cad_entity_ptr, bernstein
public cad_curve_midpoint, cad_surf_midpoint, cad_surf_center, cad_curve_cast, cad_surf_cast
contains
!------------------------------------------------------------------------------
!> Map parametric possition of a vertex to physical coordinates
!! - (u,v) -> (x,y,z)
!------------------------------------------------------------------------------
subroutine cad_vertex_eval(self,pt,u,v)
class(cad_vertex), intent(in) :: self
real(r8), intent(out) :: pt(3) !< Position vector [3]
real(r8), intent(in) :: u !< Parametric coordinate 1 (ignored)
real(r8), intent(in) :: v !< Parametric coordinate 2 (ignored)
DEBUG_STACK_PUSH
pt=self%pt
DEBUG_STACK_POP
end subroutine cad_vertex_eval
!------------------------------------------------------------------------------
!> Find the parametric representation of a boundary point in CAD representation on a vertex
!! - (x,y,z) -> (u,v)
!------------------------------------------------------------------------------
subroutine cad_vertex_find(self,pt,u,v)
class(cad_vertex), intent(in) :: self
real(r8), intent(in) :: pt(3) !< Position vector [3]
real(r8), intent(out) :: u !< Parametric coordinate 1 (ignored)
real(r8), intent(out) :: v !< Parametric coordinate 2 (ignored)
DEBUG_STACK_PUSH
u=0.d0
v=0.d0
DEBUG_STACK_POP
end subroutine cad_vertex_find
!------------------------------------------------------------------------------
!> Reflect a vertex object across a given plane
!------------------------------------------------------------------------------
subroutine cad_vertex_reflect(self,copy,tol,k)
class(cad_vertex), intent(in) :: self !< Source object to copy
class(cad_vertex), intent(out) :: copy !< Reflected copy of the source vertex
real(r8), intent(in) :: tol !< Minimum distance from plane
integer(i4), intent(in) :: k !< Coordinate index for the reflection plane
integer(i4) :: k1,k2
DEBUG_STACK_PUSH
k1=mod(k ,3)+1
k2=mod(k1,3)+1
allocate(copy%pt(3))
if(abs(self%pt(k))<=tol)then ! Object on Ref-Plane
  copy%pt=self%pt
else ! Object not on Ref-Plane
  copy%pt(k)=-self%pt(k)
  copy%pt(k1)=self%pt(k1)
  copy%pt(k2)=self%pt(k2)
end if
DEBUG_STACK_POP
end subroutine cad_vertex_reflect
!---------------------------------------------------------------------------------
!> Cast \ref bezier_cad::cad_entity "cad_entity" to \ref bezier_cad::cad_curve
!! "cad_curve"
!!
!! @result Error flag
!---------------------------------------------------------------------------------
FUNCTION cad_curve_cast(self,source) result(ierr)
class(cad_curve), pointer, intent(out) :: self !< Object of desired type, unassociated if cast fails
class(cad_entity), target, intent(in) :: source !< Source object to cast
integer(i4) :: ierr
DEBUG_STACK_PUSH
select type(source)
  class is(cad_curve)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
DEBUG_STACK_POP
end FUNCTION cad_curve_cast
!------------------------------------------------------------------------------
!> Map the parametric possition on a curve to physical coordinates.
!! - (u,v) -> (x,y,z)
!------------------------------------------------------------------------------
subroutine cad_curve_eval(self,pt,u,v)
class(cad_curve), intent(in) :: self
real(r8), intent(out) :: pt(3) !< Position vector [3]
real(r8), intent(in) :: u !< Parametric coordinate 1
real(r8), intent(in) :: v !< Parametric coordinate 2 (ignored)
integer(i4) :: i
real(r8) :: norm
DEBUG_STACK_PUSH
!---Initialize point reconstruction
pt = 0.d0
norm = 0.d0
!---Reconstruct point from CAD curve representation
do i=1,self%order ! Loop over weight points
  pt = pt + self%wt(i)*self%pt(:,i)*bernstein(u,self%order-1,i-1)
  norm = norm + self%wt(i)*bernstein(u,self%order-1,i-1)
end do
pt=pt/norm
DEBUG_STACK_POP
end subroutine cad_curve_eval
!------------------------------------------------------------------------------
!> Reflect a curve object across a given plane
!------------------------------------------------------------------------------
subroutine cad_curve_reflect(self,copy,tol,k)
class(cad_curve), intent(in) :: self !< Source object to copy
class(cad_curve), intent(out) :: copy !< Reflected copy of the source curve
real(r8), intent(in) :: tol !< Minimum distance from plane
integer(i4), intent(in) :: k !< Coordinate index for the reflection plane
integer(i4) :: i,k1,k2
DEBUG_STACK_PUSH
k1=mod(k ,3)+1
k2=mod(k1,3)+1
copy%order=self%order
allocate(copy%wt(copy%order),copy%pt(3,copy%order))
copy%wt=self%wt
do i=1,self%order
  if(abs(self%pt(k,i))<=tol)then ! Object on Ref-Plane
    copy%pt(:,i)=self%pt(:,i)
  else ! Object not on Ref-Plane
    copy%pt(k,i)=-self%pt(k,i)
    copy%pt(k1,i)=self%pt(k1,i)
    copy%pt(k2,i)=self%pt(k2,i)
  end if
end do
DEBUG_STACK_POP
end subroutine cad_curve_reflect
!------------------------------------------------------------------------------
!> Evalute a grid of points evenly spaced in parametric space
!------------------------------------------------------------------------------
subroutine cad_curve_grid(self)
class(cad_curve), intent(inout) :: self
real(r8) :: u,v
integer(i4) :: i
DEBUG_STACK_PUSH
v=0.d0
do i=1,cad_ngrid
  u=(i-1)*1.d0/(cad_ngrid-1)
  call self%eval(self%rgrid(:,i),u,v)
end do
DEBUG_STACK_POP
end subroutine cad_curve_grid
!------------------------------------------------------------------------------
!> Compute unit tangent vector for a CAD curve at a given location
!------------------------------------------------------------------------------
subroutine cad_curve_tang(self,tang,u)
class(cad_curve), intent(in) :: self
real(r8), intent(out) :: tang(3) !< Curve tangent unit vector [3]
real(r8), intent(in) :: u !< Parametric coordinate
real(r8) :: pt1(3),pt2(3),du=1.E-9_r8
DEBUG_STACK_PUSH
!---
CALL self%eval(pt1,u-du,0._r8)
CALL self%eval(pt2,u+du,0._r8)
!---
tang=pt2-pt1
tang=tang/SQRT(SUM(tang**2))
DEBUG_STACK_POP
end subroutine cad_curve_tang
!------------------------------------------------------------------------------
!> Find the parametric representation of a boundary point in CAD representation on a curve
!! - (x,y,z) -> (u,v)
!------------------------------------------------------------------------------
subroutine cad_curve_find(self,pt,u,v)
class(cad_curve), intent(in) :: self
real(r8), intent(in) :: pt(3) !< Position vector [3]
real(r8), intent(out) :: u !< Parametric coordinate 1
real(r8), intent(out) :: v !< Parametric coordinate 2 (ignored)
real(r8) :: du,un,val(3),rt(3),fp
integer :: i
DEBUG_STACK_PUSH
v=0.d0
!---Check if point is vertex 1
u=0.d0
call cad_curve_eval(self,rt,u,v)
val(2)=sqrt(sum((pt-rt)**2))
IF(val(2)<1.d-12)THEN ! Point is vertex 1
  DEBUG_STACK_POP
  RETURN
END IF
!---Check if point is vertex 2
u=1.d0
call cad_curve_eval(self,rt,u,v)
val(2)=sqrt(sum((pt-rt)**2))
IF(val(2)<1.d-12)THEN ! Point is vertex 2
  DEBUG_STACK_POP
  RETURN
END IF
!---Initialize search point
u=.5
un=.5
fp=1.d99
du=1.d-9
!---Use Newton's Method to find point
do i=1,40
  do while(.TRUE.)
    call cad_curve_eval(self,rt,un,v)
    val(2)=sqrt(sum((pt-rt)**2))
    if(val(2)<fp)exit
    !---Backtrack if necessary
    un=(un-u)*.5d0+u
    if(abs(un-u)<1.d-5)exit
  end do
  !---Update search point and error
  u=un
  fp=val(2)
  if(val(2)<1.d-6)exit
  !---Compute gradient
  call cad_curve_eval(self,rt,u-du,v)
  val(1)=sqrt(sum((pt-rt)**2))
  call cad_curve_eval(self,rt,u+du,v)
  val(3)=sqrt(sum((pt-rt)**2))
  !---Prevent divide by zero
  if(abs(val(3)-val(1))<1.d-14)then
    u=u-du/2
    cycle
  endif
  !---Update solution
  un=u-val(2)*du*2/(val(3)-val(1))
  !---Limit minimum and maximum values
  if(un>1.d0)un=1.d0
  if(un<0.d0)un=0.d0
enddo
!---Catch failure to converge
if(i>40)then
  write(*,*)self%id,u,val(2)
  write(*,*)pt
  write(*,*)rt
  call oft_abort('Did not converge','cad_curve_find',__FILE__)
endif
DEBUG_STACK_POP
end subroutine cad_curve_find
!------------------------------------------------------------------------------
!> Compute the weighted midpoint of a curve edge
!!
!! Locates the point on a given CAD curve which minimizes the weighted sum
!! of distances to 2 constraint points.
!!
!! \f[ \sum_i w_i*(r_n - p_i)^2 \f]
!------------------------------------------------------------------------------
subroutine cad_curve_midpoint(self,pt,pt1,pt2,wt1,wt2,ierr)
class(cad_curve),target, intent(in) :: self !< Curve object
real(r8), intent(inout) :: pt(3) !< Solution point [3]
real(r8), intent(in) :: pt1(3) !< Constraint point 1 [3]
real(r8), intent(in) :: pt2(3) !< Constraint point 2 [3]
real(r8), intent(in) :: wt1 !< Weight for constraint point 1
real(r8), intent(in) :: wt2 !< Weight for constraint point 2
integer(i4), intent(out) :: ierr !< Error flag
real(r8) :: u,v,val(2),rt(3)
integer :: i
!---
integer(i4), parameter :: nerr=6
integer(i4), parameter :: neq=1
real(r8) :: diag(neq),fjac(nerr,neq),qtf(neq),uv(neq)
real(r8) :: wa1(neq),wa2(neq),wa3(neq),wa4(nerr),error(nerr)
integer(i4) :: ipvt(neq)
real(r8) :: ftol,xtol,gtol,epsfcn,factor,dmin
integer :: info,ldfjac,maxfev,mode,nfev,nprint
DEBUG_STACK_PUSH
ierr=0
v=0.d0
!---Set initial guess
u=.5d0
!---Set initial guess
dmin=1.d99
do i=1,cad_ngrid
  rt=self%rgrid(:,i)
  val(2)=sqrt(sum((pt-rt)**2))
  if(val(2)<dmin)then
    dmin=val(2)
    u=(i-1)*1.d0/(cad_ngrid-1)
  end if
end do
!---Setup MINPACK parameters
mode = 1
factor = 1.d0
maxfev = 100
ftol = 1.d-6
xtol = 1.d-6
gtol = 1.d-5
epsfcn = 1.d-6
nprint = 0
ldfjac = nerr
!---Setup active points
active_curve=>self
active_endpts(1,:)=pt1; active_wts(1)=wt1
active_endpts(2,:)=pt2; active_wts(2)=wt2
uv=(/u/)
!---
CALL lmdif(cad_cmid_error,nerr,neq,uv,error, &
           ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
           nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!---
IF(info>4)THEN
  !WRITE(*,*)'Curve midpoint failed',info,nfev
  ierr=-1
ELSE
  ierr=0
END IF
u=uv(1)
v=0.d0
call self%eval(pt,u,v)
DEBUG_STACK_POP
end subroutine cad_curve_midpoint
!------------------------------------------------------------------------------
!> Evalute the error between a curve point and the current active points
!! used in a 2 point minimization.
!!
!! @note Designed to be used as the error function for minimization in
!! @ref bezier_cad::cad_curve_midpoint "cad_curve_midpoint"
!------------------------------------------------------------------------------
subroutine cad_cmid_error(m,n,uv,err,iflag)
integer(i4), intent(in) :: m !< Number of spatial dimensions (3)
integer(i4), intent(in) :: n !< Number of parametric dimensions (2)
real(r8), intent(in) :: uv(n) !< Parametric position [n]
real(r8), intent(out) :: err(m) !< Error vector between current and desired point [3]
integer(i4), intent(inout) :: iflag !< Unused flag
real(r8) :: pt(3)
DEBUG_STACK_PUSH
call active_curve%eval(pt,uv(1),0.d0)
err(1:3)=(active_endpts(1,:)-pt)*SQRT(active_wts(1))
err(4:6)=(active_endpts(2,:)-pt)*SQRT(active_wts(2))
DEBUG_STACK_POP
end subroutine cad_cmid_error
!---------------------------------------------------------------------------------
!> Cast \ref bezier_cad::cad_entity "cad_entity" to \ref bezier_cad::cad_surf
!! "cad_surf"
!!
!! @result Error flag
!---------------------------------------------------------------------------------
FUNCTION cad_surf_cast(self,source) result(ierr)
class(cad_surf), pointer, intent(out) :: self !< Object of desired type, unassociated if cast fails
class(cad_entity), target, intent(in) :: source !< Source object to cast
integer(i4) :: ierr
DEBUG_STACK_PUSH
select type(source)
  class is(cad_surf)
    self=>source
    ierr=0
  class default
    ierr=-1
end select
DEBUG_STACK_POP
end FUNCTION cad_surf_cast
!------------------------------------------------------------------------------
!> Map the parametric possition on a surface to physical coordinates.
!! - (u,v) -> (x,y,z)
!------------------------------------------------------------------------------
subroutine cad_surf_eval(self,pt,u,v)
class(cad_surf), intent(in) :: self
real(r8), intent(out) :: pt(3) !< Position vector [3]
real(r8), intent(in) :: u !< Parametric coordinate 1
real(r8), intent(in) :: v !< Parametric coordinate 2
integer(i4) :: i,j
real(r8) :: norm
DEBUG_STACK_PUSH
!---Initialize point reconstruction
pt = 0.d0
norm = 0.d0
!---Reconstruct point from CAD surface representation
do i=1,self%order(1) ! Loop over weight points
  do j=1,self%order(2)
    pt = pt + self%wt(i,j)*self%pt(:,i,j)*bernstein(u,self%order(1)-1,i-1)*bernstein(v,self%order(2)-1,j-1)
    norm = norm + self%wt(i,j)*bernstein(u,self%order(1)-1,i-1)*bernstein(v,self%order(2)-1,j-1)
  end do
end do
pt=pt/norm
DEBUG_STACK_POP
end subroutine cad_surf_eval
!------------------------------------------------------------------------------
!> Find the parametric representation of a boundary point in CAD representation on a surface.
!! - (x,y,z) -> (u,v)
!------------------------------------------------------------------------------
subroutine cad_surf_find(self,pt,u,v)
class(cad_surf), intent(in) :: self
real(r8), intent(in) :: pt(3) !< Position vector [3]
real(r8), intent(out) :: u !< Parametric coordinate 1
real(r8), intent(out) :: v !< Parametric coordinate 2
real(r8) :: du,un,vn,val(3),rt(3),grad(2),alpha,fp
integer :: i
DEBUG_STACK_PUSH
!---Check if point is vertex 1
u=0.d0; v=0.d0
call cad_surf_eval(self,rt,u,v)
val(2)=sqrt(sum((pt-rt)**2))
IF(val(2)<1.d-12)THEN ! Point is vertex 1
  DEBUG_STACK_POP
  RETURN
END IF
!---Check if point is vertex 2
u=1.d0; v=0.d0
call cad_surf_eval(self,rt,u,v)
val(2)=sqrt(sum((pt-rt)**2))
IF(val(2)<1.d-12)THEN ! Point is vertex 2
  DEBUG_STACK_POP
  RETURN
END IF
!---Check if point is vertex 3
u=1.d0; v=1.d0
call cad_surf_eval(self,rt,u,v)
val(2)=sqrt(sum((pt-rt)**2))
IF(val(2)<1.d-12)THEN ! Point is vertex 3
  DEBUG_STACK_POP
  RETURN
END IF
!---Check if point is vertex 4
u=0.d0; v=1.d0
call cad_surf_eval(self,rt,u,v)
val(2)=sqrt(sum((pt-rt)**2))
IF(val(2)<1.d-12)THEN ! Point is vertex 4
  DEBUG_STACK_POP
  RETURN
END IF
!---Use 2-D Newton's Method to find point
u=.3d0; un=.3d0
v=.7d0; vn=.7d0
fp=1.d99
du=1.d-9
do i=1,3000
  do while(.TRUE.)
    call cad_surf_eval(self,rt,un,vn)
    val(2)=sqrt(sum((pt-rt)**2))
    if(val(2)<fp)exit
    !---Backtrack if necessary
    un=(un-u)*.5d0+u
    vn=(vn-v)*.5d0+v
    if(abs(u-un)<1.d-5)exit
  end do
  !---Update search point and error
  u=un; v=vn
  if(val(2)<1.d-5)exit
  fp=val(2)
  !---Compute gradient in u
  call cad_surf_eval(self,rt,u-du,v)
  val(1)=sqrt(sum((pt-rt)**2))
  call cad_surf_eval(self,rt,u+du,v)
  val(3)=sqrt(sum((pt-rt)**2))
  grad(1)=(val(3)-val(1))/(2*du)
  !---Compute gradient in v
  call cad_surf_eval(self,rt,u,v-du)
  val(1)=sqrt(sum((pt-rt)**2))
  call cad_surf_eval(self,rt,u,v+du)
  val(3)=sqrt(sum((pt-rt)**2))
  grad(2)=(val(3)-val(1))/(2*du)
  !---Prevent divide by zero
  if(sum(grad**2)<1.d-14)then
    u=u-du/2.d0
    v=v-du/2.d0
    cycle
  end if
  !---Update solution
  un=u-val(2)*grad(1)/sum(grad**2)
  vn=v-val(2)*grad(2)/sum(grad**2)
  !---Limit minimum and maximum values
  if(un>1.d0)un=1.d0
  if(vn>1.d0)vn=1.d0
  if(un<0.d0)un=0.d0
  if(vn<0.d0)vn=0.d0
end do
!---Catch failure to converge
if(i>2000)then
  write(*,*)pt
  write(*,*)rt
  write(*,*)u,v,val(2),alpha
  call oft_abort('Did note converge','cad_surf_find',__FILE__)
end if
DEBUG_STACK_POP
end subroutine cad_surf_find
!------------------------------------------------------------------------------
!> Reflect a surface object across a given plane
!------------------------------------------------------------------------------
subroutine cad_surf_reflect(self,copy,tol,k)
class(cad_surf), intent(in) :: self !< Source surface to copy
class(cad_surf), intent(out) :: copy !< Reflected copy of the source surface
real(r8), intent(in) :: tol !< Minimum distance from plane
integer(i4), intent(in) :: k !< Coordinate index for the reflection plane
integer(i4) :: i,j,k1,k2
DEBUG_STACK_PUSH
k1=mod(k ,3)+1
k2=mod(k1,3)+1
copy%order=self%order
allocate(copy%wt(copy%order(1),copy%order(2)),copy%pt(3,copy%order(1),copy%order(2)))
copy%wt=self%wt
do i=1,self%order(1)
  do j=1,self%order(2)
    if(abs(self%pt(k,i,j))<=tol)then ! Object on Ref-Plane
      copy%pt(:,i,j)=self%pt(:,i,j)
    else ! Object not on Ref-Plane
      copy%pt(k,i,j)=-self%pt(k,i,j)
      copy%pt(k1,i,j)=self%pt(k1,i,j)
      copy%pt(k2,i,j)=self%pt(k2,i,j)
    end if
  end do
end do
DEBUG_STACK_POP
end subroutine cad_surf_reflect
!------------------------------------------------------------------------------
!> Compute unit normal vector for a CAD surface at a specified position
!------------------------------------------------------------------------------
subroutine cad_surf_norm(self,norm,u,v)
class(cad_surf), intent(in) :: self
real(r8), intent(out) :: norm(3) !< Surface normal unit vector [3]
real(r8), intent(in) :: u !< Parametric coordinate 1
real(r8), intent(in) :: v !< Parametric coordinate 2
real(r8) :: pt1(3),pt2(3),du=1.E-9_r8,ru(3),rv(3)
DEBUG_STACK_PUSH
!---
CALL self%eval(pt1,u-du,v)
CALL self%eval(pt2,u+du,v)
ru=pt2-pt1
ru=ru/SQRT(SUM(ru**2))
!---
CALL self%eval(pt1,u,v-du)
CALL self%eval(pt2,u,v+du)
rv=pt2-pt1
rv=rv/SQRT(SUM(rv**2))
!---
norm=-1.d99
DEBUG_STACK_POP
end subroutine cad_surf_norm
!------------------------------------------------------------------------------
!> Evalute a grid of points evenly spaced in parametric space
!------------------------------------------------------------------------------
subroutine cad_surf_grid(self)
class(cad_surf), intent(inout) :: self
real(r8) :: u,v
integer(i4) :: i,j
DEBUG_STACK_PUSH
do i=1,cad_ngrid
  u=(i-1)*1.d0/(cad_ngrid-1)
  do j=1,cad_ngrid
    v=(j-1)*1.d0/(cad_ngrid-1)
    call self%eval(self%rgrid(:,i,j),u,v)
  end do
end do
DEBUG_STACK_POP
end subroutine cad_surf_grid
!------------------------------------------------------------------------------
!> Compute the weighted midpoint of a surface edge
!!
!! Locates the point on a given CAD surface which minimizes the weighted sum
!! of distances to 2 constraint points.
!!
!! \f[ \sum_i w_i*(r_n - p_i)^2 \f]
!------------------------------------------------------------------------------
subroutine cad_surf_midpoint(self,pt,pt1,pt2,wt1,wt2,ierr)
class(cad_surf), target, intent(in) :: self !< Surface object
real(r8), intent(inout) :: pt(3) !< Solution point [3]
real(r8), intent(in) :: pt1(3) !< Constraint point 1 [3]
real(r8), intent(in) :: pt2(3) !< Constraint point 2 [3]
real(r8), intent(in) :: wt1 !< Weight for constraint point 1
real(r8), intent(in) :: wt2 !< Weight for constraint point 2
integer(i4), intent(out) :: ierr !< Error flag
logical :: sing
real(r8) :: u,v
!---
integer(i4), parameter :: nerr=6
integer(i4), parameter :: neq=2
real(r8) :: diag(neq),fjac(nerr,neq),qtf(neq),uv(neq)
real(r8) :: wa1(neq),wa2(neq),wa3(neq),wa4(nerr),error(nerr)
integer(i4) :: ipvt(neq)
real(r8) :: ftol,xtol,gtol,epsfcn,factor
real(r8) :: val(3),rt(3),dmin
integer :: i,j,info,ldfjac,maxfev,mode,nfev,nprint
DEBUG_STACK_PUSH
!---
sing=.FALSE.
ierr=0
!---Set initial guess
dmin=1.d99
do i=1,cad_ngrid
  do j=1,cad_ngrid
    rt=self%rgrid(:,i,j)
    val(2)=sqrt(sum((pt-rt)**2))
    if(val(2)<dmin)then
      dmin=val(2)
      u=(i-1)*1.d0/(cad_ngrid-1)
      v=(j-1)*1.d0/(cad_ngrid-1)
    end if
  end do
end do
!---Setup MINPACK parameters
mode = 1
factor = 1.d0
maxfev = 100
ftol = 1.d-6
xtol = 1.d-6
gtol = 1.d-5
epsfcn = 1.d-6
nprint = 0
ldfjac = nerr
!---Setup active points
active_surf=>self
active_endpts(1,:)=pt1; active_wts(1)=wt1
active_endpts(2,:)=pt2; active_wts(2)=wt2
uv=(/u,v/)
!---
CALL lmdif(cad_smid_error,nerr,neq,uv,error, &
           ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
           nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!---
IF(info>4)THEN
  !WRITE(*,*)'Surface midpoint failed',info,nfev
  ierr=-1
ELSE
  ierr=0
END IF
u=uv(1)
v=uv(2)
call self%eval(pt,u,v)
DEBUG_STACK_POP
end subroutine cad_surf_midpoint
!------------------------------------------------------------------------------
!> Evalute the error between a surface point and the current active points
!! used in a 2 point minimization.
!!
!! @note Designed to be used as the error function for minimization in
!! @ref bezier_cad::cad_surf_midpoint "cad_surf_midpoint"
!------------------------------------------------------------------------------
subroutine cad_smid_error(m,n,uv,err,iflag)
integer(i4), intent(in) :: m !< Number of spatial dimensions (3)
integer(i4), intent(in) :: n !< Number of parametric dimensions (2)
real(r8), intent(in) :: uv(n) !< Parametric possition [n]
real(r8), intent(out) :: err(m) !< Error vector between current and desired point [3]
integer(i4), intent(inout) :: iflag !< Unused flag
real(r8) :: pt(3)
DEBUG_STACK_PUSH
call active_surf%eval(pt,uv(1),uv(2))
err(1:3)=(active_endpts(1,:)-pt)*SQRT(active_wts(1))
err(4:6)=(active_endpts(2,:)-pt)*SQRT(active_wts(2))
DEBUG_STACK_POP
end subroutine cad_smid_error
!------------------------------------------------------------------------------
!> Compute the weighted center point of a surface triangle
!!
!! Locates the point on a given CAD surface which minimizes the weighted sum
!! of distances to 3 constraint points.
!!
!! \f[ \sum_i w_i*(r_n - p_i)^2 \f]
!------------------------------------------------------------------------------
subroutine cad_surf_center(self,pt,pt1,pt2,pt3,wt1,wt2,wt3,ierr)
class(cad_surf), target, intent(in) :: self !< Surface object
real(r8), intent(inout) :: pt(3) !< Solution point [3]
real(r8), intent(in) :: pt1(3) !< Constraint point 1 [3]
real(r8), intent(in) :: pt2(3) !< Constraint point 2 [3]
real(r8), intent(in) :: pt3(3) !< Constraint point 3 [3]
real(r8), intent(in) :: wt1 !< Weight for constraint point 1
real(r8), intent(in) :: wt2 !< Weight for constraint point 2
real(r8), intent(in) :: wt3 !< Weight for constraint point 3
integer(i4), intent(out) :: ierr !< Error flag
logical :: sing
real(r8) :: u,v
!---
integer(i4), parameter :: nerr=9
integer(i4), parameter :: neq=2
real(r8) :: diag(neq),fjac(nerr,neq),qtf(neq),uv(neq)
real(r8) :: wa1(neq),wa2(neq),wa3(neq),wa4(nerr),error(nerr)
integer(i4) :: ipvt(neq)
real(r8) :: ftol,xtol,gtol,epsfcn,factor
real(r8) :: val(3),rt(3),dmin
integer :: i,j,info,ldfjac,maxfev,mode,nfev,nprint
DEBUG_STACK_PUSH
!---
sing=.FALSE.
ierr=0
!---Set initial guess
dmin=1.d99
do i=1,cad_ngrid
  do j=1,cad_ngrid
    rt=self%rgrid(:,i,j)
    val(2)=sqrt(sum((pt-rt)**2))
    if(val(2)<dmin)then
      dmin=val(2)
      u=(i-1)*1.d0/(cad_ngrid-1)
      v=(j-1)*1.d0/(cad_ngrid-1)
    end if
  end do
end do
!u=.3
!v=.7
!---
mode = 1
factor = 1.d0
maxfev = 100
ftol = 1.d-6
xtol = 1.d-6
gtol = 1.d-5
epsfcn = 1.d-6
nprint = 0
ldfjac = nerr
!---
active_surf=>self
active_endpts(1,:)=pt1; active_wts(1)=wt1
active_endpts(2,:)=pt2; active_wts(2)=wt2
active_endpts(3,:)=pt3; active_wts(3)=wt3
uv=(/u,v/)
!---
CALL lmdif(cad_scenter_error,nerr,neq,uv,error, &
           ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
           nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!---
IF(info>4)THEN
  !WRITE(*,*)'Surface center failed',info,nfev
  ierr=-1
ELSE
  ierr=0
END IF
u=uv(1)
v=uv(2)
call self%eval(pt,u,v)
DEBUG_STACK_POP
end subroutine cad_surf_center
!------------------------------------------------------------------------------
!> Evalute the error between a surface point and the current active points
!! used in a 3 point minimization.
!!
!! @note Designed to be used as the error function for minimization in
!! @ref bezier_cad::cad_surf_center "cad_surf_center"
!------------------------------------------------------------------------------
subroutine cad_scenter_error(m,n,uv,err,iflag)
integer(i4), intent(in) :: m !< Number of spatial dimensions (3)
integer(i4), intent(in) :: n !< Number of parametric dimensions (2)
real(r8), intent(in) :: uv(n) !< Parametric possition [n]
real(r8), intent(out) :: err(m) !< Error vector between current and desired point [3]
integer(i4), intent(inout) :: iflag !< Unused flag
real(r8) :: pt(3)
DEBUG_STACK_PUSH
call active_surf%eval(pt,uv(1),uv(2))
err(1:3)=(active_endpts(1,:)-pt)*SQRT(active_wts(1))
err(4:6)=(active_endpts(2,:)-pt)*SQRT(active_wts(2))
err(7:9)=(active_endpts(3,:)-pt)*SQRT(active_wts(3))
DEBUG_STACK_POP
end subroutine cad_scenter_error
!------------------------------------------------------------------------------
!> Evaluates the bernstein polynomial.
!!
!! f(x) = \f$ B_{ij}(x) \f$
!!
!! @result \f$ B_{ij}(x) \f$
!------------------------------------------------------------------------------
real(r8) function bernstein(x,i,j)
real(r8), intent(in) :: x !< Point to evaluate
integer(i4), intent(in) :: i !< Polynomial degree
integer(i4), intent(in) :: j !< Polynomial kind
DEBUG_STACK_PUSH
bernstein=0.d0
select case(i)
  case(1) ! Degree 1
    if(j==0)then
      bernstein = 1.d0-x
    elseif(j==1)then
      bernstein = x
    endif
  case(2) ! Degree 2
    if(j==0)then
      bernstein = (1-x)**2
    elseif(j==1)then
      bernstein = 2.d0*x*(1-x)
    elseif(j==2)then
      bernstein = x**2
    endif
  case(3) ! Degree 3
    if(j==0)then
      bernstein = (1.d0-x)**3
    elseif(j==1)then
      bernstein = 3.d0*x*(1.d0-x)**2
    elseif(j==2)then
      bernstein = 3.d0*(x**2)*(1.d0-x)
    elseif(j==3)then
      bernstein = x**3
    endif
end select
DEBUG_STACK_POP
end function bernstein
end module bezier_cad
