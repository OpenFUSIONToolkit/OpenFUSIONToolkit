!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file nurbs_cad.F90
!
!> CAD utility functions and class definition for reconstruction of
!! non-uniform rational bezier CAD objects.
!!
!! This module supports internal refinement of CUBIT generated meshes while maintaining
!! the true CAD boundary. Refinement uses MINPACK to minimize the weighted sum of
!! distances between constraint points with the point constrained to the CAD object.
!! Evaluation of CAD objects is provided by an interface to the OpenNURBS library.
!!
!! @author Chris Hansen
!! @date May 2012
!! @ingroup doxy_oft_grid
!-----------------------------------------------------------------------------
MODULE nurbs_cad
#ifdef HAVE_ONURBS
USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_double, c_char
USE oft_base
IMPLICIT NONE
#include "local.h"
PRIVATE
!---
INTEGER(i4), PARAMETER :: nurbs_ngrid = 20 !< Number of grid points to use for object meshes
REAL(r8), PARAMETER :: nurbs_singstep = .1d0 !< Distance to offset guesses from singular edges
!---------------------------------------------------------------------------
! CLASS nurbs_entity
!---------------------------------------------------------------------------
!> CAD entity class
!---------------------------------------------------------------------------
TYPE, ABSTRACT :: nurbs_entity
  INTEGER(i4) :: id !< Input ID of CAD object
  INTEGER(i4) :: reflect = -1 !< Index to reflect
  REAL(r8) :: coord_scales(3) = 1.d0 !< Scaling for each coordinate
CONTAINS
  !> Physical point location
  PROCEDURE(nurbs_dummy_grid), DEFERRED :: grid
  !> Physical point location
  PROCEDURE(nurbs_dummy_find), DEFERRED :: locate
  !> Parametric point evaluation
  PROCEDURE(nurbs_dummy_eval), DEFERRED :: eval
  !> Parametric point evaluation
  PROCEDURE(nurbs_dummy_wrap), DEFERRED :: wrap
  !> Parametric point evaluation
  PROCEDURE(nurbs_dummy_unwrap), DEFERRED :: unwrap
END TYPE nurbs_entity
!---------------------------------------------------------------------------
! CLASS nurbs_curve
!---------------------------------------------------------------------------
!> CAD curve class
!---------------------------------------------------------------------------
TYPE, EXTENDS(nurbs_entity) :: nurbs_curve
  INTEGER(i4) :: cid !< Curve ID in OFT indexing
  INTEGER(i4) :: nspan !< Size of span vector
  REAL(r8), POINTER, DIMENSION(:) :: span !< Currently unused
  REAL(r8) :: domain(2) !< Extent in parametric domain
  LOGICAL :: periodic = .FALSE. !< Curve is periodic
  LOGICAL :: linear = .FALSE. !< Curve is straight
  REAL(r8) :: rgrid(3,nurbs_ngrid) !< Object mesh
CONTAINS
  !> Generate grid over curve domain
  PROCEDURE :: grid => nurbs_curve_grid
  !> Physical point location
  PROCEDURE :: locate => nurbs_curve_find
  !> Parametric point evaluation
  PROCEDURE :: eval => nurbs_curve_eval
  !> Physical point location
  PROCEDURE :: wrap => nurbs_curve_wrap
  !> Parametric point evaluation
  PROCEDURE :: unwrap => nurbs_curve_unwrap
END TYPE nurbs_curve
!---------------------------------------------------------------------------
! CLASS nurbs_surf
!---------------------------------------------------------------------------
!> CAD surface class
!---------------------------------------------------------------------------
TYPE, EXTENDS(nurbs_entity) :: nurbs_surf
  INTEGER(i4) :: sid !< Surface ID in OFT indexing
  INTEGER(i4) :: nspan(2) !< Size of span vectors
  REAL(r8), POINTER, DIMENSION(:) :: span1 !< Currently unused
  REAL(r8), POINTER, DIMENSION(:) :: span2 !< Currently unused
  REAL(r8) :: domain(2,2) !< Extent in parametric domain
  LOGICAL :: periodic(2) = .FALSE. !< Surface is periodic
  LOGICAL :: singular(2,2) = .FALSE. !< Surface edge is singular
  LOGICAL :: planar = .FALSE. !< Surface is flat
  REAL(r8) :: rgrid(3,nurbs_ngrid,nurbs_ngrid) !< Object mesh
CONTAINS
  !> Generate grid over surface domain
  PROCEDURE :: grid => nurbs_surf_grid
  !> Physical point location
  PROCEDURE :: locate => nurbs_surf_find
  !> Parametric point evaluation
  PROCEDURE :: eval => nurbs_surf_eval
  !> Physical point location
  PROCEDURE :: wrap => nurbs_surf_wrap
  !> Parametric point evaluation
  PROCEDURE :: unwrap => nurbs_surf_unwrap
END TYPE nurbs_surf
!---------------------------------------------------------------------------
! TYPE nurbs_entity_ptr
!---------------------------------------------------------------------------
!> List of CAD entities
!---------------------------------------------------------------------------
TYPE :: nurbs_entity_ptr
  CLASS(nurbs_entity), POINTER :: wo => NULL() !< OpenNURBS CAD entity
END TYPE nurbs_entity_ptr
!---------------------------------------------------------------------------
! Abstract subroutine definitions
!---------------------------------------------------------------------------
ABSTRACT INTERFACE
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_dummy_eval
!---------------------------------------------------------------------------
!> Evaluate a NURBS object
!!
!! @param[out] pt Evaluated possition [3]
!! @param[in] u Parametric coordinate 1
!! @param[in] v Parametric coordinate 2 (unused if curve)
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_dummy_eval(self,pt,u,v)
    IMPORT nurbs_entity, r8
    CLASS(nurbs_entity), INTENT(in) :: self
    REAL(r8), INTENT(out) :: pt(3)
    REAL(r8), INTENT(in) :: u
    REAL(r8), INTENT(in) :: v
  END SUBROUTINE nurbs_dummy_eval
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_dummy_grid
!---------------------------------------------------------------------------
!> Evalute a grid of points evenly spaced in parametric space
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_dummy_grid(self)
    IMPORT nurbs_entity
    CLASS(nurbs_entity), INTENT(inout) :: self
  END SUBROUTINE nurbs_dummy_grid
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_dummy_find
!---------------------------------------------------------------------------
!> Find the nearest location on a NURBS object to a given point
!!
!! @param[in] pt Evaluated possition [3]
!! @param[out] u Parametric coordinate 1
!! @param[out] v Parametric coordinate 2 (unused if curve)
!! @param[out] ierr Error status (point not found if ierr<0)
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_dummy_find(self,pt,u,v,ierr)
    IMPORT nurbs_entity, i4, r8
    CLASS(nurbs_entity), TARGET, INTENT(in) :: self
    REAL(r8), INTENT(in) :: pt(3)
    REAL(r8), INTENT(out) :: u
    REAL(r8), INTENT(out) :: v
    INTEGER(i4), INTENT(out) :: ierr
  END SUBROUTINE nurbs_dummy_find
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_dummy_wrap
!---------------------------------------------------------------------------
!> Wrap the periodic coordinates of a NURBS object into a single span
!! - Example: For a periodic domain of \f$ [0,2 \pi] (-\pi/4,0) -> (3 \pi/4,0) \f$
!!
!! @param[in,out] u Parametric coordinate 1
!! @param[in,out] v Parametric coordinate 2 (unused if curve)
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_dummy_wrap(self,u,v)
    IMPORT nurbs_entity, r8
    CLASS(nurbs_entity), INTENT(in) :: self
    REAL(r8), INTENT(inout) :: u
    REAL(r8), INTENT(inout) :: v
  END SUBROUTINE nurbs_dummy_wrap
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_dummy_unwrap
!---------------------------------------------------------------------------
!> Unwrap the periodic coordinates of a NURBS object to avoid parameter cuts
!! - Example: For a periodic domain of \f$ [0,2 \pi] (3 \pi/4,0) -> (-\pi/4,0) \f$
!!
!! @param[in,out] u1 Parametric coordinate 1 of first point
!! @param[in,out] v1 Parametric coordinate 2 of first point (unused if curve)
!! @param[in,out] u2 Parametric coordinate 1 of second point
!! @param[in,out] v2 Parametric coordinate 2 of second point (unused if curve)
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_dummy_unwrap(self,u1,u2,v1,v2)
    IMPORT nurbs_entity, r8
    CLASS(nurbs_entity), INTENT(in) :: self
    REAL(r8), INTENT(inout) :: u1
    REAL(r8), INTENT(inout) :: u2
    REAL(r8), INTENT(inout) :: v1
    REAL(r8), INTENT(inout) :: v2
  END SUBROUTINE nurbs_dummy_unwrap
END INTERFACE
!---------------------------------------------------------------------------
! Interfaces to C++ NURBS subroutines
!---------------------------------------------------------------------------
INTERFACE
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_init
!---------------------------------------------------------------------------
!> Initialize the OpenNURBS library
!!
!! @param[in,out] ierr Error status
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_init(ierr) BIND(C)
    IMPORT c_int
    INTEGER(c_int), INTENT(inout) :: ierr
  END SUBROUTINE nurbs_init
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_finalize
!---------------------------------------------------------------------------
!> Finalize the OpenNURBS library
!!
!! @param[in,out] ierr Error status
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_finalize(ierr) BIND(C)
    IMPORT c_int
    INTEGER(c_int), INTENT(inout) :: ierr
  END SUBROUTINE nurbs_finalize
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_read_in
!---------------------------------------------------------------------------
!> Read-in a OpenNURBS *.3dm geometry file set objects as active
!!
!! @param[in] filename Name of geometry file
!! @param[in,out] ierr Error status
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_read_in(filename,ierr) BIND(C)
    IMPORT c_int, c_char
    CHARACTER(KIND=c_char), INTENT(in) :: filename
    INTEGER(c_int), INTENT(inout) :: ierr
  END SUBROUTINE nurbs_read_in
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_get_count
!---------------------------------------------------------------------------
!> Get a count of OpenNURBS objects in current file
!!
!! @param[in,out] nc Number of NURBS curves
!! @param[in,out] ns Number of NURBS surfaces
!! @param[in,out] ierr Error status
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_get_count(nc,ns,ierr) BIND(C)
    IMPORT c_int
    INTEGER(c_int), INTENT(inout) :: nc
    INTEGER(c_int), INTENT(inout) :: ns
    INTEGER(c_int), INTENT(inout) :: ierr
  END SUBROUTINE nurbs_get_count
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_curve_name
!---------------------------------------------------------------------------
!> Get the name of a OpenNURBS curve
!!
!! @param[in] ind Index of curve
!! @param[in,out] name Name of curve
!! @param[in,out] ierr Error status
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_curve_name(ind,name,ierr) BIND(C)
    IMPORT c_int, c_char
    INTEGER(c_int), INTENT(in) :: ind
    CHARACTER(KIND=c_char), INTENT(inout) :: name
    INTEGER(c_int), INTENT(inout) :: ierr
  END SUBROUTINE nurbs_curve_name
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_eval_curve
!---------------------------------------------------------------------------
!> Evalute the possition along a NURBS curve
!!
!! @param[in] ind Index of curve
!! @param[in] u Parametric possition
!! @param[in,out] r Possition in physical space [3]
!! @param[in,out] ierr Error status
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_eval_curve(ind,u,r,ierr) BIND(C)
    IMPORT c_int, c_double
    INTEGER(c_int), INTENT(in) :: ind
    REAL(c_double), INTENT(in) :: u
    REAL(c_double), DIMENSION(3), INTENT(inout) :: r
    INTEGER(c_int), INTENT(inout) :: ierr
  END SUBROUTINE nurbs_eval_curve
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_curve_domain
!---------------------------------------------------------------------------
!> Get the parametric domain of a NURBS curve
!!
!! @param[in] ind Index of curve
!! @param[in,out] domain Domain of parametric variable [2]
!! @param[in,out] ierr Error status
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_curve_domain(ind,domain,ierr) BIND(C)
    IMPORT c_int, c_double
    INTEGER(c_int), INTENT(in) :: ind
    REAL(c_double), DIMENSION(2), INTENT(inout) :: domain
    INTEGER(c_int), INTENT(inout) :: ierr
  END SUBROUTINE nurbs_curve_domain
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_curve_periodic
!---------------------------------------------------------------------------
!> Test periodicity of a NURBS curve
!!
!! @param[in] ind Index of curve
!! @param[in,out] p Periodic flag (1 if true, 0 if false)
!! @param[in,out] ierr Error status
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_curve_periodic(ind,p,ierr) BIND(C)
    IMPORT c_int
    INTEGER(c_int), INTENT(in) :: ind
    INTEGER(c_int), INTENT(inout) :: p
    INTEGER(c_int), INTENT(inout) :: ierr
  END SUBROUTINE nurbs_curve_periodic
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_curve_linear
!---------------------------------------------------------------------------
!> Test if NURBS curve is straight
!!
!! @param[in] ind Index of curve
!! @param[in,out] p Straight flag (1 if true, 0 if false)
!! @param[in,out] ierr Error status
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_curve_linear(ind,p,ierr) BIND(C)
    IMPORT c_int
    INTEGER(c_int), INTENT(in) :: ind
    INTEGER(c_int), INTENT(inout) :: p
    INTEGER(c_int), INTENT(inout) :: ierr
  END SUBROUTINE nurbs_curve_linear
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_surf_name
!---------------------------------------------------------------------------
!> Get the name of a OpenNURBS surface
!!
!! @param[in] ind Index of surface
!! @param[in,out] name Name of surface
!! @param[in,out] ierr Error status
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_surf_name(ind,name,ierr) BIND(C)
    IMPORT c_int, c_char
    INTEGER(c_int), INTENT(in) :: ind
    CHARACTER(KIND=c_char), INTENT(inout) :: name
    INTEGER(c_int), INTENT(inout) :: ierr
  END SUBROUTINE nurbs_surf_name
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_eval_surf
!---------------------------------------------------------------------------
!> Evalute the possition on a NURBS surface
!!
!! @param[in] ind Index of surface
!! @param[in] u Parametric possition 1
!! @param[in] v Parametric possition 2
!! @param[in,out] r Possition in physical space [3]
!! @param[in,out] ierr Error status
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_eval_surf(ind,u,v,r,ierr) BIND(C)
    IMPORT c_int, c_double
    INTEGER(c_int), INTENT(in) :: ind
    REAL(c_double), INTENT(in) :: u
    REAL(c_double), INTENT(in) :: v
    REAL(c_double), DIMENSION(3), INTENT(inout) :: r
    INTEGER(c_int), INTENT(inout) :: ierr
  END SUBROUTINE nurbs_eval_surf
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_surf_domain
!---------------------------------------------------------------------------
!> Get the parametric domain of a NURBS surface
!!
!! @param[in] ind Index of surface
!! @param[in,out] domain1 Domain of parametric variable 1 [2]
!! @param[in,out] domain2 Domain of parametric variable 2 [2]
!! @param[in,out] ierr Error status
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_surf_domain(ind,domain1,domain2,ierr) BIND(C)
    IMPORT c_int, c_double
    INTEGER(c_int), INTENT(in) :: ind
    REAL(c_double), DIMENSION(2), INTENT(inout) :: domain1
    REAL(c_double), DIMENSION(2), INTENT(inout) :: domain2
    INTEGER(c_int), INTENT(inout) :: ierr
  END SUBROUTINE nurbs_surf_domain
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_surf_periodic
!---------------------------------------------------------------------------
!> Test periodicity of a NURBS surface
!!
!! @param[in] ind Index of surface
!! @param[in,out] p1 Periodic flag for direction 1 (1 if true, 0 if false)
!! @param[in,out] p2 Periodic flag for direction 2 (1 if true, 0 if false)
!! @param[in,out] ierr Error status
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_surf_periodic(ind,p1,p2,ierr) BIND(C)
    IMPORT c_int
    INTEGER(c_int), INTENT(in) :: ind
    INTEGER(c_int), INTENT(inout) :: p1
    INTEGER(c_int), INTENT(inout) :: p2
    INTEGER(c_int), INTENT(inout) :: ierr
  END SUBROUTINE nurbs_surf_periodic
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_surf_singular
!---------------------------------------------------------------------------
!> Test if edges of a NURBS surface are degenerate
!!
!! @param[in] ind Index of surface
!! @param[in,out] p Singular flag (1 if true, 0 if false) [4]
!! @param[in,out] ierr Error status
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_surf_singular(ind,p,ierr) BIND(C)
    IMPORT c_int
    INTEGER(c_int), INTENT(in) :: ind
    INTEGER(c_int), DIMENSION(4), INTENT(inout) :: p
    INTEGER(c_int), INTENT(inout) :: ierr
  END SUBROUTINE nurbs_surf_singular
!---------------------------------------------------------------------------
! SUBROUTINE nurbs_surf_planar
!---------------------------------------------------------------------------
!> Test if NURBS surface is planar
!!
!! @param[in] ind Index of surface
!! @param[in,out] p Planar flag (1 if true, 0 if false)
!! @param[in,out] ierr Error status
!---------------------------------------------------------------------------
  SUBROUTINE nurbs_surf_planar(ind,p,ierr) BIND(C)
    IMPORT c_int
    INTEGER(c_int), INTENT(in) :: ind
    INTEGER(c_int), INTENT(inout) :: p
    INTEGER(c_int), INTENT(inout) :: ierr
  END SUBROUTINE nurbs_surf_planar
END INTERFACE
!---------------------------------------------------------------------------
! Declare public interface
!---------------------------------------------------------------------------
PUBLIC nurbs_entity, nurbs_curve, nurbs_surf, nurbs_entity_ptr, &
  nurbs_init, nurbs_finalize, nurbs_read_in, nurbs_get_count, &
  nurbs_curve_name, nurbs_eval_curve, nurbs_curve_domain, &
  nurbs_curve_periodic, nurbs_surf_name, nurbs_eval_surf, &
  nurbs_surf_domain, nurbs_surf_periodic, nurbs_surf_singular, &
  nurbs_curve_linear, nurbs_surf_planar
!---Fitting variables
CLASS(nurbs_curve), POINTER :: active_curve => NULL() !< Active curve for MINPACK fitting
CLASS(nurbs_surf), POINTER :: active_surf => NULL() !< Active surface for MINPACK fitting
real(r8) :: active_pt(3) !< Active fit point for MINPACK fitting
real(r8) :: active_endpts(3,3) !< Active constraint points for MINPACK fitting
real(r8) :: active_wts(3) !< Active constraint weights for MINPACK fitting
!$omp threadprivate(active_pt,active_endpts,active_wts,active_curve,active_surf)
PUBLIC nurbs_surf_avg, nurbs_surf_midpoint, nurbs_curve_midpoint, nurbs_surf_center
CONTAINS
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_curve_eval
!---------------------------------------------------------------------------
!> Map the parametric possition on a curve to physical coordinates.
!! - (u,v) -> (x,y,z)
!!
!! @param[out] pt Position vector [3]
!! @param[in] u Parametric coordinate 1
!! @param[in] v Parametric coordinate 2 (ignored)
!---------------------------------------------------------------------------
subroutine nurbs_curve_eval(self,pt,u,v)
class(nurbs_curve), intent(in) :: self
real(r8), intent(out) :: pt(3)
real(r8), intent(in) :: u,v
real(r8) :: utmp
integer(i4) :: ierr
DEBUG_STACK_PUSH
utmp=u
if(self%periodic)then
  if(u>self%domain(2))utmp=u-self%domain(2)+self%domain(1)
  if(u<self%domain(1))utmp=u-self%domain(1)+self%domain(2)
end if
!---Initialize point reconstruction
call nurbs_eval_curve(self%id,utmp,pt,ierr)
pt=pt*self%coord_scales
IF(self%reflect>0)pt(self%reflect)=-pt(self%reflect)
DEBUG_STACK_POP
end subroutine nurbs_curve_eval
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_curve_find
!---------------------------------------------------------------------------
!> Find the parametric representation of a boundary point in CAD representation on a curve.
!! - (x,y,z) -> (u,v)
!!
!! @param[in] pt Position vector [3]
!! @param[out] u Parametric coordinate 1
!! @param[out] v Parametric coordinate 2 (ignored)
!---------------------------------------------------------------------------
subroutine nurbs_curve_find(self,pt,u,v,ierr)
class(nurbs_curve), target, intent(in) :: self
real(r8), intent(in) :: pt(3)
real(r8), intent(out) :: u,v
integer(i4), intent(out) :: ierr
real(r8) :: du,un,val(3),rt(3),fp,ug,dmin
integer :: i
DEBUG_STACK_PUSH
ierr=0
ug=u
v=0.d0
!---Check if point is vertex 1
u=self%domain(1)
call self%eval(rt,u,v)
val(2)=sqrt(sum((pt-rt)**2))
IF(val(2)<1.d-12)THEN
  DEBUG_STACK_POP
  return ! Point is vertex 1
END IF
!---Check if point is vertex 2
u=self%domain(2)
call self%eval(rt,u,v)
val(2)=sqrt(sum((pt-rt)**2))
if(val(2)<1.d-12)THEN
  DEBUG_STACK_POP
  return ! Point is vertex 2
END IF
!---Initialize search point
!if(ug<self%domain(1).OR.ug>self%domain(2))then
!  u=(self%domain(1)+self%domain(2))/2.d0
!else
!  u=ug
!end if
dmin=1.d99
do i=1,nurbs_ngrid
  rt=self%rgrid(:,i)
  val(2)=sqrt(sum((pt-rt)**2))
  if(val(2)<dmin)then
    dmin=val(2)
    u=(i-1)*(self%domain(2)-self%domain(1))*(1.d0/(nurbs_ngrid-1))+self%domain(1)
  end if
end do
un=u
fp=1.d99
du=1.d-9
!---Use Newton's Method to find point
do i=1,40
  do while(.TRUE.)
    call self%eval(rt,un,v)
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
  call self%eval(rt,u-du,v)
  val(1)=sqrt(sum((pt-rt)**2))
  call self%eval(rt,u+du,v)
  val(3)=sqrt(sum((pt-rt)**2))
  !---Prevent divide by zero
  if(abs(val(3)-val(1))<1.d-14)then
    u=u-du/2
    cycle
  endif
  !---Update solution
  un=u-val(2)*du*2/(val(3)-val(1))
  !---Limit minimum and maximum values
  if(.NOT.self%periodic)then
    if(un>self%domain(2))un=self%domain(2)-.01d0*(self%domain(2)-self%domain(1))
    if(un<self%domain(1))un=self%domain(1)-.01d0*(self%domain(1)-self%domain(2))
  end if
enddo
!---Catch failure to converge
if(i>40)then
  ierr=-1
endif
DEBUG_STACK_POP
end subroutine nurbs_curve_find
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_curve_midpoint
!---------------------------------------------------------------------------
!> Compute the weighted midpoint of a curve edge
!!
!! Locates the point on a given CAD curve which minimizes the weighted sum
!! of distances to 2 constraint points.
!!
!! \f[ \sum_i w_i*(r_n - p_i)^2 \f]
!!
!! @param[in] self CAD curve
!! @param[in,out] pt Solution point
!! @param[in] pt1 Constraint point 1
!! @param[in] pt2 Constraint point 2
!! @param[in] wt1 Constraint weight 1
!! @param[in] wt2 Constraint weight 2
!! @param[out] ierr Error flag
!---------------------------------------------------------------------------
subroutine nurbs_curve_midpoint(self,pt,pt1,pt2,wt1,wt2,ierr)
class(nurbs_curve),target, intent(in) :: self
real(r8), intent(inout) :: pt(3)
real(r8), intent(in) :: pt1(3)
real(r8), intent(in) :: pt2(3)
real(r8), intent(in) :: wt1
real(r8), intent(in) :: wt2
integer(i4), intent(out) :: ierr
real(r8) :: u,v,val(2),rt(3)
integer :: i
!---
integer(i4), parameter :: nerr=6
integer(i4), parameter :: neq=1
real(r8) :: diag(neq),fjac(nerr,neq),qtf(neq),uv(neq)
real(r8) :: wa1(neq),wa2(neq),wa3(neq),wa4(nerr),error(nerr)
integer(i4) :: ipvt(neq)
real(r8) :: ftol,xtol,gtol,epsfcn,factor,dmin,err_in,uv_in(neq)
integer :: info,ldfjac,maxfev,mode,nfev,nprint
DEBUG_STACK_PUSH
ierr=0
v=0.d0
!---
dmin=1.d99
do i=1,nurbs_ngrid
  rt=self%rgrid(:,i)
  val(2)=sqrt(sum((pt-rt)**2))
  if(val(2)<dmin)then
    dmin=val(2)
    u=(i-1)*(self%domain(2)-self%domain(1))*(1.d0/(nurbs_ngrid-1))+self%domain(1)
  end if
end do
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
active_curve=>self
active_endpts(1,:)=pt1; active_wts(1)=wt1
active_endpts(2,:)=pt2; active_wts(2)=wt2
uv=(/u/)
uv_in=uv
CALL nurbs_cmid_error(nerr,neq,uv,error,info)
err_in = SQRT(SUM(error**2))
!---
call lmdif(nurbs_cmid_error,nerr,neq,uv,error, &
           ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
           nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!---
IF(SQRT(SUM(error**2))>err_in)uv=uv_in
IF(info>4)THEN
  !WRITE(*,*)'Curve midpoint failed',self%cid,info,nfev
  ierr=-1
ELSE
  ierr=0
END IF
u=uv(1)
v=0.d0
call self%eval(pt,u,v)
DEBUG_STACK_POP
end subroutine nurbs_curve_midpoint
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_mid_error
!---------------------------------------------------------------------------
!> Evalute the error between a curve point and the current active points
!! used in a 2 point minimization.
!!
!! @note Designed to be used as the error function for minimization in
!! @ref nurbs_cad::nurbs_curve_midpoint "nurbs_curve_midpoint"
!!
!! @param[in] m Number of spatial dimensions (3)
!! @param[in] n Number of parametric dimensions (2)
!! @param[in] uv Parametric possition [n]
!! @param[out] err Error vector between current and desired point [3]
!! @param[in,out] iflag Unused flag
!---------------------------------------------------------------------------
subroutine nurbs_cmid_error(m,n,uv,err,iflag)
integer(i4), intent(in) :: m,n
real(r8), intent(in) :: uv(n)
real(r8), intent(out) :: err(m)
integer(i4), intent(inout) :: iflag
real(r8) :: pt(3)
DEBUG_STACK_PUSH
call active_curve%eval(pt,uv(1),0.d0)
err(1:3)=(active_endpts(1,:)-pt)*SQRT(active_wts(1))
err(4:6)=(active_endpts(2,:)-pt)*SQRT(active_wts(2))
DEBUG_STACK_POP
end subroutine nurbs_cmid_error
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_curve_unwrap
!---------------------------------------------------------------------------
!> Unwrap the coordinate of a periodic NURBS curve to avoid parameter cuts
!! - Example: For a periodic domain of \f$ [0,2 \pi] (3 \pi/4,0) -> (-\pi/4,0) \f$
!!
!! @param[in,out] u1 Parametric coordinate 1 of first point
!! @param[in,out] v1 Parametric coordinate 2 of first point (unused)
!! @param[in,out] u2 Parametric coordinate 1 of second point
!! @param[in,out] v2 Parametric coordinate 2 of second point (unused)
!---------------------------------------------------------------------------
subroutine nurbs_curve_unwrap(self,u1,u2,v1,v2)
class(nurbs_curve), intent(in) :: self
real(r8), intent(inout) :: u1,u2,v1,v2
real(r8) :: ubound,lbound,per_buff=.2d0
DEBUG_STACK_PUSH
!---Limit minimum and maximum values
if(self%periodic)then
  ubound = (1.d0-per_buff)*self%domain(2)+per_buff*self%domain(1)
  lbound = per_buff*self%domain(2)+(1.d0-per_buff)*self%domain(1)
  if(u1<lbound.AND.u2>ubound)u2=u2-self%domain(2)+self%domain(1)
  if(u2<lbound.AND.u1>ubound)u1=u1-self%domain(2)+self%domain(1)
end if
DEBUG_STACK_POP
end subroutine nurbs_curve_unwrap
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_curve_wrap
!---------------------------------------------------------------------------
!> Wrap the coordinate of a periodic NURBS curve into a single span
!! - Example: For a periodic domain of \f$ [0,2 \pi] (-\pi/4,0) -> (3 \pi/4,0) \f$
!!
!! @param[in,out] u Parametric coordinate 1
!! @param[in,out] v Parametric coordinate 2 (unused)
!---------------------------------------------------------------------------
subroutine nurbs_curve_wrap(self,u,v)
class(nurbs_curve), intent(in) :: self
real(r8), intent(inout) :: u,v
real(r8) :: ubound,lbound
DEBUG_STACK_PUSH
!---Limit minimum and maximum values
if(self%periodic)then
  ubound = .9d0*self%domain(2)+.1d0*self%domain(1)
  lbound = .1d0*self%domain(2)+.9d0*self%domain(1)
  if(u<lbound)u=ubound+lbound-u
end if
DEBUG_STACK_POP
end subroutine nurbs_curve_wrap
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_curve_grid
!---------------------------------------------------------------------------
!> Evalute a grid of points evenly spaced in parametric space
!---------------------------------------------------------------------------
subroutine nurbs_curve_grid(self)
class(nurbs_curve), intent(inout) :: self
real(r8) :: u,v
integer(i4) :: i
DEBUG_STACK_PUSH
v=0.d0
do i=1,nurbs_ngrid
  u=(i-1)*(self%domain(2)-self%domain(1))*(1.d0/(nurbs_ngrid-1))+self%domain(1)
  call self%eval(self%rgrid(:,i),u,v)
end do
DEBUG_STACK_POP
end subroutine nurbs_curve_grid
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_surf_eval
!---------------------------------------------------------------------------
!> Map the parametric possition on a surface to physical coordinates.
!! - (u,v) -> (x,y,z)
!!
!! @param[out] pt Position vector [3]
!! @param[in] u Parametric coordinate 1
!! @param[in] v Parametric coordinate 2
!---------------------------------------------------------------------------
subroutine nurbs_surf_eval(self,pt,u,v)
class(nurbs_surf), intent(in) :: self
real(r8), intent(out) :: pt(3)
real(r8), intent(in) :: u,v
real(r8) :: utmp,vtmp,span
integer(i4) :: ierr,n
DEBUG_STACK_PUSH
utmp=u
vtmp=v
if(self%periodic(1))then
  span=self%domain(2,1)-self%domain(1,1)
  if(u>self%domain(2,1))then
    utmp=u-self%domain(2,1)
    n=floor(utmp/span)
    utmp=u-n*span+self%domain(1,1)
  end if
  if(u<self%domain(1,1))then
    utmp=u-self%domain(1,1)
    n=floor(utmp/span)
    utmp=u-n*span+self%domain(1,1)
  end if
end if
if(self%periodic(2))then
  !if(v>self%domain(2,2))vtmp=v-self%domain(2,2)+self%domain(1,2)
  !if(v<self%domain(1,2))vtmp=v-self%domain(1,2)+self%domain(2,2)
  span=self%domain(2,2)-self%domain(1,2)
  if(v>self%domain(2,2))then
    vtmp=v-self%domain(2,2)
    n=floor(vtmp/span)
    vtmp=v-n*span+self%domain(1,2)
  end if
  if(v<self%domain(1,2))then
    vtmp=v-self%domain(1,2)
    n=floor(vtmp/span)
    vtmp=v-n*span+self%domain(1,2)
  end if
end if
!---Initialize point reconstruction
call nurbs_eval_surf(self%id,utmp,vtmp,pt,ierr)
pt=pt*self%coord_scales
IF(self%reflect>0)pt(self%reflect)=-pt(self%reflect)
DEBUG_STACK_POP
end subroutine nurbs_surf_eval
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_surf_find
!---------------------------------------------------------------------------
!> Find the parametric representation of a boundary point in CAD representation on a surface.
!! - (x,y,z) -> (u,v)
!!
!! @param[in] pt Position vector [3]
!! @param[out] u Parametric coordinate 1
!! @param[out] v Parametric coordinate 2
!---------------------------------------------------------------------------
subroutine nurbs_surf_find(self,pt,u,v,ierr)
class(nurbs_surf), target, intent(in) :: self
real(r8), intent(in) :: pt(3)
real(r8), intent(out) :: u,v
integer(i4), intent(out) :: ierr
logical :: sing(2)
!---
integer(i4), parameter :: nerr=3
integer(i4), parameter :: neq=2
real(r8) :: diag(neq),fjac(nerr,neq),qtf(neq),uv(neq)
real(r8) :: wa1(neq),wa2(neq),wa3(neq),wa4(nerr),error(nerr)
integer(i4) :: ipvt(neq)
real(r8) :: ftol,xtol,gtol,epsfcn,factor
real(r8) :: du,dv,un,vn,val(3),vals(4),rt(3),grad(2),fp,ug,vg,dmin
integer :: i,j,step(2),info,ldfjac,maxfev,mode,nfev,nprint
DEBUG_STACK_PUSH
!---
sing=.FALSE.
ierr=0
ug=u; vg=v
!---Check if point is vertex 1
u=self%domain(1,1); v=self%domain(1,2)
call self%eval(rt,u,v)
val(2)=sqrt(sum((pt-rt)**2))
if(val(2)<1.d-12)THEN
  DEBUG_STACK_POP
  RETURN ! Point is vertex 1
END IF
!---Check if point is vertex 2
u=self%domain(2,1); v=self%domain(1,2)
call self%eval(rt,u,v)
val(2)=sqrt(sum((pt-rt)**2))
if(val(2)<1.d-12)THEN
  DEBUG_STACK_POP
  RETURN ! Point is vertex 2
END IF
!---Check if point is vertex 3
u=self%domain(2,1); v=self%domain(2,2)
call self%eval(rt,u,v)
val(2)=sqrt(sum((pt-rt)**2))
if(val(2)<1.d-12)THEN
  DEBUG_STACK_POP
  RETURN ! Point is vertex 3
END IF
!---Check if point is vertex 4
u=self%domain(1,1); v=self%domain(2,2)
call self%eval(rt,u,v)
val(2)=sqrt(sum((pt-rt)**2))
if(val(2)<1.d-12)THEN
  DEBUG_STACK_POP
  RETURN ! Point is vertex 4
END IF
!---
dmin=1.d99
do i=1,nurbs_ngrid
  do j=1,nurbs_ngrid
    rt=self%rgrid(:,i,j)
    val(2)=sqrt(sum((pt-rt)**2))
    if(val(2)<dmin)then
      dmin=val(2)
      u=(i-1)*(self%domain(2,1)-self%domain(1,1))*(1.d0/(nurbs_ngrid-1))+self%domain(1,1)
      v=(j-1)*(self%domain(2,2)-self%domain(1,2))*(1.d0/(nurbs_ngrid-1))+self%domain(1,2)
    end if
  end do
end do
sing(1)=nurbs_surf_atsingular(self,u,v)
IF(sing(1))THEN
  u=(1.d0-nurbs_singstep)*u+nurbs_singstep*(self%domain(2,1)+self%domain(1,1))/2.d0
  v=(1.d0-nurbs_singstep)*v+nurbs_singstep*(self%domain(2,2)+self%domain(1,2))/2.d0
END IF
call self%eval(rt,u,v)
!write(*,*)0,u,v,dmin
!write(*,*)rt
un=u
vn=v
!---
mode = 1
factor = 1.d0
maxfev = 100
ftol = 1.d-6
xtol = 1.d-6
gtol = 1.d-3
epsfcn = 1.d-6
nprint = 0
ldfjac = 3

active_surf=>self
active_pt=pt
uv=(/u,v/)
call lmdif(nurbs_surf_error,3,2,uv,error, &
         ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
         nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
IF(info>4)THEN
  !WRITE(*,*)'Surface midpoint failed',self%sid,info,nfev
  ierr=-1
ELSE
  ierr=0
END IF
u=uv(1)
v=uv(2)
call self%eval(rt,u,v)
val(2)=sqrt(sum((pt-rt)**2))
DEBUG_STACK_POP
end subroutine nurbs_surf_find
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_surf_midpoint
!---------------------------------------------------------------------------
!> Compute the weighted midpoint of a surface edge
!!
!! Locates the point on a given CAD surface which minimizes the weighted sum
!! of distances to 2 constraint points.
!!
!! \f[ \sum_i w_i*(r_n - p_i)^2 \f]
!!
!! @param[in] self CAD surface
!! @param[in,out] pt Solution point
!! @param[in] pt1 Constraint point 1
!! @param[in] pt2 Constraint point 2
!! @param[in] wt1 Constraint weight 1
!! @param[in] wt2 Constraint weight 2
!! @param[out] ierr Error flag
!---------------------------------------------------------------------------
subroutine nurbs_surf_midpoint(self,pt,pt1,pt2,wt1,wt2,ierr)
class(nurbs_surf), target, intent(in) :: self
real(r8), intent(inout) :: pt(3)
real(r8), intent(in) :: pt1(3)
real(r8), intent(in) :: pt2(3)
real(r8), intent(in) :: wt1
real(r8), intent(in) :: wt2
integer(i4), intent(out) :: ierr
logical :: sing
real(r8) :: u,v
!---
integer(i4), parameter :: nerr=6
integer(i4), parameter :: neq=2
real(r8) :: diag(neq),fjac(nerr,neq),qtf(neq),uv(neq)
real(r8) :: wa1(neq),wa2(neq),wa3(neq),wa4(nerr),error(nerr)
integer(i4) :: ipvt(neq)
real(r8) :: ftol,xtol,gtol,epsfcn,factor,err_in,uv_in(neq)
real(r8) :: val(3),rt(3),dmin
integer :: i,j,info,ldfjac,maxfev,mode,nfev,nprint
DEBUG_STACK_PUSH
!---
sing=.FALSE.
ierr=0
!---
dmin=1.d99
do i=1,nurbs_ngrid
  do j=1,nurbs_ngrid
    rt=self%rgrid(:,i,j)
    val(2)=sqrt(sum((pt-rt)**2))
    if(val(2)<dmin)then
      dmin=val(2)
      u=(i-1)*(self%domain(2,1)-self%domain(1,1))*(1.d0/(nurbs_ngrid-1))+self%domain(1,1)
      v=(j-1)*(self%domain(2,2)-self%domain(1,2))*(1.d0/(nurbs_ngrid-1))+self%domain(1,2)
    end if
  end do
end do
sing=nurbs_surf_atsingular(self,u,v)
IF(sing)THEN
  u=(1.d0-nurbs_singstep)*u+nurbs_singstep*(self%domain(2,1)+self%domain(1,1))/2.d0
  v=(1.d0-nurbs_singstep)*v+nurbs_singstep*(self%domain(2,2)+self%domain(1,2))/2.d0
END IF
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
uv=(/u,v/)
uv_in=uv
CALL nurbs_smid_error(nerr,neq,uv_in,error,info)
err_in = SQRT(SUM(error**2))
!---
call lmdif(nurbs_smid_error,nerr,neq,uv,error, &
           ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
           nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!---
IF(SQRT(SUM(error**2))>err_in)uv=uv_in
IF(info>4)THEN
  !WRITE(*,*)'Surface midpoint failed',self%sid,info,nfev
  ierr=-1
ELSE
  ierr=0
END IF
!WRITE(*,*)'surf',info
u=uv(1)
v=uv(2)
call self%eval(rt,u,v)
pt=rt
DEBUG_STACK_POP
end subroutine nurbs_surf_midpoint
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_smid_error
!---------------------------------------------------------------------------
!> Evalute the error between a surface point and the current active points
!! used in a 2 point minimization.
!!
!! @note Designed to be used as the error function for minimization in
!! @ref nurbs_cad::nurbs_surf_midpoint "nurbs_surf_midpoint"
!!
!! @param[in] m Number of spatial dimensions (3)
!! @param[in] n Number of parametric dimensions (2)
!! @param[in] uv Parametric possition [n]
!! @param[out] err Error vector between current and desired point [3]
!! @param[in,out] iflag Unused flag
!---------------------------------------------------------------------------
subroutine nurbs_smid_error(m,n,uv,err,iflag)
integer(i4), intent(in) :: m,n
real(r8), intent(in) :: uv(n)
real(r8), intent(out) :: err(m)
integer(i4), intent(inout) :: iflag
real(r8) :: pt(3)
DEBUG_STACK_PUSH
call active_surf%eval(pt,uv(1),uv(2))
err(1:3)=(active_endpts(1,:)-pt)*SQRT(active_wts(1))
err(4:6)=(active_endpts(2,:)-pt)*SQRT(active_wts(2))
DEBUG_STACK_POP
end subroutine nurbs_smid_error
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_surf_center
!---------------------------------------------------------------------------
!> Compute the weighted center point of a surface triangle
!!
!! Locates the point on a given CAD surface which minimizes the weighted sum
!! of distances to 3 constraint points.
!!
!! \f[ \sum_i w_i*(r_n - p_i)^2 \f]
!!
!! @param[in] self CAD surface
!! @param[in,out] pt Solution point
!! @param[in] pt1 Constraint point 1
!! @param[in] pt2 Constraint point 2
!! @param[in] pt3 Constraint point 3
!! @param[in] wt1 Constraint weight 1
!! @param[in] wt2 Constraint weight 2
!! @param[in] wt3 Constraint weight 3
!! @param[out] ierr Error flag
!---------------------------------------------------------------------------
subroutine nurbs_surf_center(self,pt,pt1,pt2,pt3,wt1,wt2,wt3,ierr)
class(nurbs_surf), target, intent(in) :: self
real(r8), intent(inout) :: pt(3)
real(r8), intent(in) :: pt1(3)
real(r8), intent(in) :: pt2(3)
real(r8), intent(in) :: pt3(3)
real(r8), intent(in) :: wt1
real(r8), intent(in) :: wt2
real(r8), intent(in) :: wt3
integer(i4), intent(out) :: ierr
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
!---
dmin=1.d99
do i=1,nurbs_ngrid
  do j=1,nurbs_ngrid
    rt=self%rgrid(:,i,j)
    val(2)=sqrt(sum((pt-rt)**2))
    if(val(2)<dmin)then
      dmin=val(2)
      u=(i-1)*(self%domain(2,1)-self%domain(1,1))*(1.d0/(nurbs_ngrid-1))+self%domain(1,1)
      v=(j-1)*(self%domain(2,2)-self%domain(1,2))*(1.d0/(nurbs_ngrid-1))+self%domain(1,2)
    end if
  end do
end do
sing=nurbs_surf_atsingular(self,u,v)
IF(sing)THEN
  u=(1.d0-nurbs_singstep)*u+nurbs_singstep*(self%domain(2,1)+self%domain(1,1))/2.d0
  v=(1.d0-nurbs_singstep)*v+nurbs_singstep*(self%domain(2,2)+self%domain(1,2))/2.d0
END IF
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
call lmdif(nurbs_scenter_error,nerr,neq,uv,error, &
           ftol,xtol,gtol,maxfev,epsfcn,diag,mode,factor,nprint,info, &
           nfev,fjac,ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
!---
IF(info>4)THEN
  !WRITE(*,*)'Surface center point failed',self%sid,info,nfev
  ierr=-1
ELSE
  ierr=0
END IF
!WRITE(*,*)'surf',info
u=uv(1)
v=uv(2)
call self%eval(pt,u,v)
DEBUG_STACK_POP
end subroutine nurbs_surf_center
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_scenter_error
!---------------------------------------------------------------------------
!> Evalute the error between a surface point and the current active points
!! used in a 3 point minimization.
!!
!! @note Designed to be used as the error function for minimization in
!! @ref nurbs_cad::nurbs_surf_center "nurbs_surf_center"
!!
!! @param[in] m Number of spatial dimensions (3)
!! @param[in] n Number of parametric dimensions (2)
!! @param[in] uv Parametric possition [n]
!! @param[out] err Error vector between current and desired point [3]
!! @param[in,out] iflag Unused flag
!---------------------------------------------------------------------------
subroutine nurbs_scenter_error(m,n,uv,err,iflag)
integer(i4), intent(in) :: m,n
real(r8), intent(in) :: uv(n)
real(r8), intent(out) :: err(m)
integer(i4), intent(inout) :: iflag
real(r8) :: pt(3)
DEBUG_STACK_PUSH
call active_surf%eval(pt,uv(1),uv(2))
err(1:3)=(active_endpts(1,:)-pt)*SQRT(active_wts(1))
err(4:6)=(active_endpts(2,:)-pt)*SQRT(active_wts(2))
err(7:9)=(active_endpts(3,:)-pt)*SQRT(active_wts(3))
DEBUG_STACK_POP
end subroutine nurbs_scenter_error
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_surf_error
!---------------------------------------------------------------------------
!> Evalute the error between a surface point and the current active point
!! used in a minimization.
!!
!! @note Designed to be used with MINPACK for non-linear least square location
!! of nearest surface point
!!
!! @param[in] m Number of spatial dimensions (3)
!! @param[in] n Number of parametric dimensions (2)
!! @param[in] uv Parametric possition [n]
!! @param[out] err Error vector between current and desired point [3]
!! @param[in,out] iflag Unused flag
!---------------------------------------------------------------------------
subroutine nurbs_surf_error(m,n,uv,err,iflag)
integer(i4), intent(in) :: m,n
real(r8), intent(in) :: uv(n)
real(r8), intent(out) :: err(m)
integer(i4), intent(inout) :: iflag
real(r8) :: pt(3)
DEBUG_STACK_PUSH
call active_surf%eval(pt,uv(1),uv(2))
err=active_pt-pt
DEBUG_STACK_POP
end subroutine nurbs_surf_error
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_surf_avg
!---------------------------------------------------------------------------
!> Unwrap the periodic coordinates of a NURBS surface to avoid parameter cuts
!! - Example: For a periodic domain of \f$ [0,2 \pi] (3 \pi/4,0) -> (-\pi/4,0) \f$
!!
!! @param[in,out] u1 Parametric coordinate 1 of first point
!! @param[in,out] v1 Parametric coordinate 2 of first point
!! @param[in,out] u2 Parametric coordinate 1 of second point
!! @param[in,out] v2 Parametric coordinate 2 of second point
!---------------------------------------------------------------------------
subroutine nurbs_surf_avg(self,n,uin,vin,uavg,vavg)
class(nurbs_entity), TARGET, intent(in) :: self
integer(i4), intent(in) :: n
real(r8), intent(in) :: uin(n),vin(n)
real(r8), intent(out) :: uavg,vavg
logical :: uskip,vskip
integer(i4) :: nu,nv,i
DEBUG_STACK_PUSH
select type(this=>self)
  type is(nurbs_surf)
    !---
    uavg=0.d0; nu=0
    vavg=0.d0; nv=0
    DO i=1,n
      uskip=.FALSE.
      vskip=.FALSE.
      !---
      if(this%singular(2,1))then
        if(vin(i)==this%domain(2,2))uskip=.TRUE.
      end if
      if(this%singular(1,1))then
        if(vin(i)==this%domain(1,2))uskip=.TRUE.
      end if
      if(this%singular(2,2))then
        if(uin(i)==this%domain(2,1))vskip=.TRUE.
      end if
      if(this%singular(1,2))then
        if(uin(i)==this%domain(1,1))vskip=.TRUE.
      end if
      IF(.NOT.uskip)THEN
        uavg=uavg+uin(i)
        nu=nu+1
      END IF
      IF(.NOT.vskip)THEN
        vavg=vavg+vin(i)
        nv=nv+1
      END IF
    END DO
    uavg=uavg/nu
    vavg=vavg/nv
  class default
    call oft_abort('Attempted to cast non-surface to surface type.','nurbs_surf_avg',__FILE__)
end select
DEBUG_STACK_POP
end subroutine nurbs_surf_avg
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_surf_atsingular
!---------------------------------------------------------------------------
!> Unwrap the periodic coordinates of a NURBS surface to avoid parameter cuts
!! - Example: For a periodic domain of \f$ [0,2 \pi] (3 \pi/4,0) -> (-\pi/4,0) \f$
!!
!! @param[in,out] u1 Parametric coordinate 1 of first point
!! @param[in,out] v1 Parametric coordinate 2 of first point
!! @param[in,out] u2 Parametric coordinate 1 of second point
!! @param[in,out] v2 Parametric coordinate 2 of second point
!---------------------------------------------------------------------------
function nurbs_surf_atsingular(self,u,v) result(test)
class(nurbs_surf), intent(in) :: self
real(r8), intent(inout) :: u,v
logical :: test
real(r8) :: tol=1.d-8
DEBUG_STACK_PUSH
!---
test=.FALSE.
if(self%singular(1,1))then
  if(ABS(v-self%domain(1,2))<tol)test=.TRUE.
end if
if(self%singular(2,1))then
  if(ABS(v-self%domain(2,2))<tol)test=.TRUE.
end if
if(self%singular(1,2))then
  if(ABS(u-self%domain(2,1))<tol)test=.TRUE.
end if
if(self%singular(2,2))then
  if(ABS(u-self%domain(1,1))<tol)test=.TRUE.
end if
DEBUG_STACK_POP
end function nurbs_surf_atsingular
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_surf_unwrap
!---------------------------------------------------------------------------
!> Unwrap the periodic coordinates of a NURBS surface to avoid parameter cuts
!! - Example: For a periodic domain of \f$ [0,2 \pi] (3 \pi/4,0) -> (-\pi/4,0) \f$
!!
!! @param[in,out] u1 Parametric coordinate 1 of first point
!! @param[in,out] v1 Parametric coordinate 2 of first point
!! @param[in,out] u2 Parametric coordinate 1 of second point
!! @param[in,out] v2 Parametric coordinate 2 of second point
!---------------------------------------------------------------------------
subroutine nurbs_surf_unwrap(self,u1,u2,v1,v2)
class(nurbs_surf), intent(in) :: self
real(r8), intent(inout) :: u1,u2,v1,v2
real(r8) :: ubound,lbound,per_buff=.3d0,sep,sep_lim
DEBUG_STACK_PUSH
!---
if(self%singular(1,1))then
  if(v1==self%domain(1,2))u1=u2
  if(v2==self%domain(1,2))u2=u1
end if
if(self%singular(2,1))then
  if(v1==self%domain(2,2))u1=u2
  if(v2==self%domain(2,2))u2=u1
end if
if(self%singular(1,2))then
  if(u1==self%domain(2,1))v1=v2
  if(u2==self%domain(2,1))v2=v1
end if
if(self%singular(2,2))then
  if(u1==self%domain(1,1))v1=v2
  if(u2==self%domain(1,1))v2=v1
end if
!---
if(self%periodic(1))then
  ubound = (1.d0-per_buff)*self%domain(2,1)+per_buff*self%domain(1,1)
  lbound = per_buff*self%domain(2,1)+(1.d0-per_buff)*self%domain(1,1)
  sep=ABS(u1-u2)
  sep_lim=.5d0*ABS(self%domain(2,1)-self%domain(1,1))
  if(sep>sep_lim)u2=u2-self%domain(2,1)+self%domain(1,1)
  !if(u2<lbound.AND.sep>sep_lim)u1=u1-self%domain(2,1)+self%domain(1,1)
  !WRITE(*,*)'unwrap u',sep,sep_lim
end if
if(self%periodic(2))then
  ubound = (1.d0-per_buff)*self%domain(2,2)+per_buff*self%domain(1,2)
  lbound = per_buff*self%domain(2,2)+(1.d0-per_buff)*self%domain(1,2)
  sep=ABS(v1-v2)
  sep_lim=.5d0*ABS(self%domain(2,2)-self%domain(1,2))
  if(sep>sep_lim)v2=v2-self%domain(2,2)+self%domain(1,2)
  !if(v2<lbound.AND.sep>sep_lim)v1=v1-self%domain(2,2)+self%domain(1,2)
  !WRITE(*,*)'unwrap v',sep,sep_lim
end if
DEBUG_STACK_POP
end subroutine nurbs_surf_unwrap
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_surf_wrap
!---------------------------------------------------------------------------
!> Wrap the periodic coordinates of a NURBS surface into a single span
!! - Example: For a periodic domain of \f$ [0,2 \pi] (-\pi/4,0) -> (3 \pi/4,0) \f$
!!
!! @param[in,out] u Parametric coordinate 1
!! @param[in,out] v Parametric coordinate 2
!---------------------------------------------------------------------------
subroutine nurbs_surf_wrap(self,u,v)
class(nurbs_surf), intent(in) :: self
real(r8), intent(inout) :: u,v
real(r8) :: ubound,lbound
DEBUG_STACK_PUSH
!---Limit minimum and maximum values
if(self%periodic(1))then
  ubound = .9d0*self%domain(2,1)+.1d0*self%domain(1,1)
  lbound = .1d0*self%domain(2,1)+.9d0*self%domain(1,1)
  if(u<lbound)u=self%domain(2,1)+self%domain(1,1)-u
end if
if(self%periodic(2))then
  ubound = .9d0*self%domain(2,2)+.1d0*self%domain(1,2)
  lbound = .1d0*self%domain(2,2)+.9d0*self%domain(1,2)
  if(v<lbound)v=self%domain(2,2)+self%domain(1,2)-v
end if
DEBUG_STACK_POP
end subroutine nurbs_surf_wrap
!---------------------------------------------------------------------------
! SUBROUTINE: nurbs_surf_grid
!---------------------------------------------------------------------------
!> Evalute a grid of points evenly spaced in parametric space
!---------------------------------------------------------------------------
subroutine nurbs_surf_grid(self)
class(nurbs_surf), intent(inout) :: self
real(r8) :: u,v
integer(i4) :: i,j
DEBUG_STACK_PUSH
do i=1,nurbs_ngrid
  u=(i-1)*(self%domain(2,1)-self%domain(1,1))*(1.d0/(nurbs_ngrid-1))+self%domain(1,1)
  do j=1,nurbs_ngrid
    v=(j-1)*(self%domain(2,2)-self%domain(1,2))*(1.d0/(nurbs_ngrid-1))+self%domain(1,2)
    call self%eval(self%rgrid(:,i,j),u,v)
  end do
end do
DEBUG_STACK_POP
end subroutine nurbs_surf_grid
#endif
end module nurbs_cad
