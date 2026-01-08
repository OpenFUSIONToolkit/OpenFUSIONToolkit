!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_gs_profiles.F90
!
!> Flux profile definitions
!!
!! @authors Chris Hansen
!! @date March 2014
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
module oft_gs_profiles
USE oft_base
USE spline_mod
USE oft_gs, ONLY: flux_func, gs_eq, oft_indent, oft_increase_indent, &
  oft_decrease_indent
implicit none
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(flux_func) :: zero_flux_func
contains
  !> Needs docs
  procedure :: f => zero_f
  !> Needs docs
  procedure :: fp => zero_fp
  !> Needs docs
  procedure :: update => zero_update
  !> Needs docs
  procedure :: set_cofs => zero_cofs_update
  !> Needs docs
  procedure :: get_cofs => zero_cofs_get
end type zero_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(flux_func) :: flat_flux_func
contains
  !> Needs docs
  procedure :: f => flat_f
  !> Needs docs
  procedure :: fp => flat_fp
  !> Needs docs
  procedure :: update => flat_update
  !> Needs docs
  procedure :: set_cofs => flat_cofs_update
  !> Needs docs
  procedure :: get_cofs => flat_cofs_get
end type flat_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(flux_func) :: linear_flux_func
  real(8) :: alpha = 0.d0 !< Slope parameter
contains
  !> Needs docs
  procedure :: f => linear_f
  !> Needs docs
  procedure :: fp => linear_fp
  !> Needs docs
  procedure :: update => linear_update
  !> Needs docs
  procedure :: set_cofs => linear_cofs_update
  !> Needs docs
  procedure :: get_cofs => linear_cofs_get
end type linear_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(flux_func) :: poly_flux_func
  integer(4) :: deg = 0 !< Degree of polynomial function for F'
  real(8) :: c0 = 1.d0 !< Offset for non-zero edge gradient
  real(8), pointer, dimension(:) :: cofs => NULL() !< Coefficients [deg]
  logical :: zero_grad = .FALSE. !< Force zero gradient at boundary
contains
  !> Needs docs
  procedure :: f => poly_f
  !> Needs docs
  procedure :: fp => poly_fp
  !> Needs docs
  procedure :: update => poly_update
  !> Needs docs
  procedure :: set_cofs => poly_cofs_update
  !> Needs docs
  procedure :: get_cofs => poly_cofs_get
end type poly_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(flux_func) :: spline_flux_func
  INTEGER(4) :: npsi = 0 !< Needs docs
  REAL(8) :: xmin = 0.d0 !< Needs docs
  REAL(8) :: xmax = 0.d0 !< Needs docs
  REAL(8) :: f1 = 0.d0 !< Needs docs
  REAL(8) :: fn = 0.d0 !< Needs docs
  REAL(8) :: yp1 = 0.d0 !< Needs docs
  REAL(8) :: ypn = 0.d0 !< Needs docs
  TYPE(spline_type) :: func !< Needs docs
  TYPE(spline_type), POINTER, DIMENSION(:) :: fun_loc => NULL() !< Needs docs
contains
  !> Needs docs
  procedure :: f => spline_f
  !> Needs docs
  procedure :: fp => spline_fp
  !> Needs docs
  procedure :: update => spline_update
  !> Needs docs
  procedure :: set_cofs => spline_cofs_update
  !> Needs docs
  procedure :: get_cofs => spline_cofs_get
end type spline_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(flux_func) :: linterp_flux_func
  integer(4) :: npsi = 0 !< Needs docs
  real(8) :: y0 = 1.d0 !< Needs docs
  real(8), pointer, dimension(:) :: x => NULL() !< Needs docs
  real(8), pointer, dimension(:) :: yp => NULL() !< Needs docs
  real(8), pointer, dimension(:) :: y => NULL() !< Needs docs
contains
  !> Needs docs
  procedure :: f => linterp_f
  !> Needs docs
  procedure :: fp => linterp_fp
  !> Needs docs
  procedure :: fpp => linterp_fpp
  !> Needs docs
  procedure :: update => linterp_update
  !> Needs docs
  procedure :: set_cofs => linterp_cofs_update
  !> Needs docs
  procedure :: get_cofs => linterp_cofs_get
end type linterp_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(flux_func) :: twolam_flux_func
  real(8) :: alpha = 0.d0 !< Jump height
  real(8) :: sep = .5d0 !< Jump location (poloidal flux)
  real(8) :: width = 200.d0 !< Width factor for TANH
contains
  !> Needs docs
  procedure :: f => twolam_f
  !> Needs docs
  procedure :: fp => twolam_fp
  !> Needs docs
  procedure :: update => twolam_update
  !> Needs docs
  procedure :: set_cofs => twolam_cofs_update
  !> Needs docs
  procedure :: get_cofs => twolam_cofs_get
end type twolam_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(twolam_flux_func) :: stepslant_flux_func
  real(8) :: beta = 0.d0
contains
  !> Needs docs
  procedure :: f => stepslant_f
  !> Needs docs
  procedure :: fp => stepslant_fp
  !> Needs docs
  procedure :: update => stepslant_update
  !> Needs docs
  procedure :: set_cofs => stepslant_cofs_update
  !> Needs docs
  procedure :: get_cofs => stepslant_cofs_get
end type stepslant_flux_func
!------------------------------------------------------------------------------
! CLASS wesson_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(flux_func) :: wesson_flux_func
  real(8) :: gamma = 0.d0
contains
  procedure :: f => wesson_f
  procedure :: fp => wesson_fp
  procedure :: update => wesson_update
  procedure :: set_cofs => wesson_cofs_update
  procedure :: get_cofs => wesson_cofs_get
end type wesson_flux_func
contains
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function zero_f(self,psi) result(b)
class(zero_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b
b=0.d0
end function zero_f
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function zero_fp(self,psi) result(b)
class(zero_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b
b=0.d0
end function zero_fp
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine zero_update(self,gseq)
class(zero_flux_func), intent(inout) :: self
class(gs_eq), intent(inout) :: gseq
end subroutine zero_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function zero_cofs_update(self,c) result(ierr)
class(zero_flux_func), intent(inout) :: self
real(8), intent(in) :: c(:)
integer(4) :: ierr
ierr=0
end function zero_cofs_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine zero_cofs_get(self,c)
class(zero_flux_func), intent(inout) :: self
real(8), intent(out) :: c(:)
end subroutine zero_cofs_get
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_flat_f(func)
CLASS(flux_func), POINTER, INTENT(out) :: func
ALLOCATE(flat_flux_func::func)
select type(self=>func)
  type is(flat_flux_func)
    self%ncofs=0
end select
END SUBROUTINE create_flat_f
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function flat_f(self,psi) result(b)
class(flat_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b,x1,x2
IF(self%plasma_bounds(1)<-1.d98)THEN
  b=psi
  RETURN
END IF
x1=self%plasma_bounds(1)
x2=self%plasma_bounds(2)
IF(psi>x1)THEN
  b=(psi-x1)
ELSE
  b=0.d0
END IF
end function flat_f
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function flat_fp(self,psi) result(b)
class(flat_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b,x1,x2
IF(self%plasma_bounds(1)<-1.d98)THEN
  b=1.d0
  RETURN
END IF
x1=self%plasma_bounds(1)
x2=self%plasma_bounds(2)
IF(psi>x1)THEN
  b=1.d0
ELSE
  b=0.d0
END IF
end function flat_fp
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine flat_update(self,gseq)
class(flat_flux_func), intent(inout) :: self
class(gs_eq), intent(inout) :: gseq
self%plasma_bounds=gseq%plasma_bounds
end subroutine flat_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function flat_cofs_update(self,c) result(ierr)
class(flat_flux_func), intent(inout) :: self
real(8), intent(in) :: c(:)
integer(4) :: ierr
ierr=0
end function flat_cofs_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine flat_cofs_get(self,c)
class(flat_flux_func), intent(inout) :: self
real(8), intent(out) :: c(:)
end subroutine flat_cofs_get
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_linear_ff(func,alpha)
CLASS(flux_func), POINTER, INTENT(out) :: func
REAL(8), INTENT(in) :: alpha
ALLOCATE(linear_flux_func::func)
select type(self=>func)
  type is(linear_flux_func)
    self%ncofs=1
    self%alpha=alpha
end select
END SUBROUTINE create_linear_ff
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function linear_f(self,psi) result(b)
class(linear_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b,x1,x2
IF(self%plasma_bounds(1)<-1.d98)THEN
  b=self%alpha*(psi**2 - psi) + psi
  RETURN
END IF
x1=self%plasma_bounds(1)
x2=self%plasma_bounds(2)
IF(psi>x1)THEN
  b = psi*(self%alpha*(psi - 2.d0*x1 - (x2-x1)) + (x2-x1))/(x2-x1) &
  - x1*(self%alpha*(x1 - 2.d0*x1 - (x2-x1)) + (x2-x1))/(x2-x1)
ELSE
  b = 0.d0
END IF
end function linear_f
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function linear_fp(self,psi) result(b)
class(linear_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b,x1,x2
IF(self%plasma_bounds(1)<-1.d98)THEN
  b=1.d0 + self%alpha*(2.d0*psi - 1.d0)
  RETURN
END IF
x1=self%plasma_bounds(1)
x2=self%plasma_bounds(2)
IF(psi>x1)THEN
  b=1.d0 + self%alpha*(2.d0*(psi-x1)/(x2-x1) - 1.d0)
ELSE
  b=0.d0
END IF
end function linear_fp
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine linear_update(self,gseq)
class(linear_flux_func), intent(inout) :: self
class(gs_eq), intent(inout) :: gseq
self%plasma_bounds=gseq%plasma_bounds
end subroutine linear_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function linear_cofs_update(self,c) result(ierr)
class(linear_flux_func), intent(inout) :: self
real(8), intent(in) :: c(:)
integer(4) :: ierr
ierr=0
self%alpha=c(1)
end function linear_cofs_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine linear_cofs_get(self,c)
class(linear_flux_func), intent(inout) :: self
real(8), intent(out) :: c(:)
c(1)=self%alpha
end subroutine linear_cofs_get
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_poly_ff(func,ncofs,cofs,zero_fp)
CLASS(flux_func), POINTER, INTENT(out) :: func
INTEGER(4), INTENT(in) :: ncofs
REAL(8), OPTIONAL, INTENT(in) :: cofs(:)
LOGICAL, OPTIONAL, INTENT(in) :: zero_fp
INTEGER(4) :: i,offset

ALLOCATE(poly_flux_func::func)
select type(self=>func)
  type is(poly_flux_func)
  !---
  self%deg=ncofs
  self%ncofs=ncofs
  IF(PRESENT(zero_fp))self%zero_grad=zero_fp
  offset=0
  IF(self%zero_grad)THEN
    offset=1
    self%ncofs=ncofs-1
    self%deg=ncofs-1
  END IF
  ALLOCATE(self%cofs(self%ncofs))
  self%cofs=0.d0
  !---
  IF(PRESENT(cofs))THEN
    DO i=1,self%ncofs
      self%cofs(i)=cofs(offset+i)
    END DO
  END IF
  self%c0=1.d0
  IF(.NOT.self%zero_grad)THEN
    self%c0=0.d0
    DO i=1,self%ncofs
      self%c0=self%c0+self%cofs(i)/REAL((i+1),8)
    END DO
    self%c0=1.d0-self%c0
  END IF
end select

END SUBROUTINE create_poly_ff
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function poly_f(self,psi) result(b)
class(poly_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
integer(4) :: i,j,offset
real(8) :: b,x1,x2
offset=0
IF(self%plasma_bounds(1)<-1.d98)THEN
  IF(self%zero_grad)THEN
    b=(psi**2)/2.d0
  ELSE
    b=self%c0*psi
  END IF
  DO i=1,self%deg
    j=offset+i
    b=b+self%cofs(i)*(psi**(j+1))/REAL(j+1,8)
  END DO
  RETURN
END IF
x1=self%plasma_bounds(1)
x2=self%plasma_bounds(2)
IF(psi>x1)THEN
  IF(self%zero_grad)THEN
    b = (psi*(psi-2.d0*x1) + x1**2)/REAL(2.d0*(x2-x1),8)
    offset=1
  ELSE
    b=self%c0*(psi-x1)
  END IF
  DO i=1,self%deg
    j=offset+i
    b=b+self%cofs(i)*((psi-x1)**(j+1))/REAL((j+1)*(x2-x1)**j,8)
  END DO
ELSE
  b = 0.d0
END IF
end function poly_f
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function poly_fp(self,psi) result(b)
class(poly_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
integer(4) :: i,j,offset
real(8) :: b,x1,x2
offset=0
IF(self%zero_grad)THEN
  b=psi
  offset=1
ELSE
  b=self%c0
END IF
IF(self%plasma_bounds(1)<-1.d98)THEN
  DO i=1,self%deg
    j=offset+i
    b=b+self%cofs(i)*(psi**j)
  END DO
  RETURN
END IF
x1=self%plasma_bounds(1)
x2=self%plasma_bounds(2)
IF(psi>x1)THEN
  IF(self%zero_grad)THEN
    b = (psi-x1)/(x2-x1)
    offset=1
  END IF
  DO i=1,self%deg
    j=offset+i
    b=b+self%cofs(i)*((psi-x1)**j)/(x2-x1)**j
  END DO
ELSE
  b = 0.d0
END IF
end function poly_fp
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine poly_update(self,gseq)
class(poly_flux_func), intent(inout) :: self
class(gs_eq), intent(inout) :: gseq
self%plasma_bounds=gseq%plasma_bounds
end subroutine poly_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function poly_cofs_update(self,c) result(ierr)
class(poly_flux_func), intent(inout) :: self
real(8), intent(in) :: c(:)
INTEGER(4) :: i,ierr
DO i=1,self%ncofs
  self%cofs(i)=c(i)
END DO
IF(.NOT.self%zero_grad)THEN
  self%c0=0.d0
  DO i=1,self%ncofs
    self%c0=self%c0+self%cofs(i)/REAL((i+1),8)
  END DO
  self%c0=1.d0-self%c0
END IF
ierr=0
end function poly_cofs_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine poly_cofs_get(self,c)
class(poly_flux_func), intent(inout) :: self
real(8), intent(out) :: c(:)
INTEGER(4) :: i
DO i=1,self%ncofs
  c(i)=self%cofs(i)
END DO
end subroutine poly_cofs_get
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_spline_ff(func,npsi,psimin,psimax,psivals)
CLASS(flux_func), POINTER, INTENT(out) :: func
INTEGER(4), INTENT(in) :: npsi
REAL(8), OPTIONAL, INTENT(in) :: psimin
REAL(8), OPTIONAL, INTENT(in) :: psimax
REAL(8), OPTIONAL, INTENT(in) :: psivals(npsi)
REAL(8) :: psi1,psi2
REAL(8), ALLOCATABLE :: c(:)
INTEGER(4) :: i,ierr
ALLOCATE(spline_flux_func::func)
select type(self=>func)
  type is(spline_flux_func)
  !---
  self%npsi=npsi
  self%ncofs=self%npsi
  ALLOCATE(self%fun_loc(omp_get_max_threads()))
  CALL spline_alloc(self%func,self%npsi-1,1)
  DO i=1,omp_get_max_threads()
    CALL spline_alloc(self%fun_loc(i),self%npsi-1,1)
  END DO
  IF(PRESENT(psivals))THEN
    !---
    DO i=1,self%npsi
      self%func%xs(i-1)=psivals(i)
      self%func%fs(i-1,1)=1.d0*psivals(i)
    END DO
  ELSE
    !---
    psi1=0.d0
    psi2=1.d0
    IF(PRESENT(psimin))psi1=psimin
    IF(PRESENT(psimax))psi2=psimax
    !---
    DO i=1,self%npsi
      self%func%xs(i-1)=(psi2-psi1)*i/(self%npsi+1) + psi1
      self%func%fs(i-1,1)=1.d0*self%func%xs(i-1)
    END DO
  END IF
  self%xmin=self%func%xs(0)
  self%xmax=self%func%xs(self%npsi-1)
  !---
  WRITE(*,*)'Fitting',self%func%xs
  !CALL spline_fit(self%func,"extrap")
  CALL spline_fit(self%func,"not-a-knot")
  DO i=1,omp_get_max_threads()
    CALL spline_copy(self%func,self%fun_loc(i))
  END DO
  ALLOCATE(c(self%ncofs))
  CALL self%get_cofs(c)
  ierr=self%set_cofs(c)
  DEALLOCATE(c)
  WRITE(*,*)'Spline Created',self%ncofs,self%func%xs
end select

END SUBROUTINE create_spline_ff
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function spline_f(self,psi) result(b)
class(spline_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b
IF(psi<self%xmin)THEN
  b=self%yp1*psi
ELSE IF(psi>=self%xmin.AND.psi<self%xmax)THEN
  CALL spline_eval(self%fun_loc(omp_get_thread_num()+1),psi,0)
  b=self%fun_loc(omp_get_thread_num()+1)%f(1)-self%f1!self%yp1*psi
ELSE
  b=self%ypn*(psi-self%xmax)+self%fn
END IF
end function spline_f
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function spline_fp(self,psi) result(b)
class(spline_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b
IF(psi<self%xmin)THEN
  b=self%yp1
ELSE IF(psi>=self%xmin.AND.psi<self%xmax)THEN
  CALL spline_eval(self%fun_loc(omp_get_thread_num()+1),psi,1)
  b=self%fun_loc(omp_get_thread_num()+1)%f1(1)
ELSE
  b=self%ypn
END IF
end function spline_fp
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine spline_update(self,gseq)
class(spline_flux_func), intent(inout) :: self
class(gs_eq), intent(inout) :: gseq
REAL(8) :: yp1,f0
INTEGER(4) :: i
self%plasma_bounds=gseq%plasma_bounds
RETURN
CALL spline_eval(self%func,self%xmin,1)
f0=self%func%f(1)
yp1=self%func%f1(1)
!---
self%func%fs(:,1)=self%func%fs(:,1)/yp1
!CALL spline_fit(self%func,"extrap")
CALL spline_fit(self%func,"not-a-knot")
!---
CALL spline_eval(self%func,self%xmin,1)
self%yp1=self%func%f1(1)
IF(ABS(self%yp1-1.d0)>1.d-4)CALL oft_abort('Error in spline normalization','spline_update',__FILE__)
CALL spline_eval(self%func,self%xmax,1)
self%ypn=self%func%f1(1)
!---
DO i=1,omp_get_max_threads()
  CALL spline_copy(self%func,self%fun_loc(i))
END DO
end subroutine spline_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function spline_cofs_update(self,c) result(ierr)
class(spline_flux_func), intent(inout) :: self
real(8), intent(in) :: c(:)
REAL(8) :: yp1
integer(4) :: ierr
INTEGER(4) :: i
DO i=0,self%npsi-1
  self%func%fs(i,1)=c(i+1)
END DO
!CALL spline_fit(self%func,"extrap")
CALL spline_fit(self%func,"not-a-knot")
!---
CALL spline_eval(self%func,self%xmin,1)
yp1=self%func%f1(1)
!---
self%func%fs(:,1)=self%func%fs(:,1)/yp1
!CALL spline_fit(self%func,"extrap")
CALL spline_fit(self%func,"not-a-knot")
!---
CALL spline_eval(self%func,self%xmin,1)
self%f1=self%func%f(1)-self%xmin
self%yp1=self%func%f1(1)
!IF(ABS(self%yp1-1.d0)>1.d-4)CALL oft_abort('Error in spline normalization','spline_cofs_update',__FILE__)
CALL spline_eval(self%func,self%xmax,1)
self%fn=self%func%f(1)-self%f1
self%ypn=self%func%f1(1)
!---
DO i=1,omp_get_max_threads()
  CALL spline_copy(self%func,self%fun_loc(i))
END DO
ierr=0
end function spline_cofs_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine spline_cofs_get(self,c)
class(spline_flux_func), intent(inout) :: self
real(8), intent(out) :: c(:)
INTEGER(4) :: i
DO i=0,self%npsi-1
  c(i+1)=self%func%fs(i,1)
END DO
end subroutine spline_cofs_get
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_linterp_ff(func,npsi,psivals,yvals,y0)
CLASS(flux_func), POINTER, INTENT(out) :: func
INTEGER(4), INTENT(in) :: npsi
REAL(8), INTENT(in) :: psivals(npsi)
REAL(8), INTENT(in) :: yvals(npsi)
REAL(8), INTENT(in) :: y0
INTEGER(4) :: i,ierr
ALLOCATE(linterp_flux_func::func)
SELECT TYPE(self=>func)
  TYPE IS(linterp_flux_func)
  !---
  self%npsi=npsi
  self%ncofs=self%npsi
  !---
  ALLOCATE(self%x(self%npsi))
  ALLOCATE(self%yp(self%npsi))
  ALLOCATE(self%y(self%npsi))
  !---
  self%y0=y0
  DO i=1,self%npsi
    self%x(i)=psivals(i)
    self%yp(i)=yvals(i)
  END DO
  IF(self%y0<1.d-8)THEN
    self%ncofs=self%ncofs-1
    ierr=self%set_cofs(yvals(2:self%npsi))
  ELSE
    ierr=self%set_cofs(yvals)
  END IF
  IF(oft_debug_print(1))WRITE(*,*)'Linear interpolator Created',self%ncofs,self%x,self%y0
END SELECT

END SUBROUTINE create_linterp_ff
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function linterp_f(self,psi) result(b)
class(linterp_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b,psihat
integer(4) :: i
IF(self%plasma_bounds(1)<-1.d98)THEN
  psihat=psi
ELSE
  psihat=(psi-self%plasma_bounds(1))/(self%plasma_bounds(2)-self%plasma_bounds(1))
END IF
b=0.d0
IF(psihat<0.d0)RETURN
!
if(psihat<=self%x(1))then
  b = psihat*(.5d0*psihat)/self%x(1)*(self%yp(1)-self%y0) + psihat*self%y0
else if(psihat>=self%x(self%npsi))then
  b=self%y(self%npsi)+self%yp(self%npsi)*(psihat-self%x(self%npsi))
else
  do i=2,self%npsi
    if(psihat>self%x(i-1).AND.psihat<=self%x(i))then
      b = self%y(i-1) + psihat*(.5d0*psihat - self%x(i-1))/(self%x(i) - self%x(i-1))* &
          (self%yp(i) - self%yp(i-1)) + psihat*self%yp(i-1)
    end if
  end do
end if
IF(self%plasma_bounds(1)>-1.d97)THEN
  b=b*(self%plasma_bounds(2)-self%plasma_bounds(1))
END IF
end function linterp_f
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function linterp_fp(self,psi) result(b)
class(linterp_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b,x,psihat
integer(4) :: i
IF(self%plasma_bounds(1)<-1.d98)THEN
  psihat=psi
ELSE
  psihat=(psi-self%plasma_bounds(1))/(self%plasma_bounds(2)-self%plasma_bounds(1))
END IF
b=0.d0
IF(psihat<0.d0)RETURN
!
if(psihat<=self%x(1))then
  x = psihat/self%x(1)
  b = self%y0 + (self%yp(1)-self%y0)*x
else if(psihat>=self%x(self%npsi))then
  b=self%yp(self%npsi)
else
  do i=2,self%npsi
    if(psihat>self%x(i-1).AND.psihat<=self%x(i))then
      x=(psihat-self%x(i-1))/(self%x(i)-self%x(i-1))
      b = self%yp(i-1) + (self%yp(i)-self%yp(i-1))*x
    end if
  end do
end if
end function linterp_fp
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function linterp_fpp(self,psi) result(b)
class(linterp_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b,x,psihat
integer(4) :: i
IF(self%plasma_bounds(1)<-1.d98)THEN
  psihat=psi
ELSE
  psihat=(psi-self%plasma_bounds(1))/(self%plasma_bounds(2)-self%plasma_bounds(1))
END IF
b=0.d0
IF(psihat<0.d0)RETURN
!
if(psihat<=self%x(1))then
  x = 1.d0/(self%plasma_bounds(2)-self%plasma_bounds(1))/self%x(1)
  b = (self%yp(1)-self%y0)*x
else if(psihat>=self%x(self%npsi))then
  b  = 0.d0
else
  do i=2,self%npsi
    if(psihat>self%x(i-1).AND.psihat<=self%x(i))then
      x=1.d0/(self%plasma_bounds(2)-self%plasma_bounds(1))/(self%x(i)-self%x(i-1))
      b = (self%yp(i)-self%yp(i-1))*x
    end if
  end do
end if
end function linterp_fpp
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine linterp_update(self,gseq)
class(linterp_flux_func), intent(inout) :: self
class(gs_eq), intent(inout) :: gseq
self%plasma_bounds=gseq%plasma_bounds
end subroutine linterp_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function linterp_cofs_update(self,c) result(ierr)
class(linterp_flux_func), intent(inout) :: self
real(8), intent(in) :: c(:)
real(8) :: x,xp
integer(4) :: ierr
integer(4) :: i,offset
!---
offset=0
IF(self%ncofs==self%npsi-1)offset=1
DO i=1,self%ncofs
  self%yp(i+offset)=c(i)
END DO
!---
DO i=1,self%npsi
  x=self%x(i)
  !---
  IF(i==1)THEN
    self%y(i)=x*(.5d0*x)/self%x(1)*(self%yp(1)-self%y0) + x*self%y0
  ELSE
    xp=self%x(i-1)
    self%y(i-1)=self%y(i-1) - (xp*(.5d0*xp-xp)/(x-xp)*(self%yp(i)-self%yp(i-1)) &
               + xp*self%yp(i-1))
    self%y(i)=x*(.5d0*x - xp)/(x - xp)*(self%yp(i) - self%yp(i-1)) + x*self%yp(i-1)
    self%y(i)=self%y(i) + self%y(i-1)
  END IF
END DO
IF(oft_debug_print(2))THEN
  WRITE(*,'(2A)')oft_indent,'Update Linterp:'
  CALL oft_increase_indent
  WRITE(*,'(2A,100ES11.3)')oft_indent,' x  =',self%x
  WRITE(*,'(2A,100ES11.3)')oft_indent,' yp =',self%yp
  WRITE(*,'(2A,100ES11.3)')oft_indent,' y  =',self%y
  WRITE(*,'(2A,ES11.3)')oft_indent,' y0  =',self%y0
  CALL oft_decrease_indent
END IF
ierr=0
end function linterp_cofs_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine linterp_cofs_get(self,c)
class(linterp_flux_func), intent(inout) :: self
real(8), intent(out) :: c(:)
integer(4) :: i,offset
!---
offset=0
IF(self%ncofs==self%npsi-1)offset=1
DO i=1,self%ncofs
  c(i)=self%yp(i+offset)
END DO
end subroutine linterp_cofs_get
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_twolam_ff(func,sep,alpha)
CLASS(flux_func), POINTER, INTENT(out) :: func
REAL(8), INTENT(in) :: sep
REAL(8), INTENT(in) :: alpha

ALLOCATE(twolam_flux_func::func)
select type(self=>func)
  type is(twolam_flux_func)
  !---
  self%ncofs=2
  self%sep=sep
  self%alpha=alpha
end select

END SUBROUTINE create_twolam_ff
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function twolam_f(self,psi) result(b)
class(twolam_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b
b=(1.d0+self%alpha/2.d0)*psi + self%alpha*LOG(COSH(self%width*(psi-self%sep)))/2.d0/self%width &
- self%alpha*LOG(COSH(-self%width*self%sep))/2.d0/self%width
end function twolam_f
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function twolam_fp(self,psi) result(b)
class(twolam_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b
b=1.d0 + self%alpha*(1.d0 + TANH(self%width*(psi-self%sep)))/2.d0
end function twolam_fp
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine twolam_update(self,gseq)
class(twolam_flux_func), intent(inout) :: self
class(gs_eq), intent(inout) :: gseq
self%plasma_bounds=gseq%plasma_bounds
end subroutine twolam_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function twolam_cofs_update(self,c) result(ierr)
class(twolam_flux_func), intent(inout) :: self
real(8), intent(in) :: c(:)
integer(4) :: ierr
ierr=0
IF(c(1)<0.d0.OR.c(1)>1.d0)ierr=-1
IF(ABS(c(2))>1.d0)ierr=-1
self%sep=c(1)
self%alpha=c(2)
end function twolam_cofs_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine twolam_cofs_get(self,c)
class(twolam_flux_func), intent(inout) :: self
real(8), intent(out) :: c(:)
c(1)=self%sep
c(2)=self%alpha
end subroutine twolam_cofs_get
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_stepslant_ff(func,sep,alpha,beta)
CLASS(flux_func), POINTER, INTENT(out) :: func
REAL(8), INTENT(in) :: sep
REAL(8), INTENT(in) :: alpha
REAL(8), INTENT(in) :: beta

ALLOCATE(stepslant_flux_func::func)
select type(self=>func)
  type is(stepslant_flux_func)
  !---
  self%ncofs=2
  self%sep=sep
  self%alpha=alpha
  self%beta=beta
end select

END SUBROUTINE create_stepslant_ff
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function stepslant_f(self,psi) result(b)
class(stepslant_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b
if(psi<self%sep)then
  b=psi
else
  b=self%sep + (1.d0+self%alpha)*(self%sep-psi)*((self%beta-2)*self%sep-self%beta*psi)/(2*self%sep)
end if
end function stepslant_f
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function stepslant_fp(self,psi) result(b)
class(stepslant_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b
if(psi<self%sep) then
  b=1.d0
else
  b=(1.d0+self%alpha)*(1.d0+self%beta*(psi-self%sep)/self%sep)
end if
end function stepslant_fp
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine stepslant_update(self,gseq)
class(stepslant_flux_func), intent(inout) :: self
class(gs_eq), intent(inout) :: gseq
self%plasma_bounds=gseq%plasma_bounds
end subroutine stepslant_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function stepslant_cofs_update(self,c) result(ierr)
class(stepslant_flux_func), intent(inout) :: self
real(8), intent(in) :: c(:)
integer(4) :: ierr
self%sep=c(1)
self%alpha=c(2)
self%beta=c(3)
ierr=0
end function stepslant_cofs_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine stepslant_cofs_get(self,c)
class(stepslant_flux_func), intent(inout) :: self
real(8), intent(out) :: c(:)
c(1)=self%sep
c(2)=self%alpha
c(3)=self%beta
end subroutine stepslant_cofs_get
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_wesson_ff(func,ncofs,gamma)
CLASS(flux_func), POINTER, INTENT(out) :: func
INTEGER(4), INTENT(in) :: ncofs
REAL(8), INTENT(in) :: gamma

ALLOCATE(wesson_flux_func::func)
select type(self=>func)
  type is(wesson_flux_func)
  !---
  self%ncofs=ncofs
  self%gamma=gamma
end select

END SUBROUTINE create_wesson_ff
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function wesson_f(self,psi) result(b)
class(wesson_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b,x1,x2,psihat
IF(self%plasma_bounds(1)<-1.d98)THEN
  b = (psi**self%gamma)/self%gamma
  RETURN
END IF
x1=self%plasma_bounds(1)
x2=self%plasma_bounds(2)
psihat=(psi-x1)/(x2-x1)
IF(psi>x1)THEN
  b = (psihat**self%gamma)*(x2-x1)/self%gamma
ELSE
  b = 0.d0
END IF
end function wesson_f
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function wesson_fp(self,psi) result(b)
class(wesson_flux_func), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b,x1,x2,psihat
IF(self%plasma_bounds(1)<-1.d98)THEN
  b = psi**(self%gamma-1.d0)
  RETURN
END IF
x1=self%plasma_bounds(1)
x2=self%plasma_bounds(2)
psihat=(psi-x1)/(x2-x1)
IF(psi>x1)THEN
  b = psihat**(self%gamma-1.d0)
ELSE
  b = 0.d0
END IF
end function wesson_fp
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine wesson_update(self,gseq)
class(wesson_flux_func), intent(inout) :: self
class(gs_eq), intent(inout) :: gseq
self%plasma_bounds=gseq%plasma_bounds
end subroutine wesson_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function wesson_cofs_update(self,c) result(ierr)
class(wesson_flux_func), intent(inout) :: self
real(8), intent(in) :: c(:)
integer(4) :: ierr
self%gamma=c(1)
ierr=0
end function wesson_cofs_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine wesson_cofs_get(self,c)
class(wesson_flux_func), intent(inout) :: self
real(8), intent(out) :: c(:)
c(1)=self%gamma
end subroutine wesson_cofs_get
end module oft_gs_profiles
