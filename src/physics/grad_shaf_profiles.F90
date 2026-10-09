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
USE oft_io, ONLY: hdf5_read, hdf5_write, hdf5_field_exist
USE oft_base
USE spline_mod
USE oft_gs, ONLY: flux_func, gs_equil, oft_indent, oft_increase_indent, &
  oft_decrease_indent, gs_get_qprof
USE oft_gs_cutcell, ONLY: gs_qgeom_cutcell
implicit none
INTEGER(4) :: torflux_qgeom_backend = 0 !< Toroidal flux map q/F backend: 0 cut-cell, 1 tracer (testing only)
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(flux_func) :: zero_flux_func
contains
  !> Delete profile
  procedure :: delete => zero_delete
  !> Needs docs
  procedure :: copy => zero_copy
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
  !> Needs docs
  procedure :: save_hdf5 => zero_save_hdf5
  procedure :: save_txt => zero_save_txt
  !> Needs docs
  procedure :: load_hdf5 => zero_load_hdf5
  procedure :: load_txt => zero_load_txt
end type zero_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(flux_func) :: flat_flux_func
contains
  !> Delete profile
  procedure :: delete => flat_delete
  !> Needs docs
  procedure :: copy => flat_copy
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
  !> Needs docs
  procedure :: save_hdf5 => flat_save_hdf5
  procedure :: save_txt => flat_save_txt
  !> Needs docs
  procedure :: load_hdf5 => flat_load_hdf5
  procedure :: load_txt => flat_load_txt
end type flat_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(flux_func) :: poly_flux_func
  integer(4) :: deg = 0 !< Degree of polynomial function for F'
  real(8) :: c0 = 1.d0 !< Offset for non-zero edge gradient
  real(8), pointer, dimension(:) :: cofs => NULL() !< Coefficients [deg]
  logical :: zero_grad = .FALSE. !< Force zero gradient at boundary
contains
  !> Delete profile
  procedure :: delete => poly_delete
  !> Needs docs
  procedure :: copy => poly_copy
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
  !> Needs docs
  procedure :: save_hdf5 => poly_save_hdf5
  procedure :: save_txt => poly_save_txt
  !> Needs docs
  procedure :: load_hdf5 => poly_load_hdf5
  procedure :: load_txt => poly_load_txt
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
  !> Delete profile
  procedure :: delete => spline_func_delete
  !> Needs docs
  procedure :: copy => spline_func_copy
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
  !> Needs docs
  procedure :: save_hdf5 => spline_save_hdf5
  procedure :: save_txt => spline_save_txt
  !> Needs docs
  procedure :: load_hdf5 => spline_load_hdf5
  procedure :: load_txt => spline_load_txt
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
  class(flux_func), pointer :: phi_wrapper => NULL() !< Enclosing @ref phi_flux_wrapper if `x` is on \f$ \hat{\Phi} \f$
contains
  !> Delete profile
  procedure :: delete => linterp_delete
  !> Needs docs
  procedure :: copy => linterp_copy
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
  !> Needs docs
  procedure :: save_hdf5 => linterp_save_hdf5
  procedure :: save_txt => linterp_save_txt
  !> Needs docs
  procedure :: load_hdf5 => linterp_load_hdf5
  procedure :: load_txt => linterp_load_txt
end type linterp_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(linterp_flux_func) :: mlinterp_flux_func
  integer(4) :: nbasis = 0 !< Needs docs
  real(8), pointer, dimension(:,:) :: yp_basis => NULL() !< Needs docs
  real(8), pointer, dimension(:) :: weights => NULL() !< Needs docs
contains
  !> Delete profile
  procedure :: delete => mlinterp_delete
  !> Needs docs
  procedure :: copy => mlinterp_copy
  !> Needs docs
  procedure :: update => mlinterp_update
  !> Needs docs
  procedure :: set_cofs => mlinterp_cofs_update
  !> Needs docs
  procedure :: get_cofs => mlinterp_cofs_get
  !> Needs docs
  procedure :: save_hdf5 => mlinterp_save_hdf5
  procedure :: save_txt => mlinterp_save_txt
  !> Needs docs
  procedure :: load_hdf5 => mlinterp_load_hdf5
  procedure :: load_txt => mlinterp_load_txt
end type mlinterp_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(flux_func) :: wesson_flux_func
  real(8) :: gamma = 0.d0
contains
  procedure :: delete => wesson_delete
  procedure :: copy => wesson_copy
  procedure :: f => wesson_f
  procedure :: fp => wesson_fp
  procedure :: update => wesson_update
  procedure :: set_cofs => wesson_cofs_update
  procedure :: get_cofs => wesson_cofs_get
  procedure :: save_hdf5 => wesson_save_hdf5
  procedure :: save_txt => wesson_save_txt
  procedure :: load_hdf5 => wesson_load_hdf5
  procedure :: load_txt => wesson_load_txt
end type wesson_flux_func
!------------------------------------------------------------------------------
!> Map from normalized poloidal to normalized toroidal flux, \f$ \hat{\Phi}(\hat{\psi}) \f$
!!
!! Stored in the internal convention (\f$ 1-\hat{\psi} \f$, \f$ 1-\hat{\Phi} \f$: 0 at the LCFS,
!! 1 at the axis) as a cubic Hermite interpolant through flux surfaces with slopes
!! \f$ J = d\hat{\Phi}/d\hat{\psi} = q/Q \f$. One map (`gs_equil%torflux`) is shared by all
!! @ref phi_flux_wrapper profiles of an equilibrium and rebuilt from q by `update`.
!! As a flux function, `f` is \f$ 1-\hat{\Phi} \f$ and `fp` its derivative.
!------------------------------------------------------------------------------
type, extends(flux_func) :: torflux_map
  logical :: active = .FALSE. !< Any profile of the equilibrium on toroidal flux?
  logical :: settle = .TRUE. !< Settle surface placement fully on the next update
  integer(4) :: ns = 0 !< Number of surfaces (identity map if < 2)
  integer(4) :: stamp = 0 !< Incremented on each map update
  real(8) :: Q = 1.d0 !< \f$ \int_0^1 q d\hat{\psi} \f$
  real(8), pointer, dimension(:) :: psis => NULL() !< \f$ 1-\hat{\psi} \f$ of each surface (ascending)
  real(8), pointer, dimension(:) :: phihat => NULL() !< \f$ 1-\hat{\Phi} \f$ of each surface
  real(8), pointer, dimension(:) :: jac => NULL() !< \f$ d\hat{\Phi}/d\hat{\psi} = q/Q \f$ of each surface
  real(8), pointer, dimension(:) :: nodes => NULL() !< Profile nodes the surfaces were placed for
  integer(4), pointer, dimension(:) :: bin => NULL() !< First interval of each uniform \f$ 1-\hat{\psi} \f$ bin
contains
  !> Evaluate \f$ 1-\hat{\Phi} \f$ and derivatives at \f$ 1-\hat{\psi} \f$
  procedure :: eval => tmap_eval
  !> Evaluate \f$ 1-\hat{\psi} \f$ at \f$ 1-\hat{\Phi} \f$
  procedure :: inv => tmap_inv
  procedure :: delete => tmap_delete
  procedure :: copy => tmap_copy
  procedure :: f => tmap_f
  procedure :: fp => tmap_fp
  procedure :: update => tmap_update
  procedure :: set_cofs => tmap_cofs_update
  procedure :: get_cofs => tmap_cofs_get
  procedure :: save_hdf5 => tmap_save_hdf5
  procedure :: save_txt => tmap_save_txt
  procedure :: load_hdf5 => tmap_load_hdf5
  procedure :: load_txt => tmap_load_txt
end type torflux_map
!------------------------------------------------------------------------------
!> Flux function on normalized toroidal flux, evaluated through @ref torflux_map
!!
!! `phi_rep` is any profile whose coordinate is \f$ \hat{\Phi} \f$ instead of \f$ \hat{\psi} \f$:
!! - "phi_n" (`relabel=.FALSE.`): \f$ f(\psi) = Y(\hat{\Phi}) \f$ and \f$ f' = Y'(\hat{\Phi}) J \f$
!! - "phi_n_relabel" (`relabel=.TRUE.`): `phi_rep` holds \f$ dY/d\hat{\psi} \f$ on a \f$ \hat{\Phi} \f$
!!   grid, so \f$ f' = Y'(\hat{\Phi}) \f$ and `f` integrates it over \f$ \hat{\psi} \f$
!------------------------------------------------------------------------------
type, extends(flux_func) :: phi_flux_wrapper
  logical :: relabel = .FALSE. !< "phi_n_relabel" (otherwise "phi_n")
  class(flux_func), pointer :: phi_rep => NULL() !< Profile on \f$ \hat{\Phi} \f$
  class(torflux_map), pointer :: map => NULL() !< Equilibrium map (`gs_equil%torflux`)
  integer(4) :: fint_stamp = -1 !< Map stamp `fint` was built for
  real(8), pointer, dimension(:) :: fint => NULL() !< Cumulative \f$ \int f' d\hat{\psi} \f$ on map surfaces (relabel)
  real(8), pointer, dimension(:,:) :: fcof => NULL() !< Per-interval quartic for partial integrals of `fint` [4,ns-1]
contains
  procedure :: delete => phi_wrap_delete
  procedure :: copy => phi_wrap_copy
  procedure :: f => phi_wrap_f
  procedure :: fp => phi_wrap_fp
  procedure :: fpp => phi_wrap_fpp
  procedure :: update => phi_wrap_update
  procedure :: set_cofs => phi_wrap_cofs_update
  procedure :: get_cofs => phi_wrap_cofs_get
  procedure :: save_hdf5 => phi_wrap_save_hdf5
  procedure :: save_txt => phi_wrap_save_txt
  procedure :: load_hdf5 => phi_wrap_load_hdf5
  procedure :: load_txt => phi_wrap_load_txt
end type phi_flux_wrapper
!---Toroidal flux map internals
PRIVATE :: sort_unique, quad_interp, fd_deriv, tmap_find, tmap_index, tmap_lookup, torflux_qgeom, &
  phi_wrap_map, phi_wrap_norm, phi_wrap_relabel_int, phi_wrap_relabel_cof, phi_wrap_relabel_f, &
  phi_wrap_build_fint, torflux_set_bounds
contains
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine zero_save_hdf5(self,filename,path)
class(zero_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
IF(.NOT.hdf5_field_exist(filename,path//'/TYPE'))CALL hdf5_write('zero',filename,path//'/TYPE')
end subroutine zero_save_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine zero_save_txt(self,io_unit)
class(zero_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
WRITE(io_unit,*)'zero'
end subroutine zero_save_txt
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine zero_load_hdf5(self,filename,path,success)
class(zero_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
logical, intent(out) :: success
self%ndofs=0
success=.TRUE.
end subroutine zero_load_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine zero_load_txt(self,io_unit)
class(zero_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
self%ndofs=0
end subroutine zero_load_txt
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine zero_copy(self,new)
class(zero_flux_func), intent(inout) :: self
class(flux_func), pointer, intent(inout) :: new
ALLOCATE(zero_flux_func::new)
new%plasma_bounds=self%plasma_bounds
new%f_offset=self%f_offset
new%ndofs=0
end subroutine zero_copy
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine zero_delete(self)
class(zero_flux_func), intent(inout) :: self
!--- No internal memory to free
end subroutine zero_delete
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
class(gs_equil), intent(inout) :: gseq
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
!> Needs Docs
!------------------------------------------------------------------------------
subroutine flat_save_hdf5(self,filename,path)
class(flat_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
IF(.NOT.hdf5_field_exist(filename,path//'/TYPE'))CALL hdf5_write('flat',filename,path//'/TYPE')
end subroutine flat_save_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine flat_save_txt(self,io_unit)
class(flat_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
WRITE(io_unit,*)'flat'
end subroutine flat_save_txt
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine flat_load_hdf5(self,filename,path,success)
class(flat_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
logical, intent(out) :: success
self%ndofs=0
success=.TRUE.
end subroutine flat_load_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine flat_load_txt(self,io_unit)
class(flat_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
self%ndofs=0
end subroutine flat_load_txt
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_flat_f(func)
CLASS(flux_func), INTENT(inout) :: func
select type(self=>func)
  type is(flat_flux_func)
    self%ndofs=0
class default
  CALL oft_abort('Invalid flux function type in create_flat_f','create_flat_f',__FILE__)
end select
END SUBROUTINE create_flat_f
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine flat_copy(self,new)
class(flat_flux_func), intent(inout) :: self
class(flux_func), pointer, intent(inout) :: new
ALLOCATE(flat_flux_func::new)
new%plasma_bounds=self%plasma_bounds
new%f_offset=self%f_offset
new%ndofs=0
end subroutine flat_copy
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine flat_delete(self)
class(flat_flux_func), intent(inout) :: self
!--- No internal memory to free
end subroutine flat_delete
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
class(gs_equil), intent(inout) :: gseq
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
!> Needs Docs
!------------------------------------------------------------------------------
subroutine poly_save_hdf5(self,filename,path)
class(poly_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
IF(.NOT.hdf5_field_exist(filename,path//'/TYPE'))CALL hdf5_write('poly',filename,path//'/TYPE')
CALL hdf5_write(self%ndofs,filename,path//'/NCOFS')
CALL hdf5_write(self%cofs,filename,path//'/COFS')
CALL hdf5_write(self%zero_grad,filename,path//'/ZERO_GRAD')
end subroutine poly_save_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine poly_save_txt(self,io_unit)
class(poly_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
WRITE(io_unit,*)'poly'
IF(self%zero_grad)THEN
  WRITE(io_unit,*)self%ndofs+1,self%zero_grad
  WRITE(io_unit,*)0.d0,self%cofs
ELSE
  WRITE(io_unit,*)self%ndofs,self%zero_grad
  WRITE(io_unit,*)self%cofs
END IF
end subroutine poly_save_txt
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine poly_load_hdf5(self,filename,path,success)
class(poly_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
logical, intent(out) :: success
CALL hdf5_read(self%ndofs,filename,path//'/NCOFS',success=success)
IF(.NOT.success)RETURN
ALLOCATE(self%cofs(self%ndofs))
CALL hdf5_read(self%cofs,filename,path//'/COFS',success=success)
IF(.NOT.success)RETURN
CALL hdf5_read(self%zero_grad,filename,path//'/ZERO_GRAD',success=success)
IF(.NOT.success)RETURN
end subroutine poly_load_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine poly_load_txt(self,io_unit)
class(poly_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
logical :: zero_grad
integer(i4) :: ncofs
real(r8), allocatable :: cofs(:)
READ(io_unit,*)ncofs,zero_grad
ALLOCATE(cofs(ncofs))
READ(io_unit,*)cofs
CALL create_poly_ff(self,ncofs,cofs,zero_grad)
DEALLOCATE(cofs)
end subroutine poly_load_txt
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_poly_ff(func,ncofs,cofs,zero_grad)
CLASS(flux_func), INTENT(inout) :: func
INTEGER(4), INTENT(in) :: ncofs
REAL(8), INTENT(in) :: cofs(:)
LOGICAL, INTENT(in) :: zero_grad
INTEGER(4) :: i,offset

select type(self=>func)
  type is(poly_flux_func)
  !---
  self%deg=ncofs
  self%ndofs=ncofs
  self%zero_grad=zero_grad
  offset=0
  IF(self%zero_grad)THEN
    offset=1
    self%ndofs=ncofs-1
    self%deg=ncofs-1
  END IF
  ALLOCATE(self%cofs(self%ndofs))
  self%cofs=0.d0
  !---
  DO i=1,self%ndofs
    self%cofs(i)=cofs(offset+i)
  END DO
  self%c0=1.d0
  IF(.NOT.self%zero_grad)THEN
    self%c0=0.d0
    DO i=1,self%ndofs
      self%c0=self%c0+self%cofs(i)/REAL((i+1),8)
    END DO
    self%c0=1.d0-self%c0
  END IF
class default
  CALL oft_abort('Invalid flux function type in create_poly_ff','create_poly_ff',__FILE__)
end select
END SUBROUTINE create_poly_ff
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine poly_copy(self,new)
class(poly_flux_func), intent(inout) :: self
class(flux_func), pointer, intent(inout) :: new
ALLOCATE(new, MOLD=self)
SELECT TYPE(new)
  CLASS IS(poly_flux_func)
    new%plasma_bounds=self%plasma_bounds
    new%f_offset=self%f_offset
    new%deg=self%deg
    new%ndofs=self%ndofs
    new%zero_grad=self%zero_grad
    new%c0=self%c0
    ALLOCATE(new%cofs(new%ndofs))
    new%cofs=self%cofs
END SELECT
end subroutine poly_copy
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine poly_delete(self)
class(poly_flux_func), intent(inout) :: self
IF(ASSOCIATED(self%cofs))DEALLOCATE(self%cofs)
self%ndofs=0
self%deg=0
end subroutine poly_delete
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
class(gs_equil), intent(inout) :: gseq
self%plasma_bounds=gseq%plasma_bounds
end subroutine poly_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function poly_cofs_update(self,c) result(ierr)
class(poly_flux_func), intent(inout) :: self
real(8), intent(in) :: c(:)
INTEGER(4) :: i,ierr
DO i=1,self%ndofs
  self%cofs(i)=c(i)
END DO
IF(.NOT.self%zero_grad)THEN
  self%c0=0.d0
  DO i=1,self%ndofs
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
DO i=1,self%ndofs
  c(i)=self%cofs(i)
END DO
end subroutine poly_cofs_get
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine spline_save_hdf5(self,filename,path)
class(spline_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
IF(.NOT.hdf5_field_exist(filename,path//'/TYPE'))CALL hdf5_write('spline',filename,path//'/TYPE')
CALL hdf5_write(self%npsi,filename,path//'/NPSI')
CALL hdf5_write(self%func%xs(0:self%npsi-1),filename,path//'/XVALS')
CALL hdf5_write(self%func%fs(0:self%npsi-1,1),filename,path//'/YVALS')
end subroutine spline_save_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine spline_save_txt(self,io_unit)
class(spline_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
WRITE(io_unit,*)'spline'
WRITE(io_unit,*)self%npsi
WRITE(io_unit,*)self%func%xs(0:self%npsi-1)
WRITE(io_unit,*)self%func%fs(0:self%npsi-1,1)
end subroutine spline_save_txt
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine spline_load_hdf5(self,filename,path,success)
class(spline_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
logical, intent(out) :: success
integer(i4) :: npsi
real(r8), allocatable :: xvals(:),yvals(:)
CALL hdf5_read(npsi,filename,path//'/NPSI',success=success)
IF(.NOT.success)RETURN
ALLOCATE(xvals(npsi),yvals(npsi))
CALL hdf5_read(xvals,filename,path//'/XVALS',success=success)
IF(.NOT.success)RETURN
CALL hdf5_read(yvals,filename,path//'/YVALS',success=success)
IF(.NOT.success)RETURN
CALL create_spline_ff(self,npsi,xvals,yvals)
DEALLOCATE(xvals,yvals)
end subroutine spline_load_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine spline_load_txt(self,io_unit)
class(spline_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
integer(i4) :: npsi
real(r8), allocatable :: xvals(:),yvals(:)
READ(io_unit,*)npsi
ALLOCATE(xvals(npsi),yvals(npsi))
READ(io_unit,*)xvals
READ(io_unit,*)yvals
CALL create_spline_ff(self,npsi,xvals,yvals)
DEALLOCATE(xvals,yvals)
end subroutine spline_load_txt
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_spline_ff(func,npsi,psivals,yvals)
CLASS(flux_func), INTENT(inout) :: func
INTEGER(4), INTENT(in) :: npsi
REAL(8), INTENT(in) :: psivals(npsi),yvals(npsi)
REAL(8) :: psi1,psi2
REAL(8), ALLOCATABLE :: c(:)
INTEGER(4) :: i,ierr
select type(self=>func)
  type is(spline_flux_func)
  !---
  self%npsi=npsi
  self%ndofs=self%npsi
  ALLOCATE(self%fun_loc(omp_get_max_threads()))
  CALL spline_alloc(self%func,self%npsi-1,1)
  DO i=1,omp_get_max_threads()
    CALL spline_alloc(self%fun_loc(i),self%npsi-1,1)
  END DO
  !---
  DO i=1,self%npsi
    self%func%xs(i-1)=psivals(i)
    self%func%fs(i-1,1)=yvals(i)
  END DO
  self%xmin=self%func%xs(0)
  self%xmax=self%func%xs(self%npsi-1)
  !---
  ! WRITE(*,*)'Fitting',self%func%xs
  !CALL spline_fit(self%func,"extrap")
  CALL spline_fit(self%func,"not-a-knot")
  DO i=1,omp_get_max_threads()
    CALL spline_copy(self%func,self%fun_loc(i))
  END DO
  ALLOCATE(c(self%ndofs))
  CALL self%get_cofs(c)
  ierr=self%set_cofs(c)
  DEALLOCATE(c)
  ! WRITE(*,*)'Spline Created',self%ndofs,self%func%xs
class default
  CALL oft_abort('Invalid flux function type in create_spline_ff','create_spline_ff',__FILE__)
end select

END SUBROUTINE create_spline_ff
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine spline_func_copy(self,new)
class(spline_flux_func), intent(inout) :: self
class(flux_func), pointer, intent(inout) :: new
integer(i4) :: i
ALLOCATE(new, MOLD=self)
SELECT TYPE(new)
  CLASS IS(spline_flux_func)
    new%plasma_bounds=self%plasma_bounds
    new%f_offset=self%f_offset
    new%npsi=self%npsi
    new%ndofs=self%ndofs
    new%xmin=self%xmin
    new%xmax=self%xmax
    new%f1=self%f1
    new%fn=self%fn
    new%yp1=self%yp1
    new%ypn=self%ypn
    CALL spline_copy(self%func,new%func)
    ALLOCATE(new%fun_loc(omp_get_max_threads()))
    DO i=1,omp_get_max_threads()
      CALL spline_alloc(new%fun_loc(i),self%npsi-1,1)
      CALL spline_copy(new%func,new%fun_loc(i))
    END DO
END SELECT
end subroutine spline_func_copy
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine spline_func_delete(self)
class(spline_flux_func), intent(inout) :: self
integer(i4) :: i
self%npsi=0
self%ndofs=0
IF(ASSOCIATED(self%fun_loc))THEN
  CALL spline_dealloc(self%func)
  DO i=1,omp_get_max_threads()
    CALL spline_dealloc(self%fun_loc(i))
  END DO
  DEALLOCATE(self%fun_loc)
END IF
end subroutine spline_func_delete
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
class(gs_equil), intent(inout) :: gseq
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
!> Needs Docs
!------------------------------------------------------------------------------
subroutine linterp_save_hdf5(self,filename,path)
class(linterp_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
IF(.NOT.hdf5_field_exist(filename,path//'/TYPE'))CALL hdf5_write('linterp',filename,path//'/TYPE')
CALL hdf5_write(self%npsi,filename,path//'/NPSI')
CALL hdf5_write(self%x,filename,path//'/XVALS')
CALL hdf5_write(self%yp,filename,path//'/YVALS')
CALL hdf5_write(self%y0,filename,path//'/Y0')
end subroutine linterp_save_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine linterp_save_txt(self,io_unit)
class(linterp_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
WRITE(io_unit,*)'linterp'
WRITE(io_unit,*)self%npsi,self%y0
WRITE(io_unit,*)self%x
WRITE(io_unit,*)self%yp
end subroutine linterp_save_txt
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine linterp_load_hdf5(self,filename,path,success)
class(linterp_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
logical, intent(out) :: success
integer(i4) :: npsi
real(r8) :: y0
real(r8), allocatable :: xvals(:),yvals(:)
CALL hdf5_read(npsi,filename,path//'/NPSI',success=success)
IF(.NOT.success)RETURN
ALLOCATE(xvals(npsi),yvals(npsi))
CALL hdf5_read(xvals,filename,path//'/XVALS',success=success)
IF(.NOT.success)RETURN
CALL hdf5_read(yvals,filename,path//'/YVALS',success=success)
IF(.NOT.success)RETURN
CALL hdf5_read(y0,filename,path//'/Y0',success=success)
IF(.NOT.success)RETURN
CALL create_linterp_ff(self,npsi,xvals,yvals,y0)
DEALLOCATE(xvals,yvals)
end subroutine linterp_load_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine linterp_load_txt(self,io_unit)
class(linterp_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
integer(i4) :: npsi
real(r8) :: y0
real(r8), allocatable :: xvals(:),yvals(:)
READ(io_unit,*)npsi,y0
ALLOCATE(xvals(npsi),yvals(npsi))
READ(io_unit,*)xvals
READ(io_unit,*)yvals
CALL create_linterp_ff(self,npsi,xvals,yvals,y0)
DEALLOCATE(xvals,yvals)
end subroutine linterp_load_txt
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_linterp_ff(func,npsi,psivals,yvals,y0)
CLASS(flux_func), INTENT(inout) :: func
INTEGER(4), INTENT(in) :: npsi
REAL(8), INTENT(in) :: psivals(npsi)
REAL(8), INTENT(in) :: yvals(npsi)
REAL(8), INTENT(in) :: y0
INTEGER(4) :: i,ierr
SELECT TYPE(self=>func)
  TYPE IS(linterp_flux_func)
  !---
  self%npsi=npsi
  self%ndofs=self%npsi
  !---
  ALLOCATE(self%x(self%npsi))
  ALLOCATE(self%yp(self%npsi))
  ALLOCATE(self%y(self%npsi))
  !---
  ! IF(y0<1.d-8)CALL oft_abort("Support for y0 = 0 has been removed","create_linterp_ff",__FILE__)
  self%y0=y0
  DO i=1,self%npsi
    self%x(i)=psivals(i)
    self%yp(i)=yvals(i)
  END DO
  IF(self%y0<1.d-8)THEN
    self%ndofs=self%ndofs-1
    ierr=self%set_cofs(yvals(2:self%npsi))
  ELSE
    ierr=self%set_cofs(yvals)
  END IF
  IF(oft_debug_print(1))WRITE(*,*)'Linear interpolator Created',self%ndofs,self%x,self%y0
class default
  CALL oft_abort('Invalid flux function type in create_linterp_ff','create_linterp_ff',__FILE__)
END SELECT

END SUBROUTINE create_linterp_ff
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine linterp_copy(self,new)
class(linterp_flux_func), intent(inout) :: self
class(flux_func), pointer, intent(inout) :: new
ALLOCATE(new, MOLD=self)
SELECT TYPE(new)
  CLASS IS(linterp_flux_func)
    new%plasma_bounds=self%plasma_bounds
    new%f_offset=self%f_offset
    new%npsi=self%npsi
    new%ndofs=self%ndofs
    new%y0=self%y0
    ALLOCATE(new%x(new%npsi),new%yp(new%npsi),new%y(new%npsi))
    new%x=self%x
    new%yp=self%yp
    new%y=self%y
END SELECT
end subroutine linterp_copy
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine linterp_delete(self)
class(linterp_flux_func), intent(inout) :: self
self%npsi=0
self%ndofs=0
IF(ASSOCIATED(self%x))DEALLOCATE(self%x)
IF(ASSOCIATED(self%yp))DEALLOCATE(self%yp)
IF(ASSOCIATED(self%y))DEALLOCATE(self%y)
end subroutine linterp_delete
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
      exit
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
      exit
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
      exit
    end if
  end do
end if
end function linterp_fpp
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine linterp_update(self,gseq)
class(linterp_flux_func), intent(inout) :: self
class(gs_equil), intent(inout) :: gseq
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
IF(self%ndofs==self%npsi-1)offset=1
DO i=1,self%ndofs
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
IF(self%ndofs==self%npsi-1)offset=1
DO i=1,self%ndofs
  c(i)=self%yp(i+offset)
END DO
end subroutine linterp_cofs_get
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mlinterp_save_hdf5(self,filename,path)
class(mlinterp_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
IF(.NOT.hdf5_field_exist(filename,path//'/TYPE'))CALL hdf5_write('mlinterp',filename,path//'/TYPE')
CALL hdf5_write(self%nbasis,filename,path//'/NBASIS')
CALL hdf5_write(self%weights,filename,path//'/WEIGHTS')
CALL hdf5_write(self%npsi,filename,path//'/NPSI')
CALL hdf5_write(self%x,filename,path//'/XVALS')
CALL hdf5_write(self%yp_basis,filename,path//'/YVALS_BASIS')
end subroutine mlinterp_save_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mlinterp_save_txt(self,io_unit)
class(mlinterp_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
integer(4) :: i
WRITE(io_unit,*)'mlinterp'
WRITE(io_unit,*)self%npsi
WRITE(io_unit,*)self%x
WRITE(io_unit,*)self%nbasis
WRITE(io_unit,*)self%weights
DO i=1,self%nbasis
  WRITE(io_unit,*)self%yp_basis(:,i)
END DO
end subroutine mlinterp_save_txt
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mlinterp_load_hdf5(self,filename,path,success)
class(mlinterp_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
logical, intent(out) :: success
integer(i4) :: npsi,nbasis
real(r8), allocatable :: xvals(:),weights(:)
real(r8), allocatable :: yvals(:,:)
CALL hdf5_read(npsi,filename,path//'/NPSI',success=success)
IF(.NOT.success)RETURN
CALL hdf5_read(nbasis,filename,path//'/NBASIS',success=success)
IF(.NOT.success)RETURN
ALLOCATE(xvals(npsi),yvals(npsi+1,nbasis),weights(nbasis-1))
CALL hdf5_read(xvals,filename,path//'/XVALS',success=success)
IF(.NOT.success)RETURN
CALL hdf5_read(yvals,filename,path//'/YVALS_BASIS',success=success)
IF(.NOT.success)RETURN
CALL hdf5_read(weights,filename,path//'/WEIGHTS',success=success)
IF(.NOT.success)RETURN
CALL create_mlinterp_ff(self,npsi,xvals,nbasis,yvals,weights)
DEALLOCATE(xvals,yvals,weights)
end subroutine mlinterp_load_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mlinterp_load_txt(self,io_unit)
class(mlinterp_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
integer(i4) :: i,npsi,nbasis
real(r8), allocatable :: xvals(:),weights(:)
real(r8), allocatable :: yvals(:,:)
READ(io_unit,*)npsi
ALLOCATE(xvals(npsi))
READ(io_unit,*)xvals
READ(io_unit,*)nbasis
ALLOCATE(yvals(npsi+1,nbasis),weights(nbasis-1))
READ(io_unit,*)weights
DO i=1,nbasis
  READ(io_unit,*)yvals(:,i)
END DO
CALL create_mlinterp_ff(self,npsi,xvals,nbasis,yvals,weights)
DEALLOCATE(xvals,yvals,weights)
end subroutine mlinterp_load_txt
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_mlinterp_ff(func,npsi,psivals,nbasis,yvals,weights)
CLASS(flux_func), INTENT(inout) :: func
INTEGER(4), INTENT(in) :: npsi
REAL(8), INTENT(in) :: psivals(npsi)
INTEGER(4), INTENT(in) :: nbasis
REAL(8), INTENT(in) :: yvals(npsi+1,nbasis)
REAL(8), INTENT(in) :: weights(nbasis-1)
INTEGER(4) :: i,ierr
SELECT TYPE(self=>func)
  TYPE IS(mlinterp_flux_func)
  !
  self%npsi=npsi
  self%nbasis=nbasis
  self%ndofs=self%nbasis-1
  !
  ALLOCATE(self%x(self%npsi))
  ALLOCATE(self%yp_basis(self%npsi+1,nbasis))
  ALLOCATE(self%weights(nbasis-1))
  self%weights=weights
  ALLOCATE(self%yp(self%npsi))
  ALLOCATE(self%y(self%npsi))
  !---
  self%x=psivals
  self%yp_basis=yvals
  self%y0=self%yp_basis(1,1)
  self%yp=self%yp_basis(2:self%npsi+1,1)
  DO i=2,self%nbasis
    self%y0=self%y0+self%weights(i-1)*self%yp_basis(1,i)
    self%yp=self%yp+self%weights(i-1)*self%yp_basis(2:self%npsi+1,i)
  END DO
  ierr=self%set_cofs(self%weights)
  IF(oft_debug_print(1))WRITE(*,*)'Multi-linear interpolator Created',self%ndofs,self%x,self%y0
class default
  CALL oft_abort('Invalid flux function type in create_mlinterp_ff','create_mlinterp_ff',__FILE__)
END SELECT
END SUBROUTINE create_mlinterp_ff
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mlinterp_copy(self,new)
class(mlinterp_flux_func), intent(inout) :: self
class(flux_func), pointer, intent(inout) :: new
ALLOCATE(new, MOLD=self)
SELECT TYPE(new)
  CLASS IS(mlinterp_flux_func)
    new%plasma_bounds=self%plasma_bounds
    new%f_offset=self%f_offset
    new%npsi=self%npsi
    new%ndofs=self%ndofs
    new%y0=self%y0
    ALLOCATE(new%x(new%npsi),new%yp(new%npsi),new%y(new%npsi))
    new%x=self%x
    new%yp=self%yp
    new%y=self%y
    new%nbasis=self%nbasis
    ALLOCATE(new%weights(new%nbasis-1),new%yp_basis(new%npsi+1,new%nbasis))
    new%weights=self%weights
    new%yp_basis=self%yp_basis
END SELECT
end subroutine mlinterp_copy
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mlinterp_delete(self)
class(mlinterp_flux_func), intent(inout) :: self
self%npsi=0
self%ndofs=0
self%nbasis=0
IF(ASSOCIATED(self%x))DEALLOCATE(self%x)
IF(ASSOCIATED(self%yp))DEALLOCATE(self%yp)
IF(ASSOCIATED(self%y))DEALLOCATE(self%y)
IF(ASSOCIATED(self%weights))DEALLOCATE(self%weights)
IF(ASSOCIATED(self%yp_basis))DEALLOCATE(self%yp_basis)
end subroutine mlinterp_delete
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine mlinterp_update(self,gseq)
class(mlinterp_flux_func), intent(inout) :: self
class(gs_equil), intent(inout) :: gseq
self%plasma_bounds=gseq%plasma_bounds
end subroutine mlinterp_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
function mlinterp_cofs_update(self,c) result(ierr)
class(mlinterp_flux_func), intent(inout) :: self
real(8), intent(in) :: c(:)
real(8) :: x,xp
integer(4) :: ierr
integer(4) :: i
!---
self%weights=c
self%y0=self%yp_basis(1,1)
self%yp=self%yp_basis(2:self%npsi+1,1)
DO i=2,self%nbasis
  self%y0=self%y0+self%weights(i-1)*self%yp_basis(1,i)
  self%yp=self%yp+self%weights(i-1)*self%yp_basis(2:self%npsi+1,i)
END DO
! WRITE(*,*)'mLinterp set: ',self%y0
! WRITE(*,*)'  yp: ',self%yp
!---
DO i=1,self%npsi
  x=self%x(i)
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
  WRITE(*,'(2A)')oft_indent,'Update Multi-linear interpolator:'
  CALL oft_increase_indent
  WRITE(*,'(2A,100ES11.3)')oft_indent,' x  =',self%x
  WRITE(*,'(2A,100ES11.3)')oft_indent,' yp =',self%yp
  WRITE(*,'(2A,100ES11.3)')oft_indent,' y  =',self%y
  WRITE(*,'(2A,ES11.3)')oft_indent,' y0  =',self%y0
  CALL oft_decrease_indent
END IF
ierr=0
end function mlinterp_cofs_update
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine mlinterp_cofs_get(self,c)
class(mlinterp_flux_func), intent(inout) :: self
real(8), intent(out) :: c(:)
c=self%weights
end subroutine mlinterp_cofs_get
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine wesson_save_hdf5(self,filename,path)
class(wesson_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
IF(.NOT.hdf5_field_exist(filename,path//'/TYPE'))CALL hdf5_write('wesson',filename,path//'/TYPE')
CALL hdf5_write(self%ndofs,filename,path//'/NCOFS')
CALL hdf5_write(self%gamma,filename,path//'/GAMMA')
end subroutine wesson_save_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine wesson_save_txt(self,io_unit)
class(wesson_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
WRITE(io_unit,*)'wesson'
WRITE(io_unit,*)self%ndofs
WRITE(io_unit,*)self%gamma
end subroutine wesson_save_txt
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine wesson_load_hdf5(self,filename,path,success)
class(wesson_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
logical, intent(out) :: success
CALL hdf5_read(self%ndofs,filename,path//'/NCOFS',success=success)
IF(.NOT.success)RETURN
CALL hdf5_read(self%gamma,filename,path//'/GAMMA',success=success)
IF(.NOT.success)RETURN
end subroutine wesson_load_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine wesson_load_txt(self,io_unit)
class(wesson_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
READ(io_unit,*)self%ndofs
READ(io_unit,*)self%gamma
end subroutine wesson_load_txt
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_wesson_ff(func,ncofs,gamma)
CLASS(flux_func), INTENT(inout) :: func
INTEGER(4), INTENT(in) :: ncofs
REAL(8), INTENT(in) :: gamma
select type(self=>func)
  type is(wesson_flux_func)
  !---
  self%ndofs=ncofs
  self%gamma=gamma
class default
  CALL oft_abort('Invalid flux function type in create_wesson_ff','create_wesson_ff',__FILE__)
end select

END SUBROUTINE create_wesson_ff
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine wesson_copy(self,new)
class(wesson_flux_func), intent(inout) :: self
class(flux_func), pointer, intent(inout) :: new
ALLOCATE(wesson_flux_func::new)
SELECT TYPE(new)
  TYPE IS(wesson_flux_func)
    new%plasma_bounds=self%plasma_bounds
    new%f_offset=self%f_offset
    new%ndofs=self%ndofs
    new%gamma=self%gamma
END SELECT
end subroutine wesson_copy
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine wesson_delete(self)
class(wesson_flux_func), intent(inout) :: self
self%ndofs=0
self%gamma=0.d0
end subroutine wesson_delete
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
class(gs_equil), intent(inout) :: gseq
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
!------------------------------------------------------------------------------
!> Wrap `F` as a profile on normalized toroidal flux (`coord`: 1 "phi_n", 2 "phi_n_relabel")
!------------------------------------------------------------------------------
subroutine torflux_wrap(F,coord)
class(flux_func), pointer, intent(inout) :: F !< Profile on \f$ \hat{\Phi} \f$, replaced by its wrapper
integer(4), intent(in) :: coord !< Coordinate of `F`
class(phi_flux_wrapper), pointer :: wrap
IF(coord==0)RETURN
SELECT TYPE(F)
  CLASS IS(spline_flux_func)
    CALL oft_abort('Toroidal flux coordinate not supported for this profile','torflux_wrap',__FILE__)
END SELECT
ALLOCATE(wrap)
wrap%relabel=(coord==2)
wrap%phi_rep=>F
wrap%ndofs=F%ndofs
wrap%f_offset=F%f_offset
wrap%plasma_bounds=F%plasma_bounds
wrap%update_on_load=F%update_on_load
F=>wrap
end subroutine torflux_wrap
!------------------------------------------------------------------------------
!> Coordinate of `F` (0 "psi_n", 1 "phi_n", 2 "phi_n_relabel")
!------------------------------------------------------------------------------
function torflux_coord(F) result(coord)
class(flux_func), intent(in) :: F
integer(4) :: coord
coord=0
SELECT TYPE(F)
  TYPE IS(phi_flux_wrapper)
    coord=MERGE(2,1,F%relabel)
END SELECT
end function torflux_coord
!------------------------------------------------------------------------------
!> Profile `F` itself, or the profile it wraps if it is on toroidal flux
!------------------------------------------------------------------------------
function torflux_unwrap(F) result(base)
class(flux_func), pointer, intent(in) :: F
class(flux_func), pointer :: base
base=>F
IF(.NOT.ASSOCIATED(F))RETURN
SELECT TYPE(F)
  TYPE IS(phi_flux_wrapper)
    base=>F%phi_rep
END SELECT
end function torflux_unwrap
!------------------------------------------------------------------------------
!> Point a profile on toroidal flux at the equilibrium's map, creating the map if needed
!------------------------------------------------------------------------------
subroutine torflux_attach(F,gseq)
class(flux_func), pointer, intent(in) :: F !< Profile (no-op unless it is a @ref phi_flux_wrapper)
class(gs_equil), intent(inout) :: gseq !< G-S object
IF(.NOT.ASSOCIATED(F))RETURN
SELECT TYPE(F)
  TYPE IS(phi_flux_wrapper)
    IF(.NOT.ASSOCIATED(gseq%torflux))ALLOCATE(torflux_map::gseq%torflux)
    SELECT TYPE(map=>gseq%torflux)
      CLASS IS(torflux_map)
        F%map=>map
    END SELECT
    SELECT TYPE(rep=>F%phi_rep)
      CLASS IS(linterp_flux_func)
        rep%phi_wrapper=>F
    END SELECT
END SELECT
end subroutine torflux_attach
!------------------------------------------------------------------------------
!> Attach all profiles of an equilibrium to its toroidal flux map (created if needed)
!------------------------------------------------------------------------------
subroutine torflux_setup(gseq,active,resettle)
class(gs_equil), intent(inout) :: gseq !< G-S object
logical, optional, intent(out) :: active !< Any profile on toroidal flux?
logical, optional, intent(in) :: resettle !< Settle surface placement fully on the next map update (new profiles or state)
logical :: any_wrapped
any_wrapped=.FALSE.
CALL attach(gseq%I,'',.FALSE.)
CALL attach(gseq%P,'',.FALSE.)
CALL attach(gseq%I_NI,'',.FALSE.)
CALL attach(gseq%eta,'eta',.TRUE.)
CALL attach(gseq%Te,'Te',.TRUE.)
CALL attach(gseq%Ti,'Ti',.TRUE.)
CALL attach(gseq%ne,'ne',.TRUE.)
CALL attach(gseq%ni,'ni',.TRUE.)
CALL attach(gseq%Zeff,'Zeff',.TRUE.)
CALL attach(gseq%jphi_fixed,'jphi_fixed',.TRUE.)
CALL attach(gseq%jphi_saw,'jphi_saw',.TRUE.)
IF(PRESENT(active))active=any_wrapped
IF(PRESENT(resettle).AND.ASSOCIATED(gseq%torflux))THEN
  SELECT TYPE(map=>gseq%torflux)
    CLASS IS(torflux_map)
      IF(resettle)map%settle=.TRUE.
  END SELECT
END IF
CONTAINS
!---Attach a wrapped profile (value-type profiles must be relabelled)
subroutine attach(F,name,value_type)
CLASS(flux_func), POINTER, INTENT(in) :: F
CHARACTER(LEN=*), INTENT(in) :: name
LOGICAL, INTENT(in) :: value_type
IF(.NOT.ASSOCIATED(F))RETURN
SELECT TYPE(F)
  TYPE IS(phi_flux_wrapper)
    IF(value_type.AND.(.NOT.F%relabel))CALL oft_abort(name//' profile on toroidal flux must use "phi_n_relabel"', &
      'torflux_setup',__FILE__)
    any_wrapped=.TRUE.
END SELECT
CALL torflux_attach(F,gseq)
end subroutine attach
end subroutine torflux_setup
!------------------------------------------------------------------------------
!> Map \f$ \psi \f$ to the coordinate of `phi_rep`, returns `.FALSE.` if no map applies
!------------------------------------------------------------------------------
function phi_wrap_map(self,psi,psim,jac,djac) result(mapped)
class(phi_flux_wrapper), intent(inout) :: self
real(8), intent(in) :: psi !< Poloidal flux
real(8), intent(out) :: psim !< Mapped flux \f$ \psi_a + (1-\hat{\Phi})(\psi_0-\psi_a) \f$
real(8), intent(out) :: jac !< \f$ d\hat{\Phi}/d\hat{\psi} \f$
real(8), optional, intent(out) :: djac !< \f$ d^2\hat{\Phi}/d\hat{\psi}^2 \f$
logical :: mapped
real(8) :: x1,dpsi,psihat,phihat
mapped=.FALSE.
psim=psi; jac=1.d0
IF(PRESENT(djac))djac=0.d0
IF(.NOT.ASSOCIATED(self%map))RETURN
IF(self%map%ns<2)RETURN
CALL phi_wrap_norm(self,x1,dpsi)
psihat=(psi-x1)/dpsi
IF(psihat<0.d0)RETURN
CALL self%map%eval(psihat,phihat,jac,djac)
psim=x1+phihat*dpsi
mapped=.TRUE.
end function phi_wrap_map
!------------------------------------------------------------------------------
!> Set the plasma bounds of `F`, and of the profile it wraps if it is on toroidal flux
!------------------------------------------------------------------------------
subroutine torflux_set_bounds(F,bounds)
class(flux_func), intent(inout) :: F !< Profile
real(8), intent(in) :: bounds(2) !< Plasma bounds
F%plasma_bounds=bounds
SELECT TYPE(F)
  TYPE IS(phi_flux_wrapper)
    F%phi_rep%plasma_bounds=bounds
END SELECT
end subroutine torflux_set_bounds
!------------------------------------------------------------------------------
!> Normalization offset and width (unset bounds: argument is already \f$ 1-\hat{\psi} \f$)
!------------------------------------------------------------------------------
pure subroutine phi_wrap_norm(self,x1,dpsi)
class(flux_func), intent(in) :: self
real(8), intent(out) :: x1 !< Flux at LCFS
real(8), intent(out) :: dpsi !< Axis minus LCFS flux
IF(self%plasma_bounds(1)<-1.d98)THEN
  x1=0.d0; dpsi=1.d0
ELSE
  x1=self%plasma_bounds(1); dpsi=self%plasma_bounds(2)-self%plasma_bounds(1)
END IF
end subroutine phi_wrap_norm
!------------------------------------------------------------------------------
!> Evaluate the profile at \f$ \psi \f$
!------------------------------------------------------------------------------
function phi_wrap_f(self,psi) result(b)
class(phi_flux_wrapper), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b,psim,jac
IF(.NOT.phi_wrap_map(self,psi,psim,jac))THEN
  b=self%phi_rep%f(psi)
ELSE IF(.NOT.self%relabel)THEN
  b=self%phi_rep%f(psim)
ELSE
  b=phi_wrap_relabel_f(self,psi)
END IF
end function phi_wrap_f
!------------------------------------------------------------------------------
!> Evaluate the first derivative of the profile at \f$ \psi \f$
!------------------------------------------------------------------------------
function phi_wrap_fp(self,psi) result(b)
class(phi_flux_wrapper), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b,psim,jac
IF(.NOT.phi_wrap_map(self,psi,psim,jac))THEN
  b=self%phi_rep%fp(psi)
ELSE IF(.NOT.self%relabel)THEN
  b=self%phi_rep%fp(psim)*jac
ELSE
  b=self%phi_rep%fp(psim)
END IF
end function phi_wrap_fp
!------------------------------------------------------------------------------
!> Evaluate the second derivative of the profile at \f$ \psi \f$
!------------------------------------------------------------------------------
function phi_wrap_fpp(self,psi) result(b)
class(phi_flux_wrapper), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b,psim,jac,djac,x1,dpsi
IF(.NOT.phi_wrap_map(self,psi,psim,jac,djac))THEN
  b=self%phi_rep%fpp(psi)
ELSE IF(.NOT.self%relabel)THEN
  CALL phi_wrap_norm(self,x1,dpsi)
  b=self%phi_rep%fpp(psim)*jac**2 + self%phi_rep%fp(psim)*djac/dpsi
ELSE
  b=self%phi_rep%fpp(psim)*jac
END IF
end function phi_wrap_fpp
!------------------------------------------------------------------------------
!> Integral of `phi_rep%fp` over \f$ 1-\hat{\psi} \in [a,b] \f$ (4-point Gauss)
!------------------------------------------------------------------------------
function phi_wrap_relabel_int(self,a,b) result(val)
class(phi_flux_wrapper), intent(inout) :: self
real(8), intent(in) :: a,b
real(8) :: val,s,phihat,jac,x1,dpsi
real(8), parameter :: xg(4)=[-0.861136311594053d0,-0.339981043584856d0, &
  0.339981043584856d0,0.861136311594053d0]
real(8), parameter :: wg(4)=[0.347854845137454d0,0.652145154862546d0, &
  0.652145154862546d0,0.347854845137454d0]
integer(4) :: i
CALL phi_wrap_norm(self,x1,dpsi)
val=0.d0
DO i=1,4
  s=(a+b)/2.d0+xg(i)*(b-a)/2.d0
  CALL self%map%eval(s,phihat,jac)
  val=val+wg(i)*self%phi_rep%fp(x1+phihat*dpsi)
END DO
val=val*(b-a)/2.d0
end function phi_wrap_relabel_int
!------------------------------------------------------------------------------
!> Rebuild the cumulative integral table on the map surfaces (relabel only)
!------------------------------------------------------------------------------
subroutine phi_wrap_build_fint(self)
class(phi_flux_wrapper), intent(inout) :: self
integer(4) :: k
IF(ASSOCIATED(self%fint))DEALLOCATE(self%fint)
IF(ASSOCIATED(self%fcof))DEALLOCATE(self%fcof)
self%fint_stamp=-1
IF(self%map%ns<2)RETURN
ALLOCATE(self%fint(self%map%ns),self%fcof(4,self%map%ns-1))
self%fint(1)=phi_wrap_relabel_int(self,0.d0,self%map%psis(1))
DO k=2,self%map%ns
  CALL phi_wrap_relabel_cof(self,self%map%psis(k-1),self%map%psis(k),self%fcof(:,k-1))
  self%fint(k)=self%fint(k-1)+SUM(self%fcof(:,k-1))
END DO
self%fint_stamp=self%map%stamp
end subroutine phi_wrap_build_fint
!------------------------------------------------------------------------------
!> Quartic \f$ \sum_j c_j u^j \f$ (u in [0,1]) for \f$ \int_a^{a+u(b-a)} f'\,d(1-\hat{\psi}) \f$
!!
!! Exact integral of the cubic interpolating `phi_rep%fp` at the 4 Gauss points of [a,b],
!! so the full-interval value equals @ref phi_wrap_relabel_int.
!------------------------------------------------------------------------------
subroutine phi_wrap_relabel_cof(self,a,b,cof)
class(phi_flux_wrapper), intent(inout) :: self
real(8), intent(in) :: a,b
real(8), intent(out) :: cof(4)
real(8), parameter :: xg(4)=[-0.861136311594053d0,-0.339981043584856d0, &
  0.339981043584856d0,0.861136311594053d0]
real(8) :: u(4),fv(4),m(4),phihat,jac,x1,dpsi,s
integer(4) :: i,j,k
CALL phi_wrap_norm(self,x1,dpsi)
u=(1.d0+xg)/2.d0
DO i=1,4
  s=a+u(i)*(b-a)
  CALL self%map%eval(s,phihat,jac)
  fv(i)=self%phi_rep%fp(x1+phihat*dpsi)
END DO
!---Monomial coefficients of the Lagrange interpolant through (u,fv)
cof=0.d0
DO i=1,4
  m=0.d0; m(1)=1.d0
  DO j=1,4
    IF(j==i)CYCLE
    DO k=4,2,-1
      m(k)=(m(k-1)-u(j)*m(k))/(u(i)-u(j))
    END DO
    m(1)=-u(j)*m(1)/(u(i)-u(j))
  END DO
  cof=cof+fv(i)*m
END DO
!---Integrate: c_j u^j -> c_j u^(j+1)/(j+1), scaled by interval width
cof=cof*(b-a)/[1.d0,2.d0,3.d0,4.d0]
end subroutine phi_wrap_relabel_cof
!------------------------------------------------------------------------------
!> Evaluate \f$ (\psi_0-\psi_a) \int f'\,d(1-\hat{\psi}) \f$ from the LCFS to \f$ 1-\hat{\psi} \f$ (relabel)
!------------------------------------------------------------------------------
function phi_wrap_relabel_f(self,psi) result(b)
class(phi_flux_wrapper), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b,x1,dpsi,psihat,u
integer(4) :: lo
CALL phi_wrap_norm(self,x1,dpsi)
psihat=(psi-x1)/dpsi
IF(self%fint_stamp/=self%map%stamp)THEN
  !---Table unavailable, integrate over all map intervals (slow)
  b=phi_wrap_relabel_int(self,0.d0,MIN(psihat,self%map%psis(1)))
  DO lo=2,self%map%ns
    IF(psihat<=self%map%psis(lo-1))EXIT
    b=b+phi_wrap_relabel_int(self,self%map%psis(lo-1),MIN(psihat,self%map%psis(lo)))
  END DO
  IF(psihat>self%map%psis(self%map%ns))b=b+phi_wrap_relabel_int(self,self%map%psis(self%map%ns),psihat)
ELSE IF(psihat<=self%map%psis(1))THEN
  b=phi_wrap_relabel_int(self,0.d0,psihat)
ELSE
  lo=tmap_lookup(self%map,psihat)
  u=(psihat-self%map%psis(lo))/(self%map%psis(lo+1)-self%map%psis(lo))
  b=self%fint(lo)+u*(self%fcof(1,lo)+u*(self%fcof(2,lo)+u*(self%fcof(3,lo)+u*self%fcof(4,lo))))
END IF
b=b*dpsi
end function phi_wrap_relabel_f
!------------------------------------------------------------------------------
!> Update the profile on \f$ \hat{\Phi} \f$ for the new equilibrium
!------------------------------------------------------------------------------
subroutine phi_wrap_update(self,gseq)
class(phi_flux_wrapper), intent(inout) :: self
class(gs_equil), intent(inout) :: gseq
real(8) :: dpsi
CALL torflux_set_bounds(self,gseq%plasma_bounds)
dpsi=ABS(gseq%plasma_bounds(2)-gseq%plasma_bounds(1))
IF(.NOT.((dpsi>0.d0).AND.(dpsi<HUGE(dpsi))))RETURN ! No closed surfaces to map
IF(.NOT.ASSOCIATED(self%map))THEN
  IF(.NOT.ASSOCIATED(gseq%torflux))ALLOCATE(torflux_map::gseq%torflux)
  SELECT TYPE(map=>gseq%torflux)
    CLASS IS(torflux_map)
      self%map=>map
  END SELECT
END IF
CALL self%phi_rep%update(gseq)
self%ndofs=self%phi_rep%ndofs
IF(self%relabel)CALL phi_wrap_build_fint(self)
end subroutine phi_wrap_update
!------------------------------------------------------------------------------
!> Set the coefficients of `phi_rep`
!------------------------------------------------------------------------------
function phi_wrap_cofs_update(self,c) result(ierr)
class(phi_flux_wrapper), intent(inout) :: self
real(8), intent(in) :: c(:)
integer(4) :: ierr
ierr=self%phi_rep%set_cofs(c)
IF(self%relabel.AND.ASSOCIATED(self%map))CALL phi_wrap_build_fint(self)
end function phi_wrap_cofs_update
!------------------------------------------------------------------------------
!> Get the coefficients of `phi_rep`
!------------------------------------------------------------------------------
subroutine phi_wrap_cofs_get(self,c)
class(phi_flux_wrapper), intent(inout) :: self
real(8), intent(out) :: c(:)
CALL self%phi_rep%get_cofs(c)
end subroutine phi_wrap_cofs_get
!------------------------------------------------------------------------------
!> Destroy the wrapper and `phi_rep`
!------------------------------------------------------------------------------
subroutine phi_wrap_delete(self)
class(phi_flux_wrapper), intent(inout) :: self
IF(ASSOCIATED(self%phi_rep))THEN
  CALL self%phi_rep%delete()
  DEALLOCATE(self%phi_rep)
END IF
IF(ASSOCIATED(self%fint))DEALLOCATE(self%fint)
IF(ASSOCIATED(self%fcof))DEALLOCATE(self%fcof)
NULLIFY(self%map)
self%fint_stamp=-1
end subroutine phi_wrap_delete
!------------------------------------------------------------------------------
!> Copy the wrapper and `phi_rep` (the copy shares the map until @ref torflux_setup)
!------------------------------------------------------------------------------
subroutine phi_wrap_copy(self,new)
class(phi_flux_wrapper), intent(inout) :: self
class(flux_func), pointer, intent(inout) :: new
class(flux_func), pointer :: rep
CALL self%phi_rep%copy(rep)
CALL torflux_wrap(rep,MERGE(2,1,self%relabel))
SELECT TYPE(rep)
  TYPE IS(phi_flux_wrapper)
    rep%plasma_bounds=self%plasma_bounds
    rep%f_offset=self%f_offset
    rep%map=>self%map
    IF(ASSOCIATED(self%fint))THEN
      ALLOCATE(rep%fint,SOURCE=self%fint)
      ALLOCATE(rep%fcof,SOURCE=self%fcof)
      rep%fint_stamp=self%fint_stamp
    END IF
END SELECT
new=>rep
end subroutine phi_wrap_copy
!------------------------------------------------------------------------------
!> Save `phi_rep` and its coordinate to HDF5
!------------------------------------------------------------------------------
subroutine phi_wrap_save_hdf5(self,filename,path)
class(phi_flux_wrapper), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
CALL self%phi_rep%save(filename,path)
CALL hdf5_write(torflux_coord(self),filename,path//'/COORD')
end subroutine phi_wrap_save_hdf5
!------------------------------------------------------------------------------
!> Save `phi_rep` to a text file, with its coordinate after the type on the first line
!------------------------------------------------------------------------------
subroutine phi_wrap_save_txt(self,io_unit)
class(phi_flux_wrapper), intent(inout) :: self
integer, intent(in) :: io_unit
character(LEN=2**17) :: line
integer(4) :: tmp_unit,ios
logical :: first
OPEN(NEWUNIT=tmp_unit,STATUS='SCRATCH')
CALL self%phi_rep%save(tmp_unit)
REWIND(tmp_unit)
first=.TRUE.
DO
  READ(tmp_unit,'(A)',IOSTAT=ios)line
  IF(ios/=0)EXIT
  IF(LEN_TRIM(line)==LEN(line))CALL oft_abort('Profile line too long','phi_wrap_save_txt',__FILE__)
  IF(first)line=TRIM(line)//MERGE(' phi_n_relabel',' phi_n        ',self%relabel)
  WRITE(io_unit,'(A)')TRIM(line)
  first=.FALSE.
END DO
CLOSE(tmp_unit)
end subroutine phi_wrap_save_txt
!------------------------------------------------------------------------------
!> Load `phi_rep` from HDF5 (allocate it first; normally loaded unwrapped, see @ref torflux_wrap)
!------------------------------------------------------------------------------
subroutine phi_wrap_load_hdf5(self,filename,path,success)
class(phi_flux_wrapper), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
logical, intent(out) :: success
CALL self%phi_rep%load(filename,path,success)
self%ndofs=self%phi_rep%ndofs
end subroutine phi_wrap_load_hdf5
!------------------------------------------------------------------------------
!> Load `phi_rep` from a text file (allocate it first; normally loaded unwrapped, see @ref torflux_wrap)
!------------------------------------------------------------------------------
subroutine phi_wrap_load_txt(self,io_unit)
class(phi_flux_wrapper), intent(inout) :: self
integer, intent(in) :: io_unit
CALL self%phi_rep%load(io_unit)
self%ndofs=self%phi_rep%ndofs
end subroutine phi_wrap_load_txt
!------------------------------------------------------------------------------
!> Evaluate \f$ 1-\hat{\Phi} \f$, \f$ J = d\hat{\Phi}/d\hat{\psi} \f$ and \f$ dJ/d\hat{\psi} \f$ at \f$ 1-\hat{\psi} \f$
!------------------------------------------------------------------------------
subroutine tmap_eval(self,psihat,phihat,jac,djac)
class(torflux_map), intent(in) :: self
real(8), intent(in) :: psihat !< Normalized poloidal flux
real(8), intent(out) :: phihat !< Normalized toroidal flux
real(8), intent(out) :: jac !< \f$ d\hat{\Phi}/d\hat{\psi} \f$
real(8), optional, intent(out) :: djac !< \f$ d^2\hat{\Phi}/d\hat{\psi}^2 \f$
integer(4) :: lo,hi
real(8) :: h,u,pa,pb,ma,mb
IF(self%ns<2)THEN
  phihat=psihat; jac=1.d0
  IF(PRESENT(djac))djac=0.d0
  RETURN
END IF
!---Linear extrapolation outside surface range
IF(psihat<=self%psis(1))THEN
  phihat=self%phihat(1)+(psihat-self%psis(1))*self%jac(1); jac=self%jac(1)
  IF(PRESENT(djac))djac=0.d0
  RETURN
ELSE IF(psihat>=self%psis(self%ns))THEN
  phihat=self%phihat(self%ns)+(psihat-self%psis(self%ns))*self%jac(self%ns); jac=self%jac(self%ns)
  IF(PRESENT(djac))djac=0.d0
  RETURN
END IF
lo=tmap_lookup(self,psihat)
hi=lo+1
!---Cubic Hermite on [lo,hi]
h=self%psis(hi)-self%psis(lo)
u=(psihat-self%psis(lo))/h
pa=self%phihat(lo); pb=self%phihat(hi)
ma=self%jac(lo)*h; mb=self%jac(hi)*h
phihat=(2.d0*u**3-3.d0*u**2+1.d0)*pa + (u**3-2.d0*u**2+u)*ma &
  + (3.d0*u**2-2.d0*u**3)*pb + (u**3-u**2)*mb
jac=((6.d0*u**2-6.d0*u)*(pa-pb) + (3.d0*u**2-4.d0*u+1.d0)*ma + (3.d0*u**2-2.d0*u)*mb)/h
IF(PRESENT(djac))djac=((12.d0*u-6.d0)*(pa-pb) + (6.d0*u-4.d0)*ma + (6.d0*u-2.d0)*mb)/h**2
end subroutine tmap_eval
!------------------------------------------------------------------------------
!> Evaluate \f$ 1-\hat{\psi} \f$ at \f$ 1-\hat{\Phi} \f$ (inverse of @ref tmap_eval)
!------------------------------------------------------------------------------
function tmap_inv(self,phihat) result(psihat)
class(torflux_map), intent(in) :: self
real(8), intent(in) :: phihat !< Normalized toroidal flux
real(8) :: psihat
real(8) :: a,b,phi,jac
integer(4) :: lo,i
IF(self%ns<2)THEN
  psihat=phihat
  RETURN
END IF
IF(phihat<=self%phihat(1))THEN
  psihat=self%psis(1)+(phihat-self%phihat(1))/self%jac(1)
  RETURN
ELSE IF(phihat>=self%phihat(self%ns))THEN
  psihat=self%psis(self%ns)+(phihat-self%phihat(self%ns))/self%jac(self%ns)
  RETURN
END IF
lo=tmap_find(self%phihat,self%ns,phihat)
a=self%psis(lo); b=self%psis(lo+1)
psihat=a+(b-a)*(phihat-self%phihat(lo))/(self%phihat(lo+1)-self%phihat(lo))
!---Safeguarded Newton
DO i=1,30
  CALL self%eval(psihat,phi,jac)
  IF(phi>phihat)THEN
    b=psihat
  ELSE
    a=psihat
  END IF
  IF(ABS(phi-phihat)<1.d-14)EXIT
  psihat=psihat-(phi-phihat)/jac
  IF((psihat<=a).OR.(psihat>=b))psihat=(a+b)/2.d0
END DO
end function tmap_inv
!------------------------------------------------------------------------------
!> Reset map to identity
!------------------------------------------------------------------------------
subroutine tmap_delete(self)
class(torflux_map), intent(inout) :: self
IF(ASSOCIATED(self%psis))DEALLOCATE(self%psis)
IF(ASSOCIATED(self%phihat))DEALLOCATE(self%phihat)
IF(ASSOCIATED(self%jac))DEALLOCATE(self%jac)
IF(ASSOCIATED(self%nodes))DEALLOCATE(self%nodes)
IF(ASSOCIATED(self%bin))DEALLOCATE(self%bin)
self%ns=0
self%Q=1.d0
self%stamp=self%stamp+1
end subroutine tmap_delete
!------------------------------------------------------------------------------
!> Copy map (profiles of the new equilibrium attach to it in @ref torflux_setup)
!------------------------------------------------------------------------------
subroutine tmap_copy(self,new)
class(torflux_map), intent(inout) :: self
class(flux_func), pointer, intent(inout) :: new
ALLOCATE(torflux_map::new)
SELECT TYPE(new)
  TYPE IS(torflux_map)
    new%plasma_bounds=self%plasma_bounds
    new%active=self%active
    new%Q=self%Q
    new%ns=self%ns
    new%stamp=self%stamp
    IF(self%ns>0)THEN
      ALLOCATE(new%psis,SOURCE=self%psis)
      ALLOCATE(new%phihat,SOURCE=self%phihat)
      ALLOCATE(new%jac,SOURCE=self%jac)
      CALL tmap_index(new)
    END IF
    IF(ASSOCIATED(self%nodes))ALLOCATE(new%nodes,SOURCE=self%nodes)
END SELECT
end subroutine tmap_copy
!------------------------------------------------------------------------------
!> Evaluate \f$ 1-\hat{\Phi} \f$ at \f$ \psi \f$
!------------------------------------------------------------------------------
function tmap_f(self,psi) result(b)
class(torflux_map), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b,x1,dpsi,jac
CALL phi_wrap_norm(self,x1,dpsi)
CALL self%eval((psi-x1)/dpsi,b,jac)
end function tmap_f
!------------------------------------------------------------------------------
!> Evaluate \f$ d(1-\hat{\Phi})/d\psi \f$ at \f$ \psi \f$
!------------------------------------------------------------------------------
function tmap_fp(self,psi) result(b)
class(torflux_map), intent(inout) :: self
real(8), intent(in) :: psi
real(8) :: b,x1,dpsi,phihat
CALL phi_wrap_norm(self,x1,dpsi)
CALL self%eval((psi-x1)/dpsi,phihat,b)
b=b/dpsi
end function tmap_fp
!------------------------------------------------------------------------------
!> The map has no free coefficients
!------------------------------------------------------------------------------
function tmap_cofs_update(self,c) result(ierr)
class(torflux_map), intent(inout) :: self
real(8), intent(in) :: c(:)
integer(4) :: ierr
ierr=0
end function tmap_cofs_update
!------------------------------------------------------------------------------
!> The map has no free coefficients
!------------------------------------------------------------------------------
subroutine tmap_cofs_get(self,c)
class(torflux_map), intent(inout) :: self
real(8), intent(out) :: c(:)
end subroutine tmap_cofs_get
!------------------------------------------------------------------------------
!> The map is rebuilt from the equilibrium, not saved
!------------------------------------------------------------------------------
subroutine tmap_save_hdf5(self,filename,path)
class(torflux_map), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
end subroutine tmap_save_hdf5
!------------------------------------------------------------------------------
!> The map is rebuilt from the equilibrium, not saved
!------------------------------------------------------------------------------
subroutine tmap_save_txt(self,io_unit)
class(torflux_map), intent(inout) :: self
integer, intent(in) :: io_unit
end subroutine tmap_save_txt
!------------------------------------------------------------------------------
!> The map is rebuilt from the equilibrium, not loaded
!------------------------------------------------------------------------------
subroutine tmap_load_hdf5(self,filename,path,success)
class(torflux_map), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
logical, intent(out) :: success
success=.FALSE.
end subroutine tmap_load_hdf5
!------------------------------------------------------------------------------
!> The map is rebuilt from the equilibrium, not loaded
!------------------------------------------------------------------------------
subroutine tmap_load_txt(self,io_unit)
class(torflux_map), intent(inout) :: self
integer, intent(in) :: io_unit
end subroutine tmap_load_txt
!------------------------------------------------------------------------------
!> Find interval `lo` such that `x(lo) <= val < x(lo+1)` (bisection, `x` ascending)
!------------------------------------------------------------------------------
pure function tmap_find(x,n,val) result(lo)
integer(4), intent(in) :: n
real(8), intent(in) :: x(n),val
integer(4) :: lo,hi,k
lo=1; hi=n
DO WHILE(hi-lo>1)
  k=(lo+hi)/2
  IF(val<x(k))THEN
    hi=k
  ELSE
    lo=k
  END IF
END DO
end function tmap_find
!------------------------------------------------------------------------------
!> Rebuild uniform-bin lookup index for `map%psis` (call after changing surfaces)
!------------------------------------------------------------------------------
subroutine tmap_index(map)
class(torflux_map), intent(inout) :: map
integer(4) :: j,nbin
IF(ASSOCIATED(map%bin))DEALLOCATE(map%bin)
IF(map%ns<2)RETURN
nbin=4*map%ns
ALLOCATE(map%bin(nbin))
DO j=1,nbin
  map%bin(j)=tmap_find(map%psis,map%ns,REAL(j-1,8)/REAL(nbin,8))
END DO
end subroutine tmap_index
!------------------------------------------------------------------------------
!> Same result as `tmap_find` on `map%psis`, starting from the bin index
!------------------------------------------------------------------------------
function tmap_lookup(map,psihat) result(lo)
class(torflux_map), intent(in) :: map
real(8), intent(in) :: psihat
integer(4) :: lo,nbin
IF(.NOT.ASSOCIATED(map%bin))THEN
  lo=tmap_find(map%psis,map%ns,psihat)
  RETURN
END IF
nbin=SIZE(map%bin)
lo=map%bin(MIN(INT(MIN(MAX(psihat,0.d0),1.d0)*REAL(nbin,8))+1,nbin))
DO WHILE(lo<map%ns-1)
  IF(psihat<map%psis(lo+1))EXIT
  lo=lo+1
END DO
end function tmap_lookup
!------------------------------------------------------------------------------
!> Update the \f$ \hat{\psi} \rightarrow \hat{\Phi} \f$ map from the equilibrium q profile
!!
!! Surfaces are placed at the \f$ \hat{\psi} \f$ of each profile node, with gaps filled to
!! `dphi_max`, plus fixed surfaces near the LCFS and axis. After new profiles or a new state
!! (`settle`, set by @ref torflux_setup) surfaces are re-placed on the updated map until they move
!! less than `place_tol`; otherwise a single pass is made and they settle over the nonlinear iterations. The geometric factor
!! \f$ q/F \f$ is traced once per pass; consistency between F and the map (when F*F' is itself on
!! \f$ \hat{\Phi} \f$) is then obtained by fixed-point iteration without further tracing.
!! Sets `gseq%torflux_ierr`: 0 ok, 1 update failed (previous map kept), 2 no closed surfaces.
!------------------------------------------------------------------------------
subroutine tmap_update(self,gseq)
class(torflux_map), intent(inout) :: self
class(gs_equil), intent(inout) :: gseq
real(8), parameter :: pad = 1.d-3 !< Minimum offset of node surfaces from LCFS and axis
real(8), parameter :: edge_s(2) = [1.d-4,1.d-3] !< Fixed surfaces near the LCFS (edge singularity fit)
real(8), parameter :: axis_s(2) = [1.d-3,3.d-3] !< Fixed surfaces near the axis (axis extrapolation)
real(8), parameter :: dphi_max = 1.d-2 !< Maximum spacing between surfaces in \f$ \hat{\Phi} \f$
real(8), parameter :: fp_tol = 1.d-10 !< Fixed-point tolerance on \f$ \hat{\Phi} \f$
integer(4), parameter :: fp_maxits = 20
real(8), parameter :: place_tol = 1.d-8 !< Surface placement tolerance in \f$ \hat{\psi} \f$
integer(4), parameter :: place_maxits = 5
integer(4) :: i,k,n,nfill,it,ip,npass
real(8) :: dpsi,b,h,err,shift
real(8), allocatable :: phinodes(:),phitarg(:),psis(:),pnew(:),g(:),q(:),r(:),dq(:),cum(:)
gseq%torflux_ierr=0
self%plasma_bounds=gseq%plasma_bounds
CALL torflux_setup(gseq,self%active)
IF((.NOT.self%active).OR.(gseq%plasma_bounds(1)<-1.d98))THEN
  IF(self%ns>0)CALL self%delete()
  RETURN
END IF
dpsi=gseq%plasma_bounds(2)-gseq%plasma_bounds(1)
IF(.NOT.((ABS(dpsi)>0.d0).AND.(ABS(dpsi)<HUGE(dpsi))))THEN
  CALL oft_warn('Toroidal flux map update skipped, no closed flux surfaces')
  gseq%torflux_ierr=2
  RETURN
END IF
!---Target 1-Phi-hat values: [0,1] + profile nodes, gaps filled to dphi_max
ALLOCATE(phinodes(0))
CALL add_nodes(gseq%I)
CALL add_nodes(gseq%P)
CALL add_nodes(gseq%eta)
CALL add_nodes(gseq%I_NI)
CALL add_nodes(gseq%Te)
CALL add_nodes(gseq%Ti)
CALL add_nodes(gseq%ne)
CALL add_nodes(gseq%ni)
CALL add_nodes(gseq%Zeff)
CALL add_nodes(gseq%jphi_fixed)
CALL add_nodes(gseq%jphi_saw)
phinodes=[0.d0,1.d0,MIN(MAX(phinodes,0.d0),1.d0)]
CALL sort_unique(phinodes)
ALLOCATE(phitarg(0))
DO k=1,SIZE(phinodes)-1
  nfill=MAX(CEILING((phinodes(k+1)-phinodes(k))/dphi_max),1)
  phitarg=[phitarg,(phinodes(k)+(phinodes(k+1)-phinodes(k))*REAL(i,8)/REAL(nfill,8), i=0,nfill-1)]
END DO
phitarg=[phitarg,1.d0]
!---Settle surface placement for new profiles or state, otherwise one pass (settles over iterations)
npass=1
IF(self%settle)THEN
  npass=place_maxits
ELSE IF(.NOT.ASSOCIATED(self%nodes))THEN
  npass=place_maxits
ELSE IF(SIZE(self%nodes)/=SIZE(phinodes))THEN
  npass=place_maxits
ELSE IF(ANY(self%nodes/=phinodes))THEN
  npass=place_maxits
END IF
!---Place surfaces at node locations of the current map and re-trace until they settle
shift=-1.d0
DO ip=1,npass
  ALLOCATE(pnew(SIZE(phitarg)))
  DO i=1,SIZE(phitarg)
    pnew(i)=MIN(MAX(self%inv(phitarg(i)),pad),1.d0-pad)
  END DO
  pnew(1)=0.d0; pnew(SIZE(pnew))=1.d0
  pnew=[pnew,edge_s,1.d0-axis_s]
  CALL sort_unique(pnew)
  IF(ip>1)THEN
    IF(SIZE(pnew)==n)THEN
      shift=MAXVAL(ABS(pnew-psis))
      IF(shift<place_tol)EXIT
    END IF
    DEALLOCATE(psis,g,q,r,dq,cum)
  END IF
  CALL MOVE_ALLOC(pnew,psis)
  n=SIZE(psis)
  IF(n<5)THEN
    CALL map_failed('too few surfaces')
    RETURN
  END IF
  !---Trace geometric factor q/F
  ALLOCATE(g(n),q(n),r(n),dq(n),cum(n))
  g=0.d0
  CALL torflux_qgeom(gseq,n-2,psis(2:n-1),g(2:n-1))
  g=ABS(g)
  IF(.NOT.fill_failed())THEN
    CALL map_failed('surface traces failed')
    RETURN
  END IF
  !---Fixed-point iteration for F/map consistency
  IF(ASSOCIATED(gseq%I))CALL torflux_set_bounds(gseq%I,gseq%plasma_bounds)
  dpsi=gseq%plasma_bounds(2)-gseq%plasma_bounds(1)
  DO it=1,fp_maxits
    DO i=2,n-1
      q(i)=ABS(fpol(gseq%plasma_bounds(1)+psis(i)*dpsi))*g(i)
    END DO
    !---Axis: quadratic extrapolation
    q(n)=quad_interp(psis(n-3:n-1),q(n-3:n-1),1.d0)
    !---LCFS: q ~ a + b*ln(s) if diverted, linear extrapolation otherwise
    b=0.d0
    IF(gseq%diverted)THEN
      b=(q(2)-q(3))/LOG(psis(2)/psis(3))
      q(1)=q(2)
    ELSE
      q(1)=MAX(q(2)-(q(3)-q(2))*psis(2)/(psis(3)-psis(2)),0.d0)
    END IF
    !---Integrate r = q - b*ln(s) by Hermite quadrature and b*ln(s) exactly
    r(1)=q(1)
    r(2:n)=q(2:n)-b*LOG(psis(2:n))
    DO i=3,n-2
      dq(i)=fd_deriv(psis(i-1:i+1),r(i-1:i+1),2)
    END DO
    dq(2)=fd_deriv(psis(2:4),r(2:4),1)
    dq(n-1)=fd_deriv(psis(n-2:n),r(n-2:n),2)
    dq(n)=fd_deriv(psis(n-2:n),r(n-2:n),3)
    cum(1)=0.d0
    IF(gseq%diverted)THEN
      cum(2)=r(2)*psis(2)+b*xlogx(psis(2))
    ELSE
      cum(2)=psis(2)*(q(1)+q(2))/2.d0
    END IF
    DO i=2,n-1
      h=psis(i+1)-psis(i)
      cum(i+1)=cum(i)+h*(r(i)+r(i+1))/2.d0+h**2*(dq(i)-dq(i+1))/12.d0 &
        +b*(xlogx(psis(i+1))-xlogx(psis(i)))
    END DO
    IF(.NOT.((cum(n)>0.d0).AND.(cum(n)<HUGE(cum(n))).AND.ALL(q(2:n)>=0.d0)))THEN
      CALL map_failed('non-finite q')
      RETURN
    END IF
    !---Store map
    err=1.d99
    IF(self%ns==n)err=MAXVAL(ABS(cum/cum(n)-self%phihat))
    IF(self%ns/=n)THEN
      CALL self%delete()
      self%ns=n
      ALLOCATE(self%psis(n),self%phihat(n),self%jac(n))
    END IF
    self%Q=cum(n)
    self%psis=psis
    self%phihat=cum/cum(n)
    self%jac=q/cum(n)
    self%stamp=self%stamp+1
    CALL tmap_index(self)
    IF(torflux_coord(gseq%I)==0)EXIT
    CALL build_fint(gseq%I)
    IF(err<fp_tol)EXIT
  END DO
END DO
IF(oft_debug_print(1))THEN
  WRITE(*,'(2A,I4,A,I3,A,I3,A,ES11.3,A,ES11.3)')oft_indent,'Toroidal flux map: ',n,' surfaces, ', &
    MIN(ip,npass),' passes, ',it,' its, shift = ',shift,', Q = ',self%Q
END IF
!---Placement settled for these nodes
self%settle=.FALSE.
IF(ASSOCIATED(self%nodes))DEALLOCATE(self%nodes)
ALLOCATE(self%nodes,SOURCE=phinodes)
IF(ALLOCATED(pnew))DEALLOCATE(pnew)
DEALLOCATE(phinodes,phitarg,psis,g,q,r,dq,cum)
CALL build_all_fint
CONTAINS
!---Rebuild the relabel tables of all profiles against the current map
subroutine build_all_fint
CALL build_fint(gseq%I)
CALL build_fint(gseq%P)
CALL build_fint(gseq%eta)
CALL build_fint(gseq%I_NI)
CALL build_fint(gseq%Te)
CALL build_fint(gseq%Ti)
CALL build_fint(gseq%ne)
CALL build_fint(gseq%ni)
CALL build_fint(gseq%Zeff)
CALL build_fint(gseq%jphi_fixed)
CALL build_fint(gseq%jphi_saw)
end subroutine build_all_fint
!---Rebuild the relabel table of a wrapped profile
subroutine build_fint(F)
CLASS(flux_func), POINTER, INTENT(in) :: F
IF(.NOT.ASSOCIATED(F))RETURN
SELECT TYPE(F)
  TYPE IS(phi_flux_wrapper)
    IF(F%relabel.AND.ASSOCIATED(F%map))CALL phi_wrap_build_fint(F)
END SELECT
end subroutine build_fint
!---Warn and flag a failed update (previous map is kept)
subroutine map_failed(reason)
CHARACTER(LEN=*), INTENT(in) :: reason
CALL oft_warn('Toroidal flux map update failed ('//reason//'), keeping previous map')
gseq%torflux_ierr=1
CALL build_all_fint
end subroutine map_failed
!---Collect the nodes of a wrapped profile (piecewise-linear profiles only)
subroutine add_nodes(F)
CLASS(flux_func), POINTER, INTENT(in) :: F
IF(.NOT.ASSOCIATED(F))RETURN
SELECT TYPE(F)
  TYPE IS(phi_flux_wrapper)
    SELECT TYPE(rep=>F%phi_rep)
      CLASS IS(linterp_flux_func)
        phinodes=[phinodes,rep%x]
    END SELECT
END SELECT
end subroutine add_nodes
!---Poloidal current function F(psi)
function fpol(psi) result(F)
real(8), intent(in) :: psi
real(8) :: F
IF(gseq%mode==0)THEN
  F=gseq%ffp_scale*gseq%I%f(psi)+gseq%I%f_offset
ELSE
  F=SIGN(1.d0,gseq%I%f_offset)*SQRT(MAX(gseq%ffp_scale*gseq%I%f(psi)+gseq%I%f_offset**2,0.d0))
END IF
end function fpol
!---Antiderivative of ln(s)
pure function xlogx(s) result(v)
real(8), intent(in) :: s
real(8) :: v
v=s*LOG(s)-s
end function xlogx
!---Replace failed traces (g=0) by linear interpolation, .FALSE. if too few remain
function fill_failed() result(ok)
logical :: ok
integer(4) :: j,jl,jr
ok=(COUNT(g(2:n-1)>0.d0)>=3)
IF(.NOT.ok)RETURN
DO j=2,n-1
  IF(g(j)>0.d0)CYCLE
  jl=j-1
  DO WHILE(jl>=2)
    IF(g(jl)>0.d0)EXIT
    jl=jl-1
  END DO
  jr=j+1
  DO WHILE(jr<=n-1)
    IF(g(jr)>0.d0)EXIT
    jr=jr+1
  END DO
  IF(jl<2)THEN
    g(j)=g(jr)
  ELSE IF(jr>n-1)THEN
    g(j)=g(jl)
  ELSE
    g(j)=g(jl)+(g(jr)-g(jl))*(psis(j)-psis(jl))/(psis(jr)-psis(jl))
  END IF
END DO
end function fill_failed
end subroutine tmap_update
!------------------------------------------------------------------------------
!> Geometric factor \f$ g = q/F = \frac{1}{2\pi} \oint \frac{dl}{R |\nabla \psi|} \f$ on flux surfaces
!!
!! Contract: `psi_q` in internal normalized flux (0 at LCFS, 1 at axis), strictly inside (0,1);
!! `g` in the units of @ref oft_gs::gs_get_qprof with F=1 (sign ignored by caller); `g(i)=0`
!! flags a failed surface (caller interpolates).
!------------------------------------------------------------------------------
subroutine torflux_qgeom(gseq,nr,psi_q,g)
class(gs_equil), intent(inout) :: gseq !< G-S object
integer(4), intent(in) :: nr !< Number of surfaces
real(8), intent(in) :: psi_q(nr) !< Surface locations
real(8), intent(out) :: g(nr) !< Geometric factor at each surface
IF(torflux_qgeom_backend==1)THEN
  CALL gs_get_qprof(gseq,nr,psi_q,g,geom_only=.TRUE.)
ELSE
  CALL gs_qgeom_cutcell(gseq,nr,psi_q,g)
END IF
end subroutine torflux_qgeom
!------------------------------------------------------------------------------
!> Sort array ascending and remove near-duplicates
!------------------------------------------------------------------------------
subroutine sort_unique(x)
real(8), allocatable, intent(inout) :: x(:)
real(8), parameter :: tol=1.d-8
integer(4) :: i,j,m
real(8) :: tmp
DO i=2,SIZE(x)
  tmp=x(i)
  j=i-1
  DO WHILE(j>=1)
    IF(x(j)<=tmp)EXIT
    x(j+1)=x(j)
    j=j-1
  END DO
  x(j+1)=tmp
END DO
m=MIN(1,SIZE(x))
DO i=2,SIZE(x)
  IF(x(i)-x(m)>tol)THEN
    m=m+1
    x(m)=x(i)
  END IF
END DO
x=x(1:m)
end subroutine sort_unique
!------------------------------------------------------------------------------
!> Quadratic (3-point Lagrange) interpolation/extrapolation
!------------------------------------------------------------------------------
pure function quad_interp(x,y,xi) result(yi)
real(8), intent(in) :: x(3),y(3),xi
real(8) :: yi
yi=y(1)*(xi-x(2))*(xi-x(3))/((x(1)-x(2))*(x(1)-x(3))) &
  + y(2)*(xi-x(1))*(xi-x(3))/((x(2)-x(1))*(x(2)-x(3))) &
  + y(3)*(xi-x(1))*(xi-x(2))/((x(3)-x(1))*(x(3)-x(2)))
end function quad_interp
!------------------------------------------------------------------------------
!> Derivative of 3-point Lagrange interpolant at point `k` of `x`
!------------------------------------------------------------------------------
pure function fd_deriv(x,y,k) result(dy)
real(8), intent(in) :: x(3),y(3)
integer(4), intent(in) :: k
real(8) :: dy,xi
xi=x(k)
dy=y(1)*(2.d0*xi-x(2)-x(3))/((x(1)-x(2))*(x(1)-x(3))) &
  + y(2)*(2.d0*xi-x(1)-x(3))/((x(2)-x(1))*(x(2)-x(3))) &
  + y(3)*(2.d0*xi-x(1)-x(2))/((x(3)-x(1))*(x(3)-x(2)))
end function fd_deriv
end module oft_gs_profiles
