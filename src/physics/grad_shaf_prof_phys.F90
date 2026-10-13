!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file grad_shaf_prof_phys.F90
!
!> Physics-based flux functions for Grad-Sharfranov equilibrium
!!
!! @authors Chris Hansen
!! @date August 2014
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
module grad_shaf_prof_phys
USE oft_io, only: hdf5_write, hdf5_read, hdf5_field_exist
use oft_base
use oft_mesh_type, only: oft_bmesh, bmesh_findcell
USE oft_la_base, ONLY: oft_vector
USE fem_utils, ONLY: bfem_interp
use oft_lag_basis, only: oft_blag_d2eval, oft_blag_geval
USE oft_blag_operators, only: oft_lag_brinterp, oft_lag_bginterp
use oft_gs, only: gs_equil, flux_func, gs_psi2r, oft_indent, &
  oft_increase_indent, oft_decrease_indent, gsinv_interp, gs_prof_interp, &
  gs_get_qprof, gs_ani_press, gs_epsilon, gs_cell_quad
use oft_gs_cutcell, only: gs_fsa_cutcell
use tracing_2d, only: set_tracer, active_tracer, tracinginv_fs
use oft_gs_profiles, only: spline_flux_func, linterp_flux_func, linterp_copy, &
  spline_func_copy, spline_func_delete, phi_flux_wrapper, torflux_qgeom_backend
use spline_mod
USE mhd_utils, ONLY: mu0
USE oft_quadrature, ONLY: oft_quad_type
implicit none
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(gsinv_interp) :: mercierinv_interp
contains
  ! procedure :: setup => minterpinv_setup
  !> Reconstruct the gradient of a Lagrange scalar field
  procedure :: interp => minterpinv_apply
end type mercierinv_interp
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(spline_flux_func) :: mercier_flux_func
  real(8) :: rs = 0.d0
  integer(4) :: ntheta = 128
  TYPE(spline_type) :: funcp
contains
  !> Needs docs
  procedure :: copy => mercier_copy
  !> Needs docs
  procedure :: delete => mercier_delete
  !> Needs docs
  procedure :: update => mercier_update
end type mercier_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(linterp_flux_func) :: jphi_flux_func
  real(8) :: j0 = 0.d0 !< LCFS Jphi value
  real(8) :: norm_last = 1.d0 !< Last Jphi normalization factor (for Ip target)
  real(8), pointer, dimension(:) :: jphi => NULL() !< Jphi(psi) profile values
contains
  !> Needs docs
  procedure :: save_hdf5 => jphi_save_hdf5
  procedure :: save_txt => jphi_save_txt
  !> Needs docs
  procedure :: load_hdf5 => jphi_load_hdf5
  procedure :: load_txt => jphi_load_txt
  !> Needs docs
  procedure :: copy => jphi_copy
  !> Needs docs
  procedure :: delete => jphi_delete
  !> Update F*F' profile from Jphi and current equilibrium state
  procedure :: update => jphi_update
end type jphi_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(gsinv_interp) :: minbinv_interp
  real(8) :: minB = 1.d99
contains
  ! procedure :: setup => minterpinv_setup
  !> Reconstruct the gradient of a Lagrange scalar field
  procedure :: interp => minbinv_apply
end type minbinv_interp
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(spline_flux_func) :: dipole_b0_flux_func
  real(r8) :: psi_pad = 1.d-1
contains
  !> Needs docs
  procedure :: save_hdf5 => dipole_b0_save_hdf5
  procedure :: save_txt => dipole_b0_save_txt
  !> Needs docs
  procedure :: load_hdf5 => dipole_b0_load_hdf5
  procedure :: load_txt => dipole_b0_load_txt
  !> Needs docs
  procedure :: copy => dipole_b0_copy
  !> Needs docs
  procedure :: delete => dipole_b0_delete
  !> Needs docs
  procedure :: update => dipole_b0_update
end type dipole_b0_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(gs_ani_press) :: dipole_ani_press
  REAL(r8) :: a_exp = 0.d0 !< Anisotropy exponent for mirror pressure profiles
  TYPE(oft_lag_brinterp), POINTER :: psi_eval => NULL() !< Needs docs
  TYPE(oft_lag_bginterp), POINTER :: psi_geval => NULL() !< Needs docs
  CLASS(flux_func), POINTER :: B0_prof => NULL() !< Dipole minimum B profile
contains
  !> Needs docs
  procedure :: copy => dipole_ani_copy
  !> Needs docs
  procedure :: setup => dipole_ani_setup
  !> Needs docs
  procedure :: delete => dipole_ani_delete
  !> Evaluate field
  procedure :: interp => dipole_ani_apply
  !>
  procedure :: update => dipole_ani_update
end type dipole_ani_press
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(spline_flux_func) :: mirror_b0_flux_func
  real(r8) :: z_midplane = 0.d0 !< Z location of mirror midplane
contains
  !> Needs docs
  procedure :: save_hdf5 => mirror_b0_save_hdf5
  procedure :: save_txt => mirror_b0_save_txt
  !> Needs docs
  procedure :: load_hdf5 => mirror_b0_load_hdf5
  procedure :: load_txt => mirror_b0_load_txt
  !> Needs docs
  procedure :: copy => mirror_b0_copy
  !> Needs docs
  procedure :: delete => mirror_b0_delete
  !> Needs docs
  procedure :: update => mirror_b0_update
end type mirror_b0_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(gs_ani_press) :: mirror_ani_slosh
  REAL(r8) :: n_exp = 0.d0 !< Anisotropy exponent for mirror pressure profiles
  REAL(r8) :: bturn = 0.d0 !< Turning point for mirror pressure profiles (relative to B_min)
  REAL(r8) :: zthroat = 0.d0 !< Mirror peak field location
  TYPE(oft_lag_brinterp), POINTER :: psi_eval => NULL() !< Needs docs
  TYPE(oft_lag_bginterp), POINTER :: psi_geval => NULL() !< Needs docs
  CLASS(flux_func), POINTER :: B0_prof => NULL() !< Mirror minimum B profile
contains
  !> Needs docs
  procedure :: copy => mirror_slosh_copy
  !> Needs docs
  procedure :: setup => mirror_slosh_setup
  !> Needs docs
  procedure :: delete => mirror_slosh_delete
  !> Evaluate field
  procedure :: interp => mirror_slosh_apply
  !>
  procedure :: update => mirror_slosh_update
end type mirror_ani_slosh
contains
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE create_mercier_ff(func,npsi)
CLASS(flux_func), POINTER, INTENT(out) :: func
INTEGER(4), INTENT(in) :: npsi
INTEGER(4) :: i

ALLOCATE(mercier_flux_func::func)
select type(self=>func)
  type is(mercier_flux_func)
    !---
    self%npsi=npsi
    CALL spline_alloc(self%func,self%npsi,1)
    CALL spline_alloc(self%funcp,self%npsi,1)
    ALLOCATE(self%fun_loc(omp_get_max_threads()))
    DO i=1,omp_get_max_threads()
      CALL spline_alloc(self%fun_loc(i),self%npsi,1)
    END DO
end select
END SUBROUTINE create_mercier_ff
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mercier_copy(self,new)
class(mercier_flux_func), intent(inout) :: self
class(flux_func), pointer, intent(inout) :: new
CALL spline_func_copy(self,new)
SELECT TYPE(new)
  CLASS IS(mercier_flux_func)
    new%plasma_bounds=self%plasma_bounds
    new%f_offset=self%f_offset
    new%rs = self%rs
    new%ntheta = self%ntheta
    CALL spline_alloc(new%funcp,new%npsi,1)
    CALL spline_copy(self%funcp,new%funcp)
END SELECT
end subroutine mercier_copy
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mercier_delete(self)
class(mercier_flux_func), intent(inout) :: self
IF(ASSOCIATED(self%fun_loc))CALL spline_dealloc(self%funcp)
CALL spline_func_delete(self)
end subroutine mercier_delete
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mercier_update(self,gseq)
class(mercier_flux_func), intent(inout) :: self
class(gs_equil), intent(inout) :: gseq
type(mercierinv_interp), target :: field
type(oft_lag_brinterp) :: psi_int
real(8) :: I,Ip,q,qp,vp,vpp,s,a,b,pp,gop(3,3),psi_surf(1),vol,rmax
real(8) :: raxis,zaxis,f(3),pt(3),x1,x2
real(8), pointer, dimension(:) :: v
real(8), parameter :: tol=1.d-10
integer(4) :: j,k,cell
!
IF(oft_debug_print(2))THEN
  WRITE(*,'(2A)')oft_indent,'Updating Mercier Pressure:'
  CALL oft_increase_indent
END IF
psi_int%u=>gseq%psi
CALL psi_int%setup(gseq%device%fe_rep)
raxis=gseq%o_point(1)
zaxis=gseq%o_point(2)
IF(oft_debug_print(2))WRITE(*,'(2A,3ES11.3)')oft_indent,'Axis Position = ',raxis,zaxis
!---
x1=0.d0; x2=1.d0
IF(gseq%plasma_bounds(1)>-1.d98)THEN
  x1=gseq%plasma_bounds(1); x2=gseq%plasma_bounds(2)
END IF
!---Find Rmax along Zaxis
rmax=raxis
cell=0
DO j=1,100
  pt=[(gseq%device%rmax-raxis)*j/REAL(100,8)+raxis,zaxis,0.d0]
  CALL bmesh_findcell(gseq%device%mesh,cell,pt,f)
  IF( (MAXVAL(f)>1.d0+tol) .OR. (MINVAL(f)<-tol) )EXIT
  CALL psi_int%interp(cell,f,gop,psi_surf)
  IF( psi_surf(1) < x1)EXIT
  rmax=pt(1)
END DO
IF(oft_debug_print(2))WRITE(*,'(2A,ES11.3)')oft_indent,'Rmax = ',rmax
!---Trace
call set_tracer(1)
!$omp parallel private(field,gop,vol,psi_surf,I,Ip,v,q,qp,vp,vpp,s,a,b,pp,pt)
pt=[(.9d0*rmax+.1d0*raxis),zaxis,0.d0]
field%u=>gseq%psi
CALL field%setup(gseq%device%fe_rep)
active_tracer%neq=8
active_tracer%B=>field
active_tracer%maxsteps=8e4
active_tracer%tol=1.d-9
active_tracer%raxis=raxis
active_tracer%zaxis=zaxis
active_tracer%inv=.TRUE.
!$omp do schedule(dynamic,1)
do j=1,self%npsi-1
    !------------------------------------------------------------------------------
    ! Trace contour
    !------------------------------------------------------------------------------
    !psi_surf(1)=(x2-x1)*(1.d0-j/REAL(self%npsi,4))**2
    psi_surf(1)=(x2-x1)*(1.d0-j/REAL(self%npsi,4))
    psi_surf(1)=x2 - psi_surf(1)
    !$omp critical
    CALL gs_psi2r(gseq,psi_surf(1),pt)
    !$omp end critical
    call tracinginv_fs(gseq%device%mesh,pt(1:2))
    !---Exit if trace fails
    if(active_tracer%status/=1)THEN
      WRITE(*,*)'Tracer Error:',psi_surf(1),pt,active_tracer%y,active_tracer%status
      call oft_abort('Trace did not complete.','mercier_update',__FILE__)
    end if
    !------------------------------------------------------------------------------
    ! Compute Mercier Profiles
    !------------------------------------------------------------------------------
    !---Compute poloidal flux
    pt(1:2)=active_tracer%y
    call bmesh_findcell(gseq%device%mesh,active_tracer%cell,pt,active_tracer%f)
    call gseq%device%mesh%jacobian(active_tracer%cell,active_tracer%f,gop,vol)
    call psi_int%interp(active_tracer%cell,active_tracer%f,gop,psi_surf)
    !---Get flux variables
    I=gseq%ffp_scale*gseq%I%f(psi_surf(1))+gseq%I%f_offset
    Ip=gseq%ffp_scale*gseq%I%fp(psi_surf(1))
    v=>active_tracer%v
    !---Compute profile variables
    q=I*v(3)/(2*pi) ! Safety Factor (q)
    qp=(Ip*v(3)-I*v(8))/(2*pi) ! q-Shear (q')
    vp=-2*pi*v(2) ! First derivative of FS volume (V')
    vpp=2*pi*v(7) ! Second derivative of FS volume (V')
    !---Compute Mercier pressure
    s=Ip*v(3)-I*v(8) ! Shear term
    a=(I**2)*(v(5)*v(6)-v(4)**2)+v(3)*v(6) ! Linear term
    b=s*I*v(4)+v(7)*((I**2)*v(5)+v(3)) ! Quadratic term
    pp=.5d0*(sqrt(b**2+a*(s**2))-b)/a ! Mercier marginally stable pressure
    !---
    self%funcp%xs(j)=psi_surf(1)
    self%funcp%fs(j,1)=pp
    !WRITE(*,*)psi_surf(1),pp
end do
CALL field%delete
!$omp end parallel
self%funcp%xs(0)=x1
self%funcp%fs(0,1)=0.d0
!self%funcp%fs(0,1)=self%funcp%fs(1,1)
self%funcp%xs(self%npsi)=x2+.05d0*x2
self%funcp%fs(self%npsi,1)=0.d0
!---Setup Spline
CALL spline_fit(self%funcp,"extrap")
CALL spline_int(self%funcp)
!---
self%func%xs=self%funcp%xs
self%func%fs(:,1)=self%funcp%fsi(:,1)
CALL spline_fit(self%func,"extrap")
DO k=1,omp_get_max_threads()
  CALL spline_copy(self%func,self%fun_loc(k))
END DO
!
self%xmin=self%funcp%xs(0)
self%xmax=self%funcp%xs(self%npsi)
IF(oft_debug_print(2))CALL oft_decrease_indent
CALL psi_int%delete
end subroutine mercier_update
! !------------------------------------------------------------------------------
! !> Needs Docs
! !------------------------------------------------------------------------------
! subroutine minterpinv_setup(self,mesh)
! class(mercierinv_interp), intent(inout) :: self
! class(oft_bmesh), target, intent(inout) :: mesh
! NULLIFY(self%uvals)
! CALL self%u%get_local(self%uvals)
! self%mesh=>mesh
! end subroutine minterpinv_setup
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine minterpinv_apply(self,cell,f,gop,val)
class(mercierinv_interp), intent(inout) :: self
integer(4), intent(in) :: cell
real(8), intent(in) :: f(:)
real(8), intent(in) :: gop(3,3)
real(8), intent(out) :: val(:)
integer(4), allocatable :: j(:)
integer(4) :: jc
real(8) :: rop(3),d2op(6),pt(3),grad(3),d2(3),d2_tmp(6),tmp
real(8) :: g2op(6,6),s,c
real(8) :: K(6,3)
!---Get dofs
allocate(j(self%lag_rep%nce))
call self%lag_rep%ncdofs(cell,j)
!---Reconstruct gradient
call self%mesh%hessian(cell,f,g2op,K)
grad=0.d0
d2_tmp=0.d0
do jc=1,self%lag_rep%nce
  call oft_blag_geval(self%lag_rep,cell,jc,f,rop,gop)
  call oft_blag_d2eval(self%lag_rep,cell,jc,f,d2op,g2op)
  grad=grad+self%uvals(j(jc))*rop
  d2_tmp=d2_tmp+self%uvals(j(jc))*(d2op-MATMUL(g2op,MATMUL(K,rop)))
end do
d2=d2_tmp([1,2,4]) ! Map from 3D to 2D
!---Get radial position
pt=self%mesh%log2phys(cell,f)
!---
s=SIN(self%t)
c=COS(self%t)
!---
val(1)=(self%rho*(grad(1)*s-grad(2)*c))/(grad(1)*c+grad(2)*s)
val(2)=pt(1)*SQRT((self%rho**2+val(1)**2)/SUM(grad**2))
!---
val(3)=1.d0/((pt(1)+self%eps)**2)
!---
val(4)=1.d0/SUM(grad**2)
val(5)=val(4)*val(3)
val(6)=val(4)/(val(3)+self%eps)
!---
tmp=(grad(1)**2-grad(2)**2)*(d2(1)-d2(3))+4.d0*grad(1)*grad(2)*d2(2)
val(7)=(tmp*val(4)-grad(1)/(pt(1)+self%eps))*val(4)
val(8)=(tmp*val(4)+grad(1)/(pt(1)+self%eps))*val(5)
val(3:8)=val(3:8)*val(2)
deallocate(j)
end subroutine minterpinv_apply
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine jphi_save_hdf5(self,filename,path)
class(jphi_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
IF(.NOT.hdf5_field_exist(filename,path//'/TYPE'))CALL hdf5_write('jphi-linterp',filename,path//'/TYPE')
CALL hdf5_write(self%npsi,filename,path//'/NPSI')
CALL hdf5_write(self%x,filename,path//'/XVALS')
CALL hdf5_write(self%jphi,filename,path//'/YVALS')
CALL hdf5_write(self%j0,filename,path//'/J0')
end subroutine jphi_save_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine jphi_save_txt(self,io_unit)
class(jphi_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
WRITE(io_unit,*)'jphi-linterp'
WRITE(io_unit,*)self%npsi,self%j0
WRITE(io_unit,*)self%x
WRITE(io_unit,*)self%jphi
end subroutine jphi_save_txt
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine jphi_load_hdf5(self,filename,path,success)
class(jphi_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
logical, intent(out) :: success
integer(i4) :: npsi
real(r8) :: J0
real(r8), allocatable :: xvals(:),yvals(:)
CALL hdf5_read(npsi,filename,path//'/NPSI',success=success)
IF(.NOT.success)RETURN
ALLOCATE(xvals(npsi),yvals(npsi))
CALL hdf5_read(xvals,filename,path//'/XVALS',success=success)
IF(.NOT.success)RETURN
CALL hdf5_read(yvals,filename,path//'/YVALS',success=success)
IF(.NOT.success)RETURN
CALL hdf5_read(J0,filename,path//'/J0',success=success)
IF(.NOT.success)RETURN
CALL create_jphi_ff(self,npsi,xvals,yvals,J0)
DEALLOCATE(xvals,yvals)
end subroutine jphi_load_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine jphi_load_txt(self,io_unit)
class(jphi_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
integer(i4) :: npsi
real(r8) :: J0
real(r8), allocatable :: xvals(:),yvals(:)
READ(io_unit,*)npsi,J0
ALLOCATE(xvals(npsi),yvals(npsi))
READ(io_unit,*)xvals
READ(io_unit,*)yvals
CALL create_jphi_ff(self,npsi,xvals,yvals,J0)
DEALLOCATE(xvals,yvals)
end subroutine jphi_load_txt
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_jphi_ff(func,npsi,psivals,yvals,y0)
CLASS(flux_func), INTENT(inout) :: func
INTEGER(4), INTENT(in) :: npsi
REAL(8), INTENT(in) :: psivals(npsi)
REAL(8), INTENT(in) :: yvals(npsi)
REAL(8), INTENT(in) :: y0
INTEGER(4) :: i,ierr
! IF(.NOT.ASSOCIATED(func))ALLOCATE(jphi_flux_func::func)
SELECT TYPE(self=>func)
  TYPE IS(jphi_flux_func)
  !---
  self%npsi=npsi
  self%ndofs=self%npsi
  !---
  ALLOCATE(self%x(self%npsi))
  ALLOCATE(self%yp(self%npsi))
  ALLOCATE(self%y(self%npsi))
  ALLOCATE(self%jphi(self%npsi))
  !---
  self%j0=y0
  self%y0=0.d0
  DO i=1,self%npsi
    self%x(i) = psivals(i)
    self%jphi(i) = yvals(i)
    self%yp(i) = psivals(i) ! Dummy initialization
  END DO
  self%yp = self%yp/(SUM(ABS(self%yp))/REAL(self%npsi,8)) ! Consistent (hopefully) normalization
  ierr=self%set_cofs(self%yp)
  IF(oft_debug_print(1))WRITE(*,*)'Jphi linear interpolator Created',self%ndofs,self%x,self%j0
class default
  CALL oft_abort('Invalid flux function type in create_jphi_ff','create_jphi_ff',__FILE__)
END SELECT

END SUBROUTINE create_jphi_ff
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine jphi_copy(self,new)
class(jphi_flux_func), intent(inout) :: self
class(flux_func), pointer, intent(inout) :: new
CALL linterp_copy(self,new)
SELECT TYPE(new)
  CLASS IS(jphi_flux_func)
    new%plasma_bounds=self%plasma_bounds
    new%f_offset=self%f_offset
    new%j0 = self%j0
    new%norm_last = self%norm_last
    ALLOCATE(new%jphi, SOURCE=self%jphi)
END SELECT
end subroutine jphi_copy
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine jphi_delete(self)
class(jphi_flux_func), intent(inout) :: self
self%j0=0.d0
IF(ASSOCIATED(self%jphi))DEALLOCATE(self%jphi)
IF(ASSOCIATED(self%x))DEALLOCATE(self%x)
IF(ASSOCIATED(self%yp))DEALLOCATE(self%yp)
IF(ASSOCIATED(self%y))DEALLOCATE(self%y)
end subroutine jphi_delete
!------------------------------------------------------------------------------
!> Update F*F' profile from Jphi, P', and current equilibrium
!------------------------------------------------------------------------------
subroutine jphi_update(self,gseq)
class(jphi_flux_func), intent(inout) :: self
class(gs_equil), intent(inout) :: gseq
INTEGER(i4) :: i
REAL(r8) :: jphi_norm,pscale,pprime,ip0,ip1
REAL(r8), ALLOCATABLE :: jtor(:),xpsi(:),ravgs(:,:),xn(:),jn(:)
self%plasma_bounds=gseq%plasma_bounds
CALL jphi_psi_nodes(self,xpsi)
IF(gseq%mode/=1)CALL oft_abort("Jphi profile requires (F^2)' formulation","jphi_update",__FILE__)
IF(gseq%Ip_target<0.d0)CALL oft_abort("Jphi profile requires Ip target","jphi_update",__FILE__)
IF(gseq%pax_target<0.d0)CALL oft_abort("Jphi profile requires Pax target","jphi_update",__FILE__)
!---Flux surface geometry for Jphi -> F*F' mapping at the LCFS (row 0) and nodes
ALLOCATE(ravgs(0:self%npsi,3))
CALL gs_ravgs(gseq, self%npsi+1, [0.d0, xpsi], ravgs)
!---Get pressure profile
CALL gseq%P%update(gseq) ! Make sure pressure profile is up to date with EQ
IF(ASSOCIATED(gseq%P_ani))CALL oft_abort('Jphi profiles do not support anistopic pressure','jphi_update',__FILE__)
pscale=gseq%P%f(gseq%plasma_bounds(2))
pscale=gseq%pax_target/pscale
!---Update jphi normalization to match Ip target: the exact I_p (A9c) is affine in the jphi scale,
!   I_p(norm) = ip0 + norm*(ip1 - ip0); past the last node jphi is held, as F*F' is
ALLOCATE(xn(0:self%npsi),jn(0:self%npsi),jtor(0:self%npsi))
xn=[0.d0,xpsi]; jn=[self%j0,self%jphi]
CALL eval_dIdA(gseq, xn, self%npsi+1, ravgs, jn, pscale, jtor)
CALL gs_flux_int(gseq,xn,jtor,self%npsi+1,ip1,extrap=1)
jn=0.d0
CALL eval_dIdA(gseq, xn, self%npsi+1, ravgs, jn, pscale, jtor)
CALL gs_flux_int(gseq,xn,jtor,self%npsi+1,ip0,extrap=1)
DEALLOCATE(xn,jn,jtor)
IF(ABS(ip1-ip0)>0.d0)THEN
  jphi_norm=(ABS(gseq%Ip_target)-ip0)/(ip1-ip0)
ELSE
  CALL oft_warn("jphi_update: jphi profile carries no (finite) net current, keeping the last normalization")
  jphi_norm=self%norm_last
END IF
self%norm_last=jphi_norm
!---Compute updated F*F' profile ! 2.0*(jtor -  R_avg * (-pprime)) * (mu0 / one_over_R_avg)
pprime=gseq%P%fp(gseq%plasma_bounds(1))
self%y0 = 2.d0*(self%j0*jphi_norm - ravgs(0,1)*pprime*pscale)/ravgs(0,2)
DO i=1,self%npsi
  pprime=gseq%P%fp(xpsi(i)*(gseq%plasma_bounds(2)-gseq%plasma_bounds(1))+gseq%plasma_bounds(1))
  self%yp(i) = 2.d0*(self%jphi(i)*jphi_norm - ravgs(i,1)*pprime*pscale)/ravgs(i,2)
END DO
! Disable Ip matching and fix F*F' scale (matching is done here instead)
gseq%skip_targets=.TRUE.
! IF(gseq%Ip_target>0.d0)gseq%Ip_target=-gseq%Ip_target
gseq%ffp_scale=1.d0
gseq%p_scale=pscale
!---Clean up
DEALLOCATE(ravgs)
i=self%set_cofs(self%yp)
end subroutine jphi_update
!------------------------------------------------------------------------------
!> Node locations in normalized poloidal flux (`x` mapped from toroidal flux if wrapped)
!------------------------------------------------------------------------------
SUBROUTINE jphi_psi_nodes(self,xpsi)
CLASS(jphi_flux_func), INTENT(inout) :: self
REAL(r8), ALLOCATABLE, INTENT(out) :: xpsi(:) !< Node locations in \f$ 1-\hat{\psi} \f$
INTEGER(i4) :: i
ALLOCATE(xpsi(self%npsi))
xpsi=self%x
IF(.NOT.ASSOCIATED(self%phi_wrapper))RETURN
SELECT TYPE(wrap=>self%phi_wrapper)
  TYPE IS(phi_flux_wrapper)
    IF(.NOT.wrap%relabel)CALL oft_abort('Jphi profiles on toroidal flux must use "phi_n_relabel"', &
      'jphi_psi_nodes',__FILE__)
    IF(.NOT.ASSOCIATED(wrap%map))RETURN
    DO i=1,self%npsi
      xpsi(i)=wrap%map%inv(self%x(i))
    END DO
END SELECT
END SUBROUTINE jphi_psi_nodes
!---------------------------------------------------------------------------------
!> Area integral over the plasma (inside the LCFS) of a flux function given on
!! normalized psi nodes (0 at the LCFS, 1 at the axis)
!---------------------------------------------------------------------------------
SUBROUTINE gs_flux_int(self,psi_tmp,field_tmp,nvals,result,extrap)
class(gs_equil), INTENT(inout) :: self !< Pointer to TokaMaker object
REAL(r8), INTENT(in) :: psi_tmp(:) !< Normalised psi nodes (OFT convention)
REAL(r8), INTENT(in) :: field_tmp(:) !< Flux function on @p psi_tmp
INTEGER(i4), INTENT(in) :: nvals !< Number of nodes
REAL(r8), INTENT(out) :: result !< Integral over the plasma cross-section
INTEGER(i4), OPTIONAL, INTENT(in) :: extrap !< Beyond the last node: 0 skip (default), 1 hold its value (see linterp)
INTEGER(i4) :: i,m,extrap_loc
REAL(r8) :: area,psitmp(1),sgop(3,3)
TYPE(gs_prof_interp) :: prof_interp_obj
TYPE(oft_quad_type), POINTER :: quad
! DEBUG_STACK_PUSH
!---Setup
extrap_loc=0
IF(PRESENT(extrap))extrap_loc=extrap
result=0.d0
prof_interp_obj%mode=4
CALL prof_interp_obj%setup(self)
!$omp parallel do private(m,psitmp,sgop,area,quad) reduction(+:result) schedule(dynamic,16)
do i=1,self%device%mesh%nc
  IF(self%device%mesh%reg(i)/=1)CYCLE
  CALL gs_cell_quad(self,i,prof_interp_obj%psi_eval%vals,quad)
  !---Loop over quadrature points
  do m=1,quad%np
    call self%device%mesh%jacobian(i,quad%pts(:,m),sgop,area)
    call prof_interp_obj%interp(i,quad%pts(:,m),sgop,psitmp)
    IF(psitmp(1)<=0.d0)CYCLE ! Outside the plasma (mode 4 gives 0 outside the LCFS or gs_test_bounds)
    psitmp(1)=linterp(psi_tmp,field_tmp,nvals,psitmp(1),extrap_loc)
    IF(psitmp(1)>-1.d98)result = result + psitmp(1)*area*quad%wts(m)
  end do
end do
!---Global reduction and cleanup
result=oft_mpi_sum(result)
CALL prof_interp_obj%delete()
! DEBUG_STACK_POP
END SUBROUTINE gs_flux_int
!---------------------------------------------------------------------------------
!> Area integrals of several flux functions inside each node's surface, in one quadrature
!! pass (same quadrature as @ref gs_flux_int). cum(1,:) is the full plasma integral.
!---------------------------------------------------------------------------------
SUBROUTINE gs_flux_cumint(self,psi_tmp,fields,nvals,nf,cum,extrap)
CLASS(gs_equil), INTENT(inout) :: self !< TokaMaker object
INTEGER(i4), INTENT(in) :: nvals !< Number of nodes
INTEGER(i4), INTENT(in) :: nf !< Number of fields
REAL(r8), INTENT(in) :: psi_tmp(nvals) !< Nodes in normalized psi (0 at LCFS, 1 at axis), ascending
REAL(r8), INTENT(in) :: fields(nvals,nf) !< Fields on the nodes
REAL(r8), INTENT(out) :: cum(nvals,nf) !< Integral over psi >= psi_tmp(k) of each field
INTEGER(i4), OPTIONAL, INTENT(in) :: extrap !< Beyond the last node: 0 skip (default), 1 hold its value (see linterp)
INTEGER(i4), PARAMETER :: nblk_cells=64 !< Cells per block (fixed-order sums, independent of the thread count)
INTEGER(i4) :: i,m,k,ib,nblk,inds(2),extrap_loc
REAL(r8) :: area,psitmp(1),sgop(3,3),facs(2)
REAL(r8), ALLOCATABLE :: bins(:,:),bblk(:,:,:)
TYPE(gs_prof_interp) :: prof_interp_obj
TYPE(oft_quad_type), POINTER :: quad
extrap_loc=0
IF(PRESENT(extrap))extrap_loc=extrap
nblk=(self%device%mesh%nc+nblk_cells-1)/nblk_cells
ALLOCATE(bins(nvals,nf),bblk(nvals,nf,nblk))
bblk=0.d0
prof_interp_obj%mode=4
CALL prof_interp_obj%setup(self)
!$omp parallel do private(i,m,psitmp,sgop,area,inds,facs,quad) schedule(dynamic,1)
DO ib=1,nblk
  DO i=(ib-1)*nblk_cells+1,MIN(ib*nblk_cells,self%device%mesh%nc)
    IF(self%device%mesh%reg(i)/=1)CYCLE
    CALL gs_cell_quad(self,i,prof_interp_obj%psi_eval%vals,quad)
    DO m=1,quad%np
      CALL self%device%mesh%jacobian(i,quad%pts(:,m),sgop,area)
      CALL prof_interp_obj%interp(i,quad%pts(:,m),sgop,psitmp)
      IF(psitmp(1)<=0.d0)CYCLE ! outside the plasma
      CALL linterp_facs(psi_tmp,nvals,psitmp(1),inds,facs,extrap_loc)
      IF(inds(1)<=0)CYCLE
      area=area*quad%wts(m)
      ! Lower node inds(2) is the innermost surface enclosing this point
      bblk(inds(2),:,ib)=bblk(inds(2),:,ib)+(fields(inds(1),:)*facs(1)+fields(inds(2),:)*facs(2))*area
    END DO
  END DO
END DO
bins=0.d0
DO ib=1,nblk
  bins=bins+bblk(:,:,ib)
END DO
DEALLOCATE(bblk)
DO k=1,nf
  bins(:,k)=oft_mpi_sum(bins(:,k),nvals)
END DO
cum(nvals,:)=bins(nvals,:)
DO k=nvals-1,1,-1
  cum(k,:)=cum(k+1,:)+bins(k,:)
END DO
CALL prof_interp_obj%delete()
DEALLOCATE(bins)
END SUBROUTINE gs_flux_cumint
!------------------------------------------------------------------------------
!> Flux surface averages <R>, <1/R>, <1/R^2> on normalized psi nodes (0 at the LCFS, 1 at the axis)
!!
!! Cut-cell quadrature (@ref oft_gs_cutcell::gs_fsa_cutcell), tracing only surfaces it cannot
!! close, or all when `torflux_qgeom_backend=1` (testing). The axis takes its exact limit; nodes
!! at or outside the LCFS are traced on it if limited, or `psi_pad` inside it if diverted. Nodes
!! that fail both methods take the averages of the nearest node that did not.
!------------------------------------------------------------------------------
SUBROUTINE gs_ravgs(gseq, n, psi_vals, ravgs)
CLASS(gs_equil), INTENT(inout) :: gseq !< G-S object
INTEGER(i4), INTENT(in) :: n !< Number of nodes
REAL(r8), INTENT(in) :: psi_vals(n) !< Node locations
REAL(r8), INTENT(out) :: ravgs(n,3) !< <R>, <1/R>, <1/R^2> at each node
REAL(r8), PARAMETER :: psi_pad = 1.d-3
REAL(r8) :: psi_q(n), qf(n), avgs(n,4)
REAL(r8), ALLOCATABLE :: rtr(:,:), qtr(:)
LOGICAL :: ok(n)
INTEGER(i4) :: k, m
psi_q = MERGE(psi_vals, MERGE(psi_pad, 0.d0, gseq%diverted), psi_vals>0.d0)
ravgs = 0.d0
ok = .FALSE.
IF(torflux_qgeom_backend/=1)THEN
  CALL gs_fsa_cutcell(gseq, n, psi_q, qf, avgs, ok)
  ravgs = avgs(:,1:3)
END IF
WHERE(psi_q>=1.d0)
  ravgs(:,1) = gseq%o_point(1)
  ravgs(:,2) = 1.d0/gseq%o_point(1)
  ravgs(:,3) = 1.d0/gseq%o_point(1)**2
  ok = .TRUE.
END WHERE
m = COUNT(.NOT.ok)
IF(m==0)RETURN
ALLOCATE(rtr(m,4), qtr(m))
rtr = 0.d0
CALL gs_get_qprof(gseq, m, PACK(psi_q, .NOT.ok), qtr, ravgs=rtr)
DO k=1,3
  ravgs(:,k) = UNPACK(rtr(:,k), .NOT.ok, ravgs(:,k))
END DO
DEALLOCATE(rtr, qtr)
!---Nodes that failed both methods (zero averages) take the averages of the nearest node that did not
ok = (ravgs(:,2) > 0.d0)
IF(ALL(ok).OR.(.NOT.ANY(ok)))RETURN ! None: no equilibrium yet, the caller sees zeros
DO k=1,n
  IF(ok(k))CYCLE
  DO m=1,n
    IF(k-m>=1)THEN
      IF(ok(k-m))THEN
        ravgs(k,:)=ravgs(k-m,:); EXIT
      END IF
    END IF
    IF(k+m<=n)THEN
      IF(ok(k+m))THEN
        ravgs(k,:)=ravgs(k+m,:); EXIT
      END IF
    END IF
  END DO
END DO
END SUBROUTINE gs_ravgs
!------------------------------------------------------------------------------
!> Current per unit cross-section area, dI/dA = <j_phi/R>/<1/R>, from TokaMaker jphi = <j_phi> on the psi_N grid psi_vals
!! (OFT convention), eq. A5 of doc_tokamaker_current_conventions; gs_flux_int of the result is
!! exactly I_p (eq. A9c). jphi and P'*pscale must share units.
!------------------------------------------------------------------------------
SUBROUTINE eval_dIdA(gseq, psi_vals, n, ravgs, jphi, pscale, jtor)
CLASS(gs_equil), INTENT(inout) :: gseq !< Equilibrium
INTEGER(i4), INTENT(in) :: n !< Number of nodes
REAL(r8), INTENT(in) :: psi_vals(n) !< Normalised psi nodes (OFT convention)
REAL(r8), INTENT(in) :: ravgs(n,3) !< <R>, <1/R>, <1/R^2> at psi_vals (@ref gs_ravgs)
REAL(r8), INTENT(in) :: jphi(n) !< TokaMaker jphi = <j_phi> at psi_vals
REAL(r8), INTENT(in) :: pscale !< Pressure scale applied to P'
REAL(r8), INTENT(out) :: jtor(n) !< <j_phi/R>/<1/R>
INTEGER(i4) :: i
REAL(r8) :: pprime
DO i = 1, n
  pprime = pscale*gseq%P%fp(psi_vals(i)*(gseq%plasma_bounds(2)-gseq%plasma_bounds(1))+gseq%plasma_bounds(1))
  jtor(i) = (pprime + ravgs(i,3)*(jphi(i) - ravgs(i,1)*pprime)/ravgs(i,2))/ravgs(i,2)
END DO
END SUBROUTINE eval_dIdA
!---------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE create_dipole_b0_prof(func,npsi)
CLASS(flux_func), INTENT(inout) :: func
INTEGER(4), INTENT(in) :: npsi
INTEGER(4) :: i
select type(self=>func)
  type is(dipole_b0_flux_func)
    !---
    self%npsi=npsi
    CALL spline_alloc(self%func,self%npsi,1)
    ALLOCATE(self%fun_loc(omp_get_max_threads()))
    DO i=1,omp_get_max_threads()
      CALL spline_alloc(self%fun_loc(i),self%npsi,1)
    END DO
class default
  CALL oft_abort('Invalid flux function type in create_dipole_b0_prof','create_dipole_b0_prof',__FILE__)
end select
END SUBROUTINE create_dipole_b0_prof
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine dipole_b0_save_hdf5(self,filename,path)
class(dipole_b0_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
IF(.NOT.hdf5_field_exist(filename,path//'/TYPE'))CALL hdf5_write('dipole_b0',filename,path//'/TYPE')
CALL hdf5_write(self%npsi,filename,path//'/NPSI')
CALL hdf5_write(self%psi_pad,filename,path//'/PSI_PAD')
end subroutine dipole_b0_save_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine dipole_b0_save_txt(self,io_unit)
class(dipole_b0_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
WRITE(io_unit,*)'dipole_b0'
WRITE(io_unit,*)self%npsi,self%psi_pad
end subroutine dipole_b0_save_txt
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine dipole_b0_load_hdf5(self,filename,path,success)
class(dipole_b0_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
logical, intent(out) :: success
integer(i4) :: npsi
real(r8), allocatable :: xvals(:),yvals(:)
CALL hdf5_read(npsi,filename,path//'/NPSI',success=success)
IF(.NOT.success)RETURN
CALL hdf5_read(self%psi_pad,filename,path//'/PSI_PAD',success=success)
IF(.NOT.success)RETURN
CALL create_dipole_b0_prof(self,npsi)
end subroutine dipole_b0_load_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine dipole_b0_load_txt(self,io_unit)
class(dipole_b0_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
integer(i4) :: npsi
READ(io_unit,*)npsi,self%psi_pad
CALL create_dipole_b0_prof(self,npsi)
end subroutine dipole_b0_load_txt
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine dipole_b0_copy(self,new)
class(dipole_b0_flux_func), intent(inout) :: self
class(flux_func), pointer, intent(inout) :: new
CALL spline_func_copy(self,new)
SELECT TYPE(new)
  TYPE IS(dipole_b0_flux_func)
    new%plasma_bounds=self%plasma_bounds
    new%f_offset=self%f_offset
    new%psi_pad=self%psi_pad
END SELECT
end subroutine dipole_b0_copy
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine dipole_b0_delete(self)
class(dipole_b0_flux_func), intent(inout) :: self
integer(i4) :: i
CALL spline_func_delete(self)
end subroutine dipole_b0_delete
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine dipole_b0_update(self,gseq)
class(dipole_b0_flux_func), intent(inout) :: self
class(gs_equil), intent(inout) :: gseq
type(minbinv_interp), pointer :: field
type(oft_lag_brinterp) :: psi_int
real(8) :: I,Ip,q,qp,vp,vpp,s,a,b,pp,gop(3,3),psi_surf(1),vol,rmax
real(8) :: raxis,zaxis,f(3),pt(3),x1,x2,xr
real(8), pointer, dimension(:) :: v
real(8), parameter :: tol=1.d-10
integer(4) :: j,k,cell
!
IF(oft_debug_print(2))THEN
  WRITE(*,'(2A)')oft_indent,'Updating Dipole B0 profile:'
  CALL oft_increase_indent
END IF
psi_int%u=>gseq%psi
CALL psi_int%setup(gseq%device%fe_rep)
raxis=gseq%o_point(1)
zaxis=gseq%o_point(2)
IF(oft_debug_print(2))WRITE(*,'(2A,3ES11.3)')oft_indent,'Axis Position = ',raxis,zaxis
!---
x1=0.d0; x2=1.d0
IF(gseq%plasma_bounds(1)>-1.d98)THEN
  x1=gseq%plasma_bounds(1); x2=gseq%plasma_bounds(2)
END IF
xr = (x2-x1)
IF(.TRUE.)THEN
  x1 = x1 + xr*self%psi_pad
  xr = (x2-x1)
END IF
!---Find Rmax along Zaxis
rmax=raxis
cell=0
DO j=1,100
  pt=[raxis*j/REAL(100,8),0.d0,0.d0]
  CALL bmesh_findcell(gseq%device%mesh,cell,pt,f)
  IF( (MAXVAL(f)>1.d0+tol) .OR. (MINVAL(f)<-tol) )EXIT
  CALL psi_int%interp(cell,f,gop,psi_surf)
  IF( psi_surf(1) > x1)EXIT
  rmax=pt(1)
END DO
IF(oft_debug_print(2))WRITE(*,'(2A,ES11.3)')oft_indent,'Rmin = ',rmax
!---Trace
call set_tracer(1)
!$omp parallel private(field,gop,vol,psi_surf,I,Ip,v,q,qp,vp,vpp,s,a,b,pp,pt)
pt=[(.9d0*rmax+.1d0*raxis),zaxis,0.d0]
ALLOCATE(field)
field%u=>gseq%psi
CALL field%setup(gseq%device%fe_rep)
active_tracer%neq=2
active_tracer%B=>field
active_tracer%maxsteps=8e4
active_tracer%tol=1.d-9
active_tracer%raxis=raxis
active_tracer%zaxis=zaxis
active_tracer%inv=.TRUE.
!$omp do schedule(dynamic,1)
do j=1,self%npsi+1
  psi_surf(1)=(x2-x1)*(1.d0-j/REAL(self%npsi+2,4))
  psi_surf(1)=x2 - psi_surf(1)
  CALL gs_psi2r(gseq,psi_surf(1),pt,psi_int)
  field%minB=1.d99
  call tracinginv_fs(gseq%device%mesh,pt(1:2))
  !---Exit if trace fails
  if(active_tracer%status/=1)THEN
    WRITE(*,*)'Tracer Error:',psi_surf(1),pt,active_tracer%y,active_tracer%status
    call oft_abort('Trace did not complete.','dipole_b0_update',__FILE__)
  end if
  !---Get surface minB
  self%func%xs(j-1)=psi_surf(1)
  self%func%fs(j-1,1)=field%minB
  ! WRITE(*,*)psi_surf(1),field%minB
end do
CALL active_tracer%delete
CALL field%delete()
DEALLOCATE(field)
!$omp end parallel
self%xmin=self%func%xs(0)+(self%func%xs(0)-self%func%xs(1))*4.d0
self%xmax=self%func%xs(self%npsi)+(self%func%xs(1)-self%func%xs(0))*4.d0
self%yp1=0.d0
self%ypn=0.d0
!---Setup Spline
CALL spline_fit(self%func,"extrap")
DO k=1,omp_get_max_threads()
  CALL spline_copy(self%func,self%fun_loc(k))
END DO
IF(oft_debug_print(2))CALL oft_decrease_indent
CALL psi_int%delete()
end subroutine dipole_b0_update
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine minbinv_apply(self,cell,f,gop,val)
class(minbinv_interp), intent(inout) :: self
integer(4), intent(in) :: cell
real(8), intent(in) :: f(:)
real(8), intent(in) :: gop(3,3)
real(8), intent(out) :: val(:)
integer(4), allocatable :: j(:)
integer(4) :: jc
real(8) :: rop(3),d2op(6),pt(3),grad(3),d2(3),d2_tmp(6),tmp
real(8) :: g2op(6,6),s,c
real(8) :: K(6,3)
!---Get dofs
allocate(j(self%lag_rep%nce))
call self%lag_rep%ncdofs(cell,j)
!---Reconstruct gradient
call self%mesh%hessian(cell,f,g2op,K)
grad=0.d0
d2_tmp=0.d0
do jc=1,self%lag_rep%nce
  call oft_blag_geval(self%lag_rep,cell,jc,f,rop,gop)
  call oft_blag_d2eval(self%lag_rep,cell,jc,f,d2op,g2op)
  grad=grad+self%uvals(j(jc))*rop
  d2_tmp=d2_tmp+self%uvals(j(jc))*(d2op-MATMUL(g2op,MATMUL(K,rop)))
end do
d2=d2_tmp([1,2,4]) ! Map from 3D to 2D
!---Get radial position
pt=self%mesh%log2phys(cell,f)
!---
s=SIN(self%t)
c=COS(self%t)
!---
val(1)=(self%rho*(grad(1)*s-grad(2)*c))/(grad(1)*c+grad(2)*s)
val(2)=pt(1)*SQRT((self%rho**2+val(1)**2)/SUM(grad**2))
self%minB=min(self%minB,SQRT((grad(1)/pt(1))**2+(grad(2)/pt(1))**2))
deallocate(j)
end subroutine minbinv_apply
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine dipole_ani_copy(self,new,new_gs)
class(dipole_ani_press), intent(inout) :: self
class(gs_ani_press), pointer, intent(inout) :: new
class(gs_equil), target, intent(inout) :: new_gs
ALLOCATE(new, MOLD=self)
SELECT TYPE(new)
  CLASS IS(dipole_ani_press)
    new%a_exp = self%a_exp
    IF(ASSOCIATED(self%psi_eval))THEN
      CALL new%setup(new_gs)
      CALL new%update(new_gs)
    END IF
  CLASS DEFAULT
    CALL oft_abort('Allocation produced wrong type','dipole_ani_copy',__FILE__)
END SELECT
end subroutine dipole_ani_copy
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine dipole_ani_setup(self,gs)
class(dipole_ani_press), intent(inout) :: self
class(gs_equil), target, intent(inout) :: gs
self%gs=>gs
self%mesh=>gs%device%mesh
ALLOCATE(self%psi_eval,self%psi_geval)
self%psi_eval%u=>self%gs%psi
CALL self%psi_eval%setup(self%gs%device%fe_rep)
CALL self%psi_geval%shared_setup(self%psi_eval)
ALLOCATE(dipole_b0_flux_func::self%B0_prof)
CALL create_dipole_b0_prof(self%B0_prof,64)
end subroutine dipole_ani_setup
!------------------------------------------------------------------------------
!> Destroy temporary internal storage and nullify references
!------------------------------------------------------------------------------
subroutine dipole_ani_delete(self)
class(dipole_ani_press), intent(inout) :: self
INTEGER(i4) :: i
IF(ASSOCIATED(self%psi_eval))THEN
  CALL self%psi_eval%delete()
  CALL self%psi_geval%delete()
  DEALLOCATE(self%psi_eval,self%psi_geval)
END IF
IF(ASSOCIATED(self%B0_prof))THEN
  SELECT TYPE(this=>self%B0_prof)
  TYPE IS(dipole_b0_flux_func)
    CALL spline_dealloc(this%func)
    DO i=1,omp_get_max_threads()
      CALL spline_dealloc(this%fun_loc(i))
    END DO
    DEALLOCATE(this%fun_loc)
  END SELECT
END IF
NULLIFY(self%gs,self%mesh)
end subroutine dipole_ani_delete
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine dipole_ani_update(self,gseq)
class(dipole_ani_press), intent(inout) :: self
class(gs_equil), intent(inout) :: gseq
CALL self%B0_prof%update(gseq)
end subroutine dipole_ani_update
!------------------------------------------------------------------------------
!> Reconstruct a component of a Grad-Shafranov solution
!------------------------------------------------------------------------------
subroutine dipole_ani_apply(self,cell,f,gop,val)
class(dipole_ani_press), intent(inout) :: self !< Interpolation object
integer(4), intent(in) :: cell !< Cell for interpolation
real(8), intent(in) :: f(:) !< Position in cell in logical coord [3]
real(8), intent(in) :: gop(3,3) !< Logical gradient vectors at f [3,3]
real(8), intent(out) :: val(:) !< Reconstructed [p_par, p_perp] factors [2]
real(8) :: pt(3),psitmp(1),gpsitmp(3),Bp,H
pt=self%gs%device%fe_rep%mesh%log2phys(cell,f)
CALL self%psi_eval%interp(cell,f,gop,psitmp)
CALL self%psi_geval%interp(cell,f,gop,gpsitmp)
Bp = SQRT((gpsitmp(1)/(pt(1)+gs_epsilon))**2 + (gpsitmp(2)/(pt(1)+gs_epsilon))**2)
H = (self%B0_prof%f(psitmp(1))/Bp)**(2.d0*self%a_exp)
val=[H/(1.d0+2.d0*self%a_exp), H]
end subroutine dipole_ani_apply
!---------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE create_mirror_b0_prof(func,npsi)
CLASS(flux_func), INTENT(inout) :: func
INTEGER(4), INTENT(in) :: npsi
INTEGER(4) :: i
! IF(.NOT.ASSOCIATED(func))ALLOCATE(mirror_b0_flux_func::func)
select type(self=>func)
  type is(mirror_b0_flux_func)
    !---
    self%npsi=npsi
    CALL spline_alloc(self%func,self%npsi,1)
    ALLOCATE(self%fun_loc(omp_get_max_threads()))
    DO i=1,omp_get_max_threads()
      CALL spline_alloc(self%fun_loc(i),self%npsi,1)
    END DO
class default
  CALL oft_abort('Invalid flux function type in create_mirror_b0_prof','create_mirror_b0_prof',__FILE__)
end select
END SUBROUTINE create_mirror_b0_prof
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mirror_b0_save_hdf5(self,filename,path)
class(mirror_b0_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
IF(.NOT.hdf5_field_exist(filename,path//'/TYPE'))CALL hdf5_write('mirror_b0',filename,path//'/TYPE')
CALL hdf5_write(self%npsi,filename,path//'/NPSI')
CALL hdf5_write(self%z_midplane,filename,path//'/Z_MIDPLANE')
end subroutine mirror_b0_save_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mirror_b0_save_txt(self,io_unit)
class(mirror_b0_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
WRITE(io_unit,*)'mirror_b0'
WRITE(io_unit,*)self%npsi
WRITE(io_unit,*)self%z_midplane
end subroutine mirror_b0_save_txt
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mirror_b0_load_hdf5(self,filename,path,success)
class(mirror_b0_flux_func), intent(inout) :: self
character(LEN=*), intent(in) :: filename
character(LEN=*), intent(in) :: path
logical, intent(out) :: success
integer(i4) :: npsi
real(r8), allocatable :: xvals(:),yvals(:)
CALL hdf5_read(npsi,filename,path//'/NPSI',success=success)
IF(.NOT.success)RETURN
CALL hdf5_read(self%z_midplane,filename,path//'/Z_MIDPLANE',success=success)
IF(.NOT.success)RETURN
CALL create_mirror_b0_prof(self,npsi)
end subroutine mirror_b0_load_hdf5
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mirror_b0_load_txt(self,io_unit)
class(mirror_b0_flux_func), intent(inout) :: self
integer, intent(in) :: io_unit
integer(i4) :: npsi
READ(io_unit,*)npsi
READ(io_unit,*)self%z_midplane
CALL create_mirror_b0_prof(self,npsi)
end subroutine mirror_b0_load_txt
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mirror_b0_copy(self,new)
class(mirror_b0_flux_func), intent(inout) :: self
class(flux_func), pointer, intent(inout) :: new
CALL spline_func_copy(self,new)
SELECT TYPE(new)
  TYPE IS(mirror_b0_flux_func)
    new%plasma_bounds=self%plasma_bounds
    new%f_offset=self%f_offset
    new%z_midplane = self%z_midplane
END SELECT
end subroutine mirror_b0_copy
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mirror_b0_delete(self)
class(mirror_b0_flux_func), intent(inout) :: self
integer(i4) :: i
CALL spline_func_delete(self)
end subroutine mirror_b0_delete
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mirror_b0_update(self,gseq)
class(mirror_b0_flux_func), intent(inout) :: self
class(gs_equil), intent(inout) :: gseq
type(minbinv_interp), target :: field
type(oft_lag_brinterp) :: psi_int
type(oft_lag_bginterp) :: psi_gint
real(8) :: I,Ip,q,qp,vp,vpp,s,a,b,pp,gop(3,3),psi_surf(1),vol,rmax,gpsi(3)
real(8) :: raxis,zaxis,f(3),pt(3),x1,x2,xr
real(8), pointer, dimension(:) :: v
real(8), parameter :: tol=1.d-10
integer(4) :: j,k,cell
!
IF(oft_debug_print(2))THEN
  WRITE(*,'(2A)')oft_indent,'Updating Mirror B0 profile:'
  CALL oft_increase_indent
END IF
psi_int%u=>gseq%psi
CALL psi_int%setup(gseq%device%fe_rep)
CALL psi_gint%shared_setup(psi_int)
x1=0.d0; x2=1.d0
IF(gseq%plasma_bounds(1)>-1.d98)THEN
  x1=gseq%plasma_bounds(1); x2=gseq%plasma_bounds(2)
END IF
xr = (x2-x1)
!---Find Rmax along Zaxis
rmax=gseq%device%rmax
cell=0
DO j=1,100
  pt=[rmax*j/REAL(100,8),self%z_midplane,0.d0]
  CALL bmesh_findcell(gseq%device%mesh,cell,pt,f)
  IF( (MAXVAL(f)>1.d0+tol) .OR. (MINVAL(f)<-tol) )EXIT
  CALL psi_int%interp(cell,f,gop,psi_surf)
  IF( psi_surf(1) > x1)EXIT
  rmax=pt(1)
END DO
!
pt=[0.5d0*rmax,self%z_midplane,0.d0]
do j=1,self%npsi+1
  psi_surf(1)=(x2-x1)*(1.d0-j/REAL(self%npsi+2,4))
  psi_surf(1)=x2 - psi_surf(1)
  CALL gs_psi2r(gseq,psi_surf(1),pt)
  CALL bmesh_findcell(gseq%device%mesh,cell,pt,f)
  call gseq%device%mesh%jacobian(cell,f,gop,vol)
  CALL psi_gint%interp(cell,f,gop,gpsi)
  !---Get surface minB
  self%func%xs(j-1)=psi_surf(1)
  self%func%fs(j-1,1)=SQRT((gpsi(1)/(pt(1)+gs_epsilon))**2 + (gpsi(2)/(pt(1)+gs_epsilon))**2)
end do
self%yp1=0.d0
self%ypn=0.d0
self%xmin=self%func%xs(0)+(self%func%xs(0)-self%func%xs(1))*4.d0
self%xmax=self%func%xs(self%npsi)+(self%func%xs(1)-self%func%xs(0))*4.d0
!---Setup Spline
CALL spline_fit(self%func,"extrap")
DO k=1,omp_get_max_threads()
  CALL spline_copy(self%func,self%fun_loc(k))
END DO
IF(oft_debug_print(2))CALL oft_decrease_indent
CALL psi_int%delete()
CALL psi_gint%delete()
end subroutine mirror_b0_update
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mirror_slosh_copy(self,new,new_gs)
class(mirror_ani_slosh), intent(inout) :: self
class(gs_ani_press), pointer, intent(inout) :: new
class(gs_equil), target, intent(inout) :: new_gs
ALLOCATE(new, MOLD=self)
SELECT TYPE(new)
  CLASS IS(mirror_ani_slosh)
    new%n_exp = self%n_exp
    new%bturn = self%bturn
    new%zthroat = self%zthroat
    IF(ASSOCIATED(self%psi_eval))THEN
      CALL new%setup(new_gs)
      CALL new%update(new_gs)
    END IF
END SELECT
end subroutine mirror_slosh_copy
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mirror_slosh_setup(self,gs)
class(mirror_ani_slosh), intent(inout) :: self
class(gs_equil), target, intent(inout) :: gs
self%gs=>gs
self%mesh=>gs%device%mesh
ALLOCATE(self%psi_eval,self%psi_geval)
self%psi_eval%u=>self%gs%psi
CALL self%psi_eval%setup(self%gs%device%fe_rep)
CALL self%psi_geval%shared_setup(self%psi_eval)
ALLOCATE(mirror_b0_flux_func::self%B0_prof)
CALL create_mirror_b0_prof(self%B0_prof,64)
end subroutine mirror_slosh_setup
!------------------------------------------------------------------------------
!> Destroy temporary internal storage and nullify references
!------------------------------------------------------------------------------
subroutine mirror_slosh_delete(self)
class(mirror_ani_slosh), intent(inout) :: self
INTEGER(i4) :: i
IF(ASSOCIATED(self%psi_eval))THEN
  CALL self%psi_eval%delete()
  CALL self%psi_geval%delete()
  DEALLOCATE(self%psi_eval,self%psi_geval)
END IF
IF(ASSOCIATED(self%B0_prof))THEN
  SELECT TYPE(this=>self%B0_prof)
  TYPE IS(mirror_b0_flux_func)
    CALL spline_dealloc(this%func)
    DO i=1,omp_get_max_threads()
      CALL spline_dealloc(this%fun_loc(i))
    END DO
    DEALLOCATE(this%fun_loc)
  END SELECT
END IF
NULLIFY(self%gs,self%mesh)
end subroutine mirror_slosh_delete
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mirror_slosh_update(self,gseq)
class(mirror_ani_slosh), intent(inout) :: self
class(gs_equil), intent(inout) :: gseq
CALL self%B0_prof%update(gseq)
end subroutine mirror_slosh_update
!------------------------------------------------------------------------------
!> Reconstruct a component of a Grad-Shafranov solution
!------------------------------------------------------------------------------
subroutine mirror_slosh_apply(self,cell,f,gop,val)
class(mirror_ani_slosh), intent(inout) :: self !< Interpolation object
integer(4), intent(in) :: cell !< Cell for interpolation
real(8), intent(in) :: f(:) !< Position in cell in logical coord [3]
real(8), intent(in) :: gop(3,3) !< Logical gradient vectors at f [3,3]
real(8), intent(out) :: val(:) !< Reconstructed [p_par, p_perp] factors [2]
real(8) :: pt(3),psitmp(1),gpsitmp(3),b_bar,b_turn
!---
pt=self%gs%device%fe_rep%mesh%log2phys(cell,f)
CALL self%psi_eval%interp(cell,f,gop,psitmp)
CALL self%psi_geval%interp(cell,f,gop,gpsitmp)
b_turn = self%B0_prof%f(psitmp(1))*self%bturn
b_bar = SQRT((gpsitmp(1)/(pt(1)+gs_epsilon))**2 + (gpsitmp(2)/(pt(1)+gs_epsilon))**2)/b_turn
IF((b_bar<=1.d0).AND.(ABS(pt(2))<self%zthroat))THEN
  val=[b_bar/self%n_exp*(1.d0-b_bar)**self%n_exp, b_bar*b_bar*(1.d0-b_bar)**(self%n_exp-1.d0)]
ELSE
  val=[0.d0,0.d0]
END IF
end subroutine mirror_slosh_apply
end module grad_shaf_prof_phys
