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
use oft_base
use oft_mesh_type, only: oft_bmesh, bmesh_findcell
USE oft_la_base, ONLY: oft_vector
USE fem_utils, ONLY: bfem_interp
use oft_lag_basis, only: oft_blag_d2eval, oft_blag_geval
USE oft_blag_operators, only: oft_lag_brinterp, oft_lag_bginterp
use oft_gs, only: gs_equil, flux_func, gs_psi2r, gs_itor_nl, oft_indent, &
  oft_increase_indent, oft_decrease_indent, gsinv_interp, gs_prof_interp, &
  gs_get_qprof, gs_ani_press, gs_epsilon, sauter_fc
use tracing_2d, only: set_tracer, active_tracer, tracinginv_fs
use oft_gs_profiles, only: spline_flux_func, linterp_flux_func, linterp_copy, spline_func_copy
use spline_mod
USE mhd_utils, ONLY: mu0
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
  procedure :: update => mercier_update
end type mercier_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(linterp_flux_func) :: jphi_flux_func
  integer(4) :: ngeom = 50 !< Number of points in psi for <R>, <1/R> evaluation
  real(8) :: j0 = 0.d0 !< LCFS Jphi value
  real(8) :: norm_last = 1.d0 !< Last Jphi normalization factor (for Ip target)
  real(8) :: alpha_last = 1.d0 !< Alpha (input Jphi rescaling factor) from previous NL iteration
  real(8) :: rescale_last = 1.d0 !< Damped jphi_rescale from previous NL iteration, to ensure Ip target is met
  logical :: freeze_j_BS = .FALSE. !< Set .TRUE. once j_BS stagnates for 2 steps; skips Sauter call (big speedup)
  logical :: freeze_alpha = .FALSE.   !< Set .TRUE. once dalpha stagnates for 2 steps; skips alpha re-solve (speedup)
  real(8) :: djBS_tol = 1.0e-6_r8    !< Hard RMS tolerance; freeze immediately if djBS drops below this
  real(8) :: dalpha_tol = 1.0e-6_r8  !< Hard tolerance; freeze alpha immediately if dalpha drops below this
  integer(4) :: djBS_no_improve = 0   !< Consecutive steps with non-decreasing djBS
  real(8) :: djBS_min = huge(1.0d0)  !< Running minimum djBS seen so far
  integer(4) :: dalpha_no_improve = 0 !< Consecutive steps with non-decreasing dalpha
  real(8) :: dalpha_min = huge(1.0d0) !< Running minimum dalpha seen so far
  real(8), pointer, dimension(:) :: jphi => NULL() !< Jphi(psi) profile values
  real(8), pointer, dimension(:) :: jphi_total_last => NULL() !< Assembled jphi_total from previous NL iteration
  real(8), pointer, dimension(:) :: j_BS_last => NULL() !< j_BS profile from previous NL iteration (for freeze check)
  !> Update mode selector: 0 = inductive only (default), 1 = bootstrap-coupled
  integer(4) :: update_mode = 0
contains
  !> Needs docs
  procedure :: copy => jphi_copy
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
  procedure :: copy => dipole_b0_copy
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
  procedure :: copy => mirror_b0_copy
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
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE create_jphi_ff(func,npsi,psivals,yvals,y0,update_mode)
CLASS(flux_func), POINTER, INTENT(out) :: func
INTEGER(4), INTENT(in) :: npsi
REAL(8), INTENT(in) :: psivals(npsi)
REAL(8), INTENT(in) :: yvals(npsi)
REAL(8), INTENT(in) :: y0
INTEGER(4), OPTIONAL, INTENT(in) :: update_mode
INTEGER(4) :: i,ierr
ALLOCATE(jphi_flux_func::func)
SELECT TYPE(self=>func)
  TYPE IS(jphi_flux_func)
  !---
  self%npsi=npsi
  self%ncofs=self%npsi
  !---
  ALLOCATE(self%x(self%npsi))
  ALLOCATE(self%yp(self%npsi))
  ALLOCATE(self%y(self%npsi))
  ALLOCATE(self%jphi(self%npsi))
  !---
  self%j0=y0
  self%y0=0.d0
  IF (PRESENT(update_mode)) THEN
    self%update_mode = update_mode
  END IF
  DO i=1,self%npsi
    self%x(i) = psivals(i)
    self%jphi(i) = yvals(i)
    self%yp(i) = psivals(i) ! Dummy initialization
  END DO
  self%yp = self%yp/(SUM(ABS(self%yp))/REAL(self%npsi,8)) ! Consistent (hopefully) normalization
  ierr=self%set_cofs(self%yp)
  IF(oft_debug_print(1))WRITE(*,*)'Jphi linear interpolator Created',self%ncofs,self%x,self%j0
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
    new%ngeom = self%ngeom
    new%j0 = self%j0
    new%norm_last = self%norm_last
    new%alpha_last = self%alpha_last
    new%rescale_last = self%rescale_last
    new%freeze_j_BS = self%freeze_j_BS
    new%freeze_alpha = self%freeze_alpha
    new%djBS_tol = self%djBS_tol
    new%dalpha_tol = self%dalpha_tol
    new%djBS_no_improve = self%djBS_no_improve
    new%djBS_min = self%djBS_min
    new%dalpha_no_improve = self%dalpha_no_improve
    new%dalpha_min = self%dalpha_min
    new%update_mode = self%update_mode
    ALLOCATE(new%jphi, SOURCE=self%jphi)
    IF(ASSOCIATED(self%jphi_total_last)) ALLOCATE(new%jphi_total_last, SOURCE=self%jphi_total_last)
    IF(ASSOCIATED(self%j_BS_last)) ALLOCATE(new%j_BS_last, SOURCE=self%j_BS_last)
END SELECT
end subroutine jphi_copy
!------------------------------------------------------------------------------
!> Update F*F' profile — dispatches on update_mode.
!> mode=0 (default): standard jphi profile update.
!> mode=1: bootstrap-coupled update.
!------------------------------------------------------------------------------
subroutine jphi_update(self,gseq)
class(jphi_flux_func), intent(inout) :: self
class(gs_equil), intent(inout) :: gseq
IF(self%update_mode==1 .AND. gseq%Itor_target<=0.d0)THEN
  CALL jphi_bs_update(self,gseq)
ELSE
  CALL jphi_update_default(self,gseq)
END IF
end subroutine jphi_update
!------------------------------------------------------------------------------
!> Build the <R> / <1/R> spline needed for jphi -> F*F' mapping.
!> Allocates and fits R_spline on ngeom points from the current equilibrium.
!> Caller is responsible for calling spline_dealloc(R_spline) when done.
!------------------------------------------------------------------------------
SUBROUTINE build_Ravg_spline(gseq, ngeom, R_spline)
CLASS(gs_equil), INTENT(inout) :: gseq
INTEGER(i4), INTENT(in) :: ngeom
TYPE(spline_type), INTENT(out) :: R_spline
INTEGER(i4) :: i
REAL(r8), ALLOCATABLE :: ravgs(:,:), psi_q(:), qprof(:)
REAL(r8), PARAMETER :: psi_pad = 1.d-3
ALLOCATE(ravgs(ngeom,2), psi_q(ngeom), qprof(ngeom))
psi_q = [(REAL(i-1,r8)/REAL(ngeom,r8), i=1,ngeom)]
IF(gseq%diverted)THEN
  psi_q(1) = MIN(psi_q(2), psi_pad)
  CALL gs_get_qprof(gseq, ngeom, psi_q, qprof, ravgs=ravgs)
  psi_q(1) = 0.d0
  ravgs(1,1) = gseq%lim_point(1)
  ravgs(1,2) = 1.d0/gseq%lim_point(1)
ELSE
  CALL gs_get_qprof(gseq, ngeom, psi_q, qprof, ravgs=ravgs)
END IF
CALL spline_alloc(R_spline, ngeom-1, 2)
R_spline%xs(0:ngeom-2) = psi_q(1:ngeom-1); R_spline%xs(ngeom-1) = 1.d0
R_spline%fs(0:ngeom-2,1) = ravgs(1:ngeom-1,1)
R_spline%fs(ngeom-1,1) = gseq%o_point(1)
R_spline%fs(0:ngeom-2,2) = ravgs(1:ngeom-1,2)
R_spline%fs(ngeom-1,2) = 1.d0/gseq%o_point(1)
CALL spline_fit(R_spline, "extrap")
DEALLOCATE(ravgs, psi_q, qprof)
END SUBROUTINE build_Ravg_spline
!------------------------------------------------------------------------------
!> Evaluate the geometric factor qtmp(i) = <R>(psi_i) * <1/R>(psi_i) on an
!> arbitrary psi_N grid by calling spline_eval on a pre-built R_spline.
!------------------------------------------------------------------------------
SUBROUTINE eval_R_qtmp(R_spline, psi_vals, n, qtmp)
TYPE(spline_type), INTENT(inout) :: R_spline
INTEGER(i4), INTENT(in) :: n
REAL(r8), INTENT(in) :: psi_vals(n)
REAL(r8), INTENT(out) :: qtmp(n)
INTEGER(i4) :: i
DO i = 1, n
  CALL spline_eval(R_spline, psi_vals(i), 0)
  qtmp(i) = R_spline%f(1) * R_spline%f(2)
END DO
END SUBROUTINE eval_R_qtmp
!------------------------------------------------------------------------------
!> Update F*F' profile from Jphi, P', and current equilibrium (default mode=0).
!------------------------------------------------------------------------------
subroutine jphi_update_default(self,gseq)
class(jphi_flux_func), intent(inout) :: self
class(gs_equil), intent(inout) :: gseq
INTEGER(i4) :: i
REAL(r8) :: jphi_norm,pscale,pprime
REAL(r8), ALLOCATABLE :: qtmp(:)
type(spline_type) :: R_spline
self%plasma_bounds=gseq%plasma_bounds
IF(gseq%mode/=1)CALL oft_abort("Jphi profile requires (F^2)' formulation","jphi_update",__FILE__)
! IF(gseq%Itor_target<0.d0)CALL oft_abort("Jphi profile requires Ip target","jphi_update",__FILE__)
IF(gseq%pax_target<0.d0)CALL oft_abort("Jphi profile requires Pax target","jphi_update",__FILE__)
!---Get updated flux surface geometry for Jphi -> F*F' mapping
CALL build_Ravg_spline(gseq, self%ngeom, R_spline)
!---Update jphi normalization to match Ip target
IF(gseq%Itor_target>0.d0)THEN
  ALLOCATE(qtmp(self%npsi))
  CALL eval_R_qtmp(R_spline, self%x, self%npsi, qtmp)
  CALL gs_flux_int(gseq,self%x,self%jphi/qtmp,self%npsi,jphi_norm)
  DEALLOCATE(qtmp)
  jphi_norm=ABS(gseq%Itor_target)/jphi_norm
  self%norm_last=jphi_norm
ELSE
  CALL gs_itor_nl(gseq,jphi_norm)
  jphi_norm=(ABS(gseq%Itor_target)*mu0/jphi_norm + self%norm_last)/2.d0
  self%norm_last=jphi_norm
END IF
!---Get pressure profile
CALL gseq%P%update(gseq) ! Make sure pressure profile is up to date with EQ
IF(ASSOCIATED(gseq%P_ani))CALL oft_abort('Jphi profiles do not support anistopic pressure','jphi_update',__FILE__) !CALL gseq%P_ani%update(gseq)
pscale=gseq%P%f(gseq%plasma_bounds(2))
pscale=gseq%pax_target/pscale
!---Compute updated F*F' profile ! 2.0*(jtor -  R_avg * (-pprime)) * (mu0 / one_over_R_avg)
CALL spline_eval(R_spline,0.d0,0)
pprime=gseq%P%fp(gseq%plasma_bounds(1))
self%y0 = 2.d0*(self%j0*jphi_norm - R_spline%f(1)*pprime*pscale)/R_spline%f(2)
DO i=1,self%npsi
  CALL spline_eval(R_spline,self%x(i),0)
  pprime=gseq%P%fp(self%x(i)*(gseq%plasma_bounds(2)-gseq%plasma_bounds(1))+gseq%plasma_bounds(1))
  self%yp(i) = 2.d0*(self%jphi(i)*jphi_norm - R_spline%f(1)*pprime*pscale)/R_spline%f(2)
END DO
! Disable Ip matching and fix F*F' scale (matching is done here instead)
IF(gseq%Itor_target>0.d0)gseq%Itor_target=-gseq%Itor_target
gseq%ffp_scale=1.d0
!---Clean up
CALL spline_dealloc(R_spline)
i=self%set_cofs(self%yp)
end subroutine jphi_update_default
!---------------------------------------------------------------------------------
!> Update F*F' profile from inductive Jphi coupled with bootstrap current.
!>
!> Each call (one NL iteration):
!>   1. Build <R>/<1/R> spline on self%x; pre-compute qtmp = <R>*<1/R>.
!>   2. Evaluate bootstrap current j_BS on self%x (or reuse cache if frozen).
!>   3. Apply edge taper to j_BS and jphi_ind component-wise.
!>   4. Compute jphi_rescale to reconcile gs_itor_nl vs gs_flux_int.
!>   5. Solve analytically for alpha: gs_flux_int is linear in alpha, so two
!>      evaluations (alpha=0, alpha=1) give alpha = (Ip_target - Ip_lo)/(Ip_hi - Ip_lo).
!>   6. Assemble jphi_total = alpha*jphi_ind + j_BS; compute F*F' knots.
!>   7. Diagnostics (if diagnose_bs is set).
!---------------------------------------------------------------------------------
SUBROUTINE jphi_bs_update(self, gseq)
CLASS(jphi_flux_func), INTENT(inout) :: self
CLASS(gs_equil), INTENT(inout) :: gseq
INTEGER(i4) :: i
REAL(r8) :: pscale, pprime
REAL(r8), ALLOCATABLE :: qtmp(:)
TYPE(spline_type) :: R_spline
! Bootstrap arrays (on self%x grid)
REAL(r8), ALLOCATABLE :: j_BS(:)
! Optional edge-spike workspace (only allocated when isolate_edge_jBS or parameterize_jBS is set)
REAL(r8), ALLOCATABLE :: j_spike_tmp(:)
! Working arrays on self%x grid
REAL(r8), ALLOCATABLE :: jphi_total(:)
REAL(r8), ALLOCATABLE :: jphi_ind(:)  !< Tapered copy of self%jphi (= self%jphi when taper off)
! Alpha-solve scalars
REAL(r8) :: alpha, ip_target, ip_ind, ip_result_lo, ip_result_hi, dalpha
! Relative change in bootstrap current for freeze check
REAL(r8) :: djBS
! gs_itor_nl / gs_flux_int reconciliation
REAL(r8) :: itor_nl = 0.0_r8, itor_flint = 0.0_r8, jphi_rescale
CHARACTER(len=80) :: char_buf
!---
self%plasma_bounds = gseq%plasma_bounds
IF(gseq%mode/=1) &
  CALL oft_abort("Jphi-BS profile requires (F^2)' formulation", &
                 "jphi_bs_update",__FILE__)
IF(gseq%pax_target<0.d0) &
  CALL oft_abort("Jphi-BS profile requires Pax target", &
                 "jphi_bs_update",__FILE__)
IF(.NOT.ASSOCIATED(gseq%Te)) &
  CALL oft_abort("Jphi-BS profile requires Te profile", &
                 "jphi_bs_update",__FILE__)
IF(.NOT.ASSOCIATED(gseq%Ti)) &
  CALL oft_abort("Jphi-BS profile requires Ti profile", &
                 "jphi_bs_update",__FILE__)
IF(.NOT.ASSOCIATED(gseq%ne)) &
  CALL oft_abort("Jphi-BS profile requires ne profile", &
                 "jphi_bs_update",__FILE__)
IF(.NOT.ASSOCIATED(gseq%ni)) &
  CALL oft_abort("Jphi-BS profile requires ni profile", &
                 "jphi_bs_update",__FILE__)
!--- 1. Build <R>/<1/R> spline; pre-compute qtmp = <R>*<1/R> on self%x.
!   R_spline stays alive until after the F*F' loop (step 6).
ALLOCATE(qtmp(self%npsi))
CALL build_Ravg_spline(gseq, self%ngeom, R_spline)
DO i = 1, self%npsi
  CALL spline_eval(R_spline, self%x(i), 0)
  qtmp(i) = R_spline%f(1) * R_spline%f(2)
END DO
!--- 2. Bootstrap current on self%x grid.
ALLOCATE(j_BS(self%npsi))
CALL gseq%P%update(gseq)
IF(self%freeze_j_BS .AND. ASSOCIATED(self%j_BS_last)) THEN
  !--- Frozen: reuse cached j_BS.
  j_BS = self%j_BS_last
  djBS = 0.0_r8
ELSE
  !--- Not frozen: run full bootstrap calculation (Sauter).
  IF (gseq%boot_ops%isolate_edge_jBS .OR. gseq%boot_ops%parameterize_jBS) THEN
    ALLOCATE(j_spike_tmp(self%npsi))
    CALL calculate_bootstrap(gseq, self%npsi, self%x, j_BS, &
        isolate_edge_jBS=gseq%boot_ops%isolate_edge_jBS, &
        parameterize_jBS=gseq%boot_ops%parameterize_jBS, &
        scale_jBS=gseq%boot_ops%scale_jBS, &
        j_spike=j_spike_tmp)
    IF(gseq%boot_ops%diagnose_bs)THEN
      WRITE(*,'(A)') '  [diagnose_bs] i  psi_N         j_BS(bulk)[A/m2]  j_spike[A/m2]   jphi[A/m2]'
      DO i = 1, self%npsi
        WRITE(*,'(A,I4,4ES15.5)') '  ', i, self%x(i), j_BS(i), j_spike_tmp(i), self%jphi(i)
      END DO
    END IF
    j_BS = j_spike_tmp
    DEALLOCATE(j_spike_tmp)
  ELSE
    CALL calculate_bootstrap(gseq, self%npsi, self%x, j_BS)
    j_BS = j_BS * gseq%boot_ops%scale_jBS
    IF(gseq%boot_ops%diagnose_bs)THEN
      WRITE(*,'(A)') '  [diagnose_bs] i  psi_N         j_BS[A/m2]      jphi[A/m2]'
      DO i = 1, self%npsi
        WRITE(*,'(A,I4,3ES15.5)') '  ', i, self%x(i), j_BS(i), self%jphi(i)
      END DO
    END IF
  END IF
  !   calculate_bootstrap returns j_BS in raw A/m², multiply by mu0.
  j_BS = j_BS * mu0
  !--- 2a. Freeze check: freeze j_BS if RMS change drops below tol, or stagnates.
  IF(ASSOCIATED(self%j_BS_last)) THEN
    djBS = SQRT(SUM((j_BS - self%j_BS_last)**2) / REAL(self%npsi,r8)) / &
           MAX(SQRT(SUM(j_BS**2) / REAL(self%npsi,r8)), 1.0e-30_r8)
    IF(djBS < self%djBS_tol .OR. self%freeze_alpha) THEN
      self%freeze_j_BS = .TRUE.
    ELSE IF(djBS >= self%djBS_min) THEN
      self%djBS_no_improve = self%djBS_no_improve + 1
      IF(self%djBS_no_improve >= 2) THEN
        self%freeze_j_BS = .TRUE.
        IF(djBS > self%djBS_tol) THEN
          WRITE(char_buf,'(A,ES12.4,A,ES12.4,A)') &
            'Bootstrap solution convergence stalled,' // &
            ' relative change per nonlinear step = ', djBS, &
            ', above recommended tolerance (djBS_tol=', self%djBS_tol, ')'
          CALL oft_warn(TRIM(char_buf))
        END IF
      END IF
    ELSE
      self%djBS_no_improve = 0
      self%djBS_min = djBS
    END IF
  ELSE
    djBS = HUGE(1.0_r8)
  END IF
  IF(.NOT.ASSOCIATED(self%j_BS_last)) ALLOCATE(self%j_BS_last(self%npsi))
  self%j_BS_last = j_BS
END IF
!--- 3. Apply edge taper to j_BS and jphi_ind.
!   self%j_BS_last caches the un-tapered j_BS so freeze comparisons track physics.
!   taper_edge_psi0 is in standard convention (0=axis,1=LCFS);
!   threshold in OFT convention (0=LCFS,1=axis) is (1 - taper_edge_psi0).
ALLOCATE(jphi_ind(self%npsi))
jphi_ind = self%jphi
IF (gseq%boot_ops%taper_edge_jBS) THEN
  CALL apply_edge_taper(self%npsi, self%x, j_BS, &
                        1.0_r8 - gseq%boot_ops%taper_edge_psi0, &
                        gseq%boot_ops%taper_edge_shape, &
                        oft_psi_conv=.TRUE.)
  CALL apply_edge_taper(self%npsi, self%x, jphi_ind, &
                        1.0_r8 - gseq%boot_ops%taper_edge_psi0, &
                        gseq%boot_ops%taper_edge_shape, &
                        oft_psi_conv=.TRUE.)
END IF
ALLOCATE(jphi_total(self%npsi))
!--- 4. Reconcile gs_itor_nl vs gs_flux_int when Itor_target < 0.
!   No Ip target: rescale jphi_total so the integrated current matches the
!   FEM solution (gs_itor_nl) rather than the profile quadrature (gs_flux_int).
jphi_rescale = self%rescale_last
IF(gseq%Itor_target < 0.d0 .AND. .NOT. self%freeze_j_BS) THEN
  CALL gs_itor_nl(gseq, itor_nl)
  IF(ASSOCIATED(self%jphi_total_last)) THEN
    CALL gs_flux_int(gseq, self%x, self%jphi_total_last/qtmp, self%npsi, itor_flint)
    IF(ABS(itor_flint) > 0.d0) THEN
      jphi_rescale = (itor_nl/itor_flint + self%rescale_last) / 2.0_r8
    END IF
  END IF
  self%rescale_last = jphi_rescale
END IF
!--- 5. Solve analytically for alpha.
!   gs_flux_int is linear in alpha; two evaluations (alpha=0 and alpha=1) give
!   alpha = (Ip_target - Ip_lo) / (Ip_hi - Ip_lo).  Skip once frozen.
ip_target = ABS(gseq%Itor_target)/jphi_rescale
IF(self%freeze_j_BS .OR. self%freeze_alpha) THEN
  !--- Frozen: reuse last converged alpha.
  alpha = self%alpha_last
  dalpha = 0.0_r8
  ip_result_lo = 0.0_r8
  ip_result_hi = 0.0_r8
ELSE
  !--- Not yet frozen: exact linear solve for alpha.
  jphi_total = j_BS
  CALL gs_flux_int(gseq, self%x, jphi_total/qtmp, self%npsi, ip_result_lo)
  jphi_total = jphi_ind + j_BS
  CALL gs_flux_int(gseq, self%x, jphi_total/qtmp, self%npsi, ip_result_hi)
  ip_ind = ip_result_hi - ip_result_lo
  IF(ABS(ip_ind) > 0.0_r8)THEN
    alpha = (ip_target - ip_result_lo) / ip_ind
    IF(alpha < -1.0_r8 .OR. alpha > 10.0_r8) THEN
      WRITE(char_buf,'(A,ES12.4)') '[jphi_bs_update] WARNING: alpha out of expected range [-1,10]: alpha=', alpha
      CALL oft_warn(TRIM(char_buf))
    END IF
  ELSE
    alpha = self%alpha_last
  END IF
  dalpha = ABS(alpha - self%alpha_last)
  !--- Freeze alpha on hard tol OR 2 consecutive steps above running minimum; also freeze bootstrap.
  IF(dalpha < self%dalpha_tol) THEN
    self%freeze_alpha = .TRUE.
    self%freeze_j_BS = .TRUE.
  ELSE IF(dalpha >= self%dalpha_min) THEN
    self%dalpha_no_improve = self%dalpha_no_improve + 1
    IF(self%dalpha_no_improve >= 2) THEN
      self%freeze_alpha = .TRUE.
      self%freeze_j_BS = .TRUE.
      IF(dalpha > self%dalpha_tol) THEN
        WRITE(char_buf,'(A,ES12.4,A,ES12.4,A)') &
          'Alpha convergence stalled,' // &
          ' relative change per nonlinear step = ', dalpha, &
          ', above recommended tolerance (dalpha_tol=', self%dalpha_tol, ')'
        CALL oft_warn(TRIM(char_buf))
      END IF
    END IF
  ELSE
    self%dalpha_no_improve = 0
    self%dalpha_min = dalpha
  END IF
END IF
self%alpha_last = alpha
!--- 6. Assemble jphi_total and compute F*F' knots.
jphi_total = alpha * jphi_ind + j_BS
IF(.NOT.ASSOCIATED(self%jphi_total_last)) ALLOCATE(self%jphi_total_last(self%npsi))
self%jphi_total_last = jphi_total
!--- Pressure scale
IF(ASSOCIATED(gseq%P_ani)) &
  CALL oft_abort('Jphi profiles do not support anisotropic pressure', &
                 'jphi_bs_update',__FILE__)
pscale = gseq%P%f(gseq%plasma_bounds(2))
pscale = gseq%pax_target / pscale
!--- Compute updated F*F' profile
CALL spline_eval(R_spline, 0.d0, 0) ! LCFS point for y0 calculation
pprime = gseq%P%fp(gseq%plasma_bounds(1))
self%y0 = 2.d0*(self%j0*alpha - R_spline%f(1)*pprime*pscale)/R_spline%f(2)
DO i = 1, self%npsi
  CALL spline_eval(R_spline, self%x(i), 0)
  pprime = gseq%P%fp(self%x(i)*(gseq%plasma_bounds(2) - &
                                  gseq%plasma_bounds(1)) + &
                                  gseq%plasma_bounds(1))
  self%yp(i) = 2.d0*(jphi_total(i) - R_spline%f(1)*pprime*pscale)/R_spline%f(2)
END DO
CALL spline_dealloc(R_spline)
! Disable Ip matching; scale is enforced above.
IF(gseq%Itor_target>0.d0)gseq%Itor_target=-gseq%Itor_target
gseq%ffp_scale=1.d0
!--- 7. Diagnostics.
IF(gseq%boot_ops%diagnose_bs)THEN
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] ip_target   = ', ip_target
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] ip_result_lo= ', ip_result_lo
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] ip_result_hi= ', ip_result_hi
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] alpha       = ', alpha
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] dalpha      = ', dalpha
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] djBS        = ', djBS
  WRITE(*,'(A,L1)')     '  [jphi_bs_update] bs_frozen   = ', self%freeze_j_BS
  WRITE(*,'(A,L1)')     '  [jphi_bs_update] freeze_alpha= ', self%freeze_alpha
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] j_BS max    = ', MAXVAL(ABS(j_BS))
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] jphi max    = ', MAXVAL(ABS(self%jphi))
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] jphi_rescale= ', jphi_rescale
  !--- Side-by-side Ip comparison: FEM nonlinear solve vs profile flux integral
  CALL gs_itor_nl(gseq, itor_nl)
  CALL gs_flux_int(gseq, self%x, jphi_total/qtmp, self%npsi, itor_flint)
  WRITE(*,'(A)') '  [jphi_bs_update] --- Ip comparison (current jphi_total) ---'
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] Ip(gs_itor_nl)    = ', itor_nl/mu0
  WRITE(*,'(A,ES12.4)') '  [jphi_bs_update] Ip(flux_int/qtmp) = ', itor_flint/mu0
END IF
!--- Clean up
DEALLOCATE(j_BS, jphi_total, jphi_ind, qtmp)
i=self%set_cofs(self%yp)
END SUBROUTINE jphi_bs_update
!------------------------------------------------------------------------------
!> Apply a smooth edge taper to an array, zeroing it at the plasma edge.
!>
!> @param n          Number of grid points
!> @param psi        Normalised psi grid
!> @param arr        Array to taper (modified in-place)
!> @param psi0       psi value (in the grid's own convention) where taper begins
!> @param shape      Taper shape: 1=cos² (Hann), 2=quintic smoothstep, 3=cubic power
!> @param oft_psi_conv  If .TRUE., psi uses OFT internal convention (0=LCFS, 1=axis)
!>                      and psi0 should be the OFT threshold (= 1 - standard psi0).
!>                      Taper then acts on points where psi <= psi0.
!>                      If .FALSE. (default), standard convention (0=axis, 1=LCFS),
!>                      taper acts on points where psi >= psi0.
!------------------------------------------------------------------------------
SUBROUTINE apply_edge_taper(n, psi, arr, psi0, shape, oft_psi_conv)
INTEGER(i4), INTENT(in)    :: n
REAL(r8),    INTENT(in)    :: psi(n)
REAL(r8),    INTENT(inout) :: arr(n)
REAL(r8),    INTENT(in)    :: psi0
INTEGER(i4), INTENT(in)    :: shape
LOGICAL,     OPTIONAL, INTENT(in) :: oft_psi_conv
!---
INTEGER(i4) :: i
LOGICAL     :: do_oft
REAL(r8)    :: t_taper, w_taper, span
CHARACTER(len=80) :: char_buf
REAL(r8), PARAMETER :: HALF_PI = 1.5707963267948966_r8
do_oft = .FALSE.
IF (PRESENT(oft_psi_conv)) do_oft = oft_psi_conv
! Taper width: distance from threshold to the plasma edge.
! OFT convention: edge at psi=0, so width = psi0.
! Standard convention: edge at psi=1, so width = 1 - psi0.
IF (do_oft) THEN
  span = psi0
ELSE
  span = 1.0_r8 - psi0
END IF
IF (span < 1.0e-6_r8) RETURN
DO i = 1, n
  IF (do_oft) THEN
    ! OFT convention: edge at psi=0, taper region is psi in [0, psi0]
    IF (psi(i) > psi0) CYCLE
    t_taper = (psi0 - psi(i)) / span
  ELSE
    ! Standard convention: edge at psi=1, taper region is psi in [psi0, 1]
    IF (psi(i) < psi0) CYCLE
    t_taper = (psi(i) - psi0) / span
  END IF
  t_taper = MIN(MAX(t_taper, 0.0_r8), 1.0_r8)
  SELECT CASE (shape)
    CASE (1)  ! cos² / Hann
      w_taper = COS(HALF_PI * t_taper)**2
    CASE (2)  ! quintic smoothstep
      w_taper = 1.0_r8 - t_taper**3*(6.0_r8*t_taper**2 - 15.0_r8*t_taper + 10.0_r8)
    CASE (3)  ! cubic power decay
      w_taper = (1.0_r8 - t_taper)**3
    CASE DEFAULT
      WRITE(char_buf,'(A,I0,A)') 'apply_edge_taper: unknown taper shape=', shape, '; no taper applied'
      CALL oft_warn(TRIM(char_buf))
      w_taper = 1.0_r8
  END SELECT
  arr(i) = arr(i) * w_taper
END DO
END SUBROUTINE apply_edge_taper
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE gs_flux_int(self,psi_tmp,field_tmp,nvals,result)
class(gs_equil), INTENT(inout) :: self !< Pointer to TokaMaker object
REAL(r8), INTENT(in) :: psi_tmp(:) !< Needs docs
REAL(r8), INTENT(in) :: field_tmp(:) !< Needs docs
INTEGER(i4), INTENT(in) :: nvals
REAL(r8), INTENT(out) :: result !< Needs docs
INTEGER(i4) :: i,m
REAL(r8) :: area,psitmp(1),sgop(3,3)
TYPE(gs_prof_interp) :: prof_interp_obj
! DEBUG_STACK_PUSH
!---Setup
result=0.d0
prof_interp_obj%mode=4
CALL prof_interp_obj%setup(self)
!$omp parallel do private(m,psitmp,sgop,area) reduction(+:result)
do i=1,self%device%mesh%nc
  IF(self%device%mesh%reg(i)/=1)CYCLE
  !---Loop over quadrature points
  do m=1,self%device%fe_rep%quad%np
    call self%device%mesh%jacobian(i,self%device%fe_rep%quad%pts(:,m),sgop,area)
    call prof_interp_obj%interp(i,self%device%fe_rep%quad%pts(:,m),sgop,psitmp)
    psitmp(1)=linterp(psi_tmp,field_tmp,nvals,psitmp(1),0)
    IF(psitmp(1)>-1.d98)result = result + psitmp(1)*area*self%device%fe_rep%quad%wts(m)
  end do
end do
!---Global reduction and cleanup
result=oft_mpi_sum(result)
CALL prof_interp_obj%delete()
! DEBUG_STACK_POP
END SUBROUTINE gs_flux_int
!------------------------------------------------------------------------------
!> Computes the bootstrap current on a uniform psi_N grid.
!>
!> Translated from Python bootstrap.py: calculate_profiles_and_bootstrap
!> (up to and including the redl_bootstrap call).
!>
!> Fixed conventions: NRL Coulomb logs, Koh ion collisionality model,
!> Redl 2021 jboot1 form with use_sign_q=.TRUE., L34=L31.
!>
!> @param gseq    Equilibrium object (must have Te, Ti, ne, ni, Zeff set)
!> @param n_psi    Number of flux surface samples
!> @param psi_N    Normalised psi grid [0,1], arbitrary spacing
!> @param j_BS Output: average toroidal bootstrap current density [A/m^2] on psi_N grid
!------------------------------------------------------------------------------
SUBROUTINE calculate_bootstrap(gseq, n_psi, psi_N, j_BS, &
                               isolate_edge_jBS, parameterize_jBS, scale_jBS, j_spike)
CLASS(gs_equil), INTENT(inout) :: gseq
INTEGER(i4), INTENT(in) :: n_psi
REAL(r8), INTENT(in) :: psi_N(n_psi)  !< Normalised psi grid [0,1], arbitrary spacing
REAL(r8), INTENT(out) :: j_BS(n_psi)
LOGICAL,  OPTIONAL, INTENT(in)  :: isolate_edge_jBS  !< If .TRUE., isolate edge spike from core
LOGICAL,  OPTIONAL, INTENT(in)  :: parameterize_jBS  !< If .TRUE., use parametrised skew-normal fit
REAL(r8), OPTIONAL, INTENT(in)  :: scale_jBS         !< Scaling factor applied to spike profile (default 1)
REAL(r8), OPTIONAL, INTENT(out) :: j_spike(n_psi)    !< Processed spike profile [A/m^2]
!---
INTEGER(i4) :: i
REAL(r8) :: psi_abs(n_psi)
REAL(r8) :: Te(n_psi), Ti(n_psi), ne(n_psi), ni(n_psi), Zeff(n_psi)
REAL(r8) :: pe(n_psi), pi_arr(n_psi)
REAL(r8) :: f(n_psi)        !< I(psi) = R*Bt [T*m]
REAL(r8) :: fc(n_psi)       !< Circulating fraction
REAL(r8) :: ft(n_psi)       !< Trapped particle fraction
REAL(r8) :: eps(n_psi)      !< Inverse aspect ratio
REAL(r8) :: qvals(n_psi)    !< Safety factor
REAL(r8) :: R_avg(n_psi)    !< <R> [m] per flux surface
REAL(r8) :: r_avgs_saut(n_psi,3)    !< Sauter FSA: <R>, <1/R>, <a>
REAL(r8) :: modb_avgs_saut(n_psi,2) !< Sauter FSA: <|B|>, <|B|^2>
REAL(r8) :: dn_e_dpsi(n_psi), dT_e_dpsi(n_psi)
REAL(r8) :: dn_i_dpsi(n_psi), dT_i_dpsi(n_psi)
REAL(r8) :: ln_le(n_psi), ln_lii(n_psi), Z_lnLam(n_psi)
REAL(r8) :: Zavg(n_psi), Zion(n_psi)
REAL(r8) :: nu_i_star(n_psi), nu_e_star(n_psi)
REAL(r8) :: B_times_Jbs(n_psi)
REAL(r8) :: psi_range, Zdom
REAL(r8), PARAMETER :: EC = 1.602176634e-19_r8
! Locals for optional edge-spike isolation
LOGICAL  :: do_isolate, do_parametrize
REAL(r8) :: scl_jBS
REAL(r8), ALLOCATABLE :: masked_spike(:), parameterized_spike(:)
! Workspace for psi_N convention flip (OFT: 0=LCFS,1=axis → standard: 0=axis,1=LCFS)
REAL(r8), ALLOCATABLE :: psi_N_std(:), j_BS_std(:), mask_std(:), param_std(:)
INTEGER(i4) :: k_flip
!---
! Compute raw psi values from the caller-supplied psi_N grid
psi_range = gseq%plasma_bounds(2) - gseq%plasma_bounds(1)
psi_abs = gseq%plasma_bounds(1) + psi_N * psi_range
! Evaluate kinetic profiles on the psi grid; convert Te/Ti from keV to eV
DO i = 1, n_psi
  Te(i) = gseq%Te%fp(psi_N(i)) * 1000.0_r8
  Ti(i) = gseq%Ti%fp(psi_N(i)) * 1000.0_r8
  ne(i) = gseq%ne%fp(psi_N(i))
  ni(i) = gseq%ni%fp(psi_N(i))
END DO
IF(MAXVAL(Te) > 5.0e5_r8)CALL oft_warn('calculate_bootstrap: max(Te) > 500 keV — profiles should be in keV, not eV')
IF(MAXVAL(Ti) > 5.0e5_r8)CALL oft_warn('calculate_bootstrap: max(Ti) > 500 keV — profiles should be in keV, not eV')
Zeff = gseq%boot_ops%Zeff
IF(gseq%boot_ops%Zeff <= 0.0_r8)CALL oft_abort('Zeff must be set explicitly via set_boot_ops before bootstrap calculation','calculate_bootstrap',__FILE__)
! Get I(psi) = R*Bt = F(psi) profile on the psi_N grid
DO i = 1, n_psi
  IF(gseq%mode==0)THEN
    f(i) = gseq%ffp_scale * gseq%I%f(psi_abs(i)) + gseq%I%f_offset
  ELSE
    f(i) = SQRT(gseq%ffp_scale * gseq%I%f(psi_abs(i)) + gseq%I%f_offset**2) &
           + gseq%I%f_offset * (1.0_r8 - SIGN(1.0_r8, gseq%I%f_offset))
  END IF
END DO
! Get flux-surface geometry: fc, eps, q, and <R> in one tracing pass
CALL sauter_fc(gseq, n_psi, psi_N, fc, r_avgs_saut, modb_avgs_saut, qprof=qvals)
R_avg = r_avgs_saut(:,1)
! ===================================================================
ft = 1.0_r8 - fc
! Pressures [Pa]
pe = EC * ne * Te
pi_arr = EC * ni * Ti
! Gradients d/dpsi [Wb^-1] via non-uniform finite differences
! (one-sided at endpoints, central at interior points)
dT_e_dpsi(1)     = (Te(2)     - Te(1))       / (psi_abs(2)     - psi_abs(1))
dT_e_dpsi(n_psi) = (Te(n_psi) - Te(n_psi-1)) / (psi_abs(n_psi) - psi_abs(n_psi-1))
dT_i_dpsi(1)     = (Ti(2)     - Ti(1))       / (psi_abs(2)     - psi_abs(1))
dT_i_dpsi(n_psi) = (Ti(n_psi) - Ti(n_psi-1)) / (psi_abs(n_psi) - psi_abs(n_psi-1))
dn_e_dpsi(1)     = (ne(2)     - ne(1))       / (psi_abs(2)     - psi_abs(1))
dn_e_dpsi(n_psi) = (ne(n_psi) - ne(n_psi-1)) / (psi_abs(n_psi) - psi_abs(n_psi-1))
dn_i_dpsi(1)     = (ni(2)     - ni(1))       / (psi_abs(2)     - psi_abs(1))
dn_i_dpsi(n_psi) = (ni(n_psi) - ni(n_psi-1)) / (psi_abs(n_psi) - psi_abs(n_psi-1))
DO i = 2, n_psi-1
  dT_e_dpsi(i) = (Te(i+1) - Te(i-1)) / (psi_abs(i+1) - psi_abs(i-1))
  dT_i_dpsi(i) = (Ti(i+1) - Ti(i-1)) / (psi_abs(i+1) - psi_abs(i-1))
  dn_e_dpsi(i) = (ne(i+1) - ne(i-1)) / (psi_abs(i+1) - psi_abs(i-1))
  dn_i_dpsi(i) = (ni(i+1) - ni(i-1)) / (psi_abs(i+1) - psi_abs(i-1))
END DO
! In the Fortran internal convention psi increases LCFS→axis (opposite to the
! standard Sauter/Redl derivation where psi increases axis→LCFS).  Negate
! all gradients so the bootstrap formula sees the conventional sign.
dT_e_dpsi = -dT_e_dpsi
dT_i_dpsi = -dT_i_dpsi
dn_e_dpsi = -dn_e_dpsi
dn_i_dpsi = -dn_i_dpsi
! Coulomb logarithms (NRL formulary)
! Electron: ne divided by 1e6 to convert m^-3 -> cm^-3
ln_le = 23.5_r8 &
      - LOG(SQRT(ne/1.0e6_r8) * Te**(-1.25_r8)) &
      - SQRT(1.0e-5_r8 + (LOG(Te) - 2.0_r8)**2 / 16.0_r8)
! Ion: constant 30 absorbs the m^-3 -> cm^-3 conversion
Z_lnLam = MAX(ne/ni, 1.0_r8)
ln_lii = 30.0_r8 - LOG(Z_lnLam**3 * SQRT(ni) / Ti**1.5_r8)
ln_le  = MAX(ln_le,  10.0_r8)
ln_lii = MAX(ln_lii, 10.0_r8)
! Ion collisionality: Koh multi-species model
Zdom = 1.0_r8    ! dominant ion charge (deuterium)
Zavg = ne / ni
Zion = (Zdom**2 * Zavg * Zeff)**0.25_r8
nu_i_star = 4.90e-18_r8 * ABS(qvals) * R_avg * ni &
          * Zion**4 * ln_lii / (Ti**2 * eps**1.5_r8)
! Electron collisionality
nu_e_star = 6.921e-18_r8 * ABS(qvals) * R_avg * ne &
          * Zeff * ln_le / (Te**2 * eps**1.5_r8)
! Compute bootstrap current via Redl 2021
CALL redl_bootstrap(n_psi, Te, Ti, ne, ni, pe, pi_arr, Zeff, qvals, eps, ft, f, &
    dT_e_dpsi, dT_i_dpsi, dn_e_dpsi, dn_i_dpsi, &
    ln_le, ln_lii, nu_e_star, nu_i_star, B_times_Jbs)
! Convert parallel bootstrap to phi component: j_phi = B_times_Jbs * <R> / F
WHERE(ABS(f) > 0.0_r8)
  j_BS = B_times_Jbs * (R_avg / f)
ELSEWHERE
  j_BS = 0.0_r8
END WHERE
! Guard NaN (where F -> 0)
WHERE(.NOT.(ABS(j_BS) < 1.0e99_r8)) j_BS = 0.0_r8
IF(gseq%boot_ops%diagnose_bs)THEN
  WRITE(*,'(A)') '  [calculate_bootstrap] geometry & collisionality sample (i=1,mid,n):'
  WRITE(*,'(A,3ES12.4)') '    <R>      : ', r_avgs_saut(1,1), r_avgs_saut(n_psi/2,1), r_avgs_saut(n_psi,1)
  WRITE(*,'(A,3ES12.4)') '    <1/R>    : ', r_avgs_saut(1,2), r_avgs_saut(n_psi/2,2), r_avgs_saut(n_psi,2)
  WRITE(*,'(A,3ES12.4)') '    eps      : ', eps(1), eps(n_psi/2), eps(n_psi)
  WRITE(*,'(A,3ES12.4)') '    q        : ', qvals(1), qvals(n_psi/2), qvals(n_psi)
  WRITE(*,'(A,3ES12.4)') '    nu_e_star: ', nu_e_star(1), nu_e_star(n_psi/2), nu_e_star(n_psi)
  WRITE(*,'(A,3ES12.4)') '    j_BS     : ', j_BS(1), j_BS(n_psi/2), j_BS(n_psi)
END IF
END SUBROUTINE calculate_bootstrap
!------------------------------------------------------------------------------
!> Redl 2021 bootstrap current formula.
!>
!> Translates Python bootstrap.py:redl_bootstrap with fixed settings:
!>   - use_legacy_L34 = .FALSE.  (L34 = L31)
!>   - use_sign_q     = .TRUE.
!>   - formula_form   = 'jboot1'
!>   - nu_e_star and nu_i_star are passed in pre-computed (no internal fallback)
!>
!> Reference: Redl et al., Phys. Plasmas 28, 022502 (2021)
!>
!> @param n           Number of flux surfaces
!> @param Te          Electron temperature [eV]
!> @param Ti          Ion temperature [eV]
!> @param ne          Electron density [m^-3]
!> @param ni          Ion density [m^-3]
!> @param pe          Electron pressure [Pa]
!> @param pi          Ion pressure [Pa]
!> @param Zeff        Effective charge
!> @param q           Safety factor
!> @param eps         Inverse aspect ratio
!> @param fT          Trapped particle fraction
!> @param I_psi       Toroidal current function F = R*Bt [T*m]
!> @param dT_e_dpsi   d(Te)/d(psi) [eV/Wb]
!> @param dT_i_dpsi   d(Ti)/d(psi) [eV/Wb]
!> @param dn_e_dpsi   d(ne)/d(psi) [m^-3/Wb]
!> @param dn_i_dpsi   d(ni)/d(psi) [m^-3/Wb]
!> @param ln_lambda_e Electron Coulomb logarithm
!> @param ln_lambda_ii Ion Coulomb logarithm
!> @param nu_e_star   Electron collisionality (pre-computed)
!> @param nu_i_star   Ion collisionality (pre-computed)
!> @param avg_j_bootstrap_times_B Output: <j_{bs,parallel}*B> [same units as -I_psi * pe * d/dpsi]
!------------------------------------------------------------------------------
SUBROUTINE redl_bootstrap(n, Te, Ti, ne, ni, pe, pi, Zeff, q, eps, fT, I_psi, &
    dT_e_dpsi, dT_i_dpsi, dn_e_dpsi, dn_i_dpsi, &
    ln_lambda_e, ln_lambda_ii, nu_e_star, nu_i_star, avg_j_bootstrap_times_B)
INTEGER(i4), INTENT(in) :: n
REAL(r8), INTENT(in) :: Te(n), Ti(n), ne(n), ni(n)
REAL(r8), INTENT(in) :: pe(n), pi(n), Zeff(n)
REAL(r8), INTENT(in) :: q(n), eps(n), fT(n), I_psi(n)
REAL(r8), INTENT(in) :: dT_e_dpsi(n), dT_i_dpsi(n)
REAL(r8), INTENT(in) :: dn_e_dpsi(n), dn_i_dpsi(n)
REAL(r8), INTENT(in) :: ln_lambda_e(n), ln_lambda_ii(n)
REAL(r8), INTENT(in) :: nu_e_star(n), nu_i_star(n)
REAL(r8), INTENT(out) :: avg_j_bootstrap_times_B(n)
!---
REAL(r8) :: R_pe(n)
REAL(r8) :: ft_31_d1(n), ft_31_d2(n), X31(n), L31(n), L34(n), dZ(n)
REAL(r8) :: dee_2(n), dee_3(n), X32_ee(n), F32_ee(n)
REAL(r8) :: dei_2(n), dei_3(n), X32_ei(n), F32_ei(n), L32(n)
REAL(r8) :: alpha0(n), alpha(n)
REAL(r8) :: dp_dpsi(n), bra1(n), bra2(n), bra3(n)
REAL(r8), PARAMETER :: EC = 1.602176634e-19_r8
!---
R_pe = pe / (pe + pi)
! =====================================================================
! L31 (Redl Eqs. 10-11)
! =====================================================================
ft_31_d1 = (0.67_r8 * (1.0_r8 - 0.7_r8*fT) * SQRT(nu_e_star)) &
         / (0.56_r8 + 0.44_r8*Zeff)
ft_31_d2 = ((0.52_r8 + 0.086_r8*SQRT(nu_e_star)) &
          * (1.0_r8 + 0.87_r8*fT) * nu_e_star) &
         / (1.0_r8 + 1.13_r8*SQRT(MAX(Zeff - 1.0_r8, 0.0_r8)))
X31 = fT / (1.0_r8 + ft_31_d1 + ft_31_d2)
dZ = Zeff**1.2_r8 - 0.71_r8
L31 = (1.0_r8 + 0.15_r8/dZ)*X31 &
    - (0.22_r8/dZ)*X31**2 &
    + (0.01_r8/dZ)*X31**3 &
    + (0.06_r8/dZ)*X31**4
! L34 = L31 (use_legacy_L34 = .FALSE.)
L34 = L31
! =====================================================================
! L32 (Redl Eqs. 12-16)
! =====================================================================
dee_2 = (0.23_r8 * (1.0_r8 - 0.96_r8*fT) * SQRT(nu_e_star)) &
      / SQRT(Zeff)
dee_3 = (0.13_r8 * (1.0_r8 - 0.38_r8*fT) * nu_e_star / (Zeff**2)) &
      * (SQRT(1.0_r8 + 2.0_r8*SQRT(MAX(Zeff - 1.0_r8, 0.0_r8))) &
       + fT**2 * SQRT((0.075_r8 + 0.25_r8*(Zeff - 1.0_r8)**2) * nu_e_star))
X32_ee = fT / (1.0_r8 + dee_2 + dee_3)
F32_ee = ( (0.1_r8 + 0.6_r8*Zeff) &
         / (Zeff*(0.77_r8 + 0.63_r8*(1.0_r8 + (Zeff - 1.0_r8)**1.1_r8))) &
         * (X32_ee - X32_ee**4) &
         + 0.7_r8/(1.0_r8 + 0.2_r8*Zeff) &
         * (X32_ee**2 - X32_ee**4 - 1.2_r8*(X32_ee**3 - X32_ee**4)) &
         + 1.3_r8/(1.0_r8 + 0.5_r8*Zeff) * X32_ee**4 )
dei_2 = (0.87_r8 * (1.0_r8 + 0.39_r8*fT) * SQRT(nu_e_star)) &
      / (1.0_r8 + 2.95_r8*(Zeff - 1.0_r8)**2)
dei_3 = 1.53_r8 * (1.0_r8 - 0.37_r8*fT) * nu_e_star &
      * (2.0_r8 + 0.375_r8*(Zeff - 1.0_r8))
X32_ei = fT / (1.0_r8 + dei_2 + dei_3)
F32_ei = ( -(0.4_r8 + 1.93_r8*Zeff) / (Zeff*(0.8_r8 + 0.6_r8*Zeff)) &
           * (X32_ei - X32_ei**4) &
           + 5.5_r8/(1.5_r8 + 2.0_r8*Zeff) &
           * (X32_ei**2 - X32_ei**4 - 0.8_r8*(X32_ei**3 - X32_ei**4)) &
           - 1.3_r8/(1.0_r8 + 0.5_r8*Zeff) * X32_ei**4 )
L32 = F32_ee + F32_ei
! =====================================================================
! Alpha (Redl Eqs. 20-21)
! =====================================================================
alpha0 = -(0.62_r8 + 0.055_r8*(Zeff - 1.0_r8)) &
        / (0.53_r8 + 0.17_r8*(Zeff - 1.0_r8)) &
        * (1.0_r8 - fT) &
        / (1.0_r8 - (0.31_r8 - 0.065_r8*(Zeff - 1.0_r8))*fT - 0.25_r8*fT**2)
alpha = ((alpha0 + 0.7_r8*Zeff*SQRT(fT)*SQRT(nu_i_star)) &
       / (1.0_r8 + 0.18_r8*SQRT(nu_i_star)) &
       - 0.002_r8*nu_i_star**2*fT**6) &
      / (1.0_r8 + 0.004_r8*nu_i_star**2*fT**6)
! =====================================================================
! Assemble: jboot1 form with use_sign_q = .TRUE.
! =====================================================================
dp_dpsi = (ne*dT_e_dpsi + Te*dn_e_dpsi &
         + ni*dT_i_dpsi + Ti*dn_i_dpsi) * EC
bra1 = L31 * dp_dpsi / pe
bra2 = L32 * dT_e_dpsi / Te
bra3 = L34 * alpha * (1.0_r8 - R_pe) / R_pe * dT_i_dpsi / Ti
avg_j_bootstrap_times_B = -I_psi * pe * (bra1 + bra2 + bra3)
! Apply sign(q)
avg_j_bootstrap_times_B = avg_j_bootstrap_times_B * SIGN(1.0_r8, q)
END SUBROUTINE redl_bootstrap
!---------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE create_dipole_b0_prof(func,npsi)
CLASS(flux_func), POINTER, INTENT(out) :: func
INTEGER(4), INTENT(in) :: npsi
INTEGER(4) :: i
ALLOCATE(dipole_b0_flux_func::func)
select type(self=>func)
  type is(dipole_b0_flux_func)
    !---
    self%npsi=npsi
    CALL spline_alloc(self%func,self%npsi,1)
    ALLOCATE(self%fun_loc(omp_get_max_threads()))
    DO i=1,omp_get_max_threads()
      CALL spline_alloc(self%fun_loc(i),self%npsi,1)
    END DO
end select
END SUBROUTINE create_dipole_b0_prof
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
CLASS(flux_func), POINTER, INTENT(out) :: func
INTEGER(4), INTENT(in) :: npsi
INTEGER(4) :: i
ALLOCATE(mirror_b0_flux_func::func)
select type(self=>func)
  type is(mirror_b0_flux_func)
    !---
    self%npsi=npsi
    CALL spline_alloc(self%func,self%npsi,1)
    ALLOCATE(self%fun_loc(omp_get_max_threads()))
    DO i=1,omp_get_max_threads()
      CALL spline_alloc(self%fun_loc(i),self%npsi,1)
    END DO
end select
END SUBROUTINE create_mirror_b0_prof
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
