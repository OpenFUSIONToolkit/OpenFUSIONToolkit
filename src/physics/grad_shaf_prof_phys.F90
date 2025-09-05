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
USE oft_blag_operators, only: oft_lag_brinterp
use oft_gs, only: gs_eq, flux_func, gs_psi2r, gs_itor_nl, oft_indent, &
  oft_increase_indent, oft_decrease_indent, gsinv_interp, gs_prof_interp, &
  gs_get_qprof
use tracing_2d, only: set_tracer, active_tracer, tracinginv_fs
use oft_gs_profiles, only: spline_flux_func, linterp_flux_func
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
    procedure :: update => mercier_update
end type mercier_flux_func
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(linterp_flux_func) :: jphi_flux_func
  integer(4) :: ngeom = 50
  real(8), pointer, dimension(:) :: jphi => NULL() !< Needs docs
contains
  !> Needs docs
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
  integer(4) :: ntheta = 128
  real(r8) :: psi_pad = 1.d-1
  TYPE(spline_type) :: funcp
contains
    procedure :: update => dipole_b0_update
end type dipole_b0_flux_func
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
    DO i=1,omp_get_max_threads()
      CALL spline_alloc(self%fun_loc(i),self%npsi,1)
    END DO
end select
END SUBROUTINE create_mercier_ff
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine mercier_update(self,gseq)
class(mercier_flux_func), intent(inout) :: self
class(gs_eq), intent(inout) :: gseq
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
CALL psi_int%setup(gseq%fe_rep)
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
  pt=[(gseq%rmax-raxis)*j/REAL(100,8)+raxis,zaxis,0.d0]
  CALL bmesh_findcell(gseq%mesh,cell,pt,f)
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
CALL field%setup(gseq%fe_rep)
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
    call tracinginv_fs(gseq%mesh,pt(1:2))
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
    call bmesh_findcell(gseq%mesh,active_tracer%cell,pt,active_tracer%f)
    call gseq%mesh%jacobian(active_tracer%cell,active_tracer%f,gop,vol)
    call psi_int%interp(active_tracer%cell,active_tracer%f,gop,psi_surf)
    !---Get flux variables
    I=gseq%alam*gseq%I%f(psi_surf(1))+gseq%I%f_offset
    Ip=gseq%alam*gseq%I%fp(psi_surf(1))
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
SUBROUTINE create_jphi_ff(func,npsi,psivals,yvals)
CLASS(flux_func), POINTER, INTENT(out) :: func
INTEGER(4), INTENT(in) :: npsi
REAL(8), INTENT(in) :: psivals(npsi)
REAL(8), INTENT(in) :: yvals(npsi)
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
  self%y0=0.d0
  DO i=1,self%npsi
    self%x(i) = psivals(i)
    self%jphi(i) = yvals(i)
    self%yp(i) = psivals(i) ! Dummy initialization
  END DO
  self%yp = self%yp/(SUM(ABS(self%yp))/REAL(self%npsi,8)) ! Consistent (hopefully) normalization
  ierr=self%set_cofs(self%yp)
  IF(oft_debug_print(1))WRITE(*,*)'Jphi linear interpolator Created',self%ncofs,self%x,self%y0
END SELECT

END SUBROUTINE create_jphi_ff
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine jphi_update(self,gseq)
class(jphi_flux_func), intent(inout) :: self
class(gs_eq), intent(inout) :: gseq
INTEGER(i4) :: i
REAL(r8) :: jphi_norm,pscale,pprime
REAL(r8), ALLOCATABLE :: ravgs(:,:),qtmp(:),psi_q(:)
type(spline_type) :: R_spline
self%plasma_bounds=gseq%plasma_bounds
IF(gseq%mode/=1)CALL oft_abort("Jphi profile requires (F^2)' formulation","jphi_update",__FILE__)
! IF(gseq%Itor_target<0.d0)CALL oft_abort("Jphi profile requires Ip target","jphi_update",__FILE__)
IF(gseq%pax_target<0.d0)CALL oft_abort("Jphi profile requires Pax target","jphi_update",__FILE__)
!---Get updated flux surface geometry for Jphi -> F*F' mapping
ALLOCATE(ravgs(self%ngeom,3),psi_q(self%ngeom),qtmp(self%ngeom))
psi_q=[(REAL(i,8)/REAL(self%ngeom+1,8),i=1,self%ngeom)]
CALL gs_get_qprof(gseq,self%ngeom,psi_q,qtmp,ravgs=ravgs)
! WRITE(*,*)'  <R>     = ',ravgs(2,1),ravgs(self%ngeom-2,1)
! WRITE(*,*)'  <1/R>   = ',ravgs(2,2),ravgs(self%ngeom-2,2)
CALL spline_alloc(R_spline,self%ngeom-1,2)
R_spline%xs(0:self%ngeom-2)=psi_q; R_spline%xs(self%ngeom-1)=1.d0
R_spline%fs(0:self%ngeom-2,1)=ravgs(1:self%ngeom-1,1); R_spline%fs(self%ngeom-1,1)=gseq%o_point(1)
R_spline%fs(0:self%ngeom-2,2)=ravgs(1:self%ngeom-1,2); R_spline%fs(self%ngeom-1,2)=1.d0/gseq%o_point(1)
CALL spline_fit(R_spline,"extrap")
DEALLOCATE(ravgs,psi_q,qtmp)
!---Update jphi normalization to match Ip target
ALLOCATE(qtmp(self%npsi))
DO i=1,self%npsi
  CALL spline_eval(R_spline,self%x(i),0)
  qtmp(i) = R_spline%f(1)*R_spline%f(2)
END DO
CALL gs_flux_int(gseq,self%x,self%jphi/qtmp,self%npsi,jphi_norm)
jphi_norm=ABS(gseq%Itor_target)/jphi_norm
DEALLOCATE(qtmp)
!---Get pressure profile
CALL gseq%P%update(gseq) ! Make sure pressure profile is up to date with EQ
pscale=gseq%P%f(gseq%plasma_bounds(2))
pscale=gseq%pax_target/pscale
!---Compute updated F*F' profile ! 2.0*(jtor -  R_avg * (-pprime)) * (mu0 / one_over_R_avg)
DO i=1,self%npsi
  CALL spline_eval(R_spline,self%x(i),0)
  pprime=gseq%P%fp(self%x(i)*(gseq%plasma_bounds(2)-gseq%plasma_bounds(1))+gseq%plasma_bounds(1))
  self%yp(i) = 2.d0*(self%jphi(i)*jphi_norm - R_spline%f(1)*pprime*pscale)/R_spline%f(2)
END DO
! Disable Ip matching and fix F*F' scale (matching is done here instead)
IF(gseq%Itor_target>0.d0)gseq%Itor_target=-gseq%Itor_target
gseq%alam=1.d0
!---Clean up
CALL spline_dealloc(R_spline)
i=self%set_cofs(self%yp)
end subroutine jphi_update
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE gs_flux_int(self,psi_tmp,field_tmp,nvals,result)
class(gs_eq), INTENT(inout) :: self !< Pointer to TokaMaker object
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
do i=1,self%mesh%nc
  IF(self%mesh%reg(i)/=1)CYCLE
  !---Loop over quadrature points
  do m=1,self%fe_rep%quad%np
    call self%mesh%jacobian(i,self%fe_rep%quad%pts(:,m),sgop,area)
    call prof_interp_obj%interp(i,self%fe_rep%quad%pts(:,m),sgop,psitmp)
    psitmp(1)=linterp(psi_tmp,field_tmp,nvals,psitmp(1),0)
    IF(psitmp(1)>-1.d98)result = result + psitmp(1)*area*self%fe_rep%quad%wts(m)
  end do
end do
!---Global reduction and cleanup
result=oft_mpi_sum(result)
CALL prof_interp_obj%delete()
! DEBUG_STACK_POP
END SUBROUTINE gs_flux_int
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
    CALL spline_alloc(self%funcp,self%npsi,1)
    DO i=1,omp_get_max_threads()
      CALL spline_alloc(self%fun_loc(i),self%npsi,1)
    END DO
end select
END SUBROUTINE create_dipole_b0_prof
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine dipole_b0_update(self,gseq)
class(dipole_b0_flux_func), intent(inout) :: self
class(gs_eq), intent(inout) :: gseq
type(minbinv_interp), target :: field
type(oft_lag_brinterp) :: psi_int
real(8) :: I,Ip,q,qp,vp,vpp,s,a,b,pp,gop(3,3),psi_surf(1),vol,rmax
real(8) :: raxis,zaxis,f(3),pt(3),x1,x2,xr
real(8), pointer, dimension(:) :: v
real(8), parameter :: tol=1.d-10
integer(4) :: j,k,cell
!
IF(oft_debug_print(2))THEN
  WRITE(*,'(2A)')oft_indent,'Updating Mercier Pressure:'
  CALL oft_increase_indent
END IF
psi_int%u=>gseq%psi
CALL psi_int%setup(gseq%fe_rep)
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
  CALL bmesh_findcell(gseq%mesh,cell,pt,f)
  IF( (MAXVAL(f)>1.d0+tol) .OR. (MINVAL(f)<-tol) )EXIT
  CALL psi_int%interp(cell,f,gop,psi_surf)
  IF( psi_surf(1) > x1)EXIT
  rmax=pt(1)
END DO
IF(oft_debug_print(2))WRITE(*,'(2A,ES11.3)')oft_indent,'Rmin = ',rmax
!---Trace
call set_tracer(1)
!!$omp parallel private(field,gop,vol,psi_surf,I,Ip,v,q,qp,vp,vpp,s,a,b,pp,pt)
pt=[(.9d0*rmax+.1d0*raxis),zaxis,0.d0]
field%u=>gseq%psi
CALL field%setup(gseq%fe_rep)
active_tracer%neq=2
active_tracer%B=>field
active_tracer%maxsteps=8e4
active_tracer%tol=1.d-9
active_tracer%raxis=raxis
active_tracer%zaxis=zaxis
active_tracer%inv=.TRUE.
!!$omp do schedule(dynamic,1)
do j=1,self%npsi+1
  !------------------------------------------------------------------------------
  ! Trace contour
  !------------------------------------------------------------------------------
  !psi_surf(1)=(x2-x1)*(1.d0-j/REAL(self%npsi,4))**2
  psi_surf(1)=(x2-x1)*(1.d0-j/REAL(self%npsi+1,4))
  psi_surf(1)=x2 - psi_surf(1)
  !!$omp critical
  CALL gs_psi2r(gseq,psi_surf(1),pt)
  !!$omp end critical
  field%minB=1.d99
  call tracinginv_fs(gseq%mesh,pt(1:2))
  !---Exit if trace fails
  if(active_tracer%status/=1)THEN
    WRITE(*,*)'Tracer Error:',psi_surf(1),pt,active_tracer%y,active_tracer%status
    call oft_abort('Trace did not complete.','mercier_update',__FILE__)
  end if
  !---Get surface minB
  self%func%xs(j-1)=psi_surf(1)
  self%func%fs(j-1,1)=field%minB
  ! WRITE(*,*)psi_surf(1),field%minB
end do
!!$omp end parallel
! self%funcp%xs(0)=x1
! self%funcp%fs(0,1)=0.d0
!self%funcp%fs(0,1)=self%funcp%fs(1,1)
! self%funcp%xs(self%npsi)=x2!+.05d0*x2
! self%funcp%fs(self%npsi,1)=0.d0
self%yp1=self%func%fs(0,1)
self%ypn=self%func%fs(self%npsi,1)
!---Setup Spline
CALL spline_fit(self%func,"extrap")
! CALL spline_int(self%funcp)
!---
! self%func%xs=self%funcp%xs
! self%func%fs(:,1)=self%funcp%fsi(:,1)
! CALL spline_fit(self%func,"extrap")
DO k=1,omp_get_max_threads()
  CALL spline_copy(self%func,self%fun_loc(k))
END DO
!
self%xmin=self%func%xs(0)
self%xmax=self%func%xs(self%npsi)
IF(oft_debug_print(2))CALL oft_decrease_indent
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
end module grad_shaf_prof_phys
