!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file mercier.F90
!
!> Mercier criterion and DCON interface
!!
!! @authors Chris Hansen
!! @date August 2014
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
module oft_gs_mercier
use oft_base
use oft_mesh_type, only: oft_bmesh, bmesh_findcell
USE oft_la_base, ONLY: oft_vector
USE fem_utils, ONLY: bfem_interp
use oft_lag_basis, only: oft_blag_d2eval, oft_blag_geval
USE oft_blag_operators, only: oft_lag_brinterp
use oft_gs, only: gs_eq, flux_func, gs_psi2r, gs_itor_nl, oft_indent, &
  oft_increase_indent, oft_decrease_indent, gsinv_interp
use tracing_2d, only: set_tracer, active_tracer, tracinginv_fs
use oft_gs_profiles, only: spline_flux_func
use spline_mod
USE mhd_utils, ONLY: mu0
implicit none
!------------------------------------------------------------------------------
! CLASS mercierinv_interp
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
! INTERFACE mercier_flux_func
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
contains
!------------------------------------------------------------------------------
! SUBROUTINE create_mercier_ff
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
! SUBROUTINE mercier_update
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
! ! SUBROUTINE minterpinv_setup
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
! SUBROUTINE minterpinv_apply
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
end module oft_gs_mercier
