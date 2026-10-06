!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file grad_shaf_torflux.F90
!
!> Flux functions defined on normalized toroidal flux
!!
!! Implements the \f$ \hat{\psi} \rightarrow \hat{\Phi} \f$ map (built from the equilibrium
!! q profile) and the transform applied by the non-overridable `f`, `fp`, `fpp`, and
!! `set_cofs` bindings of @ref oft_gs::flux_func.
!!
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
SUBMODULE (oft_gs) oft_gs_torflux
USE oft_gs_cutcell, ONLY: gs_qgeom_cutcell
IMPLICIT NONE
CONTAINS
!------------------------------------------------------------------------------
!> Get node locations in storage coordinate (none by default)
!------------------------------------------------------------------------------
MODULE PROCEDURE flux_func_get_nodes
ALLOCATE(nodes(0))
END PROCEDURE flux_func_get_nodes
!------------------------------------------------------------------------------
!> Evaluate \f$ 1-\hat{\Phi} \f$, \f$ J = d\hat{\Phi}/d\hat{\psi} \f$ and \f$ dJ/d\hat{\psi} \f$ at \f$ 1-\hat{\psi} \f$
!------------------------------------------------------------------------------
MODULE PROCEDURE tmap_eval
integer(i4) :: lo,hi
real(r8) :: h,u,pa,pb,ma,mb
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
END PROCEDURE tmap_eval
!------------------------------------------------------------------------------
!> Evaluate \f$ 1-\hat{\psi} \f$ at \f$ 1-\hat{\Phi} \f$ (inverse of @ref tmap_eval)
!------------------------------------------------------------------------------
MODULE PROCEDURE tmap_inv
real(r8) :: a,b,phi,jac
integer(i4) :: lo,i
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
END PROCEDURE tmap_inv
!------------------------------------------------------------------------------
!> Reset map to identity
!------------------------------------------------------------------------------
MODULE PROCEDURE tmap_delete
IF(ASSOCIATED(self%psis))DEALLOCATE(self%psis)
IF(ASSOCIATED(self%phihat))DEALLOCATE(self%phihat)
IF(ASSOCIATED(self%jac))DEALLOCATE(self%jac)
IF(ASSOCIATED(self%bin))DEALLOCATE(self%bin)
self%ns=0
self%Q=1.d0
self%stamp=self%stamp+1
END PROCEDURE tmap_delete
!------------------------------------------------------------------------------
!> Find interval `lo` such that `x(lo) <= val < x(lo+1)` (bisection, `x` ascending)
!------------------------------------------------------------------------------
pure function tmap_find(x,n,val) result(lo)
integer(i4), intent(in) :: n
real(r8), intent(in) :: x(n),val
integer(i4) :: lo,hi,k
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
type(gs_torflux_map), intent(inout) :: map
integer(i4) :: j,nbin
IF(ASSOCIATED(map%bin))DEALLOCATE(map%bin)
IF(map%ns<2)RETURN
nbin=4*map%ns
ALLOCATE(map%bin(nbin))
DO j=1,nbin
  map%bin(j)=tmap_find(map%psis,map%ns,REAL(j-1,r8)/REAL(nbin,r8))
END DO
end subroutine tmap_index
!------------------------------------------------------------------------------
!> Same result as `tmap_find` on `map%psis`, starting from the bin index
!------------------------------------------------------------------------------
function tmap_lookup(map,psihat) result(lo)
type(gs_torflux_map), intent(in) :: map
real(r8), intent(in) :: psihat
integer(i4) :: lo,nbin
IF(.NOT.ASSOCIATED(map%bin))THEN
  lo=tmap_find(map%psis,map%ns,psihat)
  RETURN
END IF
nbin=SIZE(map%bin)
lo=map%bin(MIN(INT(MIN(MAX(psihat,0.d0),1.d0)*REAL(nbin,r8))+1,nbin))
DO WHILE(lo<map%ns-1)
  IF(psihat<map%psis(lo+1))EXIT
  lo=lo+1
END DO
end function tmap_lookup
!------------------------------------------------------------------------------
!> Normalization offset and width (unset bounds: argument is already \f$ 1-\hat{\psi} \f$)
!------------------------------------------------------------------------------
pure subroutine flux_func_norm(self,x1,dpsi)
class(flux_func), intent(in) :: self
real(r8), intent(out) :: x1 !< Flux at LCFS
real(r8), intent(out) :: dpsi !< Axis minus LCFS flux
IF(self%plasma_bounds(1)<-1.d98)THEN
  x1=0.d0; dpsi=1.d0
ELSE
  x1=self%plasma_bounds(1); dpsi=self%plasma_bounds(2)-self%plasma_bounds(1)
END IF
end subroutine flux_func_norm
!------------------------------------------------------------------------------
!> Map \f$ \psi \f$ to storage coordinate, returns `.FALSE.` if no transform applies
!------------------------------------------------------------------------------
function flux_func_map(self,psi,psim,jac,djac) result(mapped)
class(flux_func), intent(in) :: self
real(r8), intent(in) :: psi !< Poloidal flux
real(r8), intent(out) :: psim !< Mapped flux \f$ \psi_a + (1-\hat{\Phi})(\psi_0-\psi_a) \f$
real(r8), intent(out) :: jac !< \f$ d\hat{\Phi}/d\hat{\psi} \f$
real(r8), optional, intent(out) :: djac !< \f$ d^2\hat{\Phi}/d\hat{\psi}^2 \f$
logical :: mapped
real(r8) :: x1,dpsi,psihat,phihat
mapped=.FALSE.
psim=psi; jac=1.d0
IF(PRESENT(djac))djac=0.d0
IF(self%coord==0)RETURN
IF(.NOT.ASSOCIATED(self%tmap))RETURN
CALL flux_func_norm(self,x1,dpsi)
psihat=(psi-x1)/dpsi
IF(psihat<0.d0)RETURN
CALL self%tmap%eval(psihat,phihat,jac,djac)
psim=x1+phihat*dpsi
mapped=.TRUE.
end function flux_func_map
!------------------------------------------------------------------------------
!> Evaluate flux function at \f$ \psi \f$
!------------------------------------------------------------------------------
MODULE PROCEDURE flux_func_f
real(r8) :: psim,jac
IF(.NOT.flux_func_map(self,psi,psim,jac))THEN
  b=self%f_native(psi)
ELSE IF(self%coord==1)THEN
  b=self%f_native(psim)
ELSE
  b=flux_func_relabel_f(self,psi)
END IF
END PROCEDURE flux_func_f
!------------------------------------------------------------------------------
!> Evaluate first derivative of flux function at \f$ \psi \f$
!------------------------------------------------------------------------------
MODULE PROCEDURE flux_func_fp
real(r8) :: psim,jac
IF(.NOT.flux_func_map(self,psi,psim,jac))THEN
  b=self%fp_native(psi)
ELSE IF(self%coord==1)THEN
  b=self%fp_native(psim)*jac
ELSE
  b=self%fp_native(psim)
END IF
END PROCEDURE flux_func_fp
!------------------------------------------------------------------------------
!> Evaluate second derivative of flux function at \f$ \psi \f$
!------------------------------------------------------------------------------
MODULE PROCEDURE flux_func_fpp
real(r8) :: psim,jac,djac,x1,dpsi
IF(.NOT.flux_func_map(self,psi,psim,jac,djac))THEN
  b=self%fpp_native(psi)
ELSE IF(self%coord==1)THEN
  CALL flux_func_norm(self,x1,dpsi)
  b=self%fpp_native(psim)*jac**2 &
    + self%fp_native(psim)*djac/dpsi
ELSE
  b=self%fpp_native(psim)*jac
END IF
END PROCEDURE flux_func_fpp
!------------------------------------------------------------------------------
!> Set coefficients and refresh derived data
!------------------------------------------------------------------------------
MODULE PROCEDURE flux_func_set_cofs
ierr=self%set_cofs_native(c)
IF(self%coord==2)CALL self%build_fint()
END PROCEDURE flux_func_set_cofs
!------------------------------------------------------------------------------
!> Integral of `fp_native` over \f$ 1-\hat{\psi} \in [a,b] \f$ (4-point Gauss)
!------------------------------------------------------------------------------
function flux_func_relabel_int(self,a,b) result(val)
class(flux_func), intent(inout) :: self
real(r8), intent(in) :: a,b
real(r8) :: val,s,phihat,jac,x1,dpsi
real(r8), parameter :: xg(4)=[-0.861136311594053d0,-0.339981043584856d0, &
  0.339981043584856d0,0.861136311594053d0]
real(r8), parameter :: wg(4)=[0.347854845137454d0,0.652145154862546d0, &
  0.652145154862546d0,0.347854845137454d0]
integer(i4) :: i
CALL flux_func_norm(self,x1,dpsi)
val=0.d0
DO i=1,4
  s=(a+b)/2.d0+xg(i)*(b-a)/2.d0
  CALL self%tmap%eval(s,phihat,jac)
  val=val+wg(i)*self%fp_native(x1+phihat*dpsi)
END DO
val=val*(b-a)/2.d0
end function flux_func_relabel_int
!------------------------------------------------------------------------------
!> Rebuild cumulative integral table on map surfaces (`coord=2` only)
!------------------------------------------------------------------------------
MODULE PROCEDURE flux_func_build_fint
integer(i4) :: k
IF(ASSOCIATED(self%fint))DEALLOCATE(self%fint)
IF(ASSOCIATED(self%fcof))DEALLOCATE(self%fcof)
self%fint_stamp=-1
IF(self%coord/=2)RETURN
IF(.NOT.ASSOCIATED(self%tmap))RETURN
IF(self%tmap%ns<2)RETURN
ALLOCATE(self%fint(self%tmap%ns),self%fcof(4,self%tmap%ns-1))
self%fint(1)=flux_func_relabel_int(self,0.d0,self%tmap%psis(1))
DO k=2,self%tmap%ns
  CALL flux_func_relabel_cof(self,self%tmap%psis(k-1),self%tmap%psis(k),self%fcof(:,k-1))
  self%fint(k)=self%fint(k-1)+SUM(self%fcof(:,k-1))
END DO
self%fint_stamp=self%tmap%stamp
END PROCEDURE flux_func_build_fint
!------------------------------------------------------------------------------
!> Quartic \f$ \sum_j c_j u^j \f$ (u in [0,1]) for \f$ \int_a^{a+u(b-a)} fp\,d(1-\hat{\psi}) \f$
!!
!! Exact integral of the cubic interpolating `fp_native` at the 4 Gauss points of [a,b],
!! so the full-interval value equals @ref flux_func_relabel_int.
!------------------------------------------------------------------------------
subroutine flux_func_relabel_cof(self,a,b,cof)
class(flux_func), intent(inout) :: self
real(r8), intent(in) :: a,b
real(r8), intent(out) :: cof(4)
real(r8), parameter :: xg(4)=[-0.861136311594053d0,-0.339981043584856d0, &
  0.339981043584856d0,0.861136311594053d0]
real(r8) :: u(4),fv(4),m(4),phihat,jac,x1,dpsi,s
integer(i4) :: i,j,k
CALL flux_func_norm(self,x1,dpsi)
u=(1.d0+xg)/2.d0
DO i=1,4
  s=a+u(i)*(b-a)
  CALL self%tmap%eval(s,phihat,jac)
  fv(i)=self%fp_native(x1+phihat*dpsi)
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
end subroutine flux_func_relabel_cof
!------------------------------------------------------------------------------
!> Evaluate \f$ (\psi_0-\psi_a) \int fp\,d(1-\hat{\psi}) \f$ from the LCFS to \f$ 1-\hat{\psi} \f$ (`coord=2`)
!------------------------------------------------------------------------------
function flux_func_relabel_f(self,psi) result(b)
class(flux_func), intent(inout) :: self
real(r8), intent(in) :: psi
real(r8) :: b,x1,dpsi,psihat,u
integer(i4) :: lo
CALL flux_func_norm(self,x1,dpsi)
psihat=(psi-x1)/dpsi
IF(self%tmap%ns<2)THEN
  b=flux_func_relabel_int(self,0.d0,psihat)
ELSE IF(self%fint_stamp/=self%tmap%stamp.OR.(.NOT.ASSOCIATED(self%fint)))THEN
  !---Table unavailable, integrate over all map intervals (slow)
  b=flux_func_relabel_int(self,0.d0,MIN(psihat,self%tmap%psis(1)))
  DO lo=2,self%tmap%ns
    IF(psihat<=self%tmap%psis(lo-1))EXIT
    b=b+flux_func_relabel_int(self,self%tmap%psis(lo-1),MIN(psihat,self%tmap%psis(lo)))
  END DO
  IF(psihat>self%tmap%psis(self%tmap%ns))b=b+flux_func_relabel_int(self,self%tmap%psis(self%tmap%ns),psihat)
ELSE IF(psihat<=self%tmap%psis(1))THEN
  b=flux_func_relabel_int(self,0.d0,psihat)
ELSE
  lo=tmap_lookup(self%tmap,psihat)
  u=(psihat-self%tmap%psis(lo))/(self%tmap%psis(lo+1)-self%tmap%psis(lo))
  b=self%fint(lo)+u*(self%fcof(1,lo)+u*(self%fcof(2,lo)+u*(self%fcof(3,lo)+u*self%fcof(4,lo))))
END IF
b=b*dpsi
end function flux_func_relabel_f
!------------------------------------------------------------------------------
!> Update plasma bounds, toroidal flux map, and F*F'/P flux functions
!------------------------------------------------------------------------------
MODULE PROCEDURE gs_update_flux_funcs
integer(i4) :: map_err
CALL gs_update_bounds(equil)
CALL gs_update_torflux_map(equil,settle,map_err)
IF(PRESENT(ierr))ierr=map_err
IF(map_err==2)RETURN ! Nothing sensible to evaluate profiles on
CALL equil%I%update(equil)
CALL equil%p%update(equil)
IF(ASSOCIATED(equil%P_ani))CALL equil%P_ani%update(equil)
!---Remap all toroidal flux functions to current bounds and map
CALL remap(equil%I)
CALL remap(equil%P)
CALL remap(equil%eta)
CALL remap(equil%I_NI)
CALL remap(equil%Te)
CALL remap(equil%Ti)
CALL remap(equil%ne)
CALL remap(equil%ni)
CALL remap(equil%Zeff)
CALL remap(equil%jphi_fixed)
!---Value-type profiles cannot use derivative (phi_n) semantics
CALL check_value(equil%eta,'eta')
CALL check_value(equil%Te,'Te')
CALL check_value(equil%Ti,'Ti')
CALL check_value(equil%ne,'ne')
CALL check_value(equil%ni,'ni')
CALL check_value(equil%Zeff,'Zeff')
CALL check_value(equil%jphi_fixed,'jphi_fixed')
CONTAINS
!---Abort if a value-type profile uses coord=1
subroutine check_value(F,name)
CLASS(flux_func), POINTER, INTENT(in) :: F
CHARACTER(LEN=*), INTENT(in) :: name
IF(.NOT.ASSOCIATED(F))RETURN
IF(F%coord==1)CALL oft_abort(name//' profile on toroidal flux must use "phi_n_relabel"', &
  'gs_update_flux_funcs',__FILE__)
end subroutine check_value
!---Rebuild relabel table against current map (bounds are owned by each profile's update)
subroutine remap(F)
CLASS(flux_func), POINTER, INTENT(in) :: F
IF(.NOT.ASSOCIATED(F))RETURN
IF(F%coord==2)CALL F%build_fint()
end subroutine remap
END PROCEDURE gs_update_flux_funcs
!------------------------------------------------------------------------------
!> Copy toroidal flux map and attach it to the flux functions of `self`
!------------------------------------------------------------------------------
MODULE PROCEDURE gs_copy_torflux
CALL self%tmap%delete()
self%tmap%active=source%tmap%active
self%tmap%Q=source%tmap%Q
self%tmap%ns=source%tmap%ns
IF(source%tmap%ns>0)THEN
  ALLOCATE(self%tmap%psis,SOURCE=source%tmap%psis)
  ALLOCATE(self%tmap%phihat,SOURCE=source%tmap%phihat)
  ALLOCATE(self%tmap%jac,SOURCE=source%tmap%jac)
  CALL tmap_index(self%tmap)
END IF
CALL attach(self%I)
CALL attach(self%P)
CALL attach(self%eta)
CALL attach(self%I_NI)
CALL attach(self%Te)
CALL attach(self%Ti)
CALL attach(self%ne)
CALL attach(self%ni)
CALL attach(self%Zeff)
CALL attach(self%jphi_fixed)
CONTAINS
!---Point at copied map and rebuild relabel table
subroutine attach(F)
CLASS(flux_func), POINTER, INTENT(in) :: F
IF(.NOT.ASSOCIATED(F))RETURN
F%tmap=>self%tmap
IF(F%coord==2)CALL F%build_fint()
end subroutine attach
END PROCEDURE gs_copy_torflux
!------------------------------------------------------------------------------
!> Storage coordinate name for profile files (empty for psi_n)
!------------------------------------------------------------------------------
MODULE PROCEDURE flux_coord_name
SELECT CASE(coord)
  CASE(1)
    name='phi_n'
  CASE(2)
    name='phi_n_relabel'
  CASE DEFAULT
    name=''
END SELECT
END PROCEDURE flux_coord_name
!------------------------------------------------------------------------------
!> Update \f$ \hat{\psi} \rightarrow \hat{\Phi} \f$ map from equilibrium q profile
!!
!! Surfaces are placed at the \f$ \hat{\psi} \f$ of each profile node, with gaps filled to
!! `dphi_max`, plus fixed surfaces near the LCFS and axis. They are re-placed on the updated map
!! until they move less than `place_tol` (single pass if `settle=.FALSE.`, used inside the nonlinear loop). The geometric factor \f$ q/F \f$ is traced once;
!! consistency between F and the map (when F*F' is itself on \f$ \hat{\Phi} \f$) is then
!! obtained by fixed-point iteration without further tracing.
!------------------------------------------------------------------------------
subroutine gs_update_torflux_map(gseq,settle,ierr)
class(gs_equil), target, intent(inout) :: gseq !< G-S object
logical, optional, intent(in) :: settle !< Re-place surfaces until converged (default: `.TRUE.`)
integer(i4), intent(out) :: ierr !< 0 ok, 1 update failed (previous map kept), 2 no closed surfaces (degenerate or non-finite bounds)
real(r8), parameter :: pad = 1.d-3 !< Minimum offset of node surfaces from LCFS and axis
real(r8), parameter :: edge_s(2) = [1.d-4,1.d-3] !< Fixed surfaces near the LCFS (edge singularity fit)
real(r8), parameter :: axis_s(2) = [1.d-3,3.d-3] !< Fixed surfaces near the axis (axis extrapolation)
real(r8), parameter :: dphi_max = 1.d-2 !< Maximum spacing between surfaces in \f$ \hat{\Phi} \f$
real(r8), parameter :: fp_tol = 1.d-10 !< Fixed-point tolerance on \f$ \hat{\Phi} \f$
integer(i4), parameter :: fp_maxits = 20
real(r8), parameter :: place_tol = 1.d-8 !< Surface placement tolerance in \f$ \hat{\psi} \f$
integer(i4), parameter :: place_maxits = 5
logical :: active
integer(i4) :: i,k,n,nfill,it,ip,npass
real(r8) :: dpsi,b,h,err,shift
real(r8), allocatable :: phinodes(:),phitarg(:),psis(:),pnew(:),g(:),q(:),r(:),dq(:),cum(:)
ierr=0
!---Collect toroidal flux functions and their nodes
ALLOCATE(phinodes(0))
active=.FALSE.
CALL add_func(gseq%I)
CALL add_func(gseq%P)
CALL add_func(gseq%eta)
CALL add_func(gseq%I_NI)
CALL add_func(gseq%Te)
CALL add_func(gseq%Ti)
CALL add_func(gseq%ne)
CALL add_func(gseq%ni)
CALL add_func(gseq%Zeff)
CALL add_func(gseq%jphi_fixed)
gseq%tmap%active=active
IF((.NOT.active).OR.(gseq%plasma_bounds(1)<-1.d98))THEN
  IF(gseq%tmap%ns>0)CALL gseq%tmap%delete()
  RETURN
END IF
dpsi=gseq%plasma_bounds(2)-gseq%plasma_bounds(1)
IF(.NOT.((ABS(dpsi)>0.d0).AND.(ABS(dpsi)<HUGE(dpsi))))THEN
  CALL oft_warn('Toroidal flux map update skipped, no closed flux surfaces')
  ierr=2
  RETURN
END IF
!---Target 1-Phi-hat values: [0,1] + nodes, gaps filled to dphi_max
phinodes=[0.d0,1.d0,MIN(MAX(phinodes,0.d0),1.d0)]
CALL sort_unique(phinodes)
ALLOCATE(phitarg(0))
DO k=1,SIZE(phinodes)-1
  nfill=MAX(CEILING((phinodes(k+1)-phinodes(k))/dphi_max),1)
  phitarg=[phitarg,(phinodes(k)+(phinodes(k+1)-phinodes(k))*REAL(i,r8)/REAL(nfill,r8), i=0,nfill-1)]
END DO
phitarg=[phitarg,1.d0]
!---Place surfaces at node locations of the current map and re-trace until they settle
shift=-1.d0
npass=place_maxits
IF(PRESENT(settle))THEN
  IF(.NOT.settle)npass=1
END IF
DO ip=1,npass
  ALLOCATE(pnew(SIZE(phitarg)))
  DO i=1,SIZE(phitarg)
    pnew(i)=MIN(MAX(gseq%tmap%inv(phitarg(i)),pad),1.d0-pad)
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
  IF(ASSOCIATED(gseq%I))gseq%I%plasma_bounds=gseq%plasma_bounds
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
    IF(gseq%tmap%ns==n)err=MAXVAL(ABS(cum/cum(n)-gseq%tmap%phihat))
    IF(gseq%tmap%ns/=n)THEN
      CALL gseq%tmap%delete()
      gseq%tmap%ns=n
      ALLOCATE(gseq%tmap%psis(n),gseq%tmap%phihat(n),gseq%tmap%jac(n))
    END IF
    gseq%tmap%Q=cum(n)
    gseq%tmap%psis=psis
    gseq%tmap%phihat=cum/cum(n)
    gseq%tmap%jac=q/cum(n)
    gseq%tmap%stamp=gseq%tmap%stamp+1
    CALL tmap_index(gseq%tmap)
    IF(gseq%I%coord==0)EXIT
    IF(gseq%I%coord==2)CALL gseq%I%build_fint()
    IF(err<fp_tol)EXIT
  END DO
END DO
IF(oft_debug_print(1))THEN
  WRITE(*,'(2A,I4,A,I3,A,I3,A,ES11.3,A,ES11.3)')oft_indent,'Toroidal flux map: ',n,' surfaces, ', &
    MIN(ip,npass),' passes, ',it,' its, shift = ',shift,', Q = ',gseq%tmap%Q
END IF
IF(ALLOCATED(pnew))DEALLOCATE(pnew)
DEALLOCATE(phinodes,phitarg,psis,g,q,r,dq,cum)
CONTAINS
!---Warn and flag a failed update (previous map is kept)
subroutine map_failed(reason)
CHARACTER(LEN=*), INTENT(in) :: reason
CALL oft_warn('Toroidal flux map update failed ('//reason//'), keeping previous map')
ierr=1
end subroutine map_failed
!---Register flux function with map
subroutine add_func(F)
CLASS(flux_func), POINTER, INTENT(in) :: F
real(r8), allocatable :: nodes(:)
IF(.NOT.ASSOCIATED(F))RETURN
F%tmap=>gseq%tmap
IF(F%coord==0)RETURN
active=.TRUE.
CALL F%get_nodes(nodes)
phinodes=[phinodes,nodes]
end subroutine add_func
!---Poloidal current function F(psi)
function fpol(psi) result(F)
real(r8), intent(in) :: psi
real(r8) :: F
IF(gseq%mode==0)THEN
  F=gseq%ffp_scale*gseq%I%f(psi)+gseq%I%f_offset
ELSE
  F=SIGN(1.d0,gseq%I%f_offset)*SQRT(MAX(gseq%ffp_scale*gseq%I%f(psi)+gseq%I%f_offset**2,0.d0))
END IF
end function fpol
!---Antiderivative of ln(s)
pure function xlogx(s) result(v)
real(r8), intent(in) :: s
real(r8) :: v
v=s*LOG(s)-s
end function xlogx
!---Replace failed traces (g=0) by linear interpolation, .FALSE. if too few remain
function fill_failed() result(ok)
logical :: ok
integer(i4) :: j,jl,jr
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
end subroutine gs_update_torflux_map
!------------------------------------------------------------------------------
!> Geometric factor \f$ g = q/F = \frac{1}{2\pi} \oint \frac{dl}{R |\nabla \psi|} \f$ on flux surfaces
!!
!! Single backend entry point for the toroidal flux map. Contract:
!! - `psi_q` in internal normalized flux (0 at LCFS, 1 at axis), strictly inside (0,1)
!! - `g` in the units of @ref gs_get_qprof with F=1 (sign ignored by caller)
!! - `g(i)=0` flags a failed surface (caller interpolates)
!! - Result must vary smoothly with \f$ \psi \f$ between nonlinear iterations
!------------------------------------------------------------------------------
subroutine torflux_qgeom(gseq,nr,psi_q,g)
class(gs_equil), intent(inout) :: gseq !< G-S object
integer(i4), intent(in) :: nr !< Number of surfaces
real(r8), intent(in) :: psi_q(nr) !< Surface locations
real(r8), intent(out) :: g(nr) !< Geometric factor at each surface
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
real(r8), allocatable, intent(inout) :: x(:)
real(r8), parameter :: tol=1.d-8
integer(i4) :: i,j,m
real(r8) :: tmp
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
real(r8), intent(in) :: x(3),y(3),xi
real(r8) :: yi
yi=y(1)*(xi-x(2))*(xi-x(3))/((x(1)-x(2))*(x(1)-x(3))) &
  + y(2)*(xi-x(1))*(xi-x(3))/((x(2)-x(1))*(x(2)-x(3))) &
  + y(3)*(xi-x(1))*(xi-x(2))/((x(3)-x(1))*(x(3)-x(2)))
end function quad_interp
!------------------------------------------------------------------------------
!> Derivative of 3-point Lagrange interpolant at point `k` of `x`
!------------------------------------------------------------------------------
pure function fd_deriv(x,y,k) result(dy)
real(r8), intent(in) :: x(3),y(3)
integer(i4), intent(in) :: k
real(r8) :: dy,xi
xi=x(k)
dy=y(1)*(2.d0*xi-x(2)-x(3))/((x(1)-x(2))*(x(1)-x(3))) &
  + y(2)*(2.d0*xi-x(1)-x(3))/((x(2)-x(1))*(x(2)-x(3))) &
  + y(3)*(2.d0*xi-x(1)-x(2))/((x(3)-x(1))*(x(3)-x(2)))
end function fd_deriv
END SUBMODULE oft_gs_torflux
