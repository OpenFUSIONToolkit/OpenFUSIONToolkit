!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file grad_shaf_cutcell.F90
!
!> Mesh-based (cut-cell) evaluation of the flux-surface geometric factor \f$ q/F \f$
!!
!! Cells are split (Bernstein form) until \f$ \hat{\psi} \f$ is monotone along a direction \f$ d \f$;
!! each level-set arc is then a graph between exact edge crossings, integrated by Gauss quadrature
!! with \f$ dl/|\nabla \psi| = |J| d\xi / |\partial_d \psi| \f$.
!!
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
MODULE oft_gs_cutcell
USE oft_base
USE oft_mesh_type, ONLY: cell_is_curved
USE oft_lag_basis, ONLY: oft_blag_npos, oft_scalar_bfem
USE oft_gs, ONLY: gs_equil, gs_test_bounds
USE oft_lu, ONLY: lapack_matinv
IMPLICIT NONE
#include "local.h"
PRIVATE
!> Bernstein tables for degree p on the triangle
TYPE :: bern_tab
  INTEGER(i4) :: p = 0 !< Degree
  INTEGER(i4) :: nb = 0 !< Number of degree-p coefficients
  INTEGER(i4) :: nb1 = 0 !< Number of degree-(p-1) coefficients
  INTEGER(i4), ALLOCATABLE :: mi(:,:) !< Degree-p multi-indices [3,nb]
  INTEGER(i4), ALLOCATABLE :: mi1(:,:) !< Degree-(p-1) multi-indices [3,nb1]
  INTEGER(i4), ALLOCATABLE :: lk(:,:) !< Index of (p-a2-a3,a2,a3) [0:p,0:p]
  INTEGER(i4) :: corner(3) = 0 !< Index of vertex coefficients
  REAL(r8), ALLOCATABLE :: vinv(:,:) !< Uniform-lattice values to coefficients [nb,nb]
  REAL(r8), ALLOCATABLE :: mvinv(:,:) !< Uniform-lattice values to monomial coefficients [nb,nb]
  REAL(r8) :: cbound = 0.d0 !< Coefficients lie within nodal range extended by `cbound` times its width
END TYPE bern_tab
TYPE(bern_tab), SAVE :: bt_cache !< Tables of the last degree used (built once per degree)
INTEGER(i4), PARAMETER :: nblk_cells = 64 !< Cells per deterministic accumulation block
INTEGER(i4), PARAMETER :: max_depth = 12 !< Maximum subdivision depth
REAL(r8), PARAMETER :: kappa = 0.25d0 !< Required min/max ratio of \f$ \partial_d \hat{\psi} \f$ per sub-cell
REAL(r8), PARAMETER :: wind_tol = 1.d-2 !< Tolerance on winding number about O-point
REAL(r8), PARAMETER :: smooth_ratio = 0.9d0 !< Min/max d psi/d d ratio above which 3 Gauss points are used (else 5)
REAL(r8), PARAMETER :: xg3(3) = [-0.774596669241483377035853d0,0.d0,0.774596669241483377035853d0]
REAL(r8), PARAMETER :: wg3(3) = [0.555555555555555555555556d0,0.888888888888888888888889d0, &
  0.555555555555555555555556d0]
REAL(r8), PARAMETER :: xg5(5) = [-0.906179845938663992797627d0,-0.538469310105683091036314d0, &
  0.d0,0.538469310105683091036314d0,0.906179845938663992797627d0]
REAL(r8), PARAMETER :: wg5(5) = [0.236926885056189087514264d0,0.478628670499366468041292d0, &
  0.568888888888888888888889d0,0.478628670499366468041292d0,0.236926885056189087514264d0]
!---Cached lattice index of each cell DOF (depends only on mesh and FE order)
TYPE(oft_scalar_bfem), POINTER :: perm_rep => NULL()
INTEGER(i8) :: perm_sum = -1
INTEGER(i4), ALLOCATABLE :: perm(:,:)
PUBLIC gs_qgeom_cutcell
CONTAINS
!------------------------------------------------------------------------------
!> Geometric factor \f$ g = q/F = \frac{1}{2\pi} \oint \frac{dl}{R |\nabla \psi|} \f$ on flux surfaces
!!
!! Same contract and units as @ref oft_gs::gs_get_qprof with `geom_only`. Only arcs inside
!! the plasma mask (@ref oft_gs::gs_test_bounds) contribute; surfaces whose winding about the
!! O-point is not one return `g=0`. Summation order is fixed, so results do not depend on
!! the thread count.
!------------------------------------------------------------------------------
SUBROUTINE gs_qgeom_cutcell(gseq,nr,psi_q,g)
CLASS(gs_equil), INTENT(inout) :: gseq !< G-S object
INTEGER(i4), INTENT(in) :: nr !< Number of surfaces
REAL(r8), INTENT(in) :: psi_q(nr) !< Surface locations (0 at LCFS, 1 at axis)
REAL(r8), INTENT(out) :: g(nr) !< Geometric factor at each surface
TYPE(bern_tab) :: bt
REAL(r8), POINTER :: vals(:)
REAL(r8), ALLOCATABLE :: gblk(:,:),wblk(:,:),gloc(:),wloc(:),wind(:)
REAL(r8) :: x1,dpsi
INTEGER(i4) :: nblk,ib,cell,i,nc
CHARACTER(LEN=OFT_ERROR_SLEN) :: error_str
g=0.d0
IF(nr<1)RETURN
IF(gseq%plasma_bounds(1)<-1.d98)RETURN
x1=gseq%plasma_bounds(1)
dpsi=gseq%plasma_bounds(2)-gseq%plasma_bounds(1)
IF(dpsi==0.d0)RETURN
CALL bern_setup(bt,gseq%device%fe_rep%order)
NULLIFY(vals)
CALL gseq%psi%get_local(vals)
nc=gseq%device%fe_rep%mesh%nc
CALL perm_update(gseq%device%fe_rep,bt)
nblk=(nc+nblk_cells-1)/nblk_cells
ALLOCATE(gblk(nr,nblk),wblk(nr,nblk))
!$omp parallel private(ib,cell,gloc,wloc)
ALLOCATE(gloc(nr),wloc(nr))
!$omp do schedule(dynamic,1)
DO ib=1,nblk
  gloc=0.d0; wloc=0.d0
  DO cell=(ib-1)*nblk_cells+1,MIN(ib*nblk_cells,nc)
    IF(gseq%device%fe_rep%mesh%reg(cell)/=1)CYCLE
    CALL cell_contrib(gseq,bt,vals,x1,dpsi,nr,psi_q,cell,gloc,wloc)
  END DO
  gblk(:,ib)=gloc; wblk(:,ib)=wloc
END DO
!$omp end do
DEALLOCATE(gloc,wloc)
!$omp end parallel
!---Fixed-order reduction
ALLOCATE(wind(nr))
wind=0.d0
DO ib=1,nblk
  g=g+gblk(:,ib)
  wind=wind+wblk(:,ib)
END DO
g=g/(2.d0*pi*ABS(dpsi))
DO i=1,nr
  IF((psi_q(i)<=0.d0).OR.(psi_q(i)>=1.d0))THEN
    g(i)=0.d0
  ELSE IF(ABS(wind(i)-1.d0)>wind_tol)THEN
    WRITE(error_str,'(A,F10.6,A,ES11.3)')'gs_qgeom_cutcell: open or missing surface at psi = ', &
      1.d0-psi_q(i),', winding = ',wind(i)
    CALL oft_warn(error_str)
    g(i)=0.d0
  END IF
END DO
DEALLOCATE(vals,gblk,wblk,wind)
END SUBROUTINE gs_qgeom_cutcell
!------------------------------------------------------------------------------
!> Add arc integrals and winding of all surfaces crossing `cell`
!------------------------------------------------------------------------------
SUBROUTINE cell_contrib(gseq,bt,vals,x1,dpsi,nr,psi_q,cell,gloc,wloc)
CLASS(gs_equil), INTENT(inout) :: gseq !< G-S object
TYPE(bern_tab), INTENT(in) :: bt !< Bernstein tables
REAL(r8), INTENT(in) :: vals(:) !< Local \f$ \psi \f$ values
REAL(r8), INTENT(in) :: x1 !< LCFS flux
REAL(r8), INTENT(in) :: dpsi !< Axis minus LCFS flux
INTEGER(i4), INTENT(in) :: nr !< Number of surfaces
REAL(r8), INTENT(in) :: psi_q(nr) !< Surface locations
INTEGER(i4), INTENT(in) :: cell !< Cell index
REAL(r8), INTENT(inout) :: gloc(nr) !< Accumulated \f$ \oint dl/(R|\nabla\hat{\psi}|) \f$
REAL(r8), INTENT(inout) :: wloc(nr) !< Accumulated winding number
TYPE(oft_scalar_bfem), POINTER :: lag_rep
INTEGER(i4) :: jdofs(32),jc,k,ns,d0,i
INTEGER(i4) :: dep(3*max_depth+4)
REAL(r8) :: lat(32),cc(32),cm(32),cs(32),tri(2,3),stk(2,3,3*max_depth+4),xv(2,3),detm0
REAL(r8) :: lo,hi,pad,dvec(2),dmu(3),ainv(2,2),m(2,3),dratio
LOGICAL :: curved,act(nr),sact(nr,3*max_depth+4)
lag_rep=>gseq%device%fe_rep
!---Bernstein coefficients of 1-psi-hat in this cell
CALL lag_rep%ncdofs(cell,jdofs)
lo=(MINVAL(vals(jdofs(1:lag_rep%nce)))-x1)/dpsi
hi=(MAXVAL(vals(jdofs(1:lag_rep%nce)))-x1)/dpsi
IF(lo>hi)CALL swap(lo,hi)
pad=bt%cbound*(hi-lo)
IF(.NOT.any_hit(lo-pad,hi+pad))RETURN
DO jc=1,lag_rep%nce
  lat(perm(jc,cell))=(vals(jdofs(jc))-x1)/dpsi
END DO
cc(1:bt%nb)=MATMUL(bt%vinv,lat(1:bt%nb))
cc(bt%corner)=lat(bt%corner)
cm(1:bt%nb)=MATMUL(bt%mvinv,lat(1:bt%nb))
IF(.NOT.any_hit(MINVAL(cc(1:bt%nb)),MAXVAL(cc(1:bt%nb))))RETURN
!---Straight-cell geometry
curved=cell_is_curved(gseq%device%fe_rep%mesh,cell)
DO k=1,3
  xv(:,k)=gseq%device%fe_rep%mesh%r(1:2,gseq%device%fe_rep%mesh%lc(k,cell))
END DO
detm0=(xv(1,2)-xv(1,1))*(xv(2,3)-xv(2,1))-(xv(2,2)-xv(2,1))*(xv(1,3)-xv(1,1))
!---Recursive subdivision (depth first); each entry carries its unresolved surfaces
ns=1
stk(:,:,1)=RESHAPE([0.d0,0.d0,1.d0,0.d0,0.d0,1.d0],[2,3])
dep(1)=0
sact(:,1)=.TRUE.
DO WHILE(ns>0)
  tri=stk(:,:,ns); d0=dep(ns); act=sact(:,ns); ns=ns-1
  IF(d0==0)THEN
    cs(1:bt%nb)=cc(1:bt%nb)
  ELSE
    CALL sub_coefs(bt,cm,tri,cs)
  END IF
  lo=MINVAL(cs(1:bt%nb)); hi=MAXVAL(cs(1:bt%nb))
  act=act.AND.(psi_q>=lo).AND.(psi_q<=hi)
  IF(.NOT.ANY(act))CYCLE
  IF(tri_simple(bt,cm,cs,tri,dvec,dmu,ainv,dratio))THEN
    DO i=1,nr
      IF(.NOT.act(i))CYCLE
      IF((d0<max_depth).AND.(.NOT.edges_ok(psi_q(i))))CYCLE
      CALL tri_arc(psi_q(i),gloc(i),wloc(i))
      act(i)=.FALSE.
    END DO
  END IF
  IF((d0<max_depth).AND.ANY(act))THEN
    m(:,1)=(tri(:,1)+tri(:,2))/2.d0
    m(:,2)=(tri(:,2)+tri(:,3))/2.d0
    m(:,3)=(tri(:,3)+tri(:,1))/2.d0
    stk(:,:,ns+1)=RESHAPE([tri(:,1),m(:,1),m(:,3)],[2,3])
    stk(:,:,ns+2)=RESHAPE([m(:,1),tri(:,2),m(:,2)],[2,3])
    stk(:,:,ns+3)=RESHAPE([m(:,3),m(:,2),tri(:,3)],[2,3])
    stk(:,:,ns+4)=RESHAPE([m(:,2),m(:,3),m(:,1)],[2,3])
    dep(ns+1:ns+4)=d0+1
    DO k=1,4
      sact(:,ns+k)=act
    END DO
    ns=ns+4
  END IF
END DO
CONTAINS
!---Exchange two values
SUBROUTINE swap(a,b)
REAL(r8), INTENT(inout) :: a,b
REAL(r8) :: t
t=a; a=b; b=t
END SUBROUTINE swap
!---Does any surface lie in [a,b]?
LOGICAL FUNCTION any_hit(a,b)
REAL(r8), INTENT(in) :: a,b
any_hit=ANY((psi_q>=a).AND.(psi_q<=b))
END FUNCTION any_hit
!---At most one crossing of s per sub-cell edge (Bernstein sign changes)
LOGICAL FUNCTION edges_ok(s)
REAL(r8), INTENT(in) :: s
INTEGER(i4) :: kk,ka,kb,q,nch,al(3)
LOGICAL :: up,upl
edges_ok=.FALSE.
DO kk=1,3
  ka=MOD(kk,3)+1; kb=MOD(kk+1,3)+1
  nch=0
  DO q=0,bt%p
    al=0; al(ka)=bt%p-q; al(kb)=q
    up=(cs(bt%lk(al(2),al(3)))>=s)
    IF(q>0.AND.(up.NEQV.upl))nch=nch+1
    upl=up
  END DO
  IF(nch>1)RETURN
END DO
edges_ok=.TRUE.
END FUNCTION edges_ok
!---Integrate the arc of surface s in the current (simple) sub-cell
SUBROUTINE tri_arc(s,gacc,wacc)
REAL(r8), INTENT(in) :: s
REAL(r8), INTENT(inout) :: gacc,wacc
INTEGER(i4) :: k0,k1,k2,ig,kk,nab,ng
LOGICAL :: ab(3)
REAL(r8) :: p1(2),p2(2),e(2),xi1,xi2,xi,c(2),muc(3),elo,ehi,eta,x(2),val,gr(2)
REAL(r8) :: pt(2),detm,sum,q1(2),q2(2),dth,sig,gmid(2),xg(5),wg(5)
DO kk=1,3
  ab(kk)=(cs(bt%corner(kk))>=s)
END DO
nab=COUNT(ab)
IF((nab==0).OR.(nab==3))RETURN
DO k0=1,3
  IF(ab(k0).EQV.(nab==1))EXIT
END DO
k1=MOD(k0,3)+1; k2=MOD(k0+1,3)+1
!---Exact edge crossings
p1=edge_cross(s,k0,k1)
p2=edge_cross(s,k0,k2)
!---Plasma mask at chord midpoint
CALL cell_geom((p1+p2)/2.d0,pt,detm)
IF(.NOT.gs_test_bounds(gseq,pt))RETURN
e=[-dvec(2),dvec(1)]
xi1=DOT_PRODUCT(e,p1); xi2=DOT_PRODUCT(e,p2)
IF(xi1==xi2)RETURN
!---Gauss quadrature over xi, projecting onto the level set along d
sum=0.d0
IF(dratio>=smooth_ratio)THEN
  ng=3; xg(1:3)=xg3; wg(1:3)=wg3
ELSE
  ng=5; xg=xg5; wg=wg5
END IF
DO ig=1,ng
  xi=(xi1+xi2)/2.d0+xg(ig)*(xi2-xi1)/2.d0
  c=p1+((xi-xi1)/(xi2-xi1))*(p2-p1)
  muc(2:3)=MATMUL(ainv,c-tri(:,1))
  muc(1)=1.d0-muc(2)-muc(3)
  elo=-HUGE(1.d0); ehi=HUGE(1.d0)
  DO kk=1,3
    IF(dmu(kk)>0.d0)THEN
      elo=MAX(elo,-muc(kk)/dmu(kk))
    ELSE IF(dmu(kk)<0.d0)THEN
      ehi=MIN(ehi,-muc(kk)/dmu(kk))
    END IF
  END DO
  elo=MIN(elo,0.d0); ehi=MAX(ehi,0.d0)
  eta=line_root(bt,cm,s,c,dvec,elo,ehi,0.d0)
  x=c+eta*dvec
  CALL poly_eval(bt,cm,x,val,gr)
  CALL cell_geom(x,pt,detm)
  sum=sum+wg(ig)*ABS(detm)/(pt(1)*ABS(DOT_PRODUCT(gr,dvec)))
  IF(ig==(ng+1)/2)gmid=gr
END DO
gacc=gacc+sum*ABS(xi2-xi1)/2.d0
!---Winding about O-point, oriented with 1-psi-hat increasing to the left
CALL cell_geom(p1,q1,detm)
CALL cell_geom(p2,q2,detm)
q1=q1-gseq%o_point; q2=q2-gseq%o_point
dth=ATAN2(q1(1)*q2(2)-q1(2)*q2(1),DOT_PRODUCT(q1,q2))
sig=SIGN(1.d0,(p2(1)-p1(1))*gmid(2)-(p2(2)-p1(2))*gmid(1))*SIGN(1.d0,detm)
wacc=wacc+sig*dth/(2.d0*pi)
END SUBROUTINE tri_arc
!---Crossing of surface s on sub-cell edge ka->kb (regula falsi start)
FUNCTION edge_cross(s,ka,kb) RESULT(pc)
REAL(r8), INTENT(in) :: s
INTEGER(i4), INTENT(in) :: ka,kb
REAL(r8) :: pc(2),fa,fb,t
fa=cs(bt%corner(ka))-s; fb=cs(bt%corner(kb))-s
t=line_root(bt,cm,s,tri(:,ka),tri(:,kb)-tri(:,ka),0.d0,1.d0,fa/(fa-fb))
pc=tri(:,ka)+t*(tri(:,kb)-tri(:,ka))
END FUNCTION edge_cross
!---Physical position and logical Jacobian determinant at uv
SUBROUTINE cell_geom(uv,pt,detm)
REAL(r8), INTENT(in) :: uv(2)
REAL(r8), INTENT(out) :: pt(2),detm
REAL(r8) :: lam(3),gop(3,3),p3(3)
IF(curved)THEN
  lam=[1.d0-uv(1)-uv(2),uv(1),uv(2)]
  CALL gseq%device%fe_rep%mesh%jacobian(cell,lam,gop,detm)
  p3=gseq%device%fe_rep%mesh%log2phys(cell,lam)
  pt=p3(1:2)
ELSE
  pt=xv(:,1)+uv(1)*(xv(:,2)-xv(:,1))+uv(2)*(xv(:,3)-xv(:,1))
  detm=detm0
END IF
END SUBROUTINE cell_geom
END SUBROUTINE cell_contrib
!------------------------------------------------------------------------------
!> Root of \f$ (1-\hat{\psi})(x_0+t\,d)=s \f$ for \f$ t \in [a,b] \f$
!!
!! Newton from `t0`, falling back to safeguarded Newton-bisection on [a,b]. Returns the
!! nearer end if the bracket does not straddle the root. `cc` holds monomial coefficients.
!------------------------------------------------------------------------------
FUNCTION line_root(bt,cc,s,x0,dir,a,b,t0) RESULT(t)
TYPE(bern_tab), INTENT(in) :: bt
REAL(r8), INTENT(in) :: cc(:),s,x0(2),dir(2),a,b,t0
REAL(r8) :: t
REAL(r8) :: ta,tb,fa,fb,f,fp,gr(2),tn
INTEGER(i4) :: it
!---Newton
t=t0
DO it=1,8
  CALL poly_eval(bt,cc,x0+t*dir,f,gr)
  fp=DOT_PRODUCT(gr,dir)
  IF(fp==0.d0)EXIT
  tn=t-(f-s)/fp
  IF((tn<a).OR.(tn>b))EXIT
  IF(ABS(tn-t)<=1.d-8*(b-a))THEN
    t=tn; RETURN ! Remaining error ~ step**2
  END IF
  t=tn
END DO
!---Safeguarded fallback
ta=a; tb=b
CALL poly_eval(bt,cc,x0+ta*dir,fa,gr); fa=fa-s
CALL poly_eval(bt,cc,x0+tb*dir,fb,gr); fb=fb-s
IF(fa==0.d0)THEN
  t=ta; RETURN
ELSE IF(fb==0.d0)THEN
  t=tb; RETURN
ELSE IF(SIGN(1.d0,fa)==SIGN(1.d0,fb))THEN
  t=MERGE(ta,tb,ABS(fa)<ABS(fb)); RETURN
END IF
t=ta-fa*(tb-ta)/(fb-fa)
DO it=1,100
  CALL poly_eval(bt,cc,x0+t*dir,f,gr)
  f=f-s
  IF(f==0.d0)EXIT
  IF(SIGN(1.d0,f)==SIGN(1.d0,fa))THEN
    ta=t; fa=f
  ELSE
    tb=t; fb=f
  END IF
  fp=DOT_PRODUCT(gr,dir)
  tn=t-f/fp
  IF((fp==0.d0).OR.(tn<=ta).OR.(tn>=tb))tn=(ta+tb)/2.d0
  IF(ABS(tn-t)<=1.d-15*(b-a))THEN
    t=tn; EXIT
  END IF
  t=tn
END DO
END FUNCTION line_root
!------------------------------------------------------------------------------
!> Test if \f$ \hat{\psi} \f$ is strictly monotone (with bounded variation) along the
!! centroid gradient direction over sub-cell `tri`
!------------------------------------------------------------------------------
FUNCTION tri_simple(bt,cc,cs,tri,dvec,dmu,ainv,dratio) RESULT(simple)
TYPE(bern_tab), INTENT(in) :: bt
REAL(r8), INTENT(in) :: cc(:) !< Cell monomial coefficients
REAL(r8), INTENT(in) :: cs(:) !< Sub-cell Bernstein coefficients
REAL(r8), INTENT(in) :: tri(2,3) !< Sub-cell vertices (uv)
REAL(r8), INTENT(out) :: dvec(2) !< Monotone direction (uv, unit)
REAL(r8), INTENT(out) :: dmu(3) !< `dvec` in sub-cell barycentric increments
REAL(r8), INTENT(out) :: ainv(2,2) !< Map from uv offset to sub-cell (mu2,mu3)
REAL(r8), INTENT(out) :: dratio !< Lower bound of min/max \f$ \partial_d \hat{\psi} \f$ over the sub-cell
LOGICAL :: simple
REAL(r8) :: val,gr(2),gn,a(2,2),det,dc,dcmin,dcmax
INTEGER(i4) :: n,k,al(3)
simple=.FALSE.
dratio=0.d0
CALL poly_eval(bt,cc,SUM(tri,DIM=2)/3.d0,val,gr)
gn=SQRT(SUM(gr**2))
IF(gn<=0.d0)RETURN
dvec=gr/gn
a(:,1)=tri(:,2)-tri(:,1); a(:,2)=tri(:,3)-tri(:,1)
det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
ainv(1,1)=a(2,2)/det; ainv(1,2)=-a(1,2)/det
ainv(2,1)=-a(2,1)/det; ainv(2,2)=a(1,1)/det
dmu(2:3)=MATMUL(ainv,dvec)
dmu(1)=-dmu(2)-dmu(3)
!---Directional derivative coefficients (degree p-1)
dcmin=HUGE(1.d0); dcmax=-HUGE(1.d0)
DO n=1,bt%nb1
  dc=0.d0
  DO k=1,3
    al=bt%mi1(:,n); al(k)=al(k)+1
    dc=dc+dmu(k)*cs(bt%lk(al(2),al(3)))
  END DO
  dcmin=MIN(dcmin,dc); dcmax=MAX(dcmax,dc)
END DO
IF((dcmax<=0.d0).OR.(dcmin<=kappa*dcmax))RETURN
dratio=dcmin/dcmax
simple=.TRUE.
END FUNCTION tri_simple
!------------------------------------------------------------------------------
!> Bernstein coefficients on sub-cell `tri` of the cell polynomial (monomial coefficients `cc`)
!------------------------------------------------------------------------------
SUBROUTINE sub_coefs(bt,cc,tri,cs)
TYPE(bern_tab), INTENT(in) :: bt
REAL(r8), INTENT(in) :: cc(:),tri(2,3)
REAL(r8), INTENT(out) :: cs(:)
REAL(r8) :: lat(32),gr(2),uv(2)
INTEGER(i4) :: n
DO n=1,bt%nb
  uv=MATMUL(tri,REAL(bt%mi(:,n),r8))/REAL(bt%p,r8)
  CALL poly_eval(bt,cc,uv,lat(n),gr)
END DO
cs(1:bt%nb)=MATMUL(bt%vinv,lat(1:bt%nb))
cs(bt%corner)=lat(bt%corner)
END SUBROUTINE sub_coefs
!------------------------------------------------------------------------------
!> Evaluate cell polynomial (monomial coefficients) and its (u,v) gradient
!------------------------------------------------------------------------------
SUBROUTINE poly_eval(bt,c,uv,val,gr)
TYPE(bern_tab), INTENT(in) :: bt
REAL(r8), INTENT(in) :: c(:),uv(2)
REAL(r8), INTENT(out) :: val,gr(2)
REAL(r8) :: up(-1:4),vp(-1:4)
INTEGER(i4) :: n,i,j
up(-1)=0.d0; vp(-1)=0.d0; up(0)=1.d0; vp(0)=1.d0
DO i=1,bt%p
  up(i)=up(i-1)*uv(1); vp(i)=vp(i-1)*uv(2)
END DO
val=0.d0; gr=0.d0
DO n=1,bt%nb
  i=bt%mi(2,n); j=bt%mi(3,n)
  val=val+c(n)*up(i)*vp(j)
  gr(1)=gr(1)+c(n)*i*up(i-1)*vp(j)
  gr(2)=gr(2)+c(n)*j*up(i)*vp(j-1)
END DO
END SUBROUTINE poly_eval
!------------------------------------------------------------------------------
!> Refresh cached lattice index of each cell DOF when mesh or FE space changes
!------------------------------------------------------------------------------
SUBROUTINE perm_update(lag_rep,bt)
TYPE(oft_scalar_bfem), POINTER, INTENT(in) :: lag_rep
TYPE(bern_tab), INTENT(in) :: bt
INTEGER(i8) :: chk
INTEGER(i4) :: cell,jc
REAL(r8) :: f(3)
chk=SUM(INT(lag_rep%mesh%lc,i8))+SUM(INT(ABS(lag_rep%mesh%lce),i8)*INT(SIGN(2,lag_rep%mesh%lce)+1,i8)) &
  +INT(lag_rep%order,i8)
IF(ASSOCIATED(perm_rep,lag_rep).AND.ALLOCATED(perm).AND.(perm_sum==chk))THEN
  IF(SIZE(perm,1)==lag_rep%nce.AND.SIZE(perm,2)==lag_rep%mesh%nc)RETURN
END IF
IF(ALLOCATED(perm))DEALLOCATE(perm)
ALLOCATE(perm(lag_rep%nce,lag_rep%mesh%nc))
DO cell=1,lag_rep%mesh%nc
  DO jc=1,lag_rep%nce
    CALL oft_blag_npos(lag_rep,cell,jc,f)
    perm(jc,cell)=bt%lk(NINT(f(2)*bt%p),NINT(f(3)*bt%p))
  END DO
END DO
perm_rep=>lag_rep
perm_sum=chk
END SUBROUTINE perm_update
!------------------------------------------------------------------------------
!> Bernstein tables for degree `p` (1-4), built on first use of each degree
!------------------------------------------------------------------------------
SUBROUTINE bern_setup(bt,p)
TYPE(bern_tab), INTENT(inout) :: bt
INTEGER(i4), INTENT(in) :: p
IF(bt_cache%p/=p)THEN
  IF(ALLOCATED(bt_cache%mi))DEALLOCATE(bt_cache%mi,bt_cache%mi1,bt_cache%lk,bt_cache%vinv,bt_cache%mvinv)
  CALL bern_build(bt_cache,p)
END IF
bt=bt_cache
END SUBROUTINE bern_setup
!------------------------------------------------------------------------------
!> Build Bernstein tables for degree `p` (1-4)
!------------------------------------------------------------------------------
SUBROUTINE bern_build(bt,p)
TYPE(bern_tab), INTENT(inout) :: bt
INTEGER(i4), INTENT(in) :: p
REAL(r8) :: v(15,15),fact(0:4),mult(15),x
INTEGER(i4) :: n,m,a2,a3,k,ierr
IF((p<1).OR.(p>4))CALL oft_abort('Unsupported FE order','bern_build',__FILE__)
bt%p=p
bt%nb=(p+1)*(p+2)/2
bt%nb1=p*(p+1)/2
ALLOCATE(bt%mi(3,bt%nb),bt%mi1(3,bt%nb1),bt%lk(0:p,0:p))
fact(0)=1.d0
DO k=1,4
  fact(k)=fact(k-1)*k
END DO
bt%lk=0
n=0
DO a3=0,p
  DO a2=0,p-a3
    n=n+1
    bt%mi(:,n)=[p-a2-a3,a2,a3]
    bt%lk(a2,a3)=n
    mult(n)=fact(p)/(fact(p-a2-a3)*fact(a2)*fact(a3))
  END DO
END DO
n=0
DO a3=0,p-1
  DO a2=0,p-1-a3
    n=n+1
    bt%mi1(:,n)=[p-1-a2-a3,a2,a3]
  END DO
END DO
bt%corner=[bt%lk(0,0),bt%lk(p,0),bt%lk(0,p)]
!---Values of each basis function on the uniform lattice
DO m=1,bt%nb
  DO n=1,bt%nb
    x=mult(n)
    DO k=1,3
      IF(bt%mi(k,n)>0)x=x*(REAL(bt%mi(k,m),r8)/REAL(p,r8))**bt%mi(k,n)
    END DO
    v(m,n)=x
  END DO
END DO
ALLOCATE(bt%vinv,SOURCE=v(1:bt%nb,1:bt%nb))
CALL lapack_matinv(bt%nb,bt%vinv,ierr)
IF(ierr/=0)CALL oft_abort('Bernstein table inversion failed','bern_build',__FILE__)
!---Monomials u^a2 v^a3 on the uniform lattice
DO m=1,bt%nb
  DO n=1,bt%nb
    v(m,n)=(REAL(bt%mi(2,m),r8)/REAL(p,r8))**bt%mi(2,n)*(REAL(bt%mi(3,m),r8)/REAL(p,r8))**bt%mi(3,n)
  END DO
END DO
ALLOCATE(bt%mvinv,SOURCE=v(1:bt%nb,1:bt%nb))
CALL lapack_matinv(bt%nb,bt%mvinv,ierr)
IF(ierr/=0)CALL oft_abort('Monomial table inversion failed','bern_build',__FILE__)
bt%cbound=(MAXVAL(SUM(ABS(bt%vinv),DIM=2))-1.d0)/2.d0+1.d-12
END SUBROUTINE bern_build
END MODULE oft_gs_cutcell
