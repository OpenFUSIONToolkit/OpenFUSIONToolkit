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
PUBLIC gs_qgeom_cutcell, gs_sauter_cutcell, gs_fsa_cutcell
CONTAINS
!------------------------------------------------------------------------------
!> Geometric factor \f$ g = q/F = \frac{1}{2\pi} \oint \frac{dl}{R |\nabla \psi|} \f$ on flux surfaces
!!
!! Same contract and units as @ref oft_gs::gs_get_qprof with `geom_only`. Only arcs inside
!! the plasma mask (@ref oft_gs::gs_test_bounds) contribute; surfaces whose winding about the
!! O-point is not one return `g=0`. Summation order is fixed, so results do not depend on
!! the thread count.
!------------------------------------------------------------------------------
subroutine gs_qgeom_cutcell(gseq,nr,psi_q,g)
class(gs_equil), intent(inout) :: gseq !< G-S object
integer(i4), intent(in) :: nr !< Number of surfaces
real(r8), intent(in) :: psi_q(nr) !< Surface locations (0 at LCFS, 1 at axis)
real(r8), intent(out) :: g(nr) !< Geometric factor at each surface
real(r8) :: wind(nr),dpsi
integer(i4) :: i
CHARACTER(LEN=OFT_ERROR_SLEN) :: error_str
CALL geom_sums(gseq,nr,psi_q,g,wind,dpsi)
IF(dpsi==0.d0)RETURN
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
end subroutine gs_qgeom_cutcell
!------------------------------------------------------------------------------
!> \f$ q/F \f$, the flux-surface averages and the surface shape of @ref oft_gs::gs_get_qprof
!! (`ravgs`, `fsa_avgs`, `shape_geo`) by the same quadrature
!------------------------------------------------------------------------------
subroutine gs_fsa_cutcell(gseq,nr,psi_q,qf,ravgs,ok,fpol,fsa,shp)
class(gs_equil), intent(inout) :: gseq !< G-S object
integer(i4), intent(in) :: nr !< Number of surfaces
real(r8), intent(in) :: psi_q(nr) !< Surface locations (0 at LCFS, 1 at axis)
real(r8), intent(out) :: qf(nr) !< \f$ q/F \f$
real(r8), intent(out) :: ravgs(nr,4) !< <R>, <1/R>, <1/R^2>, dV/dpsi
logical, intent(out) :: ok(nr) !< Surface inside (0,1) and closed about the O-point
real(r8), optional, intent(in) :: fpol(nr) !< F on each surface (with `fsa`)
real(r8), optional, intent(out) :: fsa(nr,4) !< <|grad psi|>, <|grad psi|^2>, <B_p^2>, <1/B^2>
real(r8), optional, intent(out) :: shp(nr,6) !< R_min, R_max, Z_min, Z_max, R(Z_min), R(Z_max)
real(r8) :: wind(nr),racc(3,nr),dpsi
real(r8), allocatable :: facc(:,:),xacc(:,:,:)
integer(i4) :: k
IF(PRESENT(fsa).AND.(.NOT.PRESENT(fpol)))CALL oft_abort('"fsa" requires "fpol"','gs_fsa_cutcell',__FILE__)
!---Unallocated accumulators are passed as absent (F2008)
IF(PRESENT(fsa))ALLOCATE(facc(5,nr))
IF(PRESENT(shp))ALLOCATE(xacc(2,4,nr))
CALL geom_sums(gseq,nr,psi_q,qf,wind,dpsi,racc,fpol,facc,xacc)
ok=(dpsi/=0.d0).AND.(psi_q>0.d0).AND.(psi_q<1.d0).AND.(ABS(wind-1.d0)<=wind_tol)
ravgs=0.d0
WHERE(ok)
  ravgs(:,1)=racc(2,:)/racc(1,:)
  ravgs(:,2)=racc(3,:)/racc(1,:)
  ravgs(:,3)=qf/racc(1,:)
  ravgs(:,4)=-2.d0*pi*racc(1,:)/ABS(dpsi)
  qf=qf/(2.d0*pi*ABS(dpsi))
ELSEWHERE
  qf=0.d0
END WHERE
IF(PRESENT(fsa))THEN
  fsa=0.d0
  DO k=1,4
    WHERE(ok)fsa(:,k)=facc(k+1,:)/facc(1,:)
  END DO
END IF
IF(PRESENT(shp))THEN
  shp=0.d0
  WHERE(ok)
    shp(:,1)=-xacc(1,1,:); shp(:,2)=xacc(1,2,:)
    shp(:,3)=-xacc(1,3,:); shp(:,4)=xacc(1,4,:)
    shp(:,5)=xacc(2,3,:); shp(:,6)=xacc(2,4,:)
  END WHERE
END IF
end subroutine gs_fsa_cutcell
!------------------------------------------------------------------------------
!> Fixed-order sums of \f$ \oint w/R \f$, winding and optionally \f$ \oint w\,[R, R^2, 1] \f$,
!! \f$ w = dl/|\nabla \hat{\psi}| \f$; all zero with `dpsi=0` if there is no plasma
!!
!! Surface extrema (`xacc`) are reduced with a total order (@ref ext_merge), which is exact,
!! so they too do not depend on the thread count or schedule.
!------------------------------------------------------------------------------
subroutine geom_sums(gseq,nr,psi_q,g,wind,dpsi,racc,fpol,facc,xacc)
class(gs_equil), intent(inout) :: gseq !< G-S object
integer(i4), intent(in) :: nr !< Number of surfaces
real(r8), intent(in) :: psi_q(nr) !< Surface locations (0 at LCFS, 1 at axis)
real(r8), intent(out) :: g(nr) !< \f$ \oint w/R \f$
real(r8), intent(out) :: wind(nr) !< Winding number about the O-point
real(r8), intent(out) :: dpsi !< Axis minus LCFS flux
real(r8), optional, intent(out) :: racc(3,nr) !< \f$ \oint w\,[R, R^2, 1] \f$
real(r8), optional, intent(in) :: fpol(nr) !< F on each surface (with `facc`)
real(r8), optional, intent(out) :: facc(5,nr) !< \f$ \oint w R\,[1, |\nabla\psi|, |\nabla\psi|^2, B_p^2, 1/B^2] \f$
real(r8), optional, intent(out) :: xacc(2,4,nr) !< Extrema of \f$ [-R, R, -Z, Z] \f$ with the other coordinate (see @ref ext_add)
type(bern_tab) :: bt
real(r8), pointer :: vals(:)
real(r8), allocatable :: gblk(:,:),wblk(:,:),rblk(:,:,:),fblk(:,:,:),gloc(:),wloc(:),rloc(:,:),floc(:,:),xloc(:,:,:)
real(r8) :: x1
integer(i4) :: nblk,ib,cell,nc,i
g=0.d0; wind=0.d0; dpsi=0.d0
IF(PRESENT(racc))racc=0.d0
IF(PRESENT(facc))facc=0.d0
IF(PRESENT(xacc))xacc=-HUGE(1.d0)
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
ALLOCATE(gblk(nr,nblk),wblk(nr,nblk),rblk(3,nr,MERGE(nblk,0,PRESENT(racc))),fblk(5,nr,MERGE(nblk,0,PRESENT(facc))))
!---Thread-private accumulators; unallocated ones are passed to `cell_contrib` as absent (F2008)
!$omp parallel private(ib,cell,i,gloc,wloc,rloc,floc,xloc)
ALLOCATE(gloc(nr),wloc(nr))
IF(PRESENT(racc))ALLOCATE(rloc(3,nr))
IF(PRESENT(facc))ALLOCATE(floc(5,nr))
IF(PRESENT(xacc))THEN
  ALLOCATE(xloc(2,4,nr))
  xloc=-HUGE(1.d0)
END IF
!$omp do schedule(dynamic,1)
DO ib=1,nblk
  gloc=0.d0; wloc=0.d0
  IF(PRESENT(racc))rloc=0.d0
  IF(PRESENT(facc))floc=0.d0
  DO cell=(ib-1)*nblk_cells+1,MIN(ib*nblk_cells,nc)
    IF(gseq%device%fe_rep%mesh%reg(cell)/=1)CYCLE
    CALL cell_contrib(gseq,bt,vals,x1,dpsi,nr,psi_q,cell,gloc,wloc,fpol=fpol,racc=rloc,facc=floc,xacc=xloc)
  END DO
  gblk(:,ib)=gloc; wblk(:,ib)=wloc
  IF(PRESENT(racc))rblk(:,:,ib)=rloc
  IF(PRESENT(facc))fblk(:,:,ib)=floc
END DO
!$omp end do
IF(PRESENT(xacc))THEN
!$omp critical
  DO i=1,nr
    CALL ext_merge(xacc(:,:,i),xloc(:,:,i))
  END DO
!$omp end critical
END IF
DEALLOCATE(gloc,wloc)
IF(ALLOCATED(rloc))DEALLOCATE(rloc)
IF(ALLOCATED(floc))DEALLOCATE(floc)
IF(ALLOCATED(xloc))DEALLOCATE(xloc)
!$omp end parallel
!---Fixed-order reduction
DO ib=1,nblk
  g=g+gblk(:,ib)
  wind=wind+wblk(:,ib)
  IF(PRESENT(racc))racc=racc+rblk(:,:,ib)
  IF(PRESENT(facc))facc=facc+fblk(:,:,ib)
END DO
DEALLOCATE(vals,gblk,wblk,rblk,fblk)
end subroutine geom_sums
!------------------------------------------------------------------------------
!> Flux-surface sums for @ref grad_shaf_bootstrap::sauter_fc by the same quadrature
!!
!! `acc(:,i)` holds \f$ \oint w\,[R, R^2, 1, R a, R|B|, R B^2, 1/R, R f(|B|/B_{ref})] \f$ with
!! \f$ w = dl/|\nabla \hat{\psi}| \f$, \f$ a \f$ the distance to the O-point and
!! \f$ f(b)=(1-\sqrt{1-b}(1+b/2))/b^2 \f$. \f$ B_{max} \f$ is sampled at arc ends and Gauss
!! points and refined by a parabola on each arc; \f$ R_{max} \f$ and \f$ R_{min} \f$ are exact on the
!! FE surface (@ref ext_add). Results do not depend on the thread count.
!------------------------------------------------------------------------------
subroutine gs_sauter_cutcell(gseq,nr,psi_q,fpol,bref,acc,ext,ok)
class(gs_equil), intent(inout) :: gseq !< G-S object
integer(i4), intent(in) :: nr !< Number of surfaces
real(r8), intent(in) :: psi_q(nr) !< Surface locations (0 at LCFS, 1 at axis)
real(r8), intent(in) :: fpol(nr) !< F on each surface
real(r8), intent(in) :: bref(nr) !< \f$ B_{ref} \f$ for the last sum (skipped if <=0)
real(r8), intent(out) :: acc(8,nr) !< Surface sums
real(r8), intent(out) :: ext(3,nr) !< \f$ B_{max} \f$, \f$ R_{max} \f$, \f$ -R_{min} \f$
logical, intent(out) :: ok(nr) !< Surface inside (0,1) and closed about the O-point
type(bern_tab) :: bt
real(r8), pointer :: vals(:)
real(r8), allocatable :: ablk(:,:,:),wblk(:,:),gloc(:),wloc(:),aloc(:,:),eloc(:),xloc(:,:,:),wind(:)
real(r8) :: x1,dpsi,xacc(2,4,nr)
integer(i4) :: nblk,ib,cell,i,nc
acc=0.d0; ext=-HUGE(1.d0); ok=.FALSE.; xacc=-HUGE(1.d0)
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
ALLOCATE(ablk(8,nr,nblk),wblk(nr,nblk))
!$omp parallel private(ib,cell,i,gloc,wloc,aloc,eloc,xloc)
ALLOCATE(gloc(nr),wloc(nr),aloc(8,nr),eloc(nr),xloc(2,4,nr))
eloc=-HUGE(1.d0); xloc=-HUGE(1.d0)
!$omp do schedule(dynamic,1)
DO ib=1,nblk
  gloc=0.d0; wloc=0.d0; aloc=0.d0
  DO cell=(ib-1)*nblk_cells+1,MIN(ib*nblk_cells,nc)
    IF(gseq%device%fe_rep%mesh%reg(cell)/=1)CYCLE
    CALL cell_contrib(gseq,bt,vals,x1,dpsi,nr,psi_q,cell,gloc,wloc,fpol,bref,aloc,eloc,xacc=xloc)
  END DO
  ablk(:,:,ib)=aloc; wblk(:,ib)=wloc
END DO
!$omp end do
!$omp critical
ext(1,:)=MAX(ext(1,:),eloc)
DO i=1,nr
  CALL ext_merge(xacc(:,:,i),xloc(:,:,i))
END DO
!$omp end critical
DEALLOCATE(gloc,wloc,aloc,eloc,xloc)
!$omp end parallel
ext(2,:)=xacc(1,2,:); ext(3,:)=xacc(1,1,:)
ALLOCATE(wind(nr))
wind=0.d0
DO ib=1,nblk
  acc=acc+ablk(:,:,ib)
  wind=wind+wblk(:,ib)
END DO
DO i=1,nr
  ok(i)=(psi_q(i)>0.d0).AND.(psi_q(i)<1.d0).AND.(ABS(wind(i)-1.d0)<=wind_tol)
END DO
DEALLOCATE(vals,ablk,wblk,wind)
end subroutine gs_sauter_cutcell
!------------------------------------------------------------------------------
!> Add arc integrals and winding (and optionally geometric or Sauter sums) of all surfaces crossing `cell`
!------------------------------------------------------------------------------
subroutine cell_contrib(gseq,bt,vals,x1,dpsi,nr,psi_q,cell,gloc,wloc,fpol,bref,aloc,eloc,racc,facc,xacc)
class(gs_equil), intent(inout) :: gseq !< G-S object
type(bern_tab), intent(in) :: bt !< Bernstein tables
real(r8), intent(in) :: vals(:) !< Local \f$ \psi \f$ values
real(r8), intent(in) :: x1 !< LCFS flux
real(r8), intent(in) :: dpsi !< Axis minus LCFS flux
integer(i4), intent(in) :: nr !< Number of surfaces
real(r8), intent(in) :: psi_q(nr) !< Surface locations
integer(i4), intent(in) :: cell !< Cell index
real(r8), intent(inout) :: gloc(nr) !< Accumulated \f$ \oint dl/(R|\nabla\hat{\psi}|) \f$
real(r8), intent(inout) :: wloc(nr) !< Accumulated winding number
real(r8), optional, intent(in) :: fpol(nr) !< F on each surface (Sauter sums only)
real(r8), optional, intent(in) :: bref(nr) !< \f$ B_{ref} \f$ on each surface (Sauter sums only)
real(r8), optional, intent(inout) :: aloc(:,:) !< Accumulated Sauter sums [8,nr] (see @ref gs_sauter_cutcell)
real(r8), optional, intent(inout) :: eloc(:) !< Accumulated \f$ B_{max} \f$ [nr]
real(r8), optional, intent(inout) :: racc(:,:) !< Accumulated \f$ \oint w\,[R, R^2, 1] \f$ [3,nr] (see @ref geom_sums)
real(r8), optional, intent(inout) :: facc(:,:) !< Accumulated oint w R [1, |grad psi|, |grad psi|^2, B_p^2, 1/B^2] [5,nr] (needs `fpol`)
real(r8), optional, intent(inout) :: xacc(:,:,:) !< Accumulated surface extrema [2,4,nr] (see @ref ext_add)
type(oft_scalar_bfem), pointer :: lag_rep
integer(i4) :: jdofs(32),jc,k,ns,d0,i
integer(i4) :: dep(3*max_depth+4)
real(r8) :: lat(32),cc(32),cm(32),cs(32),tri(2,3),stk(2,3,3*max_depth+4),xv(2,3),detm0
real(r8) :: lo,hi,pad,dvec(2),dmu(3),ainv(2,2),m(2,3),dratio
logical :: curved,act(nr),sact(nr,3*max_depth+4)
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
      CALL tri_arc(psi_q(i),gloc(i),wloc(i),i)
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
subroutine tri_arc(s,gacc,wacc,iq)
real(r8), intent(in) :: s
real(r8), intent(inout) :: gacc,wacc
integer(i4), intent(in) :: iq
integer(i4) :: k0,k1,k2,ig,kk,nab,ng
logical :: ab(3)
real(r8) :: p1(2),p2(2),e(2),xi1,xi2,xi,x(2),val,gr(2)
real(r8) :: pt(2),detm,sum,q1(2),q2(2),dth,sig,gmid(2),xg(5),wg(5)
real(r8) :: av(8),rv(3),fv(5),gq,wq,bq,b,sx(0:6),sb(0:6),xs(0:6),sp(2,0:6),sg(2,0:6)
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
av=0.d0; rv=0.d0; fv=0.d0
IF(dratio>=smooth_ratio)THEN
  ng=3; xg(1:3)=xg3; wg(1:3)=wg3
ELSE
  ng=5; xg=xg5; wg=wg5
END IF
IF(PRESENT(aloc))THEN
  sb(0)=b_mag(p1,iq); sx(0)=xi1
  sb(ng+1)=b_mag(p2,iq); sx(ng+1)=xi2
END IF
IF(PRESENT(xacc))THEN
  CALL poly_eval(bt,cm,p1,val,gr); sg(:,0)=phys_grad(p1,gr)
  CALL cell_geom(p1,sp(:,0),detm); xs(0)=xi1
  CALL poly_eval(bt,cm,p2,val,gr); sg(:,ng+1)=phys_grad(p2,gr)
  CALL cell_geom(p2,sp(:,ng+1),detm); xs(ng+1)=xi2
END IF
DO ig=1,ng
  xi=(xi1+xi2)/2.d0+xg(ig)*(xi2-xi1)/2.d0
  x=arc_point(s,p1,p2,xi1,xi2,xi)
  CALL poly_eval(bt,cm,x,val,gr)
  CALL cell_geom(x,pt,detm)
  IF(PRESENT(xacc))THEN
    xs(ig)=xi; sp(:,ig)=pt; sg(:,ig)=phys_grad(x,gr)
  END IF
  sum=sum+wg(ig)*ABS(detm)/(pt(1)*ABS(DOT_PRODUCT(gr,dvec)))
  IF(ig==(ng+1)/2)gmid=gr
  wq=wg(ig)*ABS(detm)/ABS(DOT_PRODUCT(gr,dvec))
  IF(PRESENT(racc))rv=rv+wq*[pt(1),pt(1)**2,1.d0]
  IF(PRESENT(facc))THEN
    gq=grad_mag(x,gr)
    fv=fv+wq*pt(1)*[1.d0,gq,gq**2,(gq/pt(1))**2,pt(1)**2/(gq**2+fpol(iq)**2)]
  END IF
  IF(PRESENT(aloc))THEN
    bq=b_mag(x,iq,gr)
    sx(ig)=xi; sb(ig)=bq
    av(1)=av(1)+wq*pt(1)
    av(2)=av(2)+wq*pt(1)**2
    av(3)=av(3)+wq
    av(4)=av(4)+wq*pt(1)*SQRT((pt(1)-gseq%o_point(1))**2+(pt(2)-gseq%o_point(2))**2)
    av(5)=av(5)+wq*pt(1)*bq
    av(6)=av(6)+wq*pt(1)*bq**2
    av(7)=av(7)+wq/pt(1)
    IF(bref(iq)>0.d0)THEN
      b=MIN(1.d0,bq/bref(iq))
      av(8)=av(8)+wq*pt(1)*(1.d0-SQRT(1.d0-b)*(1.d0+b/2.d0))/b**2
    END IF
  END IF
END DO
gacc=gacc+sum*ABS(xi2-xi1)/2.d0
IF(PRESENT(racc))racc(:,iq)=racc(:,iq)+rv*ABS(xi2-xi1)/2.d0
IF(PRESENT(facc))facc(:,iq)=facc(:,iq)+fv*ABS(xi2-xi1)/2.d0
IF(PRESENT(aloc))THEN
  aloc(:,iq)=aloc(:,iq)+av*ABS(xi2-xi1)/2.d0
  eloc(iq)=MAX(eloc(iq),arc_peak(ng,sx,sb))
END IF
IF(PRESENT(xacc))CALL arc_extrema(s,p1,p2,xi1,xi2,ng,xs,sp,sg,xacc(:,:,iq))
!---Winding about O-point, oriented with 1-psi-hat increasing to the left
CALL cell_geom(p1,q1,detm)
CALL cell_geom(p2,q2,detm)
q1=q1-gseq%o_point; q2=q2-gseq%o_point
dth=ATAN2(q1(1)*q2(2)-q1(2)*q2(1),DOT_PRODUCT(q1,q2))
sig=SIGN(1.d0,(p2(1)-p1(1))*gmid(2)-(p2(2)-p1(2))*gmid(1))*SIGN(1.d0,detm)
wacc=wacc+sig*dth/(2.d0*pi)
end subroutine tri_arc
!---Point of surface s on the current arc at xi (projected from the chord p1-p2 along d)
function arc_point(s,p1,p2,xi1,xi2,xi) result(x)
real(r8), intent(in) :: s,p1(2),p2(2),xi1,xi2,xi
real(r8) :: x(2),c(2),muc(3),elo,ehi,eta
integer(i4) :: kk
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
end function arc_point
!---Extrema of R and Z on the current arc of surface s, added to `xa`
!!
!! Candidates: arc ends, Gauss points, and bracketed roots of \f$ \partial_Z \psi \f$ / \f$ \partial_R \psi \f$.
subroutine arc_extrema(s,p1,p2,xi1,xi2,ng,xs,sp,sg,xa)
real(r8), intent(in) :: s,p1(2),p2(2),xi1,xi2
integer(i4), intent(in) :: ng
real(r8), intent(in) :: xs(0:),sp(:,0:),sg(:,0:) !< Sample xi, (R,Z) and physical gradient
real(r8), intent(inout) :: xa(2,4)
real(r8) :: xl,xr,fl,fr,xm,fm,x(2),val,gr(2),gp(2),pt(2),detm
integer(i4) :: k,jg,it,side
DO k=0,ng+1
  CALL ext_add(xa,sp(:,k))
END DO
DO jg=1,2
  DO k=0,ng
    fl=sg(jg,k); fr=sg(jg,k+1)
    IF(.NOT.(((fl<0.d0).AND.(fr>0.d0)).OR.((fl>0.d0).AND.(fr<0.d0))))CYCLE
    xl=xs(k); xr=xs(k+1); side=0
    DO it=1,60
      xm=(xl*fr-xr*fl)/(fr-fl)
      x=arc_point(s,p1,p2,xi1,xi2,xm)
      CALL poly_eval(bt,cm,x,val,gr)
      gp=phys_grad(x,gr); fm=gp(jg)
      IF(fm==0.d0)EXIT
      IF((fm>0.d0).EQV.(fr>0.d0))THEN
        xr=xm; fr=fm
        IF(side==-1)fl=fl/2.d0
        side=-1
      ELSE
        xl=xm; fl=fm
        IF(side==1)fr=fr/2.d0
        side=1
      END IF
      IF(ABS(xr-xl)<=1.d-13*ABS(xi2-xi1))EXIT
    END DO
    CALL cell_geom(x,pt,detm)
    CALL ext_add(xa,pt)
  END DO
END DO
end subroutine arc_extrema
!---Physical gradient at uv, from the (u,v) gradient of 1-psi-hat
function phys_grad(uv,grp) result(gp)
real(r8), intent(in) :: uv(2),grp(2)
real(r8) :: gp(2),lam(3),gop(3,3),dj
lam=[1.d0-uv(1)-uv(2),uv(1),uv(2)]
CALL gseq%device%fe_rep%mesh%jacobian(cell,lam,gop,dj)
gp=grp(1)*gop(1:2,2)+grp(2)*gop(1:2,3)
end function phys_grad
!---|grad(psi)| at uv, from the (u,v) gradient of 1-psi-hat
function grad_mag(uv,grp) result(gq)
real(r8), intent(in) :: uv(2),grp(2)
real(r8) :: gq
gq=ABS(dpsi)*SQRT(SUM(phys_grad(uv,grp)**2))
end function grad_mag
!---|B| at uv on surface iq
function b_mag(uv,iq,grin) result(bq)
real(r8), intent(in) :: uv(2)
integer(i4), intent(in) :: iq
real(r8), optional, intent(in) :: grin(2)
real(r8) :: bq,lam(3),p3(3),grp(2),vv
IF(PRESENT(grin))THEN
  grp=grin
ELSE
  CALL poly_eval(bt,cm,uv,vv,grp)
END IF
lam=[1.d0-uv(1)-uv(2),uv(1),uv(2)]
p3=gseq%device%fe_rep%mesh%log2phys(cell,lam)
bq=SQRT(dpsi**2*SUM(phys_grad(uv,grp)**2)+fpol(iq)**2)/p3(1)
end function b_mag
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
!> Add point `pt` = (R,Z) to surface extrema `xa`
!!
!! Column k holds (key, companion) for keys \f$ [-R, R, -Z, Z] \f$ with companions
!! \f$ [Z, Z, R, R] \f$. The larger key wins and ties go to the larger companion: a total order
!! of exact comparisons, so the result does not depend on the order points are added.
!------------------------------------------------------------------------------
pure subroutine ext_add(xa,pt)
real(r8), intent(inout) :: xa(2,4)
real(r8), intent(in) :: pt(2)
CALL ext_upd(xa(:,1),-pt(1),pt(2))
CALL ext_upd(xa(:,2),pt(1),pt(2))
CALL ext_upd(xa(:,3),-pt(2),pt(1))
CALL ext_upd(xa(:,4),pt(2),pt(1))
end subroutine ext_add
!------------------------------------------------------------------------------
!> Merge surface extrema `xb` into `xa` (same order as @ref ext_add)
!------------------------------------------------------------------------------
pure subroutine ext_merge(xa,xb)
real(r8), intent(inout) :: xa(2,4)
real(r8), intent(in) :: xb(2,4)
integer(i4) :: k
DO k=1,4
  CALL ext_upd(xa(:,k),xb(1,k),xb(2,k))
END DO
end subroutine ext_merge
!------------------------------------------------------------------------------
!> Keep the larger of (`key`,`comp`) and `e` in lexicographic order
!------------------------------------------------------------------------------
pure subroutine ext_upd(e,key,comp)
real(r8), intent(inout) :: e(2)
real(r8), intent(in) :: key,comp
IF((key>e(1)).OR.((key==e(1)).AND.(comp>e(2))))THEN
  e(1)=key; e(2)=comp
END IF
end subroutine ext_upd
!------------------------------------------------------------------------------
!> Maximum of arc samples y(0:ng+1) at x(0:ng+1), refined by a parabola through the
!! largest interior sample and its neighbours
!------------------------------------------------------------------------------
pure real(r8) function arc_peak(ng,x,y)
integer(i4), intent(in) :: ng
real(r8), intent(in) :: x(0:),y(0:)
real(r8) :: d1,d2,a,b
integer(i4) :: k
k=MAXLOC(y(0:ng+1),DIM=1)-1
arc_peak=y(k)
IF((k<1).OR.(k>ng))RETURN
d1=(y(k)-y(k-1))/(x(k)-x(k-1)); d2=(y(k+1)-y(k))/(x(k+1)-x(k))
a=(d2-d1)/(x(k+1)-x(k-1))
IF(a>=0.d0)RETURN
b=d1+a*(x(k)-x(k-1)) ! slope at x(k)
arc_peak=MAX(y(k),y(k)-b**2/(4.d0*a))
end function arc_peak
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
