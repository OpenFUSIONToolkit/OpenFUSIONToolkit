!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file axi_green.F90
!
!> Object and supporting functions for axisymmetric coil sets.
!!
!! @authors Chris Hansen
!! @date March 2014
!! @ingroup doxy_oft_physics
!---------------------------------------------------------------------------
MODULE axi_green
USE oft_base
IMPLICIT NONE
PRIVATE
!------------------------------------------------------------------------------
! TYPE axi_coil_set
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
TYPE, PUBLIC :: axi_coil_set
  INTEGER(i4) :: ncoils = 0
  REAL(r8) :: curr = 0.d0
  REAL(r8), POINTER, DIMENSION(:) :: scale => NULL()
  REAL(r8), POINTER, DIMENSION(:,:) :: pt => NULL()
END TYPE axi_coil_set
REAL(r8), PARAMETER :: ROFF = 1.d-13
PUBLIC green, grad_green, decay_eigenmodes, green_brute
CONTAINS
!------------------------------------------------------------------------------
! FUNCTION rf
!------------------------------------------------------------------------------
!> Create flux function object from definition file
!!
!! Computes Carlson's elliptic integral of the first kind, RF(x; y; z). x, y, and z must be
!! nonnegative, and at most one can be zero. TINY must be at least 5 times the machine
!! undeflow limit, BIG at most one fifth the machine overflow limit.
!------------------------------------------------------------------------------
FUNCTION rf(x,y,z)
REAL(r8), INTENT(in) :: x,y,z
REAL(r8) :: rf,alamb,ave,delx,dely,delz,e2
REAL(r8) :: e3,sqrtx,sqrty,sqrtz,xt,yt,zt
REAL(r8), PARAMETER :: ERRTOL=0.08d0
REAL(r8), PARAMETER :: loc_tiny=1.d-38
REAL(r8), PARAMETER :: loc_big=3.d37
! REAL(r8), PARAMETER :: loc_tiny=TINY(0.d0)*6.d0
! REAL(r8), PARAMETER :: loc_big=HUGE(0.d0)/6.d0
REAL(r8), PARAMETER :: loc_third=1.d0/3.d0
REAL(r8), PARAMETER :: C1=1.d0/24.d0
REAL(r8), PARAMETER :: C2=.1d0
REAL(r8), PARAMETER :: C3=3.d0/44.d0
REAL(r8), PARAMETER :: C4=1.d0/14.d0
IF((MIN(x,y,z).LT.0.d0).OR.(MIN(x+y,x+z,y+z).LT.loc_tiny).OR.(MAX(x,y,z).GT.loc_big))THEN
  CALL oft_abort('invalid arguments','rf',__FILE__)
END IF
xt=x
yt=y
zt=z
!---
DO
  sqrtx=SQRT(xt)
  sqrty=SQRT(yt)
  sqrtz=SQRT(zt)
  alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
  xt=0.25d0*(xt+alamb)
  yt=0.25d0*(yt+alamb)
  zt=0.25d0*(zt+alamb)
  ave=loc_third*(xt+yt+zt)
  delx=(ave-xt)/ave
  dely=(ave-yt)/ave
  delz=(ave-zt)/ave
  IF(MAX(ABS(delx),ABS(dely),ABS(delz)).LE.ERRTOL)EXIT
END DO
e2=delx*dely-delz**2
e3=delx*dely*delz
rf=(1.d0+(C1*e2-C2-C3*e3)*e2+C4*e3)/SQRT(ave)
END FUNCTION rf
!------------------------------------------------------------------------------
! FUNCTION rd
!------------------------------------------------------------------------------
!> Create flux function object from definition file
!!
!! Computes Carlson's elliptic integral of the second kind, RD(x; y; z). x and y must be
!! nonnegative, and at most one can be zero. z must be positive. TINY must be at least twice
!! the negative 2/3 power of the machine overflow limit. BIG must be at most 0:1ERRTOL
!! times the negative 2/3 power of the machine underflow limit.
!------------------------------------------------------------------------------
FUNCTION rd(x,y,z)
REAL(r8), INTENT(in) :: x,y,z
REAL(r8) :: rd,alamb,ave,delx,dely,delz,ea,eb
REAL(r8) :: ec,ed,ee,fac,sqrtx,sqrty,sqrtz,loc_sum,xt,yt,zt
REAL(r8), PARAMETER :: ERRTOL=0.05d0
REAL(r8), PARAMETER :: loc_tiny=1.d-25
REAL(r8), PARAMETER :: loc_big=4.5d21
! REAL(r8), PARAMETER :: loc_tiny=3.d0*(HUGE(0.d0)**-2.d0/3.d0)
! REAL(r8), PARAMETER :: loc_big=ERRTOL*(TINY(0.d0)**-2.d0/3.d0)/1.5d0
REAL(r8), PARAMETER :: C1=3.d0/14.d0
REAL(r8), PARAMETER :: C2=1.d0/6.d0
REAL(r8), PARAMETER :: C3=9.d0/22.d0
REAL(r8), PARAMETER :: C4=3.d0/26.d0
REAL(r8), PARAMETER :: C5=0.25d0*C3
REAL(r8), PARAMETER :: C6=1.5d0*C4
IF((MIN(x,y).LT.0.d0).OR.(MIN(x+y,z).LT.loc_tiny).OR.(MAX(x,y,z).GT.loc_big))THEN
  CALL oft_abort('invalid arguments','rd',__FILE__)
END IF
xt=x
yt=y
zt=z
loc_sum=0.
fac=1.
!---
DO
  sqrtx=sqrt(xt)
  sqrty=sqrt(yt)
  sqrtz=sqrt(zt)
  alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
  loc_sum=loc_sum+fac/(sqrtz*(zt+alamb))
  fac=0.25d0*fac
  xt=0.25d0*(xt+alamb)
  yt=0.25d0*(yt+alamb)
  zt=0.25d0*(zt+alamb)
  ave=0.2d0*(xt+yt+3.0d0*zt)
  delx=(ave-xt)/ave
  dely=(ave-yt)/ave
  delz=(ave-zt)/ave
  IF(MAX(ABS(delx),ABS(dely),ABS(delz)).LE.ERRTOL)EXIT
END DO
ea=delx*dely
eb=delz*delz
ec=ea-eb
ed=ea-6.d0*eb
ee=ed+ec+ec
rd=3.d0*loc_sum+fac*(1.d0+ed*(-C1+C5*ed-C6*delz*ee) &
  + delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*SQRT(ave))
END FUNCTION rd
!------------------------------------------------------------------------------
! FUNCTION ellf
!------------------------------------------------------------------------------
!> Create flux function object from definition file
!!
!! Legendre elliptic integral of the 1st kind F(; k), evaluated using Carlson's function RF.
!! The argument ranges are 0 =2, 0 ksin 1.
!------------------------------------------------------------------------------
FUNCTION ellf(phi,ak)
REAL(r8), INTENT(in) :: phi,ak
REAL(r8) :: ellf,s
s=sin(phi)
ellf=s*rf(cos(phi)**2,(1.d0-s*ak)*(1.d0+s*ak),1.d0)
END FUNCTION ellf
!------------------------------------------------------------------------------
! FUNCTION elle
!------------------------------------------------------------------------------
!> Create flux function object from definition file
!!
!! Legendre elliptic integral of the 2nd kind E(; j), evaluated using Carlson's functions RD
!! and RF. The argument ranges are 0 =2, 0 ksin 1.
!------------------------------------------------------------------------------
FUNCTION elle(phi,ak)
REAL(r8), INTENT(in) :: phi,ak
REAL(r8) :: elle,cc,q,s
s=sin(phi)
cc=cos(phi)**2
q=(1.d0-s*ak)*(1.d0+s*ak)
elle=s*(rf(cc,q,1.d0)-((s*ak)**2)*rd(cc,q,1.d0)/3.d0)
END FUNCTION elle
!------------------------------------------------------------------------------
! FUNCTION green
!------------------------------------------------------------------------------
!> Create flux function object from definition file
!!
!! Legendre elliptic integral of the 2nd kind E(; k), evaluated using Carlson's functions RD
!! and RF. The argument ranges are 0    =2, 0  ksin  1.
!------------------------------------------------------------------------------
FUNCTION green(r,z,rc,zc)
REAL(r8), intent(in) :: r,z,rc,zc
REAL(r8) :: green,k,k2,ellipk,ellipe,kextrap
IF(r < ROFF)THEN
  green = 0.d0
  RETURN
END IF
k2 = 4.d0*r*rc/((r + rc)*(r + rc) + (z - zc)*(z - zc))
kextrap=-1.d0
IF(k2 > (1.d0 - ROFF))THEN
  kextrap = ((r-rc)**2 + (z-zc)**2)/(4.d0*rc**2 + 4.d0*rc*(r-rc)+ (r-rc)**2 + (z-zc)**2)
  k2 = 1.d0 - ROFF
END IF
k=SQRT(k2)
ellipk = ellf(pi/2.d0,k);
ellipe = elle(pi/2.d0,k);
green = -SQRT(r*rc/k2)*((2.d0 - k2)*ellipk-2.d0*ellipe)/(2.d0*pi)
IF(kextrap>0.d0)green=green*LOG(kextrap)/LOG(1.d0-k2)
END FUNCTION green
!------------------------------------------------------------------------------
! SUBROUTINE grad_green
!------------------------------------------------------------------------------
!> Create flux function object from definition file
!!
!! Legendre elliptic integral of the 2nd kind E(; k), evaluated using Carlson's functions RD
!! and RF. The argument ranges are 0 =2, 0 < ksin < 1.
!------------------------------------------------------------------------------
SUBROUTINE grad_green(r,z,rc,zc,Fg,Gg)
REAL(r8), intent(in) :: r,z,rc,zc
REAL(r8), INTENT(out) :: Fg,Gg(2)
REAL(r8) :: k,k2,m1,ellipk,ellipe,dke
REAL(r8), PARAMETER :: dx = 1.d-6
! Fg=green(r,z,rc,zc)
! Gg(1)=(green(r+dx,z,rc,zc)-green(r-dx,z,rc,zc))/(dx*2.d0)
! Gg(2)=(green(r,z+dx,rc,zc)-green(r,z-dx,rc,zc))/(dx*2.d0)
k2 = 4.d0*r*rc/((r + rc)*(r + rc) + (z - zc) * (z - zc))
IF(k2 > (1.0 - ROFF))k2 = 1.d0-ROFF
m1 = 1.d0 - k2
k = sqrt(k2)
ellipk = ellf(pi/2.d0, k)
ellipe = elle(pi/2.d0, k)
Fg = -sqrt(r*rc/k2)*((2.d0 - k2)*ellipk - 2.d0*ellipe)
dKE = ellipe/m1 - ellipk
Gg(1) = 0.25d0*(Fg*k2*(r + rc)/(r*rc) - sqrt(k2/(r*rc))*(2.d0*rc-k2*(r + rc))*dKE)/(2.d0*pi)
Gg(2) = 0.25d0*k2*((z - zc)/(r*rc))*(Fg + sqrt(k2*r*rc)*dKE)/(2.d0*pi)
Fg=Fg/(2.d0*pi)
END SUBROUTINE grad_green
!------------------------------------------------------------------------------
! FUNCTION green
!------------------------------------------------------------------------------
!> Create flux function object from definition file
!!
!! Legendre elliptic integral of the 2nd kind E(; k), evaluated using Carlson's functions RD
!! and RF. The argument ranges are 0    =2, 0  ksin  1.
!------------------------------------------------------------------------------
FUNCTION green_brute(r,z,rc,zc)
REAL(r8), intent(in) :: r,z,rc,zc
REAL(r8) :: green_brute,phi,x(3),xc(3),dx
INTEGER :: i
x=[r,0.d0,z]
dx=2.d0*pi*rc/REAL(360,8)
green_brute=0.d0
DO i=1,360
  phi=i*2.d0*pi/REAL(360,8)
  xc=[rc*COS(phi),rc*SIN(phi),zc]
  green_brute=green_brute-dx*COS(phi)/SQRT(SUM((x-xc)**2))
END DO
green_brute=green_brute*r/(4.d0*pi)
END FUNCTION green_brute
!------------------------------------------------------------------------------
! SUBROUTINE decay_eigenmodes
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE decay_eigenmodes(ncoils,rc,eig_val,eig_vec,eta)
INTEGER(i4), INTENT(in) :: ncoils
REAL(r8), INTENT(in) :: rc(2,ncoils)
REAL(r8), intent(out) :: eig_val(ncoils),eig_vec(ncoils,ncoils)
REAL(r8), OPTIONAL, INTENT(in) :: eta(ncoils)
!---
INTEGER(i4) :: i,j,info,N,LDVL,LDVR,LDA,LWORK
REAL(r8), ALLOCATABLE, DIMENSION(:) :: WR,WI,WORK
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: Amat,VL,VR
CHARACTER(LEN=1) :: JOBVL,JOBVR
!--- Create coupling matrix
ALLOCATE(Amat(ncoils,ncoils))
!$omp parallel do private(j)
DO i=1,ncoils
  DO j=1,ncoils
    Amat(i,j)=green(rc(1,i),rc(2,i),rc(1,j),rc(2,j))/rc(1,i)
  END DO
END DO
!---Set eta profile
IF(PRESENT(eta))THEN
  DO i=1,ncoils
    Amat(i,:)=Amat(i,:)/eta(i)
  END DO
END IF
!--- Compute eigenvalues
JOBVL = 'V'
JOBVR = 'V'
N = ncoils
LDA = ncoils
LDVL = ncoils
LDVR = ncoils
ALLOCATE(WR(N),WI(N),VL(LDVL,N),VR(LDVR,N),WORK(1))
LWORK=-1
CALL DGEEV(JOBVL, JOBVR, N, Amat, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
LWORK=INT(WORK(1),4)
DEALLOCATE(WORK)
ALLOCATE(WORK(LWORK))
CALL DGEEV(JOBVL, JOBVR, N, Amat, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
eig_val=-1.d0/WR
eig_vec=REAL(VR,8)
DEALLOCATE(WI,WR,VL,VR,WORK,Amat)
END SUBROUTINE decay_eigenmodes
END MODULE axi_green
