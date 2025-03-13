!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_lag_poly.F90
!
!> Lagrange interpolatory polynomials
!!
!! @author Chris Hansen
!! @date October 2018
!! @ingroup doxy_oft_grid
!------------------------------------------------------------------------------
MODULE oft_lag_poly
USE oft_local
IMPLICIT NONE
#include "local.h"
REAL(r8), PARAMETER :: lag_even_x1(2) = (/0.d0,1.d0/)
REAL(r8), PARAMETER :: lag_even_x2(3) = (/0.d0,1.d0,2.d0/)/2.d0
REAL(r8), PARAMETER :: lag_even_x3(4) = (/0.d0,1.d0,2.d0,3.d0/)/3.d0
REAL(r8), PARAMETER :: lag_even_x4(5) = (/0.d0,1.d0,2.d0,3.d0,4.d0/)/5.d0
CONTAINS
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function lag_1d_norm(ii,xnodes,n) result(val)
integer(i4), intent(in) :: ii,n
real(r8), intent(in) :: xnodes(n)
integer(i4) :: j
real(r8) :: val
val=1.d0
DO j=1,n
  IF(j/=ii)val=val*(xnodes(ii)-xnodes(j))
END DO
end function lag_1d_norm
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function lag_1d_unnorm(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii,n
real(r8), intent(in) :: x,xnodes(n)
integer(i4) :: j
real(r8) :: val
val=1.d0
DO j=1,n
  IF(j/=ii)val=val*(x-xnodes(j))
END DO
end function lag_1d_unnorm
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function lag_1d(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii,n
real(r8), intent(in) :: x,xnodes(n)
integer(i4) :: j
real(r8) :: val,norm
val=1.d0
norm=1.d0
DO j=1,n
  IF(j==ii)CYCLE
  val = val*(x-xnodes(j))
  norm = norm*(xnodes(ii)-xnodes(j))
END DO
val=val/norm
end function lag_1d
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function dlag_1d_unnorm(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii,n
real(r8), intent(in) :: x,xnodes(n)
integer(i4) :: i,j
real(r8) :: val,xi,vtmp
val=0.d0
xi=xnodes(ii)
DO i=1,n
  IF(i==ii)CYCLE
  vtmp=1.d0
  DO j=1,n
    IF(j==i.OR.j==ii)CYCLE
    vtmp=vtmp*(x-xnodes(j))
  END DO
  val=val+vtmp
END DO
end function dlag_1d_unnorm
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function dlag_1d(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii,n
real(r8), intent(in) :: x,xnodes(n)
integer(i4) :: i,j
real(r8) :: val,norm,vtmp
val=0.d0
norm=1.d0
DO i=1,n
  IF(i==ii)CYCLE
  norm = norm*(xnodes(ii)-xnodes(i))
  vtmp=1.d0
  DO j=1,n
    IF(j==i.OR.j==ii)CYCLE
    vtmp=vtmp*(x-xnodes(j))
  END DO
  val = val + vtmp
END DO
val = val/norm
end function dlag_1d
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function d2lag_1d(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii,n
real(r8), intent(in) :: x,xnodes(n)
integer(i4) :: i,j,k
real(r8) :: val,norm,vtmp1,vtmp2
val=0.d0
norm=1.d0
DO k=1,n
  IF(k==ii)CYCLE
  norm = norm*(xnodes(ii)-xnodes(k))
  vtmp2=0.d0
  DO i=1,n
    IF(i==k.OR.i==ii)CYCLE
    vtmp1=1.d0
    DO j=1,n
      IF(j==k.OR.j==i.OR.j==ii)CYCLE
      vtmp1=vtmp1*(x-xnodes(j))
    END DO
    vtmp2 = vtmp2 + vtmp1
  END DO
  val = val + vtmp2
END DO
val = val/norm
end function d2lag_1d
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function lag_2d(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii(2),n
real(r8), intent(in) :: x(2),xnodes(n)
integer(i4) :: i
real(r8) :: val,xi
val=1.d0
DO i=1,2
  val=val*lag_1d(ii(i),x(i),xnodes,n)
END DO
end function lag_2d
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function dlag_2d(ii,k,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii(2),k,n
real(r8), intent(in) :: x(2),xnodes(n)
integer(i4) :: i
real(r8) :: val,norm
val=1.d0
norm=1.d0
DO i=1,2
  norm=norm*lag_1d_norm(ii(i),xnodes,n)
  IF(i==k)THEN
    val=val*dlag_1d_unnorm(ii(i),x(i),xnodes,n)
  ELSE
    val=val*lag_1d_unnorm(ii(i),x(i),xnodes,n)
  END IF
END DO
val=val/norm
end function dlag_2d
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function d2lag_2d(ii,k,m,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii(2),k,m,n
real(r8), intent(in) :: x(2),xnodes(n)
integer(i4) :: i
real(r8) :: val,xi,vtmp
val=1.d0
DO i=1,2
  IF(i==k.OR.i==m)THEN
    IF(k==m)THEN
      val=val*d2lag_1d(ii(i),x(i),xnodes,n)
    ELSE
      val=val*dlag_1d(ii(i),x(i),xnodes,n)
    END IF
  ELSE
    val=val*lag_1d(ii(i),x(i),xnodes,n)
  END IF
END DO
end function d2lag_2d
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function lag_3d(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii(3),n
real(r8), intent(in) :: x(3),xnodes(n)
integer(i4) :: i
real(r8) :: val,xi
val=1.d0
DO i=1,3
  val=val*lag_1d(ii(i),x(i),xnodes,n)
END DO
end function lag_3d
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function dlag_3d(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii(3),n
real(r8), intent(in) :: x(3),xnodes(n)
integer(i4) :: i
real(r8) :: val(3),vtmp(3),norm
DO i=1,3
  norm=lag_1d_norm(ii(i),xnodes,n)
  vtmp(i)=lag_1d_unnorm(ii(i),x(i),xnodes,n)/norm
  val(i)=dlag_1d_unnorm(ii(i),x(i),xnodes,n)/norm
END DO
val(1)=vtmp(2)*vtmp(3)*val(1)
val(2)=vtmp(1)*vtmp(3)*val(2)
val(3)=vtmp(1)*vtmp(2)*val(3)
end function dlag_3d
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function d2lag_3d(ii,k,m,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii(3),k,m,n
real(r8), intent(in) :: x(3),xnodes(n)
integer(i4) :: i
real(r8) :: val
val=1.d0
! IF(k==m)val=d2lag_1d(ii(k),x(k),xnodes,n)
! DO i=1,3
!   IF(i==k.NEQV.i==m)THEN
!     val=val*dlag_1d(ii(i),x(i),xnodes,n)
!   ELSE
!     val=val*lag_1d(ii(i),x(i),xnodes,n)
!   END IF
! END DO
DO i=1,3
  IF(i==k.OR.i==m)THEN
    IF(k==m)THEN
      val=val*d2lag_1d(ii(i),x(i),xnodes,n)
    ELSE
      val=val*dlag_1d(ii(i),x(i),xnodes,n)
    END IF
  ELSE
    val=val*lag_1d(ii(i),x(i),xnodes,n)
  END IF
END DO
end function d2lag_3d
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function lag_1d_bary(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii,n
real(r8), intent(in) :: x(2),xnodes(n)
integer(i4) :: i,il(2)
real(r8) :: val
val=1.d0
il=(/ii,n-1-ii/)+1
DO i=1,2
  val=val*lag_1d(il(i), x(i), xnodes, il(i))
END DO
end function lag_1d_bary
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function dlag_1d_bary(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii,n
real(r8), intent(in) :: x(2),xnodes(n)
integer(i4) :: i,ix,il(2)
real(r8) :: val(2),norm
val=1.d0
norm=1.d0
il=(/ii,n-1-ii/)+1
DO ix=1,2
  norm=norm*lag_1d_norm(il(ix),xnodes,il(ix))
  DO i=1,2
    IF(i==ix)THEN
      val(ix)=val(ix)*dlag_1d_unnorm(il(i), x(i), xnodes, il(i))
    ELSE
      val(ix)=val(ix)*lag_1d_unnorm(il(i), x(i), xnodes, il(i))
    END IF
  END DO
END DO
val=val/norm
end function dlag_1d_bary
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function d2lag_1d_bary(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii,n
real(r8), intent(in) :: x(2),xnodes(n)
integer(i4) :: i,ix,iy,k,il(2)
real(r8) :: val(3)
val=1.d0
il=(/ii,n-1-ii/)+1
k=0
DO ix=1,2
  DO iy=ix,2
    k=k+1
    DO i=1,2
      IF(i==ix.OR.i==iy)THEN
        IF(ix==iy)THEN
          val(k)=val(k)*d2lag_1d(il(i), x(i), xnodes, il(i))
        ELSE
          val(k)=val(k)*dlag_1d(il(i), x(i), xnodes, il(i))
        END IF
      ELSE
        val(k)=val(k)*lag_1d(il(i), x(i), xnodes, il(i))
      END IF
    END DO
  END DO
END DO
end function d2lag_1d_bary
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function lag_2d_bary(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii(2),n
real(r8), intent(in) :: x(3),xnodes(n)
integer(i4) :: i,j,il(3)
real(r8) :: val
val=1.d0
il=(/ii(1),ii(2),n-1-SUM(ii)/)+1
DO i=1,3
  val=val*lag_1d(il(i), x(i), xnodes, il(i))
END DO
end function lag_2d_bary
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function dlag_2d_bary(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii(2),n
real(r8), intent(in) :: x(3),xnodes(n)
integer(i4) :: i,j,k,ix,il(3)
real(r8) :: val(3),norm
val=1.d0
norm=1.d0
il=(/ii(1),ii(2),n-1-SUM(ii)/)+1
DO ix=1,3
  norm=norm*lag_1d_norm(il(ix),xnodes,il(ix))
  DO i=1,3
    IF(i==ix)THEN
      val(ix)=val(ix)*dlag_1d_unnorm(il(i), x(i), xnodes, il(i))
    ELSE
      val(ix)=val(ix)*lag_1d_unnorm(il(i), x(i), xnodes, il(i))
    END IF
  END DO
END DO
val=val/norm
end function dlag_2d_bary
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function d2lag_2d_bary(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii(2),n
real(r8), intent(in) :: x(3),xnodes(n)
integer(i4) :: i,j,k,ix,iy,il(3)
real(r8) :: val(6)
val=1.d0
il=(/ii(1),ii(2),n-1-SUM(ii)/)+1
k=0
DO ix=1,3
  DO iy=ix,3
    k=k+1
    DO i=1,3
      IF(i==ix.OR.i==iy)THEN
        IF(ix==iy)THEN
          val(k)=val(k)*d2lag_1d(il(i), x(i), xnodes, il(i))
        ELSE
          val(k)=val(k)*dlag_1d(il(i), x(i), xnodes, il(i))
        END IF
      ELSE
        val(k)=val(k)*lag_1d(il(i), x(i), xnodes, il(i))
      END IF
    END DO
  END DO
END DO
end function d2lag_2d_bary
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function lag_3d_bary(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii(3),n
real(r8), intent(in) :: x(4),xnodes(n)
integer(i4) :: i,j,l,il(4)
real(r8) :: val
val=1.d0
il=(/ii(1),ii(2),ii(3),n-1-SUM(ii)/)+1
DO i=1,4
  val=val*lag_1d(il(i), x(i), xnodes, il(i))
END DO
end function lag_3d_bary
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function dlag_3d_bary(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii(3),n
real(r8), intent(in) :: x(4),xnodes(n)
integer(i4) :: i,j,l,ix,il(4)
real(r8) :: val(4),norm
val=1.d0
norm=1.d0
il=(/ii(1),ii(2),ii(3),n-1-SUM(ii)/)+1
DO ix=1,4
  norm=norm*lag_1d_norm(il(ix),xnodes,il(ix))
  DO i=1,4
    IF(i==ix)THEN
      val(ix)=val(ix)*dlag_1d_unnorm(il(i), x(i), xnodes, il(i))
    ELSE
      val(ix)=val(ix)*lag_1d_unnorm(il(i), x(i), xnodes, il(i))
    END IF
  END DO
END DO
val=val/norm
end function dlag_3d_bary
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function d2lag_3d_bary(ii,x,xnodes,n) result(val)
integer(i4), intent(in) :: ii(3),n
real(r8), intent(in) :: x(4),xnodes(n)
integer(i4) :: i,j,k,l,ix,iy,il(4)
real(r8) :: val(10)
val=1.d0
il=(/ii(1),ii(2),ii(3),n-1-SUM(ii)/)+1
k=0
DO ix=1,4
  DO iy=ix,4
    k=k+1
    DO i=1,4
      IF(i==ix.OR.i==iy)THEN
        IF(ix==iy)THEN
          val(k)=val(k)*d2lag_1d(il(i), x(i), xnodes, il(i))
        ELSE
          val(k)=val(k)*dlag_1d(il(i), x(i), xnodes, il(i))
        END IF
      ELSE
        val(k)=val(k)*lag_1d(il(i), x(i), xnodes, il(i))
      END IF
    END DO
  END DO
END DO
end function d2lag_3d_bary
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function hpoly_bary(ii,x) result(val)
integer(i4), intent(in) :: ii
real(r8), intent(in) :: x
real(r8) :: val
! Old def -> x * J^(0,0)_n (2*x-1)
! New def -> x * J^(1,1)_n (2*x-1)
select case(ii)
  case(0)
    val=1.d0
  case(1)
    val=x
  case(2)
    ! val=x*(2.d0*x-1.d0)
    val=x*(4.d0*x-2.d0)
  case(3)
    ! val=x*(6.d0*x**2 - 6.d0*x + 1.d0)
    val=x*(15.d0*x**2 - 15.d0*x + 3.d0)
  case(4)
    ! val=x*(20.d0*x**3 - 30.d0*x**2 + 12.d0*x - 1.d0)
    val=x*(56.d0*x**3 - 84.d0*x**2 + 36.d0*x - 4.d0)
  case(5)
    ! val=x*(70.d0*x**4 - 140.d0*x**3 + 90.d0*x**2 - 20.d0*x + 1.d0)
    val=x*(210.d0*x**4 - 420.d0*x**3 + 280.d0*x**2 - 70.d0*x + 5.d0)
  case default
    val=0.d0
end select
end function hpoly_bary
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function dhpoly_bary(ii,x) result(val)
integer(i4), intent(in) :: ii
real(r8), intent(in) :: x
real(r8) :: val
select case(ii)
  case(1)
    val=1.d0
  case(2)
    ! val=4.d0*x-1.d0
    val=8.d0*x-2.d0
  case(3)
    ! val=18.d0*x**2 - 12.d0*x + 1.d0
    val=45.d0*x**2 - 30.d0*x + 3.d0
  case(4)
    ! val=80.d0*x**3 - 90.d0*x**2 + 24.d0*x - 1.d0
    val=224.d0*x**3 - 252.d0*x**2 + 72.d0*x - 4.d0
  case(5)
    ! val=350.d0*x**4 - 560.d0*x**3 + 270.d0*x**2 - 40.d0*x + 1.d0
    val=1050.d0*x**4 - 1680.d0*x**3 + 840.d0*x**2 - 140.d0*x + 5.d0
  case default
    val=0.d0
end select
end function dhpoly_bary
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function d2hpoly_bary(ii,x) result(val)
integer(i4), intent(in) :: ii
real(r8), intent(in) :: x
real(r8) :: val
select case(ii)
  case(2)
    ! val=4.d0
    val=8.d0
  case(3)
    ! val=36.d0*x - 12.d0
    val=90.d0*x - 30.d0
  case(4)
    ! val=240.d0*x**2 - 180.d0*x + 24.d0
    val=672.d0*x**2 - 504.d0*x + 72.d0
  case(5)
    ! val=1400.d0*x**3 - 1680.d0*x**2 + 540.d0*x - 40.d0
    val=4200.d0*x**3 - 5040.d0*x**2 + 1680.d0*x - 140.d0
  case default
    val=0.d0
end select
end function d2hpoly_bary
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function hpoly_2d(inds,x) result(val)
integer(i4), intent(in) :: inds(2)
real(r8), intent(in) :: x(2)
real(r8) :: val
integer(i4) :: i
val=1.d0
DO i=1,2
  val=val*hpoly_bary(inds(i),x(i))
END DO
end function hpoly_2d
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function dhpoly_2d(inds,x) result(val)
integer(i4), intent(in) :: inds(2)
real(r8), intent(in) :: x(2)
real(r8) :: val(3)
val(2)=hpoly_bary(inds(2),x(2))
val(3)=hpoly_bary(inds(1),x(1))
!
val(1)=val(2)*val(3)
val(2)=val(2)*dhpoly_bary(inds(1),x(1))
val(3)=val(3)*dhpoly_bary(inds(2),x(2))
end function dhpoly_2d
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function d2hpoly_2d(inds,k,m,x) result(val)
integer(i4), intent(in) :: inds(2),k,m
real(r8), intent(in) :: x(2)
real(r8) :: val
integer(i4) :: i
IF(k==m)THEN
  i=2-MOD(k,2)
  val=d2hpoly_bary(inds(k),x(k))*hpoly_bary(inds(i),x(i))
ELSE
  val=dhpoly_bary(inds(1),x(1))*dhpoly_bary(inds(2),x(2))
END IF
end function d2hpoly_2d
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
subroutine hpoly_2d_grid(order,inds)
integer(i4), intent(in) :: order
integer(i4), pointer, intent(out) :: inds(:,:)
integer(i4) :: i,j,k,m
logical, ALLOCATABLE, DIMENSION(:,:,:) :: flags
IF(order<1)RETURN
ALLOCATE(inds(2,order*order))
m=0
! WRITE(*,*)'2D'
DO k=1,order
  DO i=1,k
    DO j=1,k
      IF(i<k.AND.j<k)CYCLE
      m=m+1
      inds(:,m)=(/i,j/)
      ! WRITE(*,*)m,i,j
    END DO
  END DO
END DO
! WRITE(*,*)
end subroutine hpoly_2d_grid
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function hpoly_3d(inds,x) result(val)
integer(i4), intent(in) :: inds(3)
real(r8), intent(in) :: x(3)
real(r8) :: val
integer(i4) :: i
val=1.d0
DO i=1,3
  val=val*hpoly_bary(inds(i),x(i))
END DO
end function hpoly_3d
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function dhpoly_3d(inds,x) result(val)
integer(i4), intent(in) :: inds(3)
real(r8), intent(in) :: x(3)
real(r8) :: val(4),vtmp(3)
vtmp(1)=hpoly_bary(inds(1),x(1))
vtmp(2)=hpoly_bary(inds(2),x(2))
vtmp(3)=hpoly_bary(inds(3),x(3))
!
val(1)=vtmp(1)*vtmp(2)*vtmp(3)
val(2)=vtmp(2)*vtmp(3)*dhpoly_bary(inds(1),x(1))
val(3)=vtmp(1)*vtmp(3)*dhpoly_bary(inds(2),x(2))
val(4)=vtmp(1)*vtmp(2)*dhpoly_bary(inds(3),x(3))
end function dhpoly_3d
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
pure function d2hpoly_3d(inds,k,m,x) result(val)
integer(i4), intent(in) :: inds(3),k,m
real(r8), intent(in) :: x(3)
real(r8) :: val
integer(i4) :: i
val=1.d0
! IF(k==m)val=d2hpoly_bary(inds(k),x(k))
! DO i=1,3
!   IF(i==k.NEQV.i==m)THEN
!     val=val*dhpoly_bary(inds(i),x(i))
!   ELSE
!     val=val*hpoly_bary(inds(i),x(i))
!   END IF
! END DO
DO i=1,3
  IF(i==k.OR.i==m)THEN
    IF(k==m)THEN
      val=val*d2hpoly_bary(inds(i),x(i))
    ELSE
      val=val*dhpoly_bary(inds(i),x(i))
    END IF
  ELSE
    val=val*hpoly_bary(inds(i),x(i))
  END IF
END DO
end function d2hpoly_3d
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
subroutine hpoly_3d_grid(order,inds)
integer(i4), intent(in) :: order
integer(i4), pointer, intent(out) :: inds(:,:)
integer(i4) :: i,j,k,l,m
logical, ALLOCATABLE, DIMENSION(:,:,:,:) :: flags
IF(order<1)RETURN
ALLOCATE(inds(3,order**3))
m=0
! WRITE(*,*)'3D'
DO k=1,order
  DO i=1,k
    DO j=1,k
      DO l=1,k
        IF(i<k.AND.j<k.AND.l<k)CYCLE
        m=m+1
        inds(:,m)=(/i,j,l/)
        ! WRITE(*,*)m,i,j,l
      END DO
    END DO
  END DO
END DO
! WRITE(*,*)
end subroutine hpoly_3d_grid
END MODULE oft_lag_poly
