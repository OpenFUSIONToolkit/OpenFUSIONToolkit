c-----------------------------------------------------------------------
c     file spline.f
c     fits functions to cubic splines.
c     Reference: H. Spaeth, "Spline Algorithms for Curves and Surfaces,"
c     Translated from the German by W. D. Hoskins and H. W. Sager.
c     Utilitas Mathematica Publishing Inc., Winnepeg, 1974.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. spline_mod.
c     1. spline_alloc.
c     2. spline_dealloc.
c     3. spline_fit.
c     4. spline_fac.
c     5. spline_eval.
c     6. spline_all_eval.
c     7. spline_write1.
c     8. spline_write2.
c     9. spline_int.
c     10. spline_triluf.
c     11. spline_trilus.
c     12. spline_sherman.
c     13. spline_morrison.
c     14. spline_copy.
c-----------------------------------------------------------------------
c     subprogram 0. spline_mod.
c     defines spline_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE spline_mod
      USE oft_base
      IMPLICIT NONE

      TYPE :: spline_type
      INTEGER :: mx,nqty,ix
      REAL(r8), DIMENSION(:), POINTER :: xs,f,f1,f2,f3
      REAL(r8), DIMENSION(:,:), POINTER :: fs,fs1,fsi,xpower
      REAL(r8), DIMENSION(2) :: x0
      CHARACTER(6), DIMENSION(:), POINTER :: title
      CHARACTER(6) :: name
      LOGICAL :: periodic
      END TYPE spline_type

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. spline_alloc.
c     allocates space for spline_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_alloc(spl,mx,nqty)

      INTEGER, INTENT(IN) :: mx,nqty
      TYPE(spline_type), INTENT(INOUT) :: spl
c-----------------------------------------------------------------------
c     set scalars.
c-----------------------------------------------------------------------
      spl%mx=mx
      spl%nqty=nqty
      spl%ix=0
      spl%periodic=.FALSE.
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      ALLOCATE(spl%xs(0:mx))
      ALLOCATE(spl%f(nqty))
      ALLOCATE(spl%f1(nqty))
      ALLOCATE(spl%f2(nqty))
      ALLOCATE(spl%f3(nqty))
      ALLOCATE(spl%title(0:nqty))
      ALLOCATE(spl%fs(0:mx,nqty))
      ALLOCATE(spl%fs1(0:mx,nqty))
      ALLOCATE(spl%xpower(2,nqty))
      spl%fs1=0
      spl%xpower=0
      spl%x0=0
      NULLIFY(spl%fsi)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_alloc
c-----------------------------------------------------------------------
c     subprogram 2. spline_dealloc.
c     deallocates space for spline_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_dealloc(spl)

      TYPE(spline_type), INTENT(INOUT) :: spl
c-----------------------------------------------------------------------
c     deallocate space.
c-----------------------------------------------------------------------
      DEALLOCATE(spl%xs)
      DEALLOCATE(spl%f)
      DEALLOCATE(spl%f1)
      DEALLOCATE(spl%f2)
      DEALLOCATE(spl%f3)
      DEALLOCATE(spl%title)
      DEALLOCATE(spl%fs)
      DEALLOCATE(spl%fs1)
      DEALLOCATE(spl%xpower)
      IF(ASSOCIATED(spl%fsi))DEALLOCATE(spl%fsi)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. spline_fit.
c     fits real functions to cubic splines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_fit(spl,endmode)

      TYPE(spline_type), INTENT(INOUT) :: spl
      CHARACTER(*), INTENT(IN) :: endmode

      INTEGER :: iqty,iside
      REAL(r8), DIMENSION(-1:1,0:spl%mx) :: a
      REAL(r8), DIMENSION(spl%mx) :: b
      REAL(r8), DIMENSION(4) :: cl,cr
      REAL(r8), DIMENSION(0:spl%mx) :: xfac
c-----------------------------------------------------------------------
c     extract powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         DO iqty=1,spl%nqty
            IF(spl%xpower(iside,iqty) /= 0)THEN
               xfac=1/ABS(spl%xs-spl%x0(iside))**spl%xpower(iside,iqty)
               spl%fs(:,iqty)=spl%fs(:,iqty)*xfac
            ENDIF
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     set up grid matrix.
c-----------------------------------------------------------------------
      CALL spline_fac(spl,a,b,cl,cr,endmode)
c-----------------------------------------------------------------------
c     compute first derivatives, interior.
c-----------------------------------------------------------------------
      DO iqty=1,spl%nqty
         spl%fs1(1:spl%mx-1,iqty)=
     $        3*((spl%fs(2:spl%mx,iqty)-spl%fs(1:spl%mx-1,iqty))
     $        *b(2:spl%mx)
     $        +(spl%fs(1:spl%mx-1,iqty)-spl%fs(0:spl%mx-2,iqty))
     $        *b(1:spl%mx-1))
      ENDDO
c-----------------------------------------------------------------------
c     extrapolation boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(endmode)
      CASE("extrap")
         DO iqty=1,spl%nqty
            spl%fs1(0,iqty)=SUM(cl(1:4)*spl%fs(0:3,iqty))
            spl%fs1(spl%mx,iqty)=SUM(cr(1:4)
     $           *spl%fs(spl%mx:spl%mx-3:-1,iqty))
            spl%fs1(1,iqty)=spl%fs1(1,iqty)-spl%fs1(0,iqty)
     $           /(spl%xs(1)-spl%xs(0))
            spl%fs1(spl%mx-1,iqty)=
     $           spl%fs1(spl%mx-1,iqty)-spl%fs1(spl%mx,iqty)
     $           /(spl%xs(spl%mx)-spl%xs(spl%mx-1))
         ENDDO
         CALL spline_trilus(a(:,1:spl%mx-1),spl%fs1(1:spl%mx-1,:))
c-----------------------------------------------------------------------
c     not-a-knot boundary conditions.
c-----------------------------------------------------------------------
      CASE("not-a-knot")
         spl%fs1(1,:)=spl%fs1(1,:)-(2*spl%fs(1,:)
     $        -spl%fs(0,:)-spl%fs(2,:))*2*b(1)
         spl%fs1(spl%mx-1,:)=spl%fs1(spl%mx-1,:)
     $        +(2*spl%fs(spl%mx-1,:)-spl%fs(spl%mx,:)
     $        -spl%fs(spl%mx-2,:))*2*b(spl%mx)
         CALL spline_trilus(a(:,1:spl%mx-1),spl%fs1(1:spl%mx-1,:))
         spl%fs1(0,:)=(2*(2*spl%fs(1,:)-spl%fs(0,:)-spl%fs(2,:))
     $        +(spl%fs1(1,:)+spl%fs1(2,:))*(spl%xs(2)-spl%xs(1))
     $        -spl%fs1(1,:)*(spl%xs(1)-spl%xs(0)))/(spl%xs(1)-spl%xs(0))
         spl%fs1(spl%mx,:)=
     $        (2*(spl%fs(spl%mx-2,:)+spl%fs(spl%mx,:)
     $        -2*spl%fs(spl%mx-1,:))
     $        +(spl%fs1(spl%mx-1,:)+spl%fs1(spl%mx-2,:))
     $        *(spl%xs(spl%mx-1)-spl%xs(spl%mx-2))
     $        -spl%fs1(spl%mx-1,:)
     $        *(spl%xs(spl%mx)-spl%xs(spl%mx-1)))
     $        /(spl%xs(spl%mx)-spl%xs(spl%mx-1))
c-----------------------------------------------------------------------
c     periodic boudary conditions.
c-----------------------------------------------------------------------
      CASE("periodic")
         spl%periodic=.TRUE.
         spl%fs1(0,:)=3*((spl%fs(1,:)-spl%fs(0,:))*b(1)
     $        +(spl%fs(0,:)-spl%fs(spl%mx-1,:))*b(spl%mx))
         CALL spline_morrison(a(:,0:spl%mx-1),spl%fs1(0:spl%mx-1,:))
         spl%fs1(spl%mx,:)=spl%fs1(0,:)
c-----------------------------------------------------------------------
c     unrecognized boundary condition.
c-----------------------------------------------------------------------
      CASE("none")
        !spl%fs1(0,:)=0.d0; spl%fs1(spl%mx,:)=0.d0
        CALL spline_trilus(a(:,1:spl%mx-1),spl%fs1(1:spl%mx-1,:))
        spl%fs1(0,:)=(6*(spl%fs(1,:)-spl%fs(0,:))/(spl%xs(1)-spl%xs(0))
     $        -2*spl%fs1(1,:))/4
        spl%fs1(spl%mx,:)=(6*(spl%fs(spl%mx,:)-spl%fs(spl%mx-1,:))/
     $  (spl%xs(1)-spl%xs(0))-2*spl%fs1(spl%mx-1,:))/4
c        spl%fs1(0,:)=0.d0; spl%fs1(spl%mx,:)=0.d0
      CASE DEFAULT
c         CALL program_stop("Cannot recognize endmode = "//TRIM(endmode))
         CALL oft_abort("Cannot recognize endmode = "//TRIM(endmode),
     $   "spline_fit","spline.f")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_fit
c-----------------------------------------------------------------------
c     subprogram 4. spline_fac.
c     sets up matrix for cubic spline fitting.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_fac(spl,a,b,cl,cr,endmode)

      TYPE(spline_type), INTENT(IN) :: spl
      REAL(r8), DIMENSION(-1:1,0:spl%mx), INTENT(OUT) :: a
      REAL(r8), DIMENSION(spl%mx), INTENT(OUT) :: b
      REAL(r8), DIMENSION(4), INTENT(OUT) :: cl,cr
      CHARACTER(*), INTENT(IN) :: endmode

      INTEGER :: j
c-----------------------------------------------------------------------
c     compute interior matrix.
c-----------------------------------------------------------------------
      b=1/(spl%xs(1:spl%mx)-spl%xs(0:spl%mx-1))
      DO j=1,spl%mx-1
         a(-1,j)=b(j)
         a(0,j)=2*(b(j)+b(j+1))
         a(1,j)=b(j+1)
      ENDDO
c-----------------------------------------------------------------------
c     extrapolation boundary conditions.
c-----------------------------------------------------------------------
      SELECT CASE(endmode)
      CASE("extrap")
         b=b*b
         cl(1)=(spl%xs(0)*(3*spl%xs(0)
     $        -2*(spl%xs(1)+spl%xs(2)+spl%xs(3)))
     $        +spl%xs(1)*spl%xs(2)+spl%xs(1)*spl%xs(3)
     $        +spl%xs(2)*spl%xs(3))
     $        /((spl%xs(0)-spl%xs(1))*(spl%xs(0)-spl%xs(2))
     $        *(spl%xs(0)-spl%xs(3)))
         cl(2)=((spl%xs(2)-spl%xs(0))*(spl%xs(3)-spl%xs(0)))
     $        /((spl%xs(1)-spl%xs(0))*(spl%xs(1)-spl%xs(2))
     $        *(spl%xs(1)-spl%xs(3)))
         cl(3)=((spl%xs(0)-spl%xs(1))*(spl%xs(3)-spl%xs(0)))
     $        /((spl%xs(0)-spl%xs(2))*(spl%xs(1)-spl%xs(2))
     $        *(spl%xs(3)-spl%xs(2)))
         cl(4)=((spl%xs(1)-spl%xs(0))*(spl%xs(2)-spl%xs(0)))
     $        /((spl%xs(3)-spl%xs(0))*(spl%xs(3)-spl%xs(1))
     $        *(spl%xs(3)-spl%xs(2)))
         cr(1)=(spl%xs(spl%mx)*(3*spl%xs(spl%mx)
     $        -2*(spl%xs(spl%mx-1)+spl%xs(spl%mx-2)
     $        +spl%xs(spl%mx-3)))+spl%xs(spl%mx-1)
     $        *spl%xs(spl%mx-2)+spl%xs(spl%mx-1)*spl%xs(spl%mx
     $        -3)+spl%xs(spl%mx-2)*spl%xs(spl%mx-3))
     $        /((spl%xs(spl%mx)-spl%xs(spl%mx-1))
     $        *(spl%xs(spl%mx)-spl%xs(spl%mx-2))*(spl%xs(spl%mx
     $        )-spl%xs(spl%mx-3)))
         cr(2)=((spl%xs(spl%mx-2)-spl%xs(spl%mx))
     $        *(spl%xs(spl%mx-3)-spl%xs(spl%mx)))
     $        /((spl%xs(spl%mx-1)-spl%xs(spl%mx))
     $        *(spl%xs(spl%mx-1)-spl%xs(spl%mx-2))
     $        *(spl%xs(spl%mx-1)-spl%xs(spl%mx-3)))
         cr(3)=((spl%xs(spl%mx)-spl%xs(spl%mx-1))
     $        *(spl%xs(spl%mx-3)-spl%xs(spl%mx)))
     $        /((spl%xs(spl%mx)-spl%xs(spl%mx-2))
     $        *(spl%xs(spl%mx-1)-spl%xs(spl%mx-2))
     $        *(spl%xs(spl%mx-3)-spl%xs(spl%mx-2)))
         cr(4)=((spl%xs(spl%mx-1)-spl%xs(spl%mx))
     $        *(spl%xs(spl%mx-2)-spl%xs(spl%mx)))
     $        /((spl%xs(spl%mx-3)-spl%xs(spl%mx))
     $        *(spl%xs(spl%mx-3)-spl%xs(spl%mx-1))
     $        *(spl%xs(spl%mx-3)-spl%xs(spl%mx-2)))
         CALL spline_triluf(a(:,1:spl%mx-1))
c-----------------------------------------------------------------------
c     not-a-knot boundary conditions.
c-----------------------------------------------------------------------
      CASE("not-a-knot")
         b=b*b
         a(0,1)=a(0,1)+(spl%xs(2)+spl%xs(0)-2*spl%xs(1))*b(1)
         a(1,1)=a(1,1)+(spl%xs(2)-spl%xs(1))*b(1)
         a(0,spl%mx-1)=a(0,spl%mx-1)
     $        +(2*spl%xs(spl%mx-1)-spl%xs(spl%mx-2)
     $        -spl%xs(spl%mx))*b(spl%mx)
         a(-1,spl%mx-1)=a(-1,spl%mx-1)
     $        +(spl%xs(spl%mx-1)-spl%xs(spl%mx-2))*b(spl%mx)
         CALL spline_triluf(a(:,1:spl%mx-1))
c-----------------------------------------------------------------------
c     periodic boundary conditions.
c-----------------------------------------------------------------------
      CASE("periodic")
         a(0,0:spl%mx:spl%mx)=2*(b(spl%mx)+b(1))
         a(1,0)=b(1)
         a(-1,0)=b(spl%mx)
         b=b*b
         CALL spline_sherman(a(:,0:spl%mx-1))
c-----------------------------------------------------------------------
c     unrecognized boundary condition.
c-----------------------------------------------------------------------
      CASE("none")
        !a((/-1,0,1/),0)=(/0.d0,1.d0,0.d0/)
        !a((/-1,0,1/),spl%mx)=(/0.d0,1.d0,0.d0/)
        a(1,1)=0.d0; a(-1,spl%mx-1)=0.d0
        b=b*b
        CALL spline_triluf(a(:,1:spl%mx-1))
      CASE DEFAULT
c         CALL program_stop("Cannot recognize endmode = "//TRIM(endmode))
         CALL oft_abort("Cannot recognize endmode = "//TRIM(endmode),
     $   "spline_fac","spline.f")
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_fac
c-----------------------------------------------------------------------
c     subprogram 5. spline_eval.
c     evaluates real cubic spline function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_eval(spl,x,mode)

      TYPE(spline_type), INTENT(INOUT) :: spl
      REAL(r8), INTENT(IN) :: x
      INTEGER, INTENT(IN) :: mode

      INTEGER :: iqty,iside
      REAL(r8) :: xx,d,z,z1,xfac,dx
      REAL(r8) :: g,g1,g2,g3
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      xx=x
      spl%ix=MAX(spl%ix,0)
      spl%ix=MIN(spl%ix,spl%mx-1)
c-----------------------------------------------------------------------
c     normalize interval for periodic splines.
c-----------------------------------------------------------------------
      IF(spl%periodic)THEN
         DO
            IF(xx < spl%xs(spl%mx))EXIT
            xx=xx-spl%xs(spl%mx)
         ENDDO
         DO
            IF(xx >= spl%xs(0))EXIT
            xx=xx+spl%xs(spl%mx)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find cubic spline interval.
c-----------------------------------------------------------------------
      DO
         IF(spl%ix <= 0)EXIT
         IF(xx >= spl%xs(spl%ix))EXIT
         spl%ix=spl%ix-1
      ENDDO
      DO
         IF(spl%ix >= spl%mx-1)EXIT
         IF(xx < spl%xs(spl%ix+1))EXIT
         spl%ix=spl%ix+1
      ENDDO
c-----------------------------------------------------------------------
c     evaluate offset and related quantities.
c-----------------------------------------------------------------------
      d=spl%xs(spl%ix+1)-spl%xs(spl%ix)
      z=(xx-spl%xs(spl%ix))/d
      z1=1-z
c-----------------------------------------------------------------------
c     evaluate functions.
c-----------------------------------------------------------------------
      spl%f=spl%fs(spl%ix,1:spl%nqty)*z1*z1*(3-2*z1)
     $     +spl%fs(spl%ix+1,1:spl%nqty)*z*z*(3-2*z)
     $     +d*z*z1*(spl%fs1(spl%ix,1:spl%nqty)*z1
     $     -spl%fs1(spl%ix+1,1:spl%nqty)*z)
c-----------------------------------------------------------------------
c     evaluate first derivatives.
c-----------------------------------------------------------------------
      IF(mode > 0)THEN
         spl%f1=6*(spl%fs(spl%ix+1,1:spl%nqty)
     $        -spl%fs(spl%ix,1:spl%nqty))*z*z1/d
     $        +spl%fs1(spl%ix,1:spl%nqty)*z1*(3*z1-2)
     $        +spl%fs1(spl%ix+1,1:spl%nqty)*z*(3*z-2)
      ENDIF
c-----------------------------------------------------------------------
c     evaluate second derivatives.
c-----------------------------------------------------------------------
      IF(mode > 1)THEN
         spl%f2=(6*(spl%fs(spl%ix+1,1:spl%nqty)
     $        -spl%fs(spl%ix,1:spl%nqty))*(z1-z)/d
     $        -spl%fs1(spl%ix,1:spl%nqty)*(6*z1-2)
     $        +spl%fs1(spl%ix+1,1:spl%nqty)*(6*z-2))/d
      ENDIF
c-----------------------------------------------------------------------
c     evaluate third derivatives.
c-----------------------------------------------------------------------
      IF(mode > 2)THEN
         spl%f3=(12*(spl%fs(spl%ix,1:spl%nqty)
     $        -spl%fs(spl%ix+1,1:spl%nqty))/d
     $        +6*(spl%fs1(spl%ix,1:spl%nqty)
     $        +spl%fs1(spl%ix+1,1:spl%nqty)))/(d*d)
      ENDIF
c-----------------------------------------------------------------------
c     restore powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         dx=x-spl%x0(iside)
         DO iqty=1,spl%nqty
            IF(spl%xpower(iside,iqty) == 0)CYCLE
            xfac=ABS(dx)**spl%xpower(iside,iqty)
            g=spl%f(iqty)*xfac
            IF(mode > 0)g1=(spl%f1(iqty)+spl%f(iqty)
     $           *spl%xpower(iside,iqty)/dx)*xfac
            IF(mode > 1)g2=(spl%f2(iqty)+spl%xpower(iside,iqty)/dx
     $           *(2*spl%f1(iqty)+(spl%xpower(iside,iqty)-1)
     $           *spl%f(iqty)/dx))*xfac
            IF(mode > 2)g3=(spl%f3(iqty)+spl%xpower(iside,iqty)/dx
     $           *(3*spl%f2(iqty)+(spl%xpower(iside,iqty)-1)/dx
     $           *(3*spl%f1(iqty)+(spl%xpower(iside,iqty)-2)/dx
     $           *spl%f(iqty))))*xfac
            spl%f(iqty)=g
            IF(mode > 0)spl%f1(iqty)=g1
            IF(mode > 1)spl%f2(iqty)=g2
            IF(mode > 2)spl%f3(iqty)=g3
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_eval
c-----------------------------------------------------------------------
c     subprogram 6. spline_all_eval.
c     evaluates cubic spline function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_all_eval(spl,z,f,f1,f2,f3,mode)

      TYPE(spline_type), INTENT(INOUT) :: spl
      REAL(r8), INTENT(IN) :: z
      REAL(r8), DIMENSION(spl%mx,spl%nqty), INTENT(OUT) :: f,f1,f2,f3
      INTEGER, INTENT(IN) :: mode

      INTEGER :: iqty,nqty,n,iside
      REAL(r8) :: z1
      REAL(r8), DIMENSION(spl%mx) :: d,xfac,dx
      REAL(r8), DIMENSION(spl%mx) :: g,g1,g2,g3
c-----------------------------------------------------------------------
c     evaluate offset and related quantities.
c-----------------------------------------------------------------------
      n=spl%mx
      nqty=spl%nqty
      z1=1-z
      d=spl%xs(1:n)-spl%xs(0:n-1)
c-----------------------------------------------------------------------
c     evaluate functions.
c-----------------------------------------------------------------------
      DO iqty=1,nqty
         f(:,iqty)=spl%fs(0:n-1,iqty)*z1*z1*(3-2*z1)
     $        +spl%fs(1:n,iqty)*z*z*(3-2*z)
     $        +d*z*z1*(spl%fs1(0:n-1,iqty)*z1-spl%fs1(1:n,iqty)*z)
      ENDDO
c-----------------------------------------------------------------------
c     evaluate first derivatives.
c-----------------------------------------------------------------------
      IF(mode > 0)THEN
         DO iqty=1,nqty
            f1(:,iqty)=6*(spl%fs(1:n,iqty)-spl%fs(0:n-1,iqty))*z*z1/d
     $           +spl%fs1(0:n-1,iqty)*z1*(3*z1-2)
     $           +spl%fs1(1:n,iqty)*z*(3*z-2)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     evaluate second derivatives.
c-----------------------------------------------------------------------
      IF(mode > 1)THEN
         DO iqty=1,nqty
            f2(:,iqty)=(6*(spl%fs(1:n,iqty)-spl%fs(0:n-1,iqty))*(z1-z)/d
     $           -spl%fs1(0:n-1,iqty)*(6*z1-2)
     $           +spl%fs1(1:n,iqty)*(6*z-2))/d
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     evaluate third derivatives.
c-----------------------------------------------------------------------
      IF(mode > 2)THEN
         DO iqty=1,nqty
            f3(:,iqty)=(12*(spl%fs(0:n-1,iqty)-spl%fs(1:n,iqty))/d
     $           +6*(spl%fs1(0:n-1,iqty)+spl%fs1(1:n,iqty)))/(d*d)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     restore powers.
c-----------------------------------------------------------------------
      DO iside=1,2
         dx=(spl%xs(0:spl%mx-1)+z*d(1:spl%mx))-spl%x0(iside)
         DO iqty=1,spl%nqty
            IF(spl%xpower(iside,iqty) == 0)CYCLE
            xfac=ABS(dx)**spl%xpower(iside,iqty)
            g=f(:,iqty)*xfac
            IF(mode > 0)g1=(f1(:,iqty)
     $           +f(:,iqty)*spl%xpower(iside,iqty)/dx)*xfac
            IF(mode > 1)g2=(f2(:,iqty)+spl%xpower(iside,iqty)/dx
     $           *(2*f1(:,iqty)+(spl%xpower(iside,iqty)-1)
     $           *f(:,iqty)/dx))*xfac
     $
            IF(mode > 2)g3=(f3(:,iqty)+spl%xpower(iside,iqty)/dx
     $           *(3*f2(:,iqty)+(spl%xpower(iside,iqty)-1)/dx
     $           *(3*f1(:,iqty)+(spl%xpower(iside,iqty)-2)/dx
     $           *f(:,iqty))))*xfac
            f(:,iqty)=g
            IF(mode > 0)f1(:,iqty)=g1
            IF(mode > 1)f2(:,iqty)=g2
            IF(mode > 2)f3(:,iqty)=g3
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_all_eval
c-----------------------------------------------------------------------
c     subprogram 7. spline_write1.
c     produces ascii and binary output for real cubic spline fits.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_write1(spl,out,bin,iua,iub,interp)

      TYPE(spline_type), INTENT(INOUT) :: spl
      LOGICAL, INTENT(IN) :: out,bin
      INTEGER, INTENT(IN) :: iua,iub
      LOGICAL, INTENT(IN) :: interp

      CHARACTER(30) :: format1,format2
      INTEGER :: i,j
      REAL(r8) :: x,dx
c-----------------------------------------------------------------------
c     formats.
c-----------------------------------------------------------------------
 10   FORMAT('(/3x,"ix",',i2.2,'(4x,a6,1x)/)')
 20   FORMAT('(i5,1p,',i2.2,'e11.3)')
c 30   FORMAT('(/3x,"ix",2x,"j",',i2.2,'(4x,a6,1x)/)')
c-----------------------------------------------------------------------
c     abort.
c-----------------------------------------------------------------------
      IF(.NOT.out.AND..NOT.bin)RETURN
c-----------------------------------------------------------------------
c     print node values.
c-----------------------------------------------------------------------
      IF(out)THEN
         WRITE(format1,10)spl%nqty+1
         WRITE(format2,20)spl%nqty+1
         WRITE(iua,'(/1x,a)')'node values:'
         WRITE(iua,format1)spl%title(0:spl%nqty)
      ENDIF
      DO i=0,spl%mx
         CALL spline_eval(spl,spl%xs(i),0)
         IF(out)WRITE(iua,format2)i,spl%xs(i),spl%f
         IF(bin)WRITE(iub)REAL(spl%xs(i),4),REAL(spl%f,4)
      ENDDO
      IF(out)WRITE(iua,format1)spl%title(0:spl%nqty)
      IF(bin)WRITE(iub)
      IF(.NOT. interp)RETURN
c-----------------------------------------------------------------------
c     print header for interpolated values.
c-----------------------------------------------------------------------
      IF(out)THEN
         WRITE(iua,'(/1x,a)')'interpolated values:'
         WRITE(iua,format1)spl%title(0:spl%nqty)
      ENDIF
c-----------------------------------------------------------------------
c     print interpolated values.
c-----------------------------------------------------------------------
      DO i=0,spl%mx-1
         dx=(spl%xs(i+1)-spl%xs(i))/4
         DO j=0,4
            x=spl%xs(i)+j*dx
            CALL spline_eval(spl,x,0)
            IF(out)WRITE(iua,format2)i,x,spl%f
            IF(bin)WRITE(iub)REAL(x,4),REAL(spl%f,4)
         ENDDO
         IF(out)WRITE(iua,'(1x)')
      ENDDO
c-----------------------------------------------------------------------
c     print final interpolated values.
c-----------------------------------------------------------------------
      x=spl%xs(spl%mx)
      CALL spline_eval(spl,x,0)
      IF(out)THEN
         WRITE(iua,format2)i,x,spl%f
         WRITE(iua,format1)spl%title
      ENDIF
      IF(bin)THEN
         WRITE(iub)REAL(x,4),REAL(spl%f,4)
         WRITE(iub)
         CLOSE(iub)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_write1
c-----------------------------------------------------------------------
c     subprogram 8. spline_write2.
c     produces ascii and binary output for real cubic spline fits.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_write2(spl,out,bin,iua,iub,interp)

      TYPE(spline_type), INTENT(INOUT) :: spl
      LOGICAL, INTENT(IN) :: out,bin
      INTEGER, INTENT(IN) :: iua,iub
      LOGICAL, INTENT(IN) :: interp

      CHARACTER(30) :: format1,format2
      INTEGER :: i,j,iz
      REAL(r8) :: x,dx,z
      REAL(r8), DIMENSION(spl%mx,spl%nqty) :: f
      REAL(r8), DIMENSION(0:4*spl%mx,spl%nqty) :: g
c-----------------------------------------------------------------------
c     formats.
c-----------------------------------------------------------------------
 10   FORMAT('(/4x,"i",',i2.2,'(4x,a6,1x)/)')
 20   FORMAT('(i5,1p,',i2.2,'e11.3)')
c 30   FORMAT('(/4x,"i",2x,"j",',i2.2,'(4x,a6,1x)/)')
c-----------------------------------------------------------------------
c     compute values.
c-----------------------------------------------------------------------
      z=0
      DO iz=0,3
         CALL spline_all_eval(spl,z,f,f,f,f,0)
         g(iz:4*spl%mx-1:4,:)=f
         z=z+.25_r8
      ENDDO
      CALL spline_eval(spl,spl%xs(spl%mx),0)
      g(4*spl%mx,:)=spl%f
c-----------------------------------------------------------------------
c     print node values.
c-----------------------------------------------------------------------
      IF(.NOT.out.AND..NOT.bin)RETURN
      IF(out)THEN
         WRITE(format1,10)spl%nqty+1
         WRITE(format2,20)spl%nqty+1
         WRITE(iua,'(/1x,a)')'node values:'
         WRITE(iua,format1)spl%title(0:spl%nqty)
      ENDIF
      DO i=0,spl%mx
         IF(out)WRITE(iua,format2)i,spl%xs(i),g(4*i,:)
         IF(bin)WRITE(iub)REAL(spl%xs(i),4),REAL(g(4*i,:),4)
      ENDDO
      IF(out)WRITE(iua,format1)spl%title(0:spl%nqty)
      IF(bin)WRITE(iub)
c-----------------------------------------------------------------------
c     print header for interpolated values.
c-----------------------------------------------------------------------
      IF(.NOT. interp)RETURN
      IF(out)THEN
         WRITE(iua,'(/1x,a)')'interpolated values:'
         WRITE(iua,format1)spl%title(0:spl%nqty)
      ENDIF
c-----------------------------------------------------------------------
c     print interpolated values.
c-----------------------------------------------------------------------
      DO i=0,spl%mx-1
         dx=(spl%xs(i+1)-spl%xs(i))/4
         DO j=0,4
            x=spl%xs(i)+j*dx
            IF(out)WRITE(iua,format2)i,x,g(4*i+j,:)
            IF(bin)WRITE(iub)REAL(x,4),REAL(g(4*i+j,:),4)
         ENDDO
         IF(out)WRITE(iua,'(1x)')
      ENDDO
c-----------------------------------------------------------------------
c     print final interpolated values.
c-----------------------------------------------------------------------
      IF(out)THEN
         WRITE(iua,format2)i,x,g(spl%mx*4,:)
         WRITE(iua,format1)spl%title(0:spl%nqty)
      ENDIF
      IF(bin)THEN
         WRITE(iub)REAL(x,4),REAL(g(spl%mx*4,:),4)
         WRITE(iub)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_write2
c-----------------------------------------------------------------------
c     subprogram 9. spline_int.
c     integrates real cubic splines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_int(spl)

      TYPE(spline_type), INTENT(INOUT) :: spl

      INTEGER :: ix,iqty,ig
      REAL(r8), DIMENSION(spl%mx) :: dx
      REAL(r8), DIMENSION(spl%mx,spl%nqty) :: term,f

      INTEGER, PARAMETER :: mg=4
      REAL(r8), DIMENSION(mg) :: xg=(1+(/-0.861136311594053_r8,
     $     -0.339981043584856_r8,0.339981043584856_r8,
     $     0.861136311594053_r8/))/2
      REAL(r8), DIMENSION(mg) :: wg=(/0.347854845137454_r8,
     $     0.652145154862546_r8,0.652145154862546_r8,
     $     0.347854845137454_r8/)/2
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      IF(.NOT.ASSOCIATED(spl%fsi))ALLOCATE(spl%fsi(0:spl%mx,spl%nqty))
      dx=spl%xs(1:spl%mx)-spl%xs(0:spl%mx-1)
      term=0
c-----------------------------------------------------------------------
c     compute integrals over intervals.
c-----------------------------------------------------------------------
      DO iqty=1,spl%nqty
         IF(spl%xpower(1,iqty) == 0 .AND. spl%xpower(2,iqty) == 0)THEN
            term(:,iqty)=dx/12
     $           *(6*(spl%fs(0:spl%mx-1,iqty)+spl%fs(1:spl%mx,iqty))
     $           +dx*(spl%fs1(0:spl%mx-1,iqty)-spl%fs1(1:spl%mx,iqty)))
         ELSE
            DO ig=1,mg
               CALL spline_all_eval(spl,xg(ig),f,f,f,f,0)
               term(:,iqty)=term(:,iqty)+dx*wg(ig)*f(:,iqty)
            ENDDO
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     accumulate over intervals.
c-----------------------------------------------------------------------
      spl%fsi(0,:)=0
      DO ix=1,spl%mx
         spl%fsi(ix,:)=spl%fsi(ix-1,:)+term(ix,:)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_int
c-----------------------------------------------------------------------
c     subprogram 10. spline_triluf.
c     performs tridiagonal LU factorization.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_triluf(a)

      REAL(r8), DIMENSION(-1:,:), INTENT(INOUT) :: a

      INTEGER :: i,j,k,jmin,jmax,n
c-----------------------------------------------------------------------
c     begin loop over rows and define limits.
c-----------------------------------------------------------------------
      n=SIZE(a,2)
      DO i=1,n
         jmin=MAX(1-i,-1)
         jmax=MIN(n-i,1)
c-----------------------------------------------------------------------
c     compute lower elements.
c-----------------------------------------------------------------------
         DO j=jmin,-1
            DO k=MAX(jmin,j-1),j-1
               a(j,i)=a(j,i)-a(k,i)*a(j-k,i+k)
            ENDDO
            a(j,i)=a(j,i)*a(0,i+j)
         ENDDO
c-----------------------------------------------------------------------
c     compute diagonal element
c-----------------------------------------------------------------------
         DO k=MAX(jmin,-1),-1
            a(0,i)=a(0,i)-a(k,i)*a(-k,i+k)
         ENDDO
         a(0,i)=1/a(0,i)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_triluf
c-----------------------------------------------------------------------
c     subprogram 11. spline_trilus.
c     performs tridiagonal LU solution.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_trilus(a,x)

      REAL(r8), DIMENSION(-1:,:), INTENT(IN) :: a
      REAL(r8), DIMENSION(:,:), INTENT(INOUT) :: x

      INTEGER :: i,j,n
c-----------------------------------------------------------------------
c     down sweep.
c-----------------------------------------------------------------------
      n=SIZE(a,2)
      DO i=1,n
         DO j=MAX(1-i,-1),-1
            x(i,:)=x(i,:)-a(j,i)*x(i+j,:)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     up sweep.
c-----------------------------------------------------------------------
      DO i=n,1,-1
         DO j=1,MIN(n-i,1)
            x(i,:)=x(i,:)-a(j,i)*x(i+j,:)
         ENDDO
         x(i,:)=x(i,:)*a(0,i)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_trilus
c-----------------------------------------------------------------------
c     subprogram 12. spline_sherman.
c     uses Sherman-Morrison formula to factor periodic matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_sherman(a)

      REAL(r8), DIMENSION(-1:,:), INTENT(INOUT) :: a

      INTEGER :: j,n
      REAL(r8), DIMENSION(SIZE(a,2),1) :: u
c-----------------------------------------------------------------------
c     prepare matrices.
c-----------------------------------------------------------------------
      n=SIZE(a,2)
      a(0,1)=a(0,1)-a(-1,1)
      a(0,n)=a(0,n)-a(-1,1)
      u=RESHAPE((/1._r8,(0._r8,j=2,n-1),1._r8/),SHAPE(u))
      CALL spline_triluf(a)
      CALL spline_trilus(a,u)
      a(-1,1)=a(-1,1)/(1+a(-1,1)*(u(1,1)+u(n,1)))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_sherman
c-----------------------------------------------------------------------
c     subprogram 13. spline_morrison.
c     uses Sherman-Morrison formula to solve periodic matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_morrison(a,x)

      REAL(r8), DIMENSION(-1:,:), INTENT(IN) :: a
      REAL(r8), DIMENSION(:,:), INTENT(INOUT) :: x

      INTEGER :: n
      REAL(r8), DIMENSION(SIZE(x,1),SIZE(x,2)) :: y
c-----------------------------------------------------------------------
c     solve for x.
c-----------------------------------------------------------------------
      n=SIZE(a,2)
      y=x
      CALL spline_trilus(a,y)
      x(1,:)=x(1,:)-a(-1,1)*(y(1,:)+y(n,:))
      x(n,:)=x(n,:)-a(-1,1)*(y(1,:)+y(n,:))
      CALL spline_trilus(a,x)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_morrison
c-----------------------------------------------------------------------
c     subprogram 14. spline_copy.
c     copies one spline_type to another.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE spline_copy(spl1,spl2)

      TYPE(spline_type), INTENT(IN) :: spl1
      TYPE(spline_type), INTENT(INOUT) :: spl2
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(ASSOCIATED(spl2%xs))CALL spline_dealloc(spl2)
      CALL spline_alloc(spl2,spl1%mx,spl1%nqty)
      spl2%xs=spl1%xs
      spl2%fs=spl1%fs
      spl2%fs1=spl1%fs1
      spl2%name=spl1%name
      spl2%title=spl1%title
      spl2%periodic=spl1%periodic
      spl2%xpower=spl1%xpower
      spl2%x0=spl1%x0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE spline_copy
      END MODULE spline_mod
