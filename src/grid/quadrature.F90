!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_quadrature.F90
!
!> Definition of generic quadrature type
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_grid
!------------------------------------------------------------------------------
module oft_quadrature
USE oft_local
IMPLICIT NONE
PRIVATE
!------------------------------------------------------------------------------
!> Quadrature point definition structure
!------------------------------------------------------------------------------
TYPE, PUBLIC :: oft_quad_type
  INTEGER(i4) :: dim = -1 !< Spatial dimension
  INTEGER(i4) :: order = -1 !< Order of accuracy
  INTEGER(i4) :: np = 0 !< Number of pts
  REAL(r8), POINTER, DIMENSION(:,:) :: pts => NULL() !< Evaluation points (dim+1,np)
  REAL(r8), POINTER, DIMENSION(:) :: wts => NULL() !< Evaluation weights (np)
CONTAINS
  !> Delete quadrature object
  PROCEDURE :: delete => delete_quad
END TYPE oft_quad_type
CONTAINS
!------------------------------------------------------------------------------
!> Delete quadrature object
!------------------------------------------------------------------------------
subroutine delete_quad(self)
class(oft_quad_type), intent(inout) :: self
self%dim = -1
self%order = -1
self%np = 0
if(associated(self%pts))deallocate(self%pts)
if(associated(self%wts))deallocate(self%wts)
end subroutine delete_quad    
end module oft_quadrature
!------------------------------------------------------------------------------
!> Definition of 1D, 2D, and 3D Gauss-Legendre quadrature rules
!!
!! Zeroes and weights taken from online tabluation at
!! "https://pomax.github.io/bezierinfo/legendre-gauss.html" and verified for
!! accuracy
!!
!! @authors Chris Hansen
!! @date October 2018
!! @ingroup doxy_oft_grid
!------------------------------------------------------------------------------
module oft_gauss_quadrature
USE oft_base
USE oft_quadrature
implicit none
#include "local.h"
PRIVATE
PUBLIC set_quad_1d, set_quad_2d, set_quad_3d, oft_quad_type
!
REAL(r8), PARAMETER :: gauss_a2(2) = (/-0.5773502691896257d0, 0.5773502691896257d0/)
REAL(r8), PARAMETER :: gauss_w2(2) = (/1.0000000000000000d0, 1.0000000000000000d0/)
!
REAL(r8), PARAMETER :: gauss_a4(4) = (/-0.8611363115940526d0, -0.3399810435848563d0, &
    0.3399810435848563d0, 0.8611363115940526d0/)
REAL(r8), PARAMETER :: gauss_w4(4) = (/0.3478548451374538d0, 0.6521451548625461d0, &
    0.6521451548625461d0,0.3478548451374538d0/)
!
REAL(r8), PARAMETER :: gauss_a6(6) = (/-0.9324695142031521d0, -0.6612093864662645d0, &
    -0.2386191860831969d0, 0.2386191860831969d0, 0.6612093864662645d0, 0.9324695142031521d0/)
REAL(r8), PARAMETER :: gauss_w6(6) = (/0.1713244923791704d0, 0.3607615730481386d0, &
    0.4679139345726910d0, 0.4679139345726910d0, 0.3607615730481386d0, 0.1713244923791704d0/)
!
REAL(r8), PARAMETER :: gauss_a8(8) = (/-0.9602898564975363d0, -0.7966664774136267d0, &
    -0.5255324099163290d0, -0.1834346424956498d0, 0.1834346424956498d0, &
    0.5255324099163290d0, 0.7966664774136267d0, 0.9602898564975363d0/)
REAL(r8), PARAMETER :: gauss_w8(8) = (/0.1012285362903763d0, 0.2223810344533745d0, &
    0.3137066458778873d0, 0.3626837833783620d0, 0.3626837833783620d0, &
    0.3137066458778873d0, 0.2223810344533745d0, 0.1012285362903763d0/)
!
REAL(r8), PARAMETER :: gauss_a9(9) = (/-0.9681602395076261d0, -0.8360311073266358d0, &
    -0.6133714327005904d0, -0.3242534234038089d0, 0.d0, 0.3242534234038089d0, &
    0.6133714327005904d0, 0.8360311073266358d0, 0.96816023950762610d0/)
REAL(r8), PARAMETER :: gauss_w9(9) = (/0.0812743883615744d0, 0.1806481606948574d0, &
    0.2606106964029354d0, 0.3123470770400029d0, 0.3302393550012598d0, 0.3123470770400029d0, &
    0.2606106964029354d0, 0.1806481606948574d0, 0.0812743883615744d0/)
!
REAL(r8), PARAMETER :: gauss_a11(11) = (/-0.9782286581460570d0, -0.8870625997680953d0, &
    -0.7301520055740494d0, -0.5190961292068118d0, -0.2695431559523450d0, 0.d0, &
    0.2695431559523450d0, 0.5190961292068118d0, 0.7301520055740494d0,  &
    0.8870625997680953d0, 0.9782286581460570d0/)
REAL(r8), PARAMETER :: gauss_w11(11) = (/0.0556685671161737d0, 0.1255803694649046d0, &
    0.1862902109277343d0, 0.2331937645919905d0, 0.2628045445102467d0, 0.2729250867779006d0, &
    0.2628045445102467d0, 0.2331937645919905d0, 0.1862902109277343d0,  &
    0.1255803694649046d0, 0.0556685671161737d0/)
!
REAL(r8), PARAMETER :: gauss_a13(13) = [-0.9841830547185881d0, &
    -0.9175983992229779d0, -0.8015780907333099d0, -0.6423493394403403d0, &
    -0.44849275103644687d0, -0.23045831595513483d0, 0.0d0, 0.23045831595513483d0, &
    0.44849275103644687d0, 0.6423493394403403d0, 0.8015780907333099d0, 0.9175983992229779d0, &
    0.9841830547185881d0]

REAL(r8), PARAMETER :: gauss_w13(13) = [0.040484004765315953d0, &
    0.09212149983772865d0, 0.1388735102197874d0, 0.1781459807619455d0, &
    0.2078160475368884d0, 0.2262831802628971d0, 0.2325515532308738d0, &
    0.2262831802628971d0, 0.2078160475368884d0, 0.1781459807619455d0, &
    0.1388735102197874d0, 0.09212149983772865d0, 0.040484004765315953d0]
!
REAL(r8), PARAMETER :: gauss_a15(15) = [-0.9879925180204854d0, &
    -0.937273392400706d0, -0.8482065834104272d0, -0.7244177313601701d0, &
    -0.5709721726085388d0, -0.3941513470775634d0, -0.20119409399743454d0, &
    0.0d0, 0.20119409399743454d0, 0.3941513470775634d0, 0.5709721726085388d0, &
    0.7244177313601701d0, 0.8482065834104272d0, 0.937273392400706d0, 0.9879925180204854d0]
REAL(r8), PARAMETER :: gauss_w15(15) = [0.030753241996118136d0, &
    0.07036604748810772d0, 0.10715922046717183d0, 0.13957067792615427d0, &
    0.16626920581699398d0, 0.18616100001556216d0, 0.19843148532711147d0, &
    0.2025782419255612d0, 0.19843148532711147d0, 0.18616100001556216d0, &
    0.16626920581699398d0, 0.13957067792615427d0, 0.10715922046717183d0, &
    0.07036604748810772d0, 0.030753241996118136d0]
!
REAL(r8), PARAMETER :: gauss_a17(17) = [-0.9905754753144174d0, &
    -0.9506755217687678d0, -0.8802391537269858d0, -0.7815140038968014d0, &
    -0.6576711592166907d0, -0.5126905370864769d0, -0.35123176345387636d0, &
    -0.17848418149584783d0, 0.0d0, 0.17848418149584783d0, 0.35123176345387636d0, &
    0.5126905370864769d0, 0.6576711592166907d0, 0.7815140038968014d0, &
    0.8802391537269858d0, 0.9506755217687678d0, 0.9905754753144174d0]
REAL(r8), PARAMETER :: gauss_w17(17) = [0.02414830286854747d0, &
    0.055459529373987085d0, 0.0850361483171795d0, 0.11188384719340376d0, &
    0.1351363684685257d0, 0.15404576107681056d0, 0.16800410215644998d0, &
    0.1765627053669927d0, 0.17944647035620662d0, 0.1765627053669927d0, &
    0.16800410215644998d0, 0.15404576107681056d0, 0.1351363684685257d0, &
    0.11188384719340376d0, 0.0850361483171795d0, 0.055459529373987085d0, 0.02414830286854747d0]
!
REAL(r8), PARAMETER :: gauss_a19(19) = [-0.9924068438435844d0, &
    -0.96020815213483d0, -0.9031559036148178d0, -0.8227146565371428d0, &
    -0.7209661773352294d0, -0.600545304661681d0, -0.46457074137596094d0, &
    -0.3165640999636298d0, -0.16035864564022534d0, 0.0d0, 0.16035864564022534d0, &
    0.3165640999636298d0, 0.46457074137596094d0, 0.600545304661681d0, &
    0.7209661773352294d0, 0.8227146565371428d0, 0.9031559036148178d0, &
    0.96020815213483d0, 0.9924068438435844d0]
REAL(r8), PARAMETER :: gauss_w19(19) = [0.01946178822972895d0, &
    0.044814226765700516d0, 0.06904454273764106d0, 0.09149002162244968d0, &
    0.11156664554733368d0, 0.12875396253933566d0, 0.14260670217360596d0, &
    0.152766042065859d0, 0.15896884339395378d0, 0.16105444984878306d0, &
    0.15896884339395378d0, 0.152766042065859d0, 0.14260670217360596d0, &
    0.12875396253933566d0, 0.11156664554733368d0, 0.09149002162244968d0, &
    0.06904454273764106d0, 0.044814226765700516d0, 0.01946178822972895d0]
contains
!------------------------------------------------------------------------------
!> Get 1D quadrature rule for specified order (max = 7)
!------------------------------------------------------------------------------
subroutine set_quad_1d(quad,order)
type(oft_quad_type), intent(inout) :: quad !< Quadrature rule
integer, intent(in) :: order !< Desired quadrature order
DEBUG_STACK_PUSH
! Select quadrature rule from requested order
select case(order)
    case(1)
    quad=create_quad_1d((/0.d0/), (/2.d0/), 1)
    case(2)
    quad=create_quad_1d(gauss_a2,gauss_w2,2)
    case(3)
    quad=create_quad_1d(gauss_a4,gauss_w4,4)
    case(4)
    quad=create_quad_1d(gauss_a4,gauss_w4,4)
    case(5)
    quad=create_quad_1d(gauss_a6,gauss_w6,6)
    case(6)
    quad=create_quad_1d(gauss_a6,gauss_w6,6)
    case(7)
    quad=create_quad_1d(gauss_a8,gauss_w8,8)
    case(8)
    quad=create_quad_1d(gauss_a8,gauss_w8,8)
    case(9)
    quad=create_quad_1d(gauss_a9,gauss_w9,9)
    case(10)
    quad=create_quad_1d(gauss_a11,gauss_w11,11)
    case(11)
    quad=create_quad_1d(gauss_a11,gauss_w11,11)
    case(12)
    quad=create_quad_1d(gauss_a13,gauss_w13,13)
    case(13)
    quad=create_quad_1d(gauss_a13,gauss_w13,13)
    case(14)
    quad=create_quad_1d(gauss_a15,gauss_w15,15)
    case(15)
    quad=create_quad_1d(gauss_a15,gauss_w15,15)
    case(16)
    quad=create_quad_1d(gauss_a17,gauss_w17,17)
    case(17)
    quad=create_quad_1d(gauss_a17,gauss_w17,17)
    case(18)
    quad=create_quad_1d(gauss_a19,gauss_w19,19)
    case(19)
    quad=create_quad_1d(gauss_a19,gauss_w19,19)
    case default
    CALL oft_warn('1-D Quadrature not available, using highest order = 19')
    quad=create_quad_1d(gauss_a19,gauss_w19,19)
end select
DEBUG_STACK_POP
end subroutine set_quad_1d
!------------------------------------------------------------------------------
!> Change domain of integration from [-1,1] to [0,1]
!------------------------------------------------------------------------------
subroutine quad_change_domain(x,w,n)
integer(i4), INTENT(IN) :: n !< Number of quadrature points
real(r8), intent(inout) :: x(n) !< Quadrature points
real(r8), intent(inout) :: w(n) !< Quadrature weights
integer(i4) :: i
DO i=1,n
    w(i)=w(i)/2.d0
    x(i)=x(i)/2.d0 + 0.5d0
END DO
end subroutine quad_change_domain
!------------------------------------------------------------------------------
!> Create 1D quadrature object from base Gaussian quadrature rule
!------------------------------------------------------------------------------
function create_quad_1d(x,w,n) result(quad)
integer(i4), INTENT(IN) :: n !< Number of quadrature points
real(r8), intent(in) :: x(n) !< Quadrature points
real(r8), intent(in) :: w(n) !< Quadrature weights
type(oft_quad_type) :: quad !< Quadrature rule
quad%dim=1
quad%np=n
quad%order=n
allocate(quad%pts(1,quad%np),quad%wts(quad%np))
quad%pts(1,:)=x
quad%wts=w
CALL quad_change_domain(quad%pts(1,:), quad%wts, n)
end function create_quad_1d
!------------------------------------------------------------------------------
!> Get 2D quadrature rule for specified order (max = 11)
!------------------------------------------------------------------------------
subroutine set_quad_2d(quad,order)
type(oft_quad_type), intent(inout) :: quad !< Quadrature rule
integer(i4), intent(in) :: order !< Desired quadrature order
type(oft_quad_type) :: tmp_quad
integer(i4) :: i,j,ind
DEBUG_STACK_PUSH
! Select quadrature rule from requested order
CALL set_quad_1d(tmp_quad, order)
quad%dim=2
quad%np=tmp_quad%np**2
quad%order=order
allocate(quad%pts(2,quad%np),quad%wts(quad%np))
ind=0
DO i=1,tmp_quad%np
  DO j=1,tmp_quad%np
    ind=ind+1
    quad%pts(:,ind)=(/tmp_quad%pts(1,i), tmp_quad%pts(1,j)/)
    quad%wts(ind)=tmp_quad%wts(i)*tmp_quad%wts(j)
  END DO
END DO
DEBUG_STACK_POP
end subroutine set_quad_2d
!------------------------------------------------------------------------------
!> Get 3D quadrature rule for specified order (max = 11)
!!
!! @param[in,out] quad Quadrature rule
!! @param[in] order Quadrature order
!------------------------------------------------------------------------------
subroutine set_quad_3d(quad,order)
type(oft_quad_type), intent(inout) :: quad !< Quadrature rule
integer(i4), intent(in) :: order !< Desired quadrature order
type(oft_quad_type) :: tmp_quad
integer(i4) :: i,j,k,ind
DEBUG_STACK_PUSH
! Select quadrature rule from requested order
CALL set_quad_1d(tmp_quad, order)
quad%dim=3
quad%np=tmp_quad%np**3
quad%order=order
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
ind=0
DO i=1,tmp_quad%np
  DO j=1,tmp_quad%np
    DO k=1,tmp_quad%np
      ind=ind+1
      quad%pts(:,ind)=(/tmp_quad%pts(1,i), tmp_quad%pts(1,j), tmp_quad%pts(1,k)/)
      quad%wts(ind)=tmp_quad%wts(i)*tmp_quad%wts(j)*tmp_quad%wts(k)
    END DO
  END DO
END DO
DEBUG_STACK_POP
end subroutine set_quad_3d
end module oft_gauss_quadrature
!------------------------------------------------------------------------------
!> Definition of 1D, 2D, and 3D quadrature rules
!!
!! Symmetric quadrature rules defined by Zhang et al.
!!
!! "Zhang, L., Cui, T., & Liu, H. (January 01, 2009) A set of symmetric
!! quadrature rules on triangles and tetrahedra. Journal of Computational
!! Mathematics, 27, 1, 89-96."
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_grid
!------------------------------------------------------------------------------
module oft_tet_quadrature
USE oft_base
USE oft_quadrature
implicit none
#include "local.h"
private
public set_quad_1d, set_quad_2d, set_quad_3d, oft_quad_type
contains
!------------------------------------------------------------------------------
!> Get 1D quadrature rule for specified order (max = 7)
!------------------------------------------------------------------------------
subroutine set_quad_1d(quad,order)
type(oft_quad_type), intent(inout) :: quad !< Quadrature rule
integer(i4), intent(in) :: order !< Desired quadrature order
DEBUG_STACK_PUSH
! Select quadrature rule from requested order
select case(order)
  case(1)
    quad=quad_1d_p1()
  case(2)
    quad=quad_1d_p3()
  case(3)
    quad=quad_1d_p3()
  case(4)
    quad=quad_1d_p5()
  case(5)
    quad=quad_1d_p5()
  case(6)
    quad=quad_1d_p7()
  case(7)
    quad=quad_1d_p7()
  case default
    CALL oft_warn('1-D Quadrature not available, using highest order = 7')
    quad=quad_1d_p7()
end select
DEBUG_STACK_POP
end subroutine set_quad_1d
!------------------------------------------------------------------------------
!> Set 1D quadrature rule for 1st order
!------------------------------------------------------------------------------
function quad_1d_p1() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
quad%dim=1
quad%np=1
quad%order=1
allocate(quad%pts(2,quad%np),quad%wts(quad%np))
quad%pts(:,1)=quad_p11()
quad%wts(:)=1.d0
end function quad_1d_p1
!------------------------------------------------------------------------------
!> Set 1D quadrature rule for 3rd order
!------------------------------------------------------------------------------
function quad_1d_p3() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
quad%dim=1
quad%np=2
quad%order=3
allocate(quad%pts(2,quad%np),quad%wts(quad%np))
quad%pts(:,1:2)=quad_p12(.2113248654051871177454256097490212d0)
quad%wts(:)=.5d0
end function quad_1d_p3
!------------------------------------------------------------------------------
!> Set 1D quadrature rule for 5th order
!------------------------------------------------------------------------------
function quad_1d_p5() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
quad%dim=1
quad%np=3
quad%order=5
allocate(quad%pts(2,quad%np),quad%wts(quad%np))
quad%pts(:,1:2)=quad_p12(.1127016653792583114820734600217600d0)
quad%wts(1:2)=5.d0/18.d0
quad%pts(:,3)=quad_p11()
quad%wts(3)=4.d0/9.d0
end function quad_1d_p5
!------------------------------------------------------------------------------
!> Set 1D quadrature rule for 7th order
!------------------------------------------------------------------------------
function quad_1d_p7() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
quad%dim=1
quad%np=4
quad%order=7
allocate(quad%pts(2,quad%np),quad%wts(quad%np))
quad%pts(:,1:2)=quad_p12(.0694318442029737123880267555535953d0)
quad%wts(1:2)=.1739274225687269286865319746109997d0
quad%pts(:,3:4)=quad_p12(.3300094782075718675986671204483776d0)
quad%wts(3:4)=.3260725774312730713134680253890002d0
end function quad_1d_p7
!------------------------------------------------------------------------------
!> Set 1D points to (1/2,1/2)
!------------------------------------------------------------------------------
function quad_p11() result(b)
real(r8) :: b(2) !< Evaluation points
b=1.d0/2.d0
end function
!------------------------------------------------------------------------------
!> Permute 1D points of the form (a,1-a)
!------------------------------------------------------------------------------
function quad_p12(a) result(b)
real(r8), intent(in) :: a !< Parameter 1
real(r8) :: b(2,2) !< Evaluation points from permutation [2,2]
b(:,1)=(/a,1.d0-a/)
b(:,2)=(/1.d0-a,a/)
end function
!------------------------------------------------------------------------------
!> Get 2D quadrature rule for specified order (max = 11)
!------------------------------------------------------------------------------
subroutine set_quad_2d(quad,order)
type(oft_quad_type), intent(inout) :: quad !< Quadrature rule
integer(i4), intent(in) :: order !< Desired quadrature order
DEBUG_STACK_PUSH
! Select quadrature rule from requested order
select case(order)
  case(1)
    quad=quad_2d_p1()
  case(2)
    quad=quad_2d_p3()
  case(3)
    quad=quad_2d_p3()
  case(4)
    quad=quad_2d_p4()
  case(5)
    quad=quad_2d_p5()
  case(6)
    quad=quad_2d_p6()
  case(7)
    quad=quad_2d_p7()
  case(8)
    quad=quad_2d_p8()
  case(9)
    quad=quad_2d_p9()
  case(10)
    quad=quad_2d_p10()
  case(11)
    quad=quad_2d_p11()
  case(12)
    quad=quad_2d_p12()
  case(13)
    quad=quad_2d_p14()
  case(14)
    quad=quad_2d_p14()
  case(15)
    quad=quad_2d_p16()
  case(16)
    quad=quad_2d_p16()
  case(17)
    quad=quad_2d_p18()
  case(18)
    quad=quad_2d_p18()
  case default
    CALL oft_warn('2-D Quadrature not available, using highest order = 18')
    quad=quad_2d_p18()
end select
DEBUG_STACK_POP
end subroutine set_quad_2d
!------------------------------------------------------------------------------
!> Set 2D quadrature rule for 1st order
!------------------------------------------------------------------------------
function quad_2d_p1() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(4) :: i
quad%dim=2
quad%np=1
quad%order=1
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
i=1
CALL quad_add_p21(1.d0,i,quad)
end function quad_2d_p1
!------------------------------------------------------------------------------
!> Set 2D quadrature rule for 2nd order
!------------------------------------------------------------------------------
function quad_2d_p2() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(4) :: i
quad%dim=2
quad%np=3
quad%order=2
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
i=1
CALL quad_add_p22(1.d0/3.d0,1.d0/6.d0,i,quad)
end function quad_2d_p2
!------------------------------------------------------------------------------
!> Set 2D quadrature rule for 3rd order
!------------------------------------------------------------------------------
function quad_2d_p3() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(4) :: i
quad%dim=2
quad%np=6
quad%order=3
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
i=1
CALL quad_add_p22(.2811498024409796482535143227020770d0,.1628828503958919109001618041849063d0,i,quad)
CALL quad_add_p22(.0521835308923536850798190106312564d0,.4779198835675637000000000000000000d0,i,quad)
end function quad_2d_p3
!------------------------------------------------------------------------------
!> Set 2D quadrature rule for 4th order
!------------------------------------------------------------------------------
function quad_2d_p4() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(4) :: i
quad%dim=2
quad%np=6
quad%order=4
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
i=1
CALL quad_add_p22(.2233815896780114656950070084331228d0,.4459484909159648863183292538830519d0,i,quad)
CALL quad_add_p22(.1099517436553218676383263249002105d0,.0915762135097707434595714634022015d0,i,quad)
end function quad_2d_p4
!------------------------------------------------------------------------------
!> Set 2D quadrature rule for 5th order
!------------------------------------------------------------------------------
function quad_2d_p5() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(4) :: i
quad%dim=2
quad%np=7
quad%order=5
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
i=1
CALL quad_add_p22(.1259391805448271525956839455001813d0,.1012865073234563388009873619151238d0,i,quad)
CALL quad_add_p22(.1323941527885061807376493878331519d0,.4701420641051150897704412095134476d0,i,quad)
CALL quad_add_p21(9.d0/40.d0,i,quad)
end function quad_2d_p5
!------------------------------------------------------------------------------
!> Set 2D quadrature rule for 6th order
!------------------------------------------------------------------------------
function quad_2d_p6() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(4) :: i
quad%dim=2
quad%np=12
quad%order=6
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
i=1
CALL quad_add_p22(.0508449063702068169209368091068690d0,.0630890144915022283403316028708192d0,i,quad)
CALL quad_add_p22(.1167862757263793660252896113855794d0,.2492867451709104212916385531070191d0,i,quad)
CALL quad_add_p23(.0828510756183735751935534564204425d0,.0531450498448169473532496716313981d0,.3103524510337844054166077339565522d0,i,quad)
end function quad_2d_p6
!------------------------------------------------------------------------------
!> Set 2D quadrature rule for 7th order
!------------------------------------------------------------------------------
function quad_2d_p7() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(4) :: i
quad%dim=2
quad%np=15
quad%order=7
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
i=1
CALL quad_add_p22(.0135338625156655615668230924525939d0,.0282639241560763402235960069132400d0,i,quad)
CALL quad_add_p22(.0789512544320109813765214502977033d0,.4743113232672225752752252279318165d0,i,quad)
CALL quad_add_p22(.1286079278189060745566555330895234d0,.2411433258498488102541435126703621d0,i,quad)
CALL quad_add_p23(.0561201442833753579166666287467563d0,.7612227480245238000000000000000000d0,.0462708777988089106409255939170205d0,i,quad)
end function quad_2d_p7
!------------------------------------------------------------------------------
!> Set 2D quadrature rule for 8th order
!------------------------------------------------------------------------------
function quad_2d_p8() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(4) :: i
quad%dim=2
quad%np=16
quad%order=8
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
i=1
CALL quad_add_p21(.1443156076777871682510911104890646d0,i,quad)
CALL quad_add_p22(.1032173705347182502817915502921290d0,.1705693077517602066222935014914645d0,i,quad)
CALL quad_add_p22(.0324584976231980803109259283417806d0,.0505472283170309754584235505965989d0,i,quad)
CALL quad_add_p22(.0950916342672846247938961043885843d0,.4592925882927231560288155144941693d0,i,quad)
CALL quad_add_p23(.0272303141744349942648446900739089d0,.2631128296346381134217857862846436d0,.0083947774099576053372138345392944d0,i,quad)
end function quad_2d_p8
!------------------------------------------------------------------------------
!> Set 2D quadrature rule for 9th order
!------------------------------------------------------------------------------
function quad_2d_p9() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(4) :: i
quad%dim=2
quad%np=19
quad%order=9
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
i=1
CALL quad_add_p21(.0971357962827988338192419825072886d0,i,quad)
CALL quad_add_p22(.0313347002271390705368548312872093d0,.4896825191987376277837069248361928d0,i,quad)
CALL quad_add_p22(.0255776756586980312616787985589998d0,.0447295133944527098651065899662764d0,i,quad)
CALL quad_add_p22(.0778275410047742793167393562994040d0,.4370895914929366372699303644353550d0,i,quad)
CALL quad_add_p22(.0796477389272102530328917742640453d0,.1882035356190327302409612804673356d0,i,quad)
CALL quad_add_p23(.0432835393772893772893772893772894d0,.7411985987844980206900798735234238d0,.2219629891607656956751025276931911d0,i,quad)
end function quad_2d_p9
!------------------------------------------------------------------------------
!> Set 2D quadrature rule for 10th order
!------------------------------------------------------------------------------
function quad_2d_p10() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(4) :: i
quad%dim=2
quad%np=25
quad%order=10
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
i=1
CALL quad_add_p21(.0809374287976228802571131238165019d0,i,quad)
CALL quad_add_p22(.0772985880029631216825069823803434d0,.4272731788467755380904427175154472d0,i,quad)
CALL quad_add_p22(.0784576386123717313680939208343967d0,.1830992224486750205215743848502200d0,i,quad)
CALL quad_add_p22(.0174691679959294869176071632906781d0,.4904340197011305874539712223768484d0,i,quad)
CALL quad_add_p22(.0042923741848328280304804020901319d0,.0125724455515805327313290850210413d0,i,quad)
CALL quad_add_p23(.0374688582104676429790207654850445d0,.6542686679200661406665700955876279d0,.3080460016852477000000000000000000d0,i,quad)
CALL quad_add_p23(.0269493525918799596454494795810967d0,.1228045770685592734301298174812812d0,.0333718337393047862408164417747804d0,i,quad)
end function quad_2d_p10
!------------------------------------------------------------------------------
!> Set 2D quadrature rule for 11th order
!------------------------------------------------------------------------------
function quad_2d_p11() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(4) :: i
quad%dim=2
quad%np=28
quad%order=11
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
i=1
CALL quad_add_p21(.0811779602968671595154759687498236d0,i,quad)
CALL quad_add_p22(.0123240435069094941184739010162328d0,.0309383552454307848951950149913047d0,i,quad)
CALL quad_add_p22(.0628280097444101072833394281602940d0,.4364981811341288419176152765599732d0,i,quad)
CALL quad_add_p22(.0122203790493645297552122150039379d0,.4989847637025932662879869838313909d0,i,quad)
CALL quad_add_p22(.0677013489528115099209888618232256d0,.2146881979585943366068758138782509d0,i,quad)
CALL quad_add_p22(.0402196936288516904235668896075687d0,.1136831040421133902052931562283618d0,i,quad)
CALL quad_add_p23(.0147622727177161013362930655877821d0,.8256187661648629043588062003083580d0,.1597423045918501898008607882250075d0,i,quad)
CALL quad_add_p23(.0407279964582990396603369584816179d0,.6404723101348652676770365908189668d0,.3117837157095990000000000000000000d0,i,quad)
end function quad_2d_p11
!------------------------------------------------------------------------------
!> Set 2D quadrature rule for 12th order
!------------------------------------------------------------------------------
function quad_2d_p12() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(4) :: i
quad%dim=2
quad%np=33
quad%order=12
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
i=1
CALL quad_add_p22(.0061662610515590172338664837852304d0,.0213173504532103702468569755157282d0,i,quad)
CALL quad_add_p22(.0628582242178851003542705130928825d0,.2712103850121159223459513403968947d0,i,quad)
CALL quad_add_p22(.0347961129307089429893283972949994d0,.1275761455415859246738963251542836d0,i,quad)
CALL quad_add_p22(.0436925445380384021354572625574750d0,.4397243922944602729797366234843611d0,i,quad)
CALL quad_add_p22(.0257310664404553354177909230715644d0,.4882173897738048825646620652588110d0,i,quad)
CALL quad_add_p23(.0223567732023034457118390767023200d0,.6958360867878034221416355232360725d0,.2813255809899395482481306929745527d0,i,quad)
CALL quad_add_p23(.0173162311086588923716421008110341d0,.8580140335440726305905366166261782d0,.1162519159075971412413541478426018d0,i,quad)
CALL quad_add_p23(.0403715577663809295178286992522368d0,.6089432357797878068561924377637101d0,.2757132696855141939747963460797640d0,i,quad)
end function quad_2d_p12
!------------------------------------------------------------------------------
!> Set 2D quadrature rule for 14th order
!------------------------------------------------------------------------------
function quad_2d_p14() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(4) :: i
quad%dim=2
quad%np=46
quad%order=14
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
i=1
CALL quad_add_p21(.0585962852260285941278938063477560d0,i,quad)
CALL quad_add_p22(.0017351512297252675680618638808094d0,.0099797608064584324152935295820524d0,i,quad)
CALL quad_add_p22(.0261637825586145217778288591819783d0,.4799778935211883898105528650883899d0,i,quad)
CALL quad_add_p22(.0039197292424018290965208275701454d0,.1538119591769669000000000000000000d0,i,quad)
CALL quad_add_p22(.0122473597569408660972869899262505d0,.0740234771169878100000000000000000d0,i,quad)
CALL quad_add_p22(.0281996285032579601073663071515657d0,.1303546825033300000000000000000000d0,i,quad)
CALL quad_add_p22(.0508870871859594852960348275454540d0,.2306172260266531342996053700983831d0,i,quad)
CALL quad_add_p22(.0504534399016035991910208971341189d0,.4223320834191478241144087137913939d0,i,quad)
CALL quad_add_p23(.0170636442122334512900253993849472d0,.7862373859346610033296221140330900d0, &
  .1906163600319009042461432828653034d0,i,quad)
CALL quad_add_p23(.0096834664255066004075209630934194d0,.6305521436606074416224090755688129d0, &
.3623231377435471446183267343597729d0,i,quad)
CALL quad_add_p23(.0363857559284850056220113277642717d0,.6265773298563063142335123137534265d0, &
.2907712058836674150248168174816732d0,i,quad)
CALL quad_add_p23(.0069646633735184124253997225042413d0,.9142099849296254122399670993850469d0, &
.0711657108777507625475924502924336d0,i,quad)
end function quad_2d_p14
!------------------------------------------------------------------------------
!> Set 2D quadrature rule for 16th order
!------------------------------------------------------------------------------
function quad_2d_p16() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(4) :: i
quad%dim=2
quad%np=55
quad%order=16
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
i=1
CALL quad_add_p21(.0480221886803770905518394045805199d0,i,quad)
CALL quad_add_p22(.0147091003068019271034036428618692d0,.0817949831313738726414655931188610d0,i,quad)
CALL quad_add_p22(.0295445865493192559953097267964641d0,.1653006019697796506267619329335566d0,i,quad)
CALL quad_add_p22(.0261250173510883774985975654917156d0,.4685921053494613866946028972966056d0,i,quad)
CALL quad_add_p22(.0027803873523900069750030161386621d0,.0144388134454166826141089566956602d0,i,quad)
CALL quad_add_p22(.0318217730005366495034272900559496d0,.2417842853917833534068944592932077d0,i,quad)
CALL quad_add_p22(.0086458343495096599011737341698489d0,.4953103429877699640654950868774055d0,i,quad)
CALL quad_add_p23(.0143003329044953651466164253682521d0,.6505134026613522994311446848416867d0, &
.3313997445370895565813231681825939d0,i,quad)
CALL quad_add_p23(.0278497772036008299522298734239535d0,.6040112814959970398494041030359670d0, &
.3032471627499421850415521780783469d0,i,quad)
CALL quad_add_p23(.0070416734066360975623701880892807d0,.8021682575747416636168619478116671d0, &
.1880280595212371734441821142939888d0,i,quad)
CALL quad_add_p23(.0178998382599337286017702090758108d0,.7565056064428283965511540757580608d0, &
.1835046685222968636823802774370004d0,i,quad)
CALL quad_add_p23(.0274582003843497630724700381009172d0,.4659384387141181848838107335915464d0, &
.3596459487975046000000000000000100d0,i,quad)
CALL quad_add_p23(.0072997969394317620841125440877777d0,.9063948439920415013624996618653400d0, &
.0771943712957554322825152250527139d0,i,quad)
end function quad_2d_p16
!------------------------------------------------------------------------------
!> Set 2D quadrature rule for 18th order
!------------------------------------------------------------------------------
function quad_2d_p18() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(4) :: i
quad%dim=2
quad%np=72
quad%order=18
allocate(quad%pts(3,quad%np),quad%wts(quad%np))
i=1
CALL quad_add_p22(.0139778616452860209795840079905549d0,.0732708864643828315786196714876895d0,i,quad)
CALL quad_add_p22(.0005549069792132137850684555152509d0,.0039177489832282316427840744195806d0,i,quad)
CALL quad_add_p22(.0210268138197046690284298685162450d0,.4675973189887110616515129966229624d0,i,quad)
CALL quad_add_p22(.0340182121799276997472265274182211d0,.4179162109674113120121268105139935d0,i,quad)
CALL quad_add_p23(.0279101658047749951418434740078169d0,.1653816933602894800544902692391766d0, &
.5636967056608707538051458939380737d0,i,quad)
CALL quad_add_p23(.0182146861271508661267339566206858d0,.2875008944057839899961939131396606d0, &
.2860423261392047491209581074803029d0,i,quad)
CALL quad_add_p23(.0142670236581097930775198241095567d0,.1258893143198247960170648399490380d0, &
.6960432186424611957925748602819539d0,i,quad)
CALL quad_add_p23(.0142371230906750507043127637741560d0,.0632219159465026144935750801169980d0, &
.7605455518876824326145947637978687d0,i,quad)
CALL quad_add_p23(.0192575838546747877991373836820213d0,.0789102274540205177520722103754889d0, &
.5920196312717585633226205754022254d0,i,quad)
CALL quad_add_p23(.0097051322843806411487822763323902d0,.0380580535067857143261189915962621d0, &
.6836812596359998524801240874538131d0,i,quad)
CALL quad_add_p23(.0076297881343321289957824556338534d0,.0142903521304540256499241103130749d0, &
.8517040371370558150285216534427664d0,i,quad)
CALL quad_add_p23(.0106187391363503447944635436705283d0,.0129672723432531723123416343300903d0, &
.5747324928881490288994509386896897d0,i,quad)
CALL quad_add_p23(.0057106698032758388134142143826895d0,.0076485948208408993307926288182273d0, &
.7355104408307292987031352244816406d0,i,quad)
CALL quad_add_p23(.0043268574608764182945223447328327d0,.0127104605722554679311424918135822d0, &
.9393450876437317887074042026828225d0,i,quad)
end function quad_2d_p18
!------------------------------------------------------------------------------
!> Add points to 2D quadrature rule of the form (1/3,1/3,1/3)
!------------------------------------------------------------------------------
subroutine quad_add_p21(wt,i,quad)
real(8), intent(in) :: wt !< Weight value for points
integer(4), intent(inout) :: i !< Starting index to add points
type(oft_quad_type), intent(inout) :: quad !< Quadrature rule
quad%wts(i)=wt
quad%pts(:,i)=1.d0/3.d0
i=i+1
end subroutine quad_add_p21
!------------------------------------------------------------------------------
!> Add points to 2D quadrature rule of the form (a,a,1-2*a)
!------------------------------------------------------------------------------
subroutine quad_add_p22(wt,a,i,quad)
real(8), intent(in) :: wt !< Weight value for points
real(8), intent(in) :: a !< Parameter 1
integer(4), intent(inout) :: i !< Starting index to add points
type(oft_quad_type), intent(inout) :: quad !< Quadrature rule
quad%wts(i:i+2)=wt
quad%pts(:,i)=(/a,a,1.d0-2*a/)
quad%pts(:,i+1)=(/1.d0-2*a,a,a/)
quad%pts(:,i+2)=(/a,1.d0-2*a,a/)
i=i+3
end subroutine quad_add_p22
!------------------------------------------------------------------------------
!> Add points to 2D quadrature rule of the form (a,b,1-a-b)
!------------------------------------------------------------------------------
subroutine quad_add_p23(wt,a,b,i,quad)
real(8), intent(in) :: wt !< Weight value for points
real(8), intent(in) :: a !< Parameter 1
real(8), intent(in) :: b !< Parameter 2
integer(4), intent(inout) :: i !< Starting index to add points
type(oft_quad_type), intent(inout) :: quad !< Quadrature rule
quad%wts(i:i+5)=wt
quad%pts(:,i)=(/a,b,1.d0-a-b/)
quad%pts(:,i+1)=(/1.d0-a-b,a,b/)
quad%pts(:,i+2)=(/b,1.d0-a-b,a/)
quad%pts(:,i+3)=(/a,1.d0-a-b,b/)
quad%pts(:,i+4)=(/b,a,1.d0-a-b/)
quad%pts(:,i+5)=(/1.d0-a-b,b,a/)
i=i+6
end subroutine quad_add_p23
!------------------------------------------------------------------------------
!> Get 3D quadrature rule for specified order (max = 11)
!------------------------------------------------------------------------------
subroutine set_quad_3d(quad,order)
type(oft_quad_type), intent(inout) :: quad !< Quadrature rule
integer(i4), intent(in) :: order !< Desired quadrature order
DEBUG_STACK_PUSH
! Select quadrature rule from requested order
select case(order)
  case(1)
    quad=quad_3d_p1()
  case(2)
    quad=quad_3d_p2()
  case(3)
    quad=quad_3d_p3()
  case(4)
    quad=quad_3d_p4()
  case(5)
    quad=quad_3d_p5()
  case(6)
    quad=quad_3d_p6()
  case(7)
    quad=quad_3d_p7()
  case(8)
    quad=quad_3d_p8()
  case(9)
    quad=quad_3d_p9()
  case(10)
    quad=quad_3d_p10()
  case(11)
    quad=quad_3d_p11()
  case default
    CALL oft_warn('3-D Quadrature not available, using highest order = 11')
    quad=quad_3d_p11()
end select
DEBUG_STACK_POP
end subroutine set_quad_3d
!------------------------------------------------------------------------------
!>  Set 3D quadrature rule for 1st order
!------------------------------------------------------------------------------
function quad_3d_p1() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
quad%dim=3
quad%np=1
quad%order=1
allocate(quad%pts(4,quad%np),quad%wts(quad%np))
quad%pts(:,1:1)=quad_p31()
quad%wts=1.d0
end function quad_3d_p1
!------------------------------------------------------------------------------
!>  Set 3D quadrature rule for 2nd order
!------------------------------------------------------------------------------
function quad_3d_p2() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
quad%dim=3
quad%np=4
quad%order=2
allocate(quad%pts(4,quad%np),quad%wts(quad%np))
quad%pts(:,1:4)=quad_p32(.1381966011250105151795413165634361d0)
quad%wts=1.d0/4.d0
end function quad_3d_p2
!------------------------------------------------------------------------------
!>  Set 3D quadrature rule for 3rd order
!------------------------------------------------------------------------------
function quad_3d_p3() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
quad%dim=3
quad%np=8
quad%order=3
allocate(quad%pts(4,quad%np),quad%wts(quad%np))
quad%pts(:,1:4)=quad_p32(.3280546967114266473358058199811974d0)
quad%wts(1:4)=.1385279665118621423236176983756412d0
quad%pts(:,5:8)=quad_p32(.1069522739329306827717020415706165d0)
quad%wts(5:8)=.1114720334881378576763823016243588d0
end function quad_3d_p3
!------------------------------------------------------------------------------
!>  Set 3D quadrature rule for 4th order
!------------------------------------------------------------------------------
function quad_3d_p4() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
quad%dim=3
quad%np=14
quad%order=4
allocate(quad%pts(4,quad%np),quad%wts(quad%np))
quad%pts(:,1:4)=quad_p32(.0927352503108912262865589206603214d0)
quad%wts(1:4)=.0734930431163619493435869458636788d0
quad%pts(:,5:8)=quad_p32(.3108859192633006097581474949404033d0)
quad%wts(5:8)=.1126879257180158503650149284763889d0
quad%pts(:,9:14)=quad_p33(.0455037041256496500000000000000000d0)
quad%wts(9:14)=.0425460207770814668609320837732882d0
end function quad_3d_p4
!------------------------------------------------------------------------------
!>  Set 3D quadrature rule for 5th order
!------------------------------------------------------------------------------
function quad_3d_p5() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
quad%dim=3
quad%np=14
quad%order=5
allocate(quad%pts(4,quad%np),quad%wts(quad%np))
quad%pts(:,1:4)=quad_p32(.3108859192633006097973457337634578d0)
quad%wts(1:4)=.1126879257180158507991856523332863d0
quad%pts(:,5:8)=quad_p32(.0927352503108912264023239137370306d0)
quad%wts(5:8)=.0734930431163619495437102054863275d0
quad%pts(:,9:14)=quad_p33(.0455037041256496494918805262793394d0)
quad%wts(9:14)=.0425460207770814664380694281202574d0
end function quad_3d_p5
!------------------------------------------------------------------------------
!>  Set 3D quadrature rule for 6th order
!------------------------------------------------------------------------------
function quad_3d_p6() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
quad%dim=3
quad%np=24
quad%order=6
allocate(quad%pts(4,quad%np),quad%wts(quad%np))
quad%pts(:,1:4)=quad_p32(.2146028712591520292888392193862850d0)
quad%wts(1:4)=.0399227502581674920996906275574800d0
quad%pts(:,5:8)=quad_p32(.0406739585346113531155794489564101d0)
quad%wts(5:8)=.0100772110553206429480132374459369d0
quad%pts(:,9:12)=quad_p32(.3223378901422755103439944707624921d0)
quad%wts(9:12)=.0553571815436547220951532778537260d0
quad%pts(:,13:24)=quad_p34(.0636610018750175252992355276057270d0,.6030056647916491413674311390609397d0)
quad%wts(13:24)=.0482142857142857142857142857142857d0
end function quad_3d_p6
!------------------------------------------------------------------------------
!>  Set 3D quadrature rule for 7th order
!------------------------------------------------------------------------------
function quad_3d_p7() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
quad%dim=3
quad%np=35
quad%order=7
allocate(quad%pts(4,quad%np),quad%wts(quad%np))
quad%pts(:,1:1)=quad_p31()
quad%wts(1)=.0954852894641308488605784361172264d0
quad%pts(:,2:5)=quad_p32(.3157011497782027994234299995933115d0)
quad%wts(2:5)=.0423295812099670290762861707985467d0
quad%pts(:,6:11)=quad_p33(.0504898225983963687630538229865625d0)
quad%wts(6:11)=.0318969278328575799342748240829425d0
quad%pts(:,12:23)=quad_p34(.1888338310260010477364311038545858d0,.5751716375870000234832415770223075d0)
quad%wts(12:23)=.0372071307283346213696155611914811d0
quad%pts(:,24:35)=quad_p34(.0212654725414832459888361014998199d0,.8108302410985485611181053798482324d0)
quad%wts(24:35)=.0081107708299033415661034334910965d0
end function quad_3d_p7
!------------------------------------------------------------------------------
!>  Set 3D quadrature rule for 8th order
!------------------------------------------------------------------------------
function quad_3d_p8() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
quad%dim=3
quad%np=46
quad%order=8
allocate(quad%pts(4,quad%np),quad%wts(quad%np))
quad%pts(:,1:4)=quad_p32(.0396754230703899012650713295393895d0)
quad%wts(1:4)=.0063971477799023213214514203351730d0
quad%pts(:,5:8)=quad_p32(.3144878006980963137841605626971483d0)
quad%wts(5:8)=.0401904480209661724881611584798178d0
quad%pts(:,9:12)=quad_p32(.1019866930627033000000000000000000d0)
quad%wts(9:12)=.0243079755047703211748691087719226d0
quad%pts(:,13:16)=quad_p32(.1842036969491915122759464173489092d0)
quad%wts(13:16)=.0548588924136974404669241239903914d0
quad%pts(:,17:22)=quad_p33(.0634362877545398924051412387018983d0)
quad%wts(17:22)=.0357196122340991824649509689966176d0
quad%pts(:,23:34)=quad_p34(.0216901620677280048026624826249302d0,.7199319220394659358894349533527348d0)
quad%wts(23:34)=.0071831906978525394094511052198038d0
quad%pts(:,35:46)=quad_p34(.2044800806367957142413355748727453d0,.5805771901288092241753981713906204d0)
quad%wts(35:46)=.0163721819453191175409381397561191d0
end function quad_3d_p8
!------------------------------------------------------------------------------
!>  Set 3D quadrature rule for 9th order
!------------------------------------------------------------------------------
function quad_3d_p9() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(i4) :: i,ii
quad%dim=3
quad%np=61
quad%order=9
allocate(quad%pts(4,quad%np),quad%wts(quad%np))
ii=0
!---
i=ii; ii=i+1
quad%pts(:,i+1:ii)=quad_p31()
quad%wts(i+1:ii)=.0564266931795062065887150432761254d0
!---
i=ii; ii=i+4
quad%pts(:,i+1:ii)=quad_p32(.0340221770010448664654037088787676d0)
quad%wts(i+1:ii)=.0033410950747134804029997443047177d0
!---
i=ii; ii=i+4
quad%pts(:,i+1:ii)=quad_p32(.3227703335338005253913766832549640d0)
quad%wts(i+1:ii)=.0301137547687737639073142384315749d0
!---
i=ii; ii=i+4
quad%pts(:,i+1:ii)=quad_p32(.0604570774257749300000000000000000d0)
quad%wts(i+1:ii)=.0064909609200615346357621168945686d0
!---
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.4553629909472082118003081504416430d0,.0056831773653301799061001601457447d0)
quad%wts(i+1:ii)=.0098092858682545864319687425925550d0
!---
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.1195022553938258009779737046961144d0,.4631168324784899409762244936577296d0)
quad%wts(i+1:ii)=.0281191538233654725516326174252926d0
!---
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.0280219557834011581550575066541237d0,.7252060768398674887385659542848099d0)
quad%wts(i+1:ii)=.0078945869083315007683414920096088d0
!---
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.1748330320115746157853246459722452d0,.6166825717812564045706830909795407d0)
quad%wts(i+1:ii)=.0194928120472399967169721944892460d0
end function quad_3d_p9
!------------------------------------------------------------------------------
!>  Set 3D quadrature rule for 10th order
!------------------------------------------------------------------------------
function quad_3d_p10() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(i4) :: i,ii
quad%dim=3
quad%np=81
quad%order=10
allocate(quad%pts(4,quad%np),quad%wts(quad%np))
ii=0
!---1
i=ii; ii=i+1
quad%pts(:,i+1:ii)=quad_p31()
quad%wts(i+1:ii)=.0473997735560207383847388211780511d0
!---2
i=ii; ii=i+4
quad%pts(:,i+1:ii)=quad_p32(.3122500686951886477298083186868275d0)
quad%wts(i+1:ii)=.0269370599922686998027641610048821d0
!---3
i=ii; ii=i+4
quad%pts(:,i+1:ii)=quad_p32(.1143096538573461505873711976536504d0)
quad%wts(i+1:ii)=.0098691597167933832345577354301731d0
!---4
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.0061380088247907478475937132484154d0,.9429887673452048661976305869182508d0)
quad%wts(i+1:ii)=.0003619443443392536242398783848085d0
!---5
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.0327794682164426707747210203323242d0,.3401847940871076327889879249496713d0)
quad%wts(i+1:ii)=.0101358716797557927885164701150168d0
!---6
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.4104307392189654942878978442515117d0,.1654860256196110516044901244445264d0)
quad%wts(i+1:ii)=.0113938812201952316236209348807143d0
!---7
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.0324852815648230478355149399784262d0,.1338521522120095130978284359645666d0)
quad%wts(i+1:ii)=.0065761472770359041674557402004507d0
!---8
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.1210501811455894259938950015950505d0,.4771903799042803505441064082969072d0)
quad%wts(i+1:ii)=.0257397319804560712790360122596547d0
!---9
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.1749793421839390242849492265283104d0,.6280718454753660106932760722179097d0)
quad%wts(i+1:ii)=.0129070357988619906392954302494990d0
end function quad_3d_p10
!------------------------------------------------------------------------------
!>  Set 3D quadrature rule for 11th order
!------------------------------------------------------------------------------
function quad_3d_p11() result(quad)
type(oft_quad_type) :: quad !< Quadrature rule
integer(i4) :: i,ii
quad%dim=3
quad%np=109
quad%order=10
allocate(quad%pts(4,quad%np),quad%wts(quad%np))
ii=0
!---1
i=ii; ii=i+1
quad%pts(:,i+1:ii)=quad_p31()
quad%wts(i+1:ii)=.0394321080286588635073303344912044d0
!---2
i=ii; ii=i+4
quad%pts(:,i+1:ii)=quad_p32(.1214913677765337944977023099080722d0)
quad%wts(i+1:ii)=.0156621262272791131500885627687651d0
!---3
i=ii; ii=i+4
quad%pts(:,i+1:ii)=quad_p32(.0323162591510728963539544520895810d0)
quad%wts(i+1:ii)=.0033321723749014081444092361540149d0
!---4
i=ii; ii=i+4
quad%pts(:,i+1:ii)=quad_p32(.3249261497886067978128419024144220d0)
quad%wts(i+1:ii)=.0140260774074897474374913609976924d0
!---5
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.0041483569716600120000000000000100d0,.5982659967901863502054538427761778d0)
quad%wts(i+1:ii)=.0010859075293324663068220983772355d0
!---6
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.2246246106763771414144751511649864d0,.4736622878323495714083696692020524d0)
quad%wts(i+1:ii)=.0202359604306631789111165731654084d0
!---7
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.0519050877725656967442272164426589d0,.5631447779082798987371019763030571d0)
quad%wts(i+1:ii)=.0117902148721258635368493804677018d0
!---8
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.1349301312162402042237591723429930d0,.7083588307858189538569950051271300d0)
quad%wts(i+1:ii)=.0076903149825212959011315780207389d0
!---9
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.0251911921082524729200511850653055d0,.7837195073400773754305740342999090d0)
quad%wts(i+1:ii)=.0044373057034592039047307260214396d0
!---10
i=ii; ii=i+12
quad%pts(:,i+1:ii)=quad_p34(.3653187797817336139693319800988672d0,.1346039083168658000000000000000100d0)
quad%wts(i+1:ii)=.0114295484671840404107705525985940d0
!---11
i=ii; ii=i+24
quad%pts(:,i+1:ii)=quad_p35(.5229075395099384729652169275860292d0,.1407536305436959018425391394912785d0, &
  .0097624381964526155082922803899778d0)
quad%wts(i+1:ii)=.0061856401712178114128192550838953d0
end function quad_3d_p11
!------------------------------------------------------------------------------
!> Set 3D points to (1/4,1/4,1/4,1/4)
!------------------------------------------------------------------------------
function quad_p31() result(b)
real(r8) :: b(4,1) !< Evaluation points
b=1.d0/4.d0
end function quad_p31
!------------------------------------------------------------------------------
!> Permute 3D points of the form (a,a,a,1-3*a)
!------------------------------------------------------------------------------
function quad_p32(a) result(b)
real(r8), intent(in) :: a !< Parameter 1
real(r8) :: b(4,4) !< Evaluation points from permutation
integer :: i
b(:,1)=(/a,a,a,1.d0-3.d0*a/)
do i=1,3
    b((/1,2,3,4/),i+1)=b((/2,3,4,1/),i)
end do
end function quad_p32
!------------------------------------------------------------------------------
!> Permute 3D points of the form (a,a,1/2-a,1/2-a)
!!
!! param[in] a Parameter 1
!! @result Evaluation points from permutation
!------------------------------------------------------------------------------
function quad_p33(a) result(b)
real(r8), intent(in) :: a !< Parameter 1
real(r8) :: b(4,6) !< Evaluation points from permutation
integer :: i
b(:,1)=(/a,a,.5d0-a,.5d0-a/)
do i=1,3
    b((/1,2,3,4/),i+1)=b((/2,3,4,1/),i)
end do
b(:,5)=(/a,.5d0-a,a,.5d0-a/)
b(:,6)=(/.5d0-a,a,.5d0-a,a/)
end function quad_p33
!------------------------------------------------------------------------------
!> Permute 3D points of the form (a,a,b,1-2*a-b)
!------------------------------------------------------------------------------
function quad_p34(a,b) result(c)
real(r8), intent(in) :: a !< Parameter 1
real(r8), intent(in) :: b !< Parameter 2
real(r8) :: c(4,12) !< Evaluation points from permutation
integer :: i
c(:,1)=(/a,a,b,1.d0-2.d0*a-b/)
do i=1,3
    c((/1,2,3,4/),i+1)=c((/2,3,4,1/),i)
end do
c(:,5)=(/a,b,a,1.d0-2.d0*a-b/)
do i=5,7
    c((/1,2,3,4/),i+1)=c((/2,3,4,1/),i)
end do
c(:,9)=(/b,a,a,1.d0-2.d0*a-b/)
do i=9,11
    c((/1,2,3,4/),i+1)=c((/2,3,4,1/),i)
end do
end function quad_p34
!------------------------------------------------------------------------------
!> Permute 3D points of the form (a,b,c,1-a-b-c)
!------------------------------------------------------------------------------
function quad_p35(a,b,c) result(d)
real(r8), intent(in) :: a !< Parameter 1
real(r8), intent(in) :: b !< Parameter 2
real(r8), intent(in) :: c !< Parameter 3
real(r8) :: d(4,24) !< Evaluation points from permutation
integer :: i
d(:,1)=(/a,b,c,1.d0-a-b-c/)
do i=1,3
    d((/1,2,3,4/),i+1)=d((/2,3,4,1/),i)
end do
d(:,5)=(/b,a,c,1.d0-a-b-c/)
do i=5,7
    d((/1,2,3,4/),i+1)=d((/2,3,4,1/),i)
end do
d(:,9)=(/c,b,a,1.d0-a-b-c/)
do i=9,11
    d((/1,2,3,4/),i+1)=d((/2,3,4,1/),i)
end do
d(:,13)=(/1.d0-a-b-c,b,c,a/)
do i=13,15
    d((/1,2,3,4/),i+1)=d((/2,3,4,1/),i)
end do
d(:,17)=(/a,c,b,1.d0-a-b-c/)
do i=17,19
    d((/1,2,3,4/),i+1)=d((/2,3,4,1/),i)
end do
d(:,21)=(/a,b,1.d0-a-b-c,c/)
do i=21,23
    d((/1,2,3,4/),i+1)=d((/2,3,4,1/),i)
end do
end function quad_p35
end module oft_tet_quadrature
