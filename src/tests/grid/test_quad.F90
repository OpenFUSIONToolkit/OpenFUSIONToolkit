!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_quad.F90
!
!> Regression tests for 3D quadrature rules. Test cases are constructed
!! externally for all polynomials of the form \f$ \alpha x^i + \beta y^j +
!! \gamma z^k \f$, where \f$ i+j+k \le order \f$ upto the maximum quadrature order.
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!---------------------------------------------------------------------------
PROGRAM test_quad
USE oft_base
USE oft_tet_quadrature, ONLY: oft_quad_type, set_quad_1d
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
IMPLICIT NONE
REAL(r8) :: x(3),c(3),y,y0
INTEGER(i4) :: e(3),ntests,i,ierr,io_unit1,io_unit2
TYPE(multigrid_mesh) :: mg_mesh
INTEGER(i4) :: order
NAMELIST/test_quad_options/order
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit1,FILE=oft_env%ifile)
READ(io_unit1,test_quad_options,IOSTAT=ierr)
CLOSE(io_unit1)
IF(ierr<0)CALL oft_abort('No options found in input file.','main',__FILE__)
!---Setup grid
CALL multigrid_construct(mg_mesh)
!---Check mesh type
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test vaild for cube only.','test_quad',__FILE__)
!---Init exponent and coefficients to zero
e=0
c=0.d0
!---Run 1D test cases
OPEN(NEWUNIT=io_unit1,FILE='quad_1d.tests')
OPEN(NEWUNIT=io_unit2,FILE='quad_1d.results')
READ(io_unit1,*)ntests
DO i=1,ntests
  !---Evaluate integral
  READ(io_unit1,'(1I12,2E25.17)')e(1),c(1),y0
  IF(SUM(e)>order)CYCLE
  y=lint_eval(c,e)
  WRITE(io_unit2,*)ABS((y-y0)/y0)
END DO
CLOSE(io_unit1)
CLOSE(io_unit2)
!---Run 2D test cases
OPEN(NEWUNIT=io_unit1,FILE='quad_2d.tests')
OPEN(NEWUNIT=io_unit2,FILE='quad_2d.results')
READ(io_unit1,*)ntests
DO i=1,ntests
  !---Evaluate integral
  READ(io_unit1,'(2I12,3E25.17)')e(1:2),c(1:2),y0
  IF(SUM(e)>order)CYCLE
  y=sint_eval(c,e)
  WRITE(io_unit2,*)ABS((y-y0)/y0)
END DO
CLOSE(io_unit1)
CLOSE(io_unit2)
!---Run 3D test cases
OPEN(NEWUNIT=io_unit1,FILE='quad_3d.tests')
OPEN(NEWUNIT=io_unit2,FILE='quad_3d.results')
READ(io_unit1,*)ntests
DO i=1,ntests
  !---Evaluate integral
  READ(io_unit1,'(3I12,4E25.17)')e,c,y0
  IF(SUM(e)>order)CYCLE
  y=vint_eval(c,e)
  WRITE(io_unit2,*)ABS((y-y0)/y0)
END DO
CLOSE(io_unit1)
CLOSE(io_unit2)
!---Finish run
CALL oft_finalize
!---
CONTAINS
!---------------------------------------------------------------------------
!> Evaluate test polynomial at a point in space
!---------------------------------------------------------------------------
FUNCTION poly_eval(x,c,e) RESULT(y)
REAL(r8), INTENT(in) :: x(3) !< Array of spatial locations [3]
REAL(r8), INTENT(in) :: c(3) !< Array of polynomial coefficients [3]
INTEGER(i4), INTENT(in) :: e(3) !< Array of polynomial exponents [3]
REAL(r8) :: y !< \f$ \sum c_i x_i^{e_i} \f$
INTEGER(i4) :: i
y=0.d0
DO i=1,3
  y=y+c(i)*(x(i)**e(i))
END DO
END FUNCTION poly_eval
!---------------------------------------------------------------------------
!> Evaluate integral of the domain {[0,1],[0,1],[0,1]}
!---------------------------------------------------------------------------
FUNCTION vint_eval(c,e) RESULT(y)
REAL(r8), INTENT(in) :: c(3) !< Array of polynomial coefficients [3]
INTEGER(i4), INTENT(in) :: e(3) !< Array of polynomial exponents [3]
REAL(r8) :: y !< \f$ \int \sum{i=1}^3 c_i x_i^{e_i} dV \f$
REAL(r8) :: v,goptmp(3,4),pt(3)
INTEGER(i4) :: i,m,order
TYPE(oft_quad_type) :: quad
!---
order=MAX(MAXVAL(e),1)
CALL mg_mesh%mesh%quad_rule(order,quad)
!---
y = 0.d0
DO i=1,mg_mesh%mesh%nc
  !---Straight cell test
  CALL mg_mesh%mesh%jacobian(i,quad%pts(:,1),goptmp,v)
  !---Get local reconstructed operators
  DO m=1,quad%np ! Loop over quadrature points
    pt = mg_mesh%mesh%log2phys(i,quad%pts(:,m))
    y = y + poly_eval(pt,c,e)*v*quad%wts(m)
  END DO
END DO
CALL quad%delete
END FUNCTION vint_eval
!---------------------------------------------------------------------------
!> Evaluate integral of the domain {[0,1],[0,1],1}
!---------------------------------------------------------------------------
FUNCTION sint_eval(c,e) RESULT(y)
REAL(r8), INTENT(in) :: c(3) !< Array of polynomial coefficients [3]
INTEGER(i4), INTENT(in) :: e(3) !< Array of polynomial exponents [3]
REAL(r8) :: y !< \f$ \int \sum_{i=1}^2 c_i x_i^{e_i} dV \f$
REAL(r8) :: v,goptmp(3,3),pt(3)
INTEGER(i4) :: i,j,m,order
TYPE(oft_quad_type) :: quad
!---
order=MAX(MAXVAL(e),1)
CALL mg_mesh%smesh%quad_rule(order,quad)
!---
y = 0.d0
DO i=1,mg_mesh%smesh%nc
  j=ABS(mg_mesh%smesh%parent%lf(i))
  IF(.NOT.ALL(mg_mesh%smesh%r(3,mg_mesh%smesh%lc(:,i))==1.d0))CYCLE
  !---Straight cell test
  CALL mg_mesh%smesh%jacobian(i,quad%pts(:,1),goptmp,v)
  !---Get local reconstructed operators
  DO m=1,quad%np ! Loop over quadrature points
    pt = mg_mesh%smesh%log2phys(i,quad%pts(:,m))
    y = y + poly_eval(pt,c,e)*v*quad%wts(m)
  END DO
END DO
CALL quad%delete
END FUNCTION sint_eval
!---------------------------------------------------------------------------
!> Evaluate integral of the domain {[0,1],1,1}
!---------------------------------------------------------------------------
FUNCTION lint_eval(c,e) RESULT(y)
REAL(r8), INTENT(in) :: c(3) !< Array of polynomial coefficients [3]
INTEGER(i4), INTENT(in) :: e(3) !< Array of polynomial exponents [3]
REAL(r8) :: y !< \f$ \int \sum_{i=1}^2 c_i x_i^{e_i} dV \f$
REAL(r8) :: v,goptmp(3,3),pt(3)
INTEGER(i4) :: i,j,m,order
TYPE(oft_quad_type) :: quad
!---
order=MAX(MAXVAL(e),1)
CALL set_quad_1d(quad,order)
!---
pt = 1.d0
y = 0.d0
!---Get local reconstructed operators
DO m=1,quad%np ! Loop over quadrature points
  pt(1) = quad%pts(1,m)
  y = y + poly_eval(pt,c,e)*quad%wts(m)
END DO
CALL quad%delete
END FUNCTION lint_eval
END PROGRAM test_quad
