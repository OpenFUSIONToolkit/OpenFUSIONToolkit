!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_gs_cut_cell.F90
!
!> Regression test for the cut-cell quadrature rule of TokaMaker (@ref oft_gs::gs_cut_cell_rule "gs_cut_cell_rule").
!!
!! The rule is evaluated on one reference triangle while the vertex values of the clipping function (and of up to two
!! half-plane functions) are swept through every change of the clipped region: wholly outside, cut, wholly inside,
!! the boundary passing through a vertex or along an edge, and half-planes entering and leaving the cell. For each step
!! the number of points, the sum of the weights and the integrals of a polynomial and of a smooth non-polynomial
!! function are written; test_gs_cut_cell.py checks exactness against independent clipping and continuity in the
!! swept parameter. Each transition of a sweep (a change of the point count or of a vertex value's sign) is also
!! located by bisection and the rule is written just before and just after it, for a direct check of continuity.
!!
!! @authors Daniel Burgess
!! @date November 2026
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_gs_cut_cell
USE oft_base
USE oft_quadrature, ONLY: oft_quad_type
USE oft_tet_quadrature, ONLY: set_quad_2d
USE oft_gs, ONLY: gs_cut_cell_rule, gs_cut_cell_maxpts
IMPLICIT NONE
INTEGER(i4), PARAMETER :: nsteps(2) = [1000,10000]
INTEGER(i4), PARAMETER :: nsweeps = 6
INTEGER(i4) :: io_unit,io_unit2,io_unit3,i,iorder,ires,isweep,maxpts,nplanes
INTEGER(i4), PARAMETER :: orders(2) = [5,9]
REAL(r8) :: t,dv(3),hv(3,2)
TYPE(oft_quad_type) :: quad
REAL(r8), ALLOCATABLE :: pts(:,:),wts(:)
CALL oft_init
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='gs_cut_cell.results')
  OPEN(NEWUNIT=io_unit2,FILE='gs_cut_cell_whole.results')
  OPEN(NEWUNIT=io_unit3,FILE='gs_cut_cell_jumps.results')
  DO iorder=1,2
    CALL set_quad_2d(quad,orders(iorder))
    maxpts=gs_cut_cell_maxpts(quad%np,2)
    ALLOCATE(pts(3,maxpts),wts(maxpts))
    !---Whole cell: the rule itself
    CALL write_rule(0,0,0.d0,[1.d0,2.d0,3.d0],0,RESHAPE([0.d0],[3,0]))
    DO ires=1,2
      DO i=0,nsteps(ires)
        t=REAL(i,r8)/REAL(nsteps(ires),r8)
        DO isweep=1,nsweeps
          CALL sweep_vals(isweep,t,dv,nplanes,hv)
          CALL write_rule(isweep,ires,t,dv,nplanes,hv(:,1:nplanes))
        END DO
      END DO
    END DO
    !---The rule just before and just after each transition of each sweep
    DO isweep=1,nsweeps
      CALL write_jumps(isweep)
    END DO
    DEALLOCATE(pts,wts)
    CALL quad%delete()
  END DO
  CLOSE(io_unit)
  CLOSE(io_unit2)
  CLOSE(io_unit3)
END IF
CALL oft_finalize
CONTAINS
!> Vertex values of the clipping function and the half-planes of sweep `isweep` at parameter `t`
SUBROUTINE sweep_vals(isweep,t,dv,nplanes,hv)
INTEGER(i4), INTENT(in) :: isweep
REAL(r8), INTENT(in) :: t
REAL(r8), INTENT(out) :: dv(3),hv(3,2)
INTEGER(i4), INTENT(out) :: nplanes
nplanes=0
hv=0.d0
SELECT CASE(isweep)
  CASE(1) ! Clipping function swept from wholly negative to wholly positive (crosses each vertex)
    dv=[-0.6d0+1.4d0*t+0.3d0,-0.6d0+1.4d0*t-0.1d0,-0.6d0+1.4d0*t-0.5d0]
  CASE(2) ! Boundary reaching an edge (two vertex values pass through zero together)
    dv=[0.5d0,-0.1d0+0.2d0*t,-0.1d0+0.2d0*t]
  CASE(3) ! Boundary pivoting through a vertex
    dv=[-0.1d0+0.2d0*t,0.3d0,-0.3d0]
  CASE(4) ! Half-plane swept across a cut cell
    dv=[0.4d0,0.2d0,-0.3d0]
    nplanes=1
    hv(:,1)=[-0.6d0+1.2d0*t+0.2d0,-0.6d0+1.2d0*t-0.3d0,-0.6d0+1.2d0*t+0.1d0]
  CASE(5) ! Half-plane swept across a cell wholly inside
    dv=[1.d0,1.d0,1.d0]
    nplanes=1
    hv(:,1)=[-0.6d0+1.2d0*t-0.2d0,-0.6d0+1.2d0*t+0.3d0,-0.6d0+1.2d0*t+0.05d0]
  CASE(6) ! Two half-planes and the clipping function moving together
    dv=[-0.5d0+t,0.6d0-0.2d0*t,0.1d0+0.5d0*t]
    nplanes=2
    hv=RESHAPE([0.3d0-0.8d0*t,-0.2d0+0.6d0*t,0.5d0-0.4d0*t,-0.4d0+0.9d0*t,0.2d0,0.1d0-0.5d0*t],[3,2])
END SELECT
END SUBROUTINE sweep_vals
!> Point count, sum of the weights and the two integrals of the rule for sweep `isweep` at parameter `t`
SUBROUTINE eval_sweep(isweep,t,vals,sig)
INTEGER(i4), INTENT(in) :: isweep
REAL(r8), INTENT(in) :: t
REAL(r8), INTENT(out) :: vals(4)
INTEGER(i4), INTENT(out) :: sig(10) !< Point count and signs of all vertex values (changes at each transition)
INTEGER(i4) :: npts,np,j
REAL(r8) :: dvs(3),hvs(3,2)
CALL sweep_vals(isweep,t,dvs,np,hvs)
CALL gs_cut_cell_rule(quad%np,quad%pts,quad%wts,dvs,np,hvs(:,1:np),npts,pts,wts)
vals(1)=REAL(npts,r8)
CALL rule_integrals(npts,vals(2),vals(3),vals(4))
sig=0
sig(1)=npts
sig(2:4)=MERGE(1,0,dvs>0.d0)
DO j=1,np
  sig(3*j+2:3*j+4)=MERGE(1,0,hvs(:,j)>0.d0)
END DO
END SUBROUTINE eval_sweep
!> Locate each transition of sweep `isweep` by bisection (to round-off in the parameter) and write one record:
!! order, sweep, parameter, then point count, sum of weights and integrals at the parameter -1E-12 and +1E-12
SUBROUTINE write_jumps(isweep)
INTEGER(i4), INTENT(in) :: isweep
INTEGER(i4), PARAMETER :: nscan = 10000
REAL(r8), PARAMETER :: dt = 1.d-12
INTEGER(i4) :: j,k,sig0(10),sig1(10),sigm(10)
REAL(r8) :: ta,tb,tm,vals(4,2)
CALL eval_sweep(isweep,0.d0,vals(:,1),sig1)
DO j=1,nscan
  sig0=sig1
  ta=REAL(j-1,r8)/REAL(nscan,r8)
  tb=REAL(j,r8)/REAL(nscan,r8)
  CALL eval_sweep(isweep,tb,vals(:,1),sig1)
  IF(ALL(sig0==sig1))CYCLE
  DO k=1,60
    tm=(ta+tb)/2.d0
    CALL eval_sweep(isweep,tm,vals(:,1),sigm)
    IF(ALL(sigm==sig0))THEN
      ta=tm
    ELSE
      tb=tm
    END IF
  END DO
  CALL eval_sweep(isweep,tb-dt,vals(:,1),sigm)
  CALL eval_sweep(isweep,tb+dt,vals(:,2),sigm)
  WRITE(io_unit3,'(2I8,9ES26.17)')orders(iorder),isweep,tb,vals
END DO
END SUBROUTINE write_jumps
!> Sum of the weights and the integrals of a polynomial and a smooth function over the first `npts` points
SUBROUTINE rule_integrals(npts,wsum,ipoly,ismooth)
INTEGER(i4), INTENT(in) :: npts
REAL(r8), INTENT(out) :: wsum,ipoly,ismooth
INTEGER(i4) :: m
REAL(r8) :: l(3)
ipoly=0.d0
ismooth=0.d0
DO m=1,npts
  l=pts(:,m)
  ipoly=ipoly+wts(m)*(1.d0+2.d0*l(1)-l(2)+3.d0*l(1)**2+l(1)*l(2)-2.d0*l(2)**3+(l(1)*l(2))**2+l(1)**5)
  ismooth=ismooth+wts(m)*(EXP(l(1)-0.5d0*l(2))+1.d0/(1.5d0+l(3)))
END DO
wsum=SUM(wts(1:npts))
END SUBROUTINE rule_integrals
!> Evaluate the rule and write one record: order, sweep, resolution, npts, nplanes, t, sum of weights, integrals,
!! clipping-function and half-plane values
SUBROUTINE write_rule(isweep,ires,t,dv,nplanes,hv)
INTEGER(i4), INTENT(in) :: isweep,ires,nplanes
REAL(r8), INTENT(in) :: t,dv(3),hv(:,:)
INTEGER(i4) :: npts,m
REAL(r8) :: wsum,ipoly,ismooth,hout(3,2)
CALL gs_cut_cell_rule(quad%np,quad%pts,quad%wts,dv,nplanes,hv,npts,pts,wts)
CALL rule_integrals(npts,wsum,ipoly,ismooth)
hout=0.d0
IF(nplanes>0)hout(:,1:nplanes)=hv(:,1:nplanes)
WRITE(io_unit,'(5I8,13ES26.17)')orders(iorder),isweep,ires,npts,nplanes,t,wsum,ipoly,ismooth,dv,hout
IF(isweep==0)THEN ! Points and weights of a whole cell must be the cell's rule itself
  DO m=1,npts
    WRITE(io_unit2,'(2I8,4ES26.17)')orders(iorder),m,MAXVAL(ABS(pts(:,m)-quad%pts(:,m))),ABS(wts(m)-quad%wts(m)), &
      quad%wts(m),REAL(quad%np,r8)
  END DO
END IF
END SUBROUTINE write_rule
END PROGRAM test_gs_cut_cell
