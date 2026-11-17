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
!! swept parameter.
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
INTEGER(i4) :: io_unit,io_unit2,i,iorder,ires,maxpts
INTEGER(i4), PARAMETER :: orders(2) = [5,9]
REAL(r8) :: t
TYPE(oft_quad_type) :: quad
REAL(r8), ALLOCATABLE :: pts(:,:),wts(:)
CALL oft_init
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='gs_cut_cell.results')
  OPEN(NEWUNIT=io_unit2,FILE='gs_cut_cell_whole.results')
  DO iorder=1,2
    CALL set_quad_2d(quad,orders(iorder))
    maxpts=gs_cut_cell_maxpts(quad%np,2)
    ALLOCATE(pts(3,maxpts),wts(maxpts))
    !---Whole cell: the rule itself
    CALL write_rule(0,0,0.d0,[1.d0,2.d0,3.d0],0,RESHAPE([0.d0],[3,0]))
    DO ires=1,2
      DO i=0,nsteps(ires)
        t=REAL(i,r8)/REAL(nsteps(ires),r8)
        !---1: clipping function swept from wholly negative to wholly positive (crosses each vertex)
        CALL write_rule(1,ires,t,[-0.6d0+1.4d0*t+0.3d0,-0.6d0+1.4d0*t-0.1d0,-0.6d0+1.4d0*t-0.5d0], &
          0,RESHAPE([0.d0],[3,0]))
        !---2: boundary reaching an edge (two vertex values pass through zero together)
        CALL write_rule(2,ires,t,[0.5d0,-0.1d0+0.2d0*t,-0.1d0+0.2d0*t],0,RESHAPE([0.d0],[3,0]))
        !---3: boundary pivoting through a vertex
        CALL write_rule(3,ires,t,[-0.1d0+0.2d0*t,0.3d0,-0.3d0],0,RESHAPE([0.d0],[3,0]))
        !---4: half-plane swept across a cut cell
        CALL write_rule(4,ires,t,[0.4d0,0.2d0,-0.3d0],1, &
          RESHAPE([-0.6d0+1.2d0*t+0.2d0,-0.6d0+1.2d0*t-0.3d0,-0.6d0+1.2d0*t+0.1d0],[3,1]))
        !---5: half-plane swept across a cell wholly inside
        CALL write_rule(5,ires,t,[1.d0,1.d0,1.d0],1, &
          RESHAPE([-0.6d0+1.2d0*t-0.2d0,-0.6d0+1.2d0*t+0.3d0,-0.6d0+1.2d0*t+0.05d0],[3,1]))
        !---6: two half-planes and the clipping function moving together
        CALL write_rule(6,ires,t,[-0.5d0+t,0.6d0-0.2d0*t,0.1d0+0.5d0*t],2, &
          RESHAPE([0.3d0-0.8d0*t,-0.2d0+0.6d0*t,0.5d0-0.4d0*t,-0.4d0+0.9d0*t,0.2d0,0.1d0-0.5d0*t],[3,2]))
      END DO
    END DO
    DEALLOCATE(pts,wts)
    CALL quad%delete()
  END DO
  CLOSE(io_unit)
  CLOSE(io_unit2)
END IF
CALL oft_finalize
CONTAINS
!> Evaluate the rule and write one record: order, sweep, resolution, npts, nplanes, t, sum of weights, integrals,
!! clipping-function and half-plane values
SUBROUTINE write_rule(isweep,ires,t,dv,nplanes,hv)
INTEGER(i4), INTENT(in) :: isweep,ires,nplanes
REAL(r8), INTENT(in) :: t,dv(3),hv(:,:)
INTEGER(i4) :: npts,m
REAL(r8) :: ipoly,ismooth,l(3),hout(3,2)
CALL gs_cut_cell_rule(quad%np,quad%pts,quad%wts,dv,nplanes,hv,npts,pts,wts)
ipoly=0.d0
ismooth=0.d0
DO m=1,npts
  l=pts(:,m)
  ipoly=ipoly+wts(m)*(1.d0+2.d0*l(1)-l(2)+3.d0*l(1)**2+l(1)*l(2)-2.d0*l(2)**3+(l(1)*l(2))**2+l(1)**5)
  ismooth=ismooth+wts(m)*(EXP(l(1)-0.5d0*l(2))+1.d0/(1.5d0+l(3)))
END DO
hout=0.d0
IF(nplanes>0)hout(:,1:nplanes)=hv(:,1:nplanes)
WRITE(io_unit,'(5I8,13ES26.17)')orders(iorder),isweep,ires,npts,nplanes,t,SUM(wts(1:npts)),ipoly,ismooth,dv,hout
IF(isweep==0)THEN ! Points and weights of a whole cell must be the cell's rule itself
  DO m=1,npts
    WRITE(io_unit2,'(2I8,4ES26.17)')orders(iorder),m,MAXVAL(ABS(pts(:,m)-quad%pts(:,m))),ABS(wts(m)-quad%wts(m)), &
      quad%wts(m),REAL(quad%np,r8)
  END DO
END IF
END SUBROUTINE write_rule
END PROGRAM test_gs_cut_cell
