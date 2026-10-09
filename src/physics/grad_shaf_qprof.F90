!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file grad_shaf_qprof.F90
!
!> Cut-cell evaluation of @ref oft_gs::gs_get_qprof and @ref oft_gs::gs_trace_surf outputs (submodule, so that
!! @ref oft_gs can use @ref oft_gs_cutcell, which itself uses @ref oft_gs)
!!
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
SUBMODULE (oft_gs) oft_gs_qprof
USE oft_gs_cutcell, ONLY: gs_fsa_cutcell, gs_rays_cutcell
IMPLICIT NONE
CONTAINS
!------------------------------------------------------------------------------
!> q, `ravgs`, `fsa_avgs` and `shape_geo` of @ref gs_get_qprof from the cut-cell quadrature
!------------------------------------------------------------------------------
MODULE PROCEDURE gs_qprof_cutcell
real(8) :: qf(nr),avgs(nr,4),fpol(nr),psi_surf
real(8), allocatable :: fsa(:,:),shp(:,:)
logical :: done(nr)
integer(4) :: j
DO j=1,nr
  psi_surf=psi_q(j)*(gseq%plasma_bounds(2)-gseq%plasma_bounds(1))+gseq%plasma_bounds(1)
  IF(gseq%mode==0)THEN
    fpol(j)=gseq%ffp_scale*gseq%I%f(psi_surf)+gseq%I%f_offset
  ELSE
    fpol(j)=SIGN(1.d0,gseq%I%f_offset)*SQRT(gseq%ffp_scale*gseq%I%f(psi_surf) + gseq%I%f_offset**2)
  END IF
END DO
!---Unallocated outputs are passed as absent (F2008)
IF(PRESENT(fsa_avgs))ALLOCATE(fsa(nr,4))
IF(PRESENT(shape_geo))ALLOCATE(shp(nr,6))
CALL gs_fsa_cutcell(gseq,nr,psi_q,qf,avgs,done,fpol,fsa,shp)
ok=ALL(done)
IF(.NOT.ok)RETURN
prof=fpol*qf
IF(PRESENT(ravgs))ravgs=avgs
IF(PRESENT(fsa_avgs))fsa_avgs=fsa
IF(PRESENT(shape_geo))shape_geo=shp
END PROCEDURE gs_qprof_cutcell
!------------------------------------------------------------------------------
!> Points of one surface at uniform polar angle about the O-point, starting on the outboard
!! side and running anticlockwise (as the tracer), from the cut-cell quadrature
!------------------------------------------------------------------------------
MODULE PROCEDURE gs_trace_surf_cutcell
integer(4), parameter :: nth=1000
real(8) :: theta(nth),rz(2,nth,1),qf(1)
logical :: ok(1)
integer(4) :: k
npoints=-1
NULLIFY(points)
IF(gseq%plasma_bounds(1)<-1.d98)RETURN
DO k=1,nth
  theta(k)=2.d0*pi*(k-1)/REAL(nth,8)
END DO
CALL gs_rays_cutcell(gseq,1,[psi_in],nth,theta,rz,qf,ok)
IF(.NOT.ok(1))RETURN
ALLOCATE(points(2,nth))
points=rz(:,:,1)
npoints=nth
END PROCEDURE gs_trace_surf_cutcell
END SUBMODULE oft_gs_qprof
