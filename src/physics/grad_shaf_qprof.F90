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
!> LCFS outputs of @ref gs_get_qprof (`dl`, `rbounds`, `zbounds`) for surface `psi_q1`
!!
!! The surface is sampled at `nlcfs` uniform polar angles. If the plasma is diverted and `psi_q1`
!! is near the LCFS, the points are projected onto the LCFS along rays from the O-point (as the
!! tracer path does), and the bounds are the extreme projected points. Otherwise the bounds are the
!! exact extrema of the FE surface. `dl` is the length of the closed polygon through the points.
!------------------------------------------------------------------------------
MODULE PROCEDURE gs_lcfs_cutcell
integer(4), parameter :: nlcfs=1000
real(8) :: rz(2,nlcfs,1),qf(1),xa(2,4,1),pt(2),dir(2)
logical :: okr(1)
integer(4) :: k
type(oft_lag_brinterp) :: psi_int
ok=.FALSE.
CALL gs_rays_cutcell(gseq,1,[psi_q1],nlcfs,rz,qf,okr,xa)
IF(.NOT.okr(1))RETURN
IF(gseq%diverted.AND.(psi_q1<0.05d0))THEN
  psi_int%u=>gseq%psi
  CALL psi_int%setup(gseq%device%fe_rep)
  DO k=1,nlcfs
    pt=rz(:,k,1)
    dir=(pt-gseq%o_point)/SQRT(SUM((pt-gseq%o_point)**2))
    CALL gs_psi2pt(gseq,gseq%plasma_bounds(1),pt,gseq%o_point,dir,psi_int)
    rz(:,k,1)=pt
  END DO
  CALL psi_int%delete
  rbounds(:,1)=rz(:,MINLOC(rz(1,:,1),DIM=1),1); rbounds(:,2)=rz(:,MAXLOC(rz(1,:,1),DIM=1),1)
  zbounds(:,1)=rz(:,MINLOC(rz(2,:,1),DIM=1),1); zbounds(:,2)=rz(:,MAXLOC(rz(2,:,1),DIM=1),1)
ELSE
  rbounds(:,1)=[-xa(1,1,1),xa(2,1,1)]; rbounds(:,2)=[xa(1,2,1),xa(2,2,1)]
  zbounds(:,1)=[xa(2,3,1),-xa(1,3,1)]; zbounds(:,2)=[xa(2,4,1),xa(1,4,1)]
END IF
dl=SUM(SQRT(SUM((rz(:,:,1)-CSHIFT(rz(:,:,1),1,DIM=2))**2,DIM=1)))
ok=.TRUE.
END PROCEDURE gs_lcfs_cutcell
!------------------------------------------------------------------------------
!> Points of one surface at uniform polar angle about the O-point, starting on the outboard
!! side and running anticlockwise (as the tracer), from the cut-cell quadrature
!------------------------------------------------------------------------------
MODULE PROCEDURE gs_trace_surf_cutcell
integer(4), parameter :: nth=1000
real(8) :: rz(2,nth,1),qf(1)
logical :: ok(1)
npoints=-1
NULLIFY(points)
CALL gs_rays_cutcell(gseq,1,[psi_in],nth,rz,qf,ok)
IF(.NOT.ok(1))RETURN
ALLOCATE(points(2,nth))
points=rz(:,:,1)
npoints=nth
END PROCEDURE gs_trace_surf_cutcell
END SUBMODULE oft_gs_qprof
