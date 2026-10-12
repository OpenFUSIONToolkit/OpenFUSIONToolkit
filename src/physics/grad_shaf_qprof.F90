!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file grad_shaf_qprof.F90
!
!> Cut-cell evaluation of @ref oft_gs::gs_get_qprof outputs (submodule, so that
!! @ref oft_gs can use @ref oft_gs_cutcell, which itself uses @ref oft_gs)
!!
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
SUBMODULE (oft_gs) oft_gs_qprof
USE oft_gs_cutcell, ONLY: gs_fsa_cutcell
IMPLICIT NONE
CONTAINS
!------------------------------------------------------------------------------
!> q, `ravgs`, `fsa_avgs` and `shape_geo` of @ref gs_get_qprof from the cut-cell quadrature
!------------------------------------------------------------------------------
MODULE PROCEDURE gs_qprof_cutcell
REAL(r8) :: qf(nr),avgs(nr,4),fpol(nr),psi_surf
REAL(r8), ALLOCATABLE :: fsa(:,:),shp(:,:)
INTEGER(i4) :: j
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
DO j=1,nr
  IF(.NOT.done(j))CYCLE
  prof(j)=fpol(j)*qf(j)
  IF(PRESENT(ravgs))ravgs(j,:)=avgs(j,:)
  IF(PRESENT(fsa_avgs))fsa_avgs(j,:)=fsa(j,:)
  IF(PRESENT(shape_geo))shape_geo(j,:)=shp(j,:)
END DO
END PROCEDURE gs_qprof_cutcell
END SUBMODULE oft_gs_qprof
