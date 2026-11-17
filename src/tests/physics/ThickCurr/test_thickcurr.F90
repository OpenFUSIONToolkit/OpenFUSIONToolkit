!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
! Regression test driver for the ThickCurr (thick_wall) module. Depending on
! `test_case` either a resistive diffusion eigenmode (vector potential) of a fully conducting box is
! advanced in time ('decay') or the slowest decay times of the configuration are
! computed with the eigenvalue solver ('eigs'). Results are written to
! `thickcurr.results` for comparison with analytic values by the pytest driver.
!------------------------------------------------------------------------------
MODULE test_thickcurr_helpers
USE oft_base
USE fem_utils, ONLY: fem_interp
IMPLICIT NONE
!------------------------------------------------------------------------------
! Diffusion eigenmode in a perfectly conducting box [0,L]^3, A = psi z_hat with
! psi = sin(pi*x/L)*sin(pi*y/L), so n x A = 0 (and B.n = 0) on all walls
!------------------------------------------------------------------------------
TYPE, EXTENDS(fem_interp) :: box_mode
  REAL(r8) :: L = 1.d0 ! Box side length
CONTAINS
  PROCEDURE :: interp => box_mode_interp
END TYPE box_mode
CONTAINS
SUBROUTINE box_mode_interp(self,cell,f,gop,val)
CLASS(box_mode), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: cell
REAL(r8), INTENT(in) :: f(:)
REAL(r8), INTENT(in) :: gop(3,4)
REAL(r8), INTENT(out) :: val(:)
REAL(r8) :: pt(3),k
pt=self%mesh%log2phys(cell,f)
k=pi/self%L
val(1)=0.d0
val(2)=0.d0
val(3)=SIN(k*pt(1))*SIN(k*pt(2))
END SUBROUTINE box_mode_interp
END MODULE test_thickcurr_helpers
!------------------------------------------------------------------------------
PROGRAM test_thickcurr
USE oft_base
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_setup_interp, hcurl_grad_getmop, &
  oft_hcurl_grad_project, oft_hcurl_grad_gzerop, oft_hcurl_grad_zerob
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE thick_wall, ONLY: thickcurr_ML_hcurl, thickcurr_ML_H1grad, thickcurr_ML_hcurl_grad, &
  thickcurr_minlev, thickcurr_setup, thickcurr_run, thickcurr_eigs, thickcurr_plot, &
  thickcurr_diag
USE test_thickcurr_helpers, ONLY: box_mode
IMPLICIT NONE
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_H1,ML_H1grad,ML_hcurl
TYPE(oft_ml_fem_comp_type), TARGET :: ML_hcurl_grad
TYPE(oft_hcurl_grad_gzerop) :: gzerop
TYPE(oft_hcurl_grad_zerob) :: zerob
TYPE(box_mode) :: init_field
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v,b0
INTEGER(i4) :: i,io_unit,ierr
REAL(r8) :: e0,e1,eerr,t_final,rate
REAL(r8), ALLOCATABLE :: tau(:)
!---Input options
INTEGER(i4) :: order = 2
INTEGER(i4) :: minlev = -1
INTEGER(i4) :: nev = 5
REAL(r8) :: box_size = 1.d0
REAL(r8) :: eig_shift = -1.d0
CHARACTER(LEN=8) :: test_case = 'decay'
NAMELIST/test_thickcurr_options/order,minlev,nev,box_size,eig_shift,test_case
!------------------------------------------------------------------------------
! Initialize environment and grid
!------------------------------------------------------------------------------
CALL oft_init
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_thickcurr_options,IOSTAT=ierr)
CLOSE(io_unit)
CALL multigrid_construct(mg_mesh)
IF(minlev<0)minlev=mg_mesh%mgmax+order-1
!------------------------------------------------------------------------------
! Build FE structures
!------------------------------------------------------------------------------
CALL oft_h1_setup(mg_mesh,order+1,ML_H1,minlev=minlev)
CALL h1_setup_interp(ML_H1)
CALL oft_hcurl_setup(mg_mesh,order,ML_hcurl,minlev=minlev)
CALL hcurl_setup_interp(ML_hcurl)
CALL oft_hcurl_grad_setup(ML_hcurl,ML_H1,ML_hcurl_grad,ML_H1grad,minlev)
CALL hcurl_grad_setup_interp(ML_hcurl_grad,ML_H1,create_full=.TRUE.)
thickcurr_ML_hcurl=>ML_hcurl
thickcurr_ML_H1grad=>ML_H1grad
thickcurr_ML_hcurl_grad=>ML_hcurl_grad
thickcurr_minlev=minlev
!------------------------------------------------------------------------------
! Setup ThickCurr (resistivity taken from XML input)
!------------------------------------------------------------------------------
CALL thickcurr_setup()
SELECT CASE(TRIM(test_case))
!------------------------------------------------------------------------------
! Time advance of a single eigenmode
!------------------------------------------------------------------------------
CASE('decay')
  NULLIFY(mop)
  CALL hcurl_grad_getmop(ML_hcurl_grad%current_level,mop,"zerob")
  CALL create_cg_solver(minv)
  minv%A=>mop
  minv%its=-3
  minv%atol=1.d-12
  CALL create_diag_pre(minv%pre)
  CALL ML_hcurl_grad%vec_create(u)
  CALL ML_hcurl_grad%vec_create(v)
  CALL ML_hcurl_grad%vec_create(b0)
  !---Project initial mode (tangential component on boundary is zero)
  init_field%mesh=>mg_mesh%mesh
  init_field%L=box_size
  gzerop%ML_hcurl_grad_rep=>ML_hcurl_grad
  zerob%ML_hcurl_grad_rep=>ML_hcurl_grad
  CALL oft_hcurl_grad_project(ML_hcurl_grad%current_level,init_field,v)
  CALL zerob%apply(v)
  CALL gzerop%apply(v)
  CALL u%set(0.d0)
  CALL minv%apply(u,v)
  CALL b0%add(0.d0,1.d0,u)
  CALL thickcurr_diag(u,e0)
  !---Eigenvalue solve with a different shift (time advance operators must be restored)
  IF(eig_shift>0.d0)THEN
    ALLOCATE(tau(1))
    CALL thickcurr_eigs(1,tau,shift=eig_shift)
    DEALLOCATE(tau)
  END IF
  !---Advance in time
  CALL thickcurr_run(u,t_final=t_final)
  CALL thickcurr_diag(u,e1)
  rate=-LOG(e1/e0)/(2.d0*t_final)
  !---Shape error relative to scaled initial mode
  CALL v%add(0.d0,1.d0,u)
  CALL v%add(1.d0,-SQRT(e1/e0),b0)
  CALL thickcurr_diag(v,eerr)
  IF(oft_env%head_proc)THEN
    OPEN(NEWUNIT=io_unit,FILE='thickcurr.results')
    WRITE(io_unit,*)rate
    WRITE(io_unit,*)SQRT(eerr/e1)
    CLOSE(io_unit)
  END IF
  CALL thickcurr_plot
!------------------------------------------------------------------------------
! Eigenvalue solve
!------------------------------------------------------------------------------
CASE('eigs')
  ALLOCATE(tau(nev))
  IF(eig_shift>0.d0)THEN
    CALL thickcurr_eigs(nev,tau,save_plot=.TRUE.,shift=eig_shift)
  ELSE
    CALL thickcurr_eigs(nev,tau,save_plot=.TRUE.)
  END IF
  IF(oft_env%head_proc)THEN
    OPEN(NEWUNIT=io_unit,FILE='thickcurr.results')
    DO i=1,nev
      WRITE(io_unit,*)tau(i)
    END DO
    CLOSE(io_unit)
  END IF
CASE DEFAULT
  CALL oft_abort('Unknown test case','test_thickcurr',__FILE__)
END SELECT
CALL oft_finalize
END PROGRAM test_thickcurr
