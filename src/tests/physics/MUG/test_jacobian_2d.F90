!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_jacobian_2d.F90
!
!> Regression test for the 2D xMHD module. The approximate Jacobian assembled by
!! build_approx_jacobian_comp (or _incomp) is compared with a central finite-difference derivative of the
!! nonlinear residual about a pseudo-random state, one field (column block) at a time.
!! The relative error of every row block is printed and the largest is reported.
!!
!! If `check_abort=T` a nonlinear run is started instead, so that the input checks made by
!! run_simulation can be tested. If `check_plot=T` the plotting routine is called instead, which
!! is used to test plotting when no restart files (and so no saved background fields) exist.
!!
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_jacobian_2d
USE oft_base
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
USE oft_la_base, ONLY: oft_vector
USE mhd_utils, ONLY: proton_mass
USE xmhd_2d
IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr,i,ir,ic
REAL(r8) :: rel_err,max_err,diff_norm,fd_norm,jw_norm,f0_norm
REAL(r8), POINTER, DIMENSION(:) :: vals,vals2
REAL(r8) :: scale(7),offset(7)
TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
CLASS(oft_vector), POINTER :: u,w,up,um,fp,fm,jw,f0
!---Runtime options
INTEGER(i4) :: order = 2
LOGICAL :: cyl = .FALSE.
LOGICAL :: incomp = .FALSE.
LOGICAL :: check_abort = .FALSE.
LOGICAL :: check_plot = .FALSE.
REAL(r8) :: b0(3) = 0.d0
REAL(r8) :: dt = 1.d-5
REAL(r8) :: fd_eps = 1.d-6
NAMELIST/xmhd_jac_options/order,cyl,incomp,check_abort,check_plot,b0,dt,fd_eps
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,xmhd_jac_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr/=0)CALL oft_abort('Error reading "xmhd_jac_options" from input file','test_jacobian_2d',__FILE__)
!---------------------------------------------------------------------------
! Setup grid and simulation
!---------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)
mhd_sim%cyl_flag=cyl
mhd_sim%incomp=incomp
CALL mhd_sim%setup(mg_mesh,order)
mhd_sim%chi=1.d3
mhd_sim%nu=1.d3
mhd_sim%D_diff=1.d2
mhd_sim%gamma=5.d0/3.d0
mhd_sim%eta(:,1)=50.d0
mhd_sim%eta(:,2)=80.d0
mhd_sim%den_scale=1.d19
mhd_sim%B_0=b0
mhd_sim%g_accel=1.d8 ! Constant source: must not change the Jacobian
mhd_sim%dt=dt
!---------------------------------------------------------------------------
! Pseudo-random state (n, V, T or p, psi, B_y) and vacuum flux
!---------------------------------------------------------------------------
offset=[1.d0, 0.d0, 0.d0, 0.d0, 10.d0, 0.d0, 1.d-2]
scale=[0.3d0, 1.d4, 1.d4, 1.d4, 3.d0, 1.d-3, 3.d-3]
IF(incomp)THEN
  offset(5)=0.d0
  scale(5)=1.d2
END IF
CALL mhd_sim%fe_rep%vec_create(u)
CALL mhd_sim%fe_rep%vec_create(w)
CALL mhd_sim%fe_rep%vec_create(up)
CALL mhd_sim%fe_rep%vec_create(um)
CALL mhd_sim%fe_rep%vec_create(fp)
CALL mhd_sim%fe_rep%vec_create(fm)
CALL mhd_sim%fe_rep%vec_create(jw)
CALL mhd_sim%fe_rep%vec_create(f0)
DO ic=1,7
  NULLIFY(vals)
  CALL u%get_local(vals,ic)
  DO i=1,SIZE(vals)
    vals(i)=offset(ic)+scale(ic)*(pseudo_rand(i,ic)-0.5d0)
  END DO
  CALL u%restore_local(vals,ic)
  DEALLOCATE(vals)
END DO
NULLIFY(vals)
CALL mhd_sim%psi_vac%get_local(vals)
DO i=1,SIZE(vals)
  vals(i)=1.d-3*(pseudo_rand(i,8)-0.5d0)
END DO
CALL mhd_sim%psi_vac%restore_local(vals)
DEALLOCATE(vals)
!---------------------------------------------------------------------------
! Optionally test input checks
!---------------------------------------------------------------------------
IF(check_abort)THEN
  CALL mhd_sim%u%add(0.d0,1.d0,u)
  mhd_sim%nsteps=1
  CALL mhd_sim%run_simulation()
  CALL oft_finalize
END IF
IF(check_plot)THEN
  CALL xmhd_2d_plot(mhd_sim)
  WRITE(*,'(A)')'PLOT_DONE'
  CALL oft_finalize
END IF
!---------------------------------------------------------------------------
! Compare Jacobian with finite-difference derivative of the residual
!---------------------------------------------------------------------------
mhd_sim%jac_dt=dt
mhd_sim%nlfun%dt=dt
IF(incomp)THEN
  CALL build_approx_jacobian_incomp(mhd_sim,u)
ELSE
  CALL build_approx_jacobian_comp(mhd_sim,u)
END IF
CALL mhd_sim%nlfun%apply(u,f0)
max_err=0.d0
DO ic=1,7
  !---Perturbation in field ic only
  CALL w%set(0.d0)
  NULLIFY(vals)
  CALL w%get_local(vals,ic)
  DO i=1,SIZE(vals)
    vals(i)=scale(ic)*(pseudo_rand(i,10+ic)-0.5d0)
  END DO
  CALL w%restore_local(vals,ic)
  DEALLOCATE(vals)
  !---Central difference of the residual
  CALL up%add(0.d0,1.d0,u)
  CALL up%add(1.d0,fd_eps,w)
  CALL um%add(0.d0,1.d0,u)
  CALL um%add(1.d0,-fd_eps,w)
  CALL mhd_sim%nlfun%apply(up,fp)
  CALL mhd_sim%nlfun%apply(um,fm)
  CALL fp%add(1.d0,-1.d0,fm)
  CALL fp%scale(1.d0/(2.d0*fd_eps))
  !---Jacobian product
  CALL mhd_sim%jacobian%apply(w,jw)
  !---Compare each row block
  DO ir=1,7
    NULLIFY(vals,vals2)
    CALL fp%get_local(vals,ir)
    CALL jw%get_local(vals2,ir)
    fd_norm=SQRT(SUM(vals**2))
    jw_norm=SQRT(SUM(vals2**2))
    diff_norm=SQRT(SUM((vals-vals2)**2))
    DEALLOCATE(vals,vals2)
    NULLIFY(vals)
    CALL f0%get_local(vals,ir)
    f0_norm=SQRT(SUM(vals**2))
    DEALLOCATE(vals)
    ! The floor on the denominator excludes finite-difference round-off in blocks that are (nearly) zero
    rel_err=diff_norm/(fd_norm+jw_norm+1.d-3*f0_norm+TINY(1.d0))
    WRITE(*,'(A,I2,A,I2,A,ES12.4)')'  Block (',ir,',',ic,') relative difference = ',rel_err
    max_err=MAX(max_err,rel_err)
  END DO
END DO
WRITE(*,'(A,ES12.4)')'MAX_REL_DIFF = ',max_err
CALL oft_finalize
CONTAINS
!> Deterministic pseudo-random number in [0,1) for index i and stream k
FUNCTION pseudo_rand(i,k) RESULT(val)
INTEGER(i4), INTENT(in) :: i !< Index
INTEGER(i4), INTENT(in) :: k !< Stream
REAL(r8) :: val
val=MODULO(SIN(REAL(i,8)*12.9898d0+REAL(k,8)*78.233d0)*43758.5453d0,1.d0)
END FUNCTION pseudo_rand
END PROGRAM test_jacobian_2d
