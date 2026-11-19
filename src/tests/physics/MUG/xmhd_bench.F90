!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
! Benchmark driver for the 3D xMHD module (not a regression test; driven by
! benchmarks/xmhd_qcache/run_bench.py).
!
! Sets up a nonlinear state on any mesh supported by multigrid_construct (eg. the
! built-in cube generator for tet/hex meshes, or a native curved mesh) from
! analytic fields: uniform guide field plus a divergence-free ABC perturbation,
! a divergence-free sheared flow, and modulated density and temperature(s).
! Physics options are taken from the standard `xmhd_options` group.
!
! mode='profile' calls xmhd_run(...,profile_only=.TRUE.) (timing of Jacobian
! builds, NL residual, matvecs; see xmhd_profile). mode='run' takes `nsteps`
! time steps (from xmhd_options) and then prints the xMHD profiling counters.
! mode='qcheck' calls xmhd_run(...,qcache_check=.TRUE.), which compares the NL
! residual and Jacobian implementations (xmhd_qcache=0,1) on the initial state
! (equivalence check, driven by test_MUG_qcache.py).
!---------------------------------------------------------------------------------
MODULE xmhd_bench_fields
USE oft_base
USE fem_utils, ONLY: fem_interp
IMPLICIT NONE
! Analytic initial condition (physical units)
TYPE, EXTENDS(fem_interp) :: bench_field
  CHARACTER(LEN=2) :: field = 'b' ! Field to evaluate ('b','v','n','ti','te')
  REAL(r8) :: lam(3) = [1.d0,1.d0,2.d0] ! Periods in each direction
  REAL(r8) :: amp0 = 1.d0 ! Background amplitude
  REAL(r8) :: amp1 = 0.1d0 ! Perturbation amplitude
CONTAINS
  PROCEDURE :: interp => bench_field_interp
END TYPE bench_field
CONTAINS
SUBROUTINE bench_field_interp(self,cell,f,gop,val)
CLASS(bench_field), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: cell
REAL(r8), INTENT(in) :: f(:)
REAL(r8), INTENT(in) :: gop(3,4)
REAL(r8), INTENT(out) :: val(:)
REAL(r8) :: pt(3),k(3)
pt=self%mesh%log2phys(cell,f)
k=2.d0*pi/self%lam
SELECT CASE(TRIM(self%field))
CASE('b') ! Guide field + ABC perturbation (divergence free)
  val(1)=self%amp1*(SIN(k(3)*pt(3))+COS(k(2)*pt(2)))
  val(2)=self%amp1*(SIN(k(1)*pt(1))+COS(k(3)*pt(3)))
  val(3)=self%amp0+self%amp1*(SIN(k(2)*pt(2))+COS(k(1)*pt(1)))
CASE('v') ! Sheared, divergence-free flow
  val(1)=self%amp1*SIN(k(2)*pt(2))
  val(2)=self%amp1*SIN(k(3)*pt(3))
  val(3)=self%amp1*SIN(k(1)*pt(1))
CASE('n','ti','te')
  val(1)=self%amp0*(1.d0+self%amp1*COS(k(1)*pt(1))*COS(k(2)*pt(2))*COS(k(3)*pt(3)))
CASE DEFAULT
  CALL oft_abort('Unknown field','bench_field_interp',__FILE__)
END SELECT
END SUBROUTINE bench_field_interp
! Smooth, position dependent transport shaping (used with shaping=T to exercise the shaping paths)
FUNCTION bench_shaping(pt) RESULT(val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8) :: val
val=1.d0+0.5d0*SIN(2.d0*pi*pt(1))*COS(pi*pt(3))
END FUNCTION bench_shaping
! Second shaping function, so that each transport coefficient can use a different one
FUNCTION bench_shaping2(pt) RESULT(val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8) :: val
val=1.d0+0.3d0*COS(2.d0*pi*pt(2))
END FUNCTION bench_shaping2
! Temperature dependent resistivity profile (used with res_prof=T)
FUNCTION bench_res_profile(te) RESULT(val)
REAL(r8), INTENT(in) :: te
REAL(r8) :: val
val=(20.d0/MAX(te,1.d0))**1.5d0
END FUNCTION bench_res_profile
END MODULE xmhd_bench_fields
!---------------------------------------------------------------------------------
PROGRAM xmhd_bench
USE oft_base
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_vproject, oft_lag_vgetmop, &
  oft_lag_project, oft_lag_getmop
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_setup_interp, hcurl_grad_getmop, oft_hcurl_grad_project
USE xmhd, ONLY: xmhd_run, xmhd_minlev, xmhd_taxis, xmhd_sub_fields, xmhd_prof_reps, xmhd_prof_report, &
  xmhd_ML_hcurl, xmhd_ML_H1, xmhd_ML_hcurl_grad, xmhd_ML_H1grad, xmhd_ML_lagrange, xmhd_ML_vlagrange, &
  res_shaping, visc_shaping, diff_shaping, cond_shaping, res_profile, xmhd_den_floor => den_floor, &
  xmhd_temp_floor => temp_floor, xmhd_g_accel => g_accel, xmhd_qcheck_solid, xmhd_qcheck_negsrc, &
  xmhd_qcheck_jac
USE xmhd_bench_fields, ONLY: bench_field, bench_shaping, bench_shaping2, bench_res_profile
IMPLICIT NONE
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: v
TYPE(xmhd_sub_fields) :: ic_fields
TYPE(bench_field), TARGET :: afield
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_timer) :: mytimer
INTEGER(i4) :: io_unit,ierr
REAL(r8) :: run_time
!---Options
INTEGER(i4) :: order = 2 ! Polynomial order (Lagrange; H(Curl)=order, H1=order+1)
INTEGER(i4) :: minlev = -1 ! FE minimum level (-1: finest level only, ie. no multigrid)
CHARACTER(LEN=8) :: mode = 'profile' ! 'profile', 'run' or 'qcheck'
INTEGER(i4) :: prof_reps(3) = [3,5,10] ! Profile repetitions (builds, NL applies/build, matvecs/apply)
REAL(r8) :: lam(3) = [1.d0,1.d0,2.d0] ! Perturbation periods (match the periodic box)
REAL(r8) :: b0 = 0.05d0 ! Guide field [T]
REAL(r8) :: b1 = 0.01d0 ! Field perturbation [T]
REAL(r8) :: v1 = 2.d4 ! Flow amplitude [m/s]
REAL(r8) :: n0 = 1.d19 ! Density [m^-3]
REAL(r8) :: t0 = 20.d0 ! Ion temperature [eV]
REAL(r8) :: te0 = 20.d0 ! Electron temperature [eV]
REAL(r8) :: delta = 0.1d0 ! Relative density/temperature modulation
LOGICAL :: two_temp = .TRUE. ! Evolve separate electron temperature
LOGICAL :: pm = .FALSE. ! Print solver progress
LOGICAL :: shaping = .FALSE. ! Set all four transport shaping functions (res/visc/diff/cond) to bench_shaping
LOGICAL :: shaping_mix = .FALSE. ! With shaping=T: use a second function for diff/cond_shaping
LOGICAL :: res_prof = .FALSE. ! Set a temperature dependent resistivity profile (used with eta_temp>0)
REAL(r8) :: den_floor = -1.d0 ! Density floor [m^-3] (neg to disable)
REAL(r8) :: temp_floor = -1.d0 ! Temperature floor [eV] (neg to disable)
REAL(r8) :: g_accel = 0.d0 ! Gravitational acceleration [m/s^2]
LOGICAL :: qcheck_solid = .FALSE. ! mode='qcheck': mark every fourth cell as solid
LOGICAL :: qcheck_negsrc = .FALSE. ! mode='qcheck': set a floor source pattern
LOGICAL :: jac_check = .TRUE. ! mode='qcheck': also compare the Jacobian implementations on every level
NAMELIST/xmhd_bench_options/order,minlev,mode,prof_reps,lam,b0,b1,v1,n0,t0,te0,delta,two_temp,pm,shaping, &
  shaping_mix,res_prof,den_floor,temp_floor,g_accel,qcheck_solid,qcheck_negsrc,jac_check
!------------------------------------------------------------------------------
CALL oft_init
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,xmhd_bench_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr>0)CALL oft_abort('Error reading "xmhd_bench_options"','xmhd_bench',__FILE__)
IF(oft_env%head_proc)THEN
  WRITE(*,'(A,I4,A,I4,2A)')'XMHD_BENCH order=',order,' minlev=',minlev,' mode=',TRIM(mode)
  WRITE(*,'(A,I6)')'XMHD_BENCH nthreads=',oft_env%nthreads
  WRITE(*,'(A,I6)')'XMHD_BENCH nprocs=',oft_env%nprocs
END IF
!------------------------------------------------------------------------------
! Mesh and FE spaces
!------------------------------------------------------------------------------
CALL multigrid_construct(mg_mesh)
IF((minlev<mg_mesh%mgdim).AND.(minlev/=-1))CALL oft_abort('Only polynomial MG supported','xmhd_bench',__FILE__)
IF(oft_env%head_proc)THEN
  WRITE(*,'(A,I12)')'XMHD_BENCH ncells_global=',mg_mesh%mesh%global%nc
  WRITE(*,'(A,I12)')'XMHD_BENCH ncells_local=',mg_mesh%mesh%nc
  WRITE(*,'(A,I4)')'XMHD_BENCH mesh_order=',mg_mesh%mesh%order
  WRITE(*,'(A,I4)')'XMHD_BENCH cell_np=',mg_mesh%mesh%cell_np
END IF
ALLOCATE(xmhd_ML_lagrange,xmhd_ML_vlagrange)
CALL oft_lag_setup(mg_mesh,order,xmhd_ML_lagrange,ML_vlag_obj=xmhd_ML_vlagrange,minlev=minlev)
CALL lag_setup_interp(xmhd_ML_lagrange)
ALLOCATE(xmhd_ML_H1)
CALL oft_h1_setup(mg_mesh,order+1,xmhd_ML_H1,minlev=minlev)
CALL h1_setup_interp(xmhd_ML_H1)
ALLOCATE(xmhd_ML_hcurl)
CALL oft_hcurl_setup(mg_mesh,order,xmhd_ML_hcurl,minlev=minlev)
CALL hcurl_setup_interp(xmhd_ML_hcurl)
ALLOCATE(xmhd_ML_hcurl_grad,xmhd_ML_H1grad)
CALL oft_hcurl_grad_setup(xmhd_ML_hcurl,xmhd_ML_H1,xmhd_ML_hcurl_grad,xmhd_ML_H1grad,minlev)
CALL hcurl_grad_setup_interp(xmhd_ML_hcurl_grad,xmhd_ML_H1)
afield%mesh=>mg_mesh%mesh
afield%lam=lam
!------------------------------------------------------------------------------
! Magnetic field: L2 projection onto H(Curl) + Grad(H^1)
!------------------------------------------------------------------------------
NULLIFY(mop)
CALL hcurl_grad_getmop(xmhd_ML_hcurl_grad%current_level,mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=1.d-10
CALL create_diag_pre(minv%pre)
CALL xmhd_ML_hcurl_grad%vec_create(v)
CALL xmhd_ML_hcurl_grad%vec_create(ic_fields%B)
afield%field='b'; afield%amp0=b0; afield%amp1=b1
CALL oft_hcurl_grad_project(xmhd_ML_hcurl_grad%current_level,afield,v)
CALL ic_fields%B%set(0.d0)
CALL minv%apply(ic_fields%B,v)
CALL v%delete; DEALLOCATE(v)
CALL cleanup_solver
!------------------------------------------------------------------------------
! Velocity: vector Lagrange projection
!------------------------------------------------------------------------------
CALL oft_lag_vgetmop(xmhd_ML_vlagrange%current_level,mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=1.d-10
CALL create_diag_pre(minv%pre)
CALL xmhd_ML_vlagrange%vec_create(v)
CALL xmhd_ML_vlagrange%vec_create(ic_fields%V)
afield%field='v'; afield%amp0=0.d0; afield%amp1=v1
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,afield,v)
CALL ic_fields%V%set(0.d0)
CALL minv%apply(ic_fields%V,v)
CALL v%delete; DEALLOCATE(v)
CALL cleanup_solver
!------------------------------------------------------------------------------
! Density and temperatures: scalar Lagrange projection
!------------------------------------------------------------------------------
CALL oft_lag_getmop(xmhd_ML_lagrange%current_level,mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=1.d-10
CALL create_diag_pre(minv%pre)
CALL xmhd_ML_lagrange%vec_create(v)
CALL xmhd_ML_lagrange%vec_create(ic_fields%Ne)
CALL xmhd_ML_lagrange%vec_create(ic_fields%Ti)
afield%field='n'; afield%amp0=n0; afield%amp1=delta
CALL scalar_project(ic_fields%Ne)
afield%field='ti'; afield%amp0=t0; afield%amp1=delta
CALL scalar_project(ic_fields%Ti)
IF(two_temp)THEN
  CALL xmhd_ML_lagrange%vec_create(ic_fields%Te)
  afield%field='te'; afield%amp0=te0; afield%amp1=-delta
  CALL scalar_project(ic_fields%Te)
END IF
CALL v%delete; DEALLOCATE(v)
CALL cleanup_solver
!------------------------------------------------------------------------------
! Run
!------------------------------------------------------------------------------
xmhd_minlev=minlev
xmhd_taxis=3
oft_env%pm=pm
xmhd_prof_reps=prof_reps
IF(shaping)THEN
  res_shaping=>bench_shaping
  visc_shaping=>bench_shaping
  diff_shaping=>bench_shaping
  cond_shaping=>bench_shaping
  IF(shaping_mix)THEN
    diff_shaping=>bench_shaping2
    cond_shaping=>bench_shaping2
  END IF
END IF
IF(res_prof)res_profile=>bench_res_profile
xmhd_den_floor=den_floor
xmhd_temp_floor=temp_floor
xmhd_g_accel=g_accel
xmhd_qcheck_solid=qcheck_solid
xmhd_qcheck_negsrc=qcheck_negsrc
xmhd_qcheck_jac=jac_check
CALL mytimer%tick
IF(TRIM(mode)=='profile')THEN
  CALL xmhd_run(ic_fields,profile_only=.TRUE.) ! Calls oft_finalize
ELSE IF(TRIM(mode)=='qcheck')THEN
  CALL xmhd_run(ic_fields,qcache_check=.TRUE.) ! Calls oft_finalize
ELSE
  CALL xmhd_run(ic_fields)
  run_time=mytimer%tock()
  IF(oft_env%head_proc)WRITE(*,'(A,ES12.4)')'XMHD_BENCH run_time_s=',run_time
  run_time=SQRT(ic_fields%B%dot(ic_fields%B))
  IF(oft_env%head_proc)WRITE(*,'(A,ES24.15)')'XMHD_BENCH final_norm_B=',run_time
  run_time=SQRT(ic_fields%V%dot(ic_fields%V))
  IF(oft_env%head_proc)WRITE(*,'(A,ES24.15)')'XMHD_BENCH final_norm_V=',run_time
  run_time=SQRT(ic_fields%Ti%dot(ic_fields%Ti))
  IF(oft_env%head_proc)WRITE(*,'(A,ES24.15)')'XMHD_BENCH final_norm_Ti=',run_time
  CALL xmhd_prof_report
END IF
CALL oft_finalize
CONTAINS
SUBROUTINE scalar_project(x)
CLASS(oft_vector), INTENT(inout) :: x
CALL oft_lag_project(xmhd_ML_lagrange%current_level,afield,v)
CALL x%set(0.d0)
CALL minv%apply(x,v)
END SUBROUTINE scalar_project
SUBROUTINE cleanup_solver
CALL mop%delete
DEALLOCATE(mop)
NULLIFY(mop)
CALL minv%pre%delete
CALL minv%delete
DEALLOCATE(minv)
NULLIFY(minv)
END SUBROUTINE cleanup_solver
END PROGRAM xmhd_bench
