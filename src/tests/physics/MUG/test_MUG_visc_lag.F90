!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_MUG_visc_lag.F90
!
!> Regression test for xMHD_lag (all fields on Lagrange elements) viscosity.
!! Same as test_MUG_visc. A traveling sound wave is initialized in a triply
!! periodic box with a uniform magnetic field, which sets the field direction
!! for the parallel viscosity but is not advanced. The wave is
!! advanced for one period and the velocity is projected onto the initial
!! mode to measure viscous damping, which is compared to linear theory by the
!! driver script.
!!
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_MUG_visc_lag
USE oft_base
!--Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_vproject, oft_lag_vgetmop, &
  oft_lag_getmop, oft_lag_project
!---Physics
USE oft_vector_inits, ONLY: uniform_field
USE mhd_utils, ONLY: elec_charge, proton_mass
USE xmhd_lag, ONLY: xmhd_run, xmhd_minlev, xmhd_taxis, temp_floor, xmhd_adv_b, &
  xmhd_sub_fields, xmhd_ML_lagrange, xmhd_ML_vlagrange
USE test_phys_helpers, ONLY: sound_eig, lag_vec_mode_proj
IMPLICIT NONE
!---Mass matrix solvers
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_solver), POINTER :: bminv => NULL()
CLASS(oft_matrix), POINTER :: bmop => NULL()
!---Local variables
CLASS(oft_vector), POINTER :: u,v,b,vel,n,temp
TYPE(xmhd_sub_fields) :: equil_fields
TYPE(sound_eig), TARGET :: sound_field
TYPE(uniform_field) :: b_field
TYPE(multigrid_mesh) :: mg_mesh
INTEGER(i4) :: io_unit,ierr
REAL(r8) :: T0,v_delta,proj(2,3)
INTEGER(i4) :: order = 2
INTEGER(i4) :: minlev = 1
REAL(r8) :: v_sound = 2.d4
REAL(r8) :: N0 = 1.d19
REAL(r8) :: delta = 1.d-3
REAL(r8) :: b_dir(3) = (/0.d0,0.d0,1.d0/)
NAMELIST/test_visc_options/order,minlev,delta,b_dir
!------------------------------------------------------------------------------
! Initialize enviroment
!------------------------------------------------------------------------------
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_visc_options,IOSTAT=ierr)
CLOSE(io_unit)
!------------------------------------------------------------------------------
! Setup grid
!------------------------------------------------------------------------------
CALL multigrid_construct(mg_mesh)
!------------------------------------------------------------------------------
! Build FE structures
!------------------------------------------------------------------------------
!--- Lagrange
ALLOCATE(xmhd_ML_lagrange,xmhd_ML_vlagrange)
CALL oft_lag_setup(mg_mesh,order,xmhd_ML_lagrange,ML_vlag_obj=xmhd_ML_vlagrange,minlev=minlev)
CALL lag_setup_interp(xmhd_ML_lagrange)
!------------------------------------------------------------------------------
! Set uniform B along b_dir (sets bhat for parallel viscosity, not advanced)
!------------------------------------------------------------------------------
NULLIFY(bmop)
CALL oft_lag_vgetmop(xmhd_ML_vlagrange%current_level,bmop,"none")
CALL create_cg_solver(bminv)
bminv%A=>bmop
bminv%its=-3
bminv%atol=1.d-10
CALL create_diag_pre(bminv%pre)
CALL xmhd_ML_vlagrange%vec_create(u)
CALL xmhd_ML_vlagrange%vec_create(v)
CALL xmhd_ML_vlagrange%vec_create(b)
b_field%val=b_dir/SQRT(SUM(b_dir**2))
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,b_field,v)
CALL u%set(0.d0)
CALL bminv%apply(u,v)
CALL b%add(0.d0,1.d0,u)
CALL u%delete
CALL v%delete
CALL bmop%delete
DEALLOCATE(u,v,bmop)
CALL bminv%pre%delete
CALL bminv%delete
!------------------------------------------------------------------------------
! Create Lagrange metric solver
!------------------------------------------------------------------------------
NULLIFY(mop)
CALL oft_lag_getmop(xmhd_ML_lagrange%current_level,mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=1.d-10
CALL create_diag_pre(minv%pre)
!---
CALL xmhd_ML_lagrange%vec_create(u)
CALL xmhd_ML_lagrange%vec_create(v)
CALL xmhd_ML_lagrange%vec_create(n)
CALL xmhd_ML_lagrange%vec_create(temp)
!------------------------------------------------------------------------------
! Set n and T from sound wave init
!------------------------------------------------------------------------------
sound_field%mesh=>mg_mesh%mesh
sound_field%delta=delta
sound_field%field='n'
CALL oft_lag_project(xmhd_ML_lagrange%current_level,sound_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL n%add(0.d0,N0,u)
!---
sound_field%field='t'
CALL oft_lag_project(xmhd_ML_lagrange%current_level,sound_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
T0=(v_sound**2)*3.d0*proton_mass/(5.d0*elec_charge)
CALL temp%add(0.d0,T0,u)
!---
CALL u%delete
CALL v%delete
CALL mop%delete
DEALLOCATE(u,v,mop)
CALL minv%pre%delete
CALL minv%delete
!------------------------------------------------------------------------------
! Create Lagrange vector metric solver
!------------------------------------------------------------------------------
NULLIFY(mop)
CALL oft_lag_vgetmop(xmhd_ML_vlagrange%current_level,mop,"none")
minv%A=>mop
minv%its=-3
minv%atol=1.d-10
CALL create_diag_pre(minv%pre)
!---
CALL xmhd_ML_vlagrange%vec_create(u)
CALL xmhd_ML_vlagrange%vec_create(v)
CALL xmhd_ML_vlagrange%vec_create(vel)
!------------------------------------------------------------------------------
! Set V from sound wave init
!------------------------------------------------------------------------------
sound_field%field='v'
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,sound_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
v_delta=T0*elec_charge/(proton_mass*v_sound)
CALL vel%add(0.d0,v_delta,u)
!---
CALL u%delete
CALL v%delete
CALL mop%delete
DEALLOCATE(u,v,mop)
!------------------------------------------------------------------------------
! Run simulation
!------------------------------------------------------------------------------
xmhd_minlev=minlev
xmhd_taxis=2
xmhd_adv_b=.FALSE.
oft_env%pm=.FALSE.
temp_floor=T0*1.e-2
equil_fields%B=>b
equil_fields%V=>vel
equil_fields%Ne=>n
equil_fields%Ti=>temp
CALL xmhd_run(equil_fields)
!------------------------------------------------------------------------------
! Project final velocity onto the initial mode (normalized by initial amplitude)
!------------------------------------------------------------------------------
CALL vel%scale(1.d0/(v_delta*delta))
CALL lag_vec_mode_proj(xmhd_ML_lagrange%current_level,mg_mesh%mesh,vel, &
  sound_field%k_dir,sound_field%lam,order*2,proj)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='visc.results')
  WRITE(io_unit,*)proj(:,1)
  WRITE(io_unit,*)proj(:,2)
  WRITE(io_unit,*)proj(:,3)
  CLOSE(io_unit)
END IF
!---Finalize enviroment
CALL oft_finalize
END PROGRAM test_MUG_visc_lag
