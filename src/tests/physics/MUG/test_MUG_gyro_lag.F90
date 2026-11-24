!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_MUG_gyro_lag.F90
!
!> Regression test for xMHD_lag (all fields on Lagrange elements) gyroviscosity.
!! Same as test_MUG_gyro. A linearly polarized traveling alfven wave is
!! initialized in a triply periodic box with a uniform ion temperature and
!! advanced for one half period. Gyroviscosity splits the
!! two circular polarizations, which rotates the polarization of the wave.
!! The final velocity is projected onto both polarizations for comparison
!! to linear theory by the driver script. The same driver tests the Hall
!! term, which also splits the circular polarizations, and can advance the
!! wave with either the nonlinear or the linearized model.
!!
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_MUG_gyro_lag
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
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_vproject, oft_lag_vgetmop
!---Physics
USE oft_vector_inits, ONLY: uniform_field
USE mhd_utils, ONLY: mu0, proton_mass
USE xmhd_lag, ONLY: xmhd_run, xmhd_lin_run, xmhd_minlev, xmhd_taxis, xmhd_sub_fields, &
  xmhd_ML_lagrange, xmhd_ML_vlagrange
USE test_phys_helpers, ONLY: alfven_eig, lag_vec_mode_proj
IMPLICIT NONE
!---Mass matrix solvers
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_solver), POINTER :: lag_minv => NULL()
CLASS(oft_matrix), POINTER :: lag_mop => NULL()
!---Local variables
CLASS(oft_vector), POINTER :: b,vel,u,v,den,temp,db,dvel,dden,dtemp
TYPE(xmhd_sub_fields) :: equil_fields,pert_fields
TYPE(uniform_field) :: z_field
TYPE(alfven_eig), TARGET :: alf_field
TYPE(multigrid_mesh) :: mg_mesh
INTEGER(i4) :: io_unit,ierr
REAL(r8) :: B0,proj(2,3)
INTEGER(i4) :: order = 2
INTEGER(i4) :: minlev = 1
REAL(r8) :: delta = 1.d-3
REAL(r8) :: v_alf = 1.d4
REAL(r8) :: N0 = 1.d19
REAL(r8) :: ti0 = 2.d0
LOGICAL :: linear = .FALSE.
NAMELIST/test_gyro_options/order,minlev,delta,ti0,linear
!------------------------------------------------------------------------------
! Initialize enviroment
!------------------------------------------------------------------------------
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_gyro_options,IOSTAT=ierr)
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
! Create Lagrange vector mass matrix solver (magnetic field)
!------------------------------------------------------------------------------
NULLIFY(mop)
CALL oft_lag_vgetmop(xmhd_ML_vlagrange%current_level,mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=1.d-10
CALL create_diag_pre(minv%pre)
!---
CALL xmhd_ML_vlagrange%vec_create(u)
CALL xmhd_ML_vlagrange%vec_create(v)
CALL xmhd_ML_vlagrange%vec_create(b)
CALL xmhd_ML_vlagrange%vec_create(db)
!------------------------------------------------------------------------------
! Set B = B0*(zhat + dB) from uniform field and alfven wave init (B0*zhat and
! B0*dB are kept separate for linear runs)
!------------------------------------------------------------------------------
B0=v_alf*SQRT(mu0*2.d0*N0*proton_mass)
z_field%val=(/0.d0,0.d0,1.d0/)
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,z_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL b%add(0.d0,B0,u)
!---
alf_field%mesh=>mg_mesh%mesh
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,alf_field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL db%add(0.d0,delta*B0,u)
IF(.NOT.linear)CALL b%add(1.d0,1.d0,db)
!---
CALL u%delete
CALL v%delete
CALL mop%delete
DEALLOCATE(u,v,mop)
CALL minv%pre%delete
CALL minv%delete
!------------------------------------------------------------------------------
! Create Lagrange metric solver
!------------------------------------------------------------------------------
NULLIFY(lag_mop)
CALL oft_lag_vgetmop(xmhd_ML_vlagrange%current_level,lag_mop,"none")
CALL create_cg_solver(lag_minv)
lag_minv%A=>lag_mop
lag_minv%its=-3
lag_minv%atol=1.d-10
CALL create_diag_pre(lag_minv%pre)
!---
CALL xmhd_ML_vlagrange%vec_create(u)
CALL xmhd_ML_vlagrange%vec_create(v)
CALL xmhd_ML_vlagrange%vec_create(vel)
CALL xmhd_ML_vlagrange%vec_create(dvel)
!------------------------------------------------------------------------------
! Set V from alfven wave init (zero equilibrium flow for linear runs)
!------------------------------------------------------------------------------
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,alf_field,v)
CALL u%set(0.d0)
CALL lag_minv%apply(u,v)
CALL dvel%add(0.d0,delta*v_alf,u)
CALL vel%set(0.d0)
IF(.NOT.linear)CALL vel%add(0.d0,1.d0,dvel)
!---
CALL u%delete
CALL v%delete
CALL lag_mop%delete
DEALLOCATE(u,v,lag_mop)
!------------------------------------------------------------------------------
! Run simulation
!------------------------------------------------------------------------------
xmhd_minlev=minlev
xmhd_taxis=2
oft_env%pm=.FALSE.
CALL xmhd_ML_lagrange%vec_create(den)
CALL xmhd_ML_lagrange%vec_create(temp)
CALL den%set(N0)
CALL temp%set(ti0)
equil_fields%B=>b
equil_fields%V=>vel
equil_fields%Ne=>den
equil_fields%Ti=>temp
IF(linear)THEN
  CALL xmhd_ML_lagrange%vec_create(dden)
  CALL xmhd_ML_lagrange%vec_create(dtemp)
  CALL dden%set(0.d0)
  CALL dtemp%set(0.d0)
  pert_fields%B=>db
  pert_fields%V=>dvel
  pert_fields%Ne=>dden
  pert_fields%Ti=>dtemp
  CALL xmhd_lin_run(equil_fields,pert_fields)
  !---Final velocity is the perturbation (no equilibrium flow)
  CALL vel%add(1.d0,1.d0,dvel)
ELSE
  CALL xmhd_run(equil_fields)
END IF
!------------------------------------------------------------------------------
! Project final velocity onto both polarizations (normalized by initial amplitude)
!------------------------------------------------------------------------------
CALL vel%scale(1.d0/(delta*v_alf))
CALL lag_vec_mode_proj(xmhd_ML_lagrange%current_level,mg_mesh%mesh,vel, &
  alf_field%k_dir,alf_field%lam,order*2,proj)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='gyro.results')
  WRITE(io_unit,*)proj(:,1)
  WRITE(io_unit,*)proj(:,2)
  WRITE(io_unit,*)proj(:,3)
  CLOSE(io_unit)
END IF
!---Finalize enviroment
CALL oft_finalize
END PROGRAM test_MUG_gyro_lag
