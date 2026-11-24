!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
! Regression test for xMHD_lag on a rotationally periodic (revolved) domain. A wedge about
! the x-axis is initialized with a Z-pinch equilibrium: an azimuthal field B = c (0,-z,y),
! with uniform axial current, balanced by the pressure gradient of a quadratic temperature
! profile. Both are exactly represented by quadratic Lagrange elements, so the discrete
! state is an equilibrium only if the azimuthal vectors are correctly rotated across the
! periodic boundaries. The equilibrium is advanced and the RMS velocity (normalized by
! the Alfven speed at the wall) and the relative change in B are written for the driver.
MODULE test_revolved_fields
USE oft_base
USE fem_utils, ONLY: fem_interp
IMPLICIT NONE
TYPE, EXTENDS(fem_interp) :: pinch_field
  CHARACTER(LEN=1) :: field = 'b'
  REAL(r8) :: c = 0.1d0 ! Field gradient [T/m]
  REAL(r8) :: T0 = 1.d0 ! Temperature on axis [eV]
CONTAINS
  PROCEDURE :: interp => pinch_field_interp
END TYPE pinch_field
CONTAINS
SUBROUTINE pinch_field_interp(self,cell,f,gop,val)
CLASS(pinch_field), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: cell
REAL(r8), INTENT(in) :: f(:)
REAL(r8), INTENT(in) :: gop(3,4)
REAL(r8), INTENT(out) :: val(:)
REAL(r8) :: pt(3)
pt=self%mesh%log2phys(cell,f)
SELECT CASE(self%field)
  CASE('b')
    val(1:3)=self%c*[0.d0,-pt(3),pt(2)]
  CASE('t')
    val(1)=self%T0*(1.d0-(pt(2)**2+pt(3)**2)/2.d0)
END SELECT
END SUBROUTINE pinch_field_interp
END MODULE test_revolved_fields
!------------------------------------------------------------------------------
PROGRAM test_MUG_revolved_lag
USE oft_base
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_vproject, oft_lag_vgetmop, &
  oft_lag_getmop, oft_lag_project, oft_lag_vrinterp
USE diagnostic, ONLY: vec_energy
USE mhd_utils, ONLY: mu0, proton_mass, elec_charge
USE xmhd_lag, ONLY: xmhd_run, xmhd_minlev, xmhd_sub_fields, xmhd_ML_lagrange, xmhd_ML_vlagrange
USE test_revolved_fields, ONLY: pinch_field
IMPLICIT NONE
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v,b,b_init,vel,den,temp
TYPE(xmhd_sub_fields) :: equil_fields
TYPE(pinch_field), TARGET :: init_field
TYPE(oft_lag_vrinterp), TARGET :: bfield,vfield
TYPE(multigrid_mesh) :: mg_mesh
INTEGER(i4) :: io_unit,ierr
REAL(r8) :: v_alf,vol,b_energy,db_energy,v_energy
INTEGER(i4) :: order = 2
INTEGER(i4) :: minlev = 1
REAL(r8) :: c = 0.1d0
REAL(r8) :: n0 = 1.d19
NAMELIST/test_revolved_options/order,minlev,c,n0
CALL oft_init
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_revolved_options,IOSTAT=ierr)
CLOSE(io_unit)
CALL multigrid_construct(mg_mesh)
IF(.NOT.mg_mesh%mesh%periodic%revolved)CALL oft_abort('Mesh is not revolved','test_MUG_revolved_lag',__FILE__)
ALLOCATE(xmhd_ML_lagrange,xmhd_ML_vlagrange)
CALL oft_lag_setup(mg_mesh,order,xmhd_ML_lagrange,ML_vlag_obj=xmhd_ML_vlagrange,minlev=minlev)
CALL lag_setup_interp(xmhd_ML_lagrange)
!---Pressure balance: n0 e T0 (1+te_factor)/2 = c^2/mu0 with te_factor=1
init_field%mesh=>mg_mesh%mesh
init_field%c=c
init_field%T0=c**2/(mu0*n0*elec_charge)
!---Magnetic field
NULLIFY(mop)
CALL oft_lag_vgetmop(xmhd_ML_vlagrange%current_level,mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=1.d-12
CALL create_diag_pre(minv%pre)
CALL xmhd_ML_vlagrange%vec_create(v)
CALL xmhd_ML_vlagrange%vec_create(b)
CALL xmhd_ML_vlagrange%vec_create(b_init)
CALL xmhd_ML_vlagrange%vec_create(vel)
init_field%field='b'
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,init_field,v)
CALL b%set(0.d0)
CALL minv%apply(b,v)
CALL b_init%add(0.d0,1.d0,b)
CALL vel%set(0.d0)
CALL v%delete
CALL mop%delete
DEALLOCATE(v,mop)
CALL minv%pre%delete
CALL minv%delete
!---Temperature and density
NULLIFY(mop)
CALL oft_lag_getmop(xmhd_ML_lagrange%current_level,mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=1.d-12
CALL create_diag_pre(minv%pre)
CALL xmhd_ML_lagrange%vec_create(v)
CALL xmhd_ML_lagrange%vec_create(temp)
CALL xmhd_ML_lagrange%vec_create(den)
init_field%field='t'
CALL oft_lag_project(xmhd_ML_lagrange%current_level,init_field,v)
CALL temp%set(0.d0)
CALL minv%apply(temp,v)
CALL den%set(n0)
CALL v%delete
CALL mop%delete
DEALLOCATE(v,mop)
!---Run
xmhd_minlev=minlev
oft_env%pm=.FALSE.
equil_fields%B=>b
equil_fields%V=>vel
equil_fields%Ne=>den
equil_fields%Ti=>temp
CALL xmhd_run(equil_fields)
!---RMS velocity normalized by the wall Alfven speed and relative change in B
v_alf=c/SQRT(mu0*n0*2.d0*proton_mass)
vfield%u=>vel
CALL vfield%setup(xmhd_ML_lagrange%current_level)
v_energy=vec_energy(mg_mesh%mesh,vfield,order*2)
bfield%u=>b_init
CALL bfield%setup(xmhd_ML_lagrange%current_level)
b_energy=vec_energy(mg_mesh%mesh,bfield,order*2)
CALL b%add(1.d0,-1.d0,b_init)
bfield%u=>b
CALL bfield%setup(xmhd_ML_lagrange%current_level)
db_energy=vec_energy(mg_mesh%mesh,bfield,order*2)
vol=oft_mpi_sum(SUM(mg_mesh%mesh%cv))
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='revolved.results')
  WRITE(io_unit,*)SQRT(v_energy/vol)/v_alf
  WRITE(io_unit,*)SQRT(db_energy/b_energy)
  CLOSE(io_unit)
END IF
CALL oft_finalize
END PROGRAM test_MUG_revolved_lag
