!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
! Regression test for solid (resistive wall) regions in xMHD_lag. A unit box is split into
! a plasma region (region 1, x < 1/2) and a solid region (region 2, x > 1/2). The plasma
! starts with a uniform axial field, density and temperature, and a velocity perturbation
! V = v0 sin(2 pi x) sin(pi y) y_hat in the plasma region. The system is advanced and the
! maximum velocity at Lagrange nodes in the solid region and the RMS velocity over the box
! (both normalized by v0) are written for the driver.
MODULE test_solid_fields
USE oft_base
USE fem_utils, ONLY: fem_interp
IMPLICIT NONE
TYPE, EXTENDS(fem_interp) :: solid_test_field
  REAL(r8) :: v0 = 1.d0 ! Velocity amplitude [m/s]
CONTAINS
  PROCEDURE :: interp => solid_test_field_interp
END TYPE solid_test_field
CONTAINS
SUBROUTINE solid_test_field_interp(self,cell,f,gop,val)
CLASS(solid_test_field), INTENT(INOUT) :: self
INTEGER(i4), INTENT(IN) :: cell
REAL(r8), INTENT(IN) :: f(:)
REAL(r8), INTENT(IN) :: gop(3,4)
REAL(r8), INTENT(OUT) :: val(:)
REAL(r8) :: pt(3)
pt=self%mesh%log2phys(cell,f)
val(1:3)=0.d0
IF(pt(1)<0.5d0)val(2)=self%v0*SIN(2.d0*pi*pt(1))*SIN(pi*pt(2))
END SUBROUTINE solid_test_field_interp
END MODULE test_solid_fields
!------------------------------------------------------------------------------
PROGRAM test_MUG_solid_lag
USE oft_base
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_vproject, oft_lag_vgetmop, oft_lag_vrinterp
USE diagnostic, ONLY: vec_energy
USE xmhd_lag, ONLY: xmhd_run, xmhd_minlev, xmhd_sub_fields, xmhd_ML_lagrange, xmhd_ML_vlagrange
USE test_solid_fields, ONLY: solid_test_field
IMPLICIT NONE
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: v,b,vel,den,temp
TYPE(xmhd_sub_fields) :: equil_fields
TYPE(solid_test_field), TARGET :: init_field
TYPE(oft_lag_vrinterp), TARGET :: vfield
TYPE(multigrid_mesh) :: mg_mesh
INTEGER(i4) :: i,j,io_unit,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: j_lag
LOGICAL, ALLOCATABLE, DIMENSION(:) :: solid_dof
REAL(r8) :: vol,v_energy,vmax_solid
REAL(r8), POINTER, DIMENSION(:) :: vtmp
INTEGER(i4) :: order = 2
REAL(r8) :: b0 = 0.1d0
REAL(r8) :: n0 = 1.d19
REAL(r8) :: t0 = 10.d0
REAL(r8) :: v0 = 1.d2
NAMELIST/test_solid_options/order,b0,n0,t0,v0
CALL oft_init
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_solid_options,IOSTAT=ierr)
CLOSE(io_unit)
CALL multigrid_construct(mg_mesh)
ALLOCATE(xmhd_ML_lagrange,xmhd_ML_vlagrange)
CALL oft_lag_setup(mg_mesh,order,xmhd_ML_lagrange,ML_vlag_obj=xmhd_ML_vlagrange,minlev=-1)
CALL lag_setup_interp(xmhd_ML_lagrange)
!---Uniform magnetic field and velocity perturbation
init_field%mesh=>mg_mesh%mesh
init_field%v0=v0
CALL xmhd_ML_vlagrange%vec_create(b)
CALL xmhd_ML_vlagrange%vec_create(vel)
CALL b%set(0.d0)
NULLIFY(vtmp)
CALL b%get_local(vtmp,3)
vtmp=b0
CALL b%restore_local(vtmp,3)
DEALLOCATE(vtmp)
NULLIFY(mop)
CALL oft_lag_vgetmop(xmhd_ML_vlagrange%current_level,mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=1.d-12
CALL create_diag_pre(minv%pre)
CALL xmhd_ML_vlagrange%vec_create(v)
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,init_field,v)
CALL vel%set(0.d0)
CALL minv%apply(vel,v)
CALL v%delete
CALL mop%delete
DEALLOCATE(v,mop)
CALL minv%pre%delete
CALL minv%delete
!---Uniform temperature and density
CALL xmhd_ML_lagrange%vec_create(temp)
CALL xmhd_ML_lagrange%vec_create(den)
CALL temp%set(t0)
CALL den%set(n0)
!---Run
xmhd_minlev=xmhd_ML_lagrange%nlevels
oft_env%pm=.FALSE.
equil_fields%B=>b
equil_fields%V=>vel
equil_fields%Ne=>den
equil_fields%Ti=>temp
CALL xmhd_run(equil_fields)
!---Maximum velocity at nodes in the solid region
ALLOCATE(solid_dof(xmhd_ML_lagrange%current_level%ne),j_lag(xmhd_ML_lagrange%current_level%nce))
solid_dof=.FALSE.
DO i=1,mg_mesh%mesh%nc
  IF(mg_mesh%mesh%reg(i)/=2)CYCLE
  CALL xmhd_ML_lagrange%current_level%ncdofs(i,j_lag)
  solid_dof(j_lag)=.TRUE.
END DO
vmax_solid=0.d0
NULLIFY(vtmp)
DO j=1,3
  CALL vel%get_local(vtmp,j)
  vmax_solid=MAX(vmax_solid,MAXVAL(ABS(vtmp),MASK=solid_dof))
END DO
DEALLOCATE(vtmp,solid_dof,j_lag)
vmax_solid=oft_mpi_max(vmax_solid)
!---RMS velocity over the box
vfield%u=>vel
CALL vfield%setup(xmhd_ML_lagrange%current_level)
v_energy=vec_energy(mg_mesh%mesh,vfield,order*2)
vol=oft_mpi_sum(SUM(mg_mesh%mesh%cv))
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='solid.results')
  WRITE(io_unit,*)vmax_solid/v0
  WRITE(io_unit,*)SQRT(v_energy/vol)/v0
  CLOSE(io_unit)
END IF
CALL oft_finalize
END PROGRAM test_MUG_solid_lag
