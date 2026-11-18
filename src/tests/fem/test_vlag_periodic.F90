!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_vlag_periodic.F90
!
!> Regression tests for periodic native meshes and vector Lagrange elements on them.
!!
!! A native (tet or hex) mesh is loaded with reflection and periodicity enabled, and the
!! periodicity type on the base and finest levels, the periodic point, edge and face counts on
!! the base level, and the global entity counts on the base and finest levels are written out.
!! Periodicity on the finest (refined and/or distributed) level is carried by merged global
!! indices, so its global counts must follow from the base-level ones by uniform refinement. Two analytic
!! vector fields are then projected onto vector Lagrange elements (mass matrix inverse)
!! and the L2 errors are written out:
!! - A field compatible with the mesh periodicity, which a linear field reproduces
!!   exactly (to solver tolerance) only if periodic DOFs are mapped correctly
!! - A field incompatible with the periodicity, which cannot be reproduced if the
!!   periodic DOFs are actually shared
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup testing
!------------------------------------------------------------------------------
MODULE test_vlag_periodic_fields
USE oft_base
USE fem_utils, ONLY: fem_interp
IMPLICIT NONE
! Analytic vector field used for projection and error evaluation
TYPE, EXTENDS(fem_interp) :: periodic_field
  INTEGER(i4) :: ftype = 1 ! Field selection (see periodic_field_eval)
CONTAINS
  PROCEDURE :: interp => periodic_field_interp
END TYPE periodic_field
CONTAINS
! Evaluate the selected field at a physical position
!   1: Revolved compatible, B = (1+x, y-2z, 2y+z) commutes with rotation about x
!   2: Planar compatible, B = (1+x, y-2x, 3+y) is independent of z
!   3: Incompatible with either periodicity, B = (0, 1+z, 0)
FUNCTION periodic_field_eval(ftype,pt) RESULT(val)
INTEGER(i4), INTENT(in) :: ftype
REAL(r8), INTENT(in) :: pt(3)
REAL(r8) :: val(3)
SELECT CASE(ftype)
  CASE(1)
    val=[1.d0+pt(1), pt(2)-2.d0*pt(3), 2.d0*pt(2)+pt(3)]
  CASE(2)
    val=[1.d0+pt(1), pt(2)-2.d0*pt(1), 3.d0+pt(2)]
  CASE(3)
    val=[0.d0, 1.d0+pt(3), 0.d0]
  CASE DEFAULT
    CALL oft_abort('Unknown field type','periodic_field_eval',__FILE__)
END SELECT
END FUNCTION periodic_field_eval
! Interpolate the selected field in a cell
SUBROUTINE periodic_field_interp(self,cell,f,gop,val)
CLASS(periodic_field), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: cell
REAL(r8), INTENT(in) :: f(:)
REAL(r8), INTENT(in) :: gop(3,4)
REAL(r8), INTENT(out) :: val(:)
val(1:3)=periodic_field_eval(self%ftype,self%mesh%log2phys(cell,f))
END SUBROUTINE periodic_field_interp
END MODULE test_vlag_periodic_fields
!------------------------------------------------------------------------------
PROGRAM test_vlag_periodic
USE oft_base
USE oft_mesh_type, ONLY: oft_mesh
USE oft_mesh_native, ONLY: mesh_native_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: oft_lag_vgetmop, oft_lag_vproject, oft_lag_vrinterp
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE test_vlag_periodic_fields, ONLY: periodic_field
IMPLICIT NONE
INTEGER(i4) :: order=1,ierr,io_unit,good_field
INTEGER(i4) :: nper(3)
INTEGER(i8) :: base_counts(4),fine_counts(4)
REAL(r8) :: err_good,err_bad
CLASS(oft_mesh), POINTER :: base_mesh
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange
TYPE(oft_ml_fem_comp_type), TARGET :: ML_oft_vlagrange
NAMELIST/test_lag_options/order
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_lag_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid (aborts here for unsupported periodic planes)
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_native_id)CALL oft_abort('Wrong mesh type, test for native only.','main',__FILE__)
!---Count periodic entities on the (undistributed) base mesh
base_mesh=>mg_mesh%meshes(1)
nper=0
IF(ASSOCIATED(base_mesh%periodic%lp))nper(1)=COUNT(base_mesh%periodic%lp>0)
IF(ASSOCIATED(base_mesh%periodic%le))nper(2)=COUNT(base_mesh%periodic%le>0)
IF(ASSOCIATED(base_mesh%periodic%lf))nper(3)=COUNT(base_mesh%periodic%lf>0)
!---Global counts (periodic entities merged) on the base and finest levels
base_counts=[base_mesh%global%np,base_mesh%global%ne,base_mesh%global%nf,base_mesh%global%nc]
fine_counts=[mg_mesh%mesh%global%np,mg_mesh%mesh%global%ne,mg_mesh%mesh%global%nf,mg_mesh%mesh%global%nc]
!---Project fields on the finest level
CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,ML_vlag_obj=ML_oft_vlagrange)
good_field=2
IF(mg_mesh%mesh%periodic%revolved)good_field=1
err_good=project_error(good_field)
err_bad=project_error(3)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='periodic.results')
  WRITE(io_unit,'(2I2)')MERGE(1,0,base_mesh%periodic%revolved),MERGE(1,0,mg_mesh%mesh%periodic%revolved)
  WRITE(io_unit,'(3I12)')nper
  WRITE(io_unit,'(3I12)')base_mesh%np,base_mesh%ne,base_mesh%nf
  WRITE(io_unit,'(4I12)')base_counts
  WRITE(io_unit,'(4I12)')fine_counts
  WRITE(io_unit,'(ES24.16)')err_good
  WRITE(io_unit,'(ES24.16)')err_bad
  CLOSE(io_unit)
END IF
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!---------------------------------------------------------------------------------
! Project an analytic field (L2) and return the RMS error of the projection
!---------------------------------------------------------------------------------
FUNCTION project_error(ftype) RESULT(rms_err)
INTEGER(i4), INTENT(in) :: ftype
REAL(r8) :: rms_err
INTEGER(i4) :: i,m
REAL(r8) :: vol,det,goptmp(3,4),bh(3),bexact(3),err_int,vol_int
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: mop => NULL()
TYPE(periodic_field) :: field
TYPE(oft_lag_vrinterp) :: field_h
!---Set FE level and create solver fields
CALL ML_oft_vlagrange%set_level(ML_oft_lagrange%nlevels,propogate=.TRUE.)
CALL ML_oft_vlagrange%vec_create(u)
CALL ML_oft_vlagrange%vec_create(v)
!---Mass matrix inverse
CALL oft_lag_vgetmop(ML_oft_vlagrange%current_level,mop,'none')
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
CALL create_diag_pre(minv%pre)
!---Project field
field%ftype=ftype
field%mesh=>mg_mesh%mesh
CALL oft_lag_vproject(ML_oft_lagrange%current_level,field,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
!---Integrate error
field_h%u=>u
CALL field_h%setup(ML_oft_lagrange%current_level)
err_int=0.d0
vol_int=0.d0
DO i=1,mg_mesh%mesh%nc
  DO m=1,ML_oft_lagrange%current_level%quad%np
    CALL mg_mesh%mesh%jacobian(i,ML_oft_lagrange%current_level%quad%pts(:,m),goptmp,vol)
    det=vol*ML_oft_lagrange%current_level%quad%wts(m)
    CALL field%interp(i,ML_oft_lagrange%current_level%quad%pts(:,m),goptmp,bexact)
    CALL field_h%interp(i,ML_oft_lagrange%current_level%quad%pts(:,m),goptmp,bh)
    err_int=err_int+SUM((bh-bexact)**2)*det
    vol_int=vol_int+det
  END DO
END DO
err_int=oft_mpi_sum(err_int)
vol_int=oft_mpi_sum(vol_int)
rms_err=SQRT(err_int/vol_int)
IF(oft_env%head_proc)WRITE(*,'(A,I2,A,ES12.4)')'Field ',ftype,' RMS projection error = ',rms_err
!---Cleanup
CALL field_h%delete
CALL u%delete
CALL v%delete
DEALLOCATE(u,v)
CALL mop%delete
DEALLOCATE(mop)
CALL minv%pre%delete
DEALLOCATE(minv%pre)
CALL minv%delete
DEALLOCATE(minv)
END FUNCTION project_error
END PROGRAM test_vlag_periodic
