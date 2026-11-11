!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_mapping_jac_2d.F90
!
!> Regression tests for trimesh second derivative (Hessian) transformations.
!! 2D counterpart of test_mapping_jac.F90: quadratic functions of (x,y) are
!! interpolated onto a boundary Lagrange space on a circle or (packed) square
!! triangle mesh and their Hessians, reconstructed with oft_lag_bg2interp, are
!! compared to the exact values.
!!
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_mapping_jac_2d
USE oft_base
USE oft_quadrature
USE oft_mesh_sphere, ONLY: mesh_sphere_id
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
!---
USE oft_la_base, ONLY: oft_vector
USE fem_base, ONLY: oft_ml_fem_type
USE oft_lag_basis, ONLY: oft_lag_setup, oft_blag_npos, oft_scalar_bfem, &
  oft_2D_lagrange_cast
USE oft_blag_operators, ONLY: oft_lag_brinterp, oft_lag_bg2interp
IMPLICIT NONE
INTEGER(i4) :: xi,xj,ierr,nfail,i,ntests,io_unit
REAL(r8) :: check_vec(6),tol=1.d-6
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_blagrange
INTEGER(i4) :: order
NAMELIST/test_mapping_options/order
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_mapping_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct_surf(mg_mesh)
IF((mg_mesh%smesh%cad_type/=mesh_sphere_id).AND.(mg_mesh%smesh%cad_type/=mesh_cube_id)) &
  CALL oft_abort('Wrong mesh type, test for CIRCLE or SQUARE only.','main',__FILE__)
IF(oft_env%nprocs>1)CALL oft_abort('Test is for serial meshes only.','main',__FILE__)
!---Treat the (planar) mesh as a true 2D mesh, as for TokaMaker meshes
IF(ANY(ABS(mg_mesh%smesh%r(3,:))>1.d-10))CALL oft_abort('Mesh is not planar.','main',__FILE__)
mg_mesh%smesh%dim=2
!---Setup FEM
CALL oft_lag_setup(mg_mesh,order,ML_blag_obj=ML_oft_blagrange,minlev=-1)
!---Run test cases (same reference data as the 3D test; z terms do not apply in 2D)
nfail=0
OPEN(NEWUNIT=io_unit,FILE='mapping_jac.tests')
READ(io_unit,*)ntests
DO i=1,ntests
  READ(io_unit,'(2I12,6E25.17)')xi,xj,check_vec
  IF((xi==3).OR.(xj==3))CYCLE
  nfail=nfail+check_jac2(xi,xj,check_vec,tol)
END DO
CLOSE(io_unit)
OPEN(NEWUNIT=io_unit,FILE='mapping_jac.results')
WRITE(io_unit,*)nfail
CLOSE(io_unit)
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
! Interpolate x_i*x_j (x_0 = 1) onto the Lagrange space and compare its
! reconstructed Hessian at every quadrature point to the exact result.
!------------------------------------------------------------------------------
FUNCTION check_jac2(xi,xj,check_vec,tol) RESULT(fail_count)
INTEGER(i4), INTENT(in) :: xi,xj
REAL(r8), INTENT(in) :: check_vec(6),tol
!---
INTEGER(i4) :: i,j,k,fail_count
INTEGER(i4), ALLOCATABLE :: dofs(:)
REAL(r8) :: f(3),pt(0:3),gop(3,3),v,val(6)
REAL(r8), POINTER, DIMENSION(:) :: vals
CLASS(oft_vector), POINTER :: u
CLASS(oft_scalar_bfem), POINTER :: lag_rep
TYPE(oft_lag_brinterp), TARGET :: field
TYPE(oft_lag_bg2interp) :: field_d2
TYPE(oft_quad_type), POINTER :: quad
CHARACTER(LEN=1), PARAMETER :: coords(0:3)=(/'s','x','y','z'/)
WRITE(*,*)'Testing ',coords(xi),coords(xj)
!---
fail_count=0
IF(.NOT.oft_2D_lagrange_cast(lag_rep,ML_oft_blagrange%current_level))CALL oft_abort("Invalid FE type","check_jac2",__FILE__)
quad=>lag_rep%quad
!---Set nodal values to x_i*x_j
CALL ML_oft_blagrange%vec_create(u)
NULLIFY(vals)
CALL u%get_local(vals)
ALLOCATE(dofs(lag_rep%nce))
pt=1.d0
DO i=1,mg_mesh%smesh%nc
  CALL lag_rep%ncdofs(i,dofs)
  DO k=1,lag_rep%nce
    CALL oft_blag_npos(lag_rep,i,k,f)
    pt(1:3)=mg_mesh%smesh%log2phys(i,f)
    vals(dofs(k))=pt(xi)*pt(xj)
  END DO
END DO
DEALLOCATE(dofs)
CALL u%restore_local(vals)
DEALLOCATE(vals)
!---Reconstruct Hessian at quadrature points
field%u=>u
CALL field%setup(lag_rep)
CALL field_d2%shared_setup(field)
DO i=1,mg_mesh%smesh%nc
  DO j=1,quad%np
    CALL mg_mesh%smesh%jacobian(i,quad%pts(:,j),gop,v)
    CALL field_d2%interp(i,quad%pts(:,j),gop,val)
    IF(SQRT(SUM((val-check_vec)**2))>tol)fail_count=fail_count+1
  END DO
END DO
CALL field_d2%delete
CALL field%delete
CALL u%delete
DEALLOCATE(u)
END FUNCTION check_jac2
END PROGRAM test_mapping_jac_2d
