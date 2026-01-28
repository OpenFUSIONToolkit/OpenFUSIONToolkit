!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_blag.F90
!
!> Regression tests for boundary scalar Lagrange finite elements. Tests are
!!performed on all faces of a unit cube at different polynomial orders.
!!
!! The current test cases are:
!! - Solve the Poisson equation \f$ \nabla \cdot \nabla T = 1 \f$
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_blag
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_blag_operators, ONLY: oft_blag_getlop, oft_blag_getmop, oft_blag_zeroe
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
IMPLICIT NONE
INTEGER(i4), PARAMETER :: minlev=2
INTEGER(i4) :: ierr,io_unit
TYPE(xdmf_plot_file) :: plot_file
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange,ML_oft_blagrange
TYPE(oft_ml_fem_comp_type), TARGET :: ML_oft_vlagrange
TYPE(oft_blag_zeroe), TARGET :: blag_zeroe
INTEGER(i4) :: order
NAMELIST/test_blag_options/order
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_blag_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!------------------------------------------------------------------------------
! Setup I/0
!------------------------------------------------------------------------------
CALL plot_file%setup("Test")
CALL mg_mesh%mesh%setup_io(plot_file,order)
!---
CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,ML_blag_obj=ML_oft_blagrange,minlev=minlev)
blag_zeroe%ML_lag_rep=>ML_oft_blagrange
blag_zeroe%parent_geom_flag=>ML_oft_lagrange%current_level%bc
!---Run tests
oft_env%pm=.FALSE.
CALL test_lap
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!---------------------------------------------------------------------------------
!> Solve the Poisson equation \f$ \nabla \cdot \nabla T = 1 \f$ and output
!! required iterataions and final field energy.
!---------------------------------------------------------------------------------
SUBROUTINE test_lap
!---Create solver objects
CLASS(oft_solver), POINTER :: linv => NULL()
!---Local variables
INTEGER(i4) :: i
REAL(r8) :: uu
REAL(r8), POINTER, DIMENSION(:) :: vals => NULL()
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
!---Set FE level
CALL ML_oft_lagrange%set_level(ML_oft_blagrange%nlevels)
!---Create solver fields
CALL ML_oft_blagrange%vec_create(u)
CALL ML_oft_blagrange%vec_create(v)
!---Get FE operators
CALL oft_blag_getlop(ML_oft_blagrange%current_level,lop,'edges',ML_oft_lagrange%current_level%bc)
CALL oft_blag_getmop(ML_oft_blagrange%current_level,mop)
!---Setup matrix solver
CALL create_cg_solver(linv)
linv%A=>lop
linv%its=-3
CALL create_diag_pre(linv%pre)
!---Solve
CALL u%set(1.d0)
CALL mop%apply(u,v)
CALL blag_zeroe%apply(v)
CALL u%set(0.d0)
CALL linv%apply(u,v)
CALL u%get_local(vals)
CALL mg_mesh%smesh%save_vertex_scalar(vals,plot_file,'T')
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='lagrange.results')
  WRITE(io_unit,*)linv%cits
  WRITE(io_unit,*)uu
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
DEALLOCATE(u,v)
!---Destroy matrices
CALL lop%delete
CALL mop%delete
DEALLOCATE(lop,mop)
!---Destory preconditioner
CALL linv%pre%delete
!---Destory solver
CALL linv%delete
END SUBROUTINE test_lap
END PROGRAM test_blag
