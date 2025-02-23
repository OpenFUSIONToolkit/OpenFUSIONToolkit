!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file test_hcurl_sop.F90
!
!> Regression tests for surface operators for vector H1(Curl) finite elements.
!! Tests are performed on a unit cube at different polynomial orders.
!!
!! The current test cases are:
!! - Solve the equation \f$ \nabla \times \nabla \times A = 0 \f$, with
!! \f$ \nabla \times A \times \hat{n} = \hat{x} \f$ on the boundary
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!---------------------------------------------------------------------------
program test_hcurl_sop
USE oft_base
! USE timer
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: mg_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_set_level, oft_hcurl_nlevels, &
  oft_hcurl, oft_bhcurl, oft_hcurl_eval_all
USE oft_hcurl_fields, ONLY: oft_hcurl_create
USE oft_hcurl_operators, ONLY: oft_hcurl_getkop, oft_hcurl_getwop, hcurl_zerob, &
  oft_hcurl_cinterp, oft_hcurl_bcurl
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_vector_inits, ONLY: uniform_field
USE diagnostic, ONLY: vec_energy
IMPLICIT NONE
INTEGER(i4) :: order,ierr,io_unit
NAMELIST/test_hcurl_options/order
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_hcurl_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!---
CALL oft_hcurl_setup(order)
!---Run tests
oft_env%pm=.FALSE.
CALL test_wop
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE: test_wop
!------------------------------------------------------------------------------
!> Solve the equation \f$ \nabla \times \nabla \times B = K \hat{I} \f$, where
!! \f$ K \f$ is the helicity matrix and \f$ \hat{I} \f$ is the identity vector
!! using H1(Curl) elements.
!------------------------------------------------------------------------------
SUBROUTINE test_wop
!---Create solver objects
CLASS(oft_solver), POINTER :: winv
!---Local variables
REAL(r8) :: uu
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: wop => NULL()
TYPE(oft_hcurl_cinterp) :: Bfield
TYPE(uniform_field) :: xfield
!---Set FE level
CALL oft_hcurl_set_level(oft_hcurl_nlevels)
!---Create solver fields
CALL oft_hcurl_create(u)
CALL oft_hcurl_create(v)
!---Get FE operators
CALL oft_hcurl_getwop(wop,'none')
!---Compute RHS
CALL oft_hcurl_bcurl(xfield,v)
!---Setup matrix solver
CALL create_cg_solver(winv)
winv%A=>wop
winv%its=-3
CALL create_diag_pre(winv%pre)
!---Solve
CALL u%set(0.d0)
CALL winv%apply(u,v)
!---Check results
Bfield%u=>u
CALL Bfield%setup
uu=vec_energy(Bfield,oft_hcurl%quad%order)
!---Report results
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='hcurl.results')
  WRITE(io_unit,*)winv%cits
  WRITE(io_unit,*)uu
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
DEALLOCATE(u,v)
!---Destroy matrix
CALL wop%delete
DEALLOCATE(wop)
!---Destory preconditioner
CALL winv%pre%delete
!---Destory solver
CALL winv%delete
END SUBROUTINE test_wop
end program test_hcurl_sop
