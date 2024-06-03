!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file gs_eq.F90
!
!> Classes and subroutines used for synthetic diagnostics.
!!
!! @authors Chris Hansen
!! @date April 2014
!! @ingroup doxy_oft_physics
!---------------------------------------------------------------------------
MODULE gs_eq
USE oft_base
USE oft_mesh_type, ONLY: mesh, bmesh_findcell
USE oft_mesh_local, ONLY: bmesh_local_init
USE oft_trimesh_type, ONLY: oft_trimesh
!---
USE fem_base, ONLY: oft_afem_type, oft_bfem_type
USE fem_utils, ONLY: fem_interp
USE oft_lag_basis, ONLY: oft_lag_setup_bmesh, oft_scalar_bfem
USE oft_blag_operators, ONLY: oft_lag_brinterp, oft_lag_bvrinterp
IMPLICIT NONE
#include "local.h"
PRIVATE
!---------------------------------------------------------------------------
! TYPE oft_gs_eq
!---------------------------------------------------------------------------
!> Grad-Shafranov equilibrium field
!---------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(fem_interp) :: oft_gs_eq
  INTEGER(i4) :: order = 1 !< Order of FE representation on input mesh
  INTEGER(i4) :: interp_mode = 1 !< Interpolation mode (1 -> 'B', 2 -> 'P')
  INTEGER(i4), POINTER, DIMENSION(:) :: pcell => NULL() !< Previous 2D cell
  INTEGER(i4), POINTER, DIMENSION(:) :: cell_tri => NULL() !< Current 2D cell
  REAL(r8) :: pmax = 0.d0 !< Maximum pressure
  REAL(r8), POINTER, DIMENSION(:) :: Pvals => NULL() !< Nodal pressure values
  REAL(r8), POINTER, DIMENSION(:,:) :: Bvals => NULL() !< Nodal B-Field values
  CHARACTER(LEN=OFT_HIST_SLEN) :: grid_file = 'none' !< File containing 2D grid
  CHARACTER(LEN=OFT_HIST_SLEN) :: field_file = 'none' !< File containing field data
  TYPE(oft_trimesh) :: mesh !< 2D triangular grid
  CLASS(oft_afem_type), POINTER :: lagrange => NULL() !< FE representation on 2D grid
  TYPE(oft_lag_brinterp) :: P_interp !< Pressure field interpolation object
  TYPE(oft_lag_bvrinterp) :: B_interp !< B-Field interpolation object
CONTAINS
  !> Load GS solution from file and setup interpolation
  PROCEDURE :: setup => gs_setup
  !> Reconstruct field at a given point in the 3D mesh
  PROCEDURE :: interp => gs_interp
  !> Delete reconstruction object, 2D mesh and FE representation
  PROCEDURE :: delete => gs_delete
END TYPE oft_gs_eq
CONTAINS
!---------------------------------------------------------------------------
! SUBROUTINE: gs_setup
!---------------------------------------------------------------------------
!> Setup interpolator for Grad-Shafranov fields
!!
!! Load equilibrium from file and constructs mesh and FE objects
!!
!! @note Should only be used via class \ref oft_gs_eq or children
!---------------------------------------------------------------------------
SUBROUTINE gs_setup(self)
CLASS(oft_gs_eq), INTENT(inout) :: self
INTEGER(i4) :: i,io_unit
REAL(r8) :: pmin
DEBUG_STACK_PUSH
!---Load GS grid
CALL self%mesh%load_from_file(TRIM(self%grid_file))
CALL bmesh_local_init(self%mesh)
!---Load GS field (order)
OPEN(NEWUNIT=io_unit,FILE=TRIM(self%field_file))
READ(io_unit,*)self%order
! ALLOCATE(self%lagrange)
CALL oft_lag_setup_bmesh(self%lagrange,self%mesh,self%order)
!---Load GS field (B,P)
ALLOCATE(self%Bvals(3,self%lagrange%ne),self%Pvals(self%lagrange%ne))
DO i=1,self%lagrange%ne
  READ(io_unit,*)self%Bvals(:,i),self%Pvals(i)
END DO
CLOSE(io_unit)
!---
pmin=MINVAL(self%Pvals)
self%pmax=MAXVAL(self%Pvals)
!
self%P_interp%vals=>self%Pvals
self%B_interp%vals=>self%Bvals
SELECT TYPE(this=>self%lagrange)
CLASS IS(oft_scalar_bfem)
  self%P_interp%lag_rep=>this
  self%B_interp%lag_rep=>this
END SELECT
!---
IF(.NOT.ASSOCIATED(self%pcell))THEN
  ALLOCATE(self%pcell(0:oft_env%nthreads-1),self%cell_tri(0:oft_env%nthreads-1))
  self%pcell=0
  self%cell_tri=0
END IF
!---
IF(oft_debug_print(1))THEN
  WRITE(*,*)'G-S Equilibrium loaded: ',self%field_file
  WRITE(*,*)'  Order:    ',self%order
  WRITE(*,*)'  Pressure: ',pmin,self%pmax
END IF
DEBUG_STACK_POP
END SUBROUTINE gs_setup
!---------------------------------------------------------------------------
! SUBROUTINE: gs_interp
!---------------------------------------------------------------------------
!> Reconstruct a Grad-Shafranov field
!!
!! @note Should only be used via class \ref oft_gs_eq
!!
!! @param[in] cell Cell for interpolation
!! @param[in] f Possition in cell in logical coord [4]
!! @param[in] gop Logical gradient vectors at f [3,4]
!! @param[out] val Reconstructed field at f [1]
!---------------------------------------------------------------------------
subroutine gs_interp(self,cell,f,gop,val)
class(oft_gs_eq), intent(inout) :: self
integer(i4), intent(in) :: cell
real(r8), intent(in) :: f(:)
real(r8), intent(in) :: gop(3,4)
real(r8), intent(out) :: val(:)
real(r8) :: rz(2),f_tri(3),goptmp(3,3),pt(3),b_axi(3),rhat(2),that(2)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%Bvals))CALL oft_abort('Setup has not been called!', &
'gs_interp',__FILE__)
!---Get (R,Z) position
pt=mesh%log2phys(cell,f)
rz(1)=SQRT(SUM(pt(1:2)**2))
rz(2)=pt(3)
rhat=pt(1:2)/rz(1)
that=(/-rhat(2),rhat(1)/)
!---Get logical position on TRI-MESH
IF(cell/=self%pcell(oft_tid))self%cell_tri(oft_tid)=0
self%pcell(oft_tid)=cell
CALL bmesh_findcell(self%mesh,self%cell_tri(oft_tid),rz,f_tri)
!---Evaluate requested field
SELECT CASE(self%interp_mode)
  CASE(1)
    CALL self%B_interp%interp(self%cell_tri(oft_tid),f_tri,goptmp,b_axi)
    val(1:2)=b_axi(1)*rhat + b_axi(2)*that
    val(3)=b_axi(3)
  CASE(2)
    CALL self%P_interp%interp(self%cell_tri(oft_tid),f_tri,goptmp,val)
  CASE DEFAULT
    CALL oft_abort('Invalid interpolation type','gs_interp',__FILE__)
END SELECT
DEBUG_STACK_POP
end subroutine gs_interp
!---------------------------------------------------------------------------
! SUBROUTINE: gs_delete
!---------------------------------------------------------------------------
!> Destroy temporary internal storage
!!
!! @note Should only be used via class \ref oft_lag_rinterp or children
!---------------------------------------------------------------------------
subroutine gs_delete(self)
class(oft_gs_eq), intent(inout) :: self
IF(.NOT.ASSOCIATED(self%lagrange))RETURN
DEBUG_STACK_PUSH
!---Destroy 2D mesh
CALL self%mesh%delete
!---Destroy 2D FE representation
CALL self%lagrange%delete
!---Destroy interpolation objects (destroys B,P)
CALL self%P_interp%delete
CALL self%B_interp%delete
!---Destroy thread locals
DEALLOCATE(self%pcell,self%cell_tri)
DEBUG_STACK_POP
end subroutine gs_delete
END MODULE gs_eq
