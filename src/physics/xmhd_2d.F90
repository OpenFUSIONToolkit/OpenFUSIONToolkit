!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file xmhd_2d.F90
!
!> Solve time-dependent MHD equations with lagrange basis
!---------------------------------------------------------------------------
MODULE xmhd_2d
USE oft_base
USE oft_io, ONLY: hdf5_read, hdf5_write, oft_file_exist, &
    hdf5_create_timestep, hdf5_field_exist, oft_bin_file
USE oft_quadrature
USE oft_mesh_type, ONLY: smesh, cell_is_curved
!
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_local_mat, oft_vector_ptr, &
    vector_extrapolate
USE oft_solver_utils, ONLY: create_solver_xml, create_diag_pre
USE oft_deriv_matrices, ONLY: oft_noop_matrix, oft_mf_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver
!
USE fem_composite, ONLY: oft_fem_comp_type
USE fem_utils, ONLY: fem_dirichlet_diag, fem_dirichlet_vec
USE oft_lag_basis, ONLY: oft_lag_setup, oft_blagrange, oft_blag_eval, oft_blag_geval
IMPLICIT NONE
#include "local.h"
PRIVATE

TYPE, extends(oft_noop_matrix) :: xmhd_2d_nlfun
  REAL(r8) :: dt = -1.d0 !< Time step
  REAL(r8) :: kappa_e !< Needs docs
  REAL(r8) :: kappa_i !< Needs docs
  REAL(r8) :: tau_eq !< Needs docs
  REAL(r8) :: diag_vals(2) = 0.d0 !< Needs docs
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: Ti_bc => NULL() !< Ti BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: Te_bc => NULL() !< Te BC flag
CONTAINS
  !> Apply the matrix
  PROCEDURE :: apply_real => nlfun_apply
END TYPE xmhd_2d_nlfun