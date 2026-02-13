!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!!MUG Example: Harris Sheet Reconnection    {#doc_mug2d_harris_ex}
!!============================
!!
!![TOC]
!!
!! This example demonstrates the use of the \ref xmhd2d "extended MHD" module in OFT. In
!! this example reconnection in a periodic-slab current sheet will be simulated following
!! the setup used to the [GEM Challange](https://doi.org/10.1029/1999JA900449).
!!
!!\section doc_mug2d_harris_ex_code_helper Code Walk Through
!!
!! The code consists of the main run program and subroutines for field initialization.
!!
!!\section doc_mug2d_harris_ex_code_helper Main Code
!!
!! As usual, we start with a few module imports for required functionality.
! START SOURCE
PROGRAM xmhd_circle
!---Runtime
USE oft_base
!---Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
!
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!
USE oft_blag_operators, ONLY: oft_blag_zerob, oft_blag_getmop, oft_blag_project
USE oft_scalar_inits, ONLY: poss_scalar_bfield
USE mhd_utils, ONLY: elec_charge, proton_mass
USE xmhd_2d
IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr
REAL(r8), POINTER :: vec_vals(:)
TYPE(oft_xmhd_2d_sim) :: mhd_sim
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_blag_zerob), TARGET :: blag_zerob ! setting boundary vals to zero
!---Mass matrix solver
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v
!!\subsection doc_mug2d_recon_ex_code_vars Local Variables
!! Next we define the local variables needed to initialize our case and
!! run the time-dependent solve and post-processing. Compared to previous
!! examples, we now have specialized initial condition and post-processing
!! objects.
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: nsteps = 100
INTEGER(i4) :: rst_freq = 10
REAL(r8) :: n0 = 1.d0
REAL(r8) :: velx0 = 1.d0
REAL(r8) :: vely0 = 1.d0
REAL(r8) :: velz0 = 1.d0
REAL(r8) :: t0 = 1.d0
REAL(r8) :: psi0 = 1.d0
REAL(r8) :: by0 = 1.d0
REAL(r8) :: B_0(3) = 0.d0
REAL(r8) :: chi=1.d0 !< Needs docs
REAL(r8) :: eta=1.d0 !< Needs docs
REAL(r8) :: nu=1.d0 !< Needs docs
REAL(r8) :: gamma=1.d0
REAL(r8) :: D_diff=1.d0
REAL(r8) :: den_scale=1.d19
REAL(r8) :: den_inf = 2.d-1
REAL(r8) :: lam_n = 2.d0 !< Wavelength
REAL(r8) :: lam_b = 2.d0 !< Wavelength
REAL(r8) :: L_x = 2.d0 !< Box length in x
REAL(r8) :: L_z = 2.d0 !< Box length in z
REAL(r8) :: delta = 1.d-1 !< Relative size of perturbation (<<1)
REAL(r8) :: dt = 1.d-3
LOGICAL :: pm=.FALSE.
LOGICAL :: linear=.FALSE.
LOGICAL :: cyl_flag=.FALSE.
LOGICAL :: use_mfnk=.FALSE.
NAMELIST/xmhd_options/order, chi, eta, nu, gamma, D_diff, &
dt, nsteps, rst_freq, use_mfnk, pm, n0, psi0, velx0, vely0, velz0, &
t0, by0, den_scale, delta, lam_n, lam_b, L_x, L_z, B_0, cyl_flag, linear
!!\subsection doc_mug_recon_ex_driver_setup OFT Initialization
!! See \ref doc_api_ex1 for a detailed description of calls in this section.
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,xmhd_options,IOSTAT=ierr)
CLOSE(io_unit)

!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)
CALL mhd_sim%setup(mg_mesh,order)
blag_zerob%ML_lag_rep=>ML_oft_blagrange

!---------------------------------------------------------------------------
! Set intial conditions from analytic functions
!---------------------------------------------------------------------------
!---Generate mass matrix
NULLIFY(u, v,mop,vec_vals) ! Ensure the matrix is unallocated (pointer is NULL)
CALL oft_blag_getmop(ML_oft_blagrange%current_level,mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
!---Create fields for solver
CALL ML_oft_blagrange%vec_create(u)
CALL ML_oft_blagrange%vec_create(v)

!---Set constant values

!---Project n initial condition onto scalar Lagrange basis
field_init%func=>dens_harris
field_init%mesh=>mesh
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
! CALL blag_zerob%apply(u)
CALL u%scale(n0)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'n0')
vec_vals = vec_vals / den_scale
CALL mhd_sim%u%restore_local(vec_vals,1)
CALL mhd_sim%u0%restore_local(vec_vals,1)

!---Project v_x initial condition onto scalar Lagrange basis
field_init%func=>const_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(velx0)
! CALL blag_zerob%apply(u)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'vx0')
CALL mhd_sim%u%restore_local(vec_vals,2)
IF (linear) CALL mhd_sim%u0%restore_local(vec_vals,2)

!---Project v_y initial condition onto scalar Lagrange basis
field_init%func=>const_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(vely0)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'vy0')
CALL mhd_sim%u%restore_local(vec_vals,3)
IF (linear) CALL mhd_sim%u0%restore_local(vec_vals,3)

!---Project v_z initial condition onto scalar Lagrange basis
field_init%func=>const_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(velz0)
! CALL blag_zerob%apply(u)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'vz0')
CALL mhd_sim%u%restore_local(vec_vals,4)
IF (linear) CALL mhd_sim%u0%restore_local(vec_vals,4)

!---Project T initial condition onto scalar Lagrange basis
field_init%func=>const_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(t0)
! CALL blag_zerob%apply(u)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'T0')
CALL mhd_sim%u%restore_local(vec_vals,5)
IF (linear) CALL mhd_sim%u0%restore_local(vec_vals,5)

!---Project psi equilibrium initial condition onto scalar Lagrange basis
IF (linear) THEN 
  field_init%func=>psi_harris_eq
  CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
  CALL u%set(0.d0)
  CALL minv%apply(u,v)
  ! CALL blag_zerob%apply(u)
  CALL u%scale(psi0)
  CALL u%get_local(vec_vals)
  CALL mhd_sim%u0%restore_local(vec_vals,6)
END IF

!---Project psi perturbed initial condition onto scalar Lagrange basis
field_init%func=>psi_harris_pert
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
! CALL blag_zerob%apply(u)
CALL u%scale(psi0)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'psi0')
CALL mhd_sim%u%restore_local(vec_vals,6)

!---Project by initial condition onto scalar Lagrange basis
field_init%func=>const_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(by0)
! CALL blag_zerob%apply(u)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'by0')
CALL mhd_sim%u%restore_local(vec_vals,7)
IF (linear) CALL mhd_sim%u0%restore_local(vec_vals,7)

!---Cleanup objects used for projection
CALL u%delete ! Destroy LHS vector
CALL v%delete ! Destroy RHS vector
CALL mop%delete ! Destroy mass matrix
DEALLOCATE(u,v,mop) ! Deallocate objects
CALL minv%pre%delete ! Destroy preconditioner
DEALLOCATE(minv%pre)
CALL minv%delete ! Destroy solver
DEALLOCATE(minv)

!---------------------------------------------------------------------------
! Run simulation
!---------------------------------------------------------------------------
mhd_sim%chi=chi
mhd_sim%eta=eta
mhd_sim%nu=nu
mhd_sim%gamma=gamma
mhd_sim%D_diff=D_diff
mhd_sim%dt=dt
mhd_sim%den_scale=den_scale
mhd_sim%nsteps=nsteps
mhd_sim%rst_freq=rst_freq
mhd_sim%mfnk=use_mfnk
mhd_sim%linear=linear
oft_env%pm=pm

IF (linear) THEN
  CALL mhd_sim%run_lin_simulation()
ELSE 
  CALL mhd_sim%run_simulation()
END IF

!---Finalize enviroment
CALL oft_finalize
CONTAINS
!! To set the initial conditions we define a set of functions 
!!
!! The non-uniform initial conditions for this case is given by
!! \f[ n_0 = (cosh(z / \lambda))^{-2} + n_{\infty} \f]
!! \f[ B_0 = tanh(z / \lambda) \hat{x} \f]
!!
!! The perturbation is then given by
!! \f[ \delta B = -\frac{\pi}{L_z} cos(2 \pi x / L_x) sin(\pi z / L_z) \hat{x} + \frac{2 \pi}{L_x} sin(2 \pi x / L_x) cos(\pi z / L_z) \hat{x} \f]
SUBROUTINE const_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = 1.0
END SUBROUTINE const_init

SUBROUTINE dens_harris(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val=(den_inf+1/((COSH(pt(2)/lam_n))**2))
END SUBROUTINE dens_harris

SUBROUTINE psi_harris_eq(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = -lam_b*LOG(COSH(pt(2)/lam_b))
END SUBROUTINE psi_harris_eq

SUBROUTINE psi_harris_pert(pt, val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = -lam_b*LOG(COSH(pt(2)/lam_b)) - delta*COS(2*pi*pt(1)/L_x)*COS(pi*pt(2)/L_z)
END SUBROUTINE psi_harris_pert
END PROGRAM xmhd_circle