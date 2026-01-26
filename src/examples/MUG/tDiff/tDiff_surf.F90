PROGRAM tDiff_circle
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
USE thermal_diffusion_2d
IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr
REAL(r8), POINTER :: vec_vals(:)
TYPE(oft_tdiff_sim) :: tDiff_sim
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_blag_zerob), TARGET :: blag_zerob
!---Mass matrix solver
TYPE(poss_scalar_bfield) :: field_init
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v
!---Runtime options
INTEGER(i4) :: order = 2
INTEGER(i4) :: nsteps = 100
INTEGER(i4) :: rst_freq = 10
REAL(r8) :: ti0 = 1.d0
REAL(r8) :: te0 = 4.d0
REAL(r8) :: kappa_i = 1.d0
REAL(r8) :: kappa_e = 1.d0
REAL(r8) :: tau_eq = -1.d0
REAL(r8) :: dt = 1.d-3
LOGICAL :: pm=.FALSE.
LOGICAL :: use_mfnk=.FALSE.
NAMELIST/tdiff_options/order,ti0,te0,kappa_i,kappa_e,tau_eq, &
dt,nsteps,rst_freq,use_mfnk,pm
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,tdiff_options,IOSTAT=ierr)
CLOSE(io_unit)
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)
CALL tDiff_sim%setup(mg_mesh,order)
blag_zerob%ML_lag_rep=>ML_oft_blagrange

!---------------------------------------------------------------------------
! Set intial conditions from analytic functions
!---------------------------------------------------------------------------
!---Generate mass matrix
NULLIFY(u,v,mop) ! Ensure the matrix is unallocated (pointer is NULL)
CALL oft_blag_getmop(ML_oft_blagrange%current_level,mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
!---Create fields for solver
CALL ML_oft_blagrange%vec_create(u)
CALL ML_oft_blagrange%vec_create(v)

!---Project Ti initial condition onto scalar Lagrange basis
field_init%func=>Ti_init
field_init%mesh=>mesh
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL blag_zerob%apply(u)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,tDiff_sim%xdmf_plot,'Ti0')
CALL tDiff_sim%u%restore_local(vec_vals,1)

!---Project Te initial condition onto scalar Lagrange basis
field_init%func=>Te_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL blag_zerob%apply(u)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,tDiff_sim%xdmf_plot,'Te0')
CALL tDiff_sim%u%restore_local(vec_vals,2)

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
tDiff_sim%kappa_i=kappa_i
tDiff_sim%kappa_e=kappa_e
tDiff_sim%tau_eq=tau_eq
tDiff_sim%dt=dt
tDiff_sim%nsteps=nsteps
tDiff_sim%rst_freq=rst_freq
tDiff_sim%mfnk=use_mfnk
oft_env%pm=pm
CALL tDiff_sim%run_simulation()
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!
SUBROUTINE Ti_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = Ti0*EXP(-SUM(pt**2)/2.d0)
END SUBROUTINE Ti_init
!
SUBROUTINE Te_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = Te0*(1.d0 + 1.d-1*SIN(pt(1)))
END SUBROUTINE Te_init
END PROGRAM tDiff_circle