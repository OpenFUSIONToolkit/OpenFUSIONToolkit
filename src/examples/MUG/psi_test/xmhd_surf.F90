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
REAL(r8) :: chi=1.d0 !< Needs docs
REAL(r8) :: eta=1.d0 !< Needs docs
REAL(r8) :: nu=1.d0 !< Needs docs
REAL(r8) :: gamma=1.d0
REAL(r8) :: D_diff=1.d0
REAL(r8) :: dt = 1.d-3
LOGICAL :: pm=.FALSE.
LOGICAL :: use_mfnk=.FALSE.
NAMELIST/xmhd_options/order,chi,eta,nu,gamma, D_diff, &
dt,nsteps,rst_freq,use_mfnk,pm, n0, psi0, velx0,vely0,velz0, t0, by0
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
NULLIFY(u,v,mop,vec_vals) ! Ensure the matrix is unallocated (pointer is NULL)
CALL oft_blag_getmop(ML_oft_blagrange%current_level,mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
CALL create_diag_pre(minv%pre) ! Setup Preconditioner
!---Create fields for solver
CALL ML_oft_blagrange%vec_create(u)
CALL ML_oft_blagrange%vec_create(v)

!---Project n initial condition onto scalar Lagrange basis
field_init%func=>n_init
field_init%mesh=>mesh
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
! CALL blag_zerob%apply(u)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'n0')
CALL mhd_sim%u%restore_local(vec_vals,1)

!---Project v_x initial condition onto scalar Lagrange basis
field_init%func=>velx_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL blag_zerob%apply(u)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'vx0')
CALL mhd_sim%u%restore_local(vec_vals,2)

!---Project v_y initial condition onto scalar Lagrange basis
field_init%func=>vely_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL blag_zerob%apply(u)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'vy0')
CALL mhd_sim%u%restore_local(vec_vals,3)

!---Project v_z initial condition onto scalar Lagrange basis
field_init%func=>velz_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL blag_zerob%apply(u)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'vz0')
CALL mhd_sim%u%restore_local(vec_vals,4)

!---Project T initial condition onto scalar Lagrange basis
field_init%func=>T_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL blag_zerob%apply(u)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'T0')
CALL mhd_sim%u%restore_local(vec_vals,5)

!---Project psi initial condition onto scalar Lagrange basis
field_init%func=>psi_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL blag_zerob%apply(u)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'psi0')
CALL mhd_sim%u%restore_local(vec_vals,6)
!---Project by initial condition onto scalar Lagrange basis
field_init%func=>by_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL blag_zerob%apply(u)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'by0')
CALL mhd_sim%u%restore_local(vec_vals,7)

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
mhd_sim%nsteps=nsteps
mhd_sim%rst_freq=rst_freq
mhd_sim%mfnk=use_mfnk
oft_env%pm=pm
CALL mhd_sim%run_simulation()
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!
SUBROUTINE n_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = n0*(1.d0)
END SUBROUTINE n_init
!
SUBROUTINE velx_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = -velx0*(pt(2)/sqrt((pt(1)-0.5d0)**2+pt(2)**2))
END SUBROUTINE velx_init
!
SUBROUTINE vely_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = vely0*COS(pi*sqrt((pt(1)-0.5d0)**2+pt(2)**2)/2.d0)
END SUBROUTINE vely_init
!
SUBROUTINE velz_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = velz0*((pt(1)-0.5d0)/sqrt((pt(1)-0.5d0)**2+pt(2)**2))
END SUBROUTINE velz_init
!
SUBROUTINE t_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = t0*(1.d0)
END SUBROUTINE t_init  
!
SUBROUTINE psi_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = psi0*(1.d0)
END SUBROUTINE psi_init
!
SUBROUTINE by_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = by0*(1.d0)
END SUBROUTINE by_init
END PROGRAM xmhd_circle