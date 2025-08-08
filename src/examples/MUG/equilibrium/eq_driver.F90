PROGRAM xmhd_eq
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
USE mhd_utils, ONLY: elec_charge, proton_mass, mu0
USE xmhd_2d
USE oft_io, ONLY: hdf5_field_get_sizes, hdf5_read, hdf5_field_exist
IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr, i
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
INTEGER(i4) :: ndims
INTEGER(i4) :: npoints
integer(i4), allocatable, dimension(:) :: dim_sizes
REAL(r8) :: n0 = 1.d0
REAL(r8) :: velx0 = 1.d0
REAL(r8) :: vely0 = 1.d0
REAL(r8) :: velz0 = 1.d0
REAL(r8) :: t0 = 1.d0
REAL(r8) :: psi0 = 0.d0
REAL(r8) :: by0 = 0.d0
REAL(r8) :: chi=1.E-12 !< Needs docs
REAL(r8) :: eta=1.d0 !< Needs docs
REAL(r8) :: nu=1.E-12 !< Needs docs
REAL(r8) :: gamma=1.d0
REAL(r8) :: D_diff=1.E-12
REAL(r8) :: bx0=0.d0
REAL(r8) :: bz0=0.d0
REAL(r8) :: B_0(3)=0.d0
REAL(r8) :: den_scale=1.d19
REAL(r8) :: dt = 1.d-3

REAL(r8), allocatable, dimension(:) :: psi_eq, bt_eq, p_eq, T_eq, F_eq

LOGICAL :: pm=.FALSE.
LOGICAL :: use_mfnk=.FALSE.
LOGICAL :: success
CHARACTER(LEN=16) :: filename = 'none' !< Name of input file for mesh, fix later for variable length

NAMELIST/xmhd_options/order,chi,eta,nu,gamma, D_diff, &
dt,nsteps,rst_freq,use_mfnk,pm, n0, psi0, velx0,vely0,velz0, t0, by0, den_scale, bx0, bz0
NAMELIST/equilibrium_options/filename

!------------------------------------------------------------------------------
! Initialize enviroment
!------------------------------------------------------------------------------
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,xmhd_options,IOSTAT=ierr)
CLOSE(io_unit)

OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,equilibrium_options,IOSTAT=ierr)
CLOSE(io_unit)
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)
CALL mhd_sim%setup(mg_mesh,order)
!---------------------------------------------------------------------------
! Read equilibrium from file
!---------------------------------------------------------------------------
CALL hdf5_field_get_sizes(TRIM(filename),"tokamaker/PSI",ndims,dim_sizes)
write(*,*) dim_sizes
npoints = dim_sizes(1)
ALLOCATE(psi_eq(npoints))
CALL hdf5_read(psi_eq,TRIM(filename),"tokamaker/PSI",success)

DEALLOCATE(dim_sizes)

CALL hdf5_field_get_sizes(TRIM(filename),"tokamaker/BT",ndims,dim_sizes)
npoints = dim_sizes(1)
ALLOCATE(bt_eq(npoints))
ALLOCATE(F_eq(npoints))
CALL hdf5_read(bt_eq,TRIM(filename),"tokamaker/BT",success)

CALL hdf5_field_get_sizes(TRIM(filename),"tokamaker/P",ndims,dim_sizes)
npoints = dim_sizes(1)
ALLOCATE(p_eq(npoints))
ALLOCATE(T_eq(npoints))
CALL hdf5_read(p_eq,TRIM(filename),"tokamaker/P",success)



!---------------------------------------------------------------------------
! Set intial conditions from equilibrium
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

mhd_sim%cyl_flag = .TRUE.

!---Project n initial condition onto scalar Lagrange basis
field_init%func=>const_init
field_init%mesh=>mesh
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(n0)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'n0')
mhd_sim%den_scale = den_scale
vec_vals = vec_vals / den_scale
CALL mhd_sim%u%restore_local(vec_vals,1)

!---Project v_x initial condition onto scalar Lagrange basis
field_init%func=>const_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(0.d0)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'vx0')
CALL mhd_sim%u%restore_local(vec_vals,2)

!---Project v_y initial condition onto scalar Lagrange basis
field_init%func=>const_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(0.d0)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'vy0')
CALL mhd_sim%u%restore_local(vec_vals,3)

!---Project v_z initial condition onto scalar Lagrange basis
field_init%func=>const_init
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(0.d0)
CALL u%get_local(vec_vals)
CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'vz0')
CALL mhd_sim%u%restore_local(vec_vals,4)

!---Project T initial condition onto scalar Lagrange basis
T_eq = p_eq/(2*n0*elec_charge) !have to divide in half because p from equilibrium file is for both electrons and ions
CALL mesh%save_vertex_scalar(T_eq,mhd_sim%xdmf_plot,'T0')
CALL mhd_sim%u%restore_local(T_eq,5)
! field_init%func=>const_init
! field_init%mesh=>mesh
! CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
! CALL u%set(0.d0)
! CALL minv%apply(u,v)
! CALL u%scale(2.5d3)
! CALL u%get_local(vec_vals)
! CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'T0')
! CALL mhd_sim%u%restore_local(vec_vals,5)

!---Project psi initial condition onto scalar Lagrange basis
CALL mesh%save_vertex_scalar(psi_eq,mhd_sim%xdmf_plot,'psi0')
CALL mhd_sim%u%restore_local(psi_eq,6)
! field_init%func=>const_init
! field_init%mesh=>mesh
! CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
! CALL u%set(0.d0)
! CALL minv%apply(u,v)
! CALL u%scale(0.d0)
! CALL u%get_local(vec_vals)
! CALL mesh%save_vertex_scalar(vec_vals,mhd_sim%xdmf_plot,'psi0')
! CALL mhd_sim%u%restore_local(vec_vals,6)

!---Project by initial condition onto scalar Lagrange basis
! CALL mesh%save_vertex_scalar(F_eq,mhd_sim%xdmf_plot,'by0')
! CALL mhd_sim%u%restore_local(F_eq,7)
field_init%func=>r_init
field_init%mesh=>mesh
CALL oft_blag_project(ML_oft_blagrange%current_level,field_init,v)
CALL u%set(0.d0)
CALL minv%apply(u,v)
CALL u%scale(1.d0)
CALL u%get_local(vec_vals)
F_eq = vec_vals*bt_eq
CALL mesh%save_vertex_scalar(F_eq,mhd_sim%xdmf_plot,'by0')
CALL mhd_sim%u%restore_local(F_eq,7)

!---Cleanup objects used for projection
CALL u%delete ! Destroy LHS vector
CALL v%delete ! Destroy RHS vector
CALL mop%delete ! Destroy mass matrix
DEALLOCATE(u,v,mop) ! Deallocate objects
CALL minv%pre%delete ! Destroy preconditioner
DEALLOCATE(minv%pre)
CALL minv%delete ! Destroy solver
DEALLOCATE(minv)

DEALLOCATE(bt_eq)
DEALLOCATE(F_eq)
DEALLOCATE(p_eq)
DEALLOCATE(psi_eq)
DEALLOCATE(T_eq)

!---------------------------------------------------------------------------
! Run simulation
!---------------------------------------------------------------------------
mhd_sim%chi=chi
mhd_sim%eta=eta
mhd_sim%nu=nu
mhd_sim%gamma=gamma
mhd_sim%D_diff=D_diff
mhd_sim%B_0 = B_0
mhd_sim%dt=dt
mhd_sim%nsteps=nsteps
mhd_sim%rst_freq=rst_freq
mhd_sim%mfnk=use_mfnk
oft_env%pm=pm
CALL mhd_sim%run_simulation()


!---Finalize enviroment
CALL oft_finalize
CONTAINS

SUBROUTINE const_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = 1.0
END SUBROUTINE const_init

SUBROUTINE gauss_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val =1.d0*EXP(-((pt(1)-0.4)**2 + pt(2)**2)/5.d-2)
END SUBROUTINE gauss_init

SUBROUTINE r_init(pt,val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8), INTENT(out) :: val
val = pt(1)
END SUBROUTINE r_init



END PROGRAM xmhd_eq