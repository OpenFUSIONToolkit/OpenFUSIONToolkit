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
  hdf5_field_exist, oft_bin_file, xdmf_plot_file
USE oft_quadrature
USE oft_mesh_type, ONLY: oft_bmesh, cell_is_curved
USE multigrid, ONLY: multigrid_mesh
!
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_local_mat, oft_vector_ptr, &
  vector_extrapolate
USE oft_solver_utils, ONLY: create_solver_xml, create_diag_pre
USE oft_deriv_matrices, ONLY: oft_noop_matrix, oft_mf_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_fem_comp_type
USE fem_utils, ONLY: fem_dirichlet_diag, fem_dirichlet_vec, bfem_map_flag
USE oft_lag_basis, ONLY: oft_lag_setup,oft_scalar_bfem, oft_blag_eval, oft_blag_geval, oft_2D_lagrange_cast
USE oft_blag_operators, ONLY: oft_blag_vproject, oft_blag_getmop, oft_lag_bginterp
USE mhd_utils, ONLY: mu0, elec_charge, proton_mass
USE oft_gs, ONLY: gs_epsilon
IMPLICIT NONE
#include "local.h"
#if !defined(TDIFF_RST_LEN)
#define TDIFF_RST_LEN 5
#endif
PRIVATE

TYPE, extends(oft_noop_matrix) :: xmhd_2d_nlfun
  REAL(r8) :: dt = -1.d0 !< Time step
  REAL(r8) :: chi !< Needs docs
  REAL(r8) :: eta !< Needs docs
  REAL(r8) :: nu !< Needs docs
  REAL(r8) :: gamma
  REAL(r8) :: D_diff
  REAL(r8) :: k_boltz=elec_charge
  REAL(r8) :: m_i=proton_mass
  REAL(r8) :: den_scale = 1.d19 !< Needs docs
  REAL(r8) :: diag_vals(7) = 0.d0 !< Needs docs
  REAL (r8) :: B_0(3) = 0.d0

  LOGICAL :: cyl_flag

  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: T_bc => NULL() !< T BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velx_bc => NULL() !< vel BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: vely_bc => NULL() !< vel BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velz_bc => NULL() !< vel BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: by_bc => NULL() !< by BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: psi_bc => NULL() !< psi BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: n_bc => NULL() !< n BC flag
  CONTAINS
  !> Apply the matrix
  PROCEDURE :: apply_real => nlfun_apply
END TYPE xmhd_2d_nlfun
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
TYPE, public :: oft_xmhd_2d_sim
  LOGICAL :: mfnk = .FALSE. !< Use matrix free method?
  LOGICAL :: cyl_flag = .FALSE. !Use cylindrical coordinates version
  INTEGER(i4) :: nsteps = -1 !< Needs docs
  INTEGER(i4) :: rst_base = 0 !< Needs docs
  INTEGER(i4) :: rst_freq = 10 !< Needs docs
  REAL(r8) :: dt = -1.d0 !< Needs docs
  REAL(r8) :: t = 0.d0 !< Needs docs
  ! Edited to reflect new fields
  REAL(r8) :: chi = -1.d0 !< Needs docs
  REAL(r8) :: eta = -1.d0 !< Needs docs
  REAL(r8) :: nu = -1.d0 !< Needs docs
  REAL(r8) :: gamma = -1.d0 !< Needs docs
  REAL(r8) :: D_diff = -1.d0 !< Needs docs
  REAL(r8) :: k_boltz = elec_charge !< Needs docs
  REAL(r8) :: den_scale = 1.d19 !< Needs docs
  REAL(r8) :: m_i=proton_mass
  REAL(r8) :: lin_tol = 1.d-8 !< absolute tolerance for linear solver
  REAL(r8) :: nl_tol = 1.d-5 !< Needs docs
  REAL (r8) :: B_0(3) = 0.d0
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: n_bc => NULL() !< n BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velx_bc => NULL() !< vel BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: vely_bc => NULL() !< vel BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: velz_bc => NULL() !< vel BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: T_bc => NULL() !< T BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: psi_bc => NULL() !< psi BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: by_bc => NULL() !< by BC flag
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:,:) :: jacobian_block_mask => NULL() !< Matrix block mask
  TYPE(oft_fem_comp_type), POINTER :: fe_rep => NULL() !< Finite element representation for solution field
  TYPE(xdmf_plot_file) :: xdmf_plot
  CLASS(oft_vector), POINTER :: u => NULL() !< Needs docs
  CLASS(oft_matrix), POINTER :: jacobian => NULL() !< Needs docs
  TYPE(oft_mf_matrix), POINTER :: mf_mat => NULL() !< Matrix free operator
  TYPE(xmhd_2d_nlfun), POINTER :: nlfun => NULL() !< Needs docs
  TYPE(xml_node), POINTER :: xml_root => NULL() !< XML root element
  TYPE(xml_node), POINTER :: xml_pre_def => NULL() !< XML element for preconditioner definition
  contains
  !> Setup
  PROCEDURE :: setup => setup
  !> Run simulation
  PROCEDURE :: run_simulation => run_simulation
  !> Save restart file
  PROCEDURE :: rst_save => rst_save
  !> Load restart file
  PROCEDURE :: rst_load => rst_load
END TYPE oft_xmhd_2d_sim
TYPE(oft_xmhd_2d_sim), POINTER :: current_sim => NULL()
!
CLASS(multigrid_mesh), POINTER :: mg_mesh => NULL()
CLASS(oft_bmesh), POINTER, PUBLIC :: mesh => NULL()
TYPE(oft_ml_fem_type), TARGET, PUBLIC :: ML_oft_blagrange
CLASS(oft_scalar_bfem), POINTER :: oft_blagrange => NULL()
CONTAINS
!---------------------------------------------------------------------------
!> TODO: COMPLETE EDIT FOR NEW FIELDS
!---------------------------------------------------------------------------
subroutine run_simulation(self)
class(oft_xmhd_2d_sim), target, intent(inout) :: self
type(oft_nksolver) :: nksolver
!---Solver objects
class(oft_vector), pointer :: u,v,up
type(oft_native_gmres_solver), target :: solver
type(oft_timer) :: mytimer
!---History file fields
TYPE(oft_bin_file) :: hist_file
integer(i4) :: hist_i4(3)
real(4) :: hist_r4(9)
!---
CLASS(oft_matrix), POINTER :: lmop => NULL()
CLASS(oft_solver), POINTER :: lminv => NULL()
class(oft_vector), pointer :: ux,uy,uz,v_lag
type(oft_lag_bginterp) :: grad_psi
!---Extrapolation fields
integer(i4), parameter :: maxextrap=2
integer(i4) :: nextrap
real(r8), allocatable, dimension(:) :: extrapt
type(oft_vector_ptr), allocatable, dimension(:) :: extrap_fields
!---
character(LEN=TDIFF_RST_LEN) :: rst_char
integer(i4) :: i,j,io_stat,rst_tmp,npre
real(r8) :: n_avg, u_avg(3), T_avg, psi_avg, by_avg,elapsed_time
real(r8), pointer :: plot_vals(:),plot_vec(:,:)
current_sim=>self
!---------------------------------------------------------------------------
! Create solver fields
!---------------------------------------------------------------------------
call self%fe_rep%vec_create(u)
call self%fe_rep%vec_create(up)
call self%fe_rep%vec_create(v)
!---
call oft_blagrange%vec_create(grad_psi%u)
call oft_blagrange%vec_create(ux)
call oft_blagrange%vec_create(uy)
call oft_blagrange%vec_create(uz)
call oft_blagrange%vec_create(v_lag)
NULLIFY(lmop)
call oft_blag_getmop(oft_blagrange,lmop,'none')
CALL create_cg_solver(lminv)
lminv%A=>lmop
lminv%its=-2
CALL create_diag_pre(lminv%pre)

ALLOCATE(extrap_fields(maxextrap),extrapt(maxextrap))
DO i=1,maxextrap
  CALL self%fe_rep%vec_create(extrap_fields(i)%f)
  extrapt(i)=0.d0
END DO
nextrap=0

self%t=0.d0
CALL u%add(0.d0,1.d0,self%u)
!---Create initial conditions restart file
104 FORMAT (I TDIFF_RST_LEN.TDIFF_RST_LEN)
WRITE(rst_char,104)0
CALL self%rst_save(u, self%t, self%dt, 'xmhd2d_'//rst_char//'.rst', 'U')
NULLIFY(plot_vals)
ALLOCATE(plot_vec(3,v_lag%n))
CALL self%xdmf_plot%add_timestep(self%t)
CALL self%u%get_local(plot_vals,1)
plot_vals = plot_vals*self%den_scale
CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'n')
CALL self%u%get_local(plot_vals,2)
plot_vec(1,:)=plot_vals
CALL self%u%get_local(plot_vals,3)
plot_vec(3,:)=plot_vals
CALL self%u%get_local(plot_vals,4)
plot_vec(2,:)=plot_vals
CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'V')
CALL self%u%get_local(plot_vals,5)
CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'T')
CALL self%u%get_local(plot_vals,6)
CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'psi')
CALL grad_psi%u%restore_local(plot_vals)
!------------------------------------------------------------------------------
! Project magnetic field and plot
!------------------------------------------------------------------------------
CALL grad_psi%setup(oft_blagrange)
CALL oft_blag_vproject(oft_blagrange,grad_psi,ux,uy,uz)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,ux)
CALL ux%add(0.d0,1.d0,v_lag)
CALL v_lag%set(0.d0)
CALL lminv%apply(v_lag,uy)
CALL uy%add(0.d0,1.d0,v_lag)
!
CALL uy%get_local(plot_vals)
plot_vec(1,:)=-plot_vals
CALL self%u%get_local(plot_vals,7)
plot_vec(2,:)=plot_vals
CALL ux%get_local(plot_vals)
plot_vec(3,:)=plot_vals
CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'B')


!
ALLOCATE(self%nlfun)
self%nlfun%cyl_flag = self%cyl_flag
self%nlfun%chi=self%chi
self%nlfun%eta=self%eta
self%nlfun%nu=self%nu
self%nlfun%den_scale = self%den_scale
self%nlfun%D_diff=self%D_diff
self%nlfun%gamma=self%gamma
self%nlfun%B_0 = self%B_0
self%nlfun%n_bc=>self%n_bc
self%nlfun%velx_bc=>self%velx_bc
self%nlfun%vely_bc=>self%vely_bc
self%nlfun%velz_bc=>self%velz_bc
self%nlfun%T_bc=>self%T_bc
self%nlfun%psi_bc=>self%psi_bc
self%nlfun%by_bc=>self%by_bc
!---
CALL build_approx_jacobian(self,u) ! What is this u vector for?
!---------------------------------------------------------------------------
! Setup linear solver
!---------------------------------------------------------------------------
IF(self%mfnk)THEN
  ALLOCATE(self%mf_mat)
  CALL up%set(1.d0)
  CALL self%mf_mat%setup(u,self%nlfun,up)
  self%mf_mat%b0=1.d-5
  solver%A=>self%mf_mat
ELSE
  solver%A=>self%jacobian
END IF
solver%its=400
solver%atol=self%lin_tol
solver%itplot=1
solver%nrits=20
solver%pm=oft_env%pm
NULLIFY(solver%pre)
IF(ASSOCIATED(self%xml_pre_def))THEN
  CALL create_solver_xml(solver%pre,self%xml_pre_def)
ELSE
  CALL create_diag_pre(solver%pre)
END IF
solver%pre%A=>self%jacobian

!---------------------------------------------------------------------------
! Setup nonlinear solver
!---------------------------------------------------------------------------
nksolver%A=>self%nlfun
nksolver%J_inv=>solver
nksolver%its=30
nksolver%atol=self%nl_tol
nksolver%rtol=1.d-20 ! Disable relative tolerance
IF(self%mfnk)THEN
  nksolver%J_update=>mfnk_update
  nksolver%up_freq=1
ELSE
  nksolver%J_update=>update_jacobian
  nksolver%up_freq=4
END IF

!---------------------------------------------------------------------------
! Setup history file
!---------------------------------------------------------------------------
IF(oft_env%head_proc)THEN
  CALL hist_file%setup('oft_xmhd2d.hist', desc="History file for non-linear thermal diffusion run")
  CALL hist_file%add_field('ts',   'i4', desc="Time step index")
  CALL hist_file%add_field('lits', 'i4', desc="Linear iteration count")
  CALL hist_file%add_field('nlits','i4', desc="Non-linear iteration count")
  CALL hist_file%add_field('time', 'r4', desc="Simulation time [s]")
  CALL hist_file%add_field('ti',   'r4', desc="Average ion temperature [eV]")
  CALL hist_file%add_field('te',   'r4', desc="Average electron temperature [eV]")
  CALL hist_file%add_field('stime','r4', desc="Walltime [s]")
  CALL hist_file%write_header
  CALL hist_file%open ! Open history file
END IF

!---------------------------------------------------------------------------
! Begin time stepping
!---------------------------------------------------------------------------
npre=0
DO i=1,self%nsteps
  IF(oft_env%head_proc)CALL mytimer%tick()
  self%nlfun%dt=0.d0
  CALL self%nlfun%apply(u,v)
  n_avg=self%nlfun%diag_vals(1)
  u_avg = self%nlfun%diag_vals(2:4)
  T_avg = self%nlfun%diag_vals(5)
  psi_avg = self%nlfun%diag_vals(6)
  by_avg = self%nlfun%diag_vals(7)
  self%nlfun%dt=self%dt
  npre = npre + 1
  IF((.NOT.self%mfnk).OR.MOD(npre,4)==0)THEN
    CALL update_jacobian(u)
    CALL solver%pre%update(.TRUE.)
  END IF
  !
  DO j=maxextrap,2,-1
    CALL extrap_fields(j)%f%add(0.d0,1.d0,extrap_fields(j-1)%f)
    extrapt(j)=extrapt(j-1)
  END DO
  IF(nextrap<maxextrap)nextrap=nextrap+1
  CALL extrap_fields(1)%f%add(0.d0,1.d0,u)
  extrapt(1)=self%t
  IF(i>maxextrap)CALL vector_extrapolate(extrapt,extrap_fields,nextrap,self%t+self%dt,u)
  CALL nksolver%apply(u,v)
  IF(nksolver%cits<0)CALL oft_abort("Nonlinear solve failed","run_simulation",__FILE__)
  !---------------------------------------------------------------------------
  ! Write out initial solution progress
  !---------------------------------------------------------------------------
  IF(oft_env%head_proc)THEN
    elapsed_time=mytimer%tock()
    hist_i4=[self%rst_base+i-1,nksolver%lits,nksolver%nlits]
    hist_r4=REAL([self%t,n_avg, u_avg(1), u_avg(2), u_avg(3),T_avg,psi_avg,by_avg,elapsed_time],4)
103 FORMAT(' Timestep',I8,ES14.6,2X,I4,2X,I4,F12.3,ES12.2)
    WRITE(*,103)self%rst_base+i,self%t,nksolver%lits,nksolver%nlits,elapsed_time,self%dt
    IF(oft_debug_print(1))WRITE(*,*)
    CALL hist_file%write(data_i4=hist_i4, data_r4=hist_r4)
  END IF
  !---------------------------------------------------------------------------
  ! Update timestep and save solution
  !---------------------------------------------------------------------------
  self%t=self%t+self%dt
  CALL self%u%add(0.d0,1.d0,u)
  IF(MOD(i,self%rst_freq)==0)THEN
    IF(oft_env%head_proc)CALL mytimer%tick
    !---Create restart file
    WRITE(rst_char,104)self%rst_base+i
    READ(rst_char,104,IOSTAT=io_stat)rst_tmp
    IF((io_stat/=0).OR.(rst_tmp/=self%rst_base+i))CALL oft_abort("Step count exceeds format width", "run_simulation", __FILE__)
    CALL self%rst_save(u, self%t, self%dt, 'xhmd2d_'//rst_char//'.rst', 'U')
    IF(oft_env%head_proc)THEN
      elapsed_time=mytimer%tock()
      WRITE(*,'(2X,A,F12.3)')'I/O Time = ',elapsed_time
      CALL hist_file%flush
    END IF
    !---
    CALL self%xdmf_plot%add_timestep(self%t)
    CALL self%u%get_local(plot_vals,1)
    plot_vals = plot_vals*self%den_scale
    CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'n')
    CALL self%u%get_local(plot_vals,2)
    plot_vec(1,:)=plot_vals
    CALL self%u%get_local(plot_vals,3)
    plot_vec(3,:)=plot_vals
    CALL self%u%get_local(plot_vals,4)
    plot_vec(2,:)=plot_vals
    CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'V')
    CALL self%u%get_local(plot_vals,5)
    CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'T')
    CALL self%u%get_local(plot_vals,6)
    CALL mesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'psi')
    CALL grad_psi%u%restore_local(plot_vals)
    !------------------------------------------------------------------------------
    ! Project magnetic field and plot
    !------------------------------------------------------------------------------
    CALL grad_psi%setup(oft_blagrange)
    CALL oft_blag_vproject(oft_blagrange,grad_psi,ux,uy,uz)
    CALL v_lag%set(0.d0)
    CALL lminv%apply(v_lag,ux)
    CALL ux%add(0.d0,1.d0,v_lag)
    CALL v_lag%set(0.d0)
    CALL lminv%apply(v_lag,uy)
    CALL uy%add(0.d0,1.d0,v_lag)
    !
    CALL uy%get_local(plot_vals)
    plot_vec(1,:)=-plot_vals
    CALL self%u%get_local(plot_vals,7)
    plot_vec(2,:)=plot_vals
    CALL ux%get_local(plot_vals)
    plot_vec(3,:)=plot_vals
    plot_vec(1,:) = plot_vec(1,:)
    plot_vec(2,:) = plot_vec(2,:)
    plot_vec(3,:) = plot_vec(3,:)
    CALL mesh%save_vertex_vector(plot_vec,self%xdmf_plot,'B')
  END IF
  ! IF(nksolver%lits<4)THEN
  !   self%dt=self%dt*2.d0
  !   npre=-1
  ! ELSE IF(nksolver%lits>100)THEN
  !   npre=-1
  ! END IF
END DO
CALL hist_file%close()
CALL nksolver%delete()
CALL solver%delete()
CALL u%delete()
CALL up%delete()
CALL v%delete()
DEALLOCATE(u,up,v,plot_vals)
end subroutine run_simulation
!---------------------------------------------------------------------------
!> Compute the NL error function, where we are solving F(x) = 0
!!
!! b = F(a)
!---------------------------------------------------------------------------
subroutine nlfun_apply(self,a,b)
class(xmhd_2d_nlfun), intent(inout) :: self !< NL function object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
type(oft_quad_type), pointer :: quad
LOGICAL :: curved, cyl_flag
INTEGER(i4) :: i,m,jr, k,l
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs
REAL(r8) :: k_boltz = elec_charge
REAL(r8) :: m_i=proton_mass
REAL(r8) :: chi, eta, nu, D_diff, gamma, diag_vals(7), B_0(3), test
REAL(r8) :: n, vel(3), T, psi, by, dT(3),dn(3),dpsi(3),dby(3),&
         dvel(3,3),div_vel,jac_mat(3,4), jac_det,int_factor, btmp(3), tmp1(3), coords(3)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals,T_weights_loc,n_weights_loc, &
                     psi_weights_loc, by_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: vel_weights_loc, basis_grads,res_loc
REAL(r8), POINTER, DIMENSION(:) :: n_weights,T_weights,psi_weights,by_weights, T_res, &
                              n_res, psi_res, by_res, vtmp, velx_res, vely_res, velz_res
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights
quad=>oft_blagrange%quad
NULLIFY(n_weights, vel_weights, T_weights, psi_weights, by_weights, &
n_res, velx_res, vely_res, velz_res, T_res, psi_res, by_res)
!---Get weights from solution vector
ALLOCATE(vel_weights(3,oft_blagrange%ne))
CALL a%get_local(n_weights,1)
vtmp => vel_weights(1, :)
CALL a%get_local(vtmp ,2)
vtmp => vel_weights(2, :)
CALL a%get_local(vtmp ,3)
vtmp => vel_weights(3, :)
CALL a%get_local(vtmp, 4)
CALL a%get_local(T_weights,5)
CALL a%get_local(psi_weights,6)
CALL a%get_local(by_weights,7)
!---
chi = self%chi !< Needs docs
eta = self%eta !< Needs docs
nu = self%nu !< Needs docs
gamma = self%gamma !< Needs docs
D_diff = self%D_diff !< Needs docs
B_0 = self%B_0
cyl_flag = self%cyl_flag
! might need a version of this logic for eta or suchlike
! IF(self%tau_eq>0.d0)THEN
!   tau_eq_inv=1.d0/self%tau_eq
! ELSE
!   tau_eq_inv=1.d0
! END IF
!---Zero result and get storage array
CALL b%set(0.d0)
CALL b%get_local(n_res, 1)
CALL b%get_local(velx_res, 2)
CALL b%get_local(vely_res, 3)
CALL b%get_local(velz_res, 4)
CALL b%get_local(T_res, 5)
CALL b%get_local(psi_res, 6)
CALL b%get_local(by_res, 7)
diag_vals=0.d0

!$omp parallel private(m,jr,curved,coords,cell_dofs,basis_vals,basis_grads,T_weights_loc, &
!$omp n_weights_loc,psi_weights_loc, by_weights_loc,vel_weights_loc,res_loc,jac_mat, &
!$omp jac_det,int_factor,T,n,psi,by,vel,dT,dn,dpsi,dby,dvel,div_vel,btmp, tmp1, test) reduction(+:diag_vals)
!Edit for new fields
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(T_weights_loc(oft_blagrange%nce),n_weights_loc(oft_blagrange%nce),&
        psi_weights_loc(oft_blagrange%nce), by_weights_loc(oft_blagrange%nce),&
        vel_weights_loc(3, oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce),res_loc(oft_blagrange%nce,7))
!$omp do schedule(static)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  res_loc = 0.d0 ! Zero local (cell) contribution to function
  n_weights_loc = n_weights(cell_dofs)
  vel_weights_loc = vel_weights(:, cell_dofs)
  T_weights_loc = T_weights(cell_dofs)
  psi_weights_loc = psi_weights(cell_dofs)
  by_weights_loc = by_weights(cell_dofs)
  !---------------------------------------------------------------------------
  ! Quadrature Loop
  !---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,oft_blagrange%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(oft_blagrange,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    !--Extract spatial coordinates at current point
    coords = mesh%log2phys(i,quad%pts(:,m))
    !---Reconstruct values of solution fields
    n = 0.d0; dn = 0.d0; vel = 0.d0; dvel = 0.d0
    T = 0.d0; dT = 0.d0; psi = 0.d0; dpsi=0.d0
    by = 0.d0; dby = 0.d0
    basis_grads(3, :) = basis_grads(2,:)
    basis_grads(2,:) = 0.d0
    IF (cyl_flag) THEN
      int_factor = jac_det*quad%wts(m)*coords(1)
    ELSE
      int_factor = jac_det*quad%wts(m)
    END IF
    DO jr=1,oft_blagrange%nce
      n = n + n_weights_loc(jr)*basis_vals(jr)
      vel = vel + vel_weights_loc(:, jr)*basis_vals(jr)
      T = T + T_weights_loc(jr)*basis_vals(jr)
      psi = psi + psi_weights_loc(jr)*basis_vals(jr)
      by = by + by_weights_loc(jr)*basis_vals(jr)
      ! Note this is actually $(\nabla u)^T = jac(\vec(u))$
      ! Choosing this convention to make index contractions
      ! more consistent with Fortran convention
      dn = dn + n_weights_loc(jr)*basis_grads(:,jr)
      !Note: in cylindrical coordinates, dvel is not the actual gradient of velocity
      ! but the gradient of each component in a matrix
      dvel(:, 1) = dvel(:, 1) + vel_weights_loc(:, jr)*basis_grads(1, jr)
      dvel(:, 2) = 0.d0
      dvel(:, 3) = dvel(:, 3) + vel_weights_loc(:, jr)*basis_grads(3, jr)
      dT = dT + T_weights_loc(jr)*basis_grads(:,jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads(:,jr)
      dby = dby + by_weights_loc(jr)*basis_grads(:,jr)
    END DO
    n = n * self%den_scale
    dn = dn * self%den_scale
    IF (cyl_flag) THEN
      by = by/(coords(1)+gs_epsilon)
      dby(1) = (dby(1)-by)/ (coords(1)+gs_epsilon)
      dby(3) = dby(3)/(coords(1)+gs_epsilon)
    END IF
    IF (cyl_flag) THEN
      div_vel = dvel(1,1) +vel(1)/(coords(1)+gs_epsilon) + dvel(3,3)
    ELSE
      div_vel = dvel(1,1) + dvel(3,3)
    END IF
    diag_vals = diag_vals + [n, vel(1), vel(2), vel(3), T, psi, by]*int_factor
    IF (cyl_flag) THEN
      btmp = cross_product(dpsi/coords(1), [0.d0,1.d0,0.d0]) + by*[0.d0,1.d0,0.d0] + B_0
    ELSE
      btmp = cross_product(dpsi, [0.d0,1.d0,0.d0]) + by*[0.d0,1.d0,0.d0] + B_0
    END IF
    !---Compute local function contributions
    DO jr=1,oft_blagrange%nce
      !---Diffusion
      res_loc(jr,1) = res_loc(jr, 1) &
        + basis_vals(jr)*n*int_factor &
        +self%dt*basis_vals(jr)*DOT_PRODUCT(dn, vel)*int_factor &
        +self%dt*basis_vals(jr)*n*div_vel*int_factor &
        +self%dt*D_diff*DOT_PRODUCT(dn,basis_grads(:,jr))*int_factor
      !---Momentum
      res_loc(jr, 2:4) = res_loc(jr, 2:4) &
      + basis_vals(jr)*vel*int_factor&
        + self%dt*DOT_PRODUCT(btmp,basis_grads(:,jr))*btmp*int_factor/(mu0*m_i*n) & 
        - self%dt*DOT_PRODUCT(btmp,btmp)*basis_grads(:,jr)*int_factor/(2*mu0*m_i*n) & 
        - self%dt*basis_vals(jr)*DOT_PRODUCT(dn,btmp)*btmp*int_factor/(mu0*m_i*n**2) &
        + self%dt*basis_vals(jr)*DOT_PRODUCT(btmp,btmp)*dn*int_factor/(2*mu0*m_i*n**2) &
        + self%dt*basis_vals(jr)*2.d0*k_boltz*dT*int_factor/m_i &
        + self%dt*basis_vals(jr)*2.d0*k_boltz*T*dn*int_factor/(m_i*n)
      DO k=1,3
        res_loc(jr,k+1) = res_loc(jr, k+1) &
          + basis_vals(jr)*self%dt*DOT_PRODUCT(vel,dvel(k,:))*int_factor &
          + nu*self%dt*DOT_PRODUCT(basis_grads(:,jr),dvel(k,:))*int_factor/(m_i*n) &
          - nu*basis_vals(jr)*self%dt*DOT_PRODUCT(dn,dvel(k,:))*int_factor/(m_i*n**2)
      END DO 
      IF (cyl_flag) THEN
        res_loc(jr,2) = res_loc(jr,2) &
        + self%dt*basis_vals(jr)*(btmp(2)**2-btmp(1)**2-btmp(3)**2)*int_factor/(2.d0*mu0*m_i*n*(coords(1)+gs_epsilon)) &
        + self%dt*basis_vals(jr)*nu*vel(1)*int_factor/(m_i*n*(coords(1)+gs_epsilon)**2)  &
        - self%dt*basis_vals(jr)*vel(2)**2*int_factor/(coords(1)+gs_epsilon)! &
    
        res_loc(jr,3) = res_loc(jr,3) &
        - self%dt*basis_vals(jr)*btmp(1)*btmp(2)*int_factor /(mu0*m_i*n*(coords(1)+gs_epsilon)) & !!nonzero 
        + self%dt*basis_vals(jr)*nu*vel(2)*int_factor/(m_i*n*(coords(1)+gs_epsilon)**2) &
        + self%dt*basis_vals(jr)*vel(1)*vel(2)*int_factor/(coords(1)+gs_epsilon) 
      END IF
      !---Temperature
      res_loc(jr,5) = res_loc(jr, 5) &
        + basis_vals(jr)*T*int_factor/(gamma-1) &
        + self%dt*basis_vals(jr)*DOT_PRODUCT(vel, dT)*int_factor/(gamma-1) &
        + self%dt*basis_vals(jr)*T*div_vel*int_factor & 
        + self%dt*chi*DOT_PRODUCT(dT, basis_grads(:,jr))*int_factor &
        - self%dt*chi*basis_vals(jr)*DOT_PRODUCT(dn, dT)*int_factor/n 
      !---Induction
      tmp1 = cross_product(B_0,vel)
      res_loc(jr, 6) = res_loc(jr, 6) &
        + basis_vals(jr)*psi*int_factor &
        + basis_vals(jr)*self%dt*DOT_PRODUCT(vel, dpsi)*int_factor &
        + basis_vals(jr)*self%dt*tmp1(2)*int_factor &
        + self%dt*eta*DOT_PRODUCT(basis_grads(:,jr), dpsi)*int_factor/mu0
      IF (cyl_flag) THEN
        res_loc(jr,6) = res_loc(jr, 6) &
        + basis_vals(jr)*self%dt*eta*2.d0*dpsi(1)*int_factor/(mu0*(coords(1)+gs_epsilon))
      END IF
      tmp1 = cross_product(dpsi,dvel(2, :))
      IF (cyl_flag) THEN
        res_loc(jr, 7) = res_loc(jr, 7) &
        + (basis_vals(jr)*by*int_factor &
        - basis_vals(jr)*self%dt*tmp1(2)*int_factor/coords(1) &
        + basis_vals(jr)*self%dt*DOT_PRODUCT(vel, dby)*int_factor &
        + basis_vals(jr)*self%dt*by*(dvel(1,1) + dvel(3,3))*int_factor &
        + self%dt*eta*DOT_PRODUCT(basis_grads(:,jr), dby)*int_factor/mu0 &
        + basis_vals(jr)*self%dt*eta*by*int_factor/(mu0*(coords(1)+gs_epsilon)**2))*coords(1)
      ELSE
        res_loc(jr, 7) = res_loc(jr, 7) &
        + basis_vals(jr)*by*int_factor &
        - basis_vals(jr)*self%dt*tmp1(2)*int_factor &
        + basis_vals(jr)*self%dt*DOT_PRODUCT(vel, dby)*int_factor &
        + basis_vals(jr)*self%dt*by*div_vel*int_factor &
        + self%dt*eta*DOT_PRODUCT(basis_grads(:,jr), dby)*int_factor/mu0
      END IF 
    END DO
  END DO
    !---Add local values to full vector
  !write(*,*) res_loc(1,2)
  DO jr=1,oft_blagrange%nce
    !$omp atomic
    n_res(cell_dofs(jr)) = n_res(cell_dofs(jr)) + res_loc(jr,1)/self%den_scale
    !$omp atomic
    velx_res(cell_dofs(jr)) = velx_res(cell_dofs(jr)) + res_loc(jr,2)
    !$omp atomic
    vely_res(cell_dofs(jr)) = vely_res(cell_dofs(jr)) + res_loc(jr,3)
    !$omp atomic
    velz_res(cell_dofs(jr)) = velz_res(cell_dofs(jr)) + res_loc(jr,4)
    !$omp atomic
    T_res(cell_dofs(jr)) = T_res(cell_dofs(jr)) + res_loc(jr,5)
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr,6)
    !$omp atomic
    by_res(cell_dofs(jr)) = by_res(cell_dofs(jr)) + res_loc(jr,7)
  END DO
END DO

!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads,n_weights_loc,T_weights_loc,&
          vel_weights_loc, psi_weights_loc, by_weights_loc,cell_dofs,res_loc)
!$omp end parallel
IF(oft_debug_print(2))write(*,'(4X,A)')'Applying BCs'
CALL fem_dirichlet_vec(oft_blagrange,n_weights,n_res,self%n_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(1, :),velx_res,self%velx_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(2, :),vely_res,self%vely_bc)
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(3, :),velz_res,self%velz_bc)
CALL fem_dirichlet_vec(oft_blagrange,T_weights,T_res,self%T_bc)
CALL fem_dirichlet_vec(oft_blagrange,psi_weights,psi_res,self%psi_bc)
CALL fem_dirichlet_vec(oft_blagrange,by_weights,by_res,self%by_bc)
!---Put results into full vector
CALL b%restore_local(n_res,1,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(velx_res,2,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(vely_res,3,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(velz_res,4,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(T_res,5,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(psi_res,6,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(by_res,7,add=.TRUE.)
self%diag_vals=oft_mpi_sum(diag_vals,7)
!---Cleanup remaining storage
DEALLOCATE(n_res,velx_res,vely_res, velz_res, T_res, psi_res, by_res, &
        n_weights,vel_weights, T_weights, psi_weights, by_weights)
end subroutine nlfun_apply
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine build_approx_jacobian(self,a)
class(oft_xmhd_2d_sim), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Solution for computing jacobian
LOGICAL :: curved, cyl_flag
INTEGER(i4) :: i,m,jr,jc, k,l
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs
REAL(r8) :: k_boltz=elec_charge
REAL(r8) :: m_i = proton_mass
REAL(r8) :: chi, eta, nu, D_diff, gamma, B_0(3), diag_vals(7)
REAL(r8) :: n, vel(3), T, psi, by, dT(3),dn(3),dpsi(3),dby(3),&
dvel(3,3),div_vel,jac_mat(3,4), jac_det,int_factor, btmp(3), tmp2(3), tmp3(3), coords(3)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals,n_weights_loc,T_weights_loc,&
                                    psi_weights_loc, by_weights_loc, res_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: vel_weights_loc, basis_grads
REAL(r8), POINTER, DIMENSION(:) :: n_weights,T_weights, psi_weights, by_weights, vtmp
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights
type(oft_1d_int), allocatable, dimension(:) :: iloc
class(oft_vector), pointer :: tmp
type(oft_local_mat), allocatable, dimension(:,:) :: jac_loc
integer(KIND=omp_lock_kind), allocatable, dimension(:) :: tlocks
type(oft_quad_type), pointer :: quad
quad=>oft_blagrange%quad
CALL self%jacobian%zero
NULLIFY(n_weights,vel_weights, T_weights, &
         psi_weights, by_weights, vtmp)
!---Get weights from solution vector
CALL a%get_local(n_weights,1)
ALLOCATE(vel_weights(3,oft_blagrange%ne))
vtmp => vel_weights(1, :)
CALL a%get_local(vtmp ,2)
vtmp => vel_weights(2, :)
CALL a%get_local(vtmp ,3)
vtmp => vel_weights(3, :)
CALL a%get_local(vtmp, 4)
CALL a%get_local(T_weights,5)
CALL a%get_local(psi_weights,6)
CALL a%get_local(by_weights,7)
!---
chi = self%chi !< Needs docs
eta = self%eta !< Needs docs
nu = self%nu !< Needs docs
gamma = self%gamma !< Needs docs
D_diff = self%D_diff !< Needs docs
B_0 = self%B_0
cyl_flag = self%cyl_flag
!cyl_flag = .TRUE.
! IF(self%tau_eq>0.d0)THEN
!   tau_eq_inv=1.d0/self%tau_eq
! ELSE
!   tau_eq_inv=1.d0
! END IF
!--Setup thread locks
ALLOCATE(tlocks(self%fe_rep%nfields))
DO i=1,self%fe_rep%nfields
  call omp_init_lock(tlocks(i))
END DO
!$omp parallel private(m,jr,jc,curved,coords,cell_dofs,basis_vals,basis_grads,T_weights_loc, &
!$omp n_weights_loc,vel_weights_loc, psi_weights_loc, by_weights_loc, btmp, &
!$omp n, T, vel, by, psi,jac_loc,jac_mat,jac_det,int_factor,dn, dT, dvel, div_vel, dpsi, dby,iloc, &
!$omp tmp2, tmp3) 
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(n_weights_loc(oft_blagrange%nce),vel_weights_loc(3, oft_blagrange%nce),&
        T_weights_loc(oft_blagrange%nce), psi_weights_loc(oft_blagrange%nce),&
        by_weights_loc(oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce))
ALLOCATE(jac_loc(self%fe_rep%nfields,self%fe_rep%nfields))
ALLOCATE(iloc(self%fe_rep%nfields))
DO i=1,self%fe_rep%nfields
   iloc(i)%v=>cell_dofs
END DO
CALL self%fe_rep%mat_setup_local(jac_loc, self%jacobian_block_mask)
!$omp do schedule(static)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  CALL self%fe_rep%mat_zero_local(jac_loc) ! Zero local (cell) contribution to matrix
  n_weights_loc = n_weights(cell_dofs)
  vel_weights_loc = vel_weights(:, cell_dofs)
  T_weights_loc = T_weights(cell_dofs)
  psi_weights_loc = psi_weights(cell_dofs)
  by_weights_loc = by_weights(cell_dofs)
!---------------------------------------------------------------------------
! Quadrature Loop
!---------------------------------------------------------------------------
  DO m=1,quad%np
    if(curved.OR.(m==1))call mesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,oft_blagrange%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(oft_blagrange,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
    !--Extract spatial coordinates at current point
    coords = mesh%log2phys(i,quad%pts(:,m))
    basis_grads(3, :) = basis_grads(2,:)
    basis_grads(2,:) = 0.d0
    IF (cyl_flag) THEN
      int_factor = jac_det*quad%wts(m)*coords(1)
    ELSE
      int_factor = jac_det*quad%wts(m)
    END IF
    !---Reconstruct values of solution fields
    n = 0.d0; dn = 0.d0; vel = 0.d0; dvel = 0.d0
    T = 0.d0; dT = 0.d0; psi = 0.d0; dpsi=0.d0
    by = 0.d0; dby = 0.d0
    DO jr=1,oft_blagrange%nce
      n = n + n_weights_loc(jr)*basis_vals(jr)
      vel = vel + vel_weights_loc(:, jr)*basis_vals(jr)
      T = T + T_weights_loc(jr)*basis_vals(jr)
      psi = psi + psi_weights_loc(jr)*basis_vals(jr)
      by = by + by_weights_loc(jr)*basis_vals(jr)
      dn = dn + n_weights_loc(jr)*basis_grads(:,jr)
      ! Note this is actually $(\nabla u)^T = jac(\vec(u))$
      ! Choosing this convention to make index contractions
      ! more consistent with Fortran convention
      dvel(:, 1) = dvel(:, 1) + vel_weights_loc(:, jr)*basis_grads(1, jr)
      dvel(:, 2) = 0.d0
      dvel(:, 3) = dvel(:, 3) + vel_weights_loc(:, jr)*basis_grads(3, jr)
      dT = dT + T_weights_loc(jr)*basis_grads(:,jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads(:,jr)
      dby = dby + by_weights_loc(jr)*basis_grads(:,jr)
    END DO
    n = n * self%den_scale
    dn = dn * self%den_scale
    ! Redefine so by and dby corresponds to toroidal field and its gradient (rather than F)
    IF (cyl_flag) THEN
      by = by/(coords(1)+gs_epsilon)
      dby(1) = (dby(1)-by)/ (coords(1)+gs_epsilon)
      dby(3) = dby(3)/(coords(1)+gs_epsilon)
    END IF
    IF (cyl_flag) THEN
      div_vel = dvel(1,1) +vel(1)/(coords(1)+gs_epsilon) + dvel(3,3)
    ELSE
      div_vel = dvel(1,1) + dvel(3,3)
    END IF
    diag_vals = diag_vals + [n, vel(1), vel(2), vel(3), T, psi, by]*int_factor
    IF (cyl_flag) THEN
      btmp = cross_product(dpsi/coords(1), [0.d0,1.d0,0.d0]) + by*[0.d0,1.d0,0.d0] + B_0
    ELSE
      btmp = cross_product(dpsi, [0.d0,1.d0,0.d0]) + by*[0.d0,1.d0,0.d0] + B_0
    END IF
    !---Compute local matrix contributions
    DO jr=1,oft_blagrange%nce
      DO jc=1,oft_blagrange%nce
        ! Diffusion
        ! --n, n
        jac_loc(1, 1)%m(jr,jc) = jac_loc(1, 1)%m(jr, jc) &
        + basis_vals(jr)*basis_vals(jc)*int_factor &
        + self%dt*basis_vals(jr)*DOT_PRODUCT(basis_grads(:, jc), vel)*int_factor & !delta_n*div(u)
        + self%dt*basis_vals(jr)*basis_vals(jc)*div_vel*int_factor & ! u dot grad(delta_n)
        + self%dt*D_diff*DOT_PRODUCT(basis_grads(:, jr),basis_grads(:, jc))*int_factor
        ! --n, vel
        DO l=1,3
          jac_loc(1, l+1)%m(jr,jc) = jac_loc(1, l+1)%m(jr, jc) &
          + basis_vals(jr)*self%dt*n*basis_grads(l, jc)*int_factor & ! n*div(delta_u) = n*SUM(basis_grads)?
          + basis_vals(jr)*self%dt*basis_vals(jc)*dn(l)*int_factor ! u dot grad(n)
          IF(cyl_flag .AND. l==1)THEN
            jac_loc(1, l+1)%m(jr,jc) = jac_loc(1, l+1)%m(jr, jc) &
            + basis_vals(jr)*self%dt*n*basis_vals(jc)*int_factor/(coords(1)+gs_epsilon)
          END IF
        END DO
        !--Momentum 
        !--vel,n
        DO l=1,3
          jac_loc(l+1,1)%m(jr,jc) = jac_loc(l+1,1)%m(jr,jc) &
          + self%dt*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*dT(l)*int_factor/(m_i*n) &
          + self%dt*basis_vals(jr)*2.d0*k_boltz*T*basis_grads(l,jc)*int_factor/(m_i*n) &
          - self%dt*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*n*dT(l)*int_factor/(m_i*n**2.d0) &
          - self%dt*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*dn(l)*T*int_factor/(m_i*n**2.d0) &
          - self%dt*basis_vals(jc)*DOT_PRODUCT(basis_grads(:,jr), btmp)*btmp(l)*int_factor/(mu0*m_i*n**2.d0) &
          - self%dt*basis_vals(jr)*DOT_PRODUCT(basis_grads(:,jc), btmp)*btmp(l)*int_factor/(mu0*m_i*n**2.d0) &
          + self%dt*basis_vals(jr)*2.d0*basis_vals(jc)*DOT_PRODUCT(dn, btmp)*btmp(l)*int_factor/(mu0*m_i*n**3.d0) &
          + self%dt*basis_vals(jc)*basis_grads(l,jr)*DOT_PRODUCT(btmp, btmp)*int_factor/(2.d0*mu0*m_i*n**2.d0) &
          + self%dt*basis_vals(jr)*basis_grads(l,jc)*DOT_PRODUCT(btmp, btmp)*int_factor/(2.d0*mu0*m_i*n**2.d0) &
          - self%dt*basis_vals(jr)*basis_vals(jc)*dn(l)*DOT_PRODUCT(btmp, btmp)*int_factor/(mu0*m_i*n**3.d0)
        END DO 
        DO k=1,3
          jac_loc(k+1,1)%m(jr,jc) = jac_loc(k+1,1)%m(jr,jc) &
          - self%dt*nu*basis_vals(jc)*DOT_PRODUCT(basis_grads(:,jr), dvel(k, :))*int_factor/(m_i*n**2) & !-- not sure if indexing on dvel is right here
          - self%dt*nu*basis_vals(jr)*DOT_PRODUCT(basis_grads(:,jc), dvel(k, :))*int_factor/(m_i*n**2) &
          + self%dt*nu*basis_vals(jr)*2.d0*basis_vals(jc)*DOT_PRODUCT(dn, dvel(k, :))*int_factor/(m_i*n**3)
        END DO
        IF (cyl_flag) THEN
          jac_loc(2,1)%m(jr,jc) = jac_loc(2,1)%m(jr,jc) &
          - self%dt*basis_vals(jr)*basis_vals(jc)*(btmp(2)**2-btmp(1)**2-btmp(3)**2)*int_factor/(mu0*m_i*n**2*(coords(1)+gs_epsilon)) &
          - self%dt*basis_vals(jr)*basis_vals(jc)*nu*vel(1)*int_factor/(m_i*n**2*(coords(1)+gs_epsilon)**2)
          jac_loc(3,1)%m(jr,jc) = jac_loc(3,1)%m(jr,jc) &
          + self%dt*basis_vals(jr)*basis_vals(jc)*(btmp(1)*btmp(2))*int_factor/(mu0*m_i*n**2*(coords(1)+gs_epsilon)) &
          - self%dt*basis_vals(jr)*basis_vals(jc)*nu*vel(2)*int_factor/(m_i*n**2*(coords(1)+gs_epsilon)**2)
        END IF
        !--vel, vel
        DO k=1,3
          jac_loc(k+1,k+1)%m(jr,jc)= jac_loc(k+1,k+1)%m(jr,jc) &
            + basis_vals(jr)*basis_vals(jc)*int_factor &
            + self%dt*basis_vals(jr)*DOT_PRODUCT(vel, basis_grads(:,jc))*int_factor &
            + self%dt*nu*DOT_PRODUCT(basis_grads(:,jr), basis_grads(:,jc))*int_factor/(m_i*n) &
            - self%dt*nu*basis_vals(jr)*DOT_PRODUCT(dn, basis_grads(:,jc))*int_factor/(m_i*n**2)
          DO l=1,3
            jac_loc(k+1,l+1)%m(jr,jc)= jac_loc(k+1,l+1)%m(jr,jc) &
            + self%dt*basis_vals(jr)*basis_vals(jc)*dvel(k, l)*int_factor
          END DO
        END DO
        IF (cyl_flag) THEN
          jac_loc(2,2)%m(jr,jc)= jac_loc(2,2)%m(jr,jc) &
          + self%dt*nu*basis_vals(jr)*basis_vals(jc)*int_factor/(m_i*n*(coords(1)+gs_epsilon)**2)
          jac_loc(2,3)%m(jr,jc)= jac_loc(2,3)%m(jr,jc) &
          -self%dt*basis_vals(jr)*2.d0*vel(2)*basis_vals(jc)*int_factor/(coords(1)+gs_epsilon)
          jac_loc(3,2)%m(jr,jc)= jac_loc(3,2)%m(jr,jc) &
          + self%dt*basis_vals(jr)*basis_vals(jc)*vel(2)*int_factor/(coords(1)+gs_epsilon)
          jac_loc(3,3)%m(jr,jc)= jac_loc(3,3)%m(jr,jc) &
          + self%dt*basis_vals(jr)*basis_vals(jc)*vel(1)*int_factor/(coords(1)+gs_epsilon) &
          + self%dt*nu*basis_vals(jr)*basis_vals(jc)*int_factor/(m_i*n*(coords(1)+gs_epsilon))
        END IF
        ! --vel, T
        DO l=1,3
          jac_loc(l+1,5)%m(jr,jc) = jac_loc(l+1,5)%m(jr,jc) &
          + self%dt*basis_vals(jr)*2*k_boltz*dn(l)*basis_vals(jc)*int_factor/(m_i*n) &
          + self%dt*basis_vals(jr)*2*k_boltz*basis_grads(l,jc)*int_factor/(m_i) 
        END DO
        ! --vel, psi
        IF (cyl_flag) THEN
          tmp2 = cross_product(basis_grads(:,jc), [0.d0,1.d0/(coords(1)+gs_epsilon),0.d0]) ! this is 'dB'
        ELSE
          tmp2 = cross_product(basis_grads(:,jc), [0.d0,1.d0,0.d0])
        END IF
        DO l=1,3
          jac_loc(l+1,6)%m(jr,jc) = jac_loc(l+1,6)%m(jr,jc) &
          + self%dt*DOT_PRODUCT(basis_grads(:,jr), tmp2)*btmp(l)*int_factor/(m_i*n*mu0) &
          + self%dt*DOT_PRODUCT(basis_grads(:,jr), btmp)*tmp2(l)*int_factor/(m_i*n*mu0) &
          - self%dt*basis_grads(l,jr)*DOT_PRODUCT(tmp2, btmp)*int_factor/(m_i*n*mu0) &
          - self%dt*basis_vals(jr)*DOT_PRODUCT(dn, tmp2)*btmp(l)*int_factor/(m_i*n**2*mu0) &
          - self%dt*basis_vals(jr)*DOT_PRODUCT(dn, btmp)*tmp2(l)*int_factor/(m_i*n**2*mu0) &
          + self%dt*basis_vals(jr)*dn(l)*DOT_PRODUCT(tmp2, btmp)*int_factor/(m_i*n**2*mu0) 
        END DO
        IF (cyl_flag) THEN
          jac_loc(2,6)%m(jr,jc) = jac_loc(2,6)%m(jr,jc) &
          + self%dt*basis_vals(jr)*(btmp(2)*tmp2(2)-btmp(1)*tmp2(1)-btmp(3)*tmp2(3))*int_factor/(m_i*n*mu0*(coords(1)+gs_epsilon))
          jac_loc(3,6)%m(jr,jc) = jac_loc(3,6)%m(jr,jc) &
          - self%dt*basis_vals(jr)*(btmp(1)*tmp2(2)+btmp(2)*tmp2(1))*int_factor/(m_i*n*mu0*(coords(1)+gs_epsilon))
        END IF
        ! --vel, by
        IF (cyl_flag) THEN
          tmp2 = [0.d0,basis_vals(jc)/(coords(1)+gs_epsilon),0.d0]! this is 'B_y'
        ELSE
          tmp2 = [0.d0,basis_vals(jc),0.d0] 
        END IF
        IF (cyl_flag) THEN
          DO l=1,3
            jac_loc(l+1,7)%m(jr,jc) = jac_loc(l+1,7)%m(jr,jc) &
            + (self%dt*DOT_PRODUCT(basis_grads(:,jr),tmp2)*btmp(l)*int_factor/(m_i*n*mu0) &
            + self%dt*DOT_PRODUCT(basis_grads(:,jr), btmp)*tmp2(l)*int_factor/(m_i*n*mu0) &
            - self%dt*basis_grads(l,jr)*DOT_PRODUCT(tmp2, btmp)*int_factor/(m_i*n*mu0) &
            - self%dt*basis_vals(jc)*DOT_PRODUCT(dn, tmp2)*btmp(l)*int_factor/(m_i*n**2*mu0) &
            - self%dt*basis_vals(jc)*DOT_PRODUCT(dn, btmp)*tmp2(l)*int_factor/(m_i*n**2*mu0) &
            + self%dt*basis_vals(jc)*dn(l)*DOT_PRODUCT(tmp2, btmp)*int_factor/(m_i*n**2*mu0))/(coords(1)+gs_epsilon)
          END DO
          jac_loc(2,7)%m(jr,jc) = jac_loc(2,7)%m(jr,jc) &
          + self%dt*basis_vals(jr)*(btmp(2)*tmp2(2)-btmp(1)*tmp2(1)-btmp(3)*tmp2(3))*int_factor/(m_i*n*mu0*(coords(1)+gs_epsilon)**2)
          jac_loc(3,7)%m(jr,jc) = jac_loc(3,7)%m(jr,jc) &
          - self%dt*basis_vals(jr)*(btmp(1)*tmp2(2)+btmp(2)*tmp2(1))*int_factor/(m_i*n*mu0*(coords(1)+gs_epsilon)**2)
        ELSE
          DO l=1,3
            jac_loc(l+1,7)%m(jr,jc) = jac_loc(l+1,7)%m(jr,jc) &
            + self%dt*DOT_PRODUCT(basis_grads(:,jr),tmp2)*btmp(l)*int_factor/(m_i*n*mu0) &
            + self%dt*DOT_PRODUCT(basis_grads(:,jr), btmp)*tmp2(l)*int_factor/(m_i*n*mu0) &
            - self%dt*basis_grads(l,jr)*DOT_PRODUCT(tmp2, btmp)*int_factor/(m_i*n*mu0) &
            - self%dt*basis_vals(jc)*DOT_PRODUCT(dn, tmp2)*btmp(l)*int_factor/(m_i*n**2*mu0) &
            - self%dt*basis_vals(jc)*DOT_PRODUCT(dn, btmp)*tmp2(l)*int_factor/(m_i*n**2*mu0) &
            + self%dt*basis_vals(jc)*dn(l)*DOT_PRODUCT(tmp2, btmp)*int_factor/(m_i*n**2*mu0)
          END DO
        END IF
        ! Temperature
        ! T, n
        jac_loc(5, 1)%m(jr,jc) = jac_loc(5, 1)%m(jr, jc) &  
        - self%dt*chi*basis_vals(jr)*DOT_PRODUCT(basis_grads(:, jc), dT)*int_factor/n & ! grad(delta_n) !works w/o this line
        + self%dt*chi*basis_vals(jr)*basis_vals(jc)*DOT_PRODUCT(dn, dT)*int_factor/(n**2) ! delta_n !works w/o this line
        ! T, vel
        DO l=1,3
          jac_loc(5, l+1)%m(jr,jc) = jac_loc(5, l+1)%m(jr, jc) &  
          + basis_vals(jr) * self%dt*basis_vals(jc)*dT(l)*int_factor/(gamma-1) & ! delta_u dot Delta_T
          + basis_vals(jr) * self%dt*T*basis_grads(l, jc)*int_factor ! div(delta u) = SUM(basis_grads)?
          IF(cyl_flag .AND. l==1)THEN
            jac_loc(5, l+1)%m(jr,jc) = jac_loc(5, l+1)%m(jr, jc) &
            + basis_vals(jr)*self%dt*T*basis_vals(jc)*int_factor/(coords(1) + gs_epsilon) !!! this term is a problem
          END IF
        END DO
        ! T, T
        jac_loc(5, 5)%m(jr,jc) = jac_loc(5,5)%m(jr, jc) &  
        + basis_vals(jr) * basis_vals(jc)*int_factor/(gamma-1) &! delta_T
        + basis_vals(jr) * self%dt*DOT_PRODUCT(vel, basis_grads(:, jc))*int_factor/(gamma-1) & ! nabla(dT)
        + basis_vals(jr) * self%dt*basis_vals(jc)*div_vel*int_factor & ! dT != nabla(dT)
        + self%dt * chi * DOT_PRODUCT(basis_grads(:, jc),basis_grads(:, jr))*int_factor & ! dT_Chi
        - self%dt*basis_vals(jr)*chi*DOT_PRODUCT(dn, basis_grads(:, jc))*int_factor/n 
        ! Induction
        !--psi, vel
        DO l=1,3
          tmp2 = [0.d0, 0.d0,0.d0]
          jac_loc(6,l+1)%m(jr,jc) = jac_loc(6,l+1)%m(jr,jc) &
          + basis_vals(jr)*self%dt*basis_vals(jc)*dpsi(l)*int_factor
          tmp2(l) = 1.d0
          tmp3 = cross_product(B_0, tmp2)
          jac_loc(6,l+1)%m(jr,jc) = jac_loc(6,l+1)%m(jr,jc) &
          + basis_vals(jr)*self%dt*basis_vals(jc)*tmp3(2)*int_factor
        END DO
        ! --psi, psi
        jac_loc(6, 6)%m(jr,jc) = jac_loc(6, 6)%m(jr,jc) &
        + basis_vals(jr)*basis_vals(jc)*int_factor &
        + self%dt*basis_vals(jr)*DOT_PRODUCT(vel,basis_grads(:,jc))*int_factor &
        + self%dt*eta*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*int_factor/mu0
        IF (cyl_flag) THEN
          jac_loc(6, 6)%m(jr,jc) = jac_loc(6, 6)%m(jr,jc) &
          + self%dt*basis_vals(jr)*eta*2.d0*basis_grads(1,jc)*int_factor/(mu0*(coords(1) + gs_epsilon))
        END IF
        !--by, vel
        tmp2 = cross_product(dpsi,basis_grads(:,jc))
        IF (cyl_flag) THEN
          DO l=1,3
            jac_loc(7,l+1)%m(jr,jc) = jac_loc(7, l+1)%m(jr,jc) &
            + (basis_vals(jr)*self%dt*basis_vals(jc)*dby(l)*int_factor &
            + basis_vals(jr)*self%dt*by*basis_grads(l,jc)*int_factor)*coords(1)
            IF (l==2) THEN
              jac_loc(7,l+1)%m(jr,jc) = jac_loc(7, l+1)%m(jr,jc) &
              - basis_vals(jr)*self%dt*tmp2(l)*int_factor ! Factor of r cancels out on top and bottom
            END IF
          END DO
        ELSE
          DO l=1,3
            jac_loc(7,l+1)%m(jr,jc) = jac_loc(7, l+1)%m(jr,jc) &
            + basis_vals(jr)*self%dt*basis_vals(jc)*dby(l)*int_factor &
            + basis_vals(jr)*self%dt*by*basis_grads(l,jc)*int_factor
            IF (l==2) THEN
              jac_loc(7,l+1)%m(jr,jc) = jac_loc(7, l+1)%m(jr,jc) &
              - basis_vals(jr)*self%dt*tmp2(l)*int_factor
            END IF
          END DO
        END IF
        !-- by, psi
        tmp2 = cross_product(basis_grads(:,jc),dvel(2,:))
        jac_loc(7, 6)%m(jr,jc) = jac_loc(7, 6)%m(jr,jc) &
          - basis_vals(jr)*self%dt*tmp2(2)*int_factor! Factor of r cancels out on top and bottom     
        !-- by, by
        IF (cyl_flag) THEN
          jac_loc(7, 7)%m(jr,jc) = jac_loc(7, 7)%m(jr,jc) &
          + basis_vals(jr)*basis_vals(jc)*int_factor &
          + basis_vals(jr)*self%dt*DOT_PRODUCT(basis_grads(:,jc),vel)*int_factor & 
          + basis_vals(jr)*self%dt*basis_vals(jc)*(dvel(1,1) + dvel(3,3))*int_factor & 
          + self%dt*eta*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*int_factor/mu0 &
          + basis_vals(jr)*self%dt*eta*basis_vals(jc)*int_factor/((coords(1)+gs_epsilon)**2*mu0)
        ELSE
          jac_loc(7, 7)%m(jr,jc) = jac_loc(7, 7)%m(jr,jc) &
          + basis_vals(jr)*basis_vals(jc)*int_factor &
          + basis_vals(jr)*self%dt*DOT_PRODUCT(basis_grads(:,jc),vel)*int_factor & 
          + basis_vals(jr)*self%dt*basis_vals(jc)*div_vel*int_factor & 
          + self%dt*eta*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*int_factor/mu0 
        END IF
      END DO
    END DO
  END DO
  DO jr = 1,7
    jac_loc(1, jr)%m = jac_loc(1, jr)%m / self%den_scale
    jac_loc(jr, 1)%m = jac_loc(jr, 1)%m * self%den_scale
  END DO
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%n_bc(cell_dofs),1)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%velx_bc(cell_dofs), 2)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%vely_bc(cell_dofs), 3)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%velz_bc(cell_dofs), 4)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%T_bc(cell_dofs),5)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%psi_bc(cell_dofs),6)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%by_bc(cell_dofs), 7)
  CALL self%fe_rep%mat_add_local(self%jacobian,jac_loc,iloc,tlocks)
END DO
!---Cleanup thread-local storage
CALL self%fe_rep%mat_destroy_local(jac_loc)
DEALLOCATE(basis_vals,basis_grads,T_weights_loc,vel_weights_loc, &
          n_weights_loc, psi_weights_loc, by_weights_loc, cell_dofs,jac_loc,iloc)
!$omp end parallel
!--Destroy thread locks
DO i=1,self%fe_rep%nfields
  CALL omp_destroy_lock(tlocks(i))
END DO
DEALLOCATE(tlocks)
IF(oft_debug_print(2))write(*,'(4X,A)')'Setting BCs'
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%n_bc,1)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%velx_bc,2)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%vely_bc,3)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%velz_bc,4)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%T_bc,5)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%psi_bc,6)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%by_bc,7)
!
call self%fe_rep%vec_create(tmp)
call self%jacobian%assemble(tmp)
call tmp%delete
DEALLOCATE(tmp,n_weights,vel_weights, T_weights, &
          by_weights, psi_weights)
end subroutine build_approx_jacobian

!---------------------------------------------------------------------------
!> Update Jacobian matrices on all levels with new fields
!---------------------------------------------------------------------------
subroutine mfnk_update(uin)
class(oft_vector), target, intent(inout) :: uin !< Current field
IF(oft_debug_print(1))write(*,*)'Updating 2D MUG MF-Jacobian'
CALL current_sim%mf_mat%update(uin)
END SUBROUTINE mfnk_update
!---------------------------------------------------------------------------
!> Update Jacobian matrices on all levels with new solution
!---------------------------------------------------------------------------
subroutine update_jacobian(uin)
class(oft_vector), target, intent(inout) :: uin !< Current solution
IF(oft_debug_print(1))write(*,*)'Updating 2D MUG approximate Jacobian'
CALL build_approx_jacobian(current_sim,uin)
END SUBROUTINE update_jacobian
!---------------------------------------------------------------------------
!> Setup composite FE representation and ML environment
!---------------------------------------------------------------------------
subroutine setup(self,mg_mesh_in, order)
class(oft_xmhd_2d_sim), intent(inout) :: self
CLASS(multigrid_mesh), TARGET, intent(in) :: mg_mesh_in
integer(i4), intent(in) :: order
integer(i4) :: i,ierr,io_unit
LOGICAL, ALLOCATABLE :: vert_flag(:),edge_flag(:)
mg_mesh=>mg_mesh_in
mesh=>mg_mesh%smesh
IF(ASSOCIATED(self%fe_rep))CALL oft_abort("Setup can only be called once","setup",__FILE__)
IF(ASSOCIATED(oft_blagrange))CALL oft_abort("FE space already built","setup",__FILE__)

!---Look for XML defintion elements
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml))THEN
  CALL xml_get_element(oft_env%xml,"xmhd2d",self%xml_root,ierr)
  IF(ierr==0)THEN
    !---Look for pre node
    CALL xml_get_element(self%xml_root,"pre",self%xml_pre_def,ierr)
    IF(ierr/=0)NULLIFY(self%xml_pre_def)
  ELSE
    NULLIFY(self%xml_root)
  END IF
END IF
#endif

!---Setup FE representation
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Building lagrange FE space'
CALL oft_lag_setup(mg_mesh,order,ML_blag_obj=ML_oft_blagrange,minlev=-1)
IF(.NOT.oft_2D_lagrange_cast(oft_blagrange,ML_oft_blagrange%current_level))CALL oft_abort("Invalid lagrange FE object","setup",__FILE__)
CALL self%xdmf_plot%setup("xmhd_2d")
CALL mesh%setup_io(self%xdmf_plot,order)

!---Build composite FE definition for solution field
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Creating FE type'
ALLOCATE(self%fe_rep)
self%fe_rep%nfields=7
ALLOCATE(self%fe_rep%fields(self%fe_rep%nfields))
ALLOCATE(self%fe_rep%field_tags(self%fe_rep%nfields))
self%fe_rep%fields(1)%fe=>oft_blagrange
self%fe_rep%field_tags(1)='n'
self%fe_rep%fields(2)%fe=>oft_blagrange
self%fe_rep%field_tags(2)='velx'
self%fe_rep%fields(3)%fe=>oft_blagrange
self%fe_rep%field_tags(3)='vely'
self%fe_rep%fields(4)%fe=>oft_blagrange
self%fe_rep%field_tags(4)='velz'
self%fe_rep%fields(5)%fe=>oft_blagrange
self%fe_rep%field_tags(5)='T'
self%fe_rep%fields(6)%fe=>oft_blagrange
self%fe_rep%field_tags(6)='psi'
self%fe_rep%fields(7)%fe=>oft_blagrange
self%fe_rep%field_tags(7)='by'

!---Create solution vector
CALL self%fe_rep%vec_create(self%u)
! TODO: Boundary conditions not hard coded

!---Blocks to disable field evolutions via BCs
ALLOCATE(self%n_bc(oft_blagrange%ne)); self%n_bc=.FALSE.
! ALLOCATE(self%velx_bc(oft_blagrange%ne)); self%velx_bc=.TRUE.
! ALLOCATE(self%vely_bc(oft_blagrange%ne)); self%vely_bc=.TRUE.
! ALLOCATE(self%velz_bc(oft_blagrange%ne)); self%velz_bc=.TRUE.
ALLOCATE(self%T_bc(oft_blagrange%ne)); self%T_bc=.FALSE.
! ALLOCATE(self%psi_bc(oft_blagrange%ne)); self%psi_bc=.TRUE.
! ALLOCATE(self%by_bc(oft_blagrange%ne)); self%by_bc=.TRUE.
!ALLOCATE(self%n_bc(oft_blagrange%ne)); self%n_bc=.TRUE.
!ALLOCATE(self%velx_bc(oft_blagrange%ne)); self%velx_bc=.TRUE.
!ALLOCATE(self%vely_bc(oft_blagrange%ne)); self%vely_bc=.TRUE.
!ALLOCATE(self%velz_bc(oft_blagrange%ne)); self%velz_bc=.TRUE.
!ALLOCATE(self%T_bc(oft_blagrange%ne)); self%T_bc=.TRUE.
!ALLOCATE(self%psi_bc(oft_blagrange%ne)); self%psi_bc=.TRUE.
!ALLOCATE(self%by_bc(oft_blagrange%ne)); self%by_bc=.TRUE.

! !---Set boundary conditions for pipe flow in square mesh
! ALLOCATE(vert_flag(mesh%np),edge_flag(mesh%ne))
! !---Temperature is fixed on outlet only
! ALLOCATE(self%T_bc(oft_blagrange%ne))
! vert_flag=.FALSE.; edge_flag=.FALSE.
! DO i=1,mesh%nbe
!   IF(mesh%bes(i)==2)THEN
!     edge_flag(mesh%lbe(i))=.TRUE.
!     vert_flag(mesh%le(1,mesh%lbe(i)))=.TRUE.
!     vert_flag(mesh%le(2,mesh%lbe(i)))=.TRUE.
!   END IF
! END DO
! CALL bfem_map_flag(oft_blagrange,vert_flag,edge_flag,self%T_bc)
! !---Temperature is fixed on inlet/outlet
! ALLOCATE(self%n_bc(oft_blagrange%ne))
! vert_flag=.FALSE.; edge_flag=.FALSE.
! DO i=1,mesh%nbe
!   IF(mesh%bes(i)<3)THEN
!     edge_flag(mesh%lbe(i))=.TRUE.
!     vert_flag(mesh%le(1,mesh%lbe(i)))=.TRUE.
!     vert_flag(mesh%le(2,mesh%lbe(i)))=.TRUE.
!   END IF
! END DO
! CALL bfem_map_flag(oft_blagrange,vert_flag,edge_flag,self%n_bc)
! !---Velocity outlet needs work
! ALLOCATE(self%velx_bc(oft_blagrange%ne))
! vert_flag=.FALSE.; edge_flag=.FALSE.
! DO i=1,mesh%nbe
!   IF(mesh%bes(i)/=2)THEN
!     edge_flag(mesh%lbe(i))=.TRUE.
!     vert_flag(mesh%le(1,mesh%lbe(i)))=.TRUE.
!     vert_flag(mesh%le(2,mesh%lbe(i)))=.TRUE.
!   END IF
! END DO
! CALL bfem_map_flag(oft_blagrange,vert_flag,edge_flag,self%velx_bc)
! ! self%vely_bc=>self%velx_bc
! ! self%velz_bc=>self%velx_bc
!---Clean up flags
!DEALLOCATE(vert_flag,edge_flag)

!---Set any BCs that are not yet set
IF(.NOT.ASSOCIATED(self%n_bc))self%n_bc=>oft_blagrange%global%gbe
IF(.NOT.ASSOCIATED(self%velx_bc))self%velx_bc=>oft_blagrange%global%gbe
IF(.NOT.ASSOCIATED(self%vely_bc))self%vely_bc=>oft_blagrange%global%gbe
IF(.NOT.ASSOCIATED(self%velz_bc))self%velz_bc=>oft_blagrange%global%gbe
IF(.NOT.ASSOCIATED(self%T_bc))self%T_bc=>oft_blagrange%global%gbe
IF(.NOT.ASSOCIATED(self%psi_bc))self%psi_bc=>oft_blagrange%global%gbe
IF(.NOT.ASSOCIATED(self%by_bc))self%by_bc=>oft_blagrange%global%gbe

!---Create Jacobian matrix
ALLOCATE(self%jacobian_block_mask(self%fe_rep%nfields,self%fe_rep%nfields))
self%jacobian_block_mask=1
CALL self%fe_rep%mat_create(self%jacobian,self%jacobian_block_mask)
end subroutine setup
!---------------------------------------------------------------------------
!> Save xMHD solution state to a restart file
!---------------------------------------------------------------------------
subroutine rst_save(self,u,t,dt,filename,path)
class(oft_xmhd_2d_sim), intent(inout) :: self
class(oft_vector), pointer, intent(inout) :: u !< Solution to save
real(r8), intent(in) :: t !< Current solution time
real(r8), intent(in) :: dt !< Current timestep
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to store solution vector in file
DEBUG_STACK_PUSH
CALL self%fe_rep%vec_save(u,filename,path)
IF(oft_env%head_proc)THEN
  CALL hdf5_write(t,filename,'t')
  CALL hdf5_write(dt,filename,'dt')
END IF
DEBUG_STACK_POP
end subroutine rst_save
!---------------------------------------------------------------------------
!> Load xMHD solution state from a restart file
!---------------------------------------------------------------------------
subroutine rst_load(self,u,filename,path,t,dt)
class(oft_xmhd_2d_sim), intent(inout) :: self
class(oft_vector), pointer, intent(inout) :: u !< Solution to load
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to store solution vector in file
real(r8), optional, intent(out) :: t !< Time of loaded solution
real(r8), optional, intent(out) :: dt !< Timestep at loaded time
DEBUG_STACK_PUSH
CALL self%fe_rep%vec_load(u,filename,path)
IF(PRESENT(t))CALL hdf5_read(t,filename,'t')
IF(PRESENT(dt))CALL hdf5_read(dt,filename,'dt')
DEBUG_STACK_POP
end subroutine rst_load
END MODULE xmhd_2d