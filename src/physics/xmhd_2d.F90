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
#if !defined(TDIFF_RST_LEN)
#define TDIFF_RST_LEN 5
#endif
PRIVATE

TYPE, extends(oft_noop_matrix) :: xmhd_2d_nlfun
  REAL(r8) :: dt = -1.d0 !< Time step
  REAL(r8) :: chi !< Needs docs
  REAL(r8) :: eta !< Needs docs
  REAL(r8) :: nu !< Needs docs
  REAL(r8) :: mu_0
  REAL(r8) :: gamma
  REAL(r8) :: D_diff
  REAL(r8) :: k_boltz
  REAL(r8) :: m_i
  REAL(r8) :: diag_vals(2) = 0.d0 !< Needs docs

  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: T_bc => NULL() !< T BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:,:) :: vel_bc => NULL() !< vel BC flag
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
  INTEGER(i4) :: nsteps = -1 !< Needs docs
  INTEGER(i4) :: rst_base = 0 !< Needs docs
  INTEGER(i4) :: rst_freq = 10 !< Needs docs
  REAL(r8) :: dt = -1.d0 !< Needs docs
  REAL(r8) :: t = 0.d0 !< Needs docs
  ! Edited to reflect new fields
  REAL(r8) :: chi = -1.d0 !< Needs docs
  REAL(r8) :: eta = -1.d0 !< Needs docs
  REAL(r8) :: nu = -1.d0 !< Needs docs
  REAL(r8) :: mu_0 = -1.d0 !< Needs docs
  REAL(r8) :: gamma = -1.d0 !< Needs docs
  REAL(r8) :: D_diff = -1.d0 !< Needs docs
  REAL(r8) :: k_boltz = -1.d0 !< Needs docs
  REAL(r8) :: m_i = -1.d0
  REAL(r8) :: lin_tol = 1.d-8 !< Needs docs
  REAL(r8) :: nl_tol = 1.d-5 !< Needs docs
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: n_bc => NULL() !< n BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:,:) :: vel_bc => NULL() !< vel BC flag
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
real(4) :: hist_r4(4)
!---Extrapolation fields
integer(i4), parameter :: maxextrap=2
integer(i4) :: nextrap
real(r8), allocatable, dimension(:) :: extrapt
type(oft_vector_ptr), allocatable, dimension(:) :: extrap_fields
!---
character(LEN=TDIFF_RST_LEN) :: rst_char
integer(i4) :: i,j,io_stat,rst_tmp
real(r8) :: n_avg, u_avg(3), T_avg, psi_avg, by_avg,elapsed_time 
real(r8), pointer :: plot_vals(:)
current_sim=>self
!---------------------------------------------------------------------------
! Create solver fields
!---------------------------------------------------------------------------
call self%fe_rep%vec_create(u)
call self%fe_rep%vec_create(up)
call self%fe_rep%vec_create(v)

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
CALL self%rst_save(u, self%t, self%dt, 'tDiff_'//rst_char//'.rst', 'U')
NULLIFY(plot_vals)
CALL self%xdmf_plot%add_timestep(self%t)
CALL self%u%get_local(plot_vals,1)
CALL smesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'Ti')
CALL self%u%get_local(plot_vals,2)
CALL smesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'Te')

!
ALLOCATE(self%nlfun)
self%nlfun%chi=self%chi
self%nlfun%eta=self%eta
self%nlfun%nu=self%nu
self%nlfun%D_diff=self%D_diff
self%nlfun%mu_0=self%mu_0
self%nlfun%gamma=self%gamma
self%nlfun%k_boltz=self%k_boltz
self%nlfun%m_i=self%m_i
self%nlfun%n_bc=>self%n_bc
self%nlfun%vel_bc=>self%vel_bc
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
solver%its=200
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
  CALL hist_file%setup('oft_tdiff.hist', desc="History file for non-linear thermal diffusion run")
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
DO i=1,self%nsteps
  IF(oft_env%head_proc)CALL mytimer%tick()
  self%nlfun%dt=0.d0
  CALL self%nlfun%apply(u,v)
  Ti_avg=self%nlfun%diag_vals(1)
  Te_avg=self%nlfun%diag_vals(2)
  self%nlfun%dt=self%dt
  IF((.NOT.self%mfnk).OR.MOD(i,2)==0)THEN
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
  !---------------------------------------------------------------------------
  ! Write out initial solution progress
  !---------------------------------------------------------------------------
  IF(oft_env%head_proc)THEN
    elapsed_time=mytimer%tock()
    hist_i4=[self%rst_base+i-1,nksolver%lits,nksolver%nlits]
    hist_r4=REAL([self%t,Ti_avg,Te_avg,elapsed_time],4)
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
    CALL self%rst_save(u, self%t, self%dt, 'tDiff_'//rst_char//'.rst', 'U')
    IF(oft_env%head_proc)THEN
      elapsed_time=mytimer%tock()
      WRITE(*,'(2X,A,F12.3)')'I/O Time = ',elapsed_time
      CALL hist_file%flush
    END IF
    !---
    CALL self%xdmf_plot%add_timestep(self%t)
    CALL self%u%get_local(plot_vals,1)
    CALL smesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'Ti')
    CALL self%u%get_local(plot_vals,2)
    CALL smesh%save_vertex_scalar(plot_vals,self%xdmf_plot,'Te')
  END IF
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
LOGICAL :: curved
INTEGER(i4) :: i,m,jr, k,l
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs
REAL(r8) :: chi, eta, nu, D_diff, gamma, mu_0, k_boltz, m_i, diag_vals(7)
REAL(r8) :: n, vel(3), T, psi, by, dT(3),dn(3),dpsi(3),dby(3),&
         dvel(3,3),div_vel,jac_mat(3,4), jac_det, btmp(3)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals,T_weights_loc,n_weights_loc, &
                     psi_weights_loc, by_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: vel_weights_loc, basis_grads,res_loc
REAL(r8), POINTER, DIMENSION(:) :: n_weights,T_weights,psi_weights,by_weights, T_res, &
                              n_res, psi_res, by_res, vtmp
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights, vel_res
quad=>oft_blagrange%quad
NULLIFY(n_weights, vel_weights, T_weights, psi_weights, by_weights, &
n_res, vel_res, T_res, psi_res, by_res)
!---Get weights from solution vector
CALL a%get_local(n_weights,1)
vtmp => vel_weights(1, :)
CALL a%get_local(vtmp ,2)
vtmp => vel_weights(2, :)
CALL a%get_local(vtmp ,3)
vtmp => vel_weights(3, :)
CALL a%get_local(vtmp, 4)
CALL a%get_local(psi_weights,5)
CALL a%get_local(T_weights,6)
CALL a%get_local(by_weights,7)
!---
chi = self%chi !< Needs docs
eta = self%eta !< Needs docs
nu = self%nu !< Needs docs
mu_0 = self%mu_0 !< Needs docs
gamma = self%gamma !< Needs docs
D_diff = self%D_diff !< Needs docs
k_boltz = self%k_boltz !< Needs docs
! might need a version of this logic for eta or suchlike
! IF(self%tau_eq>0.d0)THEN
!   tau_eq_inv=1.d0/self%tau_eq
! ELSE
!   tau_eq_inv=1.d0
! END IF
!---Zero result and get storage array
CALL b%set(0.d0)
CALL b%get_local(n_res, 1)
vtmp => vel_res(1, :)
CALL b%get_local(vtmp, 2)
vtmp => vel_res(2, :)
CALL b%get_local(vtmp, 3)
vtmp => vel_res(3, :)
CALL b%get_local(vtmp, 4)
CALL b%get_local(T_res, 5)
CALL b%get_local(psi_res, 6)
CALL b%get_local(by_res, 7)
diag_vals=0.d0

!$omp parallel private(m,jr,curved,cell_dofs,basis_vals,basis_grads,T_weights_loc, &
!$omp n_weights_loc,psi_weights_loc, by_weights_loc,vel_weights_loc,res_loc,jac_mat, &
!$omp jac_det,T,n,psi,by,vel,dT,dn,dpsi,dby,dvel, btmp) reduction(+:diag_vals)
!Edit for new fields
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(T_weights_loc(oft_blagrange%nce),n_weights_loc(oft_blagrange%nce),&
        psi_weights_loc(oft_blagrange%nce), by_weights_loc(oft_blagrange%nce),&
        vel_weights_loc(3, oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce),res_loc(oft_blagrange%nce,7))
!$omp do schedule(static)
DO i=1,smesh%nc
  curved=cell_is_curved(smesh,i) ! Straight cell test
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
    if(curved.OR.(m==1))call smesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,oft_blagrange%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(oft_blagrange,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
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
      dn = dn + n_weights_loc(jr)*basis_vals(jr)
      ! Note this is actually $(\nabla u)^T = jac(\vec(u))$
      ! Choosing this convention to make index contractions
      ! more consistent with Fortran convention
      dvel(:, 1) = dvel(:, 1) + vel_weights_loc(:, jr)*basis_grads(1, jr)
      dvel(:, 2) = 0.d0
      dvel(:, 3) = dvel(:, 2) + vel_weights_loc(:, jr)*basis_grads(3, jr)
      dT = dT + T_weights_loc(jr)*basis_grads(:,jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads(:,jr)
      dby = dby + by_weights_loc(jr)*basis_grads(:,jr)
    END DO
    div_vel = dvel(1,1) + dvel(2,2) + dvel(3,3)
    diag_vals = diag_vals + [T_i,T_e]*jac_det*quad%wts(m) !TODO: update this line for new fields
    btmp = cross_product(dpsi, [0,1,0]) + by*[0,1,0]
    !---Compute local function contributions
    DO jr=1,oft_blagrange%nce
      !---Diffusion
      res_loc(jr,1) = res_loc(jr, 1) &
        + basis_vals(jr)*n*jac_det*quad%wts(m) &
        +self%dt*basis_vals(jr)*DOT_PRODUCT(dn, vel)*jac_det*quad%wts(m) &
        +self%dt*basis_vals(jr)*n*div_vel*jac_det*quad%wts(m) &
        +self%dt*D_diff*DOT_PRODUCT(dn,basis_grads(:,jr))*jac_det*quad%wts(m)
      !---Momentum
      res_loc(jr, 2:4) = res_loc(jr, 2:4) &
        + basis_vals(jr)*vel*jac_det*quad%wts(m) &
        + self%dt*DOT_PRODUCT(btmp,basis_grads(:,jr))*btmp*jac_det*quad%wts(m)/(mu_0*m_i*n) &
        - self%dt*DOT_PRODUCT(btmp,btmp)*basis_grads(:,jr)*jac_det*quad%wts(m)/(2*mu_0*m_i*n) &
        - self%dt*basis_vals(jr)*DOT_PRODUCT(dn,btmp)*btmp*jac_det*quad%wts(m)/(mu_0*m_i*n**2) &
        + self%dt*basis_vals(jr)*DOT_PRODUCT(btmp,btmp)*dn*jac_det*quad%wts(m)/(2*mu_0*m_i*n**2) &
        + 2*self%dt*k_boltz*n*basis_vals(jr)*dT*jac_det*quad%wts(m)&
        + 2*self%dt*k_boltz*T*basis_vals(jr)*dn*jac_det*quad%wts(m)
      DO k=1,3
        res_loc(jr,2:4) = res_loc(jr, 2:4) &
          + basis_vals(jr)*self%dt*vel(k)*dvel(:,k)*jac_det*quad%wts(m) &
          + nu*self%dt*basis_grads(k,jr)*dvel(:,k)*jac_det*quad%wts(m)/(m_i*n) &
          - basis_vals(jr)*self%dt*dn(k)*dvel(:,k)*jac_det*quad%wts(m)/(m_i*n**2)
      END DO
      !---Temperature
      res_loc(jr,5) = res_loc(jr, 5) &
        + basis_vals(jr)*T*jac_det*quad%wts(m)/(gamma+1) &
        + self%dt*basis_vals(jr)*DOT_PRODUCT(vel, dT)*jac_det*quad%wts(m)/(gamma+1) &
        - self%dt*k_boltz*basis_vals(jr)*T*div_vel &
        + self%dt*chi*DOT_PRODUCT(dT, basis_grads(:,jr)) &
        - self%dt*chi*basis_vals(jr)*DOT_PRODUCT(dn, dT)/n
      !---Induction
      res_loc(jr, 6) = res_loc(jr, 6) &
        + basis_vals(jr)*psi*jac_det*quad%wts(m) &
        + basis_vals(jr)*self%dt*DOT_PRODUCT(vel, dpsi)*jac_det*quad%wts(m) &
        + self%dt*eta*DOT_PRODUCT(basis_grads(:,jr), dpsi)*jac_det*quad%wts(m)/mu_0
      res_loc(jr, 7) = res_loc(jr, 7) &
        + basis_vals(jr)*by*jac_det*quad%wts(m) &
        - basis_vals(jr)*self%dt*cross_product(dpsi,dvel(2, :))*jac_det*quad%wts(m) &
        + basis_vals(jr)*self%dt*(DOT_PRODUCT(dby,vel) + by*div_vel)*jac_det*quad%wts(m) &
        + basis_vals(jr)*self%dt*nu*DOT_PRODUCT(basis_grads(:,jr), dby)*jac_det*quad%wts(m)/mu_0
    END DO
  END DO
    !---Add local values to full vector
  DO jr=1,oft_blagrange%nce
    !$omp atomic
    n_res(cell_dofs(jr)) = n_res(cell_dofs(jr)) + res_loc(jr,1)
    !$omp atomic
    vel_res(cel_dofs(jr)) = vel_res(cell_dofs(jr)) + res_loc(jr,2:4)
    !$omp atomic
    T_res(cell_dofs(jr)) = T_res(cell_dofs(jr)) + res_loc(jr,5)
    !$omp atomic
    psi_res(cell_dofs(jr)) = psi_res(cell_dofs(jr)) + res_loc(jr,6)
    !$omp atomic
    by_res(cell_dofs(jr)) = by_res(cell_dofs(jr)) + res_loc(jr,7)
  END DO
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads,n_weights_loc,T_weights_loc,btmp,&
          vel_weights_loc, psi_weights_loc, by_weights_loc,cell_dofs,res_loc)
!$omp end parallel
IF(oft_debug_print(2))write(*,'(4X,A)')'Applying BCs'
CALL fem_dirichlet_vec(oft_blagrange,n_weights,n_res,self%n_bc)
!TODO: check validity of this usage
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(1, :),vel_res(1, :),self%vel_bc(1,:))
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(2, :),vel_res(2, :),self%vel_bc(2,:))
CALL fem_dirichlet_vec(oft_blagrange,vel_weights(3, :),vel_res(3, :),self%vel_bc(3,:))
CALL fem_dirichlet_vec(oft_blagrange,T_weights,T_res,self%T_bc)
CALL fem_dirichlet_vec(oft_blagrange,psi_weights,psi_res,self%psi_bc)
CALL fem_dirichlet_vec(oft_blagrange,by_weights,by_res,self%by_bc)
!---Put results into full vector
CALL b%restore_local(n_res,1,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(vel_res(1,:),2,add=.TRUE.)
CALL b%restore_local(vel_res(2,:),3,add=.TRUE.)
CALL b%restore_local(vel_res(3,:),4,add=.TRUE.)
CALL b%restore_local(T_res,5,add=.TRUE.)
CALL b%restore_local(psi_res,6,add=.TRUE.)
CALL b%restore_local(by_res,7,add=.TRUE.)
!TODO: adjust diagonal vals
self%diag_vals=oft_mpi_sum(diag_vals,2)
!---Cleanup remaining storage
DEALLOCATE(n_res,vel_res, T_res, psi_res, by_res, &
        n_weights,vel_weights, T_weights, psi_weights, by_weights)
end subroutine nlfun_apply
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine build_approx_jacobian(self,a)
class(oft_xmhd_2d_sim), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Solution for computing jacobian
LOGICAL :: curved
INTEGER(i4) :: i,m,jr,jc, k
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs
REAL(r8) :: chi, eta, nu, D_diff, gamma, mu_0, k_boltz, m_i
REAL(r8) :: n, vel(3), T, psi, by, dT(3),dn(3),dpsi(3),dby(3),&
dvel(3,3),div_vel,jac_mat(3,4), jac_det, vtmp(3), btmp(3)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals,n_weights_loc,T_weights_loc,&
                                    psi_weights_loc, by_weights_loc, res_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: vel_weights_loc, basis_grads
REAL(r8), POINTER, DIMENSION(:) :: n_weights,T_weights, psi_weights, by_weights
REAL(r8), POINTER, DIMENSION(:,:) :: vel_weights
type(oft_1d_int), allocatable, dimension(:) :: iloc
class(oft_vector), pointer :: tmp
type(oft_local_mat), allocatable, dimension(:,:) :: jac_loc
integer(KIND=omp_lock_kind), allocatable, dimension(:) :: tlocks
type(oft_quad_type), pointer :: quad
quad=>oft_blagrange%quad
CALL self%jacobian%zero
NULLIFY(n_weights,vel_weights, T_weights, &
         psi_weights, by_weights)
!---Get weights from solution vector
CALL a%get_local(n_weights,1)
vtmp => vel_weights(1, :)
CALL a%get_local(vtmp ,2)
vtmp => vel_weights(2, :)
CALL a%get_local(vtmp ,3)
vtmp => vel_weights(3, :)
CALL a%get_local(vtmp, 4)
CALL a%get_local(psi_weights,5)
CALL a%get_local(T_weights,6)
CALL a%get_local(by_weights,7)
!---
chi = self%chi !< Needs docs
eta = self%eta !< Needs docs
nu = self%nu !< Needs docs
mu_0 = self%mu_0 !< Needs docs
gamma = self%gamma !< Needs docs
D_diff = self%D_diff !< Needs docs
k_boltz = self%k_boltz !< Needs docs
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
!$omp parallel private(m,jr,jc,curved,cell_dofs,basis_vals,basis_grads,n_weights_loc&
!$omp vel_weights_loc, T_weights_loc, psi_weights_loc, by_weights_loc, btmp, &
!$omp n, T, vel, by, psi,jac_loc,jac_mat,jac_det,dn, dT, dvel, dpsi, dby,iloc)
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(n_weights_loc(oft_blagrange%nce),vel_weights_loc(3, oft_blagrange%nce),&
        T_weights_loc(oft_blagrange%nce), psi_weights_loc(oft_blagrange%nce),&
        by_weights_loc(oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce))
ALLOCATE(jac_loc(self%fe_rep%nfields,self%fe_rep%nfields))
ALLOCATE(iloc(self%fe_rep%nfields))
iloc(1)%v=>cell_dofs
iloc(2)%v=>cell_dofs
CALL self%fe_rep%mat_setup_local(jac_loc, self%jacobian_block_mask)
!$omp do schedule(static)
DO i=1,smesh%nc
  curved=cell_is_curved(smesh,i) ! Straight cell test
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
    if(curved.OR.(m==1))call smesh%jacobian(i,quad%pts(:,m),jac_mat,jac_det) ! Evaluate spatial jacobian
    !---Evaluate value and gradients of basis functions at current point
    DO jr=1,oft_blagrange%nce ! Loop over degrees of freedom
      CALL oft_blag_eval(oft_blagrange,i,jr,quad%pts(:,m),basis_vals(jr))
      CALL oft_blag_geval(oft_blagrange,i,jr,quad%pts(:,m),basis_grads(:,jr),jac_mat)
    END DO
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
      dn = dn + n_weights_loc(jr)*basis_vals(jr)
      dvel(:, 1) = dvel(:, 1) + vel_weights_loc(:, jr)*basis_grads(1, jr)
      dvel(:, 2) = 0.d0
      dvel(:, 3) = dvel(:, 2) + vel_weights_loc(:, jr)*basis_grads(3, jr)
      dT = dT + T_weights_loc(jr)*basis_grads(:,jr)
      dpsi = dpsi + psi_weights_loc(jr)*basis_grads(:,jr)
      dby = dby + by_weights_loc(jr)*basis_grads(:,jr)
    END DO
    div_vel = dvel(1,1) + dvel(2,2) + dvel(3,3)
    btmp = cross_product(dpsi, [0.d0,1.d0,0.d0]) + by*[0.d0,1.d0,0.d0]
    !---Compute local matrix contributions
    DO jr=1,oft_blagrange%nce
      DO jc=1,oft_blagrange%nce
        ! Diffusion
        ! dx = jac_det*quad%wts(m)
        jac_loc(1, 1)%m(jr,jc) = jac_loc(1, 1)%m(jr, jc) &
        + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m) &
        + self%dt*basis_vals(jr)*DOT_PRODUCT(basis_grads(:, jc), vel)*jac_det*quad%wts(m) & !delta_n*div(u)
        + self%dt*basis_vals(jr)*basis_vals(jc)*div_vel*jac_det*quad%wts(m) & ! u dot grad(delta_n)
        + self%dt*D_diff*DOT_PRODUCT(basis_grads(:, jr),basis_grads(:, jc))*jac_det*quad%wts(m)
        jac_loc(1, 2:4)%m(jr,jc) = jac_loc(1, 2:4)%m(jr, jc) &
        + basis_vals(jr)*self%dt*n*SUM(basis_grads(: jc))*jac_det*quad%wts(m) & ! n*div(delta_u) = n*SUM(basis_grads)?
        + basis_vals(jr)*self%dt*basis_vals(jc)*SUM(dn)*jac_det*quad%wts(m) ! u dot grad(n)
        !--Momentum rows
        !--Momentum,density
        jac_loc(2:4,1)%m(jr,jc) = jac_loc(2:4,1)%m(jr,jc) &
        + self%dt*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*dT*jac_det*quad%wts(m)/(n*m_i) &
        + self%dt*basis_vals(jr)*2.d0*k_boltz*T*basis_grads(:,jc)*jac_det*quad%wts(m)/(n*m_i) &
        - self%dt*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*n*dT*jac_det*quad%wts(m)/(n**2.d0*m_i) &
        - self%dt*basis_vals(jr)*2.d0*k_boltz*basis_vals(jc)*dn*T*jac_det*quad%wts(m)/(n**2.d0*m_i) &
        - self%dt*basis_vals(jc)*DOT_PRODUCT(basis_grads(:,jr), btmp)*btmp*jac_det*quad%wts(m)/(m_i*n**2.d0*mu0) &
        - self%dt*basis_vals(jr)*DOT_PRODUCT(basis_grads(:,jc), btmp)*btmp*jac_det*quad%wts(m)/(m_i*n**2.d0*mu0) &
        + self%dt*basis_vals(jr)*2.d0*basis_vals(jc)*DOT_PRODUCT(dn, btmp)*btmp*jac_det*quad%wts(m)/(m_i*n**3.d0*mu0) &
        + self%dt*basis_vals(jc)*basis_grads(:,jr)*DOT_PRODUCT(btmp, btmp)*jac_det*quad%wts(m)/(2.d0*m_i*n**2.d0*mu0) &
        + self%dt*basis_vals(jr)*basis_grads(:,jc)*DOT_PRODUCT(btmp, btmp)*jac_det*quad%wts(m)/(2.d0*m_i*n**2.d0*mu0) &
        - self%dt*basis_vals(jr)*basis_vals(jc)*dn*DOT_PRODUCT(btmp, btmp)*jac_det*quad%wts(m)/(m_i*n**3.d0*mu0) &
        DO k=1,3
          jac_loc(k+1,1)%m(jr,jc) = jac_loc(k+1,1)%m(jr,jc) &
          - self%dt*nu*basis_vals(jc)*DOT_PRODUCT(basis_grads(:,jr), dvel(:,k))*jac_det*quad%wts(m)/(m_i*n**2) & !-- not sure if indexing on dvel is right here
          - self%dt*nu*basis_vals(jr)*DOT_PRODUCT(basis_grads(:,jc), dvel(:,k))*jac_det*quad%wts(m)/(m_i*n**2) &
          + self%dt*nu*basis_vals(jr)*basis_vals(jc)*DOT_PRODUCT(dn, dvel(:,k))*jac_det*quad%wts(m)/(m_i*n**3)
        END DO
        !--Momentum, velocity
        DO k=1,3
          jac_loc(k+1,2:4)%m(jr,jc)= jac_loc(2,2:4)%m(jr,jc) &
          + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m) &
          + self%dt*basis_vals(jr)*basis_vals(jc)*dvel(k,:)*jac_det*quad%wts(m) &
          + self%dt*basis_vals(jr)*DOT_PRODUCT(vel, basis_grads(:,jc))*jac_det*quad%wts(m) &
          + self%dt*nu*DOT_PRODUCT(basis_grads(:,jr), basis_grads(:,jc))/(m_i*n) &
          - self%dt*nu*basis_vals(jr)*DOT_PRODUCT(dn, basis_grads(:,jc))/(m_i*n**2)
        END DO
        ! --Momentum, temperature
        jac_loc(2:4,5)%m(jr,jc) = jac_loc(2:4,5)%m(jr,jc) &
        + self%dt*basis_vals(jr)*2*k_boltz*dn*basis_vals(jc)*jac_det*quad%wts(m)/(m_i*n) &
        + self%dt*basis_vals(jr)*2*k_boltz*n*basis_grads(:,jc)*jac_det*quad%wts(m)/(m_i*n) &
        ! --Momentum, psi
        jac_loc(2:4,6)%m(jr,jc) = jac_loc(2:4,6)%m(jr,jc) &
        + self%dt*DOT_PRODUCT(basis_grads(:,jr), cross_product(basis_grads(:,jc), (0,1,0)))*btmp*jac_det*quad%wts(m)/(m_i*n*mu_0) &
        + self%dt*DOT_PRODUCT(basis_grads(:,jr), btmp)*cross_product(basis_grads(:,jc), (0,1,0))*jac_det*quad%wts(m)/(m_i*n*mu_0) &
        - self%dt*basis_grads(:,jr)*DOT_PRODUCT(cross_product(basis_grads(:,jc), (0,1,0)), btmp)*jac_det*quad%wts(m)/(m_i*n*mu_0) &
        - self%dt*basis_vals(jc)*DOT_PRODUCT(dn, cross_product(basis_grads(:,jc), (0,1,0)))*btmp*jac_det*quad%wts(m)/(m_i*n**2*mu_0) &
        - self%dt*basis_vals(jc)*DOT_PRODUCT(dn, btmp)*cross_product(basis_grads(:,jc), (0,1,0))*jac_det*quad%wts(m)/(m_i*n**2*mu_0) &
        + self%dt*basis_vals(jc)*dn*DOT_PRODUCT(cross_product(basis_grads(:,jc), (0,1,0)), btmp)*jac_det*quad%wts(m)/(m_i*n**2*mu_0) 
        ! --Momentum, By
        jac_loc(2:4,7)%m(jr,jc) = jac_loc(2:4,7)%m(jr,jc) &
        + self%dt*DOT_PRODUCT(basis_grads(:,jr),[0,basis_vals(jc),0])*btmp*jac_det*quad%wts(m)/(m_i*n*mu_0) &
        + self%dt*DOT_PRODUCT(basis_grads(:,jr), btmp)*[0,basis_vals(jc),0]*jac_det*quad%wts(m)/(m_i*n*mu_0) &
        - self%dt*basis_grads(:,jr)*DOT_PRODUCT([0,basis_vals(jc),0], btmp)*jac_det*quad%wts(m)/(m_i*n*mu_0) &
        - self%dt*basis_vals(jc)*DOT_PRODUCT(dn, [0,basis_vals(jc),0])*btmp*jac_det*quad%wts(m)/(m_i*n**2*mu_0) &
        - self%dt*basis_vals(jc)*DOT_PRODUCT(dn, btmp)*[0,basis_vals(jc),0]*jac_det*quad%wts(m)/(m_i*n**2*mu_0) &
        + self%dt*basis_vals(jc)*dn*DOT_PRODUCT([0,basis_vals(jc),0], btmp)*jac_det*quad%wts(m)/(m_i*n**2*mu_0) 
        ! Temperature: delta_T=basis_vals(:, jc)
        jac_loc(5, 1)%m(jr,jc) = jac_loc(5, 1)%m(jr, jc) &  
        + basis_vals(:, jr) * basis_vals(:, jc)/(gamma-1)*jac_det*quad%wts(m) &! delta_T
        + basis_vals(:, jr) * self%dt*DOT_PRODUCT(vel, basis_grads(:, jc))/(gamma-1)*jac_det*quad%wts(m) & ! nabla(dT)
        - basis_vals(:, jr) * self%dt*k_boltz*basis_vals(:, jc)*div_vel*jac_det*quad%wts(m) & ! dT != nabla(dT)
        + self%dt * chi * DOT_PRODUCT(basis_grads(:, jc),basis_grads(:, jr))*jac_det*quad%wts(m) & ! dT_Chi
        - self%dt*basis_vals(:, jr)*DOT_PRODUCT(dn, basis_grads(:, jc))/n*jac_det*quad%wts(m) 
        ! delta_u=basis_vals(:, jc)
        jac_loc(5, 2)%m(jr,jc) = jac_loc(5, 2)%m(jr, jc) &  
        + basis_vals(:, jr) * T/(gamma-1)*jac_det*quad%wts(m) &! T
        + basis_vals(:, jr) * self%dt*DOT_PRODUCT(basis_vals(:, jc), dT)/(gamma-1)*jac_det*quad%wts(m) & ! delta_u dot Delta_T
        - basis_vals(:, jr) * self%dt*k_boltz*T*SUM(basis_grads(:, jc))*jac_det*quad%wts(m) & ! div(delta u) = SUM(basis_grads)?
        + self%dt * chi * DOT_PRODUCT(dT,basis_grads(:, jr))*jac_det*quad%wts(m) & ! dT_Chi
        - self%dt*basis_vals(:, jr)*DOT_PRODUCT(dn, dT)/n*jac_det*quad%wts(m) 
        ! delta_u=basis_vals(:, jc)
        jac_loc(5, 3)%m(jr,jc) = jac_loc(5, 3)%m(jr, jc) &  
        + basis_vals(:, jr) * T/(gamma-1)*jac_det*quad%wts(m) &! T
        + basis_vals(:, jr) * self%dt*DOT_PRODUCT(vel, dT)/(gamma-1)*jac_det*quad%wts(m) & ! u dot Delta_T
        - basis_vals(:, jr) * self%dt*k_boltz*T*div_vel*jac_det*quad%wts(m) & ! div(u)
        + self%dt * chi * DOT_PRODUCT(dT,basis_grads(:, jr))*jac_det*quad%wts(m) & ! dT_Chi
        - self%dt * basis_vals(:, jr)*DOT_PRODUCT(basis_grads(:, jc), dT)/n*jac_det*quad%wts(m) ! grad(delta_n)
        + self%dt * basis_vals(:, jr)*basis_grads(:, jc)*DOT_PRODUCT(dn, dT)/(n*n)*jac_det*quad%wts(m) ! delta_n
        ! Induction
        ! Psi rows
        jac_loc(6, 1)%m(jr,jc) = jac_loc(6, 1)%m(jr,jc) &
        + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m) &
        + self%dt*basis_vals(jr)*DOT_PRODUCT(vel,basis_grads(:,jc))*jac_det*quad%wts(m) &
        + self%dt*eta*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*jac_det*quad%wts(m)/mu_0
        jac_loc(6,2)%m(jr,jc) = jac_loc(6,2)%m(jr,jc) &
        + basis_vals(jr)*self%dt*basis_vals(jc)*SUM(dpsi)*jac_det*quad%wts(m)
        ! B_y rows
        jac_loc(7, 1)%m(jr,jc) = jac_loc(7, 1)%m(jr,jc) &
        - basis_vals(jr)*self%dt*cross_product(basis_grads(:,jc),dvel(2,:))*jac_det*quad%wts(m)
        jac_loc(7, 2)%m(jr,jc) = jac_loc(7, 2)%m(jr,jc) &
        - basis_vals(jr)*self%dt*cross_product(dpsi,basis_grads(:,jc))*jac_det*quad%wts(m) &
        + basis_vals(jr)*self%dt*basis_vals(jc)*SUM(dby)*jac_det*quad%wts(m) &
        + basis_vals(jr)*self%dt*by*SUM(basis_grads(:,jc))*jac_det*quad%wts(m) &
        jac_loc(7, 3)%m(jr,jc) = jac_loc(7, 3)%m(jr,jc) &
        + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m) &
        + basis_vals(jr)*self%dt*DOT_PRODUCT(basis_grads(:,jc),vel)*jac_det*quad%wts(m) &
        + basis_vals(jr)*self%dt*basis_vals(jc)*div_vel*jac_det*quad%wts(m) &
        + self%dt*eta*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*jac_det*quad%wts(m)/mu_0
      END DO
    END DO
  END DO
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%n_bc(cell_dofs),1)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%vel_bc(cell_dofs), 2)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%T_bc(cell_dofs),3)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%psi_bc(cell_dofs),4)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%by_bc(cell_dofs), 5)
  CALL self%fe_rep%mat_add_local(self%jacobian,jac_loc,iloc,tlocks)
END DO
!---Cleanup thread-local storage
CALL self%fe_rep%mat_destroy_local(jac_loc)
DEALLOCATE(basis_vals,basis_grads,T_weights_loc,vel_weights_loc&
          n_weights_loc, psi_weights_loc, by_weights_loc, cell_dofs,jac_loc,iloc, btmp)
!$omp end parallel
!--Destroy thread locks
DO i=1,self%fe_rep%nfields
  CALL omp_destroy_lock(tlocks(i))
END DO
DEALLOCATE(tlocks)
IF(oft_debug_print(2))write(*,'(4X,A)')'Setting BCs'
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%n_bc,1)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%vel_bc(1,:),2)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%vel_bc(2,:),3)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%vel_bc(3,:),4)
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
  IF(oft_debug_print(1))write(*,*)'Updating tDiff MF-Jacobian'
  CALL current_sim%mf_mat%update(uin)
  END SUBROUTINE mfnk_update
  !---------------------------------------------------------------------------
  !> Update Jacobian matrices on all levels with new solution
  !---------------------------------------------------------------------------
  subroutine update_jacobian(uin)
  class(oft_vector), target, intent(inout) :: uin !< Current solution
  IF(oft_debug_print(1))write(*,*)'Updating tDiff approximate Jacobian'
  CALL build_approx_jacobian(current_sim,uin)
  END SUBROUTINE update_jacobian
  !---------------------------------------------------------------------------
  !> Setup composite FE representation and ML environment
  !---------------------------------------------------------------------------
  subroutine setup(self,order)
  class(oft_xmhd_2d_sim), intent(inout) :: self
  integer(i4), intent(in) :: order
  integer(i4) :: ierr,io_unit
  IF(ASSOCIATED(self%fe_rep))CALL oft_abort("Setup can only be called once","setup",__FILE__)
  IF(ASSOCIATED(oft_blagrange))CALL oft_abort("FE space already built","setup",__FILE__)
  
  !---Look for XML defintion elements
#ifdef HAVE_XML
  IF(ASSOCIATED(oft_env%xml))THEN
    CALL xml_get_element(oft_env%xml,"tdiff",self%xml_root,ierr)
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
  CALL oft_lag_setup(order, -1)
  CALL self%xdmf_plot%setup("tdiff")
  CALL smesh%setup_io(self%xdmf_plot,order)
  
  !---Build composite FE definition for solution field
  IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Creating FE type'
  ALLOCATE(self%fe_rep)
  self%fe_rep%nfields=2
  ALLOCATE(self%fe_rep%fields(self%fe_rep%nfields))
  ALLOCATE(self%fe_rep%field_tags(self%fe_rep%nfields))
  self%fe_rep%fields(1)%fe=>oft_blagrange
  self%fe_rep%field_tags(1)='Ti'
  self%fe_rep%fields(2)%fe=>oft_blagrange
  self%fe_rep%field_tags(2)='Te'
  
  !---Create solution vector
  CALL self%fe_rep%vec_create(self%u)
  
!---Set boundary conditions (Dirichlet for now)
self%n_bc=>oft_blagrange%be
self%T_bc=>oft_blagrange%be
self%vel_bc=>oft_blagrange%be
self%psi_bc=>oft_blagrange%be
self%by_bc=>oft_blagrange%be
  
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