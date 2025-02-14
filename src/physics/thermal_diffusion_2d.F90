!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file thermal_diffusion_2d.F90
!
!> Example module for modeling thermal diffusion of two species, with equilibration, in 2D
!!
!! The equation to be solved is a simple 
!!
!! \f[ \frac{\partial T_i}{\partial t} = kappa_i T_i^(5/2) \nabla^2 T_i  + \tau_eq (T_e - T_i) + Q_i \f]
!!
!! \f[ \frac{\partial T_e}{\partial t} = kappa_e T_e^(5/2) \nabla^2 T_e  + \tau_eq (T_i - T_e) + Q_e \f]
!!
!! @authors Chris Hansen
!! @date January 2025
!! @ingroup doxy_oft_physics
!---------------------------------------------------------------------------
MODULE thermal_diffusion_2d
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
#if !defined(TDIFF_RST_LEN)
#define TDIFF_RST_LEN 5
#endif
PRIVATE
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
TYPE, extends(oft_noop_matrix) :: tdiff_nlfun
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
END TYPE tdiff_nlfun
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
TYPE, public :: oft_tdiff_sim
  LOGICAL :: mfnk = .FALSE. !< Use matrix free method?
  INTEGER(i4) :: nsteps = -1 !< Needs docs
  INTEGER(i4) :: rst_base = 0 !< Needs docs
  INTEGER(i4) :: rst_freq = 10 !< Needs docs
  REAL(r8) :: dt = -1.d0 !< Needs docs
  REAL(r8) :: t = 0.d0 !< Needs docs
  REAL(r8) :: kappa_e = -1.d0 !< Needs docs
  REAL(r8) :: kappa_i = -1.d0 !< Needs docs
  REAL(r8) :: tau_eq = 0.d0 !< Needs docs
  REAL(r8) :: lin_tol = 1.d-8 !< Needs docs
  REAL(r8) :: nl_tol = 1.d-5 !< Needs docs
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: Ti_bc => NULL() !< Ti BC flag
  LOGICAL, CONTIGUOUS, POINTER, DIMENSION(:) :: Te_bc => NULL() !< Te BC flag
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:,:) :: jacobian_block_mask => NULL() !< Matrix block mask
  TYPE(oft_fem_comp_type), POINTER :: fe_rep => NULL() !< Finite element representation for solution field
  CLASS(oft_vector), POINTER :: u => NULL() !< Needs docs
  CLASS(oft_matrix), POINTER :: jacobian => NULL() !< Needs docs
  TYPE(oft_mf_matrix), POINTER :: mf_mat => NULL() !< Matrix free operator
  TYPE(tdiff_nlfun), POINTER :: nlfun => NULL() !< Needs docs
  TYPE(fox_node), POINTER :: xml_root => NULL() !< XML root element
  TYPE(fox_node), POINTER :: xml_pre_def => NULL() !< XML element for preconditioner definition
contains
  !> Setup
  PROCEDURE :: setup => setup
  !> Run simulation
  PROCEDURE :: run_simulation => run_simulation
  !> Save restart file
  PROCEDURE :: rst_save => rst_save
  !> Load restart file
  PROCEDURE :: rst_load => rst_load
END TYPE oft_tdiff_sim
TYPE(oft_tdiff_sim), POINTER :: current_sim => NULL()
CONTAINS
!---------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------
subroutine run_simulation(self)
class(oft_tdiff_sim), target, intent(inout) :: self
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
real(r8) :: Ti_avg,Te_avg,elapsed_time
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
CALL hdf5_create_timestep(self%t)
CALL self%u%get_local(plot_vals,1)
CALL smesh%save_vertex_scalar(plot_vals,'Ti')
CALL self%u%get_local(plot_vals,2)
CALL smesh%save_vertex_scalar(plot_vals,'Te')

!
ALLOCATE(self%nlfun)
self%nlfun%kappa_e=self%kappa_e
self%nlfun%kappa_i=self%kappa_i
self%nlfun%tau_eq=self%tau_eq
self%nlfun%Ti_bc=>self%Ti_bc
self%nlfun%Te_bc=>self%Te_bc

!---
CALL build_approx_jacobian(self,u)

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
    CALL hdf5_create_timestep(self%t)
    CALL self%u%get_local(plot_vals,1)
    CALL smesh%save_vertex_scalar(plot_vals,'Ti')
    CALL self%u%get_local(plot_vals,2)
    CALL smesh%save_vertex_scalar(plot_vals,'Te')
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
class(tdiff_nlfun), intent(inout) :: self !< NL function object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
type(oft_quad_type), pointer :: quad
LOGICAL :: curved
INTEGER(i4) :: i,m,jr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cell_dofs
REAL(r8) :: kappa_e,kappa_i,tau_eq_inv,diag_vals(2)
REAL(r8) :: T_i,T_e,dT_i(3),dT_e(3),jac_mat(3,4),jac_det
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals,Ti_weights_loc,Te_weights_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads,res_loc
REAL(r8), POINTER, DIMENSION(:) :: Ti_weights,Te_weights,Ti_res,Te_res
quad=>oft_blagrange%quad
NULLIFY(Ti_weights,Te_weights,Ti_res,Te_res)
!---Get weights from solution vector
CALL a%get_local(Ti_weights,1)
CALL a%get_local(Te_weights,2)
!---
kappa_e=self%kappa_e
kappa_i=self%kappa_i
IF(self%tau_eq>0.d0)THEN
  tau_eq_inv=1.d0/self%tau_eq
ELSE
  tau_eq_inv=1.d0
END IF
!---Zero result and get storage array
CALL b%set(0.d0)
CALL b%get_local(Ti_res,1)
CALL b%get_local(Te_res,2)
diag_vals=0.d0
!$omp parallel private(m,jr,curved,cell_dofs,basis_vals,basis_grads,Ti_weights_loc, &
!$omp Te_weights_loc,res_loc,jac_mat,jac_det,T_i,T_e,dT_i,dT_e) reduction(+:diag_vals)
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(Ti_weights_loc(oft_blagrange%nce),Te_weights_loc(oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce),res_loc(oft_blagrange%nce,2))
!$omp do schedule(static)
DO i=1,smesh%nc
  curved=cell_is_curved(smesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  res_loc = 0.d0 ! Zero local (cell) contribution to function
  Ti_weights_loc = Ti_weights(cell_dofs)
  Te_weights_loc = Te_weights(cell_dofs)
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
    T_i = 0.d0; dT_i = 0.d0; T_e = 0.d0; dT_e = 0.d0
    DO jr=1,oft_blagrange%nce
      T_i = T_i + Ti_weights_loc(jr)*basis_vals(jr)
      T_e = T_e + Te_weights_loc(jr)*basis_vals(jr)
      dT_i = dT_i + Ti_weights_loc(jr)*basis_grads(:,jr)
      dT_e = dT_e + Te_weights_loc(jr)*basis_grads(:,jr)
    END DO
    diag_vals = diag_vals + [T_i,T_e]*jac_det*quad%wts(m)
    !---Compute local function contributions
    DO jr=1,oft_blagrange%nce
      res_loc(jr,1) = res_loc(jr,1) &
        + basis_vals(jr)*T_i*jac_det*quad%wts(m) &
        + self%dt*kappa_i*(T_i**2.5d0)*DOT_PRODUCT(basis_grads(:,jr),dT_i)*jac_det*quad%wts(m) &
        - self%dt*tau_eq_inv*basis_vals(jr)*(T_e-T_i)*jac_det*quad%wts(m)
      res_loc(jr,2) = res_loc(jr,2) &
        + basis_vals(jr)*T_e*jac_det*quad%wts(m) &
        + self%dt*kappa_e*(T_e**2.5d0)*DOT_PRODUCT(basis_grads(:,jr),dT_e)*jac_det*quad%wts(m) &
        - self%dt*tau_eq_inv*basis_vals(jr)*(T_i-T_e)*jac_det*quad%wts(m)
    END DO
  END DO
  !---Add local values to full vector
  DO jr=1,oft_blagrange%nce
    !$omp atomic
    Ti_res(cell_dofs(jr)) = Ti_res(cell_dofs(jr)) + res_loc(jr,1)
    !$omp atomic
    Te_res(cell_dofs(jr)) = Te_res(cell_dofs(jr)) + res_loc(jr,2)
  END DO
END DO
!---Cleanup thread-local storage
DEALLOCATE(basis_vals,basis_grads,Ti_weights_loc,Te_weights_loc,cell_dofs,res_loc)
!$omp end parallel
IF(oft_debug_print(2))write(*,'(4X,A)')'Applying BCs'
CALL fem_dirichlet_vec(oft_blagrange,Ti_weights,Ti_res,self%Ti_bc)
CALL fem_dirichlet_vec(oft_blagrange,Te_weights,Te_res,self%Te_bc)
!---Put results into full vector
CALL b%restore_local(Ti_res,1,add=.TRUE.,wait=.TRUE.)
CALL b%restore_local(Te_res,2,add=.TRUE.)
self%diag_vals=oft_mpi_sum(diag_vals,2)
!---Cleanup remaining storage
DEALLOCATE(Ti_res,Te_res,Ti_weights,Te_weights)
end subroutine nlfun_apply
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine build_approx_jacobian(self,a)
class(oft_tdiff_sim), intent(inout) :: self
class(oft_vector), intent(inout) :: a !< Solution for computing jacobian
LOGICAL :: curved
INTEGER(i4) :: i,m,jr,jc
INTEGER(i4), POINTER, DIMENSION(:) :: cell_dofs
REAL(r8) :: kappa_e,kappa_i,tau_eq_inv
REAL(r8) :: T_i,T_e,dT_i(3),dT_e(3),jac_mat(3,4),jac_det
REAL(r8), ALLOCATABLE, DIMENSION(:) :: basis_vals,Ti_weights_loc,Te_weights_loc,res_loc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: basis_grads
REAL(r8), POINTER, DIMENSION(:) :: Ti_weights,Te_weights
type(oft_1d_int), allocatable, dimension(:) :: iloc
class(oft_vector), pointer :: tmp
type(oft_local_mat), allocatable, dimension(:,:) :: jac_loc
integer(KIND=omp_lock_kind), allocatable, dimension(:) :: tlocks
type(oft_quad_type), pointer :: quad
quad=>oft_blagrange%quad
CALL self%jacobian%zero
NULLIFY(Ti_weights,Te_weights)
!---Get weights from solution vector
CALL a%get_local(Ti_weights,1)
CALL a%get_local(Te_weights,2)
!---
kappa_e=self%kappa_e
kappa_i=self%kappa_i
IF(self%tau_eq>0.d0)THEN
  tau_eq_inv=1.d0/self%tau_eq
ELSE
  tau_eq_inv=1.d0
END IF
!--Setup thread locks
ALLOCATE(tlocks(self%fe_rep%nfields))
DO i=1,self%fe_rep%nfields
  call omp_init_lock(tlocks(i))
END DO
!$omp parallel private(m,jr,jc,curved,cell_dofs,basis_vals,basis_grads,Ti_weights_loc, &
!$omp Te_weights_loc,jac_loc,jac_mat,jac_det,T_i,T_e,dT_i,dT_e,iloc)
ALLOCATE(basis_vals(oft_blagrange%nce),basis_grads(3,oft_blagrange%nce))
ALLOCATE(Ti_weights_loc(oft_blagrange%nce),Te_weights_loc(oft_blagrange%nce))
ALLOCATE(cell_dofs(oft_blagrange%nce))
allocate(jac_loc(self%fe_rep%nfields,self%fe_rep%nfields))
allocate(iloc(self%fe_rep%nfields))
iloc(1)%v=>cell_dofs
iloc(2)%v=>cell_dofs
CALL self%fe_rep%mat_setup_local(jac_loc, self%jacobian_block_mask)
!$omp do schedule(static)
DO i=1,smesh%nc
  curved=cell_is_curved(smesh,i) ! Straight cell test
  call oft_blagrange%ncdofs(i,cell_dofs) ! Get global index of local DOFs
  CALL self%fe_rep%mat_zero_local(jac_loc) ! Zero local (cell) contribution to matrix
  Ti_weights_loc = Ti_weights(cell_dofs)
  Te_weights_loc = Te_weights(cell_dofs)
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
    T_i = 0.d0; dT_i = 0.d0; T_e = 0.d0; dT_e = 0.d0
    DO jr=1,oft_blagrange%nce
      T_i = T_i + Ti_weights_loc(jr)*basis_vals(jr)
      T_e = T_e + Te_weights_loc(jr)*basis_vals(jr)
      dT_i = dT_i + Ti_weights_loc(jr)*basis_grads(:,jr)
      dT_e = dT_e + Te_weights_loc(jr)*basis_grads(:,jr)
    END DO
    !---Compute local matrix contributions
    DO jr=1,oft_blagrange%nce
      DO jc=1,oft_blagrange%nce
        !---Ion rows
        jac_loc(1,1)%m(jr,jc) = jac_loc(1,1)%m(jr,jc) &
          + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m) &
          + self%dt*kappa_i*(T_i**2.5d0)*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*jac_det*quad%wts(m) &
          + self%dt*kappa_i*(T_i**1.5d0)*basis_vals(jc)*DOT_PRODUCT(basis_grads(:,jr),dT_i)*jac_det*quad%wts(m) &
          - self%dt*tau_eq_inv*basis_vals(jr)*(-basis_vals(jc))*jac_det*quad%wts(m)
        jac_loc(1,2)%m(jr,jc) = jac_loc(1,2)%m(jr,jc) &
          - self%dt*tau_eq_inv*basis_vals(jr)*(basis_vals(jc))*jac_det*quad%wts(m)
        !---Electron rows
        jac_loc(2,2)%m(jr,jc) = jac_loc(2,2)%m(jr,jc) &
          + basis_vals(jr)*basis_vals(jc)*jac_det*quad%wts(m) &
          + self%dt*kappa_e*(T_e**2.5d0)*DOT_PRODUCT(basis_grads(:,jr),basis_grads(:,jc))*jac_det*quad%wts(m) &
          + self%dt*kappa_e*(T_e**1.5d0)*basis_vals(jc)*DOT_PRODUCT(basis_grads(:,jr),dT_e)*jac_det*quad%wts(m) &
          - self%dt*tau_eq_inv*basis_vals(jr)*(-basis_vals(jc))*jac_det*quad%wts(m)
        jac_loc(2,1)%m(jr,jc) = jac_loc(2,1)%m(jr,jc) &
          - self%dt*tau_eq_inv*basis_vals(jr)*(basis_vals(jc))*jac_det*quad%wts(m)
      END DO
    END DO
  END DO
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%Ti_bc(cell_dofs),1)
  CALL self%fe_rep%mat_zero_local_rows(jac_loc,self%Te_bc(cell_dofs),2)
  CALL self%fe_rep%mat_add_local(self%jacobian,jac_loc,iloc,tlocks)
END DO
!---Cleanup thread-local storage
CALL self%fe_rep%mat_destroy_local(jac_loc)
DEALLOCATE(basis_vals,basis_grads,Ti_weights_loc,Te_weights_loc,cell_dofs,jac_loc,iloc)
!$omp end parallel
!--Destroy thread locks
DO i=1,self%fe_rep%nfields
  CALL omp_destroy_lock(tlocks(i))
END DO
DEALLOCATE(tlocks)
IF(oft_debug_print(2))write(*,'(4X,A)')'Setting BCs'
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%Ti_bc,1)
CALL fem_dirichlet_diag(oft_blagrange,self%jacobian,self%Te_bc,2)
!
call self%fe_rep%vec_create(tmp)
call self%jacobian%assemble(tmp)
call tmp%delete
DEALLOCATE(tmp,Ti_weights,Te_weights)
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
class(oft_tdiff_sim), intent(inout) :: self
integer(i4), intent(in) :: order
integer(i4) :: ierr,io_unit
#ifdef HAVE_XML
integer(i4) :: nnodes
TYPE(fox_nodelist), POINTER :: current_nodes
#endif
IF(ASSOCIATED(self%fe_rep))CALL oft_abort("Setup can only be called once","setup",__FILE__)
IF(ASSOCIATED(oft_blagrange))CALL oft_abort("FE space already built","setup",__FILE__)

!---Look for XML defintion elements
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml))THEN
  current_nodes=>fox_getElementsByTagName(oft_env%xml,"tdiff")
  nnodes=fox_getLength(current_nodes)
  IF(nnodes>0)THEN
    self%xml_root=>fox_item(current_nodes,0)
    !---Look for pre node
    current_nodes=>fox_getElementsByTagName(self%xml_root,"pre")
    nnodes=fox_getLength(current_nodes)
    IF(nnodes>0)self%xml_pre_def=>fox_item(current_nodes,0)
  END IF
END IF
#endif

!---Setup FE representation
IF(oft_debug_print(1))WRITE(*,'(2X,A)')'Building lagrange FE space'
CALL oft_lag_setup(order, -1)
CALL smesh%setup_io(order)

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
self%Ti_bc=>oft_blagrange%be
self%Te_bc=>oft_blagrange%be

!---Create Jacobian matrix
ALLOCATE(self%jacobian_block_mask(self%fe_rep%nfields,self%fe_rep%nfields))
self%jacobian_block_mask=1
CALL self%fe_rep%mat_create(self%jacobian,self%jacobian_block_mask)
end subroutine setup
!---------------------------------------------------------------------------
!> Save xMHD solution state to a restart file
!---------------------------------------------------------------------------
subroutine rst_save(self,u,t,dt,filename,path)
class(oft_tdiff_sim), intent(inout) :: self
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
class(oft_tdiff_sim), intent(inout) :: self
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
END MODULE thermal_diffusion_2d