!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (OpenFUSIONToolkit)
!---------------------------------------------------------------------------
!> @file oft_gs_td.F90
!
!> Time-dependent G-S equilibria
!!
!! @authors Chris Hansen
!! @date May 2023
!! @ingroup doxy_oft_physics
!---------------------------------------------------------------------------
MODULE oft_gs_td
USE oft_base
USE oft_sort, ONLY: sort_array, search_array
USE oft_mesh_type, ONLY: smesh, bmesh_findcell
USE oft_mesh_local_util, ONLY: mesh_local_findedge
USE oft_quadrature, ONLY: oft_quad_type
USE oft_gauss_quadrature, ONLY: set_quad_1d
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_graph_ptr, oft_vector_ptr
USE oft_la_utils, ONLY: create_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_deriv_matrices, ONLY: oft_noop_matrix, oft_mf_matrix
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_lu, ONLY: oft_lusolver
#ifdef HAVE_ARPACK
USE oft_arpack, ONLY: oft_iram_eigsolver
#endif
USE fem_utils, ONLY: bfem_map_flag
USE oft_lag_basis, ONLY: oft_blagrange, oft_blag_geval, oft_blag_eval, oft_blag_npos
USE oft_blag_operators, ONLY: blag_zerob, oft_lag_brinterp
USE axi_green, ONLY: green, grad_green
USE oft_gs, ONLY: gs_epsilon, flux_func, gs_eq, gs_update_bounds, &
    gs_test_bounds, gs_mat_create, compute_bcmat, set_bcmat, gs_zerob
USE mhd_utils, ONLY: mu0
IMPLICIT NONE
#include "local.h"
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(oft_noop_matrix) :: oft_tmaker_td
    logical :: has_plasma = .TRUE.
    integer(4) :: nrhs = 0
    real(r8) :: dt = 1.E-3_r8 !< Time step
    real(r8) :: lam_amp = 1.E-3_r8 !< Scale factor of current source term
    real(r8) :: psi_lim = 0.d0
    real(r8) :: psi_max = 1.d0
    real(r8) :: f_scale = 1.d0
    real(r8) :: p_scale = 1.d0
    real(r8) :: ip = 0.d0
    real(r8) :: estored = 0.d0
    real(r8) :: ip_target = -1.d0
    real(r8) :: ip_ratio_target = -1.d0
    real(r8) :: lim_pt(2) = 0.d0
    real(r8) :: o_point(2) = 0.d0
    logical, pointer, dimension(:) :: fe_flag => NULL()
    integer(4), pointer, dimension(:) :: rhs_list => NULL()
    real(8), pointer, dimension(:) :: eta_reg => NULL()
    real(8), pointer, dimension(:) :: curr_reg => NULL()
    real(8), pointer, dimension(:,:) :: bc_lmat => NULL()
    real(8), pointer, dimension(:,:) :: bc_bmat => NULL()
    CLASS(flux_func), POINTER :: F => NULL() !<
    CLASS(flux_func), POINTER :: P => NULL() !<
    TYPE(gs_eq), POINTER :: gs_eq => NULL()
    CLASS(oft_matrix), POINTER :: vac_op => NULL()
contains
    !> Apply the matrix
    procedure :: apply_real => apply_mfop
    !
    ! procedure :: update_lims => update_lims
end type oft_tmaker_td
type, extends(oft_noop_matrix) :: tMaker_td_mat
    integer(4), pointer, dimension(:) :: lim_nodes => NULL()
    integer(4), pointer, dimension(:) :: ax_nodes => NULL()
    real(8), pointer, dimension(:,:) :: lim_vals => NULL()
    real(8), pointer, dimension(:,:) :: ax_vals => NULL()
    class(oft_matrix), pointer :: mat => NULL()
contains
    procedure :: apply_real => apply_gs_mat
end type tMaker_td_mat
type, extends(oft_noop_matrix) :: eig_wrapper
    class(oft_matrix), pointer :: rhs_mat => NULL()
    class(oft_solver), pointer :: lhs_inv => NULL()
contains
    procedure :: apply_real => apply_wrap
end type eig_wrapper
TYPE(oft_mf_matrix), TARGET :: mfmat !< Matrix free operator
TYPE(tMaker_td_mat), TARGET :: adv_op
TYPE(oft_tmaker_td), TARGET :: tMaker_td_obj
TYPE(oft_graph_ptr) :: graphs(1,1)
! TYPE(gs_eq) :: mygs
CLASS(oft_matrix), POINTER :: mr2op,mrop,dels
CLASS(oft_solver), POINTER :: mrop_solver
type(oft_native_gmres_solver), target :: mf_solver,adv_solver
TYPE(oft_lusolver), TARGET :: adv_pre,dels_solver
CLASS(oft_vector), POINTER :: rhs1,rhs2,psi_sol,psi_tmp,tmp_vec
type(oft_nksolver) :: nksolver
CHARACTER(LEN=3) :: coil_tag
!---Extrapolation fields
integer(i4), parameter :: maxextrap = 2
integer(i4) :: nextrap
real(r8), allocatable, dimension(:) :: extrapt
type(oft_vector_ptr), allocatable, dimension(:) :: extrap_fields
real(8) :: mfop_time
CONTAINS
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine run_gs_td(eq_in,dt,nsteps,lin_tol,nl_tol)
TYPE(gs_eq), TARGET, INTENT(inout) :: eq_in
INTEGER(4), INTENT(in) :: nsteps
REAL(8), INTENT(in) :: dt,lin_tol,nl_tol
INTEGER(4) :: i,j,k,ierr,io_unit,npts,iostat
LOGICAL :: file_exists
LOGICAL, ALLOCATABLE, DIMENSION(:) :: vel_bc
REAL(8) :: err_tmp,ip_target,p_scale,f_scale,time_val,pt(2)
REAL(8), ALLOCATABLE, DIMENSION(:) :: np_count,res_in
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: pts
REAL(8), POINTER, DIMENSION(:) :: vals_out,vals2,rhs_tmp
!
tMaker_td_obj%gs_eq=>eq_in
tMaker_td_obj%ip_target=tMaker_td_obj%gs_eq%Itor_target
tMaker_td_obj%f_scale=tMaker_td_obj%gs_eq%alam
tMaker_td_obj%p_scale=tMaker_td_obj%gs_eq%pnorm
tMaker_td_obj%dt=dt
!
! CALL gs_profile_load('f_prof.in',tMaker_td_obj%gs_eq%I)
tMaker_td_obj%F=>tMaker_td_obj%gs_eq%I
tMaker_td_obj%P=>tMaker_td_obj%gs_eq%P
! CALL gs_profile_load('p_prof.in',tMaker_td_obj%gs_eq%P)
!---------------------------------------------------------------------------
! Create Solver fields
!---------------------------------------------------------------------------
NULLIFY(vals_out)
! call oft_blagrange%vec_create(psi_sol)
psi_sol=>tMaker_td_obj%gs_eq%psi
call oft_blagrange%vec_create(rhs1)
call oft_blagrange%vec_create(rhs2)
call oft_blagrange%vec_create(psi_tmp)
call oft_blagrange%vec_create(tmp_vec)
! CALL psi_sol%add(0.d0,1.d0,tMaker_td_obj%gs_eq%psi)
!---------------------------------------------------------------------------
! Create extrapolation fields (Unused)
!---------------------------------------------------------------------------
IF(maxextrap>0)THEN
  ALLOCATE(extrap_fields(maxextrap),extrapt(maxextrap))
  DO i=1,maxextrap
    CALL oft_blagrange%vec_create(extrap_fields(i)%f)
    extrapt(i)=0.d0
  END DO
  nextrap=0
END IF
!---
ALLOCATE(tMaker_td_obj%eta_reg(smesh%nreg))
tMaker_td_obj%eta_reg=-1.d0
DO i=1,tMaker_td_obj%gs_eq%ncond_regs
  j=tMaker_td_obj%gs_eq%cond_regions(i)%id
  tMaker_td_obj%eta_reg(j)=tMaker_td_obj%gs_eq%cond_regions(i)%eta
END DO
ALLOCATE(tMaker_td_obj%curr_reg(smesh%nreg))
tMaker_td_obj%curr_reg=0.d0
DO i=1,tMaker_td_obj%gs_eq%ncoil_regs
  j=tMaker_td_obj%gs_eq%coil_regions(i)%id
  tMaker_td_obj%curr_reg(j)=tMaker_td_obj%gs_eq%coil_regions(i)%curr
END DO
! Advance using MF-NK method
CALL build_vac_op(tMaker_td_obj,tMaker_td_obj%vac_op)
CALL build_jop(tMaker_td_obj,adv_op,psi_sol)
adv_solver%A=>adv_op
adv_solver%its=10
adv_solver%nrits=10
adv_solver%pre=>adv_pre
adv_pre%A=>tMaker_td_obj%vac_op
!
CALL mfmat%setup(psi_tmp,tMaker_td_obj)
! CALL mfmat%utyp%delete()
! DEALLOCATE(mfmat%utyp)
mfmat%b0=1.d-5
mf_solver%A=>mfmat
mf_solver%its=100
mf_solver%nrits=20
mf_solver%atol=lin_tol
mf_solver%itplot=1
mf_solver%pm=oft_env%pm
! CALL create_diag_pre(mf_solver%pre)
! mf_solver%pre%A=>dels
mf_solver%pre=>adv_solver
!
nksolver%A=>tMaker_td_obj
nksolver%J_inv=>mf_solver
nksolver%its=20
nksolver%atol=nl_tol
nksolver%rtol=1.d-20 ! Disable relative tolerance
nksolver%backtrack=.FALSE.
nksolver%J_update=>tMaker_td_mfnk_update
nksolver%up_freq=1
time_val=0.d0
DO i=1,nsteps
    !
    ! Advance B
    !
    CALL psi_tmp%add(0.d0,1.d0,psi_sol)
    CALL apply_rhs(tMaker_td_obj,psi_sol,rhs1)
    CALL blag_zerob(rhs1)
    ! ! Extrapolate solution (linear)
    ! DO j=maxextrap,2,-1
    !   CALL extrap_fields(j)%f%add(0.d0,1.d0,extrap_fields(j-1)%f)
    !   extrapt(j)=extrapt(j-1)
    ! END DO
    ! IF(nextrap<maxextrap)nextrap=nextrap+1
    ! CALL extrap_fields(1)%f%add(0.d0,1.d0,psi_sol)
    ! extrapt(1)=time_val

    ! !---Rebuild plasma matrix
    ! IF(MOD(i,2)==0)THEN
    !   CALL build_jop(adv_op,psi_sol)
    !   CALL adv_pre%update(.TRUE.)
    ! END IF
    DO j=1,4
        ! CALL vector_extrapolate(extrapt,extrap_fields,nextrap,time_val+tMaker_td_obj%dt,psi_sol)
        !---MFNK iteration
        ! CALL tMaker_td_mfnk_update(psi_sol)
        CALL blag_zerob(psi_sol)
        CALL nksolver%apply(psi_sol,rhs1)
        IF(nksolver%cits<0)THEN
            CALL psi_sol%add(0.d0,1.d0,psi_tmp)
            tMaker_td_obj%dt=tMaker_td_obj%dt/2.d0
            CALL build_vac_op(tMaker_td_obj,tMaker_td_obj%vac_op)
            CALL build_jop(tMaker_td_obj,adv_op,psi_sol)
            CALL adv_pre%update(.TRUE.)
            CALL apply_rhs(tMaker_td_obj,psi_sol,rhs1)
            CALL blag_zerob(rhs1)
            CYCLE
        ELSE
            EXIT
        END IF
    END DO
    IF(j>4)CALL oft_abort("Nonlinear solve failed", "run_gs_td",__FILE__)
    !
    CALL psi_sol%get_local(vals_out)
    CALL psi_tmp%add(1.d0,-1.d0,psi_sol)
    CALL gs_zerob(psi_tmp)
    err_tmp=SQRT(psi_tmp%dot(psi_tmp))
    CALL psi_tmp%add(0.d0,1.d0,psi_sol)
    CALL gs_zerob(psi_tmp)
    err_tmp=err_tmp/SQRT(psi_tmp%dot(psi_tmp))
    time_val=time_val+tMaker_td_obj%dt
    WRITE(*,'(4ES14.6,3I4)')time_val,tMaker_td_obj%ip,tMaker_td_obj%gs_eq%o_point(2),err_tmp, &
        nksolver%nlits,nksolver%lits,j
    IF(tMaker_td_obj%ip<1.d-6)THEN
        tMaker_td_obj%p_scale=0.d0
        tMaker_td_obj%f_scale=0.d0
        tMaker_td_obj%ip_target=-1.d0
    END IF
    ! !---Plot
    ! IF(MOD(i,nplot)==0)THEN
    !     CALL hdf5_create_timestep(time_val)
    !     CALL psi_sol%get_local(vals_out)
    !     CALL smesh%save_vertex_scalar(vals_out,'Psi')
    !     CALL smesh%save_vertex_scalar(vals_out-tMaker_td_obj%F%plasma_bounds(1),'Psi_p')
    ! END IF
END DO
end subroutine run_gs_td
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine setup_gs_td(eq_in,dt,lin_tol,nl_tol)
TYPE(gs_eq), TARGET, INTENT(inout) :: eq_in
REAL(8), INTENT(in) :: dt,lin_tol,nl_tol
INTEGER(4) :: i,j,k,ierr,io_unit,npts,iostat
LOGICAL :: file_exists
LOGICAL, ALLOCATABLE, DIMENSION(:) :: vel_bc
REAL(8) :: err_tmp,ip_target,p_scale,f_scale,time_val,pt(2)
REAL(8), ALLOCATABLE, DIMENSION(:) :: np_count,res_in
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: pts
REAL(8), POINTER, DIMENSION(:) :: vals_out,vals2,rhs_tmp
!
tMaker_td_obj%gs_eq=>eq_in
tMaker_td_obj%ip_target=tMaker_td_obj%gs_eq%Itor_target
tMaker_td_obj%ip_ratio_target=tMaker_td_obj%gs_eq%Ip_ratio_target
tMaker_td_obj%f_scale=tMaker_td_obj%gs_eq%alam
tMaker_td_obj%p_scale=tMaker_td_obj%gs_eq%pnorm
tMaker_td_obj%dt=dt
!
tMaker_td_obj%F=>tMaker_td_obj%gs_eq%I
tMaker_td_obj%P=>tMaker_td_obj%gs_eq%P
!---------------------------------------------------------------------------
! Create Solver fields
!---------------------------------------------------------------------------
NULLIFY(vals_out)
psi_sol=>tMaker_td_obj%gs_eq%psi
call oft_blagrange%vec_create(rhs1)
call oft_blagrange%vec_create(rhs2)
call oft_blagrange%vec_create(psi_tmp)
call oft_blagrange%vec_create(tmp_vec)
!---------------------------------------------------------------------------
! Create extrapolation fields (Unused)
!---------------------------------------------------------------------------
IF(maxextrap>0)THEN
    ALLOCATE(extrap_fields(maxextrap),extrapt(maxextrap))
    DO i=1,maxextrap
    CALL oft_blagrange%vec_create(extrap_fields(i)%f)
    extrapt(i)=0.d0
    END DO
    nextrap=0
END IF
!---
ALLOCATE(tMaker_td_obj%eta_reg(smesh%nreg))
tMaker_td_obj%eta_reg=-1.d0
DO i=1,tMaker_td_obj%gs_eq%ncond_regs
    j=tMaker_td_obj%gs_eq%cond_regions(i)%id
    tMaker_td_obj%eta_reg(j)=tMaker_td_obj%gs_eq%cond_regions(i)%eta
END DO
ALLOCATE(tMaker_td_obj%curr_reg(smesh%nreg))
tMaker_td_obj%curr_reg=0.d0
DO i=1,tMaker_td_obj%gs_eq%ncoil_regs
    j=tMaker_td_obj%gs_eq%coil_regions(i)%id
    tMaker_td_obj%curr_reg(j)=tMaker_td_obj%gs_eq%coil_regions(i)%curr
END DO
! Advance using MF-NK method
CALL build_vac_op(tMaker_td_obj,tMaker_td_obj%vac_op)
CALL build_jop(tMaker_td_obj,adv_op,psi_sol)
adv_solver%A=>adv_op
adv_solver%its=5
adv_solver%nrits=5
adv_solver%pre=>adv_pre
adv_pre%A=>tMaker_td_obj%vac_op
!
CALL mfmat%setup(psi_tmp,tMaker_td_obj)
! CALL mfmat%utyp%delete()
! DEALLOCATE(mfmat%utyp)
mfmat%b0=1.d-5
mf_solver%A=>mfmat
mf_solver%its=100
mf_solver%nrits=20
mf_solver%atol=lin_tol
mf_solver%itplot=1
mf_solver%pm=oft_env%pm
! CALL create_diag_pre(mf_solver%pre)
! mf_solver%pre%A=>dels
! mf_solver%pre=>adv_solver
mf_solver%pre=>adv_pre
!
nksolver%A=>tMaker_td_obj
nksolver%J_inv=>mf_solver
nksolver%its=20
nksolver%atol=nl_tol
nksolver%rtol=1.d-20 ! Disable relative tolerance
nksolver%backtrack=.FALSE.
nksolver%J_update=>tMaker_td_mfnk_update
nksolver%up_freq=1
end subroutine setup_gs_td
!
subroutine step_gs_td(time,dt,nl_its,lin_its,nretry)
REAL(8), INTENT(inout) :: time,dt
INTEGER(4), INTENT(out) :: nl_its,lin_its,nretry
INTEGER(4) :: i,j,k,ierr
! Update operators if the timestep has changed
IF(dt/=tMaker_td_obj%dt)THEN
    dt=ABS(dt)
    tMaker_td_obj%dt=dt
    CALL build_vac_op(tMaker_td_obj,tMaker_td_obj%vac_op)
    CALL build_jop(tMaker_td_obj,adv_op,psi_sol)
    CALL adv_pre%update(.TRUE.)
END IF
! Update coil currents (end of time step)
DO i=1,tMaker_td_obj%gs_eq%ncoil_regs
    j=tMaker_td_obj%gs_eq%coil_regions(i)%id
    tMaker_td_obj%curr_reg(j)=tMaker_td_obj%gs_eq%coil_regions(i)%curr
END DO
! Point to profiles in case they changed
tMaker_td_obj%F=>tMaker_td_obj%gs_eq%I
tMaker_td_obj%P=>tMaker_td_obj%gs_eq%P
! Update current target and sync scale factors
tMaker_td_obj%ip_target=tMaker_td_obj%gs_eq%Itor_target
tMaker_td_obj%ip_ratio_target=tMaker_td_obj%gs_eq%Ip_ratio_target
tMaker_td_obj%f_scale=tMaker_td_obj%gs_eq%alam
tMaker_td_obj%p_scale=tMaker_td_obj%gs_eq%pnorm
!
! Advance B
!
CALL psi_tmp%add(0.d0,1.d0,psi_sol)
CALL apply_rhs(tMaker_td_obj,psi_sol,rhs1)
CALL blag_zerob(rhs1)
! ! Extrapolate solution (linear)
! DO j=maxextrap,2,-1
!   CALL extrap_fields(j)%f%add(0.d0,1.d0,extrap_fields(j-1)%f)
!   extrapt(j)=extrapt(j-1)
! END DO
! IF(nextrap<maxextrap)nextrap=nextrap+1
! CALL extrap_fields(1)%f%add(0.d0,1.d0,psi_sol)
! extrapt(1)=time_val

! !---Rebuild plasma matrix
! IF(MOD(i,2)==0)THEN
!   CALL build_jop(tMaker_td_obj,adv_op,psi_sol)
!   CALL adv_pre%update(.TRUE.)
! END IF
DO j=1,4
    ! CALL vector_extrapolate(extrapt,extrap_fields,nextrap,time_val+tMaker_td_obj%dt,psi_sol)
    !---MFNK iteration
    CALL nksolver%apply(psi_sol,rhs1)
    IF(nksolver%cits<0)THEN
        CALL psi_sol%add(0.d0,1.d0,psi_tmp)
        tMaker_td_obj%dt=tMaker_td_obj%dt/2.d0
        CALL build_vac_op(tMaker_td_obj,tMaker_td_obj%vac_op)
        CALL build_jop(tMaker_td_obj,adv_op,psi_sol)
        CALL adv_pre%update(.TRUE.)
        CALL apply_rhs(tMaker_td_obj,psi_sol,rhs1)
        CALL blag_zerob(rhs1)
        CYCLE
    ELSE
        EXIT
    END IF
END DO
!
time=time+tMaker_td_obj%dt
dt=tMaker_td_obj%dt
nl_its=nksolver%nlits
lin_its=nksolver%lits
nretry=j-1
IF(j>4)THEN
    nretry=-nretry
ELSE
    tMaker_td_obj%gs_eq%alam=tMaker_td_obj%f_scale
    tMaker_td_obj%gs_eq%pnorm=tMaker_td_obj%p_scale
END IF
! WRITE(*,'(4ES14.6,3I4)')time_val,tMaker_td_obj%ip,tMaker_td_obj%gs_eq%o_point(2),err_tmp, &
! nksolver%nlits,nksolver%lits,j
! IF(tMaker_td_obj%ip<1.d-6)THEN
!     tMaker_td_obj%p_scale=0.d0
!     tMaker_td_obj%f_scale=0.d0
!     tMaker_td_obj%ip_target=-1.d0
! END IF
end subroutine step_gs_td
!
subroutine eig_gs_td(eq_in,neigs,eigs,eig_vecs,omega,include_bounds)
TYPE(gs_eq), TARGET, INTENT(inout) :: eq_in
INTEGER(4), INTENT(in) :: neigs
REAL(r8), INTENT(out) :: eigs(:,:),eig_vecs(:,:)
REAL(r8), INTENT(in) :: omega
LOGICAL, INTENT(in) :: include_bounds
#ifdef HAVE_ARPACK
CLASS(oft_matrix), POINTER :: lhs_mat,rhs_mat
CLASS(oft_vector), POINTER :: eig_vec
TYPE(oft_iram_eigsolver) :: arsolver
TYPE(eig_wrapper), TARGET :: wrap_mat
REAL(8) :: lam0
INTEGER(4) :: i,j
!
tMaker_td_obj%gs_eq=>eq_in
tMaker_td_obj%ip_target=tMaker_td_obj%gs_eq%Itor_target
tMaker_td_obj%f_scale=tMaker_td_obj%gs_eq%alam
tMaker_td_obj%p_scale=tMaker_td_obj%gs_eq%pnorm
tMaker_td_obj%F=>tMaker_td_obj%gs_eq%I
tMaker_td_obj%P=>tMaker_td_obj%gs_eq%P
ALLOCATE(tMaker_td_obj%eta_reg(smesh%nreg))
tMaker_td_obj%eta_reg=-1.d0
DO i=1,tMaker_td_obj%gs_eq%ncond_regs
    j=tMaker_td_obj%gs_eq%cond_regions(i)%id
    tMaker_td_obj%eta_reg(j)=tMaker_td_obj%gs_eq%cond_regions(i)%eta
END DO

!---Build linearized opertors
NULLIFY(lhs_mat,rhs_mat)
CALL build_linearized(tMaker_td_obj,lhs_mat,rhs_mat,tMaker_td_obj%gs_eq%psi,omega,include_bounds)
wrap_mat%rhs_mat=>rhs_mat
wrap_mat%lhs_inv=>adv_pre
adv_pre%A=>lhs_mat

!---Setup Arnoldi eig value solver
arsolver%mode=1
arsolver%A=>wrap_mat
arsolver%tol=1.E-8_r8
arsolver%nev=neigs
arsolver%which='LM'
CALL oft_blagrange%vec_create(eig_vec)
CALL arsolver%apply(eig_vec,lam0)

DO i=1,neigs
    eigs(:,i)=1.d0/arsolver%eig_val(:,i)+omega
    eig_vecs(:,i)=arsolver%eig_vec(:,i)
END DO

!---Cleanup
CALL adv_pre%delete
CALL arsolver%delete
CALL lhs_mat%delete
CALL rhs_mat%delete
CALL eig_vec%delete
DEALLOCATE(lhs_mat,rhs_mat,eig_vec,tMaker_td_obj%eta_reg)
#else
eigs=0.d0
eig_vecs=0.d0
#endif
end subroutine eig_gs_td
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine apply_rhs(self,a,b)
class(oft_tmaker_td), intent(inout) :: self !< NL operator object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,3),elapsed_time,pt(3),eta_tmp,psi_tmp,dpsi_tmp(2),eta_source,psi_lim,psi_max
real(8) :: max_tmp,lim_tmp
real(r8), allocatable :: rop(:),gop(:,:),lop(:,:),vals_loc(:)
real(r8), pointer, dimension(:) :: pol_vals,eta_vals,rhs_vals
logical :: curved
! type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
! WRITE(*,*)'Hi'
! IF(oft_debug_print(1))THEN
!   WRITE(*,'(2X,A)')'Constructing Poloidal flux time-advance operator'
!   CALL mytimer%tick()
! END IF
!---------------------------------------------------------------------------
! Get local vector values
!---------------------------------------------------------------------------
NULLIFY(pol_vals,eta_vals,rhs_vals)
CALL a%get_local(pol_vals)
! CALL self%eta%get_local(eta_vals)
CALL b%set(0.d0)
CALL b%get_local(rhs_vals)
!---------------------------------------------------------------------------
! Operator integration
!---------------------------------------------------------------------------
!$omp parallel private(j,vals_loc,rop,gop,det,curved,goptmp,m,vol,jr,jc,pt,eta_tmp,psi_tmp,dpsi_tmp,eta_source)
allocate(j(oft_blagrange%nce),vals_loc(oft_blagrange%nce)) ! Local DOF and matrix indices
allocate(rop(oft_blagrange%nce),gop(3,oft_blagrange%nce)) ! Reconstructed gradient operator
!$omp do schedule(static,1)
do i=1,oft_blagrange%mesh%nc
    IF(smesh%reg(i)==1)CYCLE
    !---Get local to global DOF mapping
    call oft_blagrange%ncdofs(i,j)
    !---Get local reconstructed operators
    vals_loc=0.d0
    do m=1,oft_blagrange%quad%np ! Loop over quadrature points
    call oft_blagrange%mesh%jacobian(i,oft_blagrange%quad%pts(:,m),goptmp,vol)
    det=vol*oft_blagrange%quad%wts(m)
    pt=oft_blagrange%mesh%log2phys(i,oft_blagrange%quad%pts(:,m))
    psi_tmp=0.d0; dpsi_tmp=0.d0; eta_tmp=0.d0; eta_source=0.d0
    do jr=1,oft_blagrange%nce ! Loop over degrees of freedom
        call oft_blag_eval(oft_blagrange,i,jr,oft_blagrange%quad%pts(:,m),rop(jr))
        psi_tmp = psi_tmp + pol_vals(j(jr))*rop(jr)
    end do
    eta_tmp=self%eta_reg(smesh%reg(i))
    IF(eta_tmp>0.d0)THEN
        eta_source=psi_tmp*det/eta_tmp/(pt(1)+gs_epsilon)
    ELSE
        eta_source=self%dt*self%curr_reg(smesh%reg(i))*det !self%eta_reg(smesh%reg(i))*det
    END IF
    do jr=1,oft_blagrange%nce
        vals_loc(jr) = vals_loc(jr) + rop(jr)*eta_source
    end do
    end do
    do jr=1,oft_blagrange%nce
    !$omp atomic
    rhs_vals(j(jr)) = rhs_vals(j(jr)) + vals_loc(jr)
    end do
end do
deallocate(j,vals_loc,rop,gop)
!$omp end parallel
DO i=1,oft_blagrange%nbe
    rhs_vals(oft_blagrange%lbe(i))=pol_vals(oft_blagrange%lbe(i))
END DO
CALL b%restore_local(rhs_vals,add=.TRUE.)
DEALLOCATE(pol_vals,rhs_vals)
!---Report time
! IF(oft_debug_print(1))THEN
!   elapsed_time=mytimer%tock()
!   WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
! END IF
DEBUG_STACK_POP
end subroutine apply_rhs
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine picard_step(self,a,b,p_scale,f_scale,ip_target)
class(oft_tmaker_td), intent(inout) :: self !< NL operator object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
real(8), intent(in) :: p_scale
real(8), intent(inout) :: f_scale
real(8), optional, intent(in) :: ip_target
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,3),elapsed_time,pt(3),eta_tmp,psi_tmp,ip_p,ip_f
real(r8) :: dpsi_tmp(2),p_source,f_source,psi_lim,psi_max,eta_source,max_tmp,lim_tmp
real(r8), allocatable :: rop(:),gop(:,:),lop(:,:),vals_loc(:,:)
real(r8), pointer, dimension(:) :: pol_vals,rhs_vals,pvals
class(oft_vector), pointer :: ptmp
logical :: curved,in_bounds
! type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
! WRITE(*,*)'Hi'
! IF(oft_debug_print(1))THEN
!   WRITE(*,'(2X,A)')'Constructing Poloidal flux time-advance operator'
!   CALL mytimer%tick()
! END IF
!---------------------------------------------------------------------------
! Get local vector values
!---------------------------------------------------------------------------
NULLIFY(pol_vals,rhs_vals,ptmp,pvals)
CALL a%get_local(pol_vals)
!---
self%gs_eq%psi=>a
CALL gs_update_bounds(self%gs_eq)
! WRITE(*,*)'Full',self%gs_eq%plasma_bounds
self%F%plasma_bounds=self%gs_eq%plasma_bounds
self%P%plasma_bounds=self%gs_eq%plasma_bounds
CALL b%set(0.d0)
CALL b%get_local(rhs_vals)
CALL b%new(ptmp)
CALL ptmp%get_local(pvals)
!---------------------------------------------------------------------------
! Operator integration
!---------------------------------------------------------------------------
ip_p=0.d0
ip_f=0.d0
!$omp parallel private(j,vals_loc,rop,gop,det,curved,goptmp,m,vol,jr,jc,pt, &
!$omp eta_tmp,psi_tmp,dpsi_tmp,p_source,f_source,eta_source,in_bounds) &
!$omp reduction(+:ip_p) reduction(+:ip_f)
allocate(j(oft_blagrange%nce),vals_loc(oft_blagrange%nce,2)) ! Local DOF and matrix indices
allocate(rop(oft_blagrange%nce),gop(3,oft_blagrange%nce)) ! Reconstructed gradient operator
!$omp do schedule(static,1)
do i=1,oft_blagrange%mesh%nc
    IF(smesh%reg(i)/=1)CYCLE
    !---Get local to global DOF mapping
    call oft_blagrange%ncdofs(i,j)
    !---Get local reconstructed operators
    vals_loc=0.d0
    do m=1,oft_blagrange%quad%np ! Loop over quadrature points
    call oft_blagrange%mesh%jacobian(i,oft_blagrange%quad%pts(:,m),goptmp,vol)
    det=vol*oft_blagrange%quad%wts(m)
    pt=oft_blagrange%mesh%log2phys(i,oft_blagrange%quad%pts(:,m))
    psi_tmp=0.d0; dpsi_tmp=0.d0; eta_tmp=0.d0; p_source=0.d0; f_source=0.d0; eta_source=0.d0
    do jr=1,oft_blagrange%nce ! Loop over degrees of freedom
        call oft_blag_eval(oft_blagrange,i,jr,oft_blagrange%quad%pts(:,m),rop(jr))
        psi_tmp = psi_tmp + pol_vals(j(jr))*rop(jr)
    end do
    !---Compute local matrix contributions
    ! IF(smesh%reg(i)==1.AND.psi_tmp>psi_lim)THEN
    ! IF(self%allow_xpoints)THEN
        in_bounds=gs_test_bounds(self%gs_eq,pt).AND.(psi_tmp>self%gs_eq%plasma_bounds(1))
    ! ELSE
    !     in_bounds=psi_tmp>psi_lim
    ! END IF
    IF(in_bounds)THEN
        ! gs_source=self%dt*self%lam_amp*(psi_tmp-psi_lim)/(psi_max-psi_lim)/(pt(1)+gs_epsilon)
        ! p_source=p_scale*pt(1)*self%P%Fp((psi_tmp-psi_lim)/(psi_max-psi_lim))
        ! f_source=0.5d0*self%F%fp((psi_tmp-psi_lim)/(psi_max-psi_lim))/(pt(1)+gs_epsilon)
        p_source=p_scale*pt(1)*self%P%Fp(psi_tmp)
        f_source=0.5d0*self%F%fp(psi_tmp)/(pt(1)+gs_epsilon)
        ip_p=ip_p+p_source*det
        ip_f=ip_f+f_source*det
    END IF
    do jr=1,oft_blagrange%nce
        vals_loc(jr,1) = vals_loc(jr,1) - rop(jr)*self%dt*p_source*det
        vals_loc(jr,2) = vals_loc(jr,2) - rop(jr)*self%dt*f_source*det
    end do
    end do
    do jr=1,oft_blagrange%nce
    !$omp atomic
    pvals(j(jr)) = pvals(j(jr)) + vals_loc(jr,1)
    !$omp atomic
    rhs_vals(j(jr)) = rhs_vals(j(jr)) + vals_loc(jr,2)
    end do
end do
deallocate(j,vals_loc,rop,gop)
!$omp end parallel
DO i=1,oft_blagrange%nbe
    rhs_vals(oft_blagrange%lbe(i))=0.d0
    pvals(oft_blagrange%lbe(i))=0.d0
END DO
CALL b%restore_local(rhs_vals,add=.TRUE.)
CALL ptmp%restore_local(pvals,add=.TRUE.)
IF(PRESENT(ip_target))THEN
    f_scale=(ip_target-ip_p)/ip_f
!  WRITE(*,*)'Ip',ip_target,ip_p,ip_f,(ip_target-ip_p)/ip_f
ELSE
!  WRITE(*,*)'Ip',(ip_p+f_scale*ip_f)/mu0
END IF
CALL b%add(f_scale,1.d0,ptmp)
DEALLOCATE(pol_vals,rhs_vals,pvals)
CALL ptmp%delete
DEALLOCATE(ptmp)
!---Report time
! IF(oft_debug_print(1))THEN
!   elapsed_time=mytimer%tock()
!   WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
! END IF
DEBUG_STACK_POP
end subroutine picard_step
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine apply_mfop(self,a,b)
class(oft_tmaker_td), intent(inout) :: self !< NL operator object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
integer(i4) :: i,m,jr,jc,cell
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,3),elapsed_time,pt(3),alam_new,psi_tmp,diag(3),f(3)
real(r8) :: dpsi_tmp(2),p_source,f_source,psi_lim,psi_max,eta_source,max_tmp,lim_tmp
real(r8), allocatable :: rop(:),gop(:,:),lop(:,:),vals_loc(:,:)
real(r8), pointer, dimension(:) :: pol_vals,rhs_vals,alam_vals,pvals
class(oft_vector), pointer :: ptmp
type(oft_lag_brinterp) :: psi_interp
logical :: curved,in_bounds
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
! WRITE(*,*)'Hi'
! IF(oft_debug_print(1))THEN
!   WRITE(*,'(2X,A)')'Constructing Poloidal flux time-advance operator'
  CALL mytimer%tick()
! END IF
!---------------------------------------------------------------------------
! Get local vector values
!---------------------------------------------------------------------------
NULLIFY(pol_vals,rhs_vals,ptmp,pvals)
CALL a%get_local(pol_vals)
!---
self%gs_eq%psi=>a
CALL gs_update_bounds(self%gs_eq)
!
self%F%plasma_bounds=self%gs_eq%plasma_bounds
self%P%plasma_bounds=self%gs_eq%plasma_bounds
CALL b%set(0.d0)
CALL b%get_local(rhs_vals)
CALL b%new(ptmp)
ALLOCATE(alam_vals(b%n))
alam_vals=0.d0
!---------------------------------------------------------------------------
! Operator integration
!---------------------------------------------------------------------------
diag=0.d0
!$omp parallel private(j,vals_loc,rop,gop,det,curved,goptmp,m,vol,jr,jc,pt, &
!$omp psi_tmp,dpsi_tmp,p_source,f_source,eta_source,in_bounds) &
!$omp reduction(+:diag)
allocate(j(oft_blagrange%nce),vals_loc(oft_blagrange%nce,2)) ! Local DOF and matrix indices
allocate(rop(oft_blagrange%nce),gop(3,oft_blagrange%nce)) ! Reconstructed gradient operator
!$omp do schedule(static,1)
!ordered
do i=1,oft_blagrange%mesh%nc
    IF(smesh%reg(i)/=1)CYCLE
    !---Get local to global DOF mapping
    call oft_blagrange%ncdofs(i,j)
    !---Get local reconstructed operators
    vals_loc=0.d0
    do m=1,oft_blagrange%quad%np ! Loop over quadrature points
        call oft_blagrange%mesh%jacobian(i,oft_blagrange%quad%pts(:,m),goptmp,vol)
        det=vol*oft_blagrange%quad%wts(m)
        pt=oft_blagrange%mesh%log2phys(i,oft_blagrange%quad%pts(:,m))
        IF(.NOT.gs_test_bounds(self%gs_eq,pt))CYCLE
        psi_tmp=0.d0!; dpsi_tmp=0.d0
        do jr=1,oft_blagrange%nce ! Loop over degrees of freedom
            call oft_blag_eval(oft_blagrange,i,jr,oft_blagrange%quad%pts(:,m),rop(jr))
            psi_tmp = psi_tmp + pol_vals(j(jr))*rop(jr)
        end do
        IF(psi_tmp<self%gs_eq%plasma_bounds(1))CYCLE
        !---Compute local matrix contributions
        p_source=self%p_scale*pt(1)*self%P%Fp(psi_tmp)
        f_source=self%f_scale*0.5d0*self%F%fp(psi_tmp)/(pt(1)+gs_epsilon)
        diag=diag+[f_source,p_source,self%p_scale*self%P%F(psi_tmp)*pt(1)]*det
        do jr=1,oft_blagrange%nce
            vals_loc(jr,1) = vals_loc(jr,1) - rop(jr)*self%dt*p_source*det
            vals_loc(jr,2) = vals_loc(jr,2) - rop(jr)*self%dt*f_source*det
        end do
    end do !eta_source=psi_tmp*det/eta_tmp/(pt(1)+gs_epsilon)
    !!$omp ordered
    do jr=1,oft_blagrange%nce
        !$omp atomic
        rhs_vals(j(jr)) = rhs_vals(j(jr)) + vals_loc(jr,1)
        !$omp atomic
        alam_vals(j(jr)) = alam_vals(j(jr)) + vals_loc(jr,2)
    end do
    !!$omp end ordered
end do
deallocate(j,vals_loc,rop,gop)
!$omp end parallel
DO i=1,oft_blagrange%nbe
    rhs_vals(oft_blagrange%lbe(i))=0.d0
    alam_vals(oft_blagrange%lbe(i))=0.d0
END DO
IF(self%ip_target>0.d0.AND.diag(1)>1.d-8)THEN
    IF(self%ip_ratio_target>-1.d98)THEN
        f_source = self%ip_target/diag(1)/(1.d0+1.d0/self%ip_ratio_target)
        p_source = self%ip_target/diag(2)/(self%ip_ratio_target+1.d0)
        rhs_vals=rhs_vals*p_source+alam_vals*f_source
        self%f_scale=f_source*self%f_scale
        diag(1)=diag(1)*f_source
        self%p_scale=p_source*self%p_scale
        diag(2)=diag(2)*p_source
    ELSE
        f_source=(self%ip_target-diag(2))/diag(1)
        rhs_vals=rhs_vals+alam_vals*f_source
        self%f_scale=f_source*self%f_scale
        diag(1)=diag(1)*f_source
    END IF
ELSE
    rhs_vals=rhs_vals+alam_vals
END IF
CALL b%restore_local(rhs_vals,add=.TRUE.)
CALL b%new(ptmp)
CALL self%vac_op%apply(a,ptmp)
CALL b%add(1.d0,1.d0,ptmp)
CALL ptmp%delete
self%ip=(diag(1)+diag(2))/mu0
self%estored=diag(2)/mu0*3.d0/2.d0
! eta_source=b%dot(b)
! WRITE(*,*)'CHK',eta_source
! WRITE(*,*)self%ip
DEALLOCATE(pol_vals,rhs_vals,ptmp,alam_vals)
!---Report time
! IF(oft_debug_print(1))THEN
  mfop_time=mfop_time+mytimer%tock()
!   WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
! END IF
DEBUG_STACK_POP
end subroutine apply_mfop
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine apply_gs_mat(self,a,b)
class(tMaker_td_mat), intent(inout) :: self !< NL operator object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
integer(4) :: i
real(8), pointer, dimension(:) :: avals,bvals
CALL self%mat%apply(a,b)
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
! bvals=bvals+self%lim_vals*avals(self%lim_node)
DO i=1,oft_blagrange%nce
    bvals=bvals+self%lim_vals(:,i)*avals(self%lim_nodes(i))
    bvals=bvals+self%ax_vals(:,i)*avals(self%ax_nodes(i))
END DO
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
DEBUG_STACK_POP
end subroutine apply_gs_mat
! !---------------------------------------------------------------------------
! !> Needs docs
! !---------------------------------------------------------------------------
! subroutine update_lims(self,a)
! class(oft_tmaker_td), intent(inout) :: self !< NL operator object
! class(oft_vector), target, intent(inout) :: a !< Source field
! integer(i4) :: i,m,jr,jc,ilim,imax,itmp
! integer(i4), allocatable :: j(:)
! real(r8) :: vol,det,goptmp(3,3),elapsed_time,pt(3),eta_tmp,psi_tmp,ip_p,ip_f
! real(r8) :: dpsi_tmp(2),p_source,f_source,psi_lim,psi_max,eta_source,max_tmp,lim_tmp
! real(r8), allocatable :: rop(:),gop(:,:),lop(:,:),vals_loc(:)
! real(r8), pointer, dimension(:) :: pol_vals,rhs_vals,pvals
! class(oft_vector), pointer :: ptmp
! logical :: curved
! ! type(oft_timer) :: mytimer
! DEBUG_STACK_PUSH
! ! WRITE(*,*)'Hi'
! ! IF(oft_debug_print(1))THEN
! !   WRITE(*,'(2X,A)')'Constructing Poloidal flux time-advance operator'
! !   CALL mytimer%tick()
! ! END IF
! !---------------------------------------------------------------------------
! ! Get local vector values
! !---------------------------------------------------------------------------
! IF(self%allow_xpoints)THEN
!     self%gs_eq%psi=>a
!     CALL gs_update_bounds(self%gs_eq)
! ELSE
! NULLIFY(pol_vals,rhs_vals,ptmp,pvals)
! CALL a%get_local(pol_vals)
! !---Get O-point and limiter
! psi_max=-1.d99; psi_lim=1.d99
! !$omp parallel private(i,j,jr,max_tmp,lim_tmp,itmp)
! allocate(j(oft_blagrange%nce))
! max_tmp=-1.d99; lim_tmp=1.d99
! itmp=0
! !$omp do
! DO i=1,smesh%nc
!     IF(smesh%reg(i)==1)THEN
!     call oft_blagrange%ncdofs(i,j)
!     do jr=1,oft_blagrange%nce
!         IF(pol_vals(j(jr))>max_tmp)THEN
!         max_tmp=pol_vals(j(jr))
!         itmp=j(jr)
!         END IF
!     end do
!     END IF
! END DO
! !$omp critical
! IF(max_tmp>psi_max)THEN
!     psi_max=max_tmp
!     imax=itmp
! END IF
! !$omp end critical
! !$omp barrier
! !$omp do
! DO i=1,self%gs_eq%nlimiter_nds
!     IF(ABS(psi_max-pol_vals(self%gs_eq%limiter_nds(i)))<ABS(psi_max-lim_tmp))THEN
!     lim_tmp=pol_vals(self%gs_eq%limiter_nds(i))
!     itmp=self%gs_eq%limiter_nds(i)
!     END IF
! END DO
! !DO i=1,smesh%nc
! !  IF(smesh%reg(i)==8)THEN
! !    call oft_blagrange%ncdofs(i,j)
! !    do jr=1,oft_blagrange%nce
! !      IF(ABS(psi_max-pol_vals(j(jr)))<ABS(psi_max-lim_tmp))lim_tmp=pol_vals(j(jr))
! !    end do
! !  END IF
! !END DO
! !$omp critical
! IF(ABS(psi_max-lim_tmp)<ABS(psi_max-psi_lim))THEN
!     psi_lim=lim_tmp
!     ilim=itmp
! END IF
! !$omp end critical
! deallocate(j)
! !$omp end parallel
! ! WRITE(*,*)'Psi',psi_max,psi_lim
! IF(ABS(psi_max-psi_lim)<1.d-8)THEN
!     psi_max=psi_lim+1.d-8
!     ! DEALLOCATE(pol_vals)
!     ! RETURN
! END IF
! DEALLOCATE(pol_vals)
! self%gs_eq%plasma_bounds=[psi_lim,psi_max]
! END IF
! ! self%psi_max=psi_max
! ! self%psi_lim=psi_lim
! self%F%plasma_bounds=self%gs_eq%plasma_bounds ![psi_lim,psi_max]
! self%P%plasma_bounds=self%gs_eq%plasma_bounds ![psi_lim,psi_max]
! !---Report time
! ! IF(oft_debug_print(1))THEN
! !   elapsed_time=mytimer%tock()
! !   WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
! ! END IF
! DEBUG_STACK_POP
! end subroutine update_lims
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
SUBROUTINE tMaker_td_mfnk_update(a)
CLASS(oft_vector), TARGET, INTENT(inout) :: a
! CALL tMaker_td_obj%update_lims(a)
CALL mfmat%update(a)
! CALL build_jop(tMaker_td_obj,adv_op,a)
!CALL adv_solver%update(.TRUE.)
END SUBROUTINE tMaker_td_mfnk_update
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine build_jop(tMaker_td_obj,mat,a)
class(oft_tmaker_td), intent(inout) :: tMaker_td_obj
class(tMaker_td_mat), intent(inout) :: mat
class(oft_vector), target, intent(inout) :: a
integer(i4) :: i,m,jr,jc,cell
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,3),elapsed_time,pt(3),eta_tmp,eta_source,gs_source,ftmp(3)
real(r8) :: psi_lim,psi_max,psi_tmp,max_tmp,lim_tmp,psi_norm
real(r8), allocatable :: rop(:),gop(:,:),lop(:,:),lim_loc(:),ax_loc(:)
real(r8), allocatable :: lim_weights(:),ax_weights(:)
logical :: curved,in_bounds
real(r8), pointer, dimension(:) :: eta_vals,pol_vals
CLASS(oft_vector), POINTER :: oft_lag_vec
type(oft_timer) :: mytimer
!CALL build_vac_op(mat)
!RETURN
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
    WRITE(*,'(2X,A)')'Constructing Toroidal flux time-advance operator'
    CALL mytimer%tick()
END IF
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat%mat))THEN
    CALL gs_mat_create(mat%mat)
    ALLOCATE(mat%lim_nodes(oft_blagrange%nce),mat%lim_vals(a%n,oft_blagrange%nce))
    ALLOCATE(mat%ax_nodes(oft_blagrange%nce),mat%ax_vals(a%n,oft_blagrange%nce))
ELSE
    CALL mat%mat%zero
    mat%lim_vals=0.d0
    mat%ax_vals=0.d0
END IF
NULLIFY(pol_vals,eta_vals)
! CALL eta_vec%get_local(eta_vals)
CALL a%get_local(pol_vals)
!---Update plasma boundary
tMaker_td_obj%gs_eq%psi=>a
CALL gs_update_bounds(tMaker_td_obj%gs_eq)
allocate(lim_weights(oft_blagrange%nce))
cell=0
CALL bmesh_findcell(smesh,cell,tMaker_td_obj%gs_eq%lim_point,ftmp)
call oft_blagrange%ncdofs(cell,mat%lim_nodes)
do jc=1,oft_blagrange%nce ! Loop over degrees of freedom
  call oft_blag_eval(oft_blagrange,cell,jc,ftmp,lim_weights(jc))
end do
allocate(ax_weights(oft_blagrange%nce))
cell=0
CALL bmesh_findcell(smesh,cell,tMaker_td_obj%gs_eq%o_point,ftmp)
call oft_blagrange%ncdofs(cell,mat%ax_nodes)
do jc=1,oft_blagrange%nce ! Loop over degrees of freedom
  call oft_blag_eval(oft_blagrange,cell,jc,ftmp,ax_weights(jc))
end do
! !---
! tMaker_td_obj%gs_eq%psi=>a !psi_sol
! CALL gs_update_bounds(tMaker_td_obj%gs_eq)
! lim_tmp=1.d99
! DO i=1,smesh%np
!     IF(SQRT(SUM((tMaker_td_obj%gs_eq%lim_point-smesh%r(1:2,i))**2))<lim_tmp)THEN
!         mat%lim_node=i
!         lim_tmp=SQRT(SUM((tMaker_td_obj%gs_eq%lim_point-smesh%r(1:2,i))**2))
!     END IF
! END DO
! WRITE(*,*)mat%lim_node
tMaker_td_obj%F%plasma_bounds=tMaker_td_obj%gs_eq%plasma_bounds
tMaker_td_obj%P%plasma_bounds=tMaker_td_obj%gs_eq%plasma_bounds
psi_norm=tMaker_td_obj%gs_eq%plasma_bounds(2)-tMaker_td_obj%gs_eq%plasma_bounds(1)
!---------------------------------------------------------------------------
! Operator integration
!---------------------------------------------------------------------------
!$omp parallel private(j,rop,gop,det,lop,curved,goptmp,m,vol,jc,jr,pt,psi_tmp, &
!$omp in_bounds,eta_tmp,eta_source,gs_source,lim_loc,ax_loc)
allocate(j(oft_blagrange%nce)) ! Local DOF and matrix indices
allocate(rop(oft_blagrange%nce),gop(3,oft_blagrange%nce)) ! Reconstructed gradient operator
allocate(lop(oft_blagrange%nce,oft_blagrange%nce),lim_loc(oft_blagrange%nce),ax_loc(oft_blagrange%nce))
!$omp do schedule(static,1)
!ordered
do i=1,oft_blagrange%mesh%nc
    ! IF(smesh%reg(i)==1)CYCLE
    !---Get local to global DOF mapping
    call oft_blagrange%ncdofs(i,j)
    !---Get local reconstructed operators
    lop=0.d0; lim_loc=0.d0; ax_loc=0.d0
    do m=1,oft_blagrange%quad%np ! Loop over quadrature points
        call oft_blagrange%mesh%jacobian(i,oft_blagrange%quad%pts(:,m),goptmp,vol)
        det=vol*oft_blagrange%quad%wts(m)
        pt=smesh%log2phys(i,oft_blagrange%quad%pts(:,m))
        eta_tmp=0.d0; psi_tmp=0.d0
        do jc=1,oft_blagrange%nce ! Loop over degrees of freedom
            call oft_blag_eval(oft_blagrange,i,jc,oft_blagrange%quad%pts(:,m),rop(jc))
            call oft_blag_geval(oft_blagrange,i,jc,oft_blagrange%quad%pts(:,m),gop(:,jc),goptmp)
            psi_tmp=psi_tmp+pol_vals(j(jc))*rop(jc)
        end do
        eta_source=0.d0; gs_source=0.d0
        IF(smesh%reg(i)==1)THEN
            in_bounds=gs_test_bounds(tMaker_td_obj%gs_eq,pt).AND.(psi_tmp>tMaker_td_obj%gs_eq%plasma_bounds(1))
            IF(in_bounds)THEN
                gs_source=tMaker_td_obj%dt*(tMaker_td_obj%p_scale*pt(1)*pt(1)*tMaker_td_obj%P%Fpp(psi_tmp) &
                + tMaker_td_obj%f_scale*0.5d0*tMaker_td_obj%F%fpp(psi_tmp))
            END IF
        ELSE IF(smesh%reg(i)>1.AND.(tMaker_td_obj%eta_reg(smesh%reg(i))>0.d0))THEN
            eta_source=1.d0/tMaker_td_obj%eta_reg(smesh%reg(i)) !eta_tmp
        END IF
        !---Compute local matrix contributions
        do jr=1,oft_blagrange%nce
            ax_loc(jr) = ax_loc(jr) + &
              rop(jr)*gs_source*(psi_tmp-tMaker_td_obj%gs_eq%plasma_bounds(1))/psi_norm*det/(pt(1)+gs_epsilon)
            lim_loc(jr) = lim_loc(jr) + &
              rop(jr)*gs_source*(1.d0-(psi_tmp-tMaker_td_obj%gs_eq%plasma_bounds(1))/psi_norm)*det/(pt(1)+gs_epsilon)
            do jc=1,oft_blagrange%nce
            lop(jr,jc) = lop(jr,jc) + (tMaker_td_obj%dt*DOT_PRODUCT(gop(1:2,jr),gop(1:2,jc)) &
                + rop(jr)*rop(jc)*(eta_source-gs_source))*det/(pt(1)+gs_epsilon)
            end do
        end do
    end do
    !---Apply bc to local matrix
    DO jr=1,oft_blagrange%nce
        IF(oft_blagrange%be(j(jr)))THEN
            lop(jr,:)=0.d0
            lim_loc(jr)=0.d0
            ax_loc(jr)=0.d0
        END IF
    END DO
    !---Add local values to global matrix
    !!$omp ordered
    call mat%mat%atomic_add_values(j,j,lop,oft_blagrange%nce,oft_blagrange%nce)
    DO jc=1,oft_blagrange%nce
    DO jr=1,oft_blagrange%nce
        !$omp atomic
        mat%lim_vals(j(jr),jc)=mat%lim_vals(j(jr),jc)+lim_loc(jr)*lim_weights(jc)
        !$omp atomic
        mat%ax_vals(j(jr),jc)=mat%ax_vals(j(jr),jc)+ax_loc(jr)*ax_weights(jc)
    END DO
    END DO
    !!$omp end ordered
end do
deallocate(j,rop,gop,lop,lim_loc,ax_loc)
!$omp end parallel
DEALLOCATE(pol_vals,lim_weights,ax_weights)
CALL set_bcmat(tMaker_td_obj%gs_eq,mat%mat)
! !---Set diagonal entries for dirichlet rows
! ALLOCATE(lop(1,1),j(1))
! lop(1,1)=1.d0
! DO i=1,oft_blagrange%nbe
!   IF(.NOT.oft_blagrange%linkage%leo(i))CYCLE
!   j=oft_blagrange%lbe(i)
!   call mat%add_values(j,j,lop,1,1)
! END DO
! DEALLOCATE(j,lop)
!---Assemble matrix
CALL oft_blagrange%vec_create(oft_lag_vec)
CALL mat%mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec)
!---Report time
IF(oft_debug_print(1))THEN
    elapsed_time=mytimer%tock()
    WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine build_jop
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine build_vac_op(tMaker_td_obj,mat)
class(oft_tmaker_td), intent(inout) :: tMaker_td_obj
class(oft_matrix), pointer, intent(inout) :: mat
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,3),elapsed_time,pt(3),eta_tmp,eta_source,gs_source
real(r8) :: psi_lim,psi_max,psi_tmp,max_tmp,lim_tmp
real(r8), allocatable :: rop(:),gop(:,:),lop(:,:)
logical :: curved,in_bounds
real(r8), pointer, dimension(:) :: eta_vals,pol_vals
CLASS(oft_vector), POINTER :: oft_lag_vec
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
    WRITE(*,'(2X,A)')'Constructing Toroidal flux time-advance operator'
    CALL mytimer%tick()
END IF
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat))THEN
    ! CALL oft_blagrange%mat_create(mat)
    CALL gs_mat_create(mat)
ELSE
    CALL mat%zero
END IF
NULLIFY(pol_vals,eta_vals)
! CALL eta_vec%get_local(eta_vals)
! CALL psi_sol%get_local(pol_vals)
!---------------------------------------------------------------------------
! Operator integration
!---------------------------------------------------------------------------
!$omp parallel private(j,rop,gop,det,lop,curved,goptmp,m,vol,jc,jr,pt,psi_tmp, &
!$omp in_bounds,eta_tmp,eta_source,gs_source)
allocate(j(oft_blagrange%nce)) ! Local DOF and matrix indices
allocate(rop(oft_blagrange%nce),gop(3,oft_blagrange%nce)) ! Reconstructed gradient operator
allocate(lop(oft_blagrange%nce,oft_blagrange%nce)) ! Local laplacian matrix
!$omp do schedule(static,1) ordered
do i=1,oft_blagrange%mesh%nc
    ! IF(smesh%reg(i)==1)CYCLE
    !---Get local to global DOF mapping
    call oft_blagrange%ncdofs(i,j)
    !---Get local reconstructed operators
    lop=0.d0
    do m=1,oft_blagrange%quad%np ! Loop over quadrature points
    call oft_blagrange%mesh%jacobian(i,oft_blagrange%quad%pts(:,m),goptmp,vol)
    det=vol*oft_blagrange%quad%wts(m)
    pt=smesh%log2phys(i,oft_blagrange%quad%pts(:,m))
    eta_tmp=0.d0; psi_tmp=0.d0
    do jc=1,oft_blagrange%nce ! Loop over degrees of freedom
        call oft_blag_eval(oft_blagrange,i,jc,oft_blagrange%quad%pts(:,m),rop(jc))
        call oft_blag_geval(oft_blagrange,i,jc,oft_blagrange%quad%pts(:,m),gop(:,jc),goptmp)
        ! psi_tmp=psi_tmp+pol_vals(j(jc))*rop(jc)
    end do
    eta_source=0.d0; gs_source=0.d0
    IF(smesh%reg(i)>1.AND.(tMaker_td_obj%eta_reg(smesh%reg(i))>0.d0))THEN
        eta_source=1.d0/tMaker_td_obj%eta_reg(smesh%reg(i)) !eta_tmp
    END IF
    !---Compute local matrix contributions
    do jr=1,oft_blagrange%nce
        do jc=1,oft_blagrange%nce
        lop(jr,jc) = lop(jr,jc) + (tMaker_td_obj%dt*DOT_PRODUCT(gop(1:2,jr),gop(1:2,jc)) &
            + rop(jr)*rop(jc)*eta_source)*det/(pt(1)+gs_epsilon)
        
        ! lop(jr,jc) = lop(jr,jc) + DOT_PRODUCT(gop(1:2,jr),gop(1:2,jc))*det/(pt(1)+gs_epsilon)
        end do
    end do
    end do
    !---Apply bc to local matrix
    DO jr=1,oft_blagrange%nce
    IF(oft_blagrange%be(j(jr)))lop(jr,:)=0.d0
    ! IF(tMaker_td_obj%fe_flag(j(jr)))lop(jr,:)=0.d0
    END DO
    !---Add local values to global matrix
    !$omp ordered
    call mat%atomic_add_values(j,j,lop,oft_blagrange%nce,oft_blagrange%nce)
    !$omp end ordered
end do
deallocate(j,rop,gop,lop)
!$omp end parallel
! DEALLOCATE(pol_vals)
!
CALL set_bcmat(tMaker_td_obj%gs_eq,mat)
! !---Set diagonal entries for dirichlet rows
! ALLOCATE(lop(1,1),j(1))
! lop(1,1)=1.d0
! DO i=1,oft_blagrange%nbe
!   IF(.NOT.oft_blagrange%linkage%leo(i))CYCLE
!   j=oft_blagrange%lbe(i)
!   call mat%add_values(j,j,lop,1,1)
! END DO
! DEALLOCATE(j,lop)
!---Assemble matrix
CALL oft_blagrange%vec_create(oft_lag_vec)
CALL mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec)
!---Report time
IF(oft_debug_print(1))THEN
    elapsed_time=mytimer%tock()
    WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
end subroutine build_vac_op
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine build_linearized(tMaker_td_obj,lhs_mat,rhs_mat,a,sigma,include_bounds)
class(oft_tmaker_td), intent(inout) :: tMaker_td_obj
class(oft_matrix), pointer, intent(inout) :: rhs_mat,lhs_mat
class(oft_vector), target, intent(inout) :: a
real(8), intent(in) :: sigma
LOGICAL, INTENT(in) :: include_bounds
integer(i4) :: i,m,jr,jc,rowtmp(1),lim_node,ax_node,cell,nnonaxi,block_max
integer(i4), allocatable :: j(:),bnd_nodes(:),node_mark(:,:)
real(r8) :: vol,det,goptmp(3,3),elapsed_time,pt(3),eta_source,gs_source,ax_tmp,ftmp(3)
real(r8) :: psi_lim,psi_max,psi_tmp,max_tmp,lim_tmp,psi_norm,lim_source,ax_source
real(r8), allocatable :: rop(:),gop(:,:),lhs_vals(:,:),rhs_vals(:,:),lim_loc(:),ax_loc(:)
real(r8), allocatable :: lim_weights(:),ax_weights(:)
logical :: curved,in_bounds
integer(i4), allocatable, dimension(:) :: dense_flag,reg_map
real(r8), pointer, dimension(:) :: eta_vals,pol_vals,lim_vals,ax_vals
real(r8), pointer, dimension(:,:) :: nonaxi_vals
type(oft_1d_int), pointer, dimension(:) :: noaxi_nodes
CLASS(oft_vector), POINTER :: oft_lag_vec
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
    WRITE(*,'(2X,A)')'Constructing Toroidal flux time-advance operator'
    CALL mytimer%tick()
END IF
!---Update plasma boundary
tMaker_td_obj%gs_eq%psi=>a
CALL gs_update_bounds(tMaker_td_obj%gs_eq)
allocate(bnd_nodes(2*oft_blagrange%nce),lim_weights(oft_blagrange%nce),ax_weights(oft_blagrange%nce))
cell=0
CALL bmesh_findcell(smesh,cell,tMaker_td_obj%gs_eq%lim_point,ftmp)
call oft_blagrange%ncdofs(cell,bnd_nodes(1:oft_blagrange%nce))
do jc=1,oft_blagrange%nce ! Loop over degrees of freedom
  call oft_blag_eval(oft_blagrange,cell,jc,ftmp,lim_weights(jc))
end do
cell=0
CALL bmesh_findcell(smesh,cell,tMaker_td_obj%gs_eq%o_point,ftmp)
call oft_blagrange%ncdofs(cell,bnd_nodes(oft_blagrange%nce+1:2*oft_blagrange%nce))
do jc=1,oft_blagrange%nce ! Loop over degrees of freedom
  call oft_blag_eval(oft_blagrange,cell,jc,ftmp,ax_weights(jc))
end do
!---
nnonaxi=0
DO jr=1,tMaker_td_obj%gs_eq%ncond_regs
    IF(tMaker_td_obj%gs_eq%cond_regions(jr)%contiguous)CYCLE
    nnonaxi=nnonaxi+1
END DO
ALLOCATE(reg_map(smesh%nreg))
IF(nnonaxi>0)THEN
    ALLOCATE(node_mark(smesh%nreg,oft_blagrange%ne+1))
    node_mark=0
    allocate(j(oft_blagrange%nce))
    DO i=1,oft_blagrange%mesh%nc
        call oft_blagrange%ncdofs(i,j)
        DO jc=1,oft_blagrange%nce
            IF(node_mark(smesh%reg(i),j(jc))/=0)CYCLE
            node_mark(smesh%reg(i),oft_blagrange%ne+1)=node_mark(smesh%reg(i),oft_blagrange%ne+1)+1
            node_mark(smesh%reg(i),j(jc))=node_mark(smesh%reg(i),oft_blagrange%ne+1)
        END DO
    END DO
    deallocate(j)
END IF
ALLOCATE(dense_flag(oft_blagrange%ne))
dense_flag=0
WRITE(*,*)nnonaxi
ALLOCATE(noaxi_nodes(nnonaxi+1))
reg_map=0
block_max=0
m=0
DO jr=1,tMaker_td_obj%gs_eq%ncond_regs
    IF(tMaker_td_obj%gs_eq%cond_regions(jr)%contiguous)CYCLE
    i=tMaker_td_obj%gs_eq%cond_regions(jr)%id
    m=m+1
    reg_map(i)=m
    noaxi_nodes(m)%n=node_mark(i,oft_blagrange%ne+1)
    ALLOCATE(noaxi_nodes(m)%v(noaxi_nodes(m)%n))
    noaxi_nodes(m)%v=0
    DO jc=1,oft_blagrange%ne
        IF(node_mark(i,jc)>0)noaxi_nodes(m)%v(node_mark(i,jc)) = jc
    END DO
    dense_flag(noaxi_nodes(m)%v)=m
    block_max=MAX(block_max,noaxi_nodes(m)%n)
END DO
m=m+1
noaxi_nodes(m)%n=oft_blagrange%nbe
noaxi_nodes(m)%v=>oft_blagrange%lbe
dense_flag(noaxi_nodes(m)%v)=m
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(lhs_mat))THEN
    ! CALL oft_blagrange%mat_create(rhs_mat)
    CALL lin_mat_create(2*oft_blagrange%nce,bnd_nodes,lhs_mat,dense_flag,noaxi_nodes)
    ALLOCATE(lim_vals(a%n),ax_vals(a%n))
    lim_vals=0.d0; ax_vals=0.d0
    ! Remove boundary points
    WHERE(dense_flag==m)
        dense_flag=0
    END WHERE
    CALL lin_mat_create(0,bnd_nodes,rhs_mat,dense_flag,noaxi_nodes)
    DEALLOCATE(dense_flag)
ELSE
    CALL oft_abort("Matrices should not be allocated","build_linearized",__FILE__)
    ! CALL rhs_mat%zero
    ! CALL lhs_mat%mat%zero
    ! lhs_mat%lim_vals=0.d0
END IF
NULLIFY(pol_vals,eta_vals)
! CALL eta_vec%get_local(eta_vals)
CALL a%get_local(pol_vals)
! WRITE(*,*)mat%lim_node
tMaker_td_obj%F%plasma_bounds=tMaker_td_obj%gs_eq%plasma_bounds
tMaker_td_obj%P%plasma_bounds=tMaker_td_obj%gs_eq%plasma_bounds
psi_norm=tMaker_td_obj%gs_eq%plasma_bounds(2)-tMaker_td_obj%gs_eq%plasma_bounds(1)
!---------------------------------------------------------------------------
! Operator integration
!---------------------------------------------------------------------------
IF(nnonaxi>0)THEN
    ALLOCATE(nonaxi_vals(block_max+1,nnonaxi))
    nonaxi_vals=0.d0
END IF
!$omp parallel private(j,rop,gop,det,lhs_vals,rhs_vals,curved,goptmp,m,vol,jc,jr,pt,psi_tmp, &
!$omp in_bounds,eta_source,gs_source,lim_loc,lim_source,ax_source,ax_loc)
allocate(j(oft_blagrange%nce)) ! Local DOF and matrix indices
allocate(rop(oft_blagrange%nce),gop(3,oft_blagrange%nce)) ! Reconstructed gradient operator
allocate(lhs_vals(oft_blagrange%nce,oft_blagrange%nce),lim_loc(oft_blagrange%nce))
allocate(rhs_vals(oft_blagrange%nce,oft_blagrange%nce),ax_loc(oft_blagrange%nce))
!$omp do schedule(static,1)
!ordered
do i=1,oft_blagrange%mesh%nc
    ! IF(smesh%reg(i)==1)CYCLE
    !---Get local to global DOF mapping
    call oft_blagrange%ncdofs(i,j)
    !---Get local reconstructed operators
    lhs_vals=0.d0; lim_loc=0.d0; ax_loc=0.d0; rhs_vals=0.d0
    do m=1,oft_blagrange%quad%np ! Loop over quadrature points
        call oft_blagrange%mesh%jacobian(i,oft_blagrange%quad%pts(:,m),goptmp,vol)
        det=vol*oft_blagrange%quad%wts(m)
        pt=smesh%log2phys(i,oft_blagrange%quad%pts(:,m))
        psi_tmp=0.d0
        do jc=1,oft_blagrange%nce ! Loop over degrees of freedom
            call oft_blag_eval(oft_blagrange,i,jc,oft_blagrange%quad%pts(:,m),rop(jc))
            call oft_blag_geval(oft_blagrange,i,jc,oft_blagrange%quad%pts(:,m),gop(:,jc),goptmp)
            psi_tmp=psi_tmp+pol_vals(j(jc))*rop(jc)
        end do
        eta_source=0.d0; gs_source=0.d0; lim_source=0.d0; ax_source=0.d0
        IF(smesh%reg(i)==1)THEN
            in_bounds=gs_test_bounds(tMaker_td_obj%gs_eq,pt).AND.(psi_tmp>tMaker_td_obj%gs_eq%plasma_bounds(1))
            IF(in_bounds)THEN
                gs_source=(tMaker_td_obj%p_scale*pt(1)*pt(1)*tMaker_td_obj%P%Fpp(psi_tmp) &
                    + tMaker_td_obj%f_scale*0.5d0*tMaker_td_obj%F%fpp(psi_tmp))
                lim_source=gs_source*(1.d0-(psi_tmp-tMaker_td_obj%gs_eq%plasma_bounds(1))/psi_norm)
                ax_source=gs_source*(psi_tmp-tMaker_td_obj%gs_eq%plasma_bounds(1))/psi_norm
            END IF
        ELSE IF(smesh%reg(i)>1.AND.(tMaker_td_obj%eta_reg(smesh%reg(i))>0.d0))THEN
            eta_source=1.d0/tMaker_td_obj%eta_reg(smesh%reg(i))
        END IF
        !---Compute local matrix contributions
        do jr=1,oft_blagrange%nce
            lim_loc(jr) = lim_loc(jr) + rop(jr)*lim_source*det/(pt(1)+gs_epsilon)
            ax_loc(jr) = ax_loc(jr) + rop(jr)*ax_source*det/(pt(1)+gs_epsilon)
            do jc=1,oft_blagrange%nce
                lhs_vals(jr,jc) = lhs_vals(jr,jc) + (DOT_PRODUCT(gop(1:2,jr),gop(1:2,jc)) &
                    - rop(jr)*rop(jc)*gs_source)*det/(pt(1)+gs_epsilon)
                rhs_vals(jr,jc) = rhs_vals(jr,jc) + rop(jr)*rop(jc)*eta_source*det/(pt(1)+gs_epsilon)
            end do
            IF(reg_map(smesh%reg(i))>0)THEN
                !$omp atomic
                nonaxi_vals(node_mark(smesh%reg(i),j(jr)),reg_map(smesh%reg(i))) = &
                  nonaxi_vals(node_mark(smesh%reg(i),j(jr)),reg_map(smesh%reg(i))) &
                  + rop(jr)*eta_source*det/(pt(1)+gs_epsilon)
            END IF
        end do
        IF(reg_map(smesh%reg(i))>0)THEN
            !$omp atomic
            nonaxi_vals(block_max+1,reg_map(smesh%reg(i))) = &
              nonaxi_vals(block_max+1,reg_map(smesh%reg(i))) + det
        END IF
    end do
    !---Apply bc to local matrix
    DO jr=1,oft_blagrange%nce
        IF(oft_blagrange%be(j(jr)))THEN
            rhs_vals(jr,:)=0.d0
            rhs_vals(:,jr)=0.d0
            lhs_vals(jr,:)=0.d0
            ! lhs_vals(:,jr)=0.d0
            lim_loc(jr)=0.d0
            ax_loc(jr)=0.d0
        END IF
    END DO
    lhs_vals=lhs_vals-sigma*rhs_vals
    !---Add local values to global matrix
    !!$omp ordered
    call rhs_mat%atomic_add_values(j,j,rhs_vals,oft_blagrange%nce,oft_blagrange%nce)
    call lhs_mat%atomic_add_values(j,j,lhs_vals,oft_blagrange%nce,oft_blagrange%nce)
    DO jr=1,oft_blagrange%nce
        !$omp atomic
        lim_vals(j(jr))=lim_vals(j(jr))+lim_loc(jr)
        !$omp atomic
        ax_vals(j(jr))=ax_vals(j(jr))+ax_loc(jr)
    END DO
    !!$omp end ordered
end do
IF(nnonaxi>0)THEN
    !$omp do schedule(static,1)
    !ordered
    do i=1,oft_blagrange%mesh%nc
        IF(reg_map(smesh%reg(i))==0)CYCLE
        !---Get local to global DOF mapping
        call oft_blagrange%ncdofs(i,j)
        !---Get local reconstructed operators
        lim_loc=0.d0
        do m=1,oft_blagrange%quad%np ! Loop over quadrature points
            call oft_blagrange%mesh%jacobian(i,oft_blagrange%quad%pts(:,m),goptmp,vol)
            det=vol*oft_blagrange%quad%wts(m)
            do jc=1,oft_blagrange%nce ! Loop over degrees of freedom
                call oft_blag_eval(oft_blagrange,i,jc,oft_blagrange%quad%pts(:,m),rop(jc))
            end do
            !---Compute local matrix contributions
            do jr=1,oft_blagrange%nce
                lim_loc(jr) = lim_loc(jr) + rop(jr)*det
            end do
        end do
        lim_loc=-lim_loc/nonaxi_vals(block_max+1,reg_map(smesh%reg(i)))
        !---Add local values to global matrix
        m=reg_map(smesh%reg(i))
        !!$omp ordered
        do jc=1,oft_blagrange%nce
            call rhs_mat%atomic_add_values(j(jc:jc),noaxi_nodes(m)%v, &
              lim_loc(jc)*nonaxi_vals(:,m),1,noaxi_nodes(m)%n)
            call lhs_mat%atomic_add_values(j(jc:jc),noaxi_nodes(m)%v, &
              -sigma*lim_loc(jc)*nonaxi_vals(:,m),1,noaxi_nodes(m)%n)
        end do
        !!$omp end ordered
    end do
END IF
deallocate(j,rop,gop,lhs_vals,lim_loc,ax_loc,rhs_vals)
!$omp end parallel
DEALLOCATE(pol_vals)
IF(nnonaxi>0)THEN
    DO i=1,nnonaxi
        DEALLOCATE(noaxi_nodes(i)%v)
    END DO
    DEALLOCATE(reg_map,noaxi_nodes,node_mark)
    DEALLOCATE(nonaxi_vals)
END IF
IF(include_bounds)THEN
  ALLOCATE(j(1),lhs_vals(1,1))
  DO jr=1,oft_blagrange%nce
    j=bnd_nodes(jr)
    DO i=1,oft_blagrange%ne
        lhs_vals=lim_vals(i)*lim_weights(jr)
        rowtmp=i
        CALL lhs_mat%add_values(rowtmp,j,lhs_vals,1,1)
    END DO
    j=bnd_nodes(oft_blagrange%nce+jr)
    DO i=1,oft_blagrange%ne
        lhs_vals=ax_vals(i)*ax_weights(jr)
        rowtmp=i
        CALL lhs_mat%add_values(rowtmp,j,lhs_vals,1,1)
    END DO
  END DO
  DEALLOCATE(j,lhs_vals)
END IF
DEALLOCATE(bnd_nodes,lim_weights,ax_weights)
DEALLOCATE(lim_vals,ax_vals)
CALL set_bcmat(tMaker_td_obj%gs_eq,lhs_mat)
!---Assemble matrix
CALL oft_blagrange%vec_create(oft_lag_vec)
CALL rhs_mat%assemble(oft_lag_vec)
CALL lhs_mat%assemble(oft_lag_vec)
CALL oft_lag_vec%delete
DEALLOCATE(oft_lag_vec)
!---Report time
IF(oft_debug_print(1))THEN
    elapsed_time=mytimer%tock()
    WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
END IF
DEBUG_STACK_POP
contains
!
subroutine lin_mat_create(nadd,nodes_add,new,dense_flag,dense_nodes)
INTEGER(4), INTENT(in) :: nadd,nodes_add(nadd),dense_flag(:)
CLASS(oft_matrix), POINTER, INTENT(out) :: new
type(oft_1d_int), intent(in) :: dense_nodes(:)
CLASS(oft_vector), POINTER :: tmp_vec
INTEGER(4) :: i,j,k,nr,iadd,offset
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: ltmp
TYPE(oft_graph_ptr) :: graphs(1,1)
DEBUG_STACK_PUSH
CALL oft_blagrange%vec_create(tmp_vec)
ALLOCATE(graphs(1,1)%g)
graphs(1,1)%g%nr=oft_blagrange%ne
graphs(1,1)%g%nrg=oft_blagrange%global%ne
graphs(1,1)%g%nc=oft_blagrange%ne
graphs(1,1)%g%ncg=oft_blagrange%global%ne
graphs(1,1)%g%nnz=0
ALLOCATE(graphs(1,1)%g%kr(oft_blagrange%ne+1))
graphs(1,1)%g%kr=0
ALLOCATE(ltmp(oft_blagrange%ne))
DO i=1,oft_blagrange%ne
    IF(dense_flag(i)/=0)THEN
        ltmp=2*oft_blagrange%ne
        offset=dense_nodes(dense_flag(i))%n
        ltmp(1:offset)=dense_nodes(dense_flag(i))%v
        DO j=oft_blagrange%kee(i),oft_blagrange%kee(i+1)-1
            ltmp(j-oft_blagrange%kee(i)+offset+1)=oft_blagrange%lee(j)
        END DO
        DO iadd=1,nadd
            ltmp(oft_blagrange%kee(i+1)-oft_blagrange%kee(i)+offset+iadd)=nodes_add(iadd)
        END DO
        CALL sort_array(ltmp,offset+nadd+oft_blagrange%kee(i+1)-oft_blagrange%kee(i))
        graphs(1,1)%g%kr(i)=1
        DO j=2,offset+nadd+oft_blagrange%kee(i+1)-oft_blagrange%kee(i)
            IF(ltmp(j)>ltmp(j-1))graphs(1,1)%g%kr(i)=graphs(1,1)%g%kr(i)+1
        END DO
    ELSE
        graphs(1,1)%g%kr(i)=oft_blagrange%kee(i+1)-oft_blagrange%kee(i)
        DO iadd=1,nadd
            j=search_array(nodes_add(iadd),oft_blagrange%lee(oft_blagrange%kee(i):oft_blagrange%kee(i+1)-1), &
                oft_blagrange%kee(i+1)-oft_blagrange%kee(i))
            IF(j==0)graphs(1,1)%g%kr(i)=graphs(1,1)%g%kr(i)+1
        END DO
    END IF
END DO
graphs(1,1)%g%nnz=SUM(graphs(1,1)%g%kr)
graphs(1,1)%g%kr(oft_blagrange%ne+1)=graphs(1,1)%g%nnz+1
DO i=oft_blagrange%ne,1,-1
    graphs(1,1)%g%kr(i)=graphs(1,1)%g%kr(i+1)-graphs(1,1)%g%kr(i)
END DO
IF(graphs(1,1)%g%kr(1)/=1)CALL oft_abort('Bad new graph setup','gs_mat_create', &
__FILE__)
ALLOCATE(graphs(1,1)%g%lc(graphs(1,1)%g%nnz))
DO i=1,oft_blagrange%ne
    IF(dense_flag(i)/=0)THEN
        ltmp=2*oft_blagrange%ne
        offset=dense_nodes(dense_flag(i))%n
        ltmp(1:offset)=dense_nodes(dense_flag(i))%v
        DO j=oft_blagrange%kee(i),oft_blagrange%kee(i+1)-1
            ltmp(j-oft_blagrange%kee(i)+offset+1)=oft_blagrange%lee(j)
        END DO
        DO iadd=1,nadd
            ltmp(oft_blagrange%kee(i+1)-oft_blagrange%kee(i)+offset+iadd)=nodes_add(iadd)
        END DO
        CALL sort_array(ltmp,offset+nadd+oft_blagrange%kee(i+1)-oft_blagrange%kee(i))
        nr=0
        graphs(1,1)%g%lc(graphs(1,1)%g%kr(i))=ltmp(1)
        DO j=2,offset+nadd+oft_blagrange%kee(i+1)-oft_blagrange%kee(i)
            IF(ltmp(j)>ltmp(j-1))THEN
            nr=nr+1
            graphs(1,1)%g%lc(graphs(1,1)%g%kr(i)+nr)=ltmp(j)
            END IF
        END DO
    ELSE
        DO j=oft_blagrange%kee(i),oft_blagrange%kee(i+1)-1
            graphs(1,1)%g%lc(graphs(1,1)%g%kr(i)+j-oft_blagrange%kee(i))=oft_blagrange%lee(j)
        END DO
        k=0
        DO iadd=1,nadd
            j=search_array(nodes_add(iadd),oft_blagrange%lee(oft_blagrange%kee(i):oft_blagrange%kee(i+1)-1), &
                oft_blagrange%kee(i+1)-oft_blagrange%kee(i))
            IF(j==0)THEN
                k=k+1
                graphs(1,1)%g%lc(graphs(1,1)%g%kr(i+1)-k)=nodes_add(iadd)
            END IF
        END DO
        IF((graphs(1,1)%g%kr(i+1)-graphs(1,1)%g%kr(i))>(oft_blagrange%kee(i+1)-oft_blagrange%kee(i)))THEN
            CALL sort_array(graphs(1,1)%g%lc(graphs(1,1)%g%kr(i):graphs(1,1)%g%kr(i+1)-1), &
                graphs(1,1)%g%kr(i+1)-graphs(1,1)%g%kr(i))
        END IF
    END IF
END DO
!---
CALL create_matrix(new,graphs,tmp_vec,tmp_vec)
CALL tmp_vec%delete
DEALLOCATE(ltmp)
DEALLOCATE(graphs(1,1)%g,tmp_vec)
DEBUG_STACK_POP
end subroutine lin_mat_create
end subroutine build_linearized
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine apply_wrap(self,a,b)
class(eig_wrapper), intent(inout) :: self !< NL operator object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
class(oft_vector), pointer :: tmp_vec
CALL oft_blagrange%vec_create(tmp_vec)
CALL self%rhs_mat%apply(a,tmp_vec)
CALL self%lhs_inv%apply(b,tmp_vec)
CALL tmp_vec%delete()
DEALLOCATE(tmp_vec)
DEBUG_STACK_POP
end subroutine apply_wrap
end module oft_gs_td