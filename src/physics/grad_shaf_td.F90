!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
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
USE oft_mesh_type, ONLY: oft_bmesh, bmesh_findcell
USE oft_mesh_local_util, ONLY: mesh_local_findedge
USE oft_quadrature, ONLY: oft_quad_type
USE oft_gauss_quadrature, ONLY: set_quad_1d
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_graph, oft_graph_ptr, oft_vector_ptr
USE oft_la_utils, ONLY: create_matrix, graph_add_dense_blocks, graph_add_full_col
USE oft_solver_base, ONLY: oft_solver
USE oft_deriv_matrices, ONLY: oft_noop_matrix, oft_mf_matrix
USE oft_native_la, ONLY: oft_native_matrix
USE oft_native_solvers, ONLY: oft_nksolver, oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE oft_lu, ONLY: oft_lusolver
#ifdef HAVE_ARPACK
USE oft_arpack, ONLY: oft_iram_eigsolver
#endif
USE fem_utils, ONLY: bfem_map_flag
USE oft_lag_basis, ONLY: oft_blag_geval, oft_blag_eval, oft_blag_npos, &
  oft_scalar_bfem
USE oft_blag_operators, ONLY: oft_lag_brinterp
USE axi_green, ONLY: green, grad_green
USE oft_gs, ONLY: gs_epsilon, flux_func, gs_eq, gs_update_bounds, &
    gs_test_bounds, gs_mat_create, compute_bcmat, set_bcmat, build_dels
USE mhd_utils, ONLY: mu0
IMPLICIT NONE
#include "local.h"
integer(i4), parameter :: maxextrap = 2
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(oft_noop_matrix) :: oft_tmaker_td_mfop
    real(r8) :: dt = 1.E-3_r8 !< Time step size [s]
    real(r8) :: f_scale = 1.d0 !< Scale factor for \f$ F*F' \f$ term
    real(r8) :: p_scale = 1.d0 !< Scale factor for \f$ P' \f$ term
    real(r8) :: ip = 0.d0 !< Plasma current at present step
    real(r8) :: estored = 0.d0 !< Stored energy at present step
    real(r8) :: ip_target = -1.d0 !< Target plasma current
    real(r8) :: ip_ratio_target = -1.d0 !< Ip ratio target
    real(8), pointer, dimension(:) :: eta_reg => NULL() !< Resistivity by region
    real(8), pointer, dimension(:) :: curr_reg => NULL() !< Coil current by region
    CLASS(flux_func), POINTER :: F => NULL() !< Flux function for \f$ F*F' \f$ term
    CLASS(flux_func), POINTER :: P => NULL() !< Flux function for \f$ P' \f$ term
    TYPE(gs_eq), POINTER :: gs_eq => NULL() !< Equilibrium object
    CLASS(oft_matrix), POINTER :: vac_op => NULL() !< Vacuum time-advance operator
contains
    !> Setup operator, allocating internal storage
    procedure :: setup => setup_mfop
    !> Update operator with new targets, etc. 
    procedure :: update => update_mfop
    !> Delete operator, deallocating internal storage
    procedure :: delete => delete_mfop
    !> Apply operator
    procedure :: apply_real => apply_mfop
    !
    ! procedure :: update_lims => update_lims
end type oft_tmaker_td_mfop
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(oft_noop_matrix) :: tMaker_td_mat
    integer(4), pointer, dimension(:) :: lim_nodes => NULL() !< List of nodes defining limiter flux
    integer(4), pointer, dimension(:) :: ax_nodes => NULL() !< List of nodes defining O-point flux
    real(8), pointer, dimension(:,:) :: lim_vals => NULL() !< Basis function values for limiter nodes
    real(8), pointer, dimension(:,:) :: ax_vals => NULL() !< Basis function values for O-point nodes
    class(oft_matrix), pointer :: mat => NULL() !< Matrix storage
contains
    !> Apply the operator
    procedure :: apply_real => apply_gs_mat
    !> Delete operator, deallocating internal storage
    procedure :: delete => delete_gs_mat
end type tMaker_td_mat
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(oft_noop_matrix) :: eig_wrapper
    class(oft_matrix), pointer :: rhs_mat => NULL() !< Needs Docs
    class(oft_solver), pointer :: lhs_inv => NULL() !< Needs Docs
contains
    !> Apply the operator
    procedure :: apply_real => apply_wrap
end type eig_wrapper
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type :: oft_tmaker_td
    CLASS(oft_vector), POINTER :: rhs => NULL() !< Temporary RHS vector
    CLASS(oft_vector), POINTER :: psi_sol => NULL() !< Current solution vector
    CLASS(oft_vector), POINTER :: psi_tmp => NULL() !< Temporary storage vector
    CLASS(oft_vector), POINTER :: tmp_vec => NULL() !< Temporary storage vector
    TYPE(oft_tmaker_td_mfop), POINTER :: mfop => NULL() !< Time-advance operator
    TYPE(oft_mf_matrix), POINTER :: mfmat => NULL() !< Matrix free Jacobian operator
    TYPE(tMaker_td_mat), POINTER :: adv_op => NULL() !< Preconditioner operator
    type(oft_native_gmres_solver), POINTER :: mf_solver => NULL() !< Outer linear solver
    type(oft_native_gmres_solver), POINTER :: pre_solver => NULL() !< Full preconditioner
    TYPE(oft_lusolver), POINTER :: vac_pre => NULL() !< Preconditioner using vacuum operator
    type(oft_nksolver) :: nksolver !< Newton-Krylov solver for time-advance
    integer(i4) :: nextrap !< Needs Docs
    real(r8), allocatable, dimension(:) :: extrapt !< Needs Docs
    type(oft_vector_ptr), allocatable, dimension(:) :: extrap_fields !< Needs Docs
    real(8) :: mfop_time !< Needs Docs
contains
    !> Apply the matrix
    procedure :: setup => setup_gs_td
    !> Needs Docs
    procedure :: delete => delete_gs_td
    !> Needs Docs
    procedure :: step => step_gs_td
end type oft_tmaker_td
TYPE(oft_tmaker_td), POINTER :: active_tMaker_td => NULL()
!$omp threadprivate(active_tMaker_td)
CONTAINS
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine setup_gs_td(self,eq_in,dt,lin_tol,nl_tol,pre_plasma)
class(oft_tmaker_td), intent(inout) :: self !< NL operator object
TYPE(gs_eq), TARGET, INTENT(inout) :: eq_in !< Needs Docs
REAL(8), INTENT(in) :: dt !< Needs Docs
REAL(8), INTENT(in) :: lin_tol !< Needs Docs
REAL(8), INTENT(in) :: nl_tol !< Needs Docs
LOGICAL, INTENT(in) :: pre_plasma !< Needs Docs
INTEGER(4) :: i,j,k,ierr,io_unit,npts,iostat
LOGICAL :: file_exists
LOGICAL, ALLOCATABLE, DIMENSION(:) :: vel_bc
REAL(8) :: err_tmp,ip_target,p_scale,f_scale,time_val,pt(2)
REAL(8), ALLOCATABLE, DIMENSION(:) :: np_count,res_in
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: pts
REAL(8), POINTER, DIMENSION(:) :: vals_out,vals2,rhs_tmp
!---
IF(ASSOCIATED(self%mfop))CALL self%delete()
!---------------------------------------------------------------------------
! Create operator
!---------------------------------------------------------------------------
ALLOCATE(self%mfop)
self%mfop%dt=dt
CALL self%mfop%setup(eq_in)
!---------------------------------------------------------------------------
! Create Solver fields
!---------------------------------------------------------------------------
NULLIFY(vals_out)
self%psi_sol=>self%mfop%gs_eq%psi
call eq_in%fe_rep%vec_create(self%rhs)
call eq_in%fe_rep%vec_create(self%psi_tmp)
call eq_in%fe_rep%vec_create(self%tmp_vec)
!---------------------------------------------------------------------------
! Create extrapolation fields (Unused)
!---------------------------------------------------------------------------
IF(maxextrap>0)THEN
    ALLOCATE(self%extrap_fields(maxextrap),self%extrapt(maxextrap))
    DO i=1,maxextrap
        CALL eq_in%fe_rep%vec_create(self%extrap_fields(i)%f)
        self%extrapt(i)=0.d0
    END DO
    self%nextrap=0
END IF
!
ALLOCATE(self%vac_pre)
self%vac_pre%A=>self%mfop%vac_op
!
ALLOCATE(self%mfmat)
CALL self%mfmat%setup(self%psi_tmp,self%mfop)
! CALL self%mfmat%utyp%delete()
! DEALLOCATE(self%mfmat%utyp)
ALLOCATE(self%mf_solver)
self%mfmat%b0=1.d-5
self%mf_solver%A=>self%mfmat
self%mf_solver%its=100
self%mf_solver%nrits=20
self%mf_solver%atol=lin_tol
self%mf_solver%itplot=1
self%mf_solver%pm=oft_env%pm
! Setup preconditioner
! CALL create_diag_pre(self%mf_solver%pre)
! self%mf_solver%pre%A=>dels
IF(pre_plasma.AND.(self%mfop%gs_eq%region_info%nnonaxi<=0))THEN
    ALLOCATE(self%adv_op)
    CALL build_jop(self%mfop,self%adv_op,self%psi_sol)
    ALLOCATE(self%pre_solver)
    self%pre_solver%A=>self%adv_op
    self%pre_solver%its=5
    self%pre_solver%nrits=5
    self%pre_solver%pre=>self%vac_pre
    !
    self%mf_solver%pre=>self%pre_solver
ELSE
    self%mf_solver%pre=>self%vac_pre
END IF
!
self%nksolver%A=>self%mfop
self%nksolver%J_inv=>self%mf_solver
self%nksolver%its=20
self%nksolver%atol=nl_tol
self%nksolver%rtol=1.d-20 ! Disable relative tolerance
self%nksolver%backtrack=.FALSE.
self%nksolver%J_update=>tMaker_td_mfnk_update
self%nksolver%up_freq=1
end subroutine setup_gs_td
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine delete_gs_td(self)
class(oft_tmaker_td), intent(inout) :: self !< NL operator object
INTEGER(4) :: i
DEBUG_STACK_PUSH
!
self%mfop_time=0.d0
!
IF(ASSOCIATED(self%mfop))THEN
    CALL self%mfop%delete()
    DEALLOCATE(self%mfop)
END IF
!
IF(ASSOCIATED(self%rhs))THEN
    CALL self%rhs%delete()
    CALL self%psi_tmp%delete()
    CALL self%tmp_vec%delete()
    DEALLOCATE(self%rhs,self%psi_tmp,self%tmp_vec)
    NULLIFY(self%psi_sol)
    !
    IF(ALLOCATED(self%extrap_fields))THEN
        DO i=1,SIZE(self%extrap_fields)
            IF(.NOT.ASSOCIATED(self%extrap_fields(i)%f))CYCLE
            CALL self%extrap_fields(i)%f%delete()
            DEALLOCATE(self%extrap_fields(i)%f)
        END DO
        DEALLOCATE(self%extrap_fields,self%extrapt)
    END IF
    !
    CALL self%mfmat%delete()
    DEALLOCATE(self%mfmat)
    IF(ASSOCIATED(self%adv_op))THEN
        CALL self%adv_op%delete()
        DEALLOCATE(self%adv_op)
    END IF
    !
    CALL self%vac_pre%delete()
    CALL self%mf_solver%delete()
    DEALLOCATE(self%mf_solver,self%vac_pre)
    IF(ASSOCIATED(self%pre_solver))THEN
        CALL self%pre_solver%delete()
        DEALLOCATE(self%pre_solver)
    END IF
    !
    CALL self%nksolver%delete()
END IF
DEBUG_STACK_POP
end subroutine delete_gs_td
!
subroutine step_gs_td(self,time,dt,nl_its,lin_its,nretry)
class(oft_tmaker_td), target, intent(inout) :: self !< NL operator object
REAL(8), INTENT(inout) :: time,dt
INTEGER(4), INTENT(out) :: nl_its,lin_its,nretry
INTEGER(4) :: i,j,k,ierr
active_tMaker_td=>self
! Update time-advance operator
CALL self%mfop%update()
! Update operators if the timestep has changed
IF(dt/=self%mfop%dt)THEN
    dt=ABS(dt)
    self%mfop%dt=dt
    CALL build_vac_op(self%mfop,self%mfop%vac_op)
    IF(ASSOCIATED(self%adv_op))CALL build_jop(self%mfop,self%adv_op,self%psi_sol)
    CALL self%vac_pre%update(.TRUE.)
END IF
!
! Advance B
!
CALL self%psi_tmp%add(0.d0,1.d0,self%psi_sol)
CALL apply_rhs(self%mfop,self%psi_sol,self%rhs)
CALL self%mfop%gs_eq%zerob_bc%apply(self%rhs)
! ! Extrapolate solution (linear)
! DO j=maxextrap,2,-1
!   CALL extrap_fields(j)%f%add(0.d0,1.d0,extrap_fields(j-1)%f)
!   extrapt(j)=extrapt(j-1)
! END DO
! IF(nextrap<maxextrap)nextrap=nextrap+1
! CALL extrap_fields(1)%f%add(0.d0,1.d0,psi_sol)
! extrapt(1)=time_val
DO j=1,4
    ! CALL vector_extrapolate(extrapt,extrap_fields,nextrap,time_val+self%dt,psi_sol)
    !---MFNK iteration
    CALL self%nksolver%apply(self%psi_sol,self%rhs)
    IF(self%nksolver%cits<0)THEN
        CALL self%psi_sol%add(0.d0,1.d0,self%psi_tmp)
        self%mfop%dt=self%mfop%dt/2.d0
        CALL build_vac_op(self%mfop,self%mfop%vac_op)
        IF(ASSOCIATED(self%adv_op))CALL build_jop(self%mfop,self%adv_op,self%psi_sol)
        CALL self%vac_pre%update(.TRUE.)
        CALL apply_rhs(self%mfop,self%psi_sol,self%rhs)
        CALL self%mfop%gs_eq%zerob_bc%apply(self%rhs)
        CYCLE
    ELSE
        EXIT
    END IF
END DO
!
time=time+self%mfop%dt
dt=self%mfop%dt
nl_its=self%nksolver%nlits
lin_its=self%nksolver%lits
nretry=j-1
IF(j>4)THEN
    nretry=-nretry
ELSE
    self%mfop%gs_eq%alam=self%mfop%f_scale
    self%mfop%gs_eq%pnorm=self%mfop%p_scale
END IF
! WRITE(*,'(4ES14.6,3I4)')time_val,self%mfop%ip,self%mfop%gs_eq%o_point(2),err_tmp, &
! nksolver%nlits,nksolver%lits,j
! IF(self%mfop%ip<1.d-6)THEN
!     self%mfop%p_scale=0.d0
!     self%mfop%f_scale=0.d0
!     self%mfop%ip_target=-1.d0
! END IF
end subroutine step_gs_td
!
subroutine eig_gs_td(eq_in,neigs,eigs,eig_vecs,omega,include_bounds,eta_plasma)
TYPE(gs_eq), TARGET, INTENT(inout) :: eq_in
INTEGER(4), INTENT(in) :: neigs
REAL(r8), INTENT(out) :: eigs(:,:),eig_vecs(:,:)
REAL(r8), INTENT(in) :: omega
LOGICAL, INTENT(in) :: include_bounds
REAL(r8), INTENT(in) :: eta_plasma
#ifdef HAVE_ARPACK
CLASS(oft_matrix), POINTER :: lhs_mat,rhs_mat
CLASS(oft_vector), POINTER :: eig_vec
TYPE(oft_lusolver), TARGET :: adv_pre
TYPE(oft_iram_eigsolver) :: arsolver
TYPE(eig_wrapper), TARGET :: wrap_mat
INTEGER(4) :: i,j
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: sort_tmp
REAL(8) :: lam0,sigma_plasma
REAL(8), ALLOCATABLE, DIMENSION(:) :: eig_tmp
type(oft_tmaker_td_mfop) :: eig_mfop !< NL operator object
!
CALL eig_mfop%setup(eq_in)

!---Build linearized opertors
NULLIFY(lhs_mat,rhs_mat)
sigma_plasma=0.d0
IF(eta_plasma>0.d0)sigma_plasma=1.d0/eta_plasma
CALL build_linearized(eig_mfop,lhs_mat,rhs_mat,eig_mfop%gs_eq%psi,omega,include_bounds,sigma_plasma)
wrap_mat%rhs_mat=>rhs_mat
wrap_mat%lhs_inv=>adv_pre
adv_pre%A=>lhs_mat

!---Setup Arnoldi eig value solver
arsolver%mode=1
arsolver%A=>wrap_mat
arsolver%tol=1.E-8_r8
arsolver%nev=neigs
arsolver%which='LM'
CALL eq_in%fe_rep%vec_create(eig_vec)
CALL arsolver%apply(eig_vec,lam0)

IF(arsolver%info>=0)THEN
  !---Sort eigenvalues
  ALLOCATE(sort_tmp(arsolver%nev),eig_tmp(arsolver%nev))
  eig_tmp=1.d0/arsolver%eig_val(1,1:arsolver%nev)+omega
  sort_tmp=[(i,i=1,arsolver%nev)]
  CALL sort_array(eig_tmp,sort_tmp,arsolver%nev)
  DEALLOCATE(eig_tmp)
  !---Copy output
  DO i=1,neigs
    j = sort_tmp(i)
    eigs(:,i)=1.d0/arsolver%eig_val(:,j)+omega
    eig_vecs(:,i)=arsolver%eig_vec(:,j)
  END DO
  DEALLOCATE(sort_tmp)
ELSE
    eigs=-1.d99
    eig_vecs=0.d0
END IF

!---Cleanup
CALL adv_pre%delete
CALL arsolver%delete
CALL lhs_mat%delete
CALL rhs_mat%delete
CALL eig_vec%delete
DEALLOCATE(lhs_mat,rhs_mat,eig_vec,eig_mfop%eta_reg)
CALL eig_mfop%delete()
#else
eigs=0.d0
eig_vecs=0.d0
#endif
end subroutine eig_gs_td
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine apply_rhs(self,a,b)
class(oft_tmaker_td_mfop), intent(inout) :: self !< NL operator object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
integer(i4) :: i,m,jr,jc
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,3),elapsed_time,pt(3),eta_tmp,psi_tmp,eta_source,psi_lim,psi_max
real(8) :: max_tmp,lim_tmp
real(r8), allocatable :: rop(:),gop(:,:),lop(:,:),vals_loc(:),reg_source(:)
real(r8), pointer, dimension(:) :: pol_vals,rhs_vals
logical :: curved
CLASS(oft_bmesh), POINTER :: mesh
CLASS(oft_scalar_bfem), POINTER :: lag_rep
DEBUG_STACK_PUSH
mesh=>self%gs_eq%mesh
lag_rep=>self%gs_eq%fe_rep
!---------------------------------------------------------------------------
! Get local vector values
!---------------------------------------------------------------------------
NULLIFY(pol_vals,rhs_vals)
CALL a%get_local(pol_vals)
CALL b%set(0.d0)
CALL b%get_local(rhs_vals)
!
ALLOCATE(reg_source(mesh%nreg))
reg_source=0.d0
IF(ASSOCIATED(self%gs_eq%region_info%nonaxi_vals))THEN
    DO i=1,mesh%nreg
        IF(self%gs_eq%region_info%reg_map(i)==0)CYCLE
        reg_source(i)=DOT_PRODUCT(pol_vals,self%gs_eq%region_info%nonaxi_vals(:,i))
    END DO
END IF
!---------------------------------------------------------------------------
! Operator integration
!---------------------------------------------------------------------------
!$omp parallel private(j,vals_loc,rop,gop,det,curved,goptmp,m,vol,jr,jc,pt,eta_tmp,psi_tmp,eta_source)
allocate(j(lag_rep%nce),vals_loc(lag_rep%nce)) ! Local DOF and matrix indices
allocate(rop(lag_rep%nce),gop(3,lag_rep%nce)) ! Reconstructed gradient operator
!$omp do schedule(static,1)
do i=1,mesh%nc
    IF(mesh%reg(i)==1)CYCLE
    !---Get local to global DOF mapping
    call lag_rep%ncdofs(i,j)
    !---Get local reconstructed operators
    vals_loc=0.d0
    do m=1,lag_rep%quad%np ! Loop over quadrature points
        call mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
        det=vol*lag_rep%quad%wts(m)
        pt=mesh%log2phys(i,lag_rep%quad%pts(:,m))
        psi_tmp=0.d0
        do jr=1,lag_rep%nce ! Loop over degrees of freedom
            call oft_blag_eval(lag_rep,i,jr,lag_rep%quad%pts(:,m),rop(jr))
            psi_tmp = psi_tmp + pol_vals(j(jr))*rop(jr)
        end do
        eta_source=0.d0
        eta_tmp=self%eta_reg(mesh%reg(i))
        IF(eta_tmp>0.d0)THEN
            eta_source=(psi_tmp/eta_tmp/(pt(1)+gs_epsilon) + reg_source(mesh%reg(i)))*det
        ELSE
            eta_source=self%dt*self%curr_reg(mesh%reg(i))*det
        END IF
        do jr=1,lag_rep%nce
            vals_loc(jr) = vals_loc(jr) + rop(jr)*eta_source
        end do
    end do
    do jr=1,lag_rep%nce
        !$omp atomic
        rhs_vals(j(jr)) = rhs_vals(j(jr)) + vals_loc(jr)
    end do
end do
deallocate(j,vals_loc,rop,gop)
!$omp end parallel
DO i=1,lag_rep%nbe
    rhs_vals(lag_rep%lbe(i))=pol_vals(lag_rep%lbe(i))
END DO
CALL b%restore_local(rhs_vals,add=.TRUE.)
DEALLOCATE(pol_vals,rhs_vals,reg_source)
DEBUG_STACK_POP
end subroutine apply_rhs
! !---------------------------------------------------------------------------
! !> Needs docs
! !---------------------------------------------------------------------------
! subroutine picard_step(self,a,b,p_scale,f_scale,ip_target)
! class(oft_tmaker_td_mfop), intent(inout) :: self !< NL operator object
! class(oft_vector), target, intent(inout) :: a !< Source field
! class(oft_vector), intent(inout) :: b !< Result of metric function
! real(8), intent(in) :: p_scale
! real(8), intent(inout) :: f_scale
! real(8), optional, intent(in) :: ip_target
! integer(i4) :: i,m,jr,jc
! integer(i4), allocatable :: j(:)
! real(r8) :: vol,det,goptmp(3,3),elapsed_time,pt(3),eta_tmp,psi_tmp,ip_p,ip_f
! real(r8) :: dpsi_tmp(2),p_source,f_source,psi_lim,psi_max,eta_source,max_tmp,lim_tmp
! real(r8), allocatable :: rop(:),gop(:,:),lop(:,:),vals_loc(:,:)
! real(r8), pointer, dimension(:) :: pol_vals,rhs_vals,pvals
! class(oft_vector), pointer :: ptmp
! logical :: curved,in_bounds
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
! NULLIFY(pol_vals,rhs_vals,ptmp,pvals)
! CALL a%get_local(pol_vals)
! !---
! self%gs_eq%psi=>a ! HERE
! CALL gs_update_bounds(self%gs_eq)
! ! WRITE(*,*)'Full',self%gs_eq%plasma_bounds
! self%F%plasma_bounds=self%gs_eq%plasma_bounds
! self%P%plasma_bounds=self%gs_eq%plasma_bounds
! CALL b%set(0.d0)
! CALL b%get_local(rhs_vals)
! CALL b%new(ptmp)
! CALL ptmp%get_local(pvals)
! !---------------------------------------------------------------------------
! ! Operator integration
! !---------------------------------------------------------------------------
! ip_p=0.d0
! ip_f=0.d0
! !$omp parallel private(j,vals_loc,rop,gop,det,curved,goptmp,m,vol,jr,jc,pt, &
! !$omp eta_tmp,psi_tmp,dpsi_tmp,p_source,f_source,eta_source,in_bounds) &
! !$omp reduction(+:ip_p) reduction(+:ip_f)
! allocate(j(oft_blagrange%nce),vals_loc(oft_blagrange%nce,2)) ! Local DOF and matrix indices
! allocate(rop(oft_blagrange%nce),gop(3,oft_blagrange%nce)) ! Reconstructed gradient operator
! !$omp do schedule(static,1)
! do i=1,oft_blagrange%mesh%nc
!     IF(oft_blagrange%mesh%reg(i)/=1)CYCLE
!     !---Get local to global DOF mapping
!     call oft_blagrange%ncdofs(i,j)
!     !---Get local reconstructed operators
!     vals_loc=0.d0
!     do m=1,oft_blagrange%quad%np ! Loop over quadrature points
!     call oft_blagrange%mesh%jacobian(i,oft_blagrange%quad%pts(:,m),goptmp,vol)
!     det=vol*oft_blagrange%quad%wts(m)
!     pt=oft_blagrange%mesh%log2phys(i,oft_blagrange%quad%pts(:,m))
!     psi_tmp=0.d0; dpsi_tmp=0.d0; eta_tmp=0.d0; p_source=0.d0; f_source=0.d0; eta_source=0.d0
!     do jr=1,oft_blagrange%nce ! Loop over degrees of freedom
!         call oft_blag_eval(oft_blagrange,i,jr,oft_blagrange%quad%pts(:,m),rop(jr))
!         psi_tmp = psi_tmp + pol_vals(j(jr))*rop(jr)
!     end do
!     !---Compute local matrix contributions
!     ! IF(oft_blagrange%mesh%reg(i)==1.AND.psi_tmp>psi_lim)THEN
!     ! IF(self%allow_xpoints)THEN
!         in_bounds=gs_test_bounds(self%gs_eq,pt).AND.(psi_tmp>self%gs_eq%plasma_bounds(1))
!     ! ELSE
!     !     in_bounds=psi_tmp>psi_lim
!     ! END IF
!     IF(in_bounds)THEN
!         ! gs_source=self%dt*self%lam_amp*(psi_tmp-psi_lim)/(psi_max-psi_lim)/(pt(1)+gs_epsilon)
!         ! p_source=p_scale*pt(1)*self%P%Fp((psi_tmp-psi_lim)/(psi_max-psi_lim))
!         ! f_source=0.5d0*self%F%fp((psi_tmp-psi_lim)/(psi_max-psi_lim))/(pt(1)+gs_epsilon)
!         p_source=p_scale*pt(1)*self%P%Fp(psi_tmp)
!         f_source=0.5d0*self%F%fp(psi_tmp)/(pt(1)+gs_epsilon)
!         ip_p=ip_p+p_source*det
!         ip_f=ip_f+f_source*det
!     END IF
!     do jr=1,oft_blagrange%nce
!         vals_loc(jr,1) = vals_loc(jr,1) - rop(jr)*self%dt*p_source*det
!         vals_loc(jr,2) = vals_loc(jr,2) - rop(jr)*self%dt*f_source*det
!     end do
!     end do
!     do jr=1,oft_blagrange%nce
!     !$omp atomic
!     pvals(j(jr)) = pvals(j(jr)) + vals_loc(jr,1)
!     !$omp atomic
!     rhs_vals(j(jr)) = rhs_vals(j(jr)) + vals_loc(jr,2)
!     end do
! end do
! deallocate(j,vals_loc,rop,gop)
! !$omp end parallel
! DO i=1,oft_blagrange%nbe
!     rhs_vals(oft_blagrange%lbe(i))=0.d0
!     pvals(oft_blagrange%lbe(i))=0.d0
! END DO
! CALL b%restore_local(rhs_vals,add=.TRUE.)
! CALL ptmp%restore_local(pvals,add=.TRUE.)
! IF(PRESENT(ip_target))THEN
!     f_scale=(ip_target-ip_p)/ip_f
! !  WRITE(*,*)'Ip',ip_target,ip_p,ip_f,(ip_target-ip_p)/ip_f
! ELSE
! !  WRITE(*,*)'Ip',(ip_p+f_scale*ip_f)/mu0
! END IF
! CALL b%add(f_scale,1.d0,ptmp)
! DEALLOCATE(pol_vals,rhs_vals,pvals)
! CALL ptmp%delete
! DEALLOCATE(ptmp)
! !---Report time
! ! IF(oft_debug_print(1))THEN
! !   elapsed_time=mytimer%tock()
! !   WRITE(*,'(4X,A,ES11.4)')'Assembly time = ',elapsed_time
! ! END IF
! DEBUG_STACK_POP
! end subroutine picard_step
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine setup_mfop(self,eq_in)
class(oft_tmaker_td_mfop), intent(inout) :: self !< NL operator object
TYPE(gs_eq), TARGET, INTENT(inout) :: eq_in
INTEGER(4) :: i,j,k
DEBUG_STACK_PUSH
self%gs_eq=>eq_in
self%ip_target=self%gs_eq%Itor_target
self%ip_ratio_target=self%gs_eq%Ip_ratio_target
self%f_scale=self%gs_eq%alam
self%p_scale=self%gs_eq%pnorm
!
self%F=>self%gs_eq%I
self%P=>self%gs_eq%P
!
ALLOCATE(self%eta_reg(eq_in%mesh%nreg))
self%eta_reg=-1.d0
DO i=1,self%gs_eq%ncond_regs
    j=self%gs_eq%cond_regions(i)%id
    self%eta_reg(j)=self%gs_eq%cond_regions(i)%eta
END DO
ALLOCATE(self%curr_reg(eq_in%mesh%nreg))
self%curr_reg=0.d0
DO i=1,self%gs_eq%ncoils
    DO k=1,self%gs_eq%ncoil_regs
        j=self%gs_eq%coil_regions(k)%id
        self%curr_reg(j)=self%curr_reg(j) &
          + (self%gs_eq%coil_currs(i) + self%gs_eq%vcontrol_val*self%gs_eq%coil_vcont(i))*self%gs_eq%coil_nturns(j,i)
    END DO
END DO
!
CALL build_vac_op(self,self%vac_op)
DEBUG_STACK_POP
end subroutine setup_mfop
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine update_mfop(self)
class(oft_tmaker_td_mfop), intent(inout) :: self !< NL operator object
INTEGER(4) :: i,j,k
DEBUG_STACK_PUSH
! Update coil currents (end of time step)
self%curr_reg=0.d0
DO i=1,self%gs_eq%ncoils
    DO k=1,self%gs_eq%ncoil_regs
        j=self%gs_eq%coil_regions(k)%id
        self%curr_reg(j)=self%curr_reg(j) &
            + (self%gs_eq%coil_currs(i) + self%gs_eq%vcontrol_val*self%gs_eq%coil_vcont(i))*self%gs_eq%coil_nturns(j,i)
    END DO
END DO
! Point to profiles in case they changed
self%F=>self%gs_eq%I
self%P=>self%gs_eq%P
! Update current target and sync scale factors
self%ip_target=self%gs_eq%Itor_target
self%ip_ratio_target=self%gs_eq%Ip_ratio_target
self%f_scale=self%gs_eq%alam
self%p_scale=self%gs_eq%pnorm
DEBUG_STACK_POP
end subroutine update_mfop
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine delete_mfop(self)
class(oft_tmaker_td_mfop), intent(inout) :: self !< NL operator object
DEBUG_STACK_PUSH
!
self%dt=-1.d0
self%ip=-1.d0
self%estored=-1.d0
self%ip_target=-1.d0
self%ip_ratio_target=-1.d0
!
IF(ASSOCIATED(self%eta_reg))DEALLOCATE(self%eta_reg)
IF(ASSOCIATED(self%curr_reg))DEALLOCATE(self%curr_reg)
! TODO: Deallocate in the future, need to check if same as GS_EQ
NULLIFY(self%F,self%P)
NULLIFY(self%gs_eq)
!
IF(ASSOCIATED(self%vac_op))THEN
    CALL self%vac_op%delete()
    DEALLOCATE(self%vac_op)
END IF
DEBUG_STACK_POP
end subroutine delete_mfop
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine apply_mfop(self,a,b)
class(oft_tmaker_td_mfop), intent(inout) :: self !< NL operator object
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
CLASS(oft_bmesh), POINTER :: mesh
CLASS(oft_scalar_bfem), POINTER :: lag_rep
DEBUG_STACK_PUSH
mesh=>self%gs_eq%mesh
lag_rep=>self%gs_eq%fe_rep
CALL mytimer%tick()
!---------------------------------------------------------------------------
! Get local vector values
!---------------------------------------------------------------------------
NULLIFY(pol_vals,rhs_vals,ptmp,pvals)
CALL a%get_local(pol_vals)
!---
self%gs_eq%psi=>a ! HERE
CALL gs_update_bounds(self%gs_eq,track_opoint=.TRUE.)
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
allocate(j(lag_rep%nce),vals_loc(lag_rep%nce,2)) ! Local DOF and matrix indices
allocate(rop(lag_rep%nce),gop(3,lag_rep%nce)) ! Reconstructed gradient operator
!$omp do schedule(static,1)
!ordered
do i=1,mesh%nc
    IF(mesh%reg(i)/=1)CYCLE
    !---Get local to global DOF mapping
    call lag_rep%ncdofs(i,j)
    !---Get local reconstructed operators
    vals_loc=0.d0
    do m=1,lag_rep%quad%np ! Loop over quadrature points
        call mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
        det=vol*lag_rep%quad%wts(m)
        pt=mesh%log2phys(i,lag_rep%quad%pts(:,m))
        IF(.NOT.gs_test_bounds(self%gs_eq,pt))CYCLE
        psi_tmp=0.d0!; dpsi_tmp=0.d0
        do jr=1,lag_rep%nce ! Loop over degrees of freedom
            call oft_blag_eval(lag_rep,i,jr,lag_rep%quad%pts(:,m),rop(jr))
            psi_tmp = psi_tmp + pol_vals(j(jr))*rop(jr)
        end do
        IF(psi_tmp<self%gs_eq%plasma_bounds(1))CYCLE
        !---Compute local matrix contributions
        p_source=self%p_scale*pt(1)*self%P%Fp(psi_tmp)
        f_source=self%f_scale*0.5d0*self%F%fp(psi_tmp)/(pt(1)+gs_epsilon)
        diag=diag+[f_source,p_source,self%p_scale*self%P%F(psi_tmp)*pt(1)]*det
        do jr=1,lag_rep%nce
            vals_loc(jr,1) = vals_loc(jr,1) - rop(jr)*self%dt*p_source*det
            vals_loc(jr,2) = vals_loc(jr,2) - rop(jr)*self%dt*f_source*det
        end do
    end do !eta_source=psi_tmp*det/eta_tmp/(pt(1)+gs_epsilon)
    !!$omp ordered
    do jr=1,lag_rep%nce
        !$omp atomic
        rhs_vals(j(jr)) = rhs_vals(j(jr)) + vals_loc(jr,1)
        !$omp atomic
        alam_vals(j(jr)) = alam_vals(j(jr)) + vals_loc(jr,2)
    end do
    !!$omp end ordered
end do
deallocate(j,vals_loc,rop,gop)
!$omp end parallel
DO i=1,lag_rep%nbe
    rhs_vals(lag_rep%lbe(i))=0.d0
    alam_vals(lag_rep%lbe(i))=0.d0
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

! self%mfop_time=self%mfop_time+mytimer%tock()
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
integer(4) :: i,n
real(8), pointer, dimension(:) :: avals,bvals
CALL self%mat%apply(a,b)
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%get_local(bvals)
! bvals=bvals+self%lim_vals*avals(self%lim_node)
n=SIZE(self%lim_vals,DIM=2)
DO i=1,n
    bvals=bvals+self%lim_vals(:,i)*avals(self%lim_nodes(i))
    bvals=bvals+self%ax_vals(:,i)*avals(self%ax_nodes(i))
END DO
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
DEBUG_STACK_POP
end subroutine apply_gs_mat
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine delete_gs_mat(self)
class(tMaker_td_mat), intent(inout) :: self !< NL operator object
IF(ASSOCIATED(self%mat))THEN
    CALL self%delete()
    DEALLOCATE(self%mat)
END IF
IF(ASSOCIATED(self%lim_nodes))DEALLOCATE(self%lim_nodes)
IF(ASSOCIATED(self%lim_vals))DEALLOCATE(self%lim_vals)
IF(ASSOCIATED(self%ax_nodes))DEALLOCATE(self%ax_nodes)
IF(ASSOCIATED(self%ax_vals))DEALLOCATE(self%ax_vals)
end subroutine delete_gs_mat
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
! DO i=1,oft_blagrange%mesh%nc
!     IF(oft_blagrange%mesh%reg(i)==1)THEN
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
! !DO i=1,oft_blagrange%mesh%nc
! !  IF(oft_blagrange%mesh%reg(i)==8)THEN
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
! CALL active_tMaker_td%mfop%update_lims(a)
CALL active_tMaker_td%mfmat%update(a)
! CALL build_jop(active_tMaker_td%mfop,adv_op,a)
!CALL active_tMaker_td%adv_solver%update(.TRUE.)
END SUBROUTINE tMaker_td_mfnk_update
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine build_jop(self,mat,a)
class(oft_tmaker_td_mfop), intent(inout) :: self
class(tMaker_td_mat), intent(inout) :: mat
class(oft_vector), target, intent(inout) :: a
integer(i4) :: i,m,jr,jc,cell
integer(i4), allocatable :: j(:)
real(r8) :: vol,det,goptmp(3,3),elapsed_time,pt(3),eta_tmp,eta_source,gs_source,ftmp(3)
real(r8) :: psi_lim,psi_max,psi_tmp,max_tmp,lim_tmp,psi_norm
real(r8), allocatable :: rop(:),gop(:,:),lop(:,:),lim_loc(:),ax_loc(:)
real(r8), allocatable :: lim_weights(:),ax_weights(:)
logical :: curved,in_bounds
real(r8), pointer, dimension(:) :: pol_vals
CLASS(oft_vector), POINTER :: oft_lag_vec
type(oft_timer) :: mytimer
CLASS(oft_bmesh), POINTER :: mesh
CLASS(oft_scalar_bfem), POINTER :: lag_rep
DEBUG_STACK_PUSH
!CALL build_vac_op(mat)
!RETURN
mesh=>self%gs_eq%mesh
lag_rep=>self%gs_eq%fe_rep
IF(oft_debug_print(1))THEN
    WRITE(*,'(2X,A)')'Constructing Toroidal flux time-advance operator'
    CALL mytimer%tick()
END IF
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(mat%mat))THEN
    CALL gs_mat_create(self%gs_eq%fe_rep,mat%mat)
    ALLOCATE(mat%lim_nodes(lag_rep%nce),mat%lim_vals(a%n,lag_rep%nce))
    ALLOCATE(mat%ax_nodes(lag_rep%nce),mat%ax_vals(a%n,lag_rep%nce))
ELSE
    CALL mat%mat%zero
    mat%lim_vals=0.d0
    mat%ax_vals=0.d0
END IF
NULLIFY(pol_vals)
CALL a%get_local(pol_vals)
!---Update plasma boundary
self%gs_eq%psi=>a
CALL gs_update_bounds(self%gs_eq,track_opoint=.TRUE.)
allocate(lim_weights(lag_rep%nce))
cell=0
CALL bmesh_findcell(mesh,cell,self%gs_eq%lim_point,ftmp)
call lag_rep%ncdofs(cell,mat%lim_nodes)
do jc=1,lag_rep%nce ! Loop over degrees of freedom
  call oft_blag_eval(lag_rep,cell,jc,ftmp,lim_weights(jc))
end do
allocate(ax_weights(lag_rep%nce))
cell=0
CALL bmesh_findcell(mesh,cell,self%gs_eq%o_point,ftmp)
call lag_rep%ncdofs(cell,mat%ax_nodes)
do jc=1,lag_rep%nce ! Loop over degrees of freedom
  call oft_blag_eval(lag_rep,cell,jc,ftmp,ax_weights(jc))
end do
! !---
! self%gs_eq%psi=>a !psi_sol
! CALL gs_update_bounds(self%gs_eq)
! lim_tmp=1.d99
! DO i=1,mesh%np
!     IF(SQRT(SUM((self%gs_eq%lim_point-mesh%r(1:2,i))**2))<lim_tmp)THEN
!         mat%lim_node=i
!         lim_tmp=SQRT(SUM((self%gs_eq%lim_point-mesh%r(1:2,i))**2))
!     END IF
! END DO
! WRITE(*,*)mat%lim_node
self%F%plasma_bounds=self%gs_eq%plasma_bounds
self%P%plasma_bounds=self%gs_eq%plasma_bounds
psi_norm=self%gs_eq%plasma_bounds(2)-self%gs_eq%plasma_bounds(1)
!---------------------------------------------------------------------------
! Operator integration
!---------------------------------------------------------------------------
!$omp parallel private(j,rop,gop,det,lop,curved,goptmp,m,vol,jc,jr,pt,psi_tmp, &
!$omp in_bounds,eta_tmp,eta_source,gs_source,lim_loc,ax_loc)
allocate(j(lag_rep%nce)) ! Local DOF and matrix indices
allocate(rop(lag_rep%nce),gop(3,lag_rep%nce)) ! Reconstructed gradient operator
allocate(lop(lag_rep%nce,lag_rep%nce),lim_loc(lag_rep%nce),ax_loc(lag_rep%nce))
!$omp do schedule(static,1)
!ordered
do i=1,mesh%nc
    ! IF(mesh%reg(i)==1)CYCLE
    !---Get local to global DOF mapping
    call lag_rep%ncdofs(i,j)
    !---Get local reconstructed operators
    lop=0.d0; lim_loc=0.d0; ax_loc=0.d0
    do m=1,lag_rep%quad%np ! Loop over quadrature points
        call mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
        det=vol*lag_rep%quad%wts(m)
        pt=mesh%log2phys(i,lag_rep%quad%pts(:,m))
        eta_tmp=0.d0; psi_tmp=0.d0
        do jc=1,lag_rep%nce ! Loop over degrees of freedom
            call oft_blag_eval(lag_rep,i,jc,lag_rep%quad%pts(:,m),rop(jc))
            call oft_blag_geval(lag_rep,i,jc,lag_rep%quad%pts(:,m),gop(:,jc),goptmp)
            psi_tmp=psi_tmp+pol_vals(j(jc))*rop(jc)
        end do
        eta_source=0.d0; gs_source=0.d0
        IF(mesh%reg(i)==1)THEN
            in_bounds=gs_test_bounds(self%gs_eq,pt).AND.(psi_tmp>self%gs_eq%plasma_bounds(1))
            IF(in_bounds)THEN
                gs_source=self%dt*(self%p_scale*pt(1)*pt(1)*self%P%Fpp(psi_tmp) &
                + self%f_scale*0.5d0*self%F%fpp(psi_tmp))
            END IF
        ELSE IF(mesh%reg(i)>1.AND.(self%eta_reg(mesh%reg(i))>0.d0))THEN
            eta_source=1.d0/self%eta_reg(mesh%reg(i)) !eta_tmp
        END IF
        !---Compute local matrix contributions
        do jr=1,lag_rep%nce
            ax_loc(jr) = ax_loc(jr) + &
              rop(jr)*gs_source*(psi_tmp-self%gs_eq%plasma_bounds(1))/psi_norm*det/(pt(1)+gs_epsilon)
            lim_loc(jr) = lim_loc(jr) + &
              rop(jr)*gs_source*(1.d0-(psi_tmp-self%gs_eq%plasma_bounds(1))/psi_norm)*det/(pt(1)+gs_epsilon)
            do jc=1,lag_rep%nce
            lop(jr,jc) = lop(jr,jc) + (self%dt*DOT_PRODUCT(gop(1:2,jr),gop(1:2,jc)) &
                + rop(jr)*rop(jc)*(eta_source-gs_source))*det/(pt(1)+gs_epsilon)
            end do
        end do
    end do
    !---Apply bc to local matrix
    DO jr=1,lag_rep%nce
        IF(lag_rep%be(j(jr)))THEN
            lop(jr,:)=0.d0
            lim_loc(jr)=0.d0
            ax_loc(jr)=0.d0
        END IF
    END DO
    !---Add local values to global matrix
    !!$omp ordered
    call mat%mat%atomic_add_values(j,j,lop,lag_rep%nce,lag_rep%nce)
    DO jc=1,lag_rep%nce
    DO jr=1,lag_rep%nce
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
CALL set_bcmat(self%gs_eq,mat%mat)
! !---Set diagonal entries for dirichlet rows
! ALLOCATE(lop(1,1),j(1))
! lop(1,1)=1.d0
! DO i=1,lag_rep%nbe
!   IF(.NOT.lag_rep%linkage%leo(i))CYCLE
!   j=lag_rep%lbe(i)
!   call mat%add_values(j,j,lop,1,1)
! END DO
! DEALLOCATE(j,lop)
!---Assemble matrix
CALL lag_rep%vec_create(oft_lag_vec)
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
subroutine build_vac_op(self,mat)
class(oft_tmaker_td_mfop), intent(inout) :: self
class(oft_matrix), pointer, intent(inout) :: mat
CALL build_dels(mat,self%gs_eq,'free',self%dt,self%dt)
end subroutine build_vac_op
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine build_linearized(self,lhs_mat,rhs_mat,a,sigma,include_bounds,sigma_plasma)
class(oft_tmaker_td_mfop), intent(inout) :: self
class(oft_matrix), pointer, intent(inout) :: rhs_mat,lhs_mat
class(oft_vector), target, intent(inout) :: a
real(8), intent(in) :: sigma
LOGICAL, INTENT(in) :: include_bounds
real(8), intent(in) :: sigma_plasma
integer(i4) :: i,m,jr,jc,rowtmp(1),lim_node,ax_node,cell,nnonaxi
integer(i4), allocatable :: j(:),bnd_nodes(:)
real(r8) :: vol,det,goptmp(3,3),elapsed_time,pt(3),eta_source,gs_source,ax_tmp,ftmp(3)
real(r8) :: psi_lim,psi_max,psi_tmp,max_tmp,lim_tmp,psi_norm,lim_source,ax_source
real(r8), allocatable :: rop(:),gop(:,:),lhs_vals(:,:),rhs_vals(:,:),lim_loc(:),ax_loc(:)
real(r8), allocatable :: lim_weights(:),ax_weights(:)
logical :: curved,in_bounds
integer(i4), allocatable, dimension(:) :: dense_flag
real(r8), pointer, dimension(:) :: eta_vals,pol_vals,lim_vals,ax_vals,nonaxi_tmp
real(r8), pointer, dimension(:,:) :: nonaxi_vals
type(oft_1d_int), pointer, dimension(:) :: bc_nodes
CLASS(oft_vector), POINTER :: oft_lag_vec
TYPE(oft_graph_ptr) :: graphs(1,1)
TYPE(oft_graph), TARGET :: graph1,graph2
type(oft_timer) :: mytimer
CLASS(oft_bmesh), POINTER :: smesh
CLASS(oft_scalar_bfem), POINTER :: lag_rep
DEBUG_STACK_PUSH
smesh=>self%gs_eq%mesh
lag_rep=>self%gs_eq%fe_rep
IF(oft_debug_print(1))THEN
    WRITE(*,'(2X,A)')'Constructing Toroidal flux time-advance operator'
    CALL mytimer%tick()
END IF
!---Update plasma boundary
self%gs_eq%psi=>a
CALL gs_update_bounds(self%gs_eq,track_opoint=.TRUE.)
allocate(bnd_nodes(2*lag_rep%nce),lim_weights(lag_rep%nce),ax_weights(lag_rep%nce))
IF(include_bounds)THEN
    cell=0
    CALL bmesh_findcell(smesh,cell,self%gs_eq%lim_point,ftmp)
    call lag_rep%ncdofs(cell,bnd_nodes(1:lag_rep%nce))
    do jc=1,lag_rep%nce ! Loop over degrees of freedom
    call oft_blag_eval(lag_rep,cell,jc,ftmp,lim_weights(jc))
    end do
    cell=0
    CALL bmesh_findcell(smesh,cell,self%gs_eq%o_point,ftmp)
    call lag_rep%ncdofs(cell,bnd_nodes(lag_rep%nce+1:2*lag_rep%nce))
    do jc=1,lag_rep%nce ! Loop over degrees of freedom
    call oft_blag_eval(lag_rep,cell,jc,ftmp,ax_weights(jc))
    end do
END IF
!---
nnonaxi=self%gs_eq%region_info%nnonaxi
!---------------------------------------------------------------------------
! Allocate matrix
!---------------------------------------------------------------------------
IF(.NOT.ASSOCIATED(lhs_mat))THEN
    CALL lag_rep%vec_create(oft_lag_vec)
    !---
    graph1%nr=lag_rep%ne
    graph1%nrg=lag_rep%global%ne
    graph1%nc=lag_rep%ne
    graph1%ncg=lag_rep%global%ne
    graph1%nnz=lag_rep%nee
    graph1%kr=>lag_rep%kee
    graph1%lc=>lag_rep%lee
    !---Add dense blocks for non-contiguous regions
    IF(nnonaxi>0)THEN
        !---Add dense blocks
        ALLOCATE(dense_flag(lag_rep%ne))
        dense_flag=0
        DO m=1,self%gs_eq%region_info%nnonaxi
            dense_flag(self%gs_eq%region_info%noaxi_nodes(m)%v)=m
        END DO
        CALL graph_add_dense_blocks(graph1,graph2,dense_flag,self%gs_eq%region_info%noaxi_nodes)
        NULLIFY(graph1%kr,graph1%lc)
        graph1%nnz=graph2%nnz
        graph1%kr=>graph2%kr
        graph1%lc=>graph2%lc
        DEALLOCATE(dense_flag)
    END IF
    !---Create matrix
    graphs(1,1)%g=>graph1
    CALL create_matrix(rhs_mat,graphs,oft_lag_vec,oft_lag_vec)
    NULLIFY(graphs(1,1)%g)
    !---Add dense block for boundary
    IF(self%gs_eq%free)THEN
        ALLOCATE(bc_nodes(1))
        bc_nodes(1)%n=lag_rep%nbe
        bc_nodes(1)%v=>lag_rep%lbe
        ALLOCATE(dense_flag(lag_rep%ne))
        dense_flag=0
        dense_flag(bc_nodes(1)%v)=1
        !---Add dense blocks
        CALL graph_add_dense_blocks(graph1,graph2,dense_flag,bc_nodes)
        NULLIFY(graph1%kr,graph1%lc)
        graph1%nnz=graph2%nnz
        graph1%kr=>graph2%kr
        graph1%lc=>graph2%lc
        DEALLOCATE(dense_flag)
    END IF
    !---Add dense blocks
    IF(include_bounds)THEN
        CALL graph_add_full_col(graph1,graph2,2*lag_rep%nce,bnd_nodes)
        NULLIFY(graph1%kr,graph1%lc)
        graph1%nnz=graph2%nnz
        graph1%kr=>graph2%kr
        graph1%lc=>graph2%lc
    END IF
    !---Create matrix
    graphs(1,1)%g=>graph1
    CALL create_matrix(lhs_mat,graphs,oft_lag_vec,oft_lag_vec)
    ALLOCATE(lim_vals(a%n),ax_vals(a%n))
    lim_vals=0.d0; ax_vals=0.d0
    NULLIFY(graphs(1,1)%g)
    !---Cleanup
    CALL oft_lag_vec%delete
    DEALLOCATE(oft_lag_vec)
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
self%F%plasma_bounds=self%gs_eq%plasma_bounds
self%P%plasma_bounds=self%gs_eq%plasma_bounds
psi_norm=self%gs_eq%plasma_bounds(2)-self%gs_eq%plasma_bounds(1)
!---------------------------------------------------------------------------
! Operator integration
!---------------------------------------------------------------------------
IF(nnonaxi>0)THEN
    ALLOCATE(nonaxi_vals(self%gs_eq%region_info%block_max+1,nnonaxi))
    nonaxi_vals=0.d0
END IF
!$omp parallel private(j,rop,gop,det,lhs_vals,rhs_vals,curved,goptmp,m,vol,jc,jr,pt,psi_tmp, &
!$omp in_bounds,eta_source,gs_source,lim_loc,lim_source,ax_source,ax_loc,nonaxi_tmp)
allocate(j(lag_rep%nce)) ! Local DOF and matrix indices
allocate(rop(lag_rep%nce),gop(3,lag_rep%nce)) ! Reconstructed gradient operator
allocate(lhs_vals(lag_rep%nce,lag_rep%nce),lim_loc(lag_rep%nce))
allocate(rhs_vals(lag_rep%nce,lag_rep%nce),ax_loc(lag_rep%nce))
IF(nnonaxi>0)allocate(nonaxi_tmp(lag_rep%nce))
!$omp do schedule(static,1)
!ordered
do i=1,lag_rep%mesh%nc
    ! IF(smesh%reg(i)==1)CYCLE
    !---Get local to global DOF mapping
    call lag_rep%ncdofs(i,j)
    !---Get local reconstructed operators
    lhs_vals=0.d0; lim_loc=0.d0; ax_loc=0.d0; rhs_vals=0.d0
    IF(nnonaxi>0)nonaxi_tmp=0.d0
    do m=1,lag_rep%quad%np ! Loop over quadrature points
        call lag_rep%mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
        det=vol*lag_rep%quad%wts(m)
        pt=smesh%log2phys(i,lag_rep%quad%pts(:,m))
        psi_tmp=0.d0
        do jc=1,lag_rep%nce ! Loop over degrees of freedom
            call oft_blag_eval(lag_rep,i,jc,lag_rep%quad%pts(:,m),rop(jc))
            call oft_blag_geval(lag_rep,i,jc,lag_rep%quad%pts(:,m),gop(:,jc),goptmp)
            psi_tmp=psi_tmp+pol_vals(j(jc))*rop(jc)
        end do
        eta_source=0.d0; gs_source=0.d0; lim_source=0.d0; ax_source=0.d0
        IF(smesh%reg(i)==1)THEN
            in_bounds=gs_test_bounds(self%gs_eq,pt).AND.(psi_tmp>self%gs_eq%plasma_bounds(1))
            IF(in_bounds)THEN
                gs_source=(self%p_scale*pt(1)*pt(1)*self%P%Fpp(psi_tmp) &
                    + self%f_scale*0.5d0*self%F%fpp(psi_tmp))
                lim_source=gs_source*(1.d0-(psi_tmp-self%gs_eq%plasma_bounds(1))/psi_norm)
                ax_source=gs_source*(psi_tmp-self%gs_eq%plasma_bounds(1))/psi_norm
            END IF
            eta_source=sigma_plasma
        ELSE IF(smesh%reg(i)>1.AND.(self%eta_reg(smesh%reg(i))>0.d0))THEN
            eta_source=1.d0/self%eta_reg(smesh%reg(i))
        END IF
        !---Compute local matrix contributions
        do jr=1,lag_rep%nce
            lim_loc(jr) = lim_loc(jr) + rop(jr)*lim_source*det/(pt(1)+gs_epsilon)
            ax_loc(jr) = ax_loc(jr) + rop(jr)*ax_source*det/(pt(1)+gs_epsilon)
            do jc=1,lag_rep%nce
                lhs_vals(jr,jc) = lhs_vals(jr,jc) + (DOT_PRODUCT(gop(1:2,jr),gop(1:2,jc)) &
                    - rop(jr)*rop(jc)*gs_source)*det/(pt(1)+gs_epsilon)
                rhs_vals(jr,jc) = rhs_vals(jr,jc) + rop(jr)*rop(jc)*eta_source*det/(pt(1)+gs_epsilon)
            end do
            IF(nnonaxi>0.AND.self%gs_eq%region_info%reg_map(smesh%reg(i))>0)THEN
                nonaxi_tmp(jr)=nonaxi_tmp(jr) + rop(jr)*eta_source*det/(pt(1)+gs_epsilon)
            END IF
        end do
        IF(nnonaxi>0.AND.self%gs_eq%region_info%reg_map(smesh%reg(i))>0)THEN
            !$omp atomic
            nonaxi_vals(self%gs_eq%region_info%block_max+1,self%gs_eq%region_info%reg_map(smesh%reg(i))) = &
              nonaxi_vals(self%gs_eq%region_info%block_max+1,self%gs_eq%region_info%reg_map(smesh%reg(i))) + det
        END IF
    end do
    !---Apply bc to local matrix
    DO jr=1,lag_rep%nce
        IF(lag_rep%be(j(jr)))THEN
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
    call rhs_mat%atomic_add_values(j,j,rhs_vals,lag_rep%nce,lag_rep%nce)
    call lhs_mat%atomic_add_values(j,j,lhs_vals,lag_rep%nce,lag_rep%nce)
    DO jr=1,lag_rep%nce
        !$omp atomic
        lim_vals(j(jr))=lim_vals(j(jr))+lim_loc(jr)
        !$omp atomic
        ax_vals(j(jr))=ax_vals(j(jr))+ax_loc(jr)
    END DO
    IF(nnonaxi>0.AND.self%gs_eq%region_info%reg_map(smesh%reg(i))>0)THEN
        DO jr=1,lag_rep%nce
          !$omp atomic
          nonaxi_vals(self%gs_eq%region_info%node_mark(smesh%reg(i),j(jr)),self%gs_eq%region_info%reg_map(smesh%reg(i))) = &
            nonaxi_vals(self%gs_eq%region_info%node_mark(smesh%reg(i),j(jr)),self%gs_eq%region_info%reg_map(smesh%reg(i))) &
            + nonaxi_tmp(jr)
        END DO
    END IF
    !!$omp end ordered
end do
IF(nnonaxi>0)THEN
    !$omp do schedule(dynamic,1)
    !ordered
    do i=1,lag_rep%mesh%nc
        IF(self%gs_eq%region_info%reg_map(smesh%reg(i))==0)CYCLE
        !---Get local to global DOF mapping
        call lag_rep%ncdofs(i,j)
        !---Get local reconstructed operators
        lim_loc=0.d0
        do m=1,lag_rep%quad%np ! Loop over quadrature points
            call lag_rep%mesh%jacobian(i,lag_rep%quad%pts(:,m),goptmp,vol)
            det=vol*lag_rep%quad%wts(m)
            do jc=1,lag_rep%nce ! Loop over degrees of freedom
                call oft_blag_eval(lag_rep,i,jc,lag_rep%quad%pts(:,m),rop(jc))
                lim_loc(jc) = lim_loc(jc) + rop(jc)*det
            end do
        end do
        lim_loc=-lim_loc/nonaxi_vals(self%gs_eq%region_info%block_max+1,self%gs_eq%region_info%reg_map(smesh%reg(i)))
        !---Add local values to global matrix
        m=self%gs_eq%region_info%reg_map(smesh%reg(i))
        !!$omp ordered
        do jc=1,lag_rep%nce
            call rhs_mat%atomic_add_values(j(jc:jc),self%gs_eq%region_info%noaxi_nodes(m)%v, &
              lim_loc(jc)*nonaxi_vals(:,m),1,self%gs_eq%region_info%noaxi_nodes(m)%n)
            call lhs_mat%atomic_add_values(j(jc:jc),self%gs_eq%region_info%noaxi_nodes(m)%v, &
              -sigma*lim_loc(jc)*nonaxi_vals(:,m),1,self%gs_eq%region_info%noaxi_nodes(m)%n)
        end do
        !!$omp end ordered
    end do
END IF
deallocate(j,rop,gop,lhs_vals,lim_loc,ax_loc,rhs_vals)
IF(nnonaxi>0)deallocate(nonaxi_tmp)
!$omp end parallel
DEALLOCATE(pol_vals)
IF(nnonaxi>0)DEALLOCATE(nonaxi_vals)
IF(include_bounds)THEN
  ALLOCATE(j(1),lhs_vals(1,1))
  DO jr=1,lag_rep%nce
    j=bnd_nodes(jr)
    DO i=1,lag_rep%ne
        lhs_vals=lim_vals(i)*lim_weights(jr)
        rowtmp=i
        CALL lhs_mat%add_values(rowtmp,j,lhs_vals,1,1)
    END DO
    j=bnd_nodes(lag_rep%nce+jr)
    DO i=1,lag_rep%ne
        lhs_vals=ax_vals(i)*ax_weights(jr)
        rowtmp=i
        CALL lhs_mat%add_values(rowtmp,j,lhs_vals,1,1)
    END DO
  END DO
  DEALLOCATE(j,lhs_vals)
END IF
DEALLOCATE(bnd_nodes,lim_weights,ax_weights)
DEALLOCATE(lim_vals,ax_vals)
CALL set_bcmat(self%gs_eq,lhs_mat)
!---Assemble matrix
CALL lag_rep%vec_create(oft_lag_vec)
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
end subroutine build_linearized
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine apply_wrap(self,a,b)
class(eig_wrapper), intent(inout) :: self !< NL operator object
class(oft_vector), target, intent(inout) :: a !< Source field
class(oft_vector), intent(inout) :: b !< Result of metric function
class(oft_vector), pointer :: tmp_vec
CALL a%new(tmp_vec)
CALL self%rhs_mat%apply(a,tmp_vec)
CALL self%lhs_inv%apply(b,tmp_vec)
CALL tmp_vec%delete()
DEALLOCATE(tmp_vec)
DEBUG_STACK_POP
end subroutine apply_wrap
end module oft_gs_td