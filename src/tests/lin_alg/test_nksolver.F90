!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
! Regression test for oft_nksolver with opt-in Eisenstat-Walker adaptive forcing.
! Solves a nonlinear reaction-diffusion problem L*u + lambda*(w*u^3) = b on the
! unit cube with homogeneous Dirichlet boundary conditions using both fixed tight
! linear tolerances (ew_forcing=.FALSE.) and Eisenstat-Walker Choice 2 forcing
! (ew_forcing=.TRUE.), checking solution equivalence, linear iteration reduction,
! and restoration of J_inv tolerances on converged and non-converged exits.
!---------------------------------------------------------------------------------
MODULE test_nksolver_ops
USE oft_base
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_deriv_matrices, ONLY: oft_noop_matrix
USE oft_solver_base, ONLY: oft_solver
USE fem_base, ONLY: oft_ml_fem_type
USE oft_lag_operators, ONLY: oft_lag_getlop
IMPLICIT NONE
PRIVATE

TYPE, PUBLIC, EXTENDS(oft_noop_matrix) :: nl_rd_op
  REAL(r8) :: lambda = 50.d0
  CLASS(oft_matrix), POINTER :: lop => NULL()
  CLASS(oft_vector), POINTER :: weight => NULL()
  CLASS(oft_vector), POINTER :: tmp => NULL()
CONTAINS
  PROCEDURE :: apply_real => nl_rd_apply
  PROCEDURE :: delete => nl_rd_delete
END TYPE nl_rd_op

REAL(r8), PUBLIC :: nl_lambda = 50.d0
TYPE(oft_ml_fem_type), POINTER, PUBLIC :: nl_ml_lag => NULL()
CLASS(oft_matrix), POINTER, PUBLIC :: nl_jmat => NULL()
CLASS(oft_solver), POINTER, PUBLIC :: nl_jinv => NULL()
CLASS(oft_vector), POINTER, PUBLIC :: nl_weight => NULL()

PUBLIC nl_update_jacobian
CONTAINS

! Evaluate F(u) = L*u + lambda*(w * u^3)
SUBROUTINE nl_rd_apply(self,a,b)
CLASS(nl_rd_op), INTENT(inout) :: self
CLASS(oft_vector), TARGET, INTENT(inout) :: a
CLASS(oft_vector), INTENT(inout) :: b
CALL self%lop%apply(a,b)
IF(.NOT.ASSOCIATED(self%tmp))CALL a%new(self%tmp)
CALL self%tmp%copy(a)
CALL self%tmp%mult(a)
CALL self%tmp%mult(a)
CALL self%tmp%mult(self%weight)
CALL b%add(1.d0,self%lambda,self%tmp)
END SUBROUTINE nl_rd_apply

! Release scratch vector
SUBROUTINE nl_rd_delete(self)
CLASS(nl_rd_op), INTENT(inout) :: self
IF(ASSOCIATED(self%tmp))THEN
  CALL self%tmp%delete
  DEALLOCATE(self%tmp)
END IF
NULLIFY(self%lop,self%weight)
END SUBROUTINE nl_rd_delete

! Assemble exact Jacobian J(u) = L + diag(3*lambda*w*u^2)
SUBROUTINE nl_update_jacobian(a)
CLASS(oft_vector), TARGET, INTENT(inout) :: a
CLASS(oft_vector), POINTER :: diag_tmp
REAL(r8), POINTER, DIMENSION(:) :: u_loc,w_loc
REAL(r8) :: dval(1,1)
INTEGER(i4) :: k,jr,ind(1)
CALL oft_lag_getlop(nl_ml_lag%current_level,nl_jmat,'zerob')
NULLIFY(u_loc,w_loc)
CALL a%get_local(u_loc)
CALL nl_weight%get_local(w_loc)
! Add diagonal cubic derivative on owned DOFs only (w is zero on Dirichlet boundary)
DO k=1,a%map(1)%nslice
  jr=a%map(1)%slice(k)
  IF(w_loc(jr)/=0.d0)THEN
    dval(1,1)=3.d0*nl_lambda*w_loc(jr)*(u_loc(jr)**2)
    ind(1)=jr
    CALL nl_jmat%add_values(ind,ind,dval,1,1)
  END IF
END DO
DEALLOCATE(u_loc,w_loc)
CALL a%new(diag_tmp)
CALL nl_jmat%assemble(diag_tmp)
CALL diag_tmp%delete
DEALLOCATE(diag_tmp)
CALL nl_jinv%update(.FALSE.)
END SUBROUTINE nl_update_jacobian

END MODULE test_nksolver_ops

PROGRAM test_nksolver
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_nksolver
USE oft_solver_utils, ONLY: create_gmres_solver, create_diag_pre
USE fem_base, ONLY: oft_ml_fem_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: oft_lag_zerob, oft_lag_getlop, oft_lag_getmop
USE test_nksolver_ops, ONLY: nl_rd_op, nl_lambda, nl_ml_lag, nl_jmat, nl_jinv, &
  nl_weight, nl_update_jacobian
IMPLICIT NONE
INTEGER(i4), PARAMETER :: order=3
INTEGER(i4) :: io_unit,nfail,nlits_fixed,lits_fixed,nlits_ew,lits_ew
REAL(r8), PARAMETER :: rtol_init=1.d-10, atol_init=1.d-14
REAL(r8) :: u_norm,diff_norm
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange
TYPE(oft_lag_zerob), TARGET :: lag_zerob
TYPE(nl_rd_op), TARGET :: nl_op
TYPE(oft_nksolver) :: nksolver
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u_fixed => NULL(), u_ew => NULL(), rhs => NULL(), diff => NULL()

CALL oft_init
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,minlev=-1)
lag_zerob%ML_lag_rep=>ML_oft_lagrange
oft_env%pm=.FALSE.

CALL ML_oft_lagrange%set_level(ML_oft_lagrange%nlevels)
CALL ML_oft_lagrange%vec_create(u_fixed)
CALL ML_oft_lagrange%vec_create(u_ew)
CALL ML_oft_lagrange%vec_create(rhs)
CALL ML_oft_lagrange%vec_create(diff)
CALL ML_oft_lagrange%vec_create(nl_weight)

! Build base FE operators
CALL oft_lag_getlop(ML_oft_lagrange%current_level,lop,'zerob')
CALL oft_lag_getmop(ML_oft_lagrange%current_level,mop,'none')
CALL nl_weight%copy(mop%D)
CALL lag_zerob%apply(nl_weight)

! Build RHS b = 100 * M * 1 (strong forcing so cubic nonlinearity requires ~6 Newton steps)
CALL u_fixed%set(100.d0)
CALL mop%apply(u_fixed,rhs)
CALL lag_zerob%apply(rhs)

! Setup nonlinear operator and Jacobian workspace
nl_lambda=50.d0
nl_op%lambda=nl_lambda
nl_op%lop=>lop
nl_op%weight=>nl_weight
nl_ml_lag=>ML_oft_lagrange
CALL oft_lag_getlop(ML_oft_lagrange%current_level,nl_jmat,'zerob')

! Setup GMRES linear solver (native or PETSc depending on use_petsc)
CALL create_gmres_solver(nl_jinv,nrits=20)
nl_jinv%A=>nl_jmat
nl_jinv%its=500
nl_jinv%rtol=rtol_init
nl_jinv%atol=atol_init
CALL create_diag_pre(nl_jinv%pre)
nl_jinv%pre%A=>nl_jmat

! Setup Newton-Krylov solver
nksolver%A=>nl_op
nksolver%J_inv=>nl_jinv
nksolver%J_update=>nl_update_jacobian
nksolver%bc=>lag_zerob
nksolver%its=30
nksolver%atol=1.d-10
nksolver%rtol=1.d-10
nksolver%up_freq=1

nfail=0

! Baseline solve with fixed tight linear tolerance (ew_forcing = .FALSE.)
nksolver%ew_forcing=.FALSE.
CALL u_fixed%set(0.d0)
CALL nksolver%apply(u_fixed,rhs)
nlits_fixed=nksolver%nlits
lits_fixed=nksolver%lits
IF(nksolver%cits<0)THEN
  IF(oft_env%head_proc)WRITE(*,*)'FAIL: Fixed-tolerance Newton solve did not converge'
  nfail=nfail+1
END IF

! Adaptive solve with Eisenstat-Walker Choice 2 forcing (ew_forcing = .TRUE.)
nksolver%ew_forcing=.TRUE.
CALL u_ew%set(0.d0)
CALL nksolver%apply(u_ew,rhs)
nlits_ew=nksolver%nlits
lits_ew=nksolver%lits
IF(nksolver%cits<0)THEN
  IF(oft_env%head_proc)WRITE(*,*)'FAIL: Eisenstat-Walker Newton solve did not converge'
  nfail=nfail+1
END IF

! Verify J_inv tolerances were restored after converged EW solve
IF(nl_jinv%rtol/=rtol_init.OR.nl_jinv%atol/=atol_init)THEN
  IF(oft_env%head_proc)WRITE(*,*)'FAIL: J_inv tolerances not restored after converged EW solve', &
    nl_jinv%rtol,nl_jinv%atol
  nfail=nfail+1
END IF

! Verify both solves agree on the nonlinear solution
CALL diff%add(0.d0,1.d0,u_ew,-1.d0,u_fixed)
u_norm=SQRT(u_fixed%dot(u_fixed))
diff_norm=SQRT(diff%dot(diff))/u_norm
IF(diff_norm>1.d-8)THEN
  IF(oft_env%head_proc)WRITE(*,*)'FAIL: Solution mismatch between fixed and EW solves:',diff_norm
  nfail=nfail+1
END IF

! Verify EW uses strictly fewer linear iterations than fixed tight tolerance
IF(lits_ew>=lits_fixed)THEN
  IF(oft_env%head_proc)WRITE(*,*)'FAIL: EW did not reduce linear iterations:',lits_ew,lits_fixed
  nfail=nfail+1
END IF

! Verify J_inv tolerances are also restored on non-converged (max-iterations) exit
nksolver%its=1
CALL diff%set(0.d0)
CALL nksolver%apply(diff,rhs)
IF(nksolver%cits>=0)THEN
  IF(oft_env%head_proc)WRITE(*,*)'FAIL: Expected max-iteration exit with its=1'
  nfail=nfail+1
END IF
IF(nl_jinv%rtol/=rtol_init.OR.nl_jinv%atol/=atol_init)THEN
  IF(oft_env%head_proc)WRITE(*,*)'FAIL: J_inv tolerances not restored after max-its exit', &
    nl_jinv%rtol,nl_jinv%atol
  nfail=nfail+1
END IF

IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='nksolver.results')
  WRITE(io_unit,*)nfail
  WRITE(io_unit,*)nlits_fixed,lits_fixed
  WRITE(io_unit,*)nlits_ew,lits_ew
  WRITE(io_unit,*)diff_norm
  CLOSE(io_unit)
END IF

! Cleanup
CALL nksolver%delete
CALL nl_op%delete
CALL nl_jinv%delete(propogate=.TRUE.)
DEALLOCATE(nl_jinv)
CALL nl_jmat%delete
CALL lop%delete
CALL mop%delete
DEALLOCATE(nl_jmat,lop,mop)
CALL u_fixed%delete
CALL u_ew%delete
CALL rhs%delete
CALL diff%delete
CALL nl_weight%delete
DEALLOCATE(u_fixed,u_ew,rhs,diff,nl_weight)
CALL oft_finalize
END PROGRAM test_nksolver
