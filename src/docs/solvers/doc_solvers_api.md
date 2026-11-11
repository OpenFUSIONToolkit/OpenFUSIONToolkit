Using Solvers in Fortran     {#doc_solvers_api}
========================

[TOC]

This page describes how to create, configure and apply solver objects in Fortran.

\section doc_solvers_usage Using Solver Objects

All linear solvers and preconditioners extend \ref oft_solver_base::oft_solver "oft_solver" and share the same
calling sequence. Preconditioners are themselves solver objects, attached to a parent solver through its `pre` field.

~~~~~~~~~~~~~~~~~~~~~{.f90}
CLASS(oft_solver), POINTER :: solver
CALL create_cg_solver(solver)    ! Create a CG solver for the active backend
solver%A=>mat                    ! Matrix to invert
solver%its=-3                    ! Set options (see Solver XML Reference)
solver%atol=1.d-10
CALL create_diag_pre(solver%pre) ! Optional preconditioner (uses solver%A)
CALL solver%apply(u,g)           ! Solve A*u = g
CALL solver%update(.TRUE.)       ! After changing A (.TRUE. if non-zero pattern changed)
CALL solver%delete(.TRUE.)       ! Free storage (.TRUE. also deletes preconditioner)
~~~~~~~~~~~~~~~~~~~~~

In `apply(u,g)` the vector `u` holds the initial guess on input and the solution on output. The vector `g` holds the
right-hand side on input and is overwritten; iterative solvers return the final residual in `g`. After a solve
`cits` holds the number of iterations taken. A negative value indicates failure, eg. breakdown or the native
GMRES solver reaching its iteration limit.

Solver objects are created with the helper routines in \ref oft_solver_utils, which select the native or PETSc
implementation according to the active backend (\ref doc_solvers_backends):

| Routine | Creates |
|---------|---------|
| \ref oft_solver_utils::create_cg_solver "create_cg_solver" | Conjugate Gradient solver |
| \ref oft_solver_utils::create_gmres_solver "create_gmres_solver" | GMRES solver |
| \ref oft_solver_utils::create_diag_pre "create_diag_pre" | Jacobi (diagonal scaling) preconditioner |
| \ref oft_solver_utils::create_ilu_pre "create_ilu_pre" | ILU(0) preconditioner |
| \ref oft_solver_utils::create_bjacobi_pre "create_bjacobi_pre" | Block Jacobi preconditioner with LU sub-solves |
| \ref oft_solver_utils::create_mlpre "create_mlpre" | Multigrid preconditioner (\ref doc_solvers_mg) |
| \ref oft_solver_utils::create_solver_xml "create_solver_xml" | Any solver from an XML definition (\ref doc_solvers_xml_ref) |

Options are set as fields of the solver object before the first call to `apply`. The options shared by all
solvers are listed in \ref doc_solvers_common, and method specific options in \ref doc_solvers_xml_ref.

\section doc_solvers_eig Eigenvalue Solvers

Eigenvalue solvers extend \ref oft_solver_base::oft_eigsolver "oft_eigsolver" and are configured in Fortran only:
 - \ref oft_native_solvers::oft_native_cg_eigsolver "oft_native_cg_eigsolver": Lowest eigenvalue of a symmetric
   positive-definite generalized system using non-linear CG
 - \ref oft_arpack::oft_irlm_eigsolver "oft_irlm_eigsolver" (symmetric) and
   \ref oft_arpack::oft_iram_eigsolver "oft_iram_eigsolver" (general): Multiple eigenvalues using
   [ARPACK](https://github.com/opencollab/arpack-ng) (requires ARPACK)

\section doc_solvers_nl Non-Linear Solver

\ref oft_native_solvers::oft_nksolver "oft_nksolver" provides a Newton iteration with backtracking for
\f$ M(u) = g \f$, using a linear solver (`J_inv`) for the Jacobian system. It is configured in Fortran only.
