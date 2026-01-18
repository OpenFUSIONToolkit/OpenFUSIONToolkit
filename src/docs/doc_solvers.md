Available Solvers     {#doc_solvers}
=================

[TOC]

This page outlines the solvers that are currently available through Open FUSION Toolkit solver objects. Additional
solvers may be used through backend implementations, but only internally controlled and supported solvers
are outlined here.

\section doc_solvers_petsc Mixing Native and PETSc Solvers

Solver interfaces have been designed to allow mixing of native and PETSc methods. Most native solvers are
written with generic API calls for linear algebra and as a result can be used when either backend is used
for vectors and matrices. There are many cases however where the PETSc implementation of a particular
solver/preconditioner is more optimal for its matrices and vectors. Additionally, there are some
methods, ex. direct methods, which are only available through PETSc. The
\ref oft_petsc_solvers::oft_petsc_pre_solver "oft_petsc_pre_solver" class is provided to allow using PETSc
solvers as preconditioners for native solvers. This class acts as a wrapper around the PETSc solver, allowing
it to be called using the standard solver calling sequence.  

\section doc_solvers_basic Basic Solvers and Preconditioners

\subsection doc_solvers_basic_jacobi Jacobi Iteration

Jacobi iterations solves a system of linear equations by iterating on the method

\f[ x^{n+1} = x^n + D^{-1}(A*x^n - y) \f]

where \f$ D \f$ is a matrix of only the diagonal elements of \f$ A \f$. Functionality is provided through
either a general native interface \ref oft_solvers::oft_diag_scale "oft_diag_scale" and an optimized PETSc
interface \ref oft_petsc_solvers::oft_petsc_diagprecond "oft_petsc_diagprecond".

**XML Spec:** `<pre type="jacobi">`

\subsection doc_solvers_basic_lu LU Factorization

This method is supported through both the native and PETSc backends through
\ref oft_lu::oft_lusolver "oft_lusolver" and
\ref oft_petsc_solvers::oft_petsc_direct_solver "oft_petsc_direct_solver" for a full solver or
\ref oft_petsc_solvers::oft_petsc_luprecond "oft_petsc_luprecond" as a PETSc preconditioner.

\warning The native solver currently only supports local solves (ie as a sub-solver to Block Jacobi)

**XML Spec:** `<solver type="lu">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `package`  | Solver package to use {super, superd, mumps} | str | No |

PETSc only options

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `type`     | Inverse type {lu, ilu} | str | No |

\subsection doc_solvers_basic_bjacobi Block Jacobi

This method is supported through both the native and PETSc backends through
\ref oft_solvers::oft_bjprecond "oft_bjprecond" and
\ref oft_petsc_solvers::oft_petsc_asprecond "oft_petsc_asprecond".

**XML Spec:** `<pre type="block_jacobi">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `nlocal`   | Number of local blocks (-1 for field based) | int | Yes |
| `groups`   | Group index for each field in multi-field vector [n] | int | Yes |
| `boverlap` | Overlap boundary rows (native only) | bool | false |

\subsection doc_solvers_basic_schwarz Additive Schwarz

This method is currently only supported with the PETSc backend through
\ref oft_petsc_solvers::oft_petsc_asprecond "oft_petsc_asprecond" as a preconditioner.

**XML Spec:** `<pre type="add_schwarz">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `nlocal`   | Number of local blocks (-1 for field based) | int | Yes |
| `overlap`  | Overlap size | int | Yes |

\section doc_solvers_krylov Krylov Solvers

\subsection doc_solvers_krylov_cg Conjugate Gradient

This method is provided through both the native and PETSc backends through \ref oft_cg::oft_cg_solver
"oft_cg_solver" and \ref oft_petsc_solvers::oft_petsc_cg_solver "oft_petsc_cg_solver".

**XML Spec:** `<solver type="cg">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `atol`     | Absolute convergence tolerance | real | Yes |
| `rtol`     | Relative convergence tolerance | real | Yes |
| `its`      | Iteration limit | int | Yes |

\subsection doc_solvers_krylov_gmres Generalized Minimum Residual Method

This method is provided through both the native and PETSc backends through \ref oft_gmres::oft_gmres_solver
"oft_gmres_solver" and \ref oft_petsc_solvers::oft_petsc_gmres_solver "oft_petsc_gmres_solver".

**XML Spec:** `<solver type="gmres">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `atol`     | Absolute convergence tolerance | real | Yes |
| `rtol`     | Relative convergence tolerance | real | Yes |
| `its`      | Iteration limit | int | Yes |
| `nrit`     | Restart size | int | Yes |
