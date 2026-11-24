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

By default the native solver performs local solves (ie as a sub-solver to Block Jacobi). With the
`superd`, `mumps`, and `strumpack` packages, the full distributed system can instead be factored across
all MPI tasks by setting `distributed` to true and using the solver directly on the global system (not
as a sub-solver to Block Jacobi). This option has no effect in serial runs or for local systems.
Distributed factorization with `strumpack` requires STRUMPACK built with MPI support (`--strumpack_mpi=1`
with `build_libs.py`), in which case the STRUMPACK options below also apply.

**XML Spec:** `<solver type="lu">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `package`  | Solver package to use: native {super, superd, umfpack, mumps, mkl, lapack, strumpack}, PETSc {petsc, super, superd, umfpack, mumps} | str | No |

Native only options

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `distributed` | Factor full distributed system across all MPI tasks (`superd`, `mumps`, and `strumpack` only) | bool | No |

Native STRUMPACK options, in a `<strumpack>` group (values of `-1` or omitted use the STRUMPACK defaults)

| Group          | Description | Type      | List |
|----------------|-------------|-----------|------|
| `comp_type`    | Compression of frontal matrices {none, hss, blr} (default: none) | str | No |
| `comp_rel_tol` | Relative compression tolerance | float | No |
| `comp_abs_tol` | Absolute compression tolerance | float | No |
| `leaf_size`    | Size of leaf/tile blocks for compression | int | No |
| `max_rank`     | Maximum rank of compressed blocks | int | No |
| `min_sep_size` | Minimum separator size for compression | int | No |
| `verbose`      | Print STRUMPACK information | bool | No |

\note With compression the STRUMPACK factorization is approximate, which is intended for use as a
preconditioner (eg. Block Jacobi sub-solver). STRUMPACK's iterative refinement is not supported. A
refactorization with an unchanged non-zero pattern reuses the reordering and MC64 matching/scaling from the
last full factorization (local and distributed solves). STRUMPACK only pivots within each front, so after large
changes in matrix values the solution can be inaccurate unless a full factorization is requested
(`update(.TRUE.)`, as for a change in the non-zero pattern), even though the pattern has not changed.

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

This method is provided as a preconditioner through both the native and PETSc backends through
\ref oft_native_solvers::oft_schwarz_precond "oft_schwarz_precond" and
\ref oft_petsc_solvers::oft_petsc_asprecond "oft_petsc_asprecond". Subdomains are grown from a
non-overlapping partition by `overlap` layers along the matrix graph and solved with the
sub-solver given by the nested `pre` element. In the native backend, with `nlocal=1` each process
has one subdomain made of its owned rows, grown across processes; with several local blocks
(`nlocal>1`, or negative for field based blocks) the blocks are grown within the owned rows of
each process. `overlap=0` is block-Jacobi.

**XML Spec:** `<pre type="add_schwarz">`

| Group      | Description | Type      | List |
|------------|-------------|-----------|------|
| `nlocal`   | Number of local blocks (-1 for field based) | int | Yes |
| `groups`   | Group index for each field in multi-field vector [n] (native only) | int | Yes |
| `overlap`  | Overlap size (native default: 1; layers of cells in cell mode) | int | Yes |
| `type`     | Combination of subdomain solutions: `restricted` (default) or `additive` (native only) | str | No |
| `overlap_mode` | How blocks are grown: `graph` (default, matrix graph) or `cell` (layers of mesh cells; native only) | str | No |
| `adjacency` | Cell contact defining a layer in cell mode: `point` (default), `edge` or `face` (native only) | str | No |
| `base_layers` | Layers of base mesh cells bounding the cells exchanged between processes in cell mode (default -1: equal to `overlap`, always sufficient; smaller values can truncate blocks at process boundaries and print a warning; native only) | int | Yes |

**Cell mode** (`<overlap_mode>cell</overlap_mode>`, native only): blocks are grown by `overlap`
layers of mesh cells instead of layers of the matrix graph. Cells are adjacent if they share at
least a point, an edge or a face (`adjacency`), and a block contains every DoF of its cells. The
cells and the number of blocks per process are defined by a cell overlap builder
(\ref fem_cell_overlap::oft_fem_cell_overlap "oft_fem_cell_overlap") that the physics module or
application attaches to the matrix (`cell_ovl`); cell mode aborts if none is attached, and the
builder's block count replaces `nlocal` (field based blocks and `boverlap` are not supported). Blocks
may cross process boundaries for any number of blocks per process: the union of the blocks of each
process is assembled as one subdomain across processes, so the setup needs a fixed number of
exchanges whatever the overlap. With `overlap=0` the blocks are the DoFs of the core cells, which
already overlap at process and block interfaces. The restricted combination takes each owned DoF
from the one block it is assigned to.

With `debug>=1` in `runtime_options`, the setup of the Schwarz subdomains (both modes) prints the
time of each phase (maximum over processes) and the number of int8 index values sent between
processes (total and maximum per process) as `Setup timing:` lines.

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
