Solver XML Reference     {#doc_solvers_xml_ref}
====================

[TOC]

This page describes how to define solvers in the XML input file and lists the options for each method. The same
options are available as fields of the corresponding Fortran classes (\ref doc_solvers_api).

\section doc_solvers_xml XML Specification

Solvers can be defined at run time in the XML input file. A definition is a tree of alternating `solver` and `pre`
elements, where the `type` attribute selects the method and options are given as child elements.

~~~~~~~~~~~~~~~~~~~~~{.xml}
<solver type="gmres">        <!-- Solver -->
  <its>20</its>
  <nrits>20</nrits>
  <pre type="block_jacobi">  <!-- Preconditioner for the parent solver -->
    <nlocal>2</nlocal>
    <solver type="lu"/>      <!-- Sub-solver used by the parent preconditioner -->
  </pre>
</solver>
~~~~~~~~~~~~~~~~~~~~~

 - A `solver` element may contain one `pre` element, and a `pre` element may contain one `solver` element (used
   by block preconditioners as the local solver).
 - The top-level element is located by the calling module, which may use a different name (eg. `<pre type="gmres">`
   for \ref xmhd). Only its `type` attribute and contents are used. See \ref doc_solvers_xml_modules.
 - Numeric options accept a single value or a comma/space separated list. A list is indexed by multigrid level
   (\ref doc_solvers_mg); outside of a multigrid preconditioner the first value is used. Options marked
   "Per-level: No" accept only a single value.
 - Logical options accept `true`/`false` (also `T`/`F` or `1`/`0`).
 - The implementation is chosen by the active backend (\ref doc_solvers_backends). For the `cg`, `gmres` and
   `sjacobi` solver types the attribute `native="true"` forces the native implementation when PETSc is in use.

\subsection doc_solvers_xml_types Available Types

| `<solver type=...>` | Method | Native class | PETSc class |
|---------------------|--------|--------------|-------------|
| `cg`      | \ref doc_solvers_cg | \ref oft_native_solvers::oft_native_cg_solver "oft_native_cg_solver" | \ref oft_petsc_solvers::oft_petsc_cg_solver "oft_petsc_cg_solver" |
| `gmres`   | \ref doc_solvers_gmres | \ref oft_native_solvers::oft_native_gmres_solver "oft_native_gmres_solver" | \ref oft_petsc_solvers::oft_petsc_gmres_solver "oft_petsc_gmres_solver" |
| `sjacobi` | \ref doc_solvers_sjacobi | \ref oft_native_solvers::oft_jblock_precond "oft_jblock_precond" | \ref oft_petsc_solvers::oft_petsc_sjacobi_solver "oft_petsc_sjacobi_solver" |
| `lu`      | \ref doc_solvers_lu | \ref oft_lu::oft_lusolver "oft_lusolver" | \ref oft_petsc_solvers::oft_petsc_direct_solver "oft_petsc_direct_solver" |
| `ilu`     | \ref doc_solvers_ilu | \ref oft_lu::oft_ilusolver "oft_ilusolver" | Not available |

| `<pre type=...>` | Method | Native class | PETSc class |
|------------------|--------|--------------|-------------|
| `jacobi`       | \ref doc_solvers_jacobi | \ref oft_native_solvers::oft_diag_scale "oft_diag_scale" | \ref oft_petsc_solvers::oft_petsc_diagprecond "oft_petsc_diagprecond" |
| `block_jacobi` | \ref doc_solvers_bjacobi | \ref oft_native_solvers::oft_bjprecond "oft_bjprecond" | \ref oft_petsc_solvers::oft_petsc_asprecond "oft_petsc_asprecond" |
| `add_schwarz`  | \ref doc_solvers_bjacobi | Not available | \ref oft_petsc_solvers::oft_petsc_asprecond "oft_petsc_asprecond" |
| `lu`           | \ref doc_solvers_lu | \ref oft_lu::oft_lusolver "oft_lusolver" | \ref oft_petsc_solvers::oft_petsc_luprecond "oft_petsc_luprecond" |
| `ilu`          | \ref doc_solvers_ilu | \ref oft_lu::oft_ilusolver "oft_ilusolver" | Not available |

Not all combinations of native and PETSc methods are valid, see \ref doc_solvers_backends_mixing.

\subsection doc_solvers_xml_modules Modules Using XML Definitions

| Module | XML element | Definition |
|--------|-------------|------------|
| \ref xmhd, \ref xmhd_lag | `<xmhd><pre>` | Preconditioner for the Newton-Krylov GMRES solver (\ref doc_solvers_mg "multigrid definition" when more than one level is used) |
| \ref xmhd_2d | `<xmhd2d><pre>` | Preconditioner for the Newton-Krylov GMRES solver |
| \ref oft_lag_operators | `<lagrange><lop>` | Multigrid preconditioner for the Laplacian |
| \ref oft_h1_operators | `<h1><lop>` | Multigrid preconditioner for the Laplacian |
| \ref oft_hcurl_operators | `<hcurl><wop>`, `<hcurl><jmlb>` | Multigrid preconditioners for the curl-curl and Jacobian operators |
| \ref oft_hcurl_grad_operators | `<hcurl_grad><mop>` | Multigrid preconditioner for the mass matrix |

\section doc_solvers_common Common Options

The following options are shared by all iterative solvers. They can be set from XML or as fields of
\ref oft_solver_base::oft_solver "oft_solver" (`pm` is Fortran only).

| Option | Description | Type | Default |
|--------|-------------|------|---------|
| `its`  | Iteration limit (see below) | int | -1 |
| `atol` | Absolute tolerance, stop when \f$ \|r\| < atol \f$ | real | 1E-14 |
| `rtol` | Relative tolerance, stop when \f$ \|r\|/\|r_0\| < rtol \f$ | real | 1E-14 |
| `pm`   | Print convergence information when called by a parent solver; top-level output follows `oft_env%pm` (Fortran only) | bool | F |

For native solvers `its` also controls the stopping criteria:
 - `its > 0`: Stop after at most `its` iterations
 - `its = 0`: Return immediately without iterating
 - `its = -1` or `-2`: Stop when the solution stagnates in single (`-1`) or double (`-2`) precision. For CG this
   is when the energy functional stops decreasing; for GMRES when \f$ \|r\|/\|u\| \f$ falls below 1E-7 or 1E-14.
 - `its < -2`: No iteration limit, stop on `atol` and `rtol` only

PETSc solvers interpret `its` differently, see \ref doc_solvers_backends_petsc.

\section doc_solvers_krylov Krylov Solvers

\subsection doc_solvers_cg Conjugate Gradient

[Conjugate Gradient](https://en.wikipedia.org/wiki/Conjugate_gradient_method) method for symmetric positive-definite
matrices. Supports the \ref doc_solvers_common "common options" only.

**XML:** `<solver type="cg">`

\subsection doc_solvers_gmres Generalized Minimum Residual

Restarted [GMRES](https://en.wikipedia.org/wiki/Generalized_minimal_residual_method) method for general matrices.
Both backends use the flexible variant (FGMRES), which allows preconditioners that change between iterations (eg.
an iterative solver or multigrid). Supports the \ref doc_solvers_common "common options" and:

**XML:** `<solver type="gmres">`

| Option  | Description | Type | Default | Per-level |
|---------|-------------|------|---------|-----------|
| `nrits` | Restart length (iterations per cycle), 0 skips the solve | int | 10 | Yes |

\section doc_solvers_basic Smoothers and Preconditioners

\subsection doc_solvers_jacobi Jacobi

Diagonal scaling, \f$ u = D^{-1} g \f$ where \f$ D \f$ is the diagonal of \f$ A \f$. No options.

**XML:** `<pre type="jacobi">`

\subsection doc_solvers_sjacobi Symmetric Jacobi

Damped Jacobi (Richardson) iteration, \f$ u^{n+1} = u^n + \omega D^{-1}(g - A u^n) \f$, designed as a multigrid
smoother. Successive calls alternate between a down sweep (initial guess ignored) and the matching up sweep, so a
single object must be used for both sweeps (`direction="both"` in \ref doc_solvers_mg).

**XML:** `<solver type="sjacobi">`

| Option | Description | Type | Default | Per-level |
|--------|-------------|------|---------|-----------|
| `its`  | Number of iterations per sweep (must be > 0) | int | -1 | Yes |
| `df`   | Damping factor \f$ \omega \f$ | real | 0.9 (PETSc), must be set (native) | Yes |

\subsection doc_solvers_bjacobi Block Jacobi and Additive Schwarz

Splits the local matrix into blocks and applies a sub-solver to each block independently. The sub-solver is defined
by the child `solver` element, usually `lu` or `ilu`. Additive Schwarz (PETSc only) extends each block with
overlapping rows.

**XML:** `<pre type="block_jacobi">` or `<pre type="add_schwarz">`

Native options (\ref oft_native_solvers::oft_bjprecond "oft_bjprecond"); a child `solver` element is required.

| Option     | Description | Type | Default | Per-level |
|------------|-------------|------|---------|-----------|
| `nlocal`   | Blocks per process: `n > 0` partitions with METIS (or matrix coloring if available), `-1` one block per field of a multi-field matrix, `< -1` uses `groups` | int | 1 | Yes |
| `groups`   | Block index for each field when `nlocal < -1` (max 10 fields) | int | - | No |
| `boverlap` | Include rows shared with other processes in the local blocks | bool | F | No |
| `loverlap` | Extend each block by one layer of neighboring rows (requires `nlocal > 1`) | bool | F | No |

PETSc options (\ref oft_petsc_solvers::oft_petsc_asprecond "oft_petsc_asprecond"); if no child `solver` element is
given an LU factorization with the default package is used.

| Option    | Description | Type | Default | Per-level |
|-----------|-------------|------|---------|-----------|
| `nlocal`  | Blocks per process: `n > 1` contiguous blocks of equal size, `-1` one block per field (forces `overlap = 0`), otherwise one block | int | 0 | Yes |
| `overlap` | Number of overlapping rows | int | 0 (`block_jacobi`), 1 (`add_schwarz`) | Yes |

\section doc_solvers_direct Direct Solvers

\subsection doc_solvers_lu LU Factorization

Direct solve using an LU factorization. The factorization is computed on the first application and reused until
`update` is called.

**XML:** `<solver type="lu">` or `<pre type="lu">`

Native options (\ref oft_lu::oft_lusolver "oft_lusolver"):

| Option    | Description | Type | Default |
|-----------|-------------|------|---------|
| `package` | Factorization package: `lapack` (dense), `super` (SuperLU), `superd` (SuperLU-DIST), `umfpack`, `mumps`, `mkl` (MKL PARDISO) | str | First available of `mkl`, `umfpack`, `mumps`, `super`, `superd` |

\warning The native solver factors the process-local matrix only. In parallel runs use it as the sub-solver of
\ref doc_solvers_bjacobi.

PETSc options for `<solver type="lu">` (\ref oft_petsc_solvers::oft_petsc_direct_solver "oft_petsc_direct_solver").
`<pre type="lu">` (\ref oft_petsc_solvers::oft_petsc_luprecond "oft_petsc_luprecond") takes no options and uses the
default package.

| Option    | Description | Type | Default |
|-----------|-------------|------|---------|
| `type`    | Factorization type: `lu` or `ilu` | str | `lu` |
| `package` | Factorization package: `petsc`, `super`, `superd`, `mumps`, `umfpack` | str | First available of `mumps`, `super`, `superd`, `umfpack`, `petsc` |

\subsection doc_solvers_ilu Incomplete LU Factorization

Zero fill-in incomplete factorization, ILU(0), of the process-local matrix. Native backend only; with PETSc use
`<solver type="lu">` with `<type>ilu</type>` instead.

**XML:** `<solver type="ilu">` or `<pre type="ilu">`

| Option    | Description | Type | Default |
|-----------|-------------|------|---------|
| `package` | Factorization package: `native` or `mkl` (MKL) | str | `mkl` if available, otherwise `native` |

\section doc_solvers_mg Multigrid Preconditioners

Geometric multigrid preconditioners (\ref oft_native_solvers::oft_ml_precond "oft_ml_precond") are built over the
mesh levels with \ref oft_solver_utils::create_mlpre "create_mlpre". Each module provides default smoothers
(eg. \ref doc_input_fem), which are replaced when the module's XML element (\ref doc_solvers_xml_modules) is present:

~~~~~~~~~~~~~~~~~~~~~{.xml}
<lop>                           <!-- Element name defined by the module -->
  <smoother direction="both">
    <solver type="sjacobi">
      <its>1 2 4</its>          <!-- Values for levels 1, 2, 3 -->
      <df>0.6</df>              <!-- Single value used on all levels -->
    </solver>
  </smoother>
  <coarse>                      <!-- Optional -->
    <solver type="lu"/>
  </coarse>
</lop>
~~~~~~~~~~~~~~~~~~~~~

| Element | Description |
|---------|-------------|
| `<smoother direction="...">` | Smoother definition: `both` uses one object for the down and up sweeps; `down` and `up` define them separately and either may be omitted. At least one is required. |
| `<coarse>` | Solver used in place of the smoothers on the coarsest level (optional) |

The smoother definition is created separately on each level, and list-valued options are indexed by the mesh level
number (1 = coarsest mesh), so lists must have at least as many entries as the finest level used.
