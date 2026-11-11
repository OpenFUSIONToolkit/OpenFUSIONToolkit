Linear Algebra Backends     {#doc_solvers_backends}
=======================

[TOC]

OFT supports two linear algebra backends: a native implementation and [PETSc](https://petsc.org). The backend
determines the type of the vectors and matrices created by OFT, and the implementation used when a solver is
created with \ref oft_solver_utils or from XML.

\section doc_solvers_backends_select Selecting a Backend

The backend is selected with the `use_petsc` runtime option (\ref doc_input_runtime), which requires OFT to be built
with PETSc. The native backend is used by default.

The available methods and their classes for each backend are listed in \ref doc_solvers_xml_types.

\section doc_solvers_backends_mixing Mixing Native and PETSc Methods

Native solvers use generic linear algebra calls and work with either backend. PETSc solvers require PETSc
preconditioners, which is always the case when building from XML. A native solver (`native="true"`) may use the PETSc
`jacobi`, `block_jacobi` and `add_schwarz` preconditioners, but not `lu` or `ilu`.

The `native="true"` XML attribute (or the `force_native` argument of the `create_*` routines) can be used to select
a native solver when PETSc is in use.

The native \ref doc_solvers_bjacobi "block Jacobi", \ref doc_solvers_lu "LU" and \ref doc_solvers_ilu "ILU"
methods require native matrices and are only available with the native backend.

\section doc_solvers_backends_petsc PETSc Backend

PETSc solvers are configured through the same options as native solvers, with the following differences:
 - `its > 0` sets the iteration limit, otherwise the PETSc default limit is used (`its = 0` skips GMRES solves).
   The special values of `its` for native solvers (\ref doc_solvers_common) are not supported.
 - Convergence information, when enabled (see `pm` in \ref doc_solvers_common), is printed by the PETSc residual monitor.

\subsection doc_solvers_backends_petsc_opts Run-Time Options

PETSc reads additional [run-time options](https://petsc.org/release/manual/other/#runtime-options) from the
command line, the `PETSC_OPTIONS` environment variable or a `.petscrc` file. These options do not change solvers
created by OFT, but can be used for diagnostics. For example:

 - `-on_error_abort`: Abort on PETSc errors. By default errors are returned through `ierr` or printed to `stderr`,
   which is not handled uniformly by OFT and can produce large output files.
