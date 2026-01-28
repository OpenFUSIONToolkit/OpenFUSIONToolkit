PETSc usage  {#doc_petsc_usage}
===========

[TOC]

This page outlines information helpful for using the Open FUSION Toolkit with PETSc as the linear algebra back end.

### Run-Time options

* `-threadcomm_nthreads` Number of threads to use per MPI task
* `-threadcomm_type` Type of threads to use, generally this will be set to `openmp`

\section doc_petsc_usage_error Error Handling

By default PETSc either returns errors via `ierr` of the calling subroutine or prints an error to
`stderr`. This is not uniformly implemented and can cause issues with large output files.

### Run-Time options

* `-on_error_abort` Make PETSc internal errors fatal

\section doc_petsc_usage_threads Threading (not recommended)

In order to use threading PETSc must be compiled with threads enabled. Currently, this is only correctly
supported using the PETSc-dev branch. In order to compile with thread capabilities add the following
commands to the configure script.

\verbatim
--with-openmp --with-pthreadclasses --with-threadcomm
\endverbatim
