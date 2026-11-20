Magnetohydrodynamics on Unstructured Grids (MUG): 3D linear/nonlinear extended MHD simulation package     {#doc_mhd_main}
================

[TOC]

\image html mug_logo.png width=25%

MUG solves the linearized or full nonlinear MHD equations in three dimensions using finite elements on
unstructured grids of tetrahedra or hexahedra.

MUG has two different implementations:
 - The primary implementation (\ref xmhd) using H(Curl) and Lagrange basis functions
 - An alternative implementation (\ref xmhd_lag) using only Lagrange basis functions and divergence cleaning

\warning Not all capability is implemented in the \ref xmhd_lag "Lagrange only" variant.

\section doc_mhd_main_ex MUG Examples
The following examples illustrate usage of MUG to run time-dependent resistive and Hall MHD simulations. 

 - \subpage doc_mug_sph_ex1
 - \subpage doc_mug_sph_ex2
 - \subpage doc_mug_recon_ex
 - \subpage doc_mug_gyro_grav_ex

\section doc_mhd_main_eqns MHD System
MUG solves the single- or two-temperature extended MHD equations: continuity, momentum (with optional
viscosity and gyroviscosity), induction with a resistive/Hall/electron-inertia Ohm's law, and
temperature evolution with anisotropic thermal conduction. The complete system advanced by the
\ref xmhd "primary implementation" is documented in:

 - \subpage doc_mug_xmhd_eqns

Most features are also supported by the \ref xmhd_lag "Lagrange only" variant.

\section doc_mhd_main_perf Performance notes
The primary implementation (\ref xmhd) evaluates the nonlinear residual and the Jacobian from a precomputed
quadrature-point cache by default (`xmhd_qcache=1` in `xmhd_options`, see \ref xmhd::xmhd_read_settings
"xmhd_read_settings" and \ref fem_quad_cache). The cache stores cell-independent reference tables of the
basis functions plus per-cell geometry, so it needs little memory (a few MB to a few tens of MB per MPI rank).
`xmhd_qcache=0` selects the original direct evaluation of the basis functions at every quadrature point, which
gives the same results to round-off. Jacobian assembly into native matrices also uses precomputed CSR positions
by default (`xmhd_jac_posmap=1`), which costs `4*sum(nce_s*nce_t)` bytes per cell and level (eg. 25 KB per cell
for tetrahedra at p=3) and can be limited with `xmhd_jac_posmap_maxmb`.

\b BLAS \b threading: the Jacobian forms its local matrices with many small BLAS `dgemm` calls from inside
OpenMP parallel regions. The BLAS library must therefore run single-threaded inside OpenMP regions, which the
OpenBLAS built by `build_libs.py` (single-threaded, `USE_THREAD=0`, by default; OpenMP, `USE_OPENMP=1`, with
`--oblas_threads`), other OpenMP builds of OpenBLAS (eg. `libopenblas0-openmp` on Debian/Ubuntu) and MKL do
automatically. With a pthreads-threaded OpenBLAS (eg. the default `libopenblas0-pthread` package on Ubuntu) set
`OPENBLAS_NUM_THREADS=1` when running MUG with OpenMP. Otherwise every small GEMM is threaded as well and the
Jacobian build stops scaling with the number of OpenMP threads.
