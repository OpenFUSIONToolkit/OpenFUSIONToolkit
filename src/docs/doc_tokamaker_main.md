TokaMaker: 2D static and time-dependent Grad-Shafranov Equilibria     {#doc_gs_main}
================

[TOC]

TokaMaker solves variants of a system of equations corresponding to Ideal MHD force balance
for an axisymmetric plasma in response to actively and passively driven currents in nearby
conductors (eg. coils and structures).

\f[
\Delta^* \psi = 
\begin{cases}
    -\frac{1}{2}\frac{\partial F^2}{\partial \psi} - \mu_0 R^2 \frac{\partial P}{\partial \psi} & \text{if } \textbf{r} \in \mathcal{P}\\
    -R \mu_0 J_{\phi} & \text{if } \textbf{r} \in \mathcal{S},\mathcal{C} \\
    0 & \text{elsewhere;}
\end{cases}
\f]
where \f$ \mathcal{P} \f$, \f$ \mathcal{S} \f$, and \f$ \mathcal{C} \f$ are axisymmetric domains corresponding to the plasma, passive conducting structures (eg. vacuum vessels), coils respectively.

TokaMaker should primarily be used through the python interface using the \ref OpenFUSIONToolkit.TokaMaker "OpenFUSIONToolkit.TokaMaker" python module
and the \ref OpenFUSIONToolkit.TokaMaker._core.TokaMaker "OpenFUSIONToolkit.TokaMaker.TokaMaker" class.

\section doc_gs_main_ex TokaMaker examples
The following examples illustrate usage of TokaMaker to compute different Grad-Shafranov equilibria. For examples of how to create
new meshes, see \ref doc_gs_main_mesh_ex.

### Fixed Boundary Equilibria
 - \subpage doc_tMaker_fixed_ex1
 - \subpage doc_tMaker_fixed_ex2

### Free Boundary Equilibria
 - \subpage doc_tMaker_ITER_ex2
 - \subpage doc_tMaker_HBT_ex2
 - \subpage doc_tMaker_ITER_ex2
 - \subpage doc_tMaker_DIIID_ex2
 - \subpage doc_tMaker_CUTE_ex2
 - \subpage doc_tMaker_LTX_ex2
 - \subpage doc_tMaker_ITER_ex3

### Equilibrium Reconstruction
 - \subpage doc_tMaker_ITER_ex4

### Other workflows/emaples
 - \subpage doc_tMaker_ITER_ex5

\section doc_gs_main_mesh Building meshes using gs_Domain
TokaMaker includes built-in meshing functionality through the \ref OpenFUSIONToolkit.TokaMaker.meshing.gs_Domain "gs_Domain" class,
which leverages the [triangle](https://pypi.org/project/triangle/) python package to generated unstructured triangular grids.
The process of geometry definition is broken into two stages:

First, logical information about the geometry is defined using \ref OpenFUSIONToolkit.TokaMaker.gs_Domain.define_region "define_region()".
This method's functionality varies based on the type of region being defined, which must be one of the following:
 * `plasma`: The region where the plasma can exist and the classic Grad-Shafranov equation with \f$ F*F' \f$ and \f$ P' \f$ are allowed. **There can only be one region of this type**
 * `vacuum`: A region where no current can flow and \f$ \nabla^* \psi = 0 \f$ is solved
 * `boundary`: A special case of the `vacuum` region, which forms the outer boundary of the computational domain. **A region of this type is required if more than one region is specified**
 * `conductor`: A region where toroidal current can flow passively (no externally applied voltage). For this type of region the resistivity should be specified with the argument `eta` in units of \f$ \omega \mathrm{-m} \f$.
 * `coil`: A region where toroidal current can flow with specified amplitude through \ref OpenFUSIONToolkit.TokaMaker.TokaMaker.set_coil_currents "set_coil_currents()" or via shape optimization \ref OpenFUSIONToolkit.TokaMaker.TokaMaker.set_coil_reg "set_coil_reg()" and \ref OpenFUSIONToolkit.TokaMaker.TokaMaker.set_isoflux "set_isoflux()"

Next, geometric information about each region defined in the prior step is performed using one of the following methods:
 * \ref OpenFUSIONToolkit.TokaMaker.meshing.gs_Domain.add_rectangle "add_rectangle()" defines the region as a simple rectangle with a specified center, width, height, and optional rotation.
 * \ref OpenFUSIONToolkit.TokaMaker.meshing.gs_Domain.add_polygon "add_polygon()" defines the region as any closed polygon from a list of (R,Z) points.
 * \ref OpenFUSIONToolkit.TokaMaker.meshing.gs_Domain.add_annulus "add_annulus()" defines the region as an annulus between two closed polygons that must not cross.
 * \ref OpenFUSIONToolkit.TokaMaker.meshing.gs_Domain.add_enclosed "add_enclosed()" defines the region by specifying a point within a space enclosed by other regions (eg. between a VV and a limiter). This is useful to avoid having to define a new curve(s) in complex configurations.

\subsection doc_gs_main_mesh_ex Mesh generation examples using gs_Domain
The following examples illustrate usage of \ref OpenFUSIONToolkit.TokaMaker.meshing.gs_Domain "gs_Domain" to create meshes for TokaMaker.

 - \subpage doc_tMaker_HBT_ex1
 - \subpage doc_tMaker_ITER_ex1
 - \subpage doc_tMaker_DIIID_ex1
 - \subpage doc_tMaker_CUTE_ex1
 - \subpage doc_tMaker_LTX_ex1
