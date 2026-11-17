TokaMaker: 2D static and time-dependent Grad-Shafranov Equilibria     {#doc_gs_main}
================

[TOC]

\image html tokamaker_logo.png width=25%

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
and the \ref OpenFUSIONToolkit.TokaMaker._core.TokaMaker "OpenFUSIONToolkit.TokaMaker.TokaMaker" and \ref OpenFUSIONToolkit.TokaMaker._core.TokaMaker_equilibrium "OpenFUSIONToolkit.TokaMaker.TokaMaker_equilibrium" classes.

\section doc_gs_main_ex TokaMaker examples
The following examples illustrate usage of TokaMaker to compute different Grad-Shafranov equilibria. For examples of how to create
new meshes, see \ref doc_gs_main_mesh_ex.

### Fixed Boundary Equilibria
 - \subpage doc_tMaker_fixed_ex1
 - \subpage doc_tMaker_fixed_ex2

### Free Boundary Equilibria
 - \subpage doc_tMaker_ITER_ex2
 - \subpage doc_tMaker_HBT_ex2
 - \subpage doc_tMaker_DIIID_ex2
 - \subpage doc_tMaker_DIIID_ex4
 - \subpage doc_tMaker_CUTE_ex2
 - \subpage doc_tMaker_LTX_ex2
 - \subpage doc_tMaker_ITER_ex3
 - \subpage doc_tMaker_MANTA_ex2
 - \subpage doc_tMaker_ARC_ex2
 - \subpage doc_tMaker_SPARC_ex2

### Equilibrium Reconstruction
 - \subpage doc_tMaker_ITER_ex4

### Dipole Equilibria
 - \subpage doc_tMaker_dipole_ex2

### Mirror Equilibria
 - \subpage doc_tMaker_wham_ex2

### Other workflows/examples
 - \subpage doc_tMaker_ITER_ex5
 - \subpage doc_tMaker_HBT_ex3
 - \subpage doc_tMaker_CUTE_ex3
 - \subpage doc_tMaker_CUTE_ex4
 - \subpage doc_tMaker_ARC_ex3
 - \subpage doc_tMaker_GPEC_ex1
 - \subpage doc_tMaker_blanket_vde_ex
 - \subpage doc_tMaker_DIIID_ex3

### Deprecated workflows
 - \subpage doc_tMaker_ITER_ex4-legacy

\section doc_gs_main_conv Conventions
 - \subpage doc_tokamaker_current_conventions

\section doc_gs_main_lcfs Plasma integrals at the boundary (cut-cell quadrature)
By default the plasma source and the plasma integrals (current, stored energy, volume, ...) are evaluated at the
quadrature points of each cell in the plasma region, and a point contributes only if it lies inside the LCFS
(\f$ \psi > \psi_{LCFS} \f$ and on the plasma side of the X-point half-planes). When the pressure gradient or the
current density is finite at the LCFS (for example pedestal profiles from reconstructions), the assembled source then
changes by a finite amount whenever a quadrature point crosses the boundary between two nonlinear iterations, and the
Picard iteration can end in a limit cycle (period 2 or 3, often with a nearly constant residual) instead of converging.

Setting `lcfs_cut_cells = True` (\ref OpenFUSIONToolkit.TokaMaker._core.tokamaker_settings "TokaMaker settings")
replaces this point test by a cut-cell quadrature. In each plasma-region cell the plasma part is the part where the
linear interpolant of \f$ \psi - \psi_{LCFS} \f$ on the cell's vertices is positive (also restricted by the X-point
half-planes); the cell's quadrature rule is mapped onto triangles covering that part, and the flux functions are
evaluated at \f$ \max(\psi, \psi_{LCFS}) \f$. A cell wholly inside the LCFS gets its standard rule, which the clipped
rule reduces to exactly, so these integrals are continuous functions of the nodal \f$ \psi \f$, of \f$ \psi_{LCFS} \f$
and of the X-point locations. The same quadrature is used for the right-hand side and the target rows of the solve,
the reported global quantities (`get_stats`, `get_globals`), the plasma current and diamagnetic flux fit constraints,
the `jphi-linterp` current normalization, the bootstrap current integrals, the plasma-coil mutual inductances and the
loop voltage, so the reported plasma current equals the current of the solved source.

The continuum problem is unchanged; only the integration over the cells crossed by the LCFS differs, by the mesh's own
boundary discretization error. Converged results therefore change slightly when the setting is switched on (on the
meshes tested 1e-5 to 1e-3 relative in stored energy and q, up to about 1e-2 in q95 on coarse boundary meshes), which
is why it is off by default. Not changed: the time-dependent and linear-stability operators (`setup_td`,
`compute_linear_stability`), point evaluations of fields and profiles (`get_field_eval`, `get_jtor_plasma`) and the
pointwise pressure fit constraint. Changes of the plasma topology remain discontinuous with either quadrature: an
X-point appearing or a change of the active X-point, a switch between diverted and limited plasmas, and the choice of
the O-point between two maxima of nearly equal flux.

```
mygs.settings.lcfs_cut_cells = True
mygs.update_settings()
```

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
 - \subpage doc_tMaker_MANTA_ex1
 - \subpage doc_tMaker_LTX_ex1
 - \subpage doc_tMaker_dipole_ex1
 - \subpage doc_tMaker_wham_ex1
 - \subpage doc_tMaker_ARC_ex1
 - \subpage doc_tMaker_SPARC_ex1
