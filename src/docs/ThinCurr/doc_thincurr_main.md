ThinCurr: 3D thin-wall electromagnetic simulation package     {#doc_tw_main}
================

[TOC]

ThinCurr solves variants of a system of equations corresponding to inductive-resistive dynamics of currents
flowing in thin conducting structures (2D sheets in 3D geometry)

\f[ \mathrm{L} \frac{\partial I}{\partial t} + \mathrm{R} I = V. \f]

ThinCurr is used through four driver programs:

 - \ref thincurr_td for time-dependent runs
 - \ref thincurr_eig for eigenvalue calculations (\f$ \mathrm{L} \frac{\partial I}{\partial t} = \lambda \mathrm{R} I \f$)
 - \ref thincurr_fr for frequency-response (frequency-domain) runs (\f$ i \omega \mathrm{L} I + \mathrm{R} I = V \f$)
 - \ref thincurr_from_mode to convert plasma mode structures to a current potential source

\section doc_tw_main_ex ThinCurr Examples
The following examples illustrate usage of ThinCurr to perform calculations using the thin-wall model.

**Python interface (recommended)**
 - \subpage doc_tCurr_plate_eig
 - \subpage doc_tCurr_cyl_td
 - \subpage doc_tCurr_torus_mode
 - \subpage doc_tCurr_torus_fr
 - \subpage doc_tCurr_regcoil
 - \subpage doc_tCurr_hodlr
 - \subpage doc_tCurr_reduction

**Command line interface**
 - \subpage doc_thincurr_ex1
 - \subpage doc_thincurr_ex2
 - \subpage doc_thincurr_ex3
 - \subpage doc_thincurr_ex4 

\section doc_tw_main_settings ThinCurr settings groups

\subsection doc_tw_main_settings_fortran Fortran option groups
Driver-specific settings groups are defined for each of the programs above (eg. `thincurr_td`), follow links for description
of available settings.

The following driver-wide settings are also available:

**Option group:** `thincurr_hodlr_options` (see \ref doc_tw_main_hodlr "HODLR")
|  Option                 |  Description  | Type [dim] |
|-------------------------|---------------|------------|
|  `target_size=1500`     |  Target size for mesh partitioning | int |
|  `aca_min_its=20`       |  Minimum nuber of ACA+ iterations to perfom (if used) | int |
|  `L_svd_tol=-1.0`       |  SVD tolerance for HODLR compression of \f$ \textrm{L} \f$ matrix (negative to disable) | int |
|  `L_aca_rel_tol=-1.0`   |  ACA tolerance (relative to SVD) for HODLR compression of \f$ \textrm{L} \f$ matrix (negative to disable) | int |
|  `B_svd_tol=T`          |  SVD tolerance for HODLR compression of B reconstruction operator (negative to disable) | int |
|  `B_aca_rel_tol=F`      |  ACA tolerance (relative to SVD) for HODLR compression of B reconstruction operator (negative to disable) | int |

\subsection doc_tw_main_settings_xml XML input settings
Settings for ThinCurr runs are contained in the `oft->thincurr` element, with the following elements:
  * `eta`: Comma separated values for surface resistivity (\f$ \eta_s = \eta / t \f$) in each region
  * `sens_mask`: Comma separated integer mask for removing regions from sensor signals (optional; default: 0 for all regions)
  * `icoils`: Definition of coils with fixed current for \ref thincurr_td and \ref thincurr_fr (see \ref doc_tw_main_filament_def)
  * `vcoils`: Definition of coils with fixed voltage (see \ref doc_tw_main_filament_def)

```xml
<oft>
  <thincurr>
    <eta></eta>
    <sens_mask></sens_mask>
    <icoils>
      ...
    </icoils>
    <vcoils>
      ...
    </vcoils>
  </thincurr>
</oft>
```


\section doc_tw_main_num Description of numerics
Currents in ThinCurr are represented in terms of a surface current (\f$ \textbf{J}_s = \int \textbf{J} dt_w = \textbf{J} * t_w \f$),
where \f$ t_w \f$ is the direction through the thickness of the wall, represents a uniform current flowing in a wall of approximately
zero thickness. To produce a convienient basis that satisfies \f$ \nabla \cdot \textbf{J}_s = 0 \f$, a surface potential is used such that
\f$ \textbf{J}_s = \nabla \phi \times \hat{\textbf{n}} \f$, where \f$ \hat{\textbf{n}} \f$ is the unit normal on the surface. \f$ \phi \f$ is
then expanded in a scalar Lagrange basis

\f[ \textbf{J}_s (\textbf{r}) = \sum_j \alpha_j \phi_j (\textbf{r}) \times \hat{\textbf{n}} \f]

on an unstructured triangular grid. Additional filament elements can also be defined to represent coils and other suitable structures.
With these preliminaries a Boundary Finite Element Method (BFEM) is then used to define the inductance \f$ \mathrm{L} \f$ and
resistance \f$ \mathrm{R} \f$ operators as

\f[ \mathrm{L}_{i,j} = \frac{\mu_0}{4 \pi} \int_{\Omega} \int_{\Omega'} \left( \nabla \phi_i (r) \times \hat{\textbf{n}} \right) \cdot \frac{\nabla \phi_j (r') \times \hat{\textbf{n}}}{|\textbf{r}-\textbf{r}'|} dA' dA \f]

and

\f[ \mathrm{R}_{i,j} = \int_{\Omega} \eta_s(r) \left( \nabla \phi_i (r) \times \hat{\textbf{n}} \right) \cdot \left( \nabla \phi_j (r) \times \hat{\textbf{n}} \right) dA, \f]

yielding the final discretized system

\f[ \mathrm{L}_{i,j} \frac{\partial \alpha_j}{\partial t} + \mathrm{R}_{i,j} \alpha_j = V_i, \f]

where \f$ \alpha_j \f$ and \f$ V_i \f$ are the finite element weights and externally-applied voltages respectively. Voltages can be directly
applied to filament elements or through the action of defined current waveforms through inductive coupling as in \f$ \mathrm{L} \f$.

\image html thincurr_pot_ex.png "Potential (shading) and current (vectors) for the first (left) and second (right) eigenvalues of a square plate"

\subsection doc_tw_main_filament Filament elements
In addition to the meshed surfaces, conducting structures can also be represented using filaments. Filaments can have fixed current
(`icoils`) or voltage (`vcoils`) waveforms as a function of time. In the latter case, resistance per unit length and effective radius
can be specified for calculating self-inductance and resistive dissipation. For coupling between filaments and other filaments and
surface elements zero cross-sectional area is assumed (ie. a simple line integral).

Note that for `icoils` the current in a given `coil_set` is fixed, but the voltage on the coil will vary according to the
mutual coupling of all coils and surfaces. Conversely, for `vcoils` the voltage is fixed and the current evolution will vary
consistent with the mutual coupling of all coils and surfaces and the coils self-inductance and resistivity.

\subsubsection doc_tw_main_filament_def Defining filaments
The `icoils` element should contain one or more `coil_set` elements each with one or more `coil` elements.
 * Individual coil sets can be masked from sensor signals using the `sens_mask` attribute (default: 0).
 * All coils in a coil_set share the same current waveform with scale factor set by the optional `scale` attribute (default: 1.0).

If two values are given for inside a `coil` element, they are treated as the R,Z position of a circular coil. If a general 3D
coil is desired, the number of points may be specified using the `npts` attribute and a corresponding list of points (1 per line)
provided in the element with comma-separated X,Y,Z positions.
```xml
<icoils>
  <coil_set sens_mask="1">
    <coil scale="1.0">1.0, 0.1</coil>
    <coil npts="5" scale="1.0">1.0, 0.0, -0.1
      0.0, 1.0, -0.1
      -1.0, 0.0, -0.1
      0.0, -1.0, -0.1
      1.0, 0.0, -0.1
    </coil>
  </coil_set>
</icoils>
```

The `vcoils` element should contain one or more `coil_set` elements each with one or more `coil` elements.
 * Individual coil sets can be masked from sensor signals using the `sens_mask` attribute (default: 0).
 * All coils in a coil_set share the same current waveform with scale factor set by the optional `scale` attribute (default: 1.0).
 * The `res_per_len` attribute is used to specify the resistance per unit length (\f$ \Omega/m \f$) of the filament
 * The `radius` attribute is used to set the effective radius of the conductor for calculation of the self-inductance

As above, if two values are given for inside a `coil` element, they are treated as the R,Z position of a circular coil. If a general 3D
coil is desired, the number of points may be specified using the `npts` attribute and a corresponding list of points (1 per line)
provided in the element with comma-separated X,Y,Z positions.
```xml
<vcoils>
  <coil_set sens_mask="1">
    <coil res_per_len="2.195E-4" radius="0.005" scale="1.0">1.1, 0.25</coil>
    <coil res_per_len="2.195E-4" radius="0.005" npts="5" scale="1.0">1.0, 0.0, -0.1
      0.0, 1.0, -0.1
      -1.0, 0.0, -0.1
      0.0, -1.0, -0.1
      1.0, 0.0, -0.1
    </coil>
  </coil_set>
</vcoils>
```


\subsection doc_tw_main_bc Boundary conditions
At the edge of a surface a boundary condition is required to ensure that no current flows normal to the boundary, which would
result in divergence or lost current. By considering the potential representation we can easily show that zero normal current is
satisfied by the condition \f$ | \nabla \phi \times \hat{\textbf{n}} \times \hat{\textbf{t}} | = | \nabla \phi \cdot \hat{\textbf{t}} | = 0\f$,
where \f$ \hat{\textbf{t}} \f$ is the unit normal in the direction tangential to the boundary. This amounts to the condition that
\f$ \phi \f$ is a constant on a given boundary. A simple solution that satisfies this boundary condition is to set the
value on all boundary nodes to a single fixed value (generally zero). However as we will see in the next section this is overly
restrictive for many cases.

\image html thincurr_elem_ex.png "Example currents associated with nodes on the boundary (red) and interior (blue) showing a non-zero divergence at the edge without boundary conditions. The boundary conditions ensure no net current flows normal to the surfaces as in the red element."

\subsection doc_tw_main_holes "Hole" elements
Note the boundary condition above requires a constant potential on a given boundary. However, for multiply connected geometry, such
as a cylinder, one or more distinct boundaries exist such that the potential can be different on each boundary. For the case of
the cylinder there are two geometrically distinct boundarys, the top and bottom circular edges. If we consider the potential formulation
we can show that the current flow across a line between two points on a given surface is given by the difference in the potential between
those points

\f[ \mathrm{I}_{l} = \int_a^b \hat{\textbf{n}} \cdot (\nabla \phi \times \hat{\textbf{n}} \times \textbf{dl}) = \int_a^b \nabla \phi \cdot \textbf{dl} = (\phi_b - \phi_a). \f]

So if point \f$ a \f$ and \f$ b \f$ exist on different boundaries then the difference in the constant for each boundary
defines the current passing between them. In the case of the cylinder this corresponds to the total current flowing azimuthally around
the cylinder.

To enable the values on boundaries to vary self-consistently in the model we define a new element that corresponds to a constant potential
on a given closed boundary loop. In ThinCurr we call these elements "holes", as they are related to the topological concept of "holes".
These elements correspond to current flowing around the boundary loop, and correspondingly to capturing magnetic flux within this loop.
This captures flux that links the mesh without crossing through any of the triangles themselves (eg. flux down the center of a cylinder).

\image html thincurr_cyl_hole.png "Example current associated with a hole element from the upper edge of a cylinder"

While one may think we should add two holes in the case of the cylinder this is not the case as it would introduce a
[gauge ambiguity](https://en.wikipedia.org/wiki/Gauge_fixing) in the system, resulting in a singular \f$ \textbf{L} \f$ matrix. We can see
this by considering the azimuthal flowing current, which as above is given by the difference in potential between the top \f$ \phi_a \f$
and bottom \f$ \phi_b \f$ of the cylinder. This means that a single current in our model is dependent on two values, indicating a redundancy
that means \f$ \textbf{L} \f$ is not full rank, or in other words singular. To avoid this we only introduce only one hole, which in turn sets
the potential on the other boundary to zero and eliminates the gauge ambiguity.

\image html thincurr_torus_holes.png "Example currents associated with the hole elements for the poloidal (blue) and toroidal (red) directions"

**Multiply connected regions**:

Hole elements are not just necessary when there is a visible boundary, but whenever there are multiple distinct topological paths on
a given surface. An example of this case is the torus, where there are two topologically distinct paths corresponding to the short (poloidal)
and long (toroidal) way around the torus. We can see this using the consideration of the difference in potentials above, as the potential
difference between any two points will always be zero and as a result no net current can flow in the poloidal or toroidal directions. To
correct this two holes must be added, corresponding to loops in each of these two directions, which makes \f$ \phi \f$
[mutlivalued](https://en.wikipedia.org/wiki/Multivalued_function) so that a difference in potential can exist even on closed paths.
This use of jumps and a single-valued potential to represent a multivalued potential is common practice in numerical methods.

\subsubsection doc_tw_main_holes_def Defining holes
Holes can be defined using the `ThinCurr_compute_holes.py` script. This script analyzes the toplogy of a given
mesh and automatically locates and defines needed hole elements. This is now the recommended way of defining holes in
ThinCurr models, although \ref doc_tw_dep_leg_holes "manual definition" is also still supported.

\subsection doc_tw_main_close "Closure" elements
As stated above the fact that the solution depends only on the gradient and not on the absolute value of \f$ \phi \f$ itself can
introduce a gauge ambiguity, resulting in redundant degrees of freedom. In the cases considered above this is resolved by keeping
the potential on at least one boundary fixed to zero. However, when a surface is fully closed, as is the case for a sphere or torus,
there are no boundaries. So to fix the gauge on these surfaces we must remove an element from the system, effectively fixing the
potential to zero at that point. We call these elements "closure" elements as they are used to close the system, making it
solveable.

\subsubsection doc_tw_main_close_def Defining closures
As with holes, closures can also be defined using the `ThinCurr_compute_holes.py` script. This is now the recommended way of defining
closures in ThinCurr models, although \ref doc_tw_dep_leg_holes "manual definition" is also still supported.

\subsection doc_tw_main_hodlr Hierarchical Off-Diagonal Low-Rank (HODLR) approximation
In the standard BFEM approach described above every element interacts with every other element, leading to a dense \f$ \textrm{L} \f$ and
growth of memory and time requirements at a rate of \f$ O(N^2) \f$, where \f$ N \f$ is the number of elements in the model. This limits
the size of models on most systems to < 20k elements on a typical laptop with 16GB of RAM and < 180k on even a large workstation with
1TB of RAM. In addition to memory limits, these models will take a very long time to compute as the number of FLOPs required will also
scale as \f$ O(N^2) \f$.

To avoid this ThinCurr models can take advantage of the structure of the \f$ \textrm{L} \f$ matrix to generate an approximation of the
full matrix with some desired tolerance. As application of these models to engineering and other studies often do not require extremely
high accuracy, this can dramatically reduce the computational requirement while retaining sufficient fidelity for a desired study. This
is possible as the \f$ \textrm{L} \f$ constructed from the BFEM exhibits a low-rank structure for interactions between elements that are
physically separated by a large distance relative to their size. If we spatial partition our mesh into regions, the interactions between
separate regions, so called off-diagonal blocks of the matrix, will possess this low-rank structure, meaning they can be well
approximated by a matrix with smaller rank (eg. truncated SVD). To improve efficiency and performance such a spatial partitioning
can be done in a nested fashion, leading to a so-called Hierarchical Off-Diagonal Low-Rank (HODLR) approximation to the matrix.

\image html thincurr_aca_scaling.png "Improvement in required memory (top) and time (bottom) for a eigenvalue calculation using ACA+, compared to the direct approach"

\subsection doc_tw_main_hodlr_part Spatial partitioning
ThinCurr currently utilizes a binary tree to partition the surface mesh used to discretize currents in the model. Initially, a
bounding box is created that encloses the full mesh. The mesh is then recursively subdivided by looking at each partition and
subdividing it in one of the three principal cartesian directions if the number of elements is greater than a given size
(`target_size`). On each level, the direction of subdivision is chosen to be the direction with largest standard deviation of
position over the contained points. Near and far interactins between mesh regions are then determined on each level, taking
into account appropriate "masking" of interactions from higher levels.

\image html thincurr_aca_part.png "Example mesh regions following partition of a ThinCurr mesh for HODLR"

\subsection doc_tw_main_hodlr_aca Adaptive Cross-Approximation (ACA+)
Once the spatial partitioning is complete the individual matrix blocks can be constructed. First, diagonal blocks are fully assembled
and stored as dense matrices. Next, low-rank approximations are computed for all off-diagonal blocks. For blocks whose corresponding
mesh regions are "close" to eachother, as determined by the ratio of the center-center separation to the circumradius of
each region, the full block is computed and then compressed using an SVD. The compression is set by truncating the SVD
at a user-specified tolerance (specified by `L_svd_tol`), based on the Frobenius norm of the singular values.

For all other blocks, which correspond to off-diagonal blocks whose corresponding mesh regions are not "close" to eachother
the Adaptive Cross-Appoximation technique is used. This technique iteratively builds a low-rank approximation
(\f$ \textrm{M} \approx \textrm{U} \textrm{V}^T \f$) by sampling rows and columns of the block and stopping once a desired
tolerance, defined relative to the SVD tolerance (`L_aca_rel_tol < 1.0`), has been reached. As ACA+ has some inherent
random variation, the tolerance specified for ACA+ should be higher than the SVD tolerance (`L_aca_rel_tol < 1.0`).
After a block's approximation by ACA+ is complete the resulting block is then recompressed using SVD for greater run
to run consistency.

\note The same techniques can also be applied to the magnetic field reconstruction operator, with corresponding settings of
`B_svd_tol` and `B_aca_rel_tol`. However, as these tolerances are absolute the magnetic versions should generally be 4-5 orders of
magnitude larger than those for the \f$ \textrm{L} \f$ matrix.