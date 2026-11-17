ThickCurr: A 3D thick-wall electromagnetic simulation package     {#doc_thickcurr_main}
================

[TOC]

ThickCurr models eddy currents in 3D volumetric ("thick") conducting structures, complementing the thin-wall
approximation used in \ref doc_tw_main "ThinCurr". The magnetic vector potential \f$ A \f$, with
\f$ B = \nabla \times A \f$ and \f$ E = -\partial A / \partial t \f$, is evolved following

\f[ \frac{1}{\eta} \frac{\partial A}{\partial t} + \nabla \times \left( \frac{1}{\mu_0} \nabla \times A \right) = 0 \f]

where \f$ \eta \f$ is the resistivity of each conducting region in \f$ \Omega \cdot m \f$ and non-conducting regions
are treated as a pseudo-vacuum with a large resistivity. The domain is bounded by a perfectly conducting wall. The
vector potential is represented using the same H(Curl) + Grad(H^1) (Nedelec) finite elements as \ref doc_mhd_main "MUG".

ThickCurr supports:
 - Computation of the slowest decaying eigenmodes (L/R time constants)
 - Time-dependent simulations from a specified initial vector potential
 - Synthetic flux loops and magnetic probes, defined using the same sensor classes as ThinCurr
 - Output of the vector potential, magnetic field, electric field, and current density for visualization

See \ref thick_wall "thick_wall" for details of the implementation and \ref OpenFUSIONToolkit.ThickCurr.ThickCurr
"ThickCurr" for the Python interface, which can also be run in parallel using MPI (eg. `mpirun -np 4 python script.py`).

\section doc_thickcurr_main_ex ThickCurr Examples
The following examples illustrate usage of ThickCurr to model eddy currents in thick conducting structures.

**Python interface**
 - \subpage doc_thickcurr_slab_ex
 - \subpage doc_thickcurr_torus_ex
