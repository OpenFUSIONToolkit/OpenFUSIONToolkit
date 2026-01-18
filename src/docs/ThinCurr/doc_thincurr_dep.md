ThinCurr: Deprecated and legacy features     {#doc_tw_dep}
================

[TOC]

\section doc_tw_dep_dep Deprecated features
There are presently no deprecated features.

\section doc_tw_dep_leg Legacy features

\subsection doc_tw_dep_leg_holes Manual definition of Holes and Closures

In addition to the \ref doc_tw_main_holes_def "automated process", holes can also be manually defined using "nodesets" in
mesh definition files that mark the nodes corresponding to a given hole. There are two cases for these definitions:

 1. If the hole is on a boundary, as in the cylinder case, only a single node is needed and ThinCurr will determine the remaining
 elements in the hole by finding a closed boundary loop.

 2. If the hole is on the interior of the mesh, as in the torus case, all the nodes that form the hole must be provided. ThinCurr
 will then generate an ordered loop from this list.

Nodesets are stored in the `mesh/NODESET****` fields in native mesh files. These can be converted from nodesets defined via the
Cubit meshing software, through the Exodus II mesh format, and `convert_cubit.py` or other mesh inputs with suitable interfacing
scripts.

Closures can also be manually defined using "sidesets" in mesh definition files that mark the triangles where an element should be removed to
avoid singularities in the system. Sidesets are stored in the `mesh/SIDESET****` fields in native mesh files. These can be
converted from nodesets defined via the Cubit meshing software, through the Exodus II mesh format, and `convert_cubit.py` or other
mesh inputs with suitable interfacing scripts.