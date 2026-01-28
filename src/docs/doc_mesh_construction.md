Mesh Construction       {#doc_mesh_construction}
=================

[TOC]

Mesh construction in OFT consist of four phases:
  1) loading, where the initial mesh is imported and any associated geometry information is loaded and linked to the internal representation.
  2) Local construction, where the local grid is refined and the mesh graph and linkage lists are created.
  3) Decomposition, where the grid is partitioned and scattered amongst the available processors.
  4) Global linking, where seam information is assembled linking common mesh entities on adjacent domains.

\section doc_mesh_construction_load Mesh Import and Loading

The specific operations performed during import are specific to the mesh type being used. However,
the general actions performed during this phase are outlined here. The initial mesh is read from file
and loaded into the internal representation. Currently, OFT does not support domain identifiers so
any block information is discarded at import. An internal geometry representation is then constructed,
which is unique to each CAD type.

In most cases the \subpage doc_mesh_native "Native mesh format" should be used, along with included
python conversion scripts to map from other formats. Some direct legacy interfaces are also available,
but should only be used after interacting with an OFT developer.
  - \subpage doc_t3d
  - \subpage doc_cubit
  - \subpage doc_gmsh

\section doc_mesh_construction_local Local Construction

The imported mesh generally contains only a cell list which defines the tetrahedra in terms of the
vertices. Finite element methods requires counts and linkage information for the remaining entities
as well. During local construction this information is assembled by traversing the mesh in different
ways to identify and link vertices, edges, faces, and cells to one another. During this phase the
boundary mesh is also assembled and linked to volume mesh. This process is the same on all levels
regardless of decomposition, all global context is contained in seam and global indices which are
assembled and updated during global linkage construction.

\section doc_mesh_construction_decomp Decomposition

In order to accommodate distributed memory environments the mesh is partitioned using the METIS library.
This library partitions the mesh by attempting to evenly balance the number of cells in each partition
while minimizing the number of edges which link partitions. This can be think of as roughly load
balancing each domain, as the work should be proportional to the number of cells, while minimizing the
communication required, which scales with the number of elements present on more than one partition. The
partitioning is done on the master process and communicated to all other processes.

\section doc_mesh_construction_global Global Linkage

Once the grid has been partitioned seam information linking adjacent processors must be constructed
and maintained. This is done by locating boundary points on adjacent domains which have the same
global index. A consistent global index is maintained across all processors by refining the global
indices created on the base mesh, which is constructed on each process. Adjacent domains are known by
analyzing the grid partitioning to locate domains which share at least one vertex.
