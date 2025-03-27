CUBIT Mesh Interface    {#doc_cubit}
====================

[TOC]

The Open FUSION Toolkit (OFT) includes an interface to the CUBIT mesh generation program. This interface allows the
import and refinement of meshes created in CUBIT. Mesh import requires that OFT is built
with NETCDF. Boundary conformal refinement requires the OpenNURBS library as well.

\section cubit_generation Generating Meshes

In order to construct the boundary representaion in OFT required for refinement associated
CAD information must be output along with the mesh. This is accomplished via the `refine parallel`
command in CUBIT. Once the mesh has been created in CUBIT use the command line to export the mesh
and CAD data with the `refine parallel` command, where "base_name" is the root name for the output files.

\verbatim
refine parallel fileroot "base_name" overwrite no_execute
\endverbatim

This command will generate two files: a mesh file `base_name.in.e` and a CAD file `base_name.3dm`
which are used by OFT.

\note Currently `refine parallel` functions correctly with CUBIT 14+ only.

\section cubit_input Input Options

In order to use a CUBIT mesh `cad_type = 2` must be specified in the `mesh_options` group. The
mesh files and some additional options are specified at runtime in the `cubit_options` group.
The options available in this group are listed below.

**Option group:** `cubit_options`
|  Option  |  Description  | Type [dim] |
|------------|------------|-----|
| `filename=""`  | Name of mesh file | str(OFT_PATH_SLEN) |
| `inpname=""`   | Name of OpenNURBS geometry file | str(OFT_PATH_SLEN) |
| `lf_file=T`    | Large file format flag | bool |
| `reflect=F`    | Reflect mesh about xy-plane | bool |
| `per_ns=-1`    | Nodeset ID for periodic surfaces | int |
| `zstretch=1.`  | Scale for z-direction | float |

\section cubit_refinement Refinement Notes

For CUBIT meshes new points are placed by minimizing the weighted sum of the distances between the
new point and its parent points, subject to the constraint that the point exist on the boundary
surface. For example if a new point is to be placed at the midpoint of a boundary edge, the process
is as follows. The parent CAD entity for that edge is first determined, either curve or surface. A
point is the found by minimizing the sum of the distance from each end point to the new point by
adjusting its position in parametric space on this CAD entity. Minimization is performed using a
Levenburg-Marquardt algorith provided by the `MINPACK` library.
