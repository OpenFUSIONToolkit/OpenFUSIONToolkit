T3D Mesh Interface    {#doc_t3d}
==================

[TOC]

The Open FUSION Toolkit includes an interface to the T3D mesh generation program. This interface allows the
import and refinement of meshes created with T3D.

\section doc_t3d_generation Generating Meshes

Docs needed

\section doc_t3d_input Input Options

In order to use a T3D mesh `cad_type = 1` must be specified in the `mesh_options` group. The
mesh files and some additional options are specified at runtime in the `t3d_options` group.
The options available in this group are listed below.

**Option group:** `t3d_options`
|  Option  |  Description  | Type [dim] |
|------------|------------|-----|
| `filename=""`    | Name of mesh file | str(40) |
| `inpname=""`     | Name of T3D input file | str(40) |
| `reflect=""`     | Reflection flag (eg. 'xy') | str(3) |
| `ref_per=F,F,F`  | Periodicity flags | bool [3] |
| `zstretch=1.`    | Scale for z-direction | float |

\section doc_t3d_input_reflect Mesh Reflection

Docs needed

\section doc_t3d_input_period Mesh Periodicity

Docs needed

\section doc_t3d_refinement Refinement Notes
