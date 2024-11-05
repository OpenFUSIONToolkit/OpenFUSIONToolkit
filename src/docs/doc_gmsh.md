GMSH Mesh Interface    {#doc_gmsh}
===================

[TOC]

The Open FUSION Toolkit includes an interface to the free [GMSH](https://www.geuz.org/gmsh) mesh generation program, pre-built binary versions of
this software are available for most operating systems. Currently this interface only supports import and non-conformal refinement
and is designed to allow the use of OFT when one of the other proprietary mesh programs is not available. In the future we
hope to update this interface to support some form of boundary conformal refinement, most likely based on a quadratic input mesh.

\section doc_gmsh_generation Generating Meshes

OFT uses files output in GMSH's `mesh` format. An example is shown below for generating a simple mesh from a CAD file `torus.stp` using
the command line interface and saving it in the file `torus.mesh`.

\verbatim
gmsh -3 -format mesh -optimize -clscale .1 -order 2 -o torus.mesh torus.stp
\endverbatim

This generates a mesh of tetrahedra with a quadratic boundary representation and a basic mesh size of `0.1`. For most simple cases these
are the only options that are required with changes to the desired element size and order only. GMSH has a number of mesh generation
options, for more information see the GMSH documentation.

\section doc_gmsh_input Input Options

In order to use a GMSH mesh `cad_type = 3` must be specified in the `mesh_options` group. The
mesh files and some additional options are specified at runtime in the `gmsh_options` group.
The options available in this group are listed below.

**Option group:** `gmsh_options`
|  Option  |  Description  | Type [dim] |
|------------|------------|-----|
| `filename=""` | Name of mesh file | str(40) |
| `order=1`     | Order of input mesh (not yet used) | int |
