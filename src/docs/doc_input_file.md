Input File Settings    {#doc_input_file}
===================

[TOC]

# Introduction {#doc_input_intro}

The Open FUSION Toolkit uses a input file to provide run-time information about different parameters, ex. the input mesh and number of refinements. By default OFT looks uses the file `oft.in` in the working directory, a different input file may be specified by passing its name as the first argument.

\verbatim
~$ oft_check_mesh myinput
\endverbatim

 Runtime options which are common across different mesh and physics types are outlined below. For further options specific to a given
 mesh type consult their documentation, which may be accessed from the main page. For options specific to a given physics module and
 driver consult the documentation for that module.

Physics modules with settings groups:
 - \ref xmhd
 - \ref xmhd_lag

# Runtime Settings {#doc_input_runtime}

 This group contains options which are common to the basic runtime environment for OFT. The `debug` variable can be used to control the
 level of output produced by debugging statements. Higher values for this variable indicated more verbose status output during execution. The
 `stack_disabled` variable can be used in conjuction with OFT's \ref doc_stack functionality to disable profiling until a specific point in
 execution. This allows you to exclude setup code from profiling and provide results only on the operations of interest.

 The `ppn` variable is used to indicated the number of processors on a given NUMA node, who share a common and lower latency communication
 link. Setting this to a value above 1 causes the decomposition to take two steps: First, the mesh is decomposed into `nprocs/ppn` domains
 using METIS. Second each domain is further decomposed into `ppn` subdomains, which are then distributed. In principle this could provide
 better communication balancing by ignoring the lower latency intra-node communication, however in practice no improvement is seen.

**Option group:** `runtime_options`
|  Option  |  Description  | Type |
|------------|------------|-----|
| `ppn=1`             | Number of processors per NUMA node | int |
| `debug=0`           | Debugging level (0-3) | int |
| `stack_disabled=F`  | Disable stack functionality | bool |
| `use_petsc=F`       | Use PETSc linear algebra backend | bool |

# Mesh settings {#doc_input_mesh}

This group contains options which are common to the OFT hierarchical mesh environment. The settings control the type of mesh, and resulting
CAD respresentation, used as well as the number of multi-grid levels to be constructed and refinement parameters. Additional mesh settings specific
to a given CAD representation are specified as part of their documentation: \ref doc_t3d, \ref doc_cubit

The number of mesh levels can be set using the `nlevels` variable. This variable sets the total number of grid levels including the transfer
level which is created by domain decomposition. The variable `nbase` is used in conjunction with this setting to create a certain number of grid
refinements before parallel decomposition is performed. For example using `nlevels = 2` and `nbase = 1` will produce 2 grid levels with no refinement,
one shared memory level of the original mesh and one level with a distributed mesh.  

\note OFT distinguishes between shared memory and distributed levels using the `nbase` setting. A run using multiple MPI tasks must have
`nlevels > nbase`, otherwise an error will be thrown at run time. For shared memory runs `nbase` and `nlevels` should always be
equal.

Two additional variables control the use of CAD geometry during grid construction. The `grid_order` variable sets the order of the geometry
representation on the finest grid level, coarse grid levels always use linear tetrahedra. This allows the representation of curved surfaces in
OFT with quadratic or cubic tetrahedra. The `fix_boundary` variable can be used to enable/disable boundary adjustment during grid refinement.
By default OFT uses CAD information to place new mesh vertices on the true surface during grid refinement, when `fix_boundary = false` vertices
created by edge division will remain at the mid point of the coarse level edge resulting in a boundary which matches the coarse mesh.

\note The highest grid order supported may be dependent on the mesh type used, see the CAD type detailed documentation for more information.

**Option group:** `mesh_options`
|  Option  |  Description  | Type |
|------------|------------|-----|
| `meshname=""`     | Name of mesh for use with I/O | str(40) |
| `cad_type=1`      | Type of mesh for import (1 = T3D, 2 = Exodus, 3 = GMSH, 91 = Sphere, 92 = Cube) | int |
| `nlevels=2`       | Total number of MG levels | int |
| `nbase=1`         | Number of base refinements before decomposition to MPI | int |
| `grid_order=1`    | Order of finest level tetrahedra (1-2) | int |
| `fix_boundary=T`  | Adjust new boundary points to CAD boundary | bool |

# Finite element smoother settings {#doc_input_fem}

Currently, OFT uses 3 additional settings groups to specify the smoother settings for common multi-grid preconditioners used throughout
the code. Going forward these blocks may be reorganized or removed as the code moves away from its initial usage toward a more general framwork.

@note The presence or structure of these settings blocks may change in the future.
Make sure you use the current documentation as a reference.

## Lagrange

This group is used to specify the common smoother settings for common Lagrange scalar
operators. Currently the isotropic and anisotropic versions of the Laplacian operator
are supported.

**Option group:** `lag_op_options`
|  Option  |  Description  | Type |
|------------|------------|-----|
| `df_lop=1`   | Smoothing factors for LAG::LOP | real [nlevels] |
| `nu_lop=1`   | Number of smoother iterations for LAG::LOP | int [nlevels] |
| `df_pdop=1`  | Smoothing factors for LAG::PDOP | real [nlevels] |
| `nu_pdop=1`  | Number of smoother iterations for LAG::PDOP | int [nlevels] |

## H^1

This group is used to specify the common smoother settings for common H^1 scalar
operators. Currently only the isotropic Laplacian operator is supported.

**Option group:** `h1_op_options`
|  Option  |  Description  | Type |
|------------|------------|-----|
| `df_lop=1`   | Smoothing factors for H^1::LOP | real [nlevels] |
| `nu_lop=1`   | Number of smoother iterations for H^1::LOP | int [nlevels] |

## H(Curl)

This group is used to specify the common smoother settings for common H(Curl) vector
operators. Currently only the \f$ \nabla \times \nabla \times \f$ operator is supported.

**Option group:** `hcurl_op_options`
|  Option  |  Description  | Type |
|------------|------------|-----|
| `df_wop=1`   | Smoothing factors for H(Curl)::WOP | real [nlevels] |
| `nu_wop=1`   | Number of smoother iterations for H(Curl)::WOP | int [nlevels] |
