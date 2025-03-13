!!Taylor Example 1: Force-free eigenmodes    {#doc_tay_ex1}
!!========================
!!
!![TOC]
!!
!!This example introduces the \ref taylor "Taylor State" physics module. The \ref taylor "Taylor"
!!module is used to compute Ideal MHD force-free states with uniform \f$ \lambda \f$. This example
!!computes the Taylor state in a oblate cylinder. The T3D interface is used to generate a mesh for
!!the cylinder and provide boundary information for the setup of a quadratic spatial mapping. This
!!example also demonstrates the steps required to \ref doc_ex3_code_project "output a vector field"
!!for plotting.
!!
!!\section doc_ex3_t3d Mesh Creation with T3D
!!
!!The T3D mesh generation code creates tetrahedral meshes from a text file describing the boundary
!!with rational Bezier splines. For this example we will demonstrate using this code to generate a
!!mesh for an oblate, tuna can, cylinder.
!!
!!\subsection doc_ex3_t3d_geom Geometry Definition
!!
!!The geometry files consists of heirarchy of geometry primitives vertices, curves, and surfaces,
!!which are defined using a defined grammar. The basics of the geometry file used are outlined here
!!for more information about the T3D input files consult the T3D manual.
!!
!!To simplify the geometry file we define only a single quadrant of a cylinder and use the \ref
!!doc_t3d_input_reflect "reflection" feature of the Open FUSION Toolkit (OFT) T3D interface to produce the full mesh.
!!
!!First we define a set of vertices that will form the 6 corners of the quadrant. Notice that 8
!!vertices are actually defined, with 2 copy vertices, to represent the degenerate curve on the
!!symmetry axis. This curve is need to define the top and bottom surfaces in the logically
!!rectangular RBS definition.
!!
!!\verbatim
!!vertex   1  xyz  .0      .0       0.
!!vertex   2  xyz  1.      .0       0.
!!vertex   3  xyz  1.      .0       1.
!!vertex   4  xyz  .0      .0       1.
!!
!!vertex  11  fixed vertex 1
!!vertex  12  xyz  .0      1.       0.
!!vertex  13  xyz  .0      1.       1.
!!vertex  14  fixed vertex 4
!!\endverbatim
!!
!!With vertices defined the straight boundary curves are now defined. These curves are simply
!!constructed by specifying the start and end vertices from the list above.
!!
!!\verbatim
!!curve   1  vertex   1   2 output yes size def
!!curve   2  vertex   2   3 output yes size def
!!curve   3  vertex   3   4 output yes size def
!!curve   4  vertex   4   1 output yes size def
!!
!!curve  11  vertex   11  12 output yes size def
!!curve  12  vertex   12  13 output yes size def
!!curve  13  vertex   13  14 output yes size def
!!\endverbatim
!!
!!The curved edges of the cylinder require a more complex definition. In this case we are using
!!second order RBS to define the curves. These splines require a weight point in addition to the
!!two end points. Notice that curves 21 and 24 are the degenerate curves corresponding the bottom
!!and top surfaces respectively.
!!
!!\verbatim
!!curve  21  order  3  vertex   1  11 output yes size def
!!polygon 1  xyz   .0      .0      0. weight 0.70710678118655
!!curve  22  order  3  vertex   2  12 output yes size def
!!polygon 1  xyz   1.      1.      0. weight 0.70710678118655
!!curve  23  order  3  vertex   3  13 output yes size def
!!polygon 1  xyz   1.      1.      1. weight 0.70710678118655
!!curve  24  order  3  vertex   4  14 output yes size def
!!polygon 1  xyz   .0      .0      1. weight 0.70710678118655
!!\endverbatim
!!
!!With the curves defined the boundary surfaces can now be defined. The outer surface, which is
!!the only curved surface, does not require a separate weight point and is simply a defined by
!!its bounding curves. The remaining surfaces are all flat and are defined as "patches" that
!!are defined by the bounding curves and a normal direction.
!!
!!\verbatim
!!surface  1 curve  2 23 12 22 output yes size def
!!
!!patch 1 normal  0 -1 0 boundary curve 1 2 3 4 output yes size def
!!patch 2 normal -1 0 0 boundary curve -11 -12 -13 -4 output yes size def
!!patch 3 normal 0 0 -1 boundary curve -1 -22 11 21 output yes size def
!!patch 4 normal 0 0 1 boundary curve -3 23 13 -24 output yes size def
!!\endverbatim
!!
!!The volume to be meshed is now defined by specifying the bounding surfaces and patches.
!!
!!\verbatim
!!region 1 boundary surface 1 boundary patch 1 2 3 4 size def
!!\endverbatim
!!
!!\subsection doc_ex3_t3d_usage Generation
!!
!!With the geometry file defined the mesh can be created using the T3D program. For this example an
!!element size of 0.2 was used. Shown below are the command used to generate the mesh and the resulting
!!mesh before reflection. The file `cyl.inp` contains the geometry definitions described above while
!!the file `cyl.t3d` contains the generated mesh.
!!
!!\verbatim
!!~$ t3d -i cyl.inp -o cyl.t3d -d .2
!!\endverbatim
!!
!!\image html ex3_result.png
!!
!!\section doc_ex3_cubit Mesh Creation with CUBIT
!!
!!If the T3D mesh generation code is not avilable the mesh can be created using the CUBIT code using
!!the script below.
!!
!!\verbatim
!!reset
!!create Cylinder height 1 radius 1
!!volume 1 scheme Tetmesh
!!set tetmesher interior points on
!!set tetmesher optimize level 3 optimize overconstrained  off sliver  off
!!set tetmesher boundary recovery  off
!!volume 1 size .2
!!mesh volume 1
!!refine parallel fileroot 'cyl' overwrite no_execute
!!\endverbatim
!!
!!\note As this method does not use reflection the resulting fields will be slightly different from the T3D grid.
!!However, these differences should be small and only noticable when comapring convergence and eigenvalues.
!!
!!\section doc_ex3_gmsh Mesh Creation with Gmsh
!!
!!If the T3D and CUBIT mesh generation codes are not avilable the mesh can be created using the Gmsh code and the
!!geometry script below.
!!
!!\verbatim
!!Coherence;
!!Point(1) = {0, 0, 0, 1.0};
!!Point(2) = {1, 0, 0, 1.0};
!!Point(3) = {0, 1, 0, 1.0};
!!Point(4) = {-1, 0, 0, 1.0};
!!Point(5) = {0, -1, 0, 1.0};
!!Circle(1) = {2, 1, 3};
!!Circle(2) = {3, 1, 4};
!!Circle(3) = {4, 1, 5};
!!Circle(4) = {5, 1, 2};
!!Line Loop(5) = {2, 3, 4, 1};
!!Plane Surface(6) = {5};
!!Extrude {0, 0, 1} {
!!  Surface{6};
!!}
!!\endverbatim
!!
!!To generate a mesh, with resolution matching the T3D example above, place the script contents in a file called
!!`cyl.geo` and run the following command.
!!
!!\verbatim
!!~$ gmsh -3 -format mesh -optimize -clscale .2 -order 2 -o cyl.mesh cyl.geo
!!\endverbatim
!!
!!\note As this method does not use reflection the resulting fields will be slightly different from the T3D grid.
!!However, these differences should be small and only noticable when comapring convergence and eigenvalues.
!!
!!\section doc_ex3_code Code Walk Through
!!
!!\subsection doc_ex3_code_inc Module Includes
!!
!!The \ref taylor "Taylor" module requires the \ref lag_group "Lagrange" and \ref hcurl_group "H(Curl)"
!!finite element representations.
! START SOURCE
PROGRAM example3
!---Runtime
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
!---Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear Algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_lop_eigs, lag_setup_interp, lag_mloptions, &
  oft_lag_vgetmop, oft_lag_vproject
!---H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup
USE oft_hcurl_operators, ONLY: oft_hcurl_cinterp, hcurl_setup_interp, &
  hcurl_mloptions
!---Taylor state
USE taylor, ONLY: taylor_minlev, taylor_hmodes, taylor_hffa, ML_oft_hcurl, &
  ML_oft_lagrange, ML_oft_vlagrange
IMPLICIT NONE
#include "local.h"
!!\subsection ex3_code_vars Variable Definitions
!!
!!Although the core solvers are contained within the \ref taylor::taylor_hmodes "taylor_hmodes"
!!subroutine an additional solver is required for the field projection. A simple diagonally scaled
!!CG solver is used for this purpose as the mass matrix is not especially stiff. We also declare
!!a field interpolation object to evalute the \ref oft_hcurl_operators::oft_hcurl_cinterp "curl"
!!of the \ref taylor::taylor_hffa "vector potential" computed by \ref taylor::taylor_hmodes
!!"taylor_hmodes".
!---Lagrange mass solver
CLASS(oft_matrix), POINTER :: lmop => NULL()
CLASS(oft_solver), POINTER :: lminv => NULL()
!---Local variables
INTEGER(i4) :: i,ierr
INTEGER(i4), PARAMETER :: order = 3
REAL(r8), POINTER, DIMENSION(:) :: vals => NULL()
REAL(r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: bvout
CLASS(oft_vector), POINTER :: u,v,check
TYPE(oft_hcurl_cinterp) :: Bfield
TYPE(xdmf_plot_file) :: plot_file
TYPE(multigrid_mesh) :: mg_mesh
!!\subsection doc_ex3_code_grid Setup Grid
!!
!!As in the previous \ref ex1 "examples" the runtime environment, grid and plotting files must be setup
!!before the FE setup begins.
!---Initialize enviroment
CALL oft_init
!---Setup grid
CALL multigrid_construct(mg_mesh)
CALL plot_file%setup("Example3")
CALL mg_mesh%mesh%setup_io(plot_file,order)
!!\subsection doc_ex3_code_fem Setup FE Types
!!
!!As in \ref ex2 "example 2" we construct the finite element space, MG vector cache, and interpolation
!!operators. In this case the setup procedure is done for each required finite element space.
!---Lagrange
CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,ML_vlag_obj=ML_oft_vlagrange)
CALL lag_setup_interp(ML_oft_lagrange)
CALL lag_mloptions
!---H(Curl) space
CALL oft_hcurl_setup(mg_mesh,order,ML_oft_hcurl)
CALL hcurl_setup_interp(ML_oft_hcurl)
CALL hcurl_mloptions(ML_oft_hcurl)
!!\subsection doc_ex3_code_taylor Compute Taylor state
!!
!!The eigenstate is now computed using the \ref taylor::taylor_hmodes "taylor_hmodes" subroutine. The
!!number of eigenstates computed is controlled by an optional parameter to this subroutine. In this
!!case we are only computing the lowest eigenvalue (Spheromak state) so the default number of modes (1)
!!is used. The \ref taylor::taylor_hmodes "taylor_hmodes" subroutine uses MG preconditioning as
!!constructed by the \ref oft_hcurl_operators::hcurl_getwop_pre "hcurl_getwop_pre" subroutine. The
!!number of levels used is set by specifying the minimum level with the variable
!!\ref taylor::taylor_minlev "taylor_minlev".
!!
!!\note For our parallel example we increase the minimum level by 1 to use distributed levels only.
taylor_minlev=1
IF(oft_env%nprocs>1)taylor_minlev=2
oft_env%pm=.TRUE.
CALL taylor_hmodes(1)
!!\subsection doc_ex3_code_project Project Solution for Plotting
!!
!!In order to output the solution for plotting the field must be converted to a nodal, Lagrange,
!!representation. This is done by projecting the solution on to a Lagrange basis by integrating
!!the field against Lagrange test functions. The result is the then multiplied by the inverse of
!!the Lagrange metric matrix to produce the projected field. This field can then be exported for
!!plotting in VisIt.
!!
!!Integration is performed by the \ref oft_lag_operators::oft_lag_vproject "oft_lag_vproject"
!!subroutine, which takes a general \ref fem_utils::oft_fem_interp "interpolation" object that is
!!used to evaluate the field. The result of \ref taylor::taylor_hmodes "taylor_hmodes" is the vector
!!potential \ref taylor::taylor_hffa "taylor_hffa" in H(Curl) form so the \ref
!!oft_hcurl_operators::oft_hcurl_cinterp "oft_hcurl_cinterp" object is used to provide evaluation of
!!\f$ B = \nabla \times A \f$.
!---Construct operator
NULLIFY(lmop)
CALL oft_lag_vgetmop(ML_oft_vlagrange%current_level,lmop,'none')
!---Setup solver
CALL create_cg_solver(lminv)
lminv%A=>lmop
lminv%its=-2
CALL create_diag_pre(lminv%pre) ! Setup Preconditioner
!---Create solver fields
CALL ML_oft_vlagrange%vec_create(u)
CALL ML_oft_vlagrange%vec_create(v)
!---Setup field interpolation
Bfield%u=>taylor_hffa(1,ML_oft_hcurl%level)%f
CALL Bfield%setup(ML_oft_hcurl%current_level)
!---Project field
CALL oft_lag_vproject(ML_oft_lagrange%current_level,Bfield,v)
CALL u%set(0.d0)
CALL lminv%apply(u,v)
!---Retrieve local values and save
ALLOCATE(bvout(3,u%n/3))
vals=>bvout(1,:)
CALL u%get_local(vals,1)
vals=>bvout(2,:)
CALL u%get_local(vals,2)
vals=>bvout(3,:)
CALL u%get_local(vals,3)
call mg_mesh%mesh%save_vertex_vector(bvout,plot_file,'B')
!---Finalize enviroment
CALL oft_finalize
END PROGRAM example3
! STOP SOURCE
!!
!!\section doc_ex3_input Input file
!!
!!Below is an input file which can be used with this example in a serial environment. For this example
!!we also introduce the use of curved tetrahedrons on the boundary. High order tetrahedra can be
!!constructed to provide a more accurate boundary using the runtime option `grid_order` in the
!!`mesh_options` group. This option causes OFT to generate a quadratic representation for all boundary
!!elements using the available CAD information.
!!
!!\verbatim
!!&runtime_options
!! ppn=1
!! debug=0
!!/
!!
!!&mesh_options
!! meshname='cylinder'
!! cad_type=1
!! nlevels=1
!! nbase=1
!! grid_order=2
!!/
!!
!!&t3d_options
!! filename='cyl.t3d'
!! inpname='cyl.inp'
!! reflect='xy'
!! ref_per=false,false,false
!!/
!!
!!&lag_op_options
!! df_lop=1.,.830,.619
!! nu_lop=64,2,1
!!/
!!
!!&hcurl_op_options
!! df_wop=.689,.390,.316
!! nu_wop=64,2,1
!!/
!!\endverbatim
!!
!!\subsection doc_ex3_input_cubit CUBIT additions
!!
!!If CUBIT was used to generate the mesh changes must be made to the `mesh_options` group and
!!the `cubit_options` group must be added as below. These changes are the same for the MPI input
!!file below as well.
!!
!!@note Only the fields whose values are changed are shown below.
!!
!!\verbatim
!!&mesh_options
!! cad_type=2
!! ...
!!/
!!
!!&cubit_options
!! filename='cyl.in.e'
!! inpname='cyl.3dm'
!! lf_file=T
!!/
!!\endverbatim
!!
!!\subsection doc_ex3_input_gmsh GMSH additions
!!
!!If GMSH was used to generate the mesh changes must be made to the `mesh_options` group and
!!the `gmsh_options` group must be added as below. These changes are the same for the MPI input
!!file below as well.
!!
!!@note Only the fields whose values are changed are shown below.
!!
!!\verbatim
!!&mesh_options
!! cad_type=3
!! ...
!!/
!!
!!&gmsh_options
!! filename='cyl.mesh'
!! order=2
!!/
!!\endverbatim
!!
!!\subsection doc_ex3_input_mpi Parallel input file
!!
!!The input file below will provide the same preconditioner as the serial example, but can
!!be run in parallel.
!!
!!\verbatim
!!&runtime_options
!! ppn=1
!! debug=0
!!/
!!
!!&mesh_options
!! meshname='cylinder'
!! cad_type=1
!! nlevels=2
!! nbase=1
!! grid_order=2
!!/
!!
!!&t3d_options
!! filename='cyl.t3d'
!! inpname='cyl.inp'
!! reflect='xy'
!! ref_per=false,false,false
!!/
!!
!!&lag_op_options
!! df_lop=0.,1.,.830,.619
!! nu_lop=0,64,2,1
!!/
!!
!!&hcurl_op_options
!! df_wop=0.,.689,.390,.316
!! nu_wop=0,64,2,1
!!/
!!\endverbatim
!!
