!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!!Marklin Example: Force-free eigenmodes    {#doc_marklin_ex1}
!!========================
!!
!![TOC]
!!
!!This example introduces the \ref taylor "Taylor State" physics module. The \ref taylor "Taylor"
!!module is used to compute Ideal MHD force-free states with uniform \f$ \lambda \f$. This example
!!computes the Taylor state in a oblate cylinder. The T3D interface is used to generate a mesh for
!!the cylinder and provide boundary information for the setup of a quadratic spatial mapping. This
!!example also demonstrates the steps required to \ref doc_marklin_ex1_code_project "output a vector field"
!!for plotting.
!!
!!\section doc_marklin_ex1_code Code Walk Through
!!
!!\subsection doc_marklin_ex1_code_inc Module Includes
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
!
USE fem_composite, ONLY: oft_ml_fem_comp_type
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_lop_eigs, lag_setup_interp, lag_mloptions, &
  oft_lag_vgetmop, oft_lag_vproject
!---H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup
USE oft_hcurl_operators, ONLY: oft_hcurl_cinterp, hcurl_setup_interp, &
  hcurl_mloptions
!---Taylor state
USE taylor, ONLY: taylor_hmodes, oft_taylor_hmodes
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
TYPE(oft_ml_fem_comp_type) :: ML_vlagrange
TYPE(oft_taylor_hmodes) :: taylor_states
!!\subsection doc_marklin_ex1_code_grid Setup Grid
!!
!!As in the previous \ref ex1 "examples" the runtime environment, grid and plotting files must be setup
!!before the FE setup begins.
!---Initialize enviroment
CALL oft_init
!---Setup grid
CALL multigrid_construct(mg_mesh)
CALL plot_file%setup("Example3")
CALL mg_mesh%mesh%setup_io(plot_file,order)
!!\subsection doc_marklin_ex1_code_fem Setup FE Types
!!
!!As in \ref ex2 "example 2" we construct the finite element space, MG vector cache, and interpolation
!!operators. In this case the setup procedure is done for each required finite element space.
!---Lagrange
ALLOCATE(taylor_states%ML_lagrange)
CALL oft_lag_setup(mg_mesh,order,taylor_states%ML_lagrange,ML_vlag_obj=ML_vlagrange)
CALL lag_setup_interp(taylor_states%ML_lagrange)
CALL lag_mloptions
!---H(Curl) space
ALLOCATE(taylor_states%ML_hcurl)
CALL oft_hcurl_setup(mg_mesh,order,taylor_states%ML_hcurl)
CALL hcurl_setup_interp(taylor_states%ML_hcurl)
CALL hcurl_mloptions(taylor_states%ML_hcurl)
!!\subsection doc_marklin_ex1_code_taylor Compute Taylor state
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
taylor_states%minlev=1
IF(oft_env%nprocs>1)taylor_states%minlev=2
oft_env%pm=.TRUE.
CALL taylor_hmodes(taylor_states,1)
!!\subsection doc_marklin_ex1_code_project Project Solution for Plotting
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
CALL oft_lag_vgetmop(ML_vlagrange%current_level,lmop,'none')
!---Setup solver
CALL create_cg_solver(lminv)
lminv%A=>lmop
lminv%its=-2
CALL create_diag_pre(lminv%pre) ! Setup Preconditioner
!---Create solver fields
CALL ML_vlagrange%vec_create(u)
CALL ML_vlagrange%vec_create(v)
!---Setup field interpolation
Bfield%u=>taylor_states%hffa(1,taylor_states%ML_hcurl%level)%f
CALL Bfield%setup(taylor_states%ML_hcurl%current_level)
!---Project field
CALL oft_lag_vproject(taylor_states%ML_lagrange%current_level,Bfield,v)
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
!!\section doc_marklin_ex1_input Input file
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
!! cad_type=0
!! nlevels=1
!! nbase=1
!! grid_order=2
!!/
!!
!!&native_mesh_options
!! filename='cyl.h5'
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
!!\subsection doc_marklin_ex1_input_mpi Parallel input file
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
!! cad_type=0
!! nlevels=2
!! nbase=1
!! grid_order=2
!!/
!!
!!&native_mesh_options
!! filename='cyl.h5'
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
!!\section doc_marklin_ex1_mesh Mesh Creation
!! A mesh file `cyl.h5` is provided with this example. Instructions to generate your
!! own mesh for the geometry using [CUBIT](https://cubit.sandia.gov/) and [GMSH](https://gmsh.info/).
!!
!!\subsection doc_marklin_ex1_cubit Meshing with CUBIT
!!
!! A suitable mesh for this example, with radius of 1m and height of 2m, can be created using
!! the CUBIT script below.
!!
!!\verbatim
!!reset
!!
!!create Cylinder height 1 radius 1
!!
!!volume 1 scheme Tetmesh
!!set tetmesher interior points on
!!set tetmesher optimize level 3 optimize overconstrained  off sliver  off
!!set tetmesher boundary recovery  off
!!volume 1 size .2
!!mesh volume 1
!!
!!set duplicate block elements off
!!block 1 add volume 1 
!!block 1 element type tetra10
!!
!!set large exodus file on
!!export Genesis  "cyl.g" overwrite block 1
!!\endverbatim
!!
!! Once complete the mesh should be converted into the native mesh format using the `convert_cubit.py` script as
!! below. The script is located in `bin` following installation or `src/utilities` in the base repo.
!!
!!\verbatim
!!~$ python convert_cubit.py --in_file=cyl.g
!!\endverbatim
!!
!!\subsection doc_marklin_ex1_gmsh Meshing with Gmsh
!!
!! If the CUBIT mesh generation codes is not avilable the mesh can be created using the Gmsh code and the
!! geometry script below.
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
!! To generate a mesh, with resolution matching the Cubit example above, place the script contents in a file called
!! `cyl.geo` and run the following command.
!!
!!\verbatim
!!~$ gmsh -3 -format mesh -optimize -clscale .2 -order 2 -o cyl.mesh cyl.geo
!!\endverbatim
!!
!! Once complete the mesh should be converted into the native mesh format using the `convert_gmsh.py` script as
!! below. The script is located in `bin` following installation or `src/utilities` in the base repo.
!!
!!\verbatim
!!~$ python convert_gmsh.py --in_file=cyl.mesh
!!\endverbatim
