!!API Example 1: Helmholtz solve    {#doc_api_ex1}
!!========================
!!
!![TOC]
!!
!!This example outlines the basic elements of a program which uses the Open FUSION Toolkit (OFT) to solve a PDE.
!!In this case computing the lowest eigenvalue of the scalar Helmholtz equation
!!\f$ ( \nabla^2 - \lambda ) u = 0 \f$.
!!
!!\section doc_ex1_code Code Walk Through
!!
!!The code consists of three basic sections, required imports and variable definitions,
!!finite element setup, and system creation and solution.
!!
!!\subsection doc_ex1_code_inc Module Includes
!!The first thing that you must do is include the OFT modules which contain
!!the required functions and variables. It is good practice to restrict the included elements
!!to only those needed. This is done using the `ONLY` clause to specifically include only
!!certain definitions. The exceptions to this practice are the \ref local and \ref oft_base
!!modules, which contain a controlled set of commonly used elements that can be safely imported
!!as a whole.
! START SOURCE
PROGRAM example1
!---Runtime
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
!---Grid
USE oft_mesh_type, ONLY: mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_native_solvers, ONLY: oft_native_cg_eigsolver
USE oft_solver_utils, ONLY: create_diag_pre
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_fields, ONLY: oft_lag_create
USE oft_lag_operators, ONLY: oft_lag_getlop, oft_lag_getmop, lag_zerob
IMPLICIT NONE
!!The first two modules import runtime and helper functions. The \ref tetmesh_local module contains
!!the main global mesh object \ref tetmesh_local::mesh "mesh" and \ref multigrid_build contains
!!the grid construction driver. The three `oft_lag_*` modules contain the finite element basis
!!and setup routines, vector definitions, and operator routines respectively for a Lagrange finite
!!element representation. \ref oft_cg contains the native Conjugate-Gradient solvers and \ref oft_io
!!contains output routines for producing plot files.
!!
!!\subsection doc_ex1_code_vars Variable Definition
!!
!!This defines the vector object used to compute the eigenvalue/vector pair as well as the matrix
!!objects for the Laplace and Mass operators. An eigensolver is also instantiated along with a
!!preconditioner which utilizes diagonal scaling.
!!
!!\note Preconditioners are not generally interoperable between different backends, ie. native vs PETSc.
!!As a result different preconditioner objects should be used for each backend and controlled using
!!preprocessor gaurds. For more information on available preconditioners see \ref doc_solvers
CLASS(oft_vector), POINTER :: u
CLASS(oft_matrix), POINTER :: lop,mop
TYPE(oft_native_cg_eigsolver) :: solver
INTEGER(i4), PARAMETER :: order = 3
REAL(r8) :: lambda
REAL(r8), POINTER, DIMENSION(:) :: vtmp => NULL()
TYPE(xdmf_plot_file) :: plot_file
!!\subsection doc_ex1_code_init Initialize Enviroment
!!
!!This call setups of the basics OFT run environment, including initializing MPI and PETSc if
!!applicable.
CALL oft_init
!!\subsection doc_ex1_code_grid Setup Grid
!!
!!This call constructs the grid levels through heirarchical refinement,
!!\ref multigrid_build::multigrid_construct "multigrid_construct".
CALL multigrid_construct
!!\subsection doc_ex1_code_hdf5 Create Output Files
!!
!!This call sets up metadata files for I/O and saves the mesh for use with solution fields output
!!later in the program. The plotting grid handles high order fields by tesselating the mesh to
!!produce additional tets using the new node points. The level of tesselation is set by
!!the first argument to \ref oft_mesh_type::oft_mesh::setup_io "mesh%setup_io".
CALL plot_file%setup("Example1")
CALL mesh%setup_io(plot_file,order)
!!\subsection doc_ex1_code_fem Setup Lagrange FE
!!
!!\ref oft_lag_basis::oft_lag_setup "oft_lag_setup" builds the Lagrange finte elements on each grid
!!and polynomial level. This create element interaction lists as well as boundary and seam information.
!!All FE index fields are encapsulated in the \ref fem_base::oft_fem_type "oft_fem_type" structure,
!!see \ref fem_base::fem_setup "fem_setup".
CALL oft_lag_setup(order)
!!\subsection doc_ex1_code_ops Setup linear system
!!
!!Solving the Helmholtz eigensystem requires the operators coresponding the general eigenvalue problem
!!\f$ Ax = \lambda Mx \f$. These operators are the distcretized Laplace operator
!!\f$ A = \int \nabla u_i^T \cdot \nabla u_j dV \f$ and the finite element mass matrix
!!\f$ M = \int u_i^T u_j dV \f$. Where \c u and \f$ u^T \f$ are the Lagrange basis and test functions
!!respectively. A vector defining the solution space is also required, to store the solution and to
!!create worker vectors. In this case a dirichlet boundary condition is used where \c u is 0 everywhere
!!on the boundary, denoted by the BC flag `'zerob'` used in matrix construction.
CALL oft_lag_create(u)
!---Create Operators
NULLIFY(lop,mop)
CALL oft_lag_getlop(lop,'zerob')
CALL oft_lag_getmop(mop,'zerob')
!!\subsection doc_ex1_code_solver Setup solver
!!
!!This section assembles the solver object by fill the required references. The solver object used here
!!is the native non-linear Conjugate-Gradient iterative method. It requires the right and left hand side
!!matrices for the generalized EV problem, \c A and \c M. The convergence tolerance is also specified
!!to be double precision convergence of the eigenvalue, see \ref oft_cg::oft_cg_eigsolver
!!"oft_cg_eigsolver". Preconditioning is supplied by diagonal scaling, which requires identifying the
!!preconditioner matrix, in this case \c A.
!!
!!The vector is then initialized with a guess solution and the boundary condition is applied to make sure
!!the initial guess is consistent with the boundary condition. Finally, the assembled solver is used
!!to compute the lowest eigenvalue/eigenvector pair.
solver%A=>lop
solver%M=>mop
solver%its=-2
CALL create_diag_pre(solver%pre) ! Setup Preconditioner
!---Compute EV
CALL u%set(1.d0)
CALL lag_zerob(u)
CALL solver%apply(u,lambda)
!!\subsection doc_ex1_code_plot Save Solution
!!
!!The solution is output to the HDF5 dump files for plotting at a later time in VisIt. In order to save the
!!field the local representation must be retrieved from the vector object using
!!\ref oft_vectors::oft_vector::get_local "get_local". These values may then be saved using the
!!\ref oft_io::hdf5_spdata "hdf5_spdata" subroutine. When the program has completed
!!\ref oft_base::oft_finalize "oft_finalize" is called to cleanup the runtime environment and terminate
!!the process. This subroutine calls any required MPI and PETSc finalize subroutines as well. After
!!completing the run, the \c build_xdmf script may be used to construct VisIt files and view the solution
!!field, saved as tag \c T.
CALL u%get_local(vtmp)
CALL mesh%save_vertex_scalar(vtmp,plot_file,'T')
DEALLOCATE(vtmp)
!---Program Stop
CALL oft_finalize
END PROGRAM example1
! STOP SOURCE
!!
!!\section doc_ex1_input Input file
!!
!!Below is an input file which can be used with this example in a serial environment. For this example
!!we have disabled additional debug output, `debug=0`, and are using the built in cube mesh type,
!!`cad_type=92`, with 4 grid refinements, `nbase=nlevels=4`.
!!
!!\verbatim
!!&runtime_options
!! ppn=1
!! debug=0
!!/
!!
!!&mesh_options
!! meshname='cube'
!! cad_type=92
!! nlevels=4
!! nbase=4
!!/
!!\endverbatim
!!
!!\subsection doc_ex1_input_mpi Parallel input file
!!
!!In order to run the code in a distributed memory environment, with MPI,
!!the mesh must be decomposed. This requires that `nlevels>nbase` as explained in
!!\ref doc_input_mesh. To perform a parallel run that is analogous to the serial run shown
!!above the input file must be modified as below. The choice of nbase is somewhat arbitrary
!!here, but generally is a function of the number of MPI tasks and the size of the base mesh.
!!\c nlevels however, must be incremented by 1 in order to account for the additional transfer
!!level create during decomposition.
!!
!!\verbatim
!!&runtime_options
!! ppn=1
!! debug=0
!!/
!!
!!&mesh_options
!! meshname='cube'
!! cad_type=92
!! nlevels=5
!! nbase=3
!!/
!!\endverbatim
!!
