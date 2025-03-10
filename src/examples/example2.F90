!!API Example 2: Multi-Grid Preconditioners    {#doc_api_ex2}
!!=====================================
!!
!![TOC]
!!
!!This program expands on \ref doc_ex1 by solving the Poisson's equation, but replaces diagonal scaling
!!with a more efficient Multi-Grid preconditioner to accelerate the solver. It also introduces how to
!!handle cases where Native ane PETSc implementations will differ. Although, PETSc matrices will work
!!with almost all Native solvers, they do not currently work with the Jacobi smoothers used here. Additionally,
!!it is desirable to use the PETSc solver implementations for any complicated nested operation, such as
!!Multi-Grid preconditioning. In order to make the code portable preprocessor gaurds are added to select
!!the appropriate implementation depending on the existence of PETSc during build.
!!
!!\section doc_ex2_code Code Walk Through
!!
!!There are many similarities between this example and the previous example. As a result only the differences
!!and elements relating directly to MG will be discussed here. Please refer to \ref ex1 for the basic program
!!outline and I/O.
!!
!!\subsection doc_ex2_code_inc Module Includes
!!
!!In order to construct MG preconditioners we must include the MG constructor \ref oft_la_utils::create_mlpre
!!"create_mlpre".
! START SOURCE
PROGRAM example2
!---Runtime
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
!---Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear Algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_mlpre
!---Lagrange FE space
USE fem_base, ONLY: oft_ml_fem_type, oft_ml_fe_vecspace
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, lag_mloptions, &
  oft_lag_getmop, oft_lag_getlop, df_lop, nu_lop, oft_lag_zerob, &
  lag_base_pop, lag_base_push
IMPLICIT NONE
!!\subsection ex2_code_vars Variable Definitions
!---Solver objects
CLASS(oft_solver), POINTER :: linv
!--- ML structures for MG-preconditioner
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_int => NULL()
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_lop => NULL()
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: levels
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: nu
REAL(r8), ALLOCATABLE, DIMENSION(:) :: df
!---Local variables
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: lop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
REAL(r8) :: uu
REAL(r8), POINTER, DIMENSION(:) :: vtmp => NULL()
INTEGER(i4) :: i,nlevels
INTEGER(i4), PARAMETER :: order = 3
TYPE(xdmf_plot_file) :: plot_file
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange,ML_oft_blagrange
TYPE(oft_ml_fem_comp_type), TARGET :: ML_oft_vlagrange
TYPE(oft_lag_zerob) :: lag_zerob
TYPE(oft_ml_fe_vecspace) :: ml_vecspace
!!\subsection doc_ex2_code_fem Setup Lagrange FE
!!
!!When MG is used it is also desirable to construct field caches for each FE representation to be used. This is done with
!!the `*_vcache` subroutines in each of the corresponding FE `*_fields` modules. These subroutines create a field cache which
!!prevents the recreation of internal field lists, which are otherwise recreated each time a new vector is created using
!!`*_create` routines. See \ref mat_vec_vconst for more information.
!---Initialize enviroment
CALL oft_init
!---Setup grid
CALL multigrid_construct(mg_mesh)
CALL plot_file%setup("Example2")
CALL mg_mesh%mesh%setup_io(plot_file,order)
!---Construct FE levels
CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,ML_oft_blagrange,ML_oft_vlagrange)
CALL lag_setup_interp(ML_oft_lagrange)
CALL lag_mloptions
lag_zerob%ML_lag_rep=>ML_oft_lagrange
!!\subsection doc_ex2_code_ml Construct ML structures
!!
!!To form a MG preconditioner the approximate matrices corresponding to each level must be
!!constructed and passed to form the preconditioner. The transfer level is skipped by setting
!!the level to the negative of the true FE level. This causes the preconditioner to skip smoothing
!!and only interpolate/inject on this level.
nlevels=ML_oft_lagrange%nlevels
ALLOCATE(ml_lop(nlevels),ml_int(nlevels-1))
ALLOCATE(df(nlevels),nu(nlevels),levels(nlevels))
DO i=1,nlevels
  CALL ML_oft_lagrange%set_level(i)
  levels(i)=i
  IF(i==ML_oft_lagrange%blevel+1)levels(i)=-levels(i)
  df(i)=df_lop(i)
  nu(i)=nu_lop(i)
  !---
  NULLIFY(ml_lop(i)%M)
  CALL oft_lag_getlop(ML_oft_lagrange%current_level,ml_lop(i)%M,'zerob')
  IF(i>1)ml_int(i-1)%M=>ML_oft_lagrange%interp_matrices(i)%m !oft_lagrange_ops%interp
END DO
CALL ML_oft_lagrange%set_level(ML_oft_lagrange%nlevels)
!!\subsection doc_ex2_code_fields Setup solver fields
!---Create solver fields
CALL ML_oft_lagrange%vec_create(u)
CALL ML_oft_lagrange%vec_create(v)
!---Get FE operators
lop=>ml_lop(nlevels)%M
CALL oft_lag_getmop(ML_oft_lagrange%current_level,mop,'none')
!!\subsection doc_ex2_code_solver Setup linear solver
CALL create_cg_solver(linv)
linv%its=-3
linv%A=>lop
!---Setup MG preconditioner
ml_vecspace%ML_FE_rep=>ML_oft_lagrange
ml_vecspace%base_pop=>lag_base_pop
ml_vecspace%base_push=>lag_base_push
CALL create_mlpre(linv%pre,ml_lop,levels,nlevels=nlevels,ml_vecspace=ml_vecspace, &
       bc=lag_zerob,stype=1,df=df,nu=nu)
!!\subsection doc_ex2_code_compute Solve linear system
CALL u%set(1.d0)
CALL mop%apply(u,v)
CALL lag_zerob%apply(v)
CALL u%set(0.d0)
CALL linv%apply(u,v)
CALL u%get_local(vtmp)
CALL mg_mesh%mesh%save_vertex_scalar(vtmp,plot_file,'T')
DEALLOCATE(vtmp)
!---Finalize enviroment
CALL oft_finalize
END PROGRAM example2
! STOP SOURCE
!!
!!\section doc_ex2_input Input file
!!
!!Below is an input file which can be used with this example in a serial environment. The
!!`lag_op_options` group has been added over the groups in \ref ex1 to set parameters used
!!by the multigrid preconditioner. In this case we are setting the smoother coefficient `df`
!!and number of iterations `nu` used by the Jacobi smoother on each level.
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
!!
!!&lag_op_options
!! df_lop=1.,1.,1.,1.,.826,.645
!! nu_lop=10,64,8,4,2,1
!!/
!!\endverbatim
!!
!!\subsection doc_ex2_input_mpi Parallel input file
!!
!!The number of grid levels changes when you got to a parallel environment, as there
!!is now a transfer level in addition to the refinement levels, see \ref doc_ex1_input_mpi "Example 1".
!!As a result you must modify the preconditioner options to take this into account.
!!The input file below will provide the same preconditioner as the serial example, but can
!!be run in parallel.
!!
!!@note The smoother coefficient and number of iterations are set to zero for the transfer
!!level. This is only for clarity of the input file however, as smoothing is skipped altogether on the
!!transfer level.
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
!!
!!&lag_op_options
!! df_lop=1.,1.,1.,0.,1.,.826,.645
!! nu_lop=10,64,8,0,4,2,1
!!/
!!\endverbatim
!!
