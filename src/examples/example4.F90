!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!!Marklin Example: Inhomogeneous Ideal MHD Equilibria    {#doc_marklin_ex2}
!!=============================================
!!
!![TOC]
!!
!! This example shows how to use the \ref taylor "Taylor State" physics module
!! and mesh cut planes to compute Inhomogeneous Ideal MHD Equilibria with uniform
!! \f$ \lambda \f$. The code computes a composite Taylor state in HIT-SI with a flux
!! ratio of 6.
!!
!!\section doc_marklin_ex2_code Code Walk Through
!!
!!\subsection doc_marklin_ex2_code_inc Module Includes
!!
!!
! START SOURCE
PROGRAM example4
!---Runtime
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
!---Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear Algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre\
!
USE fem_composite, ONLY: oft_ml_fem_comp_type
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_lop_eigs, lag_setup_interp, lag_mloptions, &
  oft_lag_vgetmop, oft_lag_vproject
!---H1 FE space (Grad(H^1) subspace)
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_mloptions, h1_setup_interp
!---Full H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: oft_hcurl_cinterp, hcurl_setup_interp, &
  hcurl_mloptions
!---Taylor state
USE taylor, ONLY: taylor_hmodes, oft_taylor_rinterp, taylor_vacuum, &
  taylor_injectors, oft_taylor_hmodes, oft_taylor_ifield, taylor_tag_size
!---Tracing
USE tracing, ONLY: oft_tracer, create_tracer, tracing_poincare
IMPLICIT NONE
#include "local.h"
!!\subsection ex4_code_vars Variable Definitions
!!
!!
!---Lagrange mass solver
CLASS(oft_matrix), POINTER :: lmop => NULL()
CLASS(oft_solver), POINTER :: lminv => NULL()
!---Fixed constants
INTEGER(i4), PARAMETER :: order = 2
INTEGER(i4), PARAMETER :: nh=2
REAL(r8), PARAMETER :: fr = 6
!---Local variables
INTEGER(i4) :: i,npts,ierr
REAL(r8) :: fluxes(nh),hcpc(3,nh),hcpv(3,nh)
CHARACTER(LEN=taylor_tag_size) :: htags(nh)
REAL(r8), POINTER, DIMENSION(:) :: vals => NULL()
REAL(r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: bvout,pts
CLASS(oft_vector), POINTER :: u,v,check
TYPE(oft_taylor_rinterp), TARGET :: Bfield
CLASS(oft_tracer), POINTER :: tracer
TYPE(xdmf_plot_file) :: plot_file
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_comp_type) :: ML_vlagrange
TYPE(oft_taylor_hmodes) :: hmodes
TYPE(oft_taylor_ifield) :: ff_obj
!!\subsection doc_marklin_ex2_code_grid Setup Grid
!!
!!As in the previous \ref ex1 "examples" the runtime environment, grid and plotting
!!files must be setup before the FE setup begins.
!---Initialize enviroment
CALL oft_init
!---Setup grid
CALL multigrid_construct(mg_mesh)
CALL plot_file%setup("Example4")
CALL mg_mesh%mesh%setup_io(plot_file,order)
!!\subsection doc_marklin_ex2_code_fem Setup FE Types
!!
!!As in \ref ex2 "example 2" we construct the finite element space, MG vector cache, and interpolation
!!operators. In this case the setup procedure is done for each required finite element space.
!--- Lagrange
ALLOCATE(hmodes%ML_lagrange)
ff_obj%ML_lagrange=>hmodes%ML_lagrange
CALL oft_lag_setup(mg_mesh,order,hmodes%ML_lagrange,ML_vlag_obj=ML_vlagrange)
CALL lag_setup_interp(hmodes%ML_lagrange)
CALL lag_mloptions
!--- Grad(H^1) subspace
ALLOCATE(ff_obj%ML_H1)
CALL oft_h1_setup(mg_mesh,order+1,ff_obj%ML_H1)
CALL h1_setup_interp(ff_obj%ML_H1)
CALL h1_mloptions
!--- H(Curl) subspace
ALLOCATE(hmodes%ML_hcurl)
ff_obj%ML_hcurl=>hmodes%ML_hcurl
CALL oft_hcurl_setup(mg_mesh,order,hmodes%ML_hcurl)
CALL hcurl_setup_interp(hmodes%ML_hcurl)
CALL hcurl_mloptions(hmodes%ML_hcurl)
!--- Full H(Curl) space
ALLOCATE(ff_obj%ML_hcurl_grad,ff_obj%ML_h1grad)
CALL oft_hcurl_grad_setup(hmodes%ML_hcurl,ff_obj%ML_H1,ff_obj%ML_hcurl_grad,ff_obj%ML_h1grad)
!!\subsection doc_marklin_ex2_code_taylor Compute Taylor state
!!
!!For composite Taylor states the lowest eigenmode is used used in addition to the injector fields. This
!!is possible in HIT-SI due to decoupling of the injector vacuum field and lowest eigenmode. The lowest
!!eigenmode is computed as in \ref ex3 "example 3" using \ref taylor::taylor_hmodes "taylor_hmodes".
hmodes%minlev=1
IF(oft_env%nprocs>1)hmodes%minlev=2
ff_obj%minlev=hmodes%minlev
CALL taylor_hmodes(hmodes,1)
!!\subsection doc_marklin_ex2_code_jumps Injector Cut Planes
!!
!!The ideal MHD equilibrium solvers with inhomogeneous source terms require vacuum fields to use as the
!!external source term. For HIT-SI the \ref taylor::taylor_vacuum "taylor_vacuum" can be used to compute
!!vacuum fields for the injectors from cut planes in the mesh. Cut planes are identified by specifying
!!the cut plane center (`hcpc`) and normal direction (`hcpv`) for each jump. The jump is also localized
!!to a circular region around the center point, defined by the magnitude of `hcpv`, such that \f$ \left( \vec{r}
!!- \vec{r}_{hcpc} \right) \cdot \vec{r}_{hcpv} \leq 1 \f$. Each jump can also be given a character identifier
!!which is used here to indicate the naming convention on HIT-SI.
!--- X-Injector
hcpc(:,1)=(/.0,.0,-.6/)
hcpv(:,1)=(/-5.,.0,.0/)
htags(1)='Xinj'
!--- Y-Injector
hcpc(:,2)=(/.0,.0,.6/)
hcpv(:,2)=(/.0,-5.,.0/)
htags(2)='Yinj'
!!\subsection doc_marklin_ex2_code_inj Compute Inhomogeneous state
!!
!!With cut planes specified, the injector equilibrium can then be compute using the \ref taylor::taylor_vacuum
!!"taylor_vacuum" and \ref taylor::taylor_injectors "taylor_injectors" subroutines. \ref taylor::taylor_injectors
!!"taylor_injectors" computes the plasma response using the vacuum fields computed for the cut planes with unit
!!fluxes for a specified \f$ \lambda \f$. Since we are interested in composite Taylor states \f$ \lambda \f$ must
!!be equal to the value for the lowest homogeneous eigenmode (Taylor state).
!!
!!Once the plasma response and vacuum fields have been computed the composite field can be projected for plotting. The
!!full field is defined as \f$\vec{B} = \vec{B}_{hffa} + flux^i \left( \vec{B}^i_{hvac} + \vec{B}^i_{gffa} \right)\f$,
!!where \f$ \vec{B}_{hffa} \f$ is the field for the lowest eigenmode, \f$ \vec{B}^i_{hvac} \f$ is the vacuum field for
!!the injector, and \f$ \vec{B}^i_{gffa} \f$ is the plasma component of the inhomogeneous equilibrium field. The
!!interpolation object \ref taylor::oft_taylor_rinterp "oft_taylor_rinterp" is designed to support this type of field
!!and is populated once the subfields are computed.
CALL ff_obj%setup(nh,hcpc,hcpv,htags)
CALL taylor_vacuum(ff_obj,hmodes=hmodes)
CALL taylor_injectors(ff_obj,hmodes,hmodes%hlam(1,hmodes%ML_hcurl%level))
!---Setup field interpolation object
fluxes=(/1.d0,0.d0/)
CALL ff_obj%ML_hcurl_grad%vec_create(Bfield%uvac)
DO i=1,nh
  CALL Bfield%uvac%add(1.d0,fluxes(i),ff_obj%hvac(i)%f)
END DO
CALL hmodes%ML_hcurl%vec_create(Bfield%ua)
CALL Bfield%ua%add(0.d0,fr/hmodes%htor(1,hmodes%ML_hcurl%level),hmodes%hffa(1,hmodes%ML_hcurl%level)%f)
DO i=1,nh
  CALL Bfield%ua%add(1.d0,fluxes(i),ff_obj%gffa(i)%f)
END DO
!!\subsection doc_marklin_ex2_code_project Project Solution for Plotting
!!
!!With the interpolation operator defined the field can be projected to a Lagrange basis for plotting as in
!!\ref ex3 "example 3".
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
CALL Bfield%setup(hmodes%ML_hcurl%current_level,ff_obj%ML_H1%current_level)
!---Project field
CALL oft_lag_vproject(hmodes%ML_lagrange%current_level,Bfield,v)
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
!!\subsection doc_marklin_ex2_code_poincare Create Poincare section
!!
!! Poincare sections can be created using the subroutine \ref tracing::tracing_poincare
!! "tracing_poincare". This subroutine requires a tracing object, which is used to advance
!! the streamline ODE and a list of launch points for each streamline. The tracing object
!! is used as a template to spawn additional tracers that are advanced in parallel. Tracing
!! objects can be created using the \ref tracing::create_tracer "create_tracer" subroutine.
!! Once the tracer has been created basic properties are set including the tracing tolerance,
!! maximum number of steps, maximum number of domain transfers and the field interpolation
!! object.
!!
!! For this example a line of points in the xz-plane are used, positioned at the radius of
!! the magnetic axis of the Taylor state. The resulting puncture list is stored in the file
!! 'Example4.poin'.
!!
!! \note Poincare tracing requires the Open FUSION Toolkit (OFT) to be built with OpenMP enabled. A single thread
!! may still be used for execution by setting the environment variable "OMP_NUM_THREADS=1",
!! and a second thread will be used only during tracing.
!---Setup tracer
CALL create_tracer(tracer,1)
tracer%tol=1.d-9
tracer%maxsteps=1E5
tracer%maxtrans=2E3
tracer%B=>Bfield
!---Create launch point list
npts=20
ALLOCATE(pts(3,npts))
DO i=1,npts
  pts(:,i)=(/0.d0,.33d0,-.3d0 + .6d0*i/REAL(npts,8)/)
END DO
!---Perform tracing
CALL tracing_poincare(tracer,pts,npts,'Example4.poin')
!---Finalize enviroment
CALL oft_finalize
END PROGRAM example4
! STOP SOURCE
!!
!!\section doc_marklin_ex2_input Input file
!!
!!Below is an input file which can be used with this example in a serial environment.
!!
!!\verbatim
!!&runtime_options
!! ppn=1
!! debug=0
!!/
!!
!!&mesh_options
!! meshname='hitsi'
!! cad_type=1
!! nlevels=1
!! nbase=1
!! grid_order=2
!!/
!!
!!&t3d_options
!! filename='hitsi_rs.t3d'
!! inpname='hitsi_rs.inp'
!! reflect='xy'
!!/
!!
!!&lag_op_options
!! df_lop=1.,.837,.613
!! nu_lop=64,2,1
!!/
!!
!!&h1_op_options
!! df_lop=.98,.564,.441,.363
!! nu_lop=64,4,2,1
!!/
!!
!!&hcurl_op_options
!! df_wop=.680,.384,.321
!! nu_wop=64,2,1
!!/
!!\endverbatim
!!
!!\subsection doc_marklin_ex2_input_mpi Parallel input file
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
!! meshname='hitsi'
!! cad_type=1
!! nlevels=2
!! nbase=1
!! grid_order=2
!!/
!!
!!&t3d_options
!! filename='hitsi_rs.t3d'
!! inpname='hitsi_rs.inp'
!! reflect='xy'
!!/
!!
!!&lag_op_options
!! df_lop=0.,1.,.837,.613
!! nu_lop=0,64,2,1
!!/
!!
!!&h1_op_options
!! df_lop=0.,.98,.564,.441,.363
!! nu_lop=0,64,4,2,1
!!/
!!
!!&hcurl_op_options
!! df_wop=0.,.680,.384,.321
!! nu_wop=0,64,2,1
!!/
!!\endverbatim
!!
