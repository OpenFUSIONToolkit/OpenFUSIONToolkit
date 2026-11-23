!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!!MUG Example: Cylindrical pinch kink {#doc_mug_cyl_kink}
!!============================
!!
!![TOC]
!!
!! This example demonstrates the use of the \ref xmhd "extended MHD" module in the Open FUSION Toolkit (OFT). In
!! this example we simulate an external kink mode in a cylindrical screw pinch.
!!
!!\section doc_mug_cyl_kink_code Code Walk Through
!!
!!The code consists of three basic sections, required imports and variable definitions,
!!finite element setup, and system creation and solution.
!!
!!\subsection doc_mug_cyl_kink_code_helper Helper Module
! START SOURCE
MODULE cyl_pinch
USE oft_base
USE oft_io, ONLY: hdf5_read, hdf5_write
USE oft_mesh_type, ONLY: oft_mesh
USE oft_la_base, ONLY: oft_vector
USE fem_base, ONLY: oft_fem_type
USE fem_utils, ONLY: fem_interp
USE oft_lag_operators, ONLY: oft_lag_vrinterp
USE oft_hcurl_grad_operators, ONLY: oft_hcurl_grad_rinterp
USE diagnostic, ONLY: field_probe
USE mhd_utils, ONLY: mu0
USE xmhd, ONLY: xmhd_sub_fields, oft_xmhd_probe, xmhd_probe_apply, xmhd_ML_lagrange, xmhd_ML_hcurl_grad
IMPLICIT NONE
!! Needs docs
!!
type, extends(fem_interp) :: cyl_pinch_init
  character(LEN=40) :: filename = 'none'
  character(LEN=1) :: field = 'n' !< Field component to initialize
  integer(4) :: npts = 0
  real(8), pointer, dimension(:,:) :: vals => NULL()
  integer(4) :: n = -1 !< Toroidal mode number
  integer(4) :: m = 2 !< Poloidal mode number
  real(8) :: phase0 = 0.d0 !< Mode phase
  real(8) :: rc = 1.0d0 !< Center of circulation layer
  real(8) :: rwidth = .2d0 !< Width of circulation layer
  real(8) :: rw = 2.0d0 !
  real(8) :: L = 2.0d0 !
contains
  !> Reconstruct fields"
  procedure :: setup => cyl_pinch_setup
  !> Reconstruct fields"
  procedure :: interp => cyl_pinch_apply
end type cyl_pinch_init
!! Needs docs
!!
TYPE, EXTENDS(oft_xmhd_probe) :: cyl_probe
  LOGICAL :: initialized = .FALSE. !< Flag to indicate setup has been called
  INTEGER(4) :: npts = 360
  TYPE(oft_hcurl_grad_rinterp), POINTER :: Bfield => NULL() !< Magnetic field interpolation class
  TYPE(oft_lag_vrinterp), POINTER :: Vfield => NULL()
  TYPE(field_probe) :: b_probes !< Synthetic surface magnetic probes
  TYPE(field_probe) :: v_probes !< Synthetic surface velocity probes
CONTAINS
  !> Extract probe signals
  PROCEDURE :: apply => cyl_probe_apply
END TYPE cyl_probe
REAL(8) :: rmax
CONTAINS
!! Needs docs
!!
SUBROUTINE cyl_pinch_setup(self,mesh)
CLASS(cyl_pinch_init), INTENT(inout) :: self
class(oft_mesh), target, intent(in) :: mesh
INTEGER(4) :: i,io_unit,iedge(1)
REAL(8) :: pmin,zmin,zmax
self%mesh=>mesh
!
zmin=MINVAL(self%mesh%r(3,:))
zmax=MAXVAL(self%mesh%r(3,:))
rmax=MAX(MAXVAL(self%mesh%r(1,:)),MAXVAL(self%mesh%r(2,:)))
zmin=-oft_mpi_max(-zmin)
zmax=oft_mpi_max(zmax)
rmax=oft_mpi_max(rmax)
self%L=(zmax-zmin)
!---Load GS field (order)
OPEN(NEWUNIT=io_unit,FILE=TRIM(self%filename))
READ(io_unit,*)self%npts
ALLOCATE(self%vals(4,self%npts))
DO i=1,self%npts
  READ(io_unit,*)self%vals(:,i)
!  WRITE(*,*)self%vals(:,i)
END DO
CLOSE(io_unit)
!---
IF(oft_env%head_proc)THEN
  iedge = MAXLOC(ABS(self%vals(2,:)))
  WRITE(*,'(2A)')'Cyl-Pinch Equilibrium loaded: ',self%filename
  WRITE(*,'(A,ES11.3)')'  L         =  ',self%L
  WRITE(*,'(A,ES11.3)')'  R_w       =  ',rmax
  WRITE(*,'(A,ES11.3)')'  B0        =  ',self%vals(3,self%npts)
  WRITE(*,'(A,ES11.3)')'  Ip        =  ',self%vals(2,iedge(1))*2.d0*pi*self%vals(1,iedge(1))/mu0
  WRITE(*,'(A,I8)')    '  # of points: ',self%npts
END IF
END SUBROUTINE cyl_pinch_setup
!! Needs docs
!!
subroutine cyl_pinch_apply(self,cell,f,gop,val)
class(cyl_pinch_init), intent(inout) :: self
integer(4), intent(in) :: cell
real(8), intent(in) :: f(:)
real(8), intent(in) :: gop(3,4)
real(8), intent(out) :: val(:)
real(8) :: pt(3),rad,theta,Bp,Bz,den,rhat(3),that(3),phase
! Map logical position to physical coordinates
pt=self%mesh%log2phys(cell,f)
rad = SQRT(pt(1)**2 + pt(2)**2)
theta = ATAN2(pt(2),pt(1))
Bp = linterp(self%vals(1,:),self%vals(2,:),self%npts,rad)
Bz = linterp(self%vals(1,:),self%vals(3,:),self%npts,rad)
den = linterp(self%vals(1,:),self%vals(4,:),self%npts,rad)
rhat = (/COS(theta),SIN(theta),0.d0/)
that = (/-SIN(theta),COS(theta),0.d0/)
! Return requested field evaluated at "(cell,f) -> pt"
SELECT CASE(self%field)
  CASE('n') ! Density
    val=den !1.d0
    IF(den==0.d0)val=1.d0
  CASE('t') ! Temperature
    val=1.d0
  CASE('v') ! Velocity
    rad = rad*self%rw/rmax
    ! val=0.d0 !1.d0
    !---Compute field
    phase = self%m*theta + self%n*2.d0*pi*(pt(3)+self%L/2.d0)/self%L + self%phase0
    ! val = rhat*COS(phase)*rad*EXP(-((rad-self%rc)**2)/self%rwidth)!/self%rwidth !(rad-self%rc)/self%rwidth!SIN((rad-self%rc)*2.d0*pi/self%rwidth)
    ! val = rhat*COS(phase)*EXP(-((rad-self%rc)**2)*2.d0/self%rwidth)
    ! val = rhat*COS(phase)*COS((rad-self%rc)*pi/self%rwidth+pi/2.d0)*self%rwidth/pi &
    !   + that*SIN(phase)*SIN((rad-self%rc)*pi/self%rwidth+pi/2.d0)/REAL(self%m,8)
    val = rhat*COS(phase)*(1.2d0-rad)*rad
    ! IF(ABS(rad-self%rc)>self%rwidth)val=0.d0
    IF(rad>1.2d0)val=0.d0
  CASE('b') ! Magnetic field
    val = (/-Bp*SIN(theta),Bp*COS(theta),Bz/)
  CASE('d') ! dB
    rad = rad*self%rw/rmax
    phase = self%m*theta + self%n*2.d0*pi*(pt(3)+self%L/2.d0)/self%L + self%phase0
    IF(rad<=1.d0)THEN
      val = rad**(self%m-1)*COS(phase)*rhat &
          - rad**(self%m-1)*SIN(phase)*that
      val = -REAL(self%m,8)*val
    ELSE
      val = (self%rw**self%m*rad**(-self%m-1) - self%rw**(-self%m)*rad**(self%m-1))*COS(phase)*rhat &
          + (self%rw**self%m*rad**(-self%m-1) + self%rw**(-self%m)*rad**(self%m-1))*SIN(phase)*that
      val = -REAL(self%m,8)*val/(self%rw**self%m - self%rw**(-self%m))
    END IF
  CASE DEFAULT
    CALL oft_abort('Unknown field component','cyl_pinch_apply',__FILE__)
END SELECT
end subroutine cyl_pinch_apply
!! Needs docs
!!
SUBROUTINE cyl_probe_apply(self,sub_fields,t)
CLASS(cyl_probe), intent(inout) :: self
type(xmhd_sub_fields), intent(inout) :: sub_fields ! Needs docs
REAL(8), intent(in) :: t
!---Local variables
INTEGER(4) :: i,npts
REAL(8) :: theta
REAL(8), ALLOCATABLE :: pts(:,:)
CHARACTER(LEN=40) :: file_tmp
!---------------------------------------------------------------------------
! Setup if necessary
!---------------------------------------------------------------------------
IF(.NOT.self%initialized)THEN
  !---
  npts=self%npts
  ALLOCATE(pts(3,npts))
  DO i=1,npts
    theta=i*2.d0*pi/REAL(npts,8)
    ! pts(:,i)=[1.0d0*COS(theta),1.0d0*SIN(theta),0.d0]
    pts(:,i)=[1.0d0*COS(theta),1.0d0*SIN(theta),0.d0]*rmax/2.d0
  END DO
  !---Setup V-probes
  ALLOCATE(self%Vfield)
  SELECT TYPE(this=>xmhd_ML_lagrange%current_level)
    CLASS IS(oft_fem_type)
      self%v_probes%mesh=>this%mesh
      self%b_probes%mesh=>this%mesh
  END SELECT
  CALL self%v_probes%setup(npts,pts)
  file_tmp='v_probes.hist'
  CALL self%v_probes%setup_save(file_tmp,'V')
  self%v_probes%B=>self%Vfield
  !---Setup B-probes
  ALLOCATE(self%Bfield)
  CALL self%b_probes%setup(npts,pts)
  file_tmp='b_probes.hist'
  CALL self%b_probes%setup_save(file_tmp,'B')
  self%b_probes%B=>self%Bfield
  !---
  self%initialized=.TRUE.
  DEALLOCATE(pts)
END IF
!---------------------------------------------------------------------------
! Setup interpolators
!---------------------------------------------------------------------------
self%Bfield%u=>sub_fields%B
CALL self%Bfield%setup(xmhd_ML_hcurl_grad%current_level)
self%Vfield%u=>sub_fields%V
CALL self%Vfield%setup(xmhd_ML_lagrange%current_level)
!---------------------------------------------------------------------------
! Sample signals
!---------------------------------------------------------------------------
CALL self%b_probes%save(t)
CALL self%v_probes%save(t)
END SUBROUTINE cyl_probe_apply
!! Needs docs
!!
FUNCTION vac_diff_shape(pt) result(val)
REAL(r8), INTENT(in) :: pt(3)
REAL(r8) :: val,rad
rad = magnitude(pt(1:2))*2.d0/rmax
! val = (1.d0+1.d3*(1.d0+TANH(50.d0*(magnitude(pt(1:2))-1.d0)))/2.d0)
val = 1.d0 + 1.d3*MAX(MIN((rad-1.0d0)/0.1d0,1.d0),0.d0)
END FUNCTION vac_diff_shape
END MODULE cyl_pinch
!!\subsection doc_mug_cyl_kink_code_prog Main Program
!! Driver program for MHD simulation of an external kink in a cylindrical pinch
!!
!! Runtime options are set in the main input file using the group
!! `cyl_pinch_options` group (see below).
PROGRAM cyl_xmhd
USE oft_base
USE oft_io, ONLY: xdmf_plot_file
!--Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre, create_bjacobi_pre, &
  create_ilu_pre, create_gmres_solver
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_getmop, oft_lag_project, &
  oft_lag_vgetmop, oft_lag_vproject, oft_vlag_zerob
!---H1(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
!---H1 FE space (Grad(H^1) subspace)
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp
!---Full H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_setup_interp, oft_hcurl_grad_rinterp, &
  oft_hcurl_grad_project, hcurl_grad_getmop, oft_hcurl_grad_gzerop
!---Physics
USE diagnostic, ONLY: scal_energy, vec_energy
USE mhd_utils, ONLY: elec_charge, proton_mass
USE xmhd, ONLY: xmhd_lin_run, xmhd_plot, xmhd_minlev, xmhd_bnorm_force, xmhd_taxis, &
  vel_scale, temp_floor, den_scale, den_floor, xmhd_sub_fields, &
  xmhd_ML_hcurl, xmhd_ML_H1, xmhd_ML_hcurl_grad, xmhd_ML_H1grad, xmhd_ML_lagrange, xmhd_ML_vlagrange, &
  res_shaping
!---Self
USE cyl_pinch, ONLY: cyl_pinch_init, cyl_probe, vac_diff_shape
IMPLICIT NONE
!---Metric solver
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
!---Local variables
TYPE(multigrid_mesh) :: mg_mesh
TYPE(xdmf_plot_file) :: plot_file
TYPE(xmhd_sub_fields) :: ic_fields,pert_fields
CLASS(oft_vector), POINTER :: u,v
TYPE(oft_hcurl_grad_rinterp), TARGET :: bfield
TYPE(oft_hcurl_grad_gzerop) :: hcurl_grad_gzerop
TYPE(oft_vlag_zerob) :: vlag_zerob
TYPE(cyl_pinch_init), TARGET :: init_field
TYPE(cyl_probe), TARGET :: probes
REAL(8), POINTER :: vals(:),vec_vals(:,:)
LOGICAL :: file_exist
LOGICAL :: pm = .FALSE.
REAL(8) :: V0 = 0.d0
INTEGER(4) :: ierr,io_unit,i
INTEGER(4) :: order = 2
INTEGER(4) :: minlev = 1
INTEGER(4) :: mm = 1
REAL(8) :: B0 = 1.d0
REAL(8) :: N0 = 1.d19
REAL(8) :: T0 = 1.d2
REAL(8) :: L0 = 1.d0
LOGICAL :: plot_run = .FALSE.
LOGICAL :: view_only = .FALSE.
CHARACTER(LEN=40) :: eq_file = 'none'
NAMELIST/cyl_pinch_options/order,minlev,B0,N0,T0,V0,L0,mm,plot_run,view_only,eq_file,pm
!---------------------------------------------------------------------------
! Initialize enviroment
!---------------------------------------------------------------------------
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,cyl_pinch_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(TRIM(eq_file)=='none')CALL oft_abort("Must specify equilibrium file","main",__FILE__)
!---------------------------------------------------------------------------
! Setup grid
!---------------------------------------------------------------------------
CALL multigrid_construct(mg_mesh,[0.d0,1.d0,0.d0])
IF((minlev<mg_mesh%mgdim).AND.(minlev/=-1))CALL oft_abort("Only polynomial MG currently supported for xMHD","main",__FILE__)
!---------------------------------------------------------------------------
! Build FE structures
!---------------------------------------------------------------------------
!--- Lagrange
ALLOCATE(xmhd_ML_lagrange,xmhd_ML_vlagrange)
CALL oft_lag_setup(mg_mesh,order,xmhd_ML_lagrange,ML_vlag_obj=xmhd_ML_vlagrange,minlev=minlev)
CALL lag_setup_interp(xmhd_ML_lagrange)
vlag_zerob%ML_vlag_rep=>xmhd_ML_vlagrange
!--- Grad(H^1) subspace
ALLOCATE(xmhd_ML_H1)
CALL oft_h1_setup(mg_mesh,order+1,xmhd_ML_H1,minlev=minlev)
CALL h1_setup_interp(xmhd_ML_H1)
!--- H(Curl) subspace
ALLOCATE(xmhd_ML_hcurl)
CALL oft_hcurl_setup(mg_mesh,order,xmhd_ML_hcurl,minlev=minlev)
CALL hcurl_setup_interp(xmhd_ML_hcurl)
!--- Full H(Curl) space
ALLOCATE(xmhd_ML_hcurl_grad,xmhd_ML_H1grad)
CALL oft_hcurl_grad_setup(xmhd_ML_hcurl,xmhd_ML_H1,xmhd_ML_hcurl_grad,xmhd_ML_H1grad,minlev=minlev)
CALL hcurl_grad_setup_interp(xmhd_ML_hcurl_grad,xmhd_ML_H1)
hcurl_grad_gzerop%ML_hcurl_grad_rep=>xmhd_ML_hcurl_grad
!---
init_field%filename=eq_file
! init_field%L=L0
init_field%m=mm
CALL init_field%setup(mg_mesh%mesh)
!---------------------------------------------------------------------------
! If plot run, just run and finish
!---------------------------------------------------------------------------
IF(plot_run)THEN
  CALL xmhd_plot(probes)
  CALL oft_finalize
END IF
!---Setup I/0
CALL plot_file%setup("cyl_pinch")
CALL mg_mesh%mesh%setup_io(plot_file,order)
!---------------------------------------------------------------------------
! Create scalar Lagrange mass matrix and solver
!---------------------------------------------------------------------------
!---Generate mass matrix
NULLIFY(mop) ! Ensure the matrix is unallocated (pointer is NULL)
CALL oft_lag_getmop(xmhd_ML_lagrange%current_level,mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
CALL create_diag_pre(minv%pre)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
! minv%pre=>minv_pre ! Set jacobi preconditioning (diagonal scaling)
!---Create fields for solver
CALL xmhd_ML_lagrange%vec_create(u) ! LHS solver field
CALL xmhd_ML_lagrange%vec_create(v) ! RHS solver field
!--------------------------------------------------------------------------
! Set up output arrays for plotting
!--------------------------------------------------------------------------
NULLIFY(vals)
ALLOCATE(vec_vals(3,u%n))
!---------------------------------------------------------------------------
! Set intial density from analytic definition
!---------------------------------------------------------------------------
!---Compute projection
init_field%field='n' ! Set desired field to poject
CALL oft_lag_project(xmhd_ML_lagrange%current_level,init_field,v) ! Project onto Lagrange basis
!---Retrieve field
CALL u%set(0.d0) ! Set initial guess for linear solve
CALL minv%apply(u,v) ! Invert mass matrix to retrieve field
!---Create density field and set values
CALL xmhd_ML_lagrange%vec_create(ic_fields%ne)
CALL ic_fields%ne%add(0.d0,1.d0,u)
! CALL ic_fields%ne%scale(N0) ! Scale to desired value
!---Plot field
CALL ic_fields%ne%get_local(vals) ! Fetch local values
CALL mg_mesh%mesh%save_vertex_scalar(vals,plot_file,'N0') ! Add field to plotting file
!---------------------------------------------------------------------------
! Set intial temperature from analytic definition
!---------------------------------------------------------------------------
! !---Compute projection
! init_field%field='t' ! Set desired field to poject
! CALL oft_lag_project(xmhd_ML_lagrange%current_level,init_field,v) ! Project onto Lagrange basis
! !---Retrieve field
! CALL u%set(0.d0) ! Set initial guess for linear solve
! CALL minv%apply(u,v) ! Invert mass matrix to retrieve field
! !---Create temperature field and set values
CALL xmhd_ML_lagrange%vec_create(ic_fields%Ti)
! CALL ic_fields%Ti%add(0.d0,1.d0,u)
! CALL ic_fields%Ti%scale(T0) ! Scale to desired value
CALL ic_fields%Ti%set(T0) ! Scale to desired value
!---Plot field
CALL ic_fields%Ti%get_local(vals) ! Fetch local values
CALL mg_mesh%mesh%save_vertex_scalar(vals,plot_file,'T0') ! Add field to plotting file
!---------------------------------------------------------------------------
! Delete fields, matrices and solvers for scalar Lagrange projection
!
! Note: Destory resests an object to a clean state so it may be used again for
! a different case. When an object is polymorphic (ie. a class) it must also be
! deallocated so that the pointer (object type) is in a clean state for reuse.
!---------------------------------------------------------------------------
CALL u%delete ! Destroy LHS vector
CALL v%delete ! Destroy RHS vector
CALL mop%delete ! Destroy mass matrix
DEALLOCATE(u,v,mop) ! Deallocate objects
! CALL minv%pre%delete ! Destory jacobi preconditioner
CALL minv%delete(.TRUE.) ! Destroy CG solver
DEALLOCATE(minv)
!---------------------------------------------------------------------------
! Create Nedelec H1 mass matrix and solver
!---------------------------------------------------------------------------
!---Generate mass matrix
NULLIFY(mop) ! Ensure the matrix is unallocated (pointer is NULL)
CALL hcurl_grad_getmop(xmhd_ML_hcurl_grad%current_level,mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_gmres_solver(minv,nrits=100)
CALL create_gmres_solver(minv%pre,nrits=2)
minv%pre%its=2
CALL create_bjacobi_pre(minv%pre%pre,5)
! CALL create_ilu_pre(minv%pre%pre%pre)
! CALL create_diag_pre(minv%pre%pre%pre)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
minv%atol=1.d-11
! ! minv%pre=>minv_pre ! Set jacobi preconditioning (diagonal scaling)
! CALL create_cg_solver(minv)
! CALL create_bjacobi_pre(minv%pre,1)
! CALL create_diag_pre(minv%pre)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
!---Create fields for solver
CALL xmhd_ML_hcurl_grad%vec_create(u) ! LHS solver field
CALL xmhd_ML_hcurl_grad%vec_create(v) ! RHS solver field
!---------------------------------------------------------------------------
! Set intial magnetic field from analytic definition
!---------------------------------------------------------------------------
INQUIRE(FILE="cyl_gs.rst", EXIST=file_exist)
IF(file_exist)THEN
  CALL xmhd_ML_hcurl_grad%current_level%vec_load(u, 'cyl_gs.rst', 'U')
ELSE
  init_field%field='b' ! Set desired field to poject
  CALL oft_hcurl_grad_project(xmhd_ML_hcurl_grad%current_level,init_field,v) ! Project onto Nedelec H1 basis
  CALL hcurl_grad_gzerop%apply(v) ! Zero out redundant vertex degrees of freedom
  !---Retrieve field
  CALL u%set(0.d0) ! Set initial guess for linear solve
  CALL minv%apply(u,v) ! Invert mass matrix to retrieve field
  CALL xmhd_ML_hcurl_grad%current_level%vec_save(u, 'cyl_gs.rst', 'U')
END IF
!---Create magnetic field and set values
CALL xmhd_ML_hcurl_grad%vec_create(ic_fields%B)
CALL ic_fields%B%add(0.d0,1.d0,u)
CALL ic_fields%B%scale(B0) ! Scale to desired value
! Note: plotting for the magnetic field is done later, see line ~342
IF(V0<0.d0)THEN
  !---------------------------------------------------------------------------
  ! Set intial magnetic field from analytic definition
  !---------------------------------------------------------------------------
  INQUIRE(FILE="cyl_db.rst", EXIST=file_exist)
  IF(file_exist)THEN
    CALL xmhd_ML_hcurl_grad%current_level%vec_load(u, 'cyl_db.rst', 'U')
  ELSE
    init_field%field='d' ! Set desired field to poject
    CALL oft_hcurl_grad_project(xmhd_ML_hcurl_grad%current_level,init_field,v) ! Project onto Nedelec H1 basis
    CALL hcurl_grad_gzerop%apply(v) ! Zero out redundant vertex degrees of freedom
    !---Retrieve field
    CALL u%set(0.d0) ! Set initial guess for linear solve
    CALL minv%apply(u,v) ! Invert mass matrix to retrieve field
    CALL xmhd_ML_hcurl_grad%current_level%vec_save(u, 'cyl_db.rst', 'U')
  END IF
  !---Create magnetic field and set values
  CALL xmhd_ML_hcurl_grad%vec_create(pert_fields%B)
  CALL pert_fields%B%add(0.d0,1.d0,u)
  CALL pert_fields%B%scale(1.d-3*V0) ! Scale to desired value
  ! Note: plotting for the magnetic field is done later, see line ~342
ELSE
  CALL xmhd_ML_hcurl_grad%vec_create(pert_fields%B)
END IF
!---------------------------------------------------------------------------
! Delete fields, matrices and solvers for Nedelec H1 projection
!---------------------------------------------------------------------------
CALL u%delete ! Destroy LHS vector
CALL v%delete ! Destroy RHS vector
CALL mop%delete ! Destroy mass matrix
DEALLOCATE(u,v,mop) ! Deallocate objects
! CALL minv%pre%pre%pre%delete
! CALL minv%pre%pre%delete
! DEALLOCATE(minv%pre%pre)
! CALL minv%pre%delete ! Destory jacobi preconditioner
! DEALLOCATE(minv%pre)
CALL minv%delete(.TRUE.) ! Destroy CG solver
DEALLOCATE(minv)
!---------------------------------------------------------------------------
! Create vector Lagrange mass matrix and solver
!---------------------------------------------------------------------------
!---Generate mass matrix
NULLIFY(mop) ! Ensure the matrix is unallocated (pointer is NULL)
CALL oft_lag_vgetmop(xmhd_ML_vlagrange%current_level,mop,"none") ! Construct mass matrix with "none" BC
!---Setup linear solver
CALL create_cg_solver(minv)
CALL create_diag_pre(minv%pre)
minv%A=>mop ! Set matrix to be solved
minv%its=-2 ! Set convergence type (in this case "full" CG convergence)
! minv%pre=>minv_pre ! Set jacobi preconditioning (diagonal scaling)
!---Create fields for solver
CALL xmhd_ML_vlagrange%vec_create(u) ! LHS solver field
CALL xmhd_ML_vlagrange%vec_create(v) ! RHS solver field
!---------------------------------------------------------------------------
! Set intial velocity field from analytic definition
!---------------------------------------------------------------------------
IF(V0>0.d0)THEN
  init_field%field='v' ! Set desired field to poject
  CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,init_field,v) ! Project onto vector Lagrange basis
  !---Retrieve field
  CALL u%set(0.d0) ! Set initial guess for linear solve
  CALL minv%apply(u,v) ! Invert mass matrix to retrieve field
  CALL vlag_zerob%apply(u)
  !---Create velocity field and set values
  CALL xmhd_ML_vlagrange%vec_create(pert_fields%V)
  CALL pert_fields%V%add(0.d0,1.d0,u)
  CALL pert_fields%V%scale(V0) ! Scale to desired value
ELSE
  CALL xmhd_ML_vlagrange%vec_create(pert_fields%V)
END IF
!---Plot field, this time we must extract each direction into a single 2D
!   array for plotting. This is done by using the "iblock" argument to "get_local".
!   Also, notice that we preallocated the "vec_vals" array.
DO i=1,3
  vals=>vec_vals(i,:) ! Set retrieval array to a slice of the full array
  CALL pert_fields%V%get_local(vals,i) ! Fetch local values for a given direction
END DO
CALL mg_mesh%mesh%save_vertex_vector(vec_vals,plot_file,'V0') ! Add field to plotting file
!---------------------------------------------------------------------------
! Project B for plotting
!
! As the magnetic field is not represented in a nodal basis set it cannot be
! plotting directly. So we first project the field onto a vector Lagrange
! basis and then plot the resulting field.
!---------------------------------------------------------------------------
!---Setup interpolation object for H1 field
bfield%u=>ic_fields%B ! Set target field to previously computed magnetic field
CALL bfield%setup(xmhd_ML_hcurl_grad%current_level) ! Setup interpolation object
!--- Project onto vector Lagrange basis
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,bfield,v)
!---Retrieve field
CALL u%set(0.d0) ! Set initial guess for linear solve
CALL minv%apply(u,v) ! Invert mass matrix to retrieve field
!---Plot field
DO i=1,3
  vals=>vec_vals(i,:) ! Set retrieval array to a slice of the full array
  CALL u%get_local(vals,i) ! Fetch local values for a given direction
END DO
CALL mg_mesh%mesh%save_vertex_vector(vec_vals,plot_file,'B0') ! Add field to plotting file
!---Setup interpolation object for H1 field
bfield%u=>pert_fields%B ! Set target field to previously computed magnetic field
CALL bfield%setup(xmhd_ML_hcurl_grad%current_level) ! Setup interpolation object
!--- Project onto vector Lagrange basis
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,bfield,v)
!---Retrieve field
CALL u%set(0.d0) ! Set initial guess for linear solve
CALL minv%apply(u,v) ! Invert mass matrix to retrieve field
!---Plot field
DO i=1,3
  vals=>vec_vals(i,:) ! Set retrieval array to a slice of the full array
  CALL u%get_local(vals,i) ! Fetch local values for a given direction
END DO
CALL mg_mesh%mesh%save_vertex_vector(vec_vals,plot_file,'dB') ! Add field to plotting file
!---------------------------------------------------------------------------
! Delete fields, matrices and solvers for vector Lagrange projection
!---------------------------------------------------------------------------
CALL u%delete ! Destroy LHS vector
CALL v%delete ! Destroy RHS vector
CALL mop%delete ! Destroy mass matrix
DEALLOCATE(u,v,mop) ! Deallocate objects
! CALL minv%pre%delete ! Destory jacobi preconditioner
CALL minv%delete(.TRUE.) ! Destroy CG solver
DEALLOCATE(minv)
IF(view_only)CALL oft_finalize
!---------------------------------------------------------------------------
! Run simulation
!---------------------------------------------------------------------------
xmhd_minlev=minlev  ! Set minimum level for multigrid preconditioning
xmhd_taxis=3        ! Set axis of rotation for toroidal flux/current
temp_floor=T0*1.d-2 ! Set temperature floor
den_floor=N0*1.d-2  ! Set density floor
den_scale=N0        ! Set density equation scale factor
!---Set custom diffusivity shaping in vacuum
res_shaping=>vac_diff_shape
!---Run simulation
! CALL b%add(1.d0,1.d0,db)
! CALL xmhd_run(b,vel,n,temp)
CALL xmhd_ML_vlagrange%vec_create(ic_fields%V)
CALL xmhd_ML_lagrange%vec_create(pert_fields%Ne)
CALL xmhd_ML_lagrange%vec_create(pert_fields%Ti)
oft_env%pm=pm
CALL xmhd_lin_run(ic_fields,pert_fields)!,escale=2.d-6)
!---Finalize enviroment
CALL oft_finalize
END PROGRAM cyl_xmhd
! STOP SOURCE
!!
!!\section doc_mug_cyl_kink_input Input file
!!
!! Below is an input file which can be used with this example in a parallel environment.
!! This example should only be run with multiple processes.
!! Some annotation of the options is provided inline below, for more information on the
!! available options in the \c xmhd_options group see \ref xmhd::xmhd_plot "xmhd_plot".
!!
! OFT_DOC_INCLUDE: oft.in
!!
!!\subsection doc_mug_cyl_kink_input_solver Solver specification
!!
!! Time dependent MHD solvers are accelerated significantly by the use of
!! a more sophisticated preconditioner than the default method. Below is
!! an example `oft_in.xml` file that constructs an appropriate MG preconditioner.
!! Currently, this preconditioner method is the suggest preconditioner for all
!! time-dependent MHD solves.
!!
!! This solver can be used by specifying both the FORTRAN input and XML input files
!! to the executable as below.
!!
!!\verbatim
!!~$ ./MUG_cyl_pinch oft.in oft_in.xml
!!\endverbatim
!!
! OFT_DOC_INCLUDE: oft_in.xml 
!!
!!\section doc_mug_cyl_kink_mesh Mesh Creation
!! A mesh file `MUG-cyl_pinch.h5` is provided with this example. A suitable mesh for this example,
!! with radius of 0.2m and height of 0.1m, can be created using the following [CUBIT](https://cubit.sandia.gov/) script.
!!
! OFT_DOC_INCLUDE: cyl_pinch-no_wall.jou
!!
!! Once complete, the mesh should be converted into the native mesh format using the `convert_cubit.py` script as
!! below. The script is located in `bin` following installation or `src/utilities` in the base repo.
!!
!!\verbatim
!!~$ python convert_cubit.py --in_file=MUG-cyl_pinch.g --periodic_nodeset=1
!!\endverbatim
!!
!! @note The option `--periodic_nodeset=1` indicates that the surface corresponding to the first, and only, nodeset
!! should be treated as periodic following reflection