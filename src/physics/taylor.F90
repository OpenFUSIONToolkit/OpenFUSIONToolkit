!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file taylor.f90
!
!> @defgroup doxy_oft_physics Physics
!! Physics modules for the Open FUSION Toolkit
!
!> Subroutines and fields for Taylor state calculations using
!! mimetic operators.
!! - Force-Free eigenmodes
!! - Vacuum fields for geometries with cut planes
!! - Inhomogeneous force free states from vacuum fields
!!
!! @author Chris Hansen
!! @date June 2010
!! @ingroup doxy_oft_physics
!---------------------------------------------------------------------------
module taylor
USE oft_base
USE oft_io
USE oft_mesh_type, ONLY: oft_mesh
!---
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, oft_matrix, oft_matrix_ptr
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_native_cg_solver, oft_native_cg_eigsolver, &
  oft_native_gmres_solver, oft_ml_precond
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---
USE fem_base, ONLY: oft_fem_type, oft_afem_type, oft_ml_fem_type
USE fem_utils, ONLY: fem_interp
USE fem_composite, ONLY: oft_fem_comp_type, oft_ml_fem_comp_type
USE oft_lag_operators, ONLY: oft_lag_zerob, lag_getlop_pre, oft_lag_getlop
USE oft_h1_basis, ONLY: oft_h1_geval_all, oft_h1_fem
USE oft_h1_operators, ONLY: oft_h1_zerogrnd, h1_getlop_pre, oft_h1_getlop
USE oft_hcurl_basis, ONLY: oft_hcurl_eval_all, oft_hcurl_ceval_all, &
  oft_hcurl_get_cgops, oft_hcurl_fem
USE oft_hcurl_operators, ONLY: oft_hcurl_cinterp, oft_hcurl_orthog, &
  oft_hcurl_divout, hcurl_getwop_pre, oft_hcurl_zerob, oft_hcurl_getmop, oft_hcurl_getkop, &
  oft_hcurl_getwop, oft_hcurl_getjmlb, hcurl_getjmlb_pre
USE oft_hcurl_grad_operators, ONLY: oft_hcurl_grad_divout, hcurl_grad_getmop, hcurl_grad_mc
!---
USE diagnostic, ONLY: tfluxfun
implicit none
#include "local.h"
!---------------------------------------------------------------------------
!> Interpolate a Taylor state field
!!
!! Taylor state fields consist of a gradient component, the vacuum field, and
!! a curl component, the plasma field.
!---------------------------------------------------------------------------
type, extends(fem_interp) :: oft_taylor_rinterp
  class(oft_vector), pointer :: ua => NULL() !< Plasma vector potential
  class(oft_vector), pointer :: uvac => NULL() !< Vacuum magnectic field
  real(r8), pointer, dimension(:) :: vac_grad => NULL() !< Local vacuum field (gradient)
  real(r8), pointer, dimension(:) :: vac_curl => NULL() !< Local vacuum field (curl)
  real(r8), pointer, dimension(:) :: acurl => NULL() !< Local vector potential
  class(oft_h1_fem), pointer :: grad_rep => NULL() !< Grad(H^1) FE representation
  class(oft_hcurl_fem), pointer :: curl_rep => NULL() !< H(Curl) FE representation
contains
  !> Retrieve local values for interpolation
  generic :: setup => setup1, setup2
  procedure :: setup1 => taylor_rinterp_setup1
  procedure :: setup2 => taylor_rinterp_setup2
  !> Reconstruct field
  procedure :: interp => taylor_rinterp
end type oft_taylor_rinterp
!---------------------------------------------------------------------------
!> Force-free eigenmode object
!---------------------------------------------------------------------------
type, public :: oft_taylor_eigs
  integer(i4) :: minlev = -1 !< Lowest FE level for MG solvers
  integer(i4) :: nm = 0 !< Number of force-free fields to be computed
  real(r8), pointer, dimension(:,:) :: hlam => NULL() !< Homogeneous force-free lambdas
  real(r8), pointer, dimension(:,:) :: htor => NULL() !< Homogeneous force-free toroidal fluxes
  type(oft_vector_ptr), pointer, dimension(:,:) :: hffa => NULL() !< Homogeneous force-free fields
  type(oft_hcurl_orthog), pointer :: orthog => NULL()
  TYPE(oft_ml_fem_type), POINTER :: ML_hcurl => NULL()
  TYPE(oft_ml_fem_type), POINTER :: ML_lagrange => NULL()
CONTAINS
  PROCEDURE :: setup => eigs_setup
end type oft_taylor_eigs
integer(i4), parameter :: taylor_tag_size = 4 !< Size of  jump planes character tags
!---------------------------------------------------------------------------
!> Force-free vacuum source object
!---------------------------------------------------------------------------
type, public :: oft_taylor_inhomo
  integer(i4) :: minlev = -1 !< Lowest FE level for MG solvers
  integer(i4) :: nh = 2 !< Number of jump planes in current geometry
  real(r8)  :: jtol = 1.d-6 !< Tolerance for identifying edges on jump plane
  real(r8), pointer, dimension(:,:)  :: hcpc => NULL() !< Center points of jump planes
  real(r8), pointer, dimension(:,:)  :: hcpv => NULL() !< Normal vectors for jump planes
  character(LEN=taylor_tag_size), pointer, dimension(:) :: htag => NULL() !< Injector names
  type(oft_vector_ptr), pointer, dimension(:,:) :: hvac => NULL() !< Vacuum magnetic fields
  type(oft_vector_ptr), pointer, dimension(:,:) :: hcur => NULL() !< Inhomogeneous source fields
  type(oft_vector_ptr), pointer, dimension(:,:) :: gffa => NULL() !< Inhomogeneous force-free fields
  TYPE(oft_ml_fem_type), POINTER :: ML_lagrange => NULL()
  TYPE(oft_ml_fem_type), POINTER :: ML_h1 => NULL()
  TYPE(oft_ml_fem_type), POINTER :: ML_h1grad => NULL()
  TYPE(oft_ml_fem_type), POINTER :: ML_hcurl => NULL()
  TYPE(oft_ml_fem_comp_type), POINTER :: ML_hcurl_grad => NULL()
CONTAINS
  PROCEDURE :: setup => ff_setup
end type oft_taylor_inhomo
contains
!---------------------------------------------------------------------------
!> Setup eigenmodes object
!---------------------------------------------------------------------------
subroutine eigs_setup(self,ML_hcurl,ML_lagrange,minlev)
class(oft_taylor_eigs), intent(inout) :: self
TYPE(oft_ml_fem_type), TARGET, INTENT(in) :: ML_hcurl
TYPE(oft_ml_fem_type), TARGET, INTENT(in) :: ML_lagrange
integer(i4), optional, intent(in) :: minlev
self%ML_hcurl=>ML_hcurl
self%ML_lagrange=>ML_lagrange
!
self%minlev=-1
IF(PRESENT(minlev))self%minlev=minlev
IF(self%minlev<0)self%minlev=self%ML_hcurl%level
end subroutine eigs_setup
!---------------------------------------------------------------------------
!> Setup eigenmodes object
!---------------------------------------------------------------------------
subroutine ff_setup(self,ML_hcurl,ML_h1,ML_hcurl_grad,ML_h1grad,ML_lagrange,minlev)
class(oft_taylor_inhomo), intent(inout) :: self
TYPE(oft_ml_fem_type), TARGET, INTENT(in) :: ML_hcurl
TYPE(oft_ml_fem_type), TARGET, INTENT(in) :: ML_h1
TYPE(oft_ml_fem_comp_type), TARGET, INTENT(in) :: ML_hcurl_grad
TYPE(oft_ml_fem_type), TARGET, INTENT(in) :: ML_h1grad
TYPE(oft_ml_fem_type), TARGET, INTENT(in) :: ML_lagrange
integer(i4), optional, intent(in) :: minlev
self%ML_hcurl=>ML_hcurl
self%ML_h1=>ML_h1
self%ML_hcurl_grad=>ML_hcurl_grad
self%ML_h1grad=>ML_h1grad
self%ML_lagrange=>ML_lagrange
!
self%minlev=-1
IF(PRESENT(minlev))self%minlev=minlev
IF(self%minlev<0)self%minlev=self%ML_hcurl%level
end subroutine ff_setup
!---------------------------------------------------------------------------
!> Compute 'taylor_nm' Force-Free eignemodes.
!!
!! @note When `taylor_rst=.TRUE.` one restart files will be generated for
!! each computed mode on each MG level `hffa_*.rst`
!---------------------------------------------------------------------------
subroutine taylor_hmodes(self,nm,rst_filename)
type(oft_taylor_eigs), intent(inout) :: self
integer(i4), optional, intent(in) :: nm !< Number of modes to compute (optional: 1)
character(LEN=OFT_PATH_SLEN), optional, intent(in) :: rst_filename
class(oft_vector), pointer :: u,tmp
!--- Taylor eigenvalue solver
TYPE(oft_native_cg_eigsolver) :: eigsolver
!--- Divergence cleaner
CLASS(oft_solver), POINTER :: linv => NULL()
TYPE(oft_hcurl_divout) :: divout
!--- ML structures for MG-preconditioner
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_wop => NULL()
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_lop => NULL()
INTEGER(i4) :: nlevels_wop,nlevels_lop
!--- Matrix pointers
class(oft_matrix), pointer :: kop => NULL()
class(oft_matrix), pointer :: wop => NULL()
class(oft_matrix), pointer :: lop => NULL()
!--- Local variables
type(oft_hcurl_cinterp) :: Bfield
real(r8) :: alam,elapsed_time
integer(i4) :: i,j,k,ierr
character(LEN=16) :: field_name
logical :: save_rst
type(oft_timer) :: mytimer
type(oft_hcurl_zerob), target :: hcurl_zerob
type(oft_lag_zerob), target :: lag_zerob
!---
DEBUG_STACK_PUSH
IF(PRESENT(nm))THEN
  IF(nm<=0.OR.nm>20)CALL oft_abort("Invalid number of modes requested.", "taylor_hmodes", __FILE__)
  self%nm=nm
ELSE
  self%nm=1
END IF
save_rst=.FALSE.
IF(PRESENT(rst_filename))save_rst=.TRUE.
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(A)')'============================'
  WRITE(*,'(A)')'Starting calculation of Taylor states'
  WRITE(*,'(A)')'============================'
  WRITE(*,*)
  CALL mytimer%tick
END IF
!---Allocate storage
ALLOCATE(self%hffa(self%nm,self%ML_hcurl%level))
ALLOCATE(self%hlam(self%nm,self%ML_hcurl%level))
ALLOCATE(self%htor(self%nm,self%ML_hcurl%level))
!---Setup orthogonalization
ALLOCATE(self%orthog)
self%orthog%ML_hcurl_rep=>self%ML_hcurl
self%orthog%orthog=>self%hffa
!---------------------------------------------------------------------------
! Create ML Matrices
!---------------------------------------------------------------------------
nlevels_wop=self%ML_hcurl%nlevels-self%minlev+1
nlevels_lop=self%ML_lagrange%ml_mesh%mgdim-self%minlev+1
ALLOCATE(ml_wop(nlevels_wop),ml_lop(nlevels_lop))
DO i=1,nlevels_wop
  CALL self%ML_hcurl%set_level(self%minlev+i-1)
  NULLIFY(ml_wop(i)%M)
  CALL oft_hcurl_getwop(self%ML_hcurl%current_level,ml_wop(i)%M,'zerob')
  IF(self%ML_hcurl%level<=self%ML_lagrange%ml_mesh%mgdim)THEN
    CALL self%ML_lagrange%set_level(self%ML_hcurl%level)
    NULLIFY(ml_lop(i)%M)
    CALL oft_lag_getlop(self%ML_lagrange%current_level,ml_lop(i)%M,"zerob")
  END IF
END DO
CALL self%ML_hcurl%set_level(self%ML_hcurl%level)
!---------------------------------------------------------------------------
! Loop over desired number of modes
!---------------------------------------------------------------------------
hcurl_zerob%ML_hcurl_rep=>self%ML_hcurl
lag_zerob%ML_lag_rep=>self%ML_lagrange
do k=1,self%nm
  self%orthog%nm=k-1
  call self%ML_hcurl%set_level(self%minlev)
  !---Loop over levels for each mode
  do i=self%minlev,self%ML_hcurl%nlevels
    !---Build level field
    call self%ML_hcurl%set_level(i)
    call self%ML_hcurl%vec_create(self%hffa(k,i)%f,i)
    !---Alias to general field
    u=>self%hffa(k,i)%f
    if((i>self%minlev).AND.(self%ML_hcurl%level==self%ML_hcurl%blevel+1))cycle
    !---Setup Solver
    wop=>ml_wop(i-self%minlev+1)%M
    self%orthog%wop=>wop
    NULLIFY(kop)
    CALL oft_hcurl_getkop(self%ML_hcurl%current_level,kop,"zerob")
!---------------------------------------------------------------------------
! Create eigenvalue solver
!---------------------------------------------------------------------------
    eigsolver%A=>wop
    eigsolver%M=>kop
    eigsolver%its=-2
    eigsolver%bc=>hcurl_zerob
    eigsolver%orthog=>self%orthog
    IF(k>1)THEN
      eigsolver%nrestarts=4
      IF(i==self%minlev)THEN
        eigsolver%ninner=1000
      ELSE
        eigsolver%ninner=100
      END IF
    END IF
!---------------------------------------------------------------------------
! Setup preconditioner
!---------------------------------------------------------------------------
    if(i==self%minlev)then ! Lowest level uses diag precond
      !---Setup Preconditioner
      CALL create_diag_pre(eigsolver%pre)
      call u%set(1.d0,random=.TRUE.) ! Initialize guess
    else ! Higher levels use MG
      CALL hcurl_getwop_pre(self%ML_hcurl,eigsolver%pre,ml_wop,nlevels=i-self%minlev+1)
      SELECT TYPE(this=>eigsolver%pre)
        CLASS IS(oft_ml_precond)
          call self%ML_hcurl%set_level(i-1)
          call this%ml_vecspace%interp(self%hffa(k,i-1)%f,u)
        CLASS DEFAULT
          CALL oft_abort("Invalid ML preconditioner","self%hmodes",__FILE__)
      END SELECT
      ! call hcurl_interp(self%hffa(k,i-1)%f,u)
      if(k/=1)call u%set(1.d0,random=.TRUE.)
    end if
!---------------------------------------------------------------------------
! Initialize guess
!---------------------------------------------------------------------------
    ! if(i==self%minlev)then
    !   call u%set(1.d0,random=.TRUE.)
    ! else
    !   call hcurl_interp(self%hffa(k,i-1)%f,u)
    !   if(k/=1)call u%set(1.d0,random=.TRUE.)
    !   if(ML_oft_hcurl%level==ML_oft_hcurl%blevel+1)cycle
    ! end if
    if(save_rst)then
      if(oft_file_exist(rst_filename))then
        WRITE(field_name,'(A6,A2,A2,I2.2,A2,I2.2)')'hffa_g',self%ML_hcurl%ml_mesh%rlevel,'_p',self%ML_hcurl%current_level%order,'_m',k
        WRITE(*,*)field_name
        IF(hdf5_field_exist(rst_filename,field_name))CALL self%ML_hcurl%current_level%vec_load(u,rst_filename,field_name)
      end if
    end if
!---------------------------------------------------------------------------
! Solve
!---------------------------------------------------------------------------
    CALL hcurl_zerob%apply(u) ! Apply BC
    CALL eigsolver%apply(u,self%hlam(k,i))
    CALL eigsolver%pre%delete
    DEALLOCATE(eigsolver%pre)
    CALL eigsolver%delete
    CALL kop%delete
    DEALLOCATE(kop)
    !---
    Bfield%u=>u
    CALL Bfield%setup(self%ML_hcurl%current_level)
    self%htor(k,i) = tfluxfun(Bfield%mesh,Bfield,self%ML_hcurl%current_level%quad%order)
    CALL Bfield%delete
!---------------------------------------------------------------------------
! Create divergence cleaner
!---------------------------------------------------------------------------
    CALL self%ML_lagrange%set_level(MIN(self%ML_hcurl%level,self%ML_lagrange%ml_mesh%mgdim))
    IF(nlevels_lop<=0)THEN
      NULLIFY(lop)
      CALL oft_lag_getlop(self%ML_lagrange%current_level,lop,"zerob")
    ELSE
      lop=>ml_lop(self%ML_lagrange%level-self%minlev+1)%M
    END IF
    if(self%ML_lagrange%level<=self%minlev)then ! Lowest level uses diag precond
      CALL create_cg_solver(linv)
      CALL create_diag_pre(linv%pre)
    else ! Higher levels use MG
      CALL create_cg_solver(linv, force_native=.TRUE.)
      CALL lag_getlop_pre(self%ML_lagrange,linv%pre,ml_lop,nlevels=self%ML_lagrange%level-self%minlev+1)
    end if
    linv%A=>lop
    linv%its=-2
    CALL divout%setup(self%ML_hcurl,self%ML_lagrange,bc='zero',solver=linv)
    ! divout%ML_hcurl_rep=>ML_oft_hcurl
    ! divout%ML_lag_rep=>ML_oft_lagrange
    ! divout%solver=>linv
    ! divout%bc=>lag_zerob
    divout%pm=oft_env%pm
    CALL divout%apply(u)
    CALL linv%pre%delete
    DEALLOCATE(linv%pre)
    CALL linv%delete
    DEALLOCATE(linv)
    IF(nlevels_lop<=0)CALL lop%delete
!---------------------------------------------------------------------------
! Write restart files
!---------------------------------------------------------------------------
    if(save_rst)then
      if(.NOT.oft_file_exist(rst_filename))then
        CALL oft_mpi_barrier(ierr)
        WRITE(field_name,'(A6,A2,A2,I2.2,A2,I2.2)')'hffa_g',self%ML_hcurl%ml_mesh%rlevel,'_p',self%ML_hcurl%current_level%order,'_m',k
        CALL self%ML_hcurl%current_level%vec_save(u,rst_filename,field_name)
      end if
    end if
  end do ! End level loop
end do ! End mode loop
if(oft_env%head_proc)then
  DO i=self%minlev,self%ML_hcurl%nlevels
    WRITE(*,'(2X,A,I3)')'Level =',i
    DO k=1,self%nm
      WRITE(*,'(4X,A,I4,A,ES14.6)')'Mode =',k,'   Lambda = ',self%hlam(k,i)
    END DO
  END DO
  elapsed_time=mytimer%tock()
  WRITE(*,*)
  WRITE(*,'(2X,A,F12.3)')'Time Elapsed = ',elapsed_time
end if
!---------------------------------------------------------------------------
! Deallocate ML matrices
!---------------------------------------------------------------------------
DO i=1,nlevels_wop
  CALL ml_wop(i)%M%delete
  DEALLOCATE(ml_wop(i)%M)
  IF(i<nlevels_lop)THEN
    CALL ml_lop(i)%M%delete
    DEALLOCATE(ml_lop(i)%M)
  END IF
END DO
DEALLOCATE(ml_wop,ml_lop)
CALL self%ML_lagrange%set_level(self%ML_lagrange%nlevels)
DEBUG_STACK_POP
end subroutine taylor_hmodes
!---------------------------------------------------------------------------
!> Generate vacuum fields for a geometry with cut planes.
!!
!! @note When `taylor_rst=.TRUE.` two restart files will be generated for
!! each jump plane `hvac_*.rst` and `hcur_*.rst`
!---------------------------------------------------------------------------
subroutine taylor_vacuum(self,nh,hcpc,hcpv,htags,energy,hmodes,rst_filename)
type(oft_taylor_inhomo), intent(inout) :: self
integer(i4), intent(in) :: nh !< Number of jump planes
real(r8), intent(in) :: hcpc(3,nh) !< Jump plane center possitions [3,nh]
real(r8), intent(in) :: hcpv(3,nh) !< Jump plane normal vectors [3,nh]
character(LEN=taylor_tag_size), optional, intent(in) :: htags(nh) !< Names for each jump plane [LEN=taylor_tag_size,nh] (optional)
real(r8), optional, intent(out) :: energy(nh) !< Vacuum energy for each jump plan (optional)
type(oft_taylor_eigs), optional, intent(inout) :: hmodes
character(LEN=OFT_PATH_SLEN), optional, intent(in) :: rst_filename
!---H(Curl) full space divergence cleaner
CLASS(oft_solver), POINTER :: linv => NULL()
TYPE(oft_hcurl_grad_divout) :: divout
!---WOP solver
CLASS(oft_solver), POINTER :: winv => NULL()
!---H(Curl) subspace divergence cleaner
CLASS(oft_solver), POINTER :: linv_lag => NULL()
TYPE(oft_hcurl_divout), TARGET :: hcurl_divout
!--- ML structures for MG-preconditioner
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_int => NULL()
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_lop => NULL()
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_wop => NULL()
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: levels,nu
REAL(r8), ALLOCATABLE, DIMENSION(:) :: df
!---
class(oft_matrix), pointer :: lop => NULL()
class(oft_matrix), pointer :: lop_lag => NULL()
class(oft_matrix), pointer :: mop => NULL()
class(oft_matrix), pointer :: mop_hcurl => NULL()
class(oft_matrix), pointer :: wop => NULL()
!---
class(oft_vector), pointer :: u,b,tmp
real(r8), pointer, dimension(:) :: vals => NULL()
real(r8) :: alam,venergy
integer(i4) :: i,j,k,nlevels,ierr
logical :: have_rst,save_rst
character(LEN=16) :: field_name
type(oft_hcurl_zerob), target :: hcurl_zerob
type(oft_lag_zerob), target :: lag_zerob
DEBUG_STACK_PUSH
NULLIFY(lop,lop_lag,mop,mop_hcurl,wop)
NULLIFY(ml_int,ml_lop,ml_wop)
!---Create taylor module variables
self%nh=nh
ALLOCATE(self%hcpc(3,self%nh),self%hcpv(3,self%nh))
self%hcpc=hcpc
self%hcpv=hcpv
ALLOCATE(self%htag(self%nh))
IF(PRESENT(htags))THEN
  self%htag=htags
ELSE
  DO i=1,self%nh
    WRITE(self%htag(i),'(A3,I1)')'INJ',i
  END DO
END IF
hcurl_zerob%ML_hcurl_rep=>self%ML_hcurl
lag_zerob%ML_lag_rep=>self%ML_lagrange
!---Loop over cut planes
save_rst=.FALSE.
IF(PRESENT(rst_filename))save_rst=.TRUE.
IF(save_rst)THEN
  have_rst=oft_file_exist(TRIM(rst_filename))
  DO i=1,self%nh
    WRITE(field_name,'(A6,A2,A2,I2.2,A2,I2.2)')'hvac_g',self%ML_hcurl%ml_mesh%rlevel,'_p',self%ML_hcurl%current_level%order,'_h',i
    have_rst=have_rst.AND.hdf5_field_exist(rst_filename,field_name)
    WRITE(field_name,'(A6,A2,A2,I2.2,A2,I2.2)')'hcur_g',self%ML_hcurl%ml_mesh%rlevel,'_p',self%ML_hcurl%current_level%order,'_h',i
    have_rst=have_rst.AND.hdf5_field_exist(rst_filename,field_name)
  END DO
ELSE
  have_rst=.FALSE.
END IF
IF(.NOT.have_rst)THEN
!---------------------------------------------------------------------------
! Setup H^1::LOP preconditioner
!---------------------------------------------------------------------------
  if(self%minlev==self%ML_h1%nlevels-1)then ! Lowest level uses diag precond
    CALL oft_h1_getlop(self%ML_h1%current_level,lop,'grnd')
    CALL create_cg_solver(linv)
    CALL create_diag_pre(linv%pre)
  else ! Nested levels use MG
    CALL create_cg_solver(linv, force_native=.TRUE.)
    CALL h1_getlop_pre(self%ML_h1,linv%pre,ml_lop,'grnd',nlevels=self%ML_h1%nlevels-self%minlev+1)
      lop=>ml_lop(self%ML_h1%nlevels-self%minlev+1)%M
  end if
!---------------------------------------------------------------------------
! Create divergence cleaner
!---------------------------------------------------------------------------
  linv%A=>lop
  linv%its=-2
  ! divout%solver=>linv
  CALL divout%setup(self%ML_hcurl_grad,'grnd',solver=linv)
ELSE
  CALL create_cg_solver(linv)
  CALL create_diag_pre(linv%pre)
END IF
!---------------------------------------------------------------------------
! Setup H(Curl)::WOP preconditioner
!---------------------------------------------------------------------------
CALL create_cg_solver(winv, force_native=.TRUE.)
IF((self%minlev==self%ML_hcurl%level).OR.have_rst)THEN ! Lowest level uses diag precond
  NULLIFY(wop)
  CALL oft_hcurl_getwop(self%ML_hcurl%current_level,wop,'zerob')
  CALL create_diag_pre(winv%pre)
ELSE ! Nested levels use MG
  CALL hcurl_getwop_pre(self%ML_hcurl,winv%pre,ml_wop,nlevels=self%ML_hcurl%level-self%minlev+1)
  wop=>ml_wop(self%ML_hcurl%level-self%minlev+1)%M
END IF
!---------------------------------------------------------------------------
! Create H(Curl)::WOP solver
!---------------------------------------------------------------------------
winv%A=>wop
winv%its=-3
winv%atol=1.d-9
SELECT TYPE(this=>winv)
CLASS IS(oft_native_cg_solver)
    this%cleaner=>hcurl_divout
  CLASS DEFAULT
    CALL oft_abort('Error allocating winv solver', 'taylor_vacuum', __FILE__)
END SELECT
!---------------------------------------------------------------------------
! Create H(Curl) divergence cleaner
!---------------------------------------------------------------------------
CALL self%ML_lagrange%set_level(MIN(self%ML_hcurl%level,self%ML_lagrange%ml_mesh%mgdim))
CALL oft_lag_getlop(self%ML_lagrange%current_level,lop_lag,"zerob")
CALL create_cg_solver(linv_lag)
linv_lag%A=>lop_lag
!linv_lag%its=-2
CALL create_diag_pre(linv_lag%pre)
linv_lag%its=40
hcurl_divout%solver=>linv_lag
! hcurl_divout%bc=>lag_zerob
CALL hcurl_divout%setup(self%ML_hcurl,self%ML_lagrange,bc='zero',solver=linv_lag)
hcurl_divout%app_freq=2
CALL oft_hcurl_getmop(self%ML_hcurl%current_level,mop_hcurl,'zerob')
hcurl_divout%mop=>mop_hcurl
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
NULLIFY(tmp)
!---Get H(Curl) + Grad(H^1) mass matrix
CALL hcurl_grad_getmop(self%ML_hcurl_grad%current_level,mop,'none')
!---Allocate vacuum and current field containers
ALLOCATE(self%hvac(self%nh,self%ML_hcurl_grad%nlevels))
ALLOCATE(self%hcur(self%nh,self%ML_hcurl_grad%nlevels))
!---Create temporary H(Curl) vector
CALL self%ML_hcurl%vec_create(b)
!---Loop over cut planes
DO i=1,self%nh
!---------------------------------------------------------------------------
! Compute vacuum fields
!---------------------------------------------------------------------------
  !---Setup level fields
  CALL self%ML_hcurl_grad%vec_create(self%hvac(i,self%ML_hcurl_grad%level)%f)
  u=>self%hvac(i,self%ML_hcurl_grad%level)%f
  IF(.NOT.ASSOCIATED(tmp))CALL u%new(tmp)
  have_rst=.FALSE.
  IF(save_rst)THEN
    IF(oft_file_exist(rst_filename))THEN
      WRITE(field_name,'(A6,A2,A2,I2.2,A2,I2.2)')'hvac_g',self%ML_hcurl%ml_mesh%rlevel,'_p',self%ML_hcurl%current_level%order,'_h',i
      IF(hdf5_field_exist(rst_filename,field_name))THEN
        CALL self%ML_hcurl_grad%current_level%vec_load(u,rst_filename,field_name)
        have_rst=.TRUE.
      END IF
    END IF
  END IF
  !---Compute jump and solve
  IF(.NOT.have_rst)THEN
    CALL hcurl_grad_mc(self%ML_hcurl%ml_mesh%mesh,u,self%hcpc(:,i),self%hcpv(:,i),self%jtol)
    venergy=u%dot(u)
    IF(venergy<1.d-12)CALL oft_abort('Plane does not intersect mesh.','taylor_vacuum',__FILE__)
    divout%pm=oft_env%pm
    CALL divout%apply(u)
  END IF
!---------------------------------------------------------------------------
! Write restart file
!---------------------------------------------------------------------------
  IF(save_rst)THEN
    IF(.NOT.oft_file_exist(rst_filename))THEN
      WRITE(field_name,'(A6,A2,A2,I2.2,A2,I2.2)')'hvac_g',self%ML_hcurl%ml_mesh%rlevel,'_p',self%ML_hcurl%current_level%order,'_h',i
      CALL oft_mpi_barrier(ierr)
      CALL self%ML_hcurl_grad%current_level%vec_save(u,rst_filename,field_name)
    END IF
  END IF
  !---Compute field energy
  CALL mop%apply(u,tmp)
  venergy=u%dot(tmp)
  CALL u%scale(1.d0/venergy)
  IF(oft_env%head_proc)THEN
    WRITE(*,*)'Injector      = ',self%htag(i)
    WRITE(*,*)'Vacuum Energy = ',1.d0/venergy
  END IF
  IF(PRESENT(energy))energy(i)=1.d0/venergy
!---------------------------------------------------------------------------
! Compute current field
!---------------------------------------------------------------------------
  !---Use vacuum field as source term
  call tmp%scale(1.d0/venergy)
  call self%ML_hcurl%vec_create(self%hcur(i,self%ML_hcurl%level)%f)
  u=>self%hcur(i,self%ML_hcurl%level)%f
  have_rst=.FALSE.
  IF(save_rst)THEN
    IF(oft_file_exist(rst_filename))THEN
      WRITE(field_name,'(A6,A2,A2,I2.2,A2,I2.2)')'hcur_g',self%ML_hcurl%ml_mesh%rlevel,'_p',self%ML_hcurl%current_level%order,'_h',i
      IF(hdf5_field_exist(rst_filename,field_name))THEN
        CALL self%ML_hcurl%current_level%vec_load(u,rst_filename,field_name)
        have_rst=.TRUE.
      END IF
    END IF
  END IF
  IF(.NOT.have_rst)THEN
    !---Copy H(Curl) subpace into new vector
    CALL b%set(0.d0)
    NULLIFY(vals)
    CALL tmp%get_slice(vals,1)
    CALL b%restore_slice(vals)
    !---Compute current field
    hcurl_divout%pm=.FALSE.
    hcurl_divout%mop=>mop_hcurl
    CALL hcurl_zerob%apply(b)
    CALL winv%apply(u,b)
  END IF
  !---Clean divergence
  NULLIFY(hcurl_divout%mop)
  hcurl_divout%app_freq=1
  hcurl_divout%pm=.TRUE.
  CALL hcurl_divout%apply(u)
  !---
  IF(PRESENT(hmodes))THEN
    DO j=1,hmodes%nm
      CALL wop%apply(hmodes%hffa(j,self%ML_hcurl%level)%f,b)
      venergy = u%dot(b)
      IF(oft_env%head_proc)WRITE(*,'(A,I3,A,E10.3)')'Mode ',j,' Coupling = ',venergy
    END DO
  END IF
!---------------------------------------------------------------------------
! Write restart file
!---------------------------------------------------------------------------
  IF(save_rst)THEN
    IF(.NOT.oft_file_exist(rst_filename))THEN
      CALL oft_mpi_barrier(ierr)
      WRITE(field_name,'(A6,A2,A2,I2.2,A2,I2.2)')'hcur_g',self%ML_hcurl%ml_mesh%rlevel,'_p',self%ML_hcurl%current_level%order,'_h',i
      CALL self%ML_hcurl%current_level%vec_save(u,rst_filename,field_name)
    END IF
  END IF
END DO
!---
CALL tmp%delete
CALL b%delete
NULLIFY(tmp,b)
!---
CALL mop%delete
CALL mop_hcurl%delete
CALL lop_lag%delete
IF(ASSOCIATED(ml_wop))THEN
  DO i=1,SIZE(ml_wop)
    CALL ml_wop(i)%M%delete
    DEALLOCATE(ml_wop(i)%M)
  END DO
  DEALLOCATE(ml_wop)
END IF
IF(ASSOCIATED(ml_lop))THEN
  DO i=1,SIZE(ml_lop)
    CALL ml_lop(i)%M%delete
    DEALLOCATE(ml_lop(i)%M)
  END DO
  DEALLOCATE(ml_lop)
END IF
CALL linv%pre%delete
DEALLOCATE(linv%pre)
CALL linv%delete
DEALLOCATE(linv)
CALL winv%pre%delete
DEALLOCATE(winv%pre)
CALL winv%delete
DEALLOCATE(winv)
CALL linv_lag%pre%delete
DEALLOCATE(linv_lag%pre)
CALL linv_lag%delete
DEALLOCATE(linv_lag)
CALL self%ML_lagrange%set_level(self%ML_lagrange%nlevels)
DEBUG_STACK_POP
end subroutine taylor_vacuum
!---------------------------------------------------------------------------
!> Compute force-free plasma response to external fields generated by @ref
!! taylor::taylor_injectors "taylor_injectors"
!!
!! @note When `taylor_rst=.TRUE.` one restart files will be generated for
!! each jump plane `gffa_*.rst`
!---------------------------------------------------------------------------
subroutine taylor_injectors(self,hmodes,lambda,rst_filename)
type(oft_taylor_inhomo), intent(inout) :: self
type(oft_taylor_eigs), intent(inout) :: hmodes
real(r8), intent(in) :: lambda !< Desired lambda for force-free state
character(LEN=OFT_PATH_SLEN), optional, intent(in) :: rst_filename
CLASS(oft_matrix), POINTER :: lop_lag => NULL()
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_jmlb => NULL()
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_wop => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_matrix), POINTER :: kop => NULL()
CLASS(oft_matrix), POINTER :: wop => NULL()
CLASS(oft_matrix), POINTER :: jmlb_mat => NULL()
!---JMLB solver
TYPE(oft_native_gmres_solver), TARGET :: jmlb_inv
!---H(Curl) Divergence cleaner
CLASS(oft_solver), POINTER :: linv_lag => NULL()
TYPE(oft_hcurl_divout), TARGET :: hcurl_divout
!---
class(oft_vector), pointer :: u,b,tmp
real(r8) :: venergy,lam_file
integer(i4) :: i,k,ierr
logical :: do_orthog,have_rst,save_rst
character(LEN=16) :: field_name
type(oft_hcurl_zerob), target :: hcurl_zerob
type(oft_lag_zerob), target :: lag_zerob
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%hcpc))CALL oft_abort("Vacuum fields not available", &
"taylor_injectors", __FILE__)
!
save_rst=.FALSE.
IF(PRESENT(rst_filename))save_rst=.TRUE.
IF(save_rst)THEN
  have_rst=oft_file_exist(rst_filename)
  DO i=1,self%nh
    WRITE(field_name,'(A6,A2,A2,I2.2,A2,I2.2)')'gffa_g',self%ML_hcurl%ml_mesh%rlevel,'_p',self%ML_hcurl%current_level%order,'_h',i
    have_rst=have_rst.AND.hdf5_field_exist(rst_filename,field_name)
  END DO
ELSE
  have_rst=.FALSE.
END IF
hcurl_zerob%ML_hcurl_rep=>self%ML_hcurl
lag_zerob%ML_lag_rep=>self%ML_lagrange
NULLIFY(mop,kop,wop,jmlb_mat,ml_jmlb)
IF(.NOT.have_rst)THEN
  !---Orthogonalize (if within 5% of Taylor state)
  do_orthog=.FALSE.
  IF(hmodes%nm>0)THEN
    IF(ABS((lambda-hmodes%hlam(1,self%ML_hcurl%level))/hmodes%hlam(1,self%ML_hcurl%level))<5.d-2)THEN
      CALL oft_hcurl_getwop(self%ML_hcurl%current_level,wop,'zerob')
      hmodes%orthog%nm=1
      hmodes%orthog%wop=>wop
      do_orthog=.TRUE.
    END IF
  END IF
!---------------------------------------------------------------------------
! Setup H(Curl)::WOP preconditioner
!---------------------------------------------------------------------------
  IF(self%minlev==self%ML_hcurl%nlevels)THEN ! Lowest level uses diag precond
    CALL oft_hcurl_getjmlb(self%ML_hcurl%current_level,jmlb_mat,lambda,'zerob')
    CALL create_diag_pre(jmlb_inv%pre)
  ELSE ! Nested levels use MG
    CALL hcurl_getjmlb_pre(self%ML_hcurl,jmlb_inv%pre,ml_jmlb,lambda,nlevels=self%ML_hcurl%level-self%minlev+1)
    jmlb_mat=>ml_jmlb(self%ML_hcurl%level-self%minlev+1)%M
  END IF
!---------------------------------------------------------------------------
! Create H(Curl)::WOP solver
!---------------------------------------------------------------------------
  jmlb_inv%A=>jmlb_mat
  jmlb_inv%atol=1.d-8
  jmlb_inv%its=-3
  jmlb_inv%nrits=20
  jmlb_inv%itplot=1
  !jmlb_inv%bc=>hcurl_zerob
END IF
!---------------------------------------------------------------------------
! Create H(Curl) divergence cleaner
!---------------------------------------------------------------------------
CALL self%ML_lagrange%set_level(MIN(self%ML_hcurl%level,self%ML_lagrange%ml_mesh%mgdim))
CALL oft_lag_getlop(self%ML_lagrange%current_level,lop_lag,"zerob")
CALL create_cg_solver(linv_lag)
linv_lag%A=>lop_lag
linv_lag%its=-2
CALL create_diag_pre(linv_lag%pre)
linv_lag%its=40
hcurl_divout%solver=>linv_lag
CALL hcurl_divout%setup(self%ML_hcurl,self%ML_lagrange,bc='zero',solver=linv_lag)
! hcurl_divout%bc=>lag_zerob
hcurl_divout%app_freq=10
CALL oft_hcurl_getmop(self%ML_hcurl%current_level,mop,'zerob')
hcurl_divout%mop=>mop
!jmlb_inv%cleaner=>hcurl_divout
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
call self%ML_hcurl_grad%set_level(self%ML_hcurl_grad%nlevels,propogate=.TRUE.)
!---
call self%ML_hcurl%vec_create(tmp)
allocate(self%gffa(self%nh,self%ML_hcurl_grad%nlevels))
IF(.NOT.have_rst)CALL oft_hcurl_getkop(self%ML_hcurl%current_level,kop,'zerob')
!---
do i=1,self%nh
  !---
  CALL self%ML_hcurl%vec_create(self%gffa(i,self%ML_hcurl_grad%level)%f)
  u=>self%gffa(i,self%ML_hcurl_grad%level)%f
  have_rst=.FALSE.
  IF(save_rst)THEN
    IF(oft_file_exist(rst_filename))THEN
      lam_file=-1.d99
      IF(hdf5_field_exist(rst_filename,'lambda'))CALL hdf5_read(lam_file,rst_filename,'lambda')
      IF(ABS(lam_file-lambda)<1.d-5)THEN
        WRITE(field_name,'(A6,A2,A2,I2.2,A2,I2.2)')'gffa_g',self%ML_hcurl%ml_mesh%rlevel,'_p',self%ML_hcurl%current_level%order,'_h',i
        IF(hdf5_field_exist(rst_filename,field_name))THEN
          CALL self%ML_hcurl%current_level%vec_load(u,rst_filename,field_name)
          have_rst=.TRUE.
        END IF
      END IF
    END IF
  END IF
  IF(.NOT.have_rst)THEN
    CALL u%add(0.d0,1.d0,self%hcur(i,self%ML_hcurl_grad%level)%f)
    IF(do_orthog)CALL hmodes%orthog%apply(u)
    CALL kop%apply(u,tmp)
    !---Solve
    b=>tmp
    CALL hcurl_zerob%apply(b)
    CALL u%set(0.d0)
    CALL jmlb_inv%apply(u,b)
  END IF
  !---Clean divergence
  NULLIFY(hcurl_divout%mop)
  hcurl_divout%app_freq=1
  hcurl_divout%pm=.TRUE.
  CALL hcurl_divout%apply(u)
!---------------------------------------------------------------------------
! Write restart file
!---------------------------------------------------------------------------
  IF(save_rst)THEN
    IF(.NOT.oft_file_exist(rst_filename))THEN
      WRITE(field_name,'(A6,A2,A2,I2.2,A2,I2.2)')'gffa_g',self%ML_hcurl%ml_mesh%rlevel,'_p',self%ML_hcurl%current_level%order,'_h',i
      CALL oft_mpi_barrier(ierr)
      CALL self%ML_hcurl%current_level%vec_save(u,rst_filename,field_name)
      IF(oft_env%head_proc)CALL hdf5_write(lambda,rst_filename,'lambda')
    END IF
  END IF
  !---
  CALL u%add(lambda**2,lambda,self%hcur(i,self%ML_hcurl_grad%level)%f)
end do
!---
CALL tmp%delete
NULLIFY(tmp,b)
!---
CALL mop%delete
IF(ASSOCIATED(kop))THEN
  CALL kop%delete
  DEALLOCATE(kop)
END IF
CALL lop_lag%delete
IF(ASSOCIATED(ml_jmlb))THEN
  DO i=1,SIZE(ml_jmlb)
    CALL ml_jmlb(i)%M%delete
    DEALLOCATE(ml_jmlb(i)%M)
  END DO
  DEALLOCATE(ml_jmlb)
ELSE
  IF(ASSOCIATED(ml_jmlb))THEN
    CALL jmlb_mat%delete
    DEALLOCATE(jmlb_mat)
  END IF
  IF(ASSOCIATED(wop))THEN
    CALL wop%delete
    DEALLOCATE(wop)
  END IF
END IF
!---
IF(.NOT.have_rst)CALL jmlb_inv%pre%delete
CALL self%ML_lagrange%set_level(self%ML_lagrange%nlevels)
DEBUG_STACK_POP
end subroutine taylor_injectors
!---------------------------------------------------------------------------
!> Compute force-free plasma response to external fields generated by @ref
!! taylor::taylor_vacuum "taylor_vacuum"
!---------------------------------------------------------------------------
subroutine taylor_injector_single(self,hmodes,lambda,fluxes,gffa)
type(oft_taylor_inhomo), intent(inout) :: self
type(oft_taylor_eigs), intent(inout) :: hmodes
real(r8), intent(in) :: lambda !< Desired lambda for force-free state
real(r8), intent(in) :: fluxes(:) !< Flux for each injector
class(oft_vector), pointer, intent(inout) :: gffa !< Plasma component (non-vacuum) of injector field
CLASS(oft_matrix), POINTER :: lop_lag => NULL()
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_jmlb => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_matrix), POINTER :: kop => NULL()
CLASS(oft_matrix), POINTER :: wop => NULL()
CLASS(oft_matrix), POINTER :: jmlb_mat => NULL()
!---JMLB solver
TYPE(oft_native_gmres_solver), TARGET :: jmlb_inv
!---H(Curl) Divergence cleaner
CLASS(oft_solver), POINTER :: linv_lag => NULL()
TYPE(oft_hcurl_divout), TARGET :: hcurl_divout
!---
class(oft_vector), pointer :: tmp
real(r8) :: venergy,lam_file
integer(i4) :: i,k
character(2) :: pnum,mnum
character(40) :: filename
type(oft_hcurl_zerob), target :: hcurl_zerob
type(oft_lag_zerob), target :: lag_zerob
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%hcpc))CALL oft_abort("Vacuum fields not available", &
"taylor_injector_single", __FILE__)
NULLIFY(mop,kop,wop,ml_jmlb)
! IF(taylor_minlev<0)taylor_minlev=taylor_ML_hcurl%nlevels
!---------------------------------------------------------------------------
! Setup H(Curl)::WOP preconditioner
!---------------------------------------------------------------------------
IF(self%minlev==self%ML_hcurl%nlevels)THEN ! Lowest level uses diag precond
  CALL oft_hcurl_getjmlb(self%ML_hcurl%current_level,jmlb_mat,lambda,'zerob')
  CALL create_diag_pre(jmlb_inv%pre)
ELSE ! Nested levels use MG
  CALL hcurl_getjmlb_pre(self%ML_hcurl,jmlb_inv%pre,ml_jmlb,lambda,nlevels=1)
  jmlb_mat=>ml_jmlb(1)%M
END IF
hcurl_zerob%ML_hcurl_rep=>self%ML_hcurl
lag_zerob%ML_lag_rep=>self%ML_lagrange
!---------------------------------------------------------------------------
! Create H(Curl)::WOP solver
!---------------------------------------------------------------------------
jmlb_inv%A=>jmlb_mat
jmlb_inv%atol=1.d-7
jmlb_inv%its=-3
jmlb_inv%nrits=20
jmlb_inv%itplot=1
!jmlb_inv%bc=>hcurl_zerob
!---------------------------------------------------------------------------
! Create H(Curl) divergence cleaner
!---------------------------------------------------------------------------
CALL self%ML_lagrange%set_level(MIN(self%ML_hcurl%level,self%ML_lagrange%ml_mesh%mgdim))
CALL oft_lag_getlop(self%ML_lagrange%current_level,lop_lag,"zerob")
CALL create_cg_solver(linv_lag)
linv_lag%A=>lop_lag
linv_lag%its=-2
CALL create_diag_pre(linv_lag%pre)
linv_lag%its=40
hcurl_divout%solver=>linv_lag
CALL hcurl_divout%setup(self%ML_hcurl,self%ML_lagrange,bc='zero',solver=linv_lag)
! hcurl_divout%bc=>lag_zerob
hcurl_divout%app_freq=10
CALL oft_hcurl_getmop(self%ML_hcurl%current_level,mop,'zerob')
hcurl_divout%mop=>mop
!jmlb_inv%cleaner=>hcurl_divout
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
call self%ML_hcurl_grad%set_level(self%ML_hcurl_grad%nlevels,propogate=.TRUE.)
!---
call self%ML_hcurl%vec_create(gffa)
call self%ML_hcurl%vec_create(tmp)
!---
CALL gffa%set(0.d0)
do i=1,self%nh
  CALL gffa%add(1.d0,fluxes(i),self%hcur(i,self%ML_hcurl_grad%level)%f)
end do
!---Orthogonalize (if within 5% of Taylor state)
IF(hmodes%nm>0)THEN
  IF(ABS((lambda-hmodes%hlam(1,self%ML_hcurl%level))/hmodes%hlam(1,self%ML_hcurl%level))<5.d-2)THEN
    CALL oft_hcurl_getwop(self%ML_hcurl%current_level,wop,'zerob')
    ! orthog%orthog=>self%hffa
    hmodes%orthog%nm=1
    hmodes%orthog%wop=>wop
    CALL hmodes%orthog%apply(gffa)
    CALL wop%delete
    DEALLOCATE(wop)
  END IF
END IF
!---Get H(Curl) helicity matrix
CALL oft_hcurl_getkop(self%ML_hcurl%current_level,kop,'zerob')
CALL kop%apply(gffa,tmp)
CALL kop%delete
DEALLOCATE(kop)
!---Solve
CALL hcurl_zerob%apply(tmp)
CALL gffa%set(0.d0)
CALL jmlb_inv%apply(gffa,tmp)
!---Clean divergence
NULLIFY(hcurl_divout%mop)
hcurl_divout%app_freq=1
hcurl_divout%pm=.TRUE.
CALL hcurl_divout%apply(gffa)
!---
CALL gffa%scale(lambda**2)
do i=1,self%nh
  CALL gffa%add(1.d0,fluxes(i)*lambda,self%hcur(i,self%ML_hcurl_grad%level)%f)
end do
!---
CALL tmp%delete
DEALLOCATE(tmp)
!---
CALL lop_lag%delete
IF(ASSOCIATED(ml_jmlb))THEN
  DO i=1,SIZE(ml_jmlb)
    CALL ml_jmlb(i)%M%delete
    DEALLOCATE(ml_jmlb(i)%M)
  END DO
  DEALLOCATE(ml_jmlb)
ELSE
  CALL jmlb_mat%delete
  DEALLOCATE(jmlb_mat)
END IF
!---
CALL jmlb_inv%pre%delete
CALL self%ML_lagrange%set_level(self%ML_lagrange%nlevels)
DEBUG_STACK_POP
end subroutine taylor_injector_single
!---------------------------------------------------------------------------
!> Setup interpolator for composite Taylor state fields
!!
!! Fetches local representation used for interpolation from vector object
!!
!! @note Should only be used via class \ref oft_taylor_rinterp or children
!---------------------------------------------------------------------------
subroutine taylor_rinterp_setup1(self,hcurl_grad_rep)
class(oft_taylor_rinterp), intent(inout) :: self
class(oft_fem_comp_type), target, intent(inout) :: hcurl_grad_rep
DEBUG_STACK_PUSH
!---Get local slice
CALL self%ua%get_local(self%acurl)
CALL self%uvac%get_local(self%vac_curl,1)
CALL self%uvac%get_local(self%vac_grad,2)
SELECT TYPE(this=>hcurl_grad_rep%fields(1)%fe)
  CLASS IS(oft_hcurl_fem)
    self%curl_rep=>this
    self%mesh=>this%mesh
  CLASS DEFAULT
    CALL oft_abort("Invalid HCurl space","taylor_rinterp_setup1",__FILE__)
END SELECT
SELECT TYPE(this=>hcurl_grad_rep%fields(2)%fe)
  CLASS IS(oft_h1_fem)
    self%grad_rep=>this
  CLASS DEFAULT
    CALL oft_abort("Invalid HGrad space","taylor_rinterp_setup1",__FILE__)
END SELECT
DEBUG_STACK_POP
end subroutine taylor_rinterp_setup1
!---------------------------------------------------------------------------
!> Setup interpolator for composite Taylor state fields
!!
!! Fetches local representation used for interpolation from vector object
!!
!! @note Should only be used via class \ref oft_taylor_rinterp or children
!---------------------------------------------------------------------------
subroutine taylor_rinterp_setup2(self,hcurl_rep,hgrad_rep)
class(oft_taylor_rinterp), intent(inout) :: self
class(oft_afem_type), target, intent(inout) :: hcurl_rep
class(oft_afem_type), target, intent(inout) :: hgrad_rep
DEBUG_STACK_PUSH
!---Get local slice
CALL self%ua%get_local(self%acurl)
CALL self%uvac%get_local(self%vac_curl,1)
CALL self%uvac%get_local(self%vac_grad,2)
SELECT TYPE(hcurl_rep)
  CLASS IS(oft_hcurl_fem)
    self%curl_rep=>hcurl_rep
    self%mesh=>hcurl_rep%mesh
  CLASS DEFAULT
    CALL oft_abort("Invalid HCurl space","taylor_rinterp_setup2",__FILE__)
END SELECT
SELECT TYPE(hgrad_rep)
  CLASS IS(oft_h1_fem)
    self%grad_rep=>hgrad_rep
  CLASS DEFAULT
    CALL oft_abort("Invalid HGrad space","taylor_rinterp_setup2",__FILE__)
END SELECT
DEBUG_STACK_POP
end subroutine taylor_rinterp_setup2
!---------------------------------------------------------------------------
!> Reconstruct a composite Taylor state field
!---------------------------------------------------------------------------
subroutine taylor_rinterp(self,cell,f,gop,val)
class(oft_taylor_rinterp), intent(inout) :: self
integer(i4), intent(in) :: cell !< Cell for interpolation
real(r8), intent(in) :: f(:) !< Position in cell in logical coord [4]
real(r8), intent(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
real(r8), intent(out) :: val(:) !< Reconstructed field at f [3]
integer(i4), allocatable :: j(:)
integer(i4) :: jc
real(r8) :: cgop(3,6)
real(r8), allocatable :: rop(:,:)
DEBUG_STACK_PUSH
!---
IF(.NOT.ASSOCIATED(self%vac_grad))CALL oft_abort('Setup has not been called!','taylor_rinterp',__FILE__)
val=0.d0
!---Get curl dofs
allocate(j(self%curl_rep%nce),rop(3,self%curl_rep%nce))
call self%curl_rep%ncdofs(cell,j) ! get DOFs
CALL oft_hcurl_get_cgops(gop,cgop)
!---Reconstruct field
call oft_hcurl_eval_all(self%curl_rep,cell,f,rop,gop)
do jc=1,self%curl_rep%nce
  val=val+self%vac_curl(j(jc))*rop(:,jc)
end do
call oft_hcurl_ceval_all(self%curl_rep,cell,f,rop,cgop)
do jc=1,self%curl_rep%nce
  val=val+self%acurl(j(jc))*rop(:,jc)
end do
deallocate(j,rop)
!---Get curl dofs
allocate(j(self%grad_rep%nce),rop(3,self%grad_rep%nce))
call self%grad_rep%ncdofs(cell,j) ! get DOFs
!---Reconstruct field
call oft_h1_geval_all(self%grad_rep,cell,f,rop,gop)
do jc=1,self%grad_rep%nce
  val=val+self%vac_grad(j(jc))*rop(:,jc)
end do
deallocate(j,rop)
DEBUG_STACK_POP
end subroutine taylor_rinterp
end module taylor
