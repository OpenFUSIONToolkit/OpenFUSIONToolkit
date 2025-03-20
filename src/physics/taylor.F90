!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Interpolate a Taylor state field
!!
!! Taylor state fields consist of a gradient component, the vacuum field, and
!! a curl component, the plasma field.
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Force-free, uniform \f$ \lambda \f$ eigenmode object
!!
!! Used to compute/store force-free solutions to the eigenproblem \[$ \nabla \times \nabla \times A = \lambda \nabla \times A \]$.
!------------------------------------------------------------------------------
type :: oft_taylor_hmodes
  integer(i4) :: minlev = -1 !< Lowest FE level for multi-level solvers (`-1` indicates single level solve)
  integer(i4) :: nm = 0 !< Number of force-free fields computed (see @ref taylor_hmodes)
  integer(i4) :: htor_axis = 3 !< Index of coordinate to use as axis for toroidal flux calculation
  real(r8), pointer, dimension(:,:) :: hlam => NULL() !< Lambda values for each mode/level [nm,nlevels]
  real(r8), pointer, dimension(:,:) :: htor => NULL() !< Toroidal flux for each mode/level [nm,nlevels]
  type(oft_vector_ptr), pointer, dimension(:,:) :: hffa => NULL() !< Vector potential for each mode/level [nm,nlevels]
  type(oft_hcurl_orthog), pointer :: orthog => NULL() !< Orthogonalization operator
  TYPE(oft_ml_fem_type), POINTER :: ML_hcurl => NULL() !< Multi-level H(Curl) FE object
  TYPE(oft_ml_fem_type), POINTER :: ML_lagrange => NULL() !< Multi-level Lagrange FE object
CONTAINS
  !> Setup object before solution
  PROCEDURE :: setup => hmodes_setup
  !> Destory/reset object
  PROCEDURE :: delete => hmodes_delete
end type oft_taylor_hmodes
integer(i4), parameter :: taylor_tag_size = 4 !< Size of  jump planes character tags
!------------------------------------------------------------------------------
!> Force-free, uniform \f$ \lambda \f$ field object (inhomogeneous BCs)
!!
!! Used to compute/store fields of the form \[$ J_p = \lambda (B_p + B_v), \]$
!! where \f$ B_v \f$ is a vacuum field and \f$ B_p \f$ is the plasma response
!! to make the full field force-free for the given value of \f$ \lambda \f$.
!------------------------------------------------------------------------------
type :: oft_taylor_ifield
  integer(i4) :: minlev = -1 !< Lowest FE level for MG solvers (`-1` indicates single level solve)
  integer(i4) :: nh = -1 !< Number of handles (jump planes) in current geometry
  real(r8) :: jtol = 1.d-6 !< Tolerance for identifying edges on jump plane
  real(r8) :: lambda = 0.d0 !< Lambda value for solutions in `gffa`
  real(r8), pointer, dimension(:,:)  :: hcpc => NULL() !< Center points of handles (jump planes)
  real(r8), pointer, dimension(:,:)  :: hcpv => NULL() !< Normal vectors for handles (jump planes)
  character(LEN=taylor_tag_size), pointer, dimension(:) :: htag => NULL() !< Handle names
  type(oft_vector_ptr), pointer, dimension(:) :: hvac => NULL() !< Vacuum magnetic fields
  type(oft_vector_ptr), pointer, dimension(:) :: hcur => NULL() !< Inhomogeneous source fields
  type(oft_vector_ptr), pointer, dimension(:) :: gffa => NULL() !< Vector potential for plasma fields at specified \f$ \lambda \f$
  TYPE(oft_ml_fem_type), POINTER :: ML_lagrange => NULL() !< Multi-level Lagrange FE object
  TYPE(oft_ml_fem_type), POINTER :: ML_h1 => NULL() !< Multi-level H^1 FE object
  TYPE(oft_ml_fem_type), POINTER :: ML_h1grad => NULL() !< Multi-level grad(H^1) FE object
  TYPE(oft_ml_fem_type), POINTER :: ML_hcurl => NULL() !< Multi-level H(Curl) FE object
  TYPE(oft_ml_fem_comp_type), POINTER :: ML_hcurl_grad => NULL() !< Multi-level H(Curl) + grad(H^1) FE object
CONTAINS
  !> Setup object before solution
  PROCEDURE :: setup => ff_setup
  !> Setup object before solution
  PROCEDURE :: delete => ff_delete
end type oft_taylor_ifield
contains
!------------------------------------------------------------------------------
!> Setup eigenmodes object
!------------------------------------------------------------------------------
subroutine hmodes_setup(self,ML_hcurl,ML_lagrange,minlev,htor_axis)
CLASS(oft_taylor_hmodes), INTENT(inout) :: self !< Force-free eigenmode object
TYPE(oft_ml_fem_type), OPTIONAL, TARGET, INTENT(in) :: ML_hcurl !< Multi-level H(Curl) FE representation
TYPE(oft_ml_fem_type), OPTIONAL, TARGET, INTENT(in) :: ML_lagrange !< Multi-level Lagrange FE representation
INTEGER(i4), OPTIONAL, INTENT(in) :: minlev
INTEGER(i4), OPTIONAL, INTENT(in) :: htor_axis
IF(PRESENT(ML_hcurl))THEN
  self%ML_hcurl=>ML_hcurl
ELSE
  IF(.NOT.ASSOCIATED(self%ML_hcurl))CALL oft_abort("No H(Curl) FE representation","hmodes_setup",__FILE__)
END IF
IF(PRESENT(ML_lagrange))THEN
  self%ML_lagrange=>ML_lagrange
ELSE
  IF(.NOT.ASSOCIATED(self%ML_lagrange))CALL oft_abort("No Lagrange FE representation","hmodes_setup",__FILE__)
END IF
!
IF(PRESENT(minlev))self%minlev=minlev
IF(self%minlev<0)self%minlev=self%ML_hcurl%level
IF(PRESENT(htor_axis))THEN
  IF((htor_axis<1).OR.(htor_axis>3))CALL oft_abort("Invalid value for 'htor_axis'","hmodes_setup",__FILE__)
  self%htor_axis=htor_axis
END IF
end subroutine hmodes_setup
!------------------------------------------------------------------------------
!> Setup eigenmodes object
!------------------------------------------------------------------------------
subroutine hmodes_delete(self,storage_only)
CLASS(oft_taylor_hmodes), INTENT(inout) :: self !< Force-free eigenmode object
LOGICAL, OPTIONAL, INTENT(in) :: storage_only !< Only reset storage, but do not clear references
LOGICAL :: do_nullify
INTEGER(i4) :: i,j
do_nullify=.TRUE.
IF(PRESENT(storage_only))do_nullify=.NOT.storage_only
!---Deallocate fields
IF(ASSOCIATED(self%hffa))THEN
  DO i=1,SIZE(self%hffa,DIM=1)
    DO j=1,SIZE(self%hffa,DIM=2)
      IF(ASSOCIATED(self%hffa(i,j)%f))THEN
        CALL self%hffa(i,j)%f%delete()
        DEALLOCATE(self%hffa(i,j)%f)
      END IF
    END DO
  END DO
  DEALLOCATE(self%hffa)
END IF
IF(ASSOCIATED(self%hlam))DEALLOCATE(self%hlam)
IF(ASSOCIATED(self%htor))DEALLOCATE(self%htor)
IF(ASSOCIATED(self%orthog))THEN
  CALL self%orthog%delete()
  DEALLOCATE(self%orthog)
END IF
!---Nullify pointers and reset defaults
IF(do_nullify)THEN
  NULLIFY(self%ML_hcurl,self%ML_lagrange)
  self%minlev=-1
  self%htor_axis=3
END IF
end subroutine hmodes_delete
!------------------------------------------------------------------------------
!> Setup eigenmodes object
!------------------------------------------------------------------------------
subroutine ff_setup(self,nh,hcpc,hcpv,htags,ML_hcurl,ML_h1,ML_hcurl_grad,ML_h1grad,ML_lagrange,minlev)
CLASS(oft_taylor_ifield), INTENT(inout) :: self !< Force-free field object
INTEGER(i4), INTENT(in) :: nh !< Number of jump planes
REAL(r8), INTENT(in) :: hcpc(3,nh) !< Jump plane center possitions [3,nh]
REAL(r8), INTENT(in) :: hcpv(3,nh) !< Jump plane normal vectors [3,nh]
CHARACTER(LEN=taylor_tag_size), OPTIONAL, INTENT(in) :: htags(nh) !< Names for each jump plane [LEN=taylor_tag_size,nh] (optional)
TYPE(oft_ml_fem_type), OPTIONAL, TARGET, INTENT(in) :: ML_hcurl !< Multi-level H(Curl) FE representation
TYPE(oft_ml_fem_type), OPTIONAL, TARGET, INTENT(in) :: ML_h1 !< Multi-level H^1 FE representation
TYPE(oft_ml_fem_comp_type), OPTIONAL, TARGET, INTENT(in) :: ML_hcurl_grad !< Multi-level H(Curl) + grad(H^1) FE representation
TYPE(oft_ml_fem_type), OPTIONAL, TARGET, INTENT(in) :: ML_h1grad !< Multi-level grad(H^1) FE representation
TYPE(oft_ml_fem_type), OPTIONAL, TARGET, INTENT(in) :: ML_lagrange !< Multi-level Lagrange FE representation
INTEGER(i4), OPTIONAL, INTENT(in) :: minlev
integer(i4) :: i
IF(PRESENT(ML_hcurl))THEN
  self%ML_hcurl=>ML_hcurl
ELSE
  IF(.NOT.ASSOCIATED(self%ML_hcurl))CALL oft_abort("No H(Curl) FE representation","ff_setup",__FILE__)
END IF
IF(PRESENT(ML_h1))THEN
  self%ML_h1=>ML_h1
ELSE
  IF(.NOT.ASSOCIATED(self%ML_h1))CALL oft_abort("No H^1 FE representation","ff_setup",__FILE__)
END IF
IF(PRESENT(ML_hcurl_grad))THEN
  self%ML_hcurl_grad=>ML_hcurl_grad
ELSE
  IF(.NOT.ASSOCIATED(self%ML_hcurl_grad))CALL oft_abort("No H(Curl) + grad(H^1) FE representation","ff_setup",__FILE__)
END IF
IF(PRESENT(ML_h1grad))THEN
  self%ML_h1grad=>ML_h1grad
ELSE
  IF(.NOT.ASSOCIATED(self%ML_h1grad))CALL oft_abort("No grad(H^1) FE representation","ff_setup",__FILE__)
END IF
IF(PRESENT(ML_lagrange))THEN
  self%ML_lagrange=>ML_lagrange
ELSE
  IF(.NOT.ASSOCIATED(self%ML_lagrange))CALL oft_abort("No Lagrange FE representation","ff_setup",__FILE__)
END IF
!
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
!
IF(PRESENT(minlev))self%minlev=minlev
IF(self%minlev<0)self%minlev=self%ML_hcurl%level
end subroutine ff_setup
!------------------------------------------------------------------------------
!> Setup eigenmodes object
!------------------------------------------------------------------------------
subroutine ff_delete(self,storage_only)
CLASS(oft_taylor_ifield), INTENT(inout) :: self !< Force-free eigenmode object
LOGICAL, OPTIONAL, INTENT(in) :: storage_only !< Only reset storage, but do not clear references
LOGICAL :: do_nullify
INTEGER(i4) :: i,j
do_nullify=.TRUE.
IF(PRESENT(storage_only))do_nullify=.NOT.storage_only
!---Deallocate fields
IF(ASSOCIATED(self%hvac))THEN
  DO i=1,SIZE(self%hvac)
    IF(ASSOCIATED(self%hvac(i)%f))THEN
      CALL self%hvac(i)%f%delete()
      DEALLOCATE(self%hvac(i)%f)
    END IF
  END DO
  DEALLOCATE(self%hvac)
END IF
IF(ASSOCIATED(self%hcur))THEN
  DO i=1,SIZE(self%hcur)
    IF(ASSOCIATED(self%hcur(i)%f))THEN
      CALL self%hcur(i)%f%delete()
      DEALLOCATE(self%hcur(i)%f)
    END IF
  END DO
  DEALLOCATE(self%hcur)
END IF
IF(ASSOCIATED(self%gffa))THEN
  DO i=1,SIZE(self%gffa)
    IF(ASSOCIATED(self%gffa(i)%f))THEN
      CALL self%gffa(i)%f%delete()
      DEALLOCATE(self%gffa(i)%f)
    END IF
  END DO
  DEALLOCATE(self%gffa)
END IF
IF(ASSOCIATED(self%hcpc))DEALLOCATE(self%hcpc)
IF(ASSOCIATED(self%hcpv))DEALLOCATE(self%hcpv)
IF(ASSOCIATED(self%htag))DEALLOCATE(self%htag)
!---Nullify pointers and reset defaults
IF(do_nullify)THEN
  NULLIFY(self%ML_hcurl,self%ML_lagrange)
  NULLIFY(self%ML_h1,self%ML_hcurl_grad,self%ML_h1grad)
  self%minlev=-1
  self%nh=0
END IF
end subroutine ff_delete
!------------------------------------------------------------------------------
!> Compute 'taylor_nm' Force-Free eignemodes.
!------------------------------------------------------------------------------
subroutine taylor_hmodes(self,nm,rst_filename)
type(oft_taylor_hmodes), intent(inout) :: self !< Force-free eigenmode object
integer(i4), optional, intent(in) :: nm !< Number of modes to compute (optional: 1)
character(LEN=*), optional, intent(in) :: rst_filename !< File name to store/load restart information
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
logical :: save_rst,rst_append
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
rst_append=.FALSE.
IF(PRESENT(rst_filename))THEN
  IF(TRIM(rst_filename)/='')save_rst=.TRUE.
END IF
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(A)')'============================'
  WRITE(*,'(A)')'Starting calculation of Taylor states'
  WRITE(*,'(A)')'============================'
  WRITE(*,*)
  CALL mytimer%tick
END IF
!---Deallocate storage if exists
CALL self%delete(storage_only=.TRUE.)
!---Allocate storage
ALLOCATE(self%hffa(self%nm,self%ML_hcurl%level))
ALLOCATE(self%hlam(self%nm,self%ML_hcurl%level))
ALLOCATE(self%htor(self%nm,self%ML_hcurl%level))
!---Setup orthogonalization
ALLOCATE(self%orthog)
self%orthog%ML_hcurl_rep=>self%ML_hcurl
self%orthog%orthog=>self%hffa
!------------------------------------------------------------------------------
! Create ML Matrices
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
! Loop over desired number of modes
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
! Create eigenvalue solver
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
! Setup preconditioner
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
! Initialize guess
!------------------------------------------------------------------------------
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
        IF(hdf5_field_exist(rst_filename,field_name))THEN
          CALL self%ML_hcurl%current_level%vec_load(u,rst_filename,field_name,ierr)
          IF(ierr/=0)WRITE(*,*)'Error reading field, skipping'
        END IF
      end if
    end if
!------------------------------------------------------------------------------
! Solve
!------------------------------------------------------------------------------
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
    self%htor(k,i) = tfluxfun(Bfield%mesh,Bfield,self%ML_hcurl%current_level%quad%order,self%htor_axis)
    CALL Bfield%delete
!------------------------------------------------------------------------------
! Create divergence cleaner
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
! Write restart files
!------------------------------------------------------------------------------
    if(save_rst)then
      CALL oft_mpi_barrier(ierr)
      WRITE(field_name,'(A6,A2,A2,I2.2,A2,I2.2)')'hffa_g',self%ML_hcurl%ml_mesh%rlevel,'_p',self%ML_hcurl%current_level%order,'_m',k
      CALL self%ML_hcurl%current_level%vec_save(u,rst_filename,field_name,append=rst_append)
      rst_append=.TRUE.
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
!------------------------------------------------------------------------------
! Deallocate ML matrices
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Generate vacuum fields for specified handle (cut planes)
!------------------------------------------------------------------------------
subroutine taylor_vacuum(self,energy,hmodes,rst_filename)
type(oft_taylor_ifield), intent(inout) :: self !< Force-free field object
real(r8), optional, intent(out) :: energy(:) !< Vacuum energy for each handle (optional)
type(oft_taylor_hmodes), optional, intent(inout) :: hmodes !< Force-free eigenmode object (optional)
character(LEN=*), optional, intent(in) :: rst_filename !< File name to store/load restart information (optional)
!---H(Curl) full space divergence cleaner
CLASS(oft_solver), POINTER :: linv => NULL()
TYPE(oft_hcurl_grad_divout) :: divout
!--- ML structures for MG-preconditioner
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_lop => NULL()
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: levels,nu
REAL(r8), ALLOCATABLE, DIMENSION(:) :: df
!---
class(oft_matrix), pointer :: lop => NULL()
class(oft_matrix), pointer :: mop => NULL()
!---
class(oft_vector), pointer :: u,b,tmp
real(r8), pointer, dimension(:) :: vals => NULL()
real(r8) :: alam,venergy
integer(i4) :: i,j,k,nlevels,ierr
logical :: have_rst,save_rst
character(LEN=16) :: field_name
DEBUG_STACK_PUSH
NULLIFY(lop,mop)
NULLIFY(ml_lop)
!---Loop over cut planes
save_rst=.FALSE.
IF(PRESENT(rst_filename))save_rst=.TRUE.
IF(save_rst)THEN
  have_rst=oft_file_exist(TRIM(rst_filename))
  DO i=1,self%nh
    WRITE(field_name,'(A6,A2,A2,I2.2,A2,I2.2)')'hvac_g',self%ML_hcurl%ml_mesh%rlevel,'_p',self%ML_hcurl%current_level%order,'_h',i
    have_rst=have_rst.AND.hdf5_field_exist(rst_filename,field_name)
  END DO
ELSE
  have_rst=.FALSE.
END IF
IF(.NOT.have_rst)THEN
!------------------------------------------------------------------------------
! Setup H^1::LOP preconditioner
!------------------------------------------------------------------------------
  if(self%minlev==self%ML_h1%nlevels-1)then ! Lowest level uses diag precond
    CALL oft_h1_getlop(self%ML_h1%current_level,lop,'grnd')
    CALL create_cg_solver(linv)
    CALL create_diag_pre(linv%pre)
  else ! Nested levels use MG
    CALL create_cg_solver(linv, force_native=.TRUE.)
    CALL h1_getlop_pre(self%ML_h1,linv%pre,ml_lop,'grnd',nlevels=self%ML_h1%nlevels-self%minlev+1)
      lop=>ml_lop(self%ML_h1%nlevels-self%minlev+1)%M
  end if
!------------------------------------------------------------------------------
! Create divergence cleaner
!------------------------------------------------------------------------------
  linv%A=>lop
  linv%its=-2
  ! divout%solver=>linv
  CALL divout%setup(self%ML_hcurl_grad,'grnd',solver=linv)
ELSE
  CALL create_cg_solver(linv)
  CALL create_diag_pre(linv%pre)
END IF
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
NULLIFY(tmp)
!---Get H(Curl) + Grad(H^1) mass matrix
CALL hcurl_grad_getmop(self%ML_hcurl_grad%current_level,mop,'none')
!---Allocate vacuum and current field containers
ALLOCATE(self%hvac(self%nh))
!---Create temporary H(Curl) vector
CALL self%ML_hcurl%vec_create(b)
!---Loop over cut planes
DO i=1,self%nh
!------------------------------------------------------------------------------
! Compute vacuum fields
!------------------------------------------------------------------------------
  !---Setup level fields
  CALL self%ML_hcurl_grad%vec_create(self%hvac(i)%f)
  u=>self%hvac(i)%f
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
!------------------------------------------------------------------------------
! Write restart file
!------------------------------------------------------------------------------
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
END DO
!---
CALL tmp%delete
CALL b%delete
NULLIFY(tmp,b)
!---
CALL mop%delete
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
CALL self%ML_lagrange%set_level(self%ML_lagrange%nlevels)
DEBUG_STACK_POP
end subroutine taylor_vacuum
!------------------------------------------------------------------------------
!> Generate vector potential whose corresponding current matches the
!! vacuum fields stored in @ref taylor::taylor_hvac
!------------------------------------------------------------------------------
subroutine taylor_vac_curr(self,hmodes,rst_filename)
type(oft_taylor_ifield), intent(inout) :: self !< Force-free field object
type(oft_taylor_hmodes), intent(inout) :: hmodes !< Force-free eigenmode object
character(LEN=*), optional, intent(in) :: rst_filename !< File name to store/load restart information (optional)
!---WOP solver
CLASS(oft_solver), POINTER :: winv => NULL()
!---H(Curl) subspace divergence cleaner
CLASS(oft_solver), POINTER :: linv_lag => NULL()
TYPE(oft_hcurl_divout), TARGET :: hcurl_divout
!--- ML structures for MG-preconditioner
TYPE(oft_matrix_ptr), POINTER, DIMENSION(:) :: ml_wop => NULL()
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: levels,nu
REAL(r8), ALLOCATABLE, DIMENSION(:) :: df
!---
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
NULLIFY(lop_lag,mop,mop_hcurl,wop)
NULLIFY(ml_wop)
!---Create taylor module variables
hcurl_zerob%ML_hcurl_rep=>self%ML_hcurl
lag_zerob%ML_lag_rep=>self%ML_lagrange
!---Loop over cut planes
save_rst=.FALSE.
IF(PRESENT(rst_filename))save_rst=.TRUE.
IF(save_rst)THEN
  have_rst=oft_file_exist(TRIM(rst_filename))
  DO i=1,self%nh
    WRITE(field_name,'(A6,A2,A2,I2.2,A2,I2.2)')'hcur_g',self%ML_hcurl%ml_mesh%rlevel,'_p',self%ML_hcurl%current_level%order,'_h',i
    have_rst=have_rst.AND.hdf5_field_exist(rst_filename,field_name)
  END DO
ELSE
  have_rst=.FALSE.
END IF
!------------------------------------------------------------------------------
! Setup H(Curl)::WOP preconditioner
!------------------------------------------------------------------------------
CALL create_cg_solver(winv, force_native=.TRUE.)
IF((self%minlev==self%ML_hcurl%level).OR.have_rst)THEN ! Lowest level uses diag precond
  NULLIFY(wop)
  CALL oft_hcurl_getwop(self%ML_hcurl%current_level,wop,'zerob')
  CALL create_diag_pre(winv%pre)
ELSE ! Nested levels use MG
  CALL hcurl_getwop_pre(self%ML_hcurl,winv%pre,ml_wop,nlevels=self%ML_hcurl%level-self%minlev+1)
  wop=>ml_wop(self%ML_hcurl%level-self%minlev+1)%M
END IF
!------------------------------------------------------------------------------
! Create H(Curl)::WOP solver
!------------------------------------------------------------------------------
winv%A=>wop
winv%its=-3
winv%atol=1.d-9
SELECT TYPE(this=>winv)
CLASS IS(oft_native_cg_solver)
    this%cleaner=>hcurl_divout
  CLASS DEFAULT
    CALL oft_abort('Error allocating winv solver', 'taylor_vac_curr', __FILE__)
END SELECT
!------------------------------------------------------------------------------
! Create H(Curl) divergence cleaner
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
NULLIFY(tmp)
!---Get H(Curl) + Grad(H^1) mass matrix
CALL hcurl_grad_getmop(self%ML_hcurl_grad%current_level,mop,'none')
!---Allocate vacuum and current field containers
ALLOCATE(self%hcur(self%nh))
!---Create temporary H(Curl) vector
CALL self%ML_hcurl%vec_create(b)
!---Loop over cut planes
DO i=1,self%nh
!------------------------------------------------------------------------------
! Compute current field
!------------------------------------------------------------------------------
  !---Use vacuum field as source term
  u=>self%hvac(i)%f
  IF(.NOT.ASSOCIATED(tmp))CALL u%new(tmp)
  CALL mop%apply(u,tmp)
  call self%ML_hcurl%vec_create(self%hcur(i)%f)
  u=>self%hcur(i)%f
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
  DO j=1,hmodes%nm
    CALL wop%apply(hmodes%hffa(j,self%ML_hcurl%level)%f,b)
    venergy = u%dot(b)
    IF(oft_env%head_proc)WRITE(*,'(A,I3,A,E10.3)')'Mode ',j,' Coupling = ',venergy
  END DO
!------------------------------------------------------------------------------
! Write restart file
!------------------------------------------------------------------------------
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
end subroutine taylor_vac_curr
!------------------------------------------------------------------------------
!> Compute force-free plasma response to external fields generated by @ref
!! taylor::taylor_injectors "taylor_injectors"
!------------------------------------------------------------------------------
subroutine taylor_injectors(self,hmodes,lambda,rst_filename)
type(oft_taylor_ifield), intent(inout) :: self !< Force-free field object
type(oft_taylor_hmodes), intent(inout) :: hmodes !< Force-free eigenmode object
real(r8), intent(in) :: lambda !< Desired lambda for force-free state
character(LEN=*), optional, intent(in) :: rst_filename !< File name to store/load restart information (optional)
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
IF(.NOT.ASSOCIATED(self%hcur))CALL taylor_vac_curr(self,hmodes,rst_filename)
WRITE(*,*)'Inj curr done'
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
!------------------------------------------------------------------------------
! Setup H(Curl)::WOP preconditioner
!------------------------------------------------------------------------------
  IF(self%minlev==self%ML_hcurl%nlevels)THEN ! Lowest level uses diag precond
    CALL oft_hcurl_getjmlb(self%ML_hcurl%current_level,jmlb_mat,lambda,'zerob')
    CALL create_diag_pre(jmlb_inv%pre)
  ELSE ! Nested levels use MG
    CALL hcurl_getjmlb_pre(self%ML_hcurl,jmlb_inv%pre,ml_jmlb,lambda,nlevels=self%ML_hcurl%level-self%minlev+1)
    jmlb_mat=>ml_jmlb(self%ML_hcurl%level-self%minlev+1)%M
  END IF
!------------------------------------------------------------------------------
! Create H(Curl)::WOP solver
!------------------------------------------------------------------------------
  jmlb_inv%A=>jmlb_mat
  jmlb_inv%atol=1.d-8
  jmlb_inv%its=-3
  jmlb_inv%nrits=20
  jmlb_inv%itplot=1
  !jmlb_inv%bc=>hcurl_zerob
END IF
!------------------------------------------------------------------------------
! Create H(Curl) divergence cleaner
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
call self%ML_hcurl_grad%set_level(self%ML_hcurl_grad%nlevels,propogate=.TRUE.)
!---
call self%ML_hcurl%vec_create(tmp)
allocate(self%gffa(self%nh))
IF(.NOT.have_rst)CALL oft_hcurl_getkop(self%ML_hcurl%current_level,kop,'zerob')
!---
do i=1,self%nh
  !---
  CALL self%ML_hcurl%vec_create(self%gffa(i)%f)
  u=>self%gffa(i)%f
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
    CALL u%add(0.d0,1.d0,self%hcur(i)%f)
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
!------------------------------------------------------------------------------
! Write restart file
!------------------------------------------------------------------------------
  IF(save_rst)THEN
    IF(.NOT.oft_file_exist(rst_filename))THEN
      WRITE(field_name,'(A6,A2,A2,I2.2,A2,I2.2)')'gffa_g',self%ML_hcurl%ml_mesh%rlevel,'_p',self%ML_hcurl%current_level%order,'_h',i
      CALL oft_mpi_barrier(ierr)
      CALL self%ML_hcurl%current_level%vec_save(u,rst_filename,field_name)
      IF(oft_env%head_proc)CALL hdf5_write(lambda,rst_filename,'lambda')
    END IF
  END IF
  !---
  CALL u%add(lambda**2,lambda,self%hcur(i)%f)
end do
self%lambda=lambda
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
!------------------------------------------------------------------------------
!> Compute force-free plasma response to external fields generated by @ref
!! taylor::taylor_vacuum "taylor_vacuum"
!------------------------------------------------------------------------------
subroutine taylor_injector_single(self,hmodes,lambda,fluxes,gffa)
type(oft_taylor_ifield), intent(inout) :: self !< Force-free field object
type(oft_taylor_hmodes), intent(inout) :: hmodes !< Force-free eigenmode object
real(r8), intent(in) :: lambda !< Desired lambda for force-free state
real(r8), intent(in) :: fluxes(:) !< Flux for each handle
class(oft_vector), pointer, intent(inout) :: gffa !< Plasma component (non-vacuum) of handle field
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
!------------------------------------------------------------------------------
! Setup H(Curl)::WOP preconditioner
!------------------------------------------------------------------------------
IF(self%minlev==self%ML_hcurl%nlevels)THEN ! Lowest level uses diag precond
  CALL oft_hcurl_getjmlb(self%ML_hcurl%current_level,jmlb_mat,lambda,'zerob')
  CALL create_diag_pre(jmlb_inv%pre)
ELSE ! Nested levels use MG
  CALL hcurl_getjmlb_pre(self%ML_hcurl,jmlb_inv%pre,ml_jmlb,lambda,nlevels=1)
  jmlb_mat=>ml_jmlb(1)%M
END IF
hcurl_zerob%ML_hcurl_rep=>self%ML_hcurl
lag_zerob%ML_lag_rep=>self%ML_lagrange
!------------------------------------------------------------------------------
! Create H(Curl)::WOP solver
!------------------------------------------------------------------------------
jmlb_inv%A=>jmlb_mat
jmlb_inv%atol=1.d-7
jmlb_inv%its=-3
jmlb_inv%nrits=20
jmlb_inv%itplot=1
!jmlb_inv%bc=>hcurl_zerob
!------------------------------------------------------------------------------
! Create H(Curl) divergence cleaner
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
call self%ML_hcurl_grad%set_level(self%ML_hcurl_grad%nlevels,propogate=.TRUE.)
!---
call self%ML_hcurl%vec_create(gffa)
call self%ML_hcurl%vec_create(tmp)
!---
CALL gffa%set(0.d0)
do i=1,self%nh
  CALL gffa%add(1.d0,fluxes(i),self%hcur(i)%f)
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
  CALL gffa%add(1.d0,fluxes(i)*lambda,self%hcur(i)%f)
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
!------------------------------------------------------------------------------
!> Setup interpolator for composite Taylor state fields
!!
!! Fetches local representation used for interpolation from vector object
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Setup interpolator for composite Taylor state fields
!!
!! Fetches local representation used for interpolation from vector object
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Reconstruct a composite Taylor state field
!------------------------------------------------------------------------------
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
