!---------------------------------------------------------------------------
!> @file marklin_f.F90
!
!> Fortran part of Python wrapper for Marklin force-free ideal MHD equilibrium functionality
!!
!! @authors Chris Hansen
!! @date May 2023
!! @ingroup doxy_oft_python
!---------------------------------------------------------------------------
MODULE marklin_f
USE iso_c_binding, ONLY: c_int, c_double, c_char, c_loc, c_null_char, c_ptr, &
    c_f_pointer, c_bool, c_null_ptr
!---Base
USE oft_base
USE oft_io, ONLY: hdf5_create_file
!--Grid
USE oft_mesh_type, ONLY: mesh, mesh_findcell
USE oft_mesh_native, ONLY: r_mem, lc_mem, reg_mem
USE multigrid_build, ONLY: multigrid_construct
!---Linear Algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---
USE fem_utils, ONLY: fem_interp
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup, oft_lagrange_nlevels
USE oft_lag_fields, ONLY: oft_lag_vcreate
USE oft_lag_operators, ONLY: lag_lop_eigs, lag_setup_interp, lag_mloptions, &
    oft_lag_vgetmop, oft_lag_vproject
!---H1(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl, oft_hcurl_setup, oft_hcurl_level, &
  oft_hcurl_minlev
USE oft_hcurl_operators, ONLY: oft_hcurl_cinterp, hcurl_setup_interp, &
    hcurl_mloptions
!---H1(Grad) FE space
USE oft_h0_basis, ONLY: oft_h0_setup
USE oft_h0_operators, ONLY: h0_setup_interp, oft_h0_getlop, h0_zerogrnd
!---H1 FE space
USE oft_h1_basis, ONLY: oft_h1_setup, oft_h1_nlevels
USE oft_h1_fields, ONLY: oft_h1_create
USE oft_h1_operators, ONLY: oft_h1_divout, h1_zeroi, h1_mc, h1curl_zerob, &
  h1_setup_interp, oft_h1_rinterp
!---Taylor state
USE taylor, ONLY: taylor_minlev, taylor_hmodes, taylor_hffa, taylor_nm, &
  taylor_rst
USE mhd_utils, ONLY: mu0
!---Wrappers
USE oft_base_f, ONLY: copy_string, copy_string_rev
IMPLICIT NONE
#include "local.h"
!
integer(i4), POINTER :: lc_plot(:,:) !< Needs docs
integer(i4), POINTER :: reg_plot(:) !< Needs docs
real(r8), POINTER :: r_plot(:,:) !< Needs docs
CONTAINS
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_setup(order,nmodes,minlev,save_rst,error_str) BIND(C,NAME="marklin_setup")
INTEGER(KIND=c_int), VALUE, INTENT(in) :: order !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: nmodes !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: minlev !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: save_rst !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(80) !< Needs docs
!---Lagrange mass solver
CLASS(oft_matrix), POINTER :: lmop => NULL()
CLASS(oft_solver), POINTER :: lminv => NULL()
!---Local variables
INTEGER(i4) :: i,io_unit,ierr
REAL(r8), POINTER, DIMENSION(:) :: vals => NULL()
REAL(r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: bvout
CLASS(oft_vector), POINTER :: u,v,check
TYPE(oft_hcurl_cinterp) :: Bfield
CHARACTER(LEN=3) :: pltnum
!---Clear error flag
CALL copy_string('',error_str)
!---Setup grid
! CALL mesh%setup_io(order)
! smesh%tess_order=order
!---Lagrange
CALL oft_lag_setup(order,minlev)
CALL lag_setup_interp
CALL lag_mloptions
!---H1(Curl) subspace
CALL oft_hcurl_setup(order,minlev)
CALL hcurl_setup_interp
CALL hcurl_mloptions
!---Compute modes
IF(minlev<0)THEN
  taylor_minlev=oft_hcurl_level
ELSE
  taylor_minlev=minlev
  IF(oft_env%nprocs>1)taylor_minlev=MAX(oft_env%nbase+1,minlev)
END IF
oft_env%pm=.TRUE.
taylor_rst=save_rst
CALL taylor_hmodes(nmodes)
END SUBROUTINE marklin_setup
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_save_visit(error_str) BIND(C,NAME="marklin_save_visit")
CHARACTER(KIND=c_char), INTENT(out) :: error_str(80) !< Needs docs
!---Lagrange mass solver
CLASS(oft_matrix), POINTER :: lmop => NULL()
CLASS(oft_solver), POINTER :: lminv => NULL()
!---Local variables
INTEGER(i4) :: i,io_unit,ierr
REAL(r8), POINTER, DIMENSION(:) :: vals => NULL()
REAL(r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: bvout
CLASS(oft_vector), POINTER :: u,v,check
TYPE(oft_hcurl_cinterp) :: Bfield
CHARACTER(LEN=3) :: pltnum
CALL mesh%setup_io(oft_hcurl%order)
!---Construct operator
NULLIFY(lmop)
CALL oft_lag_vgetmop(lmop,'none')
!---Setup solver
CALL create_cg_solver(lminv)
CALL create_diag_pre(lminv%pre)
lminv%A=>lmop
lminv%its=-2
!---Create solver fields
CALL oft_lag_vcreate(u)
CALL oft_lag_vcreate(v)
ALLOCATE(bvout(3,u%n/3))
!---Save modes
DO i=1,taylor_nm
  WRITE(pltnum,'(I3.3)')i
  ! CALL oft_hcurl%vec_save(taylor_hffa(i,oft_hcurl_level)%f, &
  !                         'taylor.rst','A_'//pltnum, append=(i/=1))
  !---Setup field interpolation
  Bfield%u=>taylor_hffa(i,oft_hcurl_level)%f
  CALL Bfield%setup
  !---Project field
  CALL oft_lag_vproject(Bfield,v)
  CALL u%set(0.d0)
  CALL lminv%apply(u,v)
  !---Retrieve local values and save
  vals=>bvout(1,:)
  CALL u%get_local(vals,1)
  vals=>bvout(2,:)
  CALL u%get_local(vals,2)
  vals=>bvout(3,:)
  CALL u%get_local(vals,3)
  call mesh%save_vertex_vector(bvout,'B_'//pltnum)
END DO
END SUBROUTINE marklin_save_visit
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_get_aint(imode,int_obj,error_str) BIND(C,NAME="marklin_get_aint")
INTEGER(KIND=c_int), VALUE, INTENT(in) :: imode !< Needs docs
TYPE(c_ptr), INTENT(out) :: int_obj !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(80) !< Needs docs
TYPE(oft_h1_rinterp), POINTER :: interp_obj
CLASS(oft_solver), POINTER :: linv => NULL()
TYPE(oft_h1_divout) :: divout
CLASS(oft_matrix), POINTER :: lop => NULL()
REAL(r8), POINTER, DIMENSION(:) :: tmp => NULL()
IF(oft_h1_nlevels==0)THEN
  !---H1(Grad) subspace
  CALL oft_h0_setup(oft_hcurl%order+1,oft_hcurl_minlev+1)
  CALL h0_setup_interp
  !---H1 full space
  CALL oft_h1_setup(oft_hcurl%order,oft_hcurl_minlev)
  CALL h1_setup_interp
END IF
!---------------------------------------------------------------------------
! Create divergence cleaner
!---------------------------------------------------------------------------
NULLIFY(lop,tmp)
CALL oft_h0_getlop(lop,"grnd")
CALL create_cg_solver(linv)
linv%A=>lop
linv%its=-2
CALL create_diag_pre(linv%pre) ! Setup Preconditioner
divout%solver=>linv
divout%bc=>h0_zerogrnd
divout%keep_boundary=.TRUE.
!---------------------------------------------------------------------------
! Setup initial conditions
!---------------------------------------------------------------------------
ALLOCATE(interp_obj)
CALL oft_h1_create(interp_obj%u)
CALL taylor_hffa(1,oft_hcurl_level)%f%get_local(tmp)
CALL interp_obj%u%restore_local(tmp,1)
WRITE(*,*)'Setting gauge'
divout%pm=.TRUE.
CALL divout%apply(interp_obj%u)
CALL interp_obj%setup()
int_obj=C_LOC(interp_obj)
!---Cleanup
CALL divout%delete()
END SUBROUTINE marklin_get_aint
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_get_bint(imode,int_obj,error_str) BIND(C,NAME="marklin_get_bint")
INTEGER(KIND=c_int), VALUE, INTENT(in) :: imode !< Needs docs
TYPE(c_ptr), INTENT(out) :: int_obj !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(80) !< Needs docs
TYPE(oft_hcurl_cinterp), POINTER :: interp_obj
ALLOCATE(interp_obj)
interp_obj%u=>taylor_hffa(imode,oft_hcurl_level)%f
CALL interp_obj%setup()
int_obj=C_LOC(interp_obj)
END SUBROUTINE marklin_get_bint
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_apply_int(int_obj,int_type,pt,fbary_tol,cell,field) BIND(C,NAME="marklin_apply_int")
TYPE(c_ptr), VALUE, INTENT(in) :: int_obj !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: int_type !< Needs docs
REAL(c_double), INTENT(in) :: pt(3) !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: fbary_tol !< Needs docs
INTEGER(c_int), INTENT(inout) :: cell !< Needs docs
REAL(c_double), INTENT(out) :: field(3) !< Needs docs
TYPE(oft_h1_rinterp), POINTER :: ainterp_obj
TYPE(oft_hcurl_cinterp), POINTER :: binterp_obj
REAL(8) :: f(4),goptmp(3,4),vol,fmin,fmax
IF(int_type<0)THEN
  SELECT CASE(int_type)
  CASE(1)
    CALL c_f_pointer(int_obj, ainterp_obj)
    CALL ainterp_obj%delete()
  CASE(2)
    CALL c_f_pointer(int_obj, binterp_obj)
    CALL binterp_obj%delete()
  END SELECT
  RETURN
END IF
call mesh_findcell(mesh,cell,pt,f)
IF(cell==0)RETURN
fmin=MINVAL(f); fmax=MAXVAL(f)
IF(( fmax>1.d0+fbary_tol ).OR.( fmin<-fbary_tol ))THEN
  cell=-ABS(cell)
  RETURN
END IF
CALL mesh%jacobian(cell,f,goptmp,vol)
SELECT CASE(int_type)
  CASE(1)
    CALL c_f_pointer(int_obj, ainterp_obj)
    CALL ainterp_obj%interp(cell,f,goptmp,field)
  CASE(2)
    CALL c_f_pointer(int_obj, binterp_obj)
    CALL binterp_obj%interp(cell,f,goptmp,field)
  CASE DEFAULT
    cell=-(mesh%np+1)
END SELECT
END SUBROUTINE marklin_apply_int
END MODULE marklin_f