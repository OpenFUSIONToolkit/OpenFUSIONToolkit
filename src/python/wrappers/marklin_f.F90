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
USE oft_lag_basis, ONLY: oft_lag_setup, ML_oft_lagrange
USE oft_lag_fields, ONLY: oft_lag_vcreate, oft_lag_create
USE oft_lag_operators, ONLY: lag_lop_eigs, lag_setup_interp, lag_mloptions, &
    oft_lag_vgetmop, oft_lag_vproject, oft_lag_getpdop, oft_lag_getmop
!---H1(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl, oft_hcurl_setup, oft_hcurl_level, &
  oft_hcurl_minlev, ML_oft_hcurl
USE oft_hcurl_operators, ONLY: oft_hcurl_cinterp, hcurl_setup_interp, &
    hcurl_mloptions
!---H1(Grad) FE space
USE oft_h0_basis, ONLY: oft_h0_setup, ML_oft_h0
USE oft_h0_operators, ONLY: h0_setup_interp, oft_h0_getlop, h0_zerogrnd, &
  h0_mloptions
!---H1 FE space
USE oft_h1_basis, ONLY: oft_h1_setup, oft_h1_nlevels, ML_oft_h1, oft_h1_level
USE oft_h1_fields, ONLY: oft_h1_create
USE oft_h1_operators, ONLY: oft_h1_divout, h1_zeroi, h1_mc, h1curl_zerob, &
  h1_setup_interp, oft_h1_rinterp
!---Taylor state
USE taylor, ONLY: taylor_minlev, taylor_hmodes, taylor_hffa, taylor_nm, &
  taylor_rst, taylor_vacuum, oft_taylor_rinterp, taylor_nh, taylor_hvac
USE mhd_utils, ONLY: mu0
!---Wrappers
USE oft_base_f, ONLY: copy_string, copy_string_rev
IMPLICIT NONE
!
integer(i4), POINTER :: lc_plot(:,:) !< Needs docs
integer(i4), POINTER :: reg_plot(:) !< Needs docs
real(r8), POINTER :: r_plot(:,:) !< Needs docs
CONTAINS
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_compute_eigs(order,minlev,nmodes,save_rst,error_str) BIND(C,NAME="marklin_compute_eigs")
INTEGER(KIND=c_int), VALUE, INTENT(in) :: order !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: minlev !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: nmodes !< Needs docs
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
!---Lagrange
IF(ML_oft_lagrange%nlevels==0)THEN
  CALL oft_lag_setup(order,minlev)
  IF(minlev>0)THEN
    CALL lag_setup_interp
    CALL lag_mloptions
  END IF
END IF
!---H1(Curl) subspace
IF(ML_oft_hcurl%nlevels==0)THEN
  CALL oft_hcurl_setup(order,minlev)
  IF(minlev>0)THEN
    CALL hcurl_setup_interp
    CALL hcurl_mloptions
  END IF
END IF
!---
IF(ML_oft_h0%nlevels==0)THEN
  CALL oft_h0_setup(order+1,minlev)
  IF(minlev>0)THEN
    CALL h0_setup_interp
    CALL h0_mloptions
  END IF
END IF
!---
IF(ML_oft_h1%nlevels==0)CALL oft_h1_setup(order,minlev)
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
END SUBROUTINE marklin_compute_eigs
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_compute_vac(order,minlev,nh,hcpc,hcpv,save_rst,error_str) BIND(C,NAME="marklin_compute_vac")
INTEGER(KIND=c_int), VALUE, INTENT(in) :: order !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: minlev !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: nh !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: hcpc !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: hcpv !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: save_rst !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(80) !< Needs docs
!---Lagrange mass solver
CLASS(oft_matrix), POINTER :: lmop => NULL()
CLASS(oft_solver), POINTER :: lminv => NULL()
!---Local variables
INTEGER(i4) :: i,io_unit,ierr
REAL(r8), POINTER, DIMENSION(:) :: vals => NULL()
REAL(r8), POINTER, DIMENSION(:,:) :: hcpc_tmp,hcpv_tmp
REAL(r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: bvout
CLASS(oft_vector), POINTER :: u,v,check
TYPE(oft_hcurl_cinterp) :: Bfield
CHARACTER(LEN=3) :: pltnum
!---Clear error flag
CALL copy_string('',error_str)
!---Lagrange
IF(ML_oft_lagrange%nlevels==0)THEN
  CALL oft_lag_setup(order,minlev)
  IF(minlev>0)THEN
    CALL lag_setup_interp
    CALL lag_mloptions
  END IF
END IF
!---H1(Curl) subspace
IF(ML_oft_hcurl%nlevels==0)THEN
  CALL oft_hcurl_setup(order,minlev)
  IF(minlev>0)THEN
    CALL hcurl_setup_interp
    CALL hcurl_mloptions
  END IF
END IF
!---
IF(ML_oft_h0%nlevels==0)THEN
  CALL oft_h0_setup(order+1,minlev)
  IF(minlev>0)THEN
    CALL h0_setup_interp
    CALL h0_mloptions
  END IF
END IF
!---
IF(ML_oft_h1%nlevels==0)CALL oft_h1_setup(order,minlev)
!---
IF(minlev<0)THEN
  taylor_minlev=oft_hcurl_level
ELSE
  taylor_minlev=minlev
  IF(oft_env%nprocs>1)taylor_minlev=MAX(oft_env%nbase+1,minlev)
END IF
oft_env%pm=.TRUE.
taylor_rst=save_rst
CALL c_f_pointer(hcpc, hcpc_tmp, [3,nh])
CALL c_f_pointer(hcpv, hcpv_tmp, [3,nh])
CALL taylor_vacuum(nh,hcpc_tmp,hcpv_tmp)
END SUBROUTINE marklin_compute_vac
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_compute_pardiff(int_obj,int_type,k_perp,error_str) BIND(C,NAME="marklin_compute_pardiff")
TYPE(c_ptr), VALUE, INTENT(in) :: int_obj !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: int_type !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: k_perp !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(80) !< Needs docs
!---Lagrange mass solver
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_matrix), POINTER :: pdop => NULL()
CLASS(oft_solver), POINTER :: pdinv => NULL()
!---Local variables
INTEGER(i4) :: i,io_unit,ierr
REAL(r8), POINTER, DIMENSION(:) :: vals => NULL()
REAL(r8), POINTER, DIMENSION(:,:) :: hcpc_tmp,hcpv_tmp
CLASS(oft_vector), POINTER :: u,v,check
TYPE(oft_hcurl_cinterp) :: Bfield
CHARACTER(LEN=3) :: pltnum
TYPE(oft_h1_rinterp), POINTER :: ainterp_obj
TYPE(oft_taylor_rinterp), POINTER :: binterp_obj
!---Clear error flag
CALL copy_string('',error_str)
!
NULLIFY(pdop,mop,vals)
SELECT CASE(int_type)
  CASE(1)
    CALL c_f_pointer(int_obj, ainterp_obj)
    CALL oft_lag_getpdop(pdop,ainterp_obj,'zerob',k_perp)
  CASE(2)
    CALL c_f_pointer(int_obj, binterp_obj)
    CALL oft_lag_getpdop(pdop,binterp_obj,'zerob',k_perp)
  CASE DEFAULT
    CALL copy_string('Invalid interpolation type',error_str)
    RETURN
END SELECT
!---Setup solver
CALL create_cg_solver(pdinv)
CALL create_diag_pre(pdinv%pre)
pdinv%A=>pdop
pdinv%its=-2
!---Create solver fields
CALL oft_lag_create(u)
CALL oft_lag_create(v)
CALL oft_lag_getmop(mop,'none')
CALL u%set(1.d0)
CALL mop%apply(u,v)
!
CALL u%set(0.d0)
CALL pdinv%apply(u,v)
!
IF(mesh%tess_order==0)CALL mesh%setup_io(oft_hcurl%order)
CALL u%get_local(vals)
CALL mesh%save_vertex_scalar(vals,'T')
!
CALL mop%delete()
CALL pdinv%pre%delete()
CALL pdinv%delete()
CALL pdop%delete()
CALL u%delete()
CALL v%delete()
DEALLOCATE(u,v,mop,pdop,pdinv)
END SUBROUTINE marklin_compute_pardiff
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
TYPE(oft_h1_rinterp) :: Bvfield
CHARACTER(LEN=3) :: pltnum
IF(mesh%tess_order==0)CALL mesh%setup_io(oft_hcurl%order)
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
  call mesh%save_vertex_vector(bvout,'Bh_'//pltnum)
END DO
IF(taylor_nm>0)CALL Bfield%delete()
!---Save vacuum fields
DO i=1,taylor_nh
  WRITE(pltnum,'(I3.3)')i
  !---Setup field interpolation
  Bvfield%u=>taylor_hvac(i,ML_oft_h1%level)%f
  CALL Bvfield%setup
  !---Project field
  CALL oft_lag_vproject(Bvfield,v)
  CALL u%set(0.d0)
  CALL lminv%apply(u,v)
  !---Retrieve local values and save
  vals=>bvout(1,:)
  CALL u%get_local(vals,1)
  vals=>bvout(2,:)
  CALL u%get_local(vals,2)
  vals=>bvout(3,:)
  CALL u%get_local(vals,3)
  call mesh%save_vertex_vector(bvout,'Bv_'//pltnum)
END DO
IF(taylor_nh>0)CALL Bvfield%delete()
CALL u%delete
CALL v%delete
DEALLOCATE(u,v,bvout)
END SUBROUTINE marklin_save_visit
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_get_aint(hmode_facs,int_obj,error_str) BIND(C,NAME="marklin_get_aint")
TYPE(c_ptr), VALUE, INTENT(in) :: hmode_facs !< Needs docs
TYPE(c_ptr), INTENT(out) :: int_obj !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(80) !< Needs docs
INTEGER(4) :: i
CLASS(oft_vector), POINTER :: utmp
TYPE(oft_h1_rinterp), POINTER :: interp_obj
CLASS(oft_solver), POINTER :: linv => NULL()
TYPE(oft_h1_divout) :: divout
CLASS(oft_matrix), POINTER :: lop => NULL()
REAL(r8), POINTER, DIMENSION(:) :: tmp,facs_tmp
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
CALL oft_h1_create(utmp)
CALL c_f_pointer(hmode_facs, facs_tmp, [taylor_nm])
DO i=1,taylor_nm
  CALL taylor_hffa(i,oft_hcurl_level)%f%get_local(tmp)
  CALL utmp%restore_local(tmp,1)
  CALL interp_obj%u%add(1.d0,facs_tmp(i),utmp)
END DO
CALL utmp%delete()
DEALLOCATE(tmp,utmp)
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
SUBROUTINE marklin_get_bint(hmode_facs,vac_facs,int_obj,error_str) BIND(C,NAME="marklin_get_bint")
TYPE(c_ptr), VALUE, INTENT(in) :: hmode_facs !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: vac_facs !< Needs docs
TYPE(c_ptr), INTENT(out) :: int_obj !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(80) !< Needs docs
INTEGER(4) :: i
REAL(r8), POINTER, DIMENSION(:) :: hmode_tmp,vac_tmp
TYPE(oft_taylor_rinterp), POINTER :: interp_obj
ALLOCATE(interp_obj)
CALL oft_h1_create(interp_obj%uvac)
CALL c_f_pointer(vac_facs, vac_tmp, [taylor_nh])
DO i=1,taylor_nh
  CALL interp_obj%uvac%add(1.d0,vac_tmp(i),taylor_hvac(i,ML_oft_h1%level)%f)
END DO
CALL ML_oft_hcurl%vec_create(interp_obj%ua)
CALL c_f_pointer(hmode_facs, hmode_tmp, [taylor_nm])
DO i=1,taylor_nm
  CALL interp_obj%ua%add(1.d0,hmode_tmp(i),taylor_hffa(i,oft_hcurl_level)%f)
END DO
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
TYPE(oft_taylor_rinterp), POINTER :: binterp_obj
REAL(8) :: f(4),goptmp(3,4),vol,fmin,fmax
IF(int_type<0)THEN
  SELECT CASE(int_type)
  CASE(1)
    CALL c_f_pointer(int_obj, ainterp_obj)
    CALL ainterp_obj%u%delete()
    CALL ainterp_obj%delete()
  CASE(2)
    CALL c_f_pointer(int_obj, binterp_obj)
    CALL binterp_obj%uvac%delete()
    CALL binterp_obj%ua%delete()
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