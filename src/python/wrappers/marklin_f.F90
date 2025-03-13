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
    c_f_pointer, c_bool, c_null_ptr, c_associated
!---Base
USE oft_base
USE oft_io, ONLY: hdf5_create_file, xdmf_plot_file
!--Grid
USE oft_mesh_type, ONLY: mesh_findcell
USE oft_mesh_native, ONLY: r_mem, lc_mem, reg_mem
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear Algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE fem_utils, ONLY: fem_interp
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_lop_eigs, lag_setup_interp, lag_mloptions, &
    oft_lag_vgetmop, oft_lag_vproject
!---H1 FE space (Grad(H^1) subspace)
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp, oft_h1_getlop, oft_h1_zerogrnd, &
  oft_h1_zerob
!---Full H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: oft_hcurl_cinterp, hcurl_setup_interp, &
    hcurl_mloptions
USE oft_hcurl_grad_operators, ONLY: oft_hcurl_grad_divout, oft_hcurl_grad_zeroi, hcurl_grad_mc, oft_hcurl_grad_czerob, &
  hcurl_grad_setup_interp, oft_hcurl_grad_rinterp
!---Taylor state
USE taylor, ONLY: taylor_minlev, taylor_hmodes, taylor_hffa, taylor_nm, &
  taylor_rst, taylor_hlam, ML_oft_hcurl, ML_oft_h1, &
  ML_hcurl_grad, ML_h1grad, ML_oft_lagrange, ML_oft_vlagrange
USE mhd_utils, ONLY: mu0
!---Wrappers
USE oft_base_f, ONLY: copy_string, copy_string_rev
IMPLICIT NONE
#include "local.h"
TYPE :: marklin_obj
  INTEGER(i4) :: order = 1
  INTEGER(i4) :: minlev = 1
  INTEGER(i4), POINTER, DIMENSION(:) :: reg_plot => NULL() !< Needs docs
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lc_plot => NULL() !< Needs docs
  REAL(r8), POINTER, DIMENSION(:,:) :: r_plot => NULL() !< Needs docs
  TYPE(xdmf_plot_file) :: xdmf_plot
  TYPE(multigrid_mesh), POINTER :: ml_mesh => NULL()
  ! TYPE(oft_ml_fem_type), POINTER :: ML_oft_lagrange => NULL()
  ! TYPE(oft_ml_fem_type), POINTER :: ML_oft_h1 => NULL()
  ! TYPE(oft_ml_fem_type), POINTER :: ML_oft_hgrad => NULL()
  ! TYPE(oft_ml_fem_type), POINTER :: ML_oft_hcurl => NULL()
  ! TYPE(oft_ml_fem_comp_type), POINTER :: ML_oft_h1 => NULL()
  ! TYPE(oft_ml_fem_comp_type), POINTER :: ML_oft_vlagrange => NULL()
END TYPE marklin_obj
CONTAINS
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
FUNCTION marklin_ccast(marklin_cptr,marklin_fptr,error_str) RESULT(success)
TYPE(c_ptr), INTENT(in) :: marklin_cptr !< C pointer to TokaMaker object
TYPE(marklin_obj), POINTER, INTENT(out) :: marklin_fptr
CHARACTER(KIND=c_char), OPTIONAL, INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
LOGICAL :: success
!---Clear error flag
IF(PRESENT(error_str))CALL copy_string('',error_str)
IF(.NOT.c_associated(marklin_cptr))THEN
  IF(PRESENT(error_str))CALL copy_string('TokaMaker object not associated',error_str)
  success=.FALSE.
  RETURN
END IF
CALL c_f_pointer(marklin_cptr,marklin_fptr)
success=.TRUE.
END FUNCTION marklin_ccast
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_setup(marklin_ptr,mesh_ptr,order,minlev,error_str) BIND(C,NAME="marklin_setup")
TYPE(c_ptr), INTENT(out) :: marklin_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: mesh_ptr !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: order !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: minlev !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
!---Local variables
INTEGER(i4) :: i,io_unit,ierr
TYPE(marklin_obj), POINTER :: self
!---Clear error flag
CALL copy_string('',error_str)
IF(.NOT.c_associated(mesh_ptr))THEN
  CALL copy_string('Mesh object not associated',error_str)
  RETURN
END IF
ALLOCATE(self)
self%order=order
self%minlev=minlev
CALL c_f_pointer(mesh_ptr,self%ml_mesh)
!---Lagrange
CALL oft_lag_setup(self%ml_mesh,self%order,ML_oft_lagrange,ML_vlag_obj=ML_oft_vlagrange,minlev=self%minlev)
!---H(Curl) subspace
CALL oft_hcurl_setup(self%ml_mesh,self%order,ML_oft_hcurl,minlev=self%minlev)
!---Compute modes
IF(self%minlev<0)THEN
  taylor_minlev=ML_oft_hcurl%level
ELSE
  taylor_minlev=self%minlev
  IF(oft_env%nprocs>1)taylor_minlev=MAX(oft_env%nbase+1,self%minlev)
END IF
IF(taylor_minlev<ML_oft_hcurl%level)THEN
  CALL lag_setup_interp(ML_oft_lagrange)
  CALL lag_mloptions
  CALL hcurl_setup_interp(ML_oft_hcurl)
  CALL hcurl_mloptions(ML_oft_hcurl)
END IF
!
marklin_ptr=C_LOC(self)
END SUBROUTINE marklin_setup
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_compute(marklin_ptr,nmodes,save_rst,eig_vals,error_str) BIND(C,NAME="marklin_compute")
TYPE(c_ptr), VALUE, INTENT(in) :: marklin_ptr !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: nmodes !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: save_rst !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: eig_vals !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
!---Lagrange mass solver
CLASS(oft_matrix), POINTER :: lmop => NULL()
CLASS(oft_solver), POINTER :: lminv => NULL()
!---Local variables
INTEGER(i4) :: i,io_unit,ierr
REAL(r8), POINTER, DIMENSION(:) :: vals_tmp => NULL()
REAL(r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: bvout
CLASS(oft_vector), POINTER :: u,v,check
TYPE(oft_hcurl_cinterp) :: Bfield
TYPE(marklin_obj), POINTER :: self
CHARACTER(LEN=3) :: pltnum
IF(.NOT.marklin_ccast(marklin_ptr,self,error_str))RETURN
IF(taylor_nm>0)THEN
  CALL copy_string('Eigenstates already computed',error_str)
  RETURN
END IF
oft_env%pm=.TRUE.
taylor_rst=save_rst
CALL taylor_hmodes(nmodes)
CALL c_f_pointer(eig_vals, vals_tmp, [nmodes])
vals_tmp=taylor_hlam(:,ML_oft_hcurl%level)
END SUBROUTINE marklin_compute
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_setup_io(marklin_ptr,basepath,error_str) BIND(C,NAME="marklin_setup_io")
TYPE(c_ptr), VALUE, INTENT(in) :: marklin_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: basepath(OFT_PATH_SLEN) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
TYPE(marklin_obj), POINTER :: self
CHARACTER(LEN=OFT_PATH_SLEN) :: pathprefix = ''
IF(.NOT.marklin_ccast(marklin_ptr,self,error_str))RETURN
CALL copy_string_rev(basepath,pathprefix)
!---Setup I/0
IF(TRIM(pathprefix)/='')THEN
  CALL self%xdmf_plot%setup('Marklin',pathprefix)
  CALL self%ml_mesh%mesh%setup_io(self%xdmf_plot,ML_oft_hcurl%current_level%order)
ELSE
  CALL self%xdmf_plot%setup('Marklin')
  CALL self%ml_mesh%mesh%setup_io(self%xdmf_plot,ML_oft_hcurl%current_level%order)
END IF
END SUBROUTINE marklin_setup_io
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_save_visit(marklin_ptr,int_obj,int_type,key,error_str) BIND(C,NAME="marklin_save_visit")
TYPE(c_ptr), VALUE, INTENT(in) :: marklin_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: int_obj !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: int_type !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: key(OFT_SLEN) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
!---Lagrange mass solver
CLASS(oft_matrix), POINTER :: lmop => NULL()
CLASS(oft_solver), POINTER :: lminv => NULL()
!---Local variables
INTEGER(i4) :: i
REAL(r8), POINTER, DIMENSION(:) :: vals => NULL()
REAL(r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: bvout
TYPE(marklin_obj), POINTER :: self
CLASS(oft_vector), POINTER :: u,v
TYPE(oft_hcurl_grad_rinterp), POINTER :: ainterp_obj
TYPE(oft_hcurl_cinterp), POINTER :: binterp_obj
CHARACTER(LEN=80) :: name_tmp = ''
IF(.NOT.marklin_ccast(marklin_ptr,self,error_str))RETURN
CALL copy_string_rev(key,name_tmp)
!---Construct operator
NULLIFY(lmop)
CALL oft_lag_vgetmop(ML_oft_vlagrange%current_level,lmop,'none')
!---Setup solver
CALL create_cg_solver(lminv)
CALL create_diag_pre(lminv%pre)
lminv%A=>lmop
lminv%its=-2
!---Create solver fields
CALL ML_oft_vlagrange%vec_create(u)
CALL ML_oft_vlagrange%vec_create(v)
ALLOCATE(bvout(3,u%n/3))
!---Project field onto plotting mesh
SELECT CASE(int_type)
CASE(1)
  CALL c_f_pointer(int_obj, ainterp_obj)
  CALL oft_lag_vproject(ML_oft_lagrange%current_level,ainterp_obj,v)
CASE(2)
  CALL c_f_pointer(int_obj, binterp_obj)
  CALL oft_lag_vproject(ML_oft_lagrange%current_level,binterp_obj,v)
END SELECT
CALL u%set(0.d0)
CALL lminv%apply(u,v)
!---Retrieve local values and save
vals=>bvout(1,:)
CALL u%get_local(vals,1)
vals=>bvout(2,:)
CALL u%get_local(vals,2)
vals=>bvout(3,:)
CALL u%get_local(vals,3)
call self%ml_mesh%mesh%save_vertex_vector(bvout,self%xdmf_plot,TRIM(name_tmp))
!---Cleanup
CALL lminv%pre%delete
DEALLOCATE(lminv%pre)
CALL u%delete()
CALL v%delete()
CALL lmop%delete()
CALL lminv%delete()
DEALLOCATE(u,v,lmop,lminv)
END SUBROUTINE marklin_save_visit
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_get_aint(marklin_ptr,imode,int_obj,zero_norm,error_str) BIND(C,NAME="marklin_get_aint")
TYPE(c_ptr), VALUE, INTENT(in) :: marklin_ptr !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: imode !< Needs docs
TYPE(c_ptr), INTENT(out) :: int_obj !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: zero_norm !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
TYPE(marklin_obj), POINTER :: self
TYPE(oft_hcurl_grad_rinterp), POINTER :: interp_obj
CLASS(oft_solver), POINTER :: linv => NULL()
TYPE(oft_hcurl_grad_divout) :: divout
CLASS(oft_matrix), POINTER :: lop => NULL()
REAL(r8), POINTER, DIMENSION(:) :: tmp => NULL()
TYPE(oft_h1_zerogrnd), TARGET :: h1_zerogrnd
TYPE(oft_h1_zerob), TARGET :: h1_zerob
IF(.NOT.marklin_ccast(marklin_ptr,self,error_str))RETURN
IF(ML_hcurl_grad%nlevels==0)THEN
  !---Grad(H^1) subspace
  CALL oft_h1_setup(self%ml_mesh,ML_oft_hcurl%current_level%order+1,ML_oft_h1,minlev=ML_oft_hcurl%minlev+1)
  CALL h1_setup_interp(ML_oft_h1)
  !---Full H(Curl) + Grad(H^1) space
  CALL oft_hcurl_grad_setup(ML_oft_hcurl,ML_oft_h1,ML_hcurl_grad,ML_h1grad,ML_oft_hcurl%minlev)
  CALL hcurl_grad_setup_interp(ML_hcurl_grad,ML_oft_h1)
END IF
!---------------------------------------------------------------------------
! Create divergence cleaner
!---------------------------------------------------------------------------
NULLIFY(lop,tmp)
IF(zero_norm)THEN
  CALL oft_h1_getlop(ML_oft_h1%current_level,lop,"grnd")
ELSE
  CALL oft_h1_getlop(ML_oft_h1%current_level,lop,"zerob")
END IF
CALL create_cg_solver(linv)
linv%A=>lop
linv%its=-2
CALL create_diag_pre(linv%pre) ! Setup Preconditioner
CALL divout%setup(ML_hcurl_grad,'none',solver=linv)
IF(zero_norm)THEN
  h1_zerogrnd%ML_H1_rep=>ML_h1grad
  divout%bc=>h1_zerogrnd
  divout%keep_boundary=.TRUE.
ELSE
  h1_zerob%ML_H1_rep=>ML_h1grad
  divout%bc=>h1_zerob
END IF
!---------------------------------------------------------------------------
! Setup initial conditions
!---------------------------------------------------------------------------
ALLOCATE(interp_obj)
CALL ML_hcurl_grad%vec_create(interp_obj%u)
CALL taylor_hffa(1,ML_oft_hcurl%level)%f%get_local(tmp)
CALL interp_obj%u%restore_local(tmp,1)
IF(zero_norm)WRITE(*,*)'Setting gauge'
divout%pm=.TRUE.
CALL divout%apply(interp_obj%u)
CALL interp_obj%setup(ML_oft_hcurl%current_level,ML_oft_h1%current_level)
int_obj=C_LOC(interp_obj)
!---Cleanup
CALL lop%delete()
DEALLOCATE(lop)
CALL linv%pre%delete()
DEALLOCATE(linv%pre)
CALL linv%delete()
DEALLOCATE(linv,tmp)
CALL divout%bc%delete()
NULLIFY(divout%bc)
CALL divout%delete()
END SUBROUTINE marklin_get_aint
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_get_bint(marklin_ptr,imode,int_obj,error_str) BIND(C,NAME="marklin_get_bint")
TYPE(c_ptr), VALUE, INTENT(in) :: marklin_ptr !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: imode !< Needs docs
TYPE(c_ptr), INTENT(out) :: int_obj !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
TYPE(marklin_obj), POINTER :: self
TYPE(oft_hcurl_cinterp), POINTER :: interp_obj
IF(.NOT.marklin_ccast(marklin_ptr,self,error_str))RETURN
ALLOCATE(interp_obj)
interp_obj%u=>taylor_hffa(imode,ML_oft_hcurl%level)%f
CALL interp_obj%setup(ML_oft_hcurl%current_level)
int_obj=C_LOC(interp_obj)
END SUBROUTINE marklin_get_bint
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_apply_int(marklin_ptr,int_obj,int_type,pt,fbary_tol,cell,field) BIND(C,NAME="marklin_apply_int")
TYPE(c_ptr), VALUE, INTENT(in) :: marklin_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: int_obj !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: int_type !< Needs docs
REAL(c_double), INTENT(in) :: pt(3) !< Needs docs
REAL(c_double), VALUE, INTENT(in) :: fbary_tol !< Needs docs
INTEGER(c_int), INTENT(inout) :: cell !< Needs docs
REAL(c_double), INTENT(out) :: field(3) !< Needs docs
TYPE(marklin_obj), POINTER :: self
TYPE(oft_hcurl_grad_rinterp), POINTER :: ainterp_obj
TYPE(oft_hcurl_cinterp), POINTER :: binterp_obj
REAL(8) :: f(4),goptmp(3,4),vol,fmin,fmax
IF(.NOT.c_associated(marklin_ptr))THEN
  cell=-1
  RETURN
END IF
CALL c_f_pointer(marklin_ptr,self)
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
call mesh_findcell(self%ml_mesh%mesh,cell,pt,f)
IF(cell==0)RETURN
fmin=MINVAL(f); fmax=MAXVAL(f)
IF(( fmax>1.d0+fbary_tol ).OR.( fmin<-fbary_tol ))THEN
  cell=-ABS(cell)
  RETURN
END IF
CALL self%ml_mesh%mesh%jacobian(cell,f,goptmp,vol)
SELECT CASE(int_type)
  CASE(1)
    CALL c_f_pointer(int_obj, ainterp_obj)
    CALL ainterp_obj%interp(cell,f,goptmp,field)
  CASE(2)
    CALL c_f_pointer(int_obj, binterp_obj)
    CALL binterp_obj%interp(cell,f,goptmp,field)
  CASE DEFAULT
    cell=-(self%ml_mesh%mesh%np+1)
END SELECT
END SUBROUTINE marklin_apply_int
END MODULE marklin_f