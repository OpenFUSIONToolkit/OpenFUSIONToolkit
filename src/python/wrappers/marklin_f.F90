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
USE multigrid, ONLY: multigrid_mesh, multigrid_reset
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
    oft_lag_vgetmop, oft_lag_vproject, oft_lag_getpdop, oft_lag_getmop
!---H1 FE space (Grad(H^1) subspace)
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp, oft_h1_getlop, oft_h1_zerogrnd, &
  oft_h1_zerob, h1_mloptions
!---Full H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: oft_hcurl_cinterp, hcurl_setup_interp, &
    hcurl_mloptions
USE oft_hcurl_grad_operators, ONLY: oft_hcurl_grad_divout, oft_hcurl_grad_zeroi, hcurl_grad_mc, oft_hcurl_grad_czerob, &
  hcurl_grad_setup_interp, oft_hcurl_grad_rinterp, hcurl_grad_mloptions
!---Taylor state
USE taylor, ONLY: oft_taylor_hmodes, oft_taylor_ifield, taylor_hmodes, oft_taylor_rinterp, &
  taylor_vacuum, taylor_injectors
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
  TYPE(oft_taylor_hmodes) :: eig_obj
  TYPE(oft_taylor_ifield) :: ff_obj
  TYPE(xdmf_plot_file) :: xdmf_plot
  TYPE(multigrid_mesh), POINTER :: ml_mesh => NULL()
  TYPE(oft_ml_fem_type), POINTER :: ML_lagrange => NULL()
  TYPE(oft_ml_fem_type), POINTER :: ML_h1 => NULL()
  TYPE(oft_ml_fem_type), POINTER :: ML_h1grad => NULL()
  TYPE(oft_ml_fem_type), POINTER :: ML_hcurl => NULL()
  TYPE(oft_ml_fem_comp_type), POINTER :: ML_hcurl_grad => NULL()
  TYPE(oft_ml_fem_comp_type), POINTER :: ML_vlagrange => NULL()
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
ALLOCATE(self%ML_lagrange,self%ML_vlagrange)
CALL oft_lag_setup(self%ml_mesh,self%order,self%ML_lagrange,ML_vlag_obj=self%ML_vlagrange,minlev=self%minlev)
!---H(Curl) subspace
ALLOCATE(self%ML_hcurl)
CALL oft_hcurl_setup(self%ml_mesh,self%order,self%ML_hcurl,minlev=self%minlev)
!---Compute modes
IF(self%minlev<0)THEN
  self%minlev=self%ML_hcurl%level
ELSE
  IF(oft_env%nprocs>1)self%minlev=MAX(oft_env%nbase+1,self%minlev)
END IF
IF(self%minlev<self%ML_hcurl%level)THEN
  CALL lag_setup_interp(self%ML_lagrange)
  CALL lag_mloptions
  CALL hcurl_setup_interp(self%ML_hcurl)
  CALL hcurl_mloptions(self%ML_hcurl)
END IF
!
marklin_ptr=C_LOC(self)
END SUBROUTINE marklin_setup
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_destroy(marklin_ptr,error_str) BIND(C,NAME="marklin_destroy")
TYPE(c_ptr), VALUE, INTENT(in) :: marklin_ptr !< Pointer to Marklin object
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error string (empty if no error)
INTEGER(4) :: i,ierr,io_unit,npts,iostat
REAL(8) :: theta
LOGICAL :: file_exists
real(r8), POINTER :: vals_tmp(:)
TYPE(marklin_obj), POINTER :: self
IF(.NOT.marklin_ccast(marklin_ptr,self,error_str))RETURN
!---Destroy objects
IF(ASSOCIATED(self%r_plot))DEALLOCATE(self%r_plot)
IF(ASSOCIATED(self%lc_plot))DEALLOCATE(self%lc_plot)
IF(ASSOCIATED(self%reg_plot))DEALLOCATE(self%reg_plot)
CALL self%eig_obj%delete()
CALL self%ff_obj%delete()
IF(ASSOCIATED(self%ML_lagrange))THEN
  CALL self%ML_lagrange%delete()
  DEALLOCATE(self%ML_lagrange)
END IF
IF(ASSOCIATED(self%ML_h1))THEN
  CALL self%ML_h1%delete()
  DEALLOCATE(self%ML_h1)
END IF
IF(ASSOCIATED(self%ML_hcurl))THEN
  CALL self%ML_hcurl%delete()
  DEALLOCATE(self%ML_hcurl)
END IF
IF(ASSOCIATED(self%ML_hcurl_grad))THEN
  CALL self%ML_hcurl_grad%delete()
  DEALLOCATE(self%ML_hcurl_grad)
END IF
! IF(ASSOCIATED(self%ML_h1grad))THEN
!   CALL self%ML_h1grad%delete()
!   DEALLOCATE(self%ML_h1grad)
! END IF
IF(ASSOCIATED(self%ml_mesh))THEN
  CALL multigrid_reset(self%ml_mesh)
  DEALLOCATE(self%ml_mesh)
END IF
DEALLOCATE(self)
END SUBROUTINE marklin_destroy
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_compute_eig(marklin_ptr,nmodes,eig_vals,cache_file,error_str) BIND(C,NAME="marklin_compute_eig")
TYPE(c_ptr), VALUE, INTENT(in) :: marklin_ptr !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: nmodes !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: eig_vals !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: cache_file(OFT_PATH_SLEN) !< Needs docs
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
CHARACTER(LEN=OFT_PATH_SLEN) :: rst_filename = ''
IF(.NOT.marklin_ccast(marklin_ptr,self,error_str))RETURN
IF(self%eig_obj%nm>0)THEN
  CALL copy_string('Eigenstates already computed',error_str)
  RETURN
END IF
oft_env%pm=.TRUE.
CALL self%eig_obj%setup(self%ML_hcurl,self%ML_lagrange,self%minlev)
CALL copy_string_rev(cache_file,rst_filename)
IF(TRIM(rst_filename)=='')THEN
  CALL taylor_hmodes(self%eig_obj,nmodes)
ELSE
  CALL taylor_hmodes(self%eig_obj,nmodes,rst_filename)
END IF
CALL c_f_pointer(eig_vals, vals_tmp, [nmodes])
vals_tmp=self%eig_obj%hlam(:,self%ML_hcurl%level)
END SUBROUTINE marklin_compute_eig
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_compute_vac(marklin_ptr,nh,hcpc,hcpv,cache_file,error_str) BIND(C,NAME="marklin_compute_vac")
TYPE(c_ptr), VALUE, INTENT(in) :: marklin_ptr !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: nh !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: hcpc !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: hcpv !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: cache_file(OFT_PATH_SLEN) !< Needs docs
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
TYPE(marklin_obj), POINTER :: self
CHARACTER(LEN=OFT_PATH_SLEN) :: rst_filename = ''
!---Clear error flag
CALL copy_string('',error_str)
IF(.NOT.c_associated(marklin_ptr))THEN
  CALL copy_string('Marklin object not associated',error_str)
  RETURN
END IF
CALL c_f_pointer(marklin_ptr,self)
!---Build FE spaces if not yet built
IF(.NOT.ASSOCIATED(self%ML_hcurl_grad))THEN
  ALLOCATE(self%ML_h1,self%ML_hcurl_grad,self%ML_h1grad)
  !---Grad(H^1) subspace
  CALL oft_h1_setup(self%ml_mesh,self%ML_hcurl%current_level%order+1,self%ML_h1,minlev=self%ML_hcurl%minlev)
  !---Full H(Curl) + Grad(H^1) space
  CALL oft_hcurl_grad_setup(self%ML_hcurl,self%ML_h1,self%ML_hcurl_grad,self%ML_h1grad,self%ML_hcurl%minlev)
  !---MG setup
  IF(self%minlev<self%ML_hcurl%level)THEN
    CALL h1_setup_interp(self%ML_h1)
    CALL h1_mloptions()
    CALL hcurl_grad_setup_interp(self%ML_hcurl_grad,self%ML_h1)
    CALL hcurl_grad_mloptions()
  END IF
END IF
!---
oft_env%pm=.TRUE.
CALL c_f_pointer(hcpc, hcpc_tmp, [3,nh])
CALL c_f_pointer(hcpv, hcpv_tmp, [3,nh])
CALL self%ff_obj%setup(nh,hcpc_tmp,hcpv_tmp,ML_hcurl=self%ML_hcurl,ML_h1=self%ML_h1,ML_hcurl_grad=self%ML_hcurl_grad, &
  ML_h1grad=self%ML_h1grad,ML_lagrange=self%ML_lagrange,minlev=self%minlev)
CALL copy_string_rev(cache_file,rst_filename)
IF(TRIM(rst_filename)=='')THEN
  CALL taylor_vacuum(self%ff_obj)
ELSE
  CALL taylor_vacuum(self%ff_obj,rst_filename=rst_filename)
END IF
END SUBROUTINE marklin_compute_vac
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_compute_pardiff(marklin_ptr,int_obj,int_type,k_perp,error_str) BIND(C,NAME="marklin_compute_pardiff")
TYPE(c_ptr), VALUE, INTENT(in) :: marklin_ptr !< Needs docs
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
CHARACTER(LEN=3) :: pltnum
TYPE(marklin_obj), POINTER :: self
CLASS(fem_interp), POINTER :: interp_obj
TYPE(oft_hcurl_grad_rinterp), POINTER :: ainterp_obj
TYPE(oft_taylor_rinterp), POINTER :: binterp_obj
!---Clear error flag
CALL copy_string('',error_str)
IF(.NOT.c_associated(marklin_ptr))THEN
  CALL copy_string('Marklin object not associated',error_str)
  RETURN
END IF
CALL c_f_pointer(marklin_ptr,self)
!---Check that interpolator is allocated
IF(.NOT.C_ASSOCIATED(int_obj))THEN
  CALL copy_string('Interpolation object not associated',error_str)
  RETURN
END IF
!
NULLIFY(pdop,mop,vals)
SELECT CASE(int_type)
  CASE(1)
    CALL c_f_pointer(int_obj, ainterp_obj)
    interp_obj=>ainterp_obj
  CASE(2)
    CALL c_f_pointer(int_obj, binterp_obj)
    interp_obj=>binterp_obj
  CASE DEFAULT
    CALL copy_string('Invalid interpolation type',error_str)
    RETURN
END SELECT
CALL oft_lag_getpdop(self%ML_lagrange%current_level,pdop,interp_obj,'zerob',k_perp)
!---Setup solver
CALL create_cg_solver(pdinv)
CALL create_diag_pre(pdinv%pre)
pdinv%A=>pdop
pdinv%its=-2
!---Create solver fields
CALL self%ML_lagrange%vec_create(u)
CALL self%ML_lagrange%vec_create(v)
CALL oft_lag_getmop(self%ML_lagrange%current_level,mop,'none')
CALL u%set(1.d0)
CALL mop%apply(u,v)
!
CALL u%set(0.d0)
CALL pdinv%apply(u,v)
!
CALL u%get_local(vals)
CALL self%ml_mesh%mesh%save_vertex_scalar(vals,self%xdmf_plot,'T')
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
  CALL self%xdmf_plot%setup('marklin',pathprefix)
  CALL self%ml_mesh%mesh%setup_io(self%xdmf_plot,self%ML_hcurl%current_level%order)
ELSE
  CALL self%xdmf_plot%setup('marklin')
  CALL self%ml_mesh%mesh%setup_io(self%xdmf_plot,self%ML_hcurl%current_level%order)
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
TYPE(oft_taylor_rinterp), POINTER :: binterp_obj
CHARACTER(LEN=80) :: name_tmp = ''
IF(.NOT.marklin_ccast(marklin_ptr,self,error_str))RETURN
CALL copy_string_rev(key,name_tmp)
!---Construct operator
NULLIFY(lmop)
CALL oft_lag_vgetmop(self%ML_vlagrange%current_level,lmop,'none')
!---Setup solver
CALL create_cg_solver(lminv)
CALL create_diag_pre(lminv%pre)
lminv%A=>lmop
lminv%its=-2
!---Create solver fields
CALL self%ML_vlagrange%vec_create(u)
CALL self%ML_vlagrange%vec_create(v)
ALLOCATE(bvout(3,u%n/3))
!---Project field onto plotting mesh
SELECT CASE(int_type)
CASE(1)
  CALL c_f_pointer(int_obj, ainterp_obj)
  CALL oft_lag_vproject(self%ML_lagrange%current_level,ainterp_obj,v)
CASE(2)
  CALL c_f_pointer(int_obj, binterp_obj)
  CALL oft_lag_vproject(self%ML_lagrange%current_level,binterp_obj,v)
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
SUBROUTINE marklin_get_aint(marklin_ptr,hmode_facs,int_obj,zero_norm,error_str) BIND(C,NAME="marklin_get_aint")
TYPE(c_ptr), VALUE, INTENT(in) :: marklin_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: hmode_facs !< Needs docs
TYPE(c_ptr), INTENT(out) :: int_obj !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: zero_norm !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
INTEGER(4) :: i
CLASS(oft_vector), POINTER :: utmp
TYPE(marklin_obj), POINTER :: self
TYPE(oft_hcurl_grad_rinterp), POINTER :: interp_obj
CLASS(oft_solver), POINTER :: linv => NULL()
TYPE(oft_hcurl_grad_divout) :: divout
CLASS(oft_matrix), POINTER :: lop => NULL()
REAL(r8), POINTER, DIMENSION(:) :: facs_tmp
REAL(r8), POINTER, DIMENSION(:) :: tmp => NULL()
TYPE(oft_h1_zerogrnd), TARGET :: h1_zerogrnd
TYPE(oft_h1_zerob), TARGET :: h1_zerob
IF(.NOT.marklin_ccast(marklin_ptr,self,error_str))RETURN
IF(.NOT.ASSOCIATED(self%ML_hcurl_grad))THEN
  ALLOCATE(self%ML_h1,self%ML_hcurl_grad,self%ML_h1grad)
  !---Grad(H^1) subspace
  CALL oft_h1_setup(self%ml_mesh,self%ML_hcurl%current_level%order+1,self%ML_h1,minlev=self%ML_hcurl%minlev+1)
  CALL h1_setup_interp(self%ML_h1)
  !---Full H(Curl) + Grad(H^1) space
  CALL oft_hcurl_grad_setup(self%ML_hcurl,self%ML_h1,self%ML_hcurl_grad,self%ML_h1grad,self%ML_hcurl%minlev)
  CALL hcurl_grad_setup_interp(self%ML_hcurl_grad,self%ML_h1)
END IF
!---------------------------------------------------------------------------
! Create divergence cleaner
!---------------------------------------------------------------------------
NULLIFY(lop,tmp)
IF(zero_norm)THEN
  CALL oft_h1_getlop(self%ML_h1%current_level,lop,"grnd")
ELSE
  CALL oft_h1_getlop(self%ML_h1%current_level,lop,"zerob")
END IF
CALL create_cg_solver(linv)
linv%A=>lop
linv%its=-2
CALL create_diag_pre(linv%pre) ! Setup Preconditioner
CALL divout%setup(self%ML_hcurl_grad,'none',solver=linv)
IF(zero_norm)THEN
  h1_zerogrnd%ML_H1_rep=>self%ML_h1grad
  divout%bc=>h1_zerogrnd
  divout%keep_boundary=.TRUE.
ELSE
  h1_zerob%ML_H1_rep=>self%ML_h1grad
  divout%bc=>h1_zerob
END IF
!---------------------------------------------------------------------------
! Setup initial conditions
!---------------------------------------------------------------------------
ALLOCATE(interp_obj)
CALL self%ML_hcurl_grad%vec_create(interp_obj%u)
CALL self%ML_hcurl_grad%vec_create(utmp)
CALL c_f_pointer(hmode_facs, facs_tmp, [self%eig_obj%nm])
DO i=1,self%eig_obj%nm
  CALL self%eig_obj%hffa(i,self%ML_hcurl%level)%f%get_local(tmp)
  CALL utmp%restore_local(tmp,1)
  CALL interp_obj%u%add(1.d0,facs_tmp(i),utmp)
END DO
CALL utmp%delete()
DEALLOCATE(tmp,utmp)
IF(zero_norm)WRITE(*,*)'Setting gauge'
divout%pm=.TRUE.
CALL divout%apply(interp_obj%u)
CALL interp_obj%setup(self%ML_hcurl%current_level,self%ML_h1%current_level)
int_obj=C_LOC(interp_obj)
!---Cleanup
CALL lop%delete()
DEALLOCATE(lop)
CALL linv%pre%delete()
DEALLOCATE(linv%pre)
CALL linv%delete()
DEALLOCATE(linv)
CALL divout%bc%delete()
NULLIFY(divout%bc)
CALL divout%delete()
END SUBROUTINE marklin_get_aint
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE marklin_get_bint(marklin_ptr,hmode_facs,vac_facs,int_obj,error_str) BIND(C,NAME="marklin_get_bint")
TYPE(c_ptr), VALUE, INTENT(in) :: marklin_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: hmode_facs !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: vac_facs !< Needs docs
TYPE(c_ptr), INTENT(out) :: int_obj !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
INTEGER(4) :: i
REAL(r8), POINTER, DIMENSION(:) :: hmode_tmp,vac_tmp
TYPE(marklin_obj), POINTER :: self
TYPE(oft_taylor_rinterp), POINTER :: interp_obj
IF(.NOT.marklin_ccast(marklin_ptr,self,error_str))RETURN
ALLOCATE(interp_obj)
CALL self%ML_hcurl_grad%vec_create(interp_obj%uvac)
CALL c_f_pointer(vac_facs, vac_tmp, [self%ff_obj%nh])
DO i=1,self%ff_obj%nh
  CALL interp_obj%uvac%add(1.d0,vac_tmp(i),self%ff_obj%hvac(i)%f)
END DO
CALL self%ML_hcurl%vec_create(interp_obj%ua)
CALL c_f_pointer(hmode_facs, hmode_tmp, [self%eig_obj%nm])
DO i=1,self%eig_obj%nm
  CALL interp_obj%ua%add(1.d0,hmode_tmp(i),self%eig_obj%hffa(i,self%ML_hcurl%level)%f)
END DO
CALL interp_obj%setup(self%ML_hcurl_grad%current_level)
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
TYPE(oft_taylor_rinterp), POINTER :: binterp_obj
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