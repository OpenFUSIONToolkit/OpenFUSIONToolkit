!---------------------------------------------------------------------------
!> @file thincurr_f.F90
!
!> Fortran part of Python wrapper for ThinCurr thin-wall eddy current functionality
!!
!! @authors Chris Hansen
!! @date March 2024
!! @ingroup doxy_oft_python
!---------------------------------------------------------------------------
MODULE thincurr_f
USE iso_c_binding, ONLY: c_int, c_double, c_char, c_loc, c_null_char, c_ptr, &
    c_f_pointer, c_bool, c_null_ptr, c_associated
!---Base
USE oft_base
USE spline_mod
USE oft_io, ONLY: hdf5_create_file, hdf5_field_get_sizes, hdf5_read, hdf5_field_exist
!--Grid
USE oft_trimesh_type, ONLY: oft_trimesh
USE oft_mesh_native, ONLY: r_mem, lc_mem, reg_mem, native_read_nodesets, native_read_sidesets
!---Linear Algebra
USE oft_la_base, ONLY: oft_vector
!---
USE fem_utils, ONLY: fem_interp
USE thin_wall, ONLY: tw_type, tw_save_pfield, tw_compute_LmatDirect, tw_compute_Rmat, &
  tw_compute_Ael2dr, tw_sensors, tw_compute_mutuals, tw_load_sensors, tw_compute_Lmat_MF
USE thin_wall_hodlr, ONLY: oft_tw_hodlr_op
USE thin_wall_solvers, ONLY: lr_eigenmodes_arpack, lr_eigenmodes_direct, frequency_response
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
SUBROUTINE thincurr_setup(mesh_file,np,r_loc,nc,lc_loc,reg_loc,pmap_loc,tw_ptr,sizes,error_str,xml_ptr) BIND(C,NAME="thincurr_setup")
CHARACTER(KIND=c_char), INTENT(in) :: mesh_file(80) !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: r_loc !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: lc_loc !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: reg_loc !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: pmap_loc !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: np !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: nc !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: sizes !< Needs docs
TYPE(c_ptr), INTENT(out) :: tw_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(200) !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: xml_ptr !< Needs docs
!
LOGICAL :: success,is_2d
INTEGER(4) :: ndims,ierr
integer(i4), allocatable, dimension(:) :: dim_sizes
INTEGER(i4), POINTER, DIMENSION(:) :: sizes_tmp,pmap_tmp
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: rtmp
CHARACTER(LEN=OFT_PATH_SLEN) :: filename = 'none'
TYPE(tw_type), POINTER :: tw_obj
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_nsets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_ssets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: hole_nsets => NULL()
#ifdef HAVE_XML
TYPE(fox_node), POINTER :: xml_node
#endif
CALL copy_string('',error_str)
CALL copy_string_rev(mesh_file,filename)
CALL c_f_pointer(sizes, sizes_tmp, [8])
IF(TRIM(filename)=='none')THEN
  CALL copy_string('Array-based loading not yet supported',error_str)
  RETURN
  ! CALL c_f_pointer(pmap_loc, pmap_tmp, [1])
  ! IF(pmap_tmp(1)/=-1)THEN
  !   CALL c_f_pointer(pmap_loc, pmap_tmp, [tw_obj%mesh%np])
  !   ALLOCATE(tw_obj%pmap(tw_obj%mesh%np))
  !   tw_obj%pmap=pmap_tmp
  ! END IF
ELSE
  !---------------------------------------------------------------------------
  ! Load model from file
  !---------------------------------------------------------------------------
  ALLOCATE(tw_obj)
  ALLOCATE(oft_trimesh::tw_obj%mesh)
  CALL tw_obj%mesh%setup(-1,.FALSE.)
  CALL hdf5_field_get_sizes(TRIM(filename),"mesh/R",ndims,dim_sizes)
  tw_obj%mesh%np=dim_sizes(2)
  is_2d=(dim_sizes(1)==2)
  DEALLOCATE(dim_sizes)
  CALL hdf5_field_get_sizes(TRIM(filename),"mesh/LC",ndims,dim_sizes)
  tw_obj%mesh%nc=dim_sizes(2)
  DEALLOCATE(dim_sizes)
  ALLOCATE(tw_obj%mesh%r(3,tw_obj%mesh%np))
  IF(is_2d)THEN
    tw_obj%mesh%r=0.d0
    ALLOCATE(rtmp(2,tw_obj%mesh%np))
    CALL hdf5_read(rtmp,TRIM(filename),"mesh/R",success)
    tw_obj%mesh%r(1:2,:)=rtmp
    DEALLOCATE(rtmp)
  ELSE
    CALL hdf5_read(tw_obj%mesh%r,TRIM(filename),"mesh/R",success)
  END IF
  ALLOCATE(tw_obj%mesh%lc(3,tw_obj%mesh%nc))
  CALL hdf5_read(tw_obj%mesh%lc,TRIM(filename),"mesh/LC",success)
  ALLOCATE(tw_obj%mesh%reg(tw_obj%mesh%nc))
  tw_obj%mesh%reg=1.d0
  !
  IF(hdf5_field_exist(TRIM(filename),'thincurr/periodicity/pmap'))THEN
    ALLOCATE(tw_obj%pmap(tw_obj%mesh%np))
    CALL hdf5_read(tw_obj%pmap,TRIM(filename),'thincurr/periodicity/pmap',success)
  END IF
  !
  NULLIFY(mesh_nsets,mesh_ssets,hole_nsets)
  CALL native_read_nodesets(mesh_nsets,native_filename=TRIM(filename))
  CALL native_read_sidesets(mesh_ssets,native_filename=TRIM(filename))
  IF(ASSOCIATED(mesh_ssets))THEN
    IF(mesh_ssets(1)%n>0)THEN
      tw_obj%nclosures=mesh_ssets(1)%n
      ALLOCATE(tw_obj%closures(tw_obj%nclosures))
      tw_obj%closures=mesh_ssets(1)%v
    END IF
  END IF
  hole_nsets=>mesh_nsets
END IF
IF(c_associated(xml_ptr))THEN
#ifdef HAVE_XML
  CALL c_f_pointer(xml_ptr, xml_node)
  CALL xml_get_element(xml_node,"thincurr",tw_obj%xml,ierr)
  IF(ierr/=0)THEN
    CALL copy_string('Error getting ThinCurr XML node',error_str)
    RETURN
  END IF
#else
  CALL copy_string('OFT not compiled with XML support',error_str)
  RETURN
#endif
END IF
CALL tw_obj%setup(hole_nsets)
!
tw_ptr=C_LOC(tw_obj)
sizes_tmp=[tw_obj%mesh%np,tw_obj%mesh%ne,tw_obj%mesh%nc,tw_obj%np_active, &
  tw_obj%nholes,tw_obj%n_vcoils,tw_obj%nelems,tw_obj%n_icoils]
END SUBROUTINE thincurr_setup
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_setup_io(tw_ptr,basepath,save_debug,error_str) BIND(C,NAME="thincurr_setup_io")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: basepath(80) !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: save_debug !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(200) !< Needs docs
!
CHARACTER(LEN=OFT_PATH_SLEN) :: pathprefix = ''
TYPE(tw_type), POINTER :: tw_obj
CALL copy_string('',error_str)
CALL c_f_pointer(tw_ptr, tw_obj)
CALL copy_string_rev(basepath,pathprefix)
!---Setup I/0
IF(TRIM(pathprefix)/='')THEN
  CALL tw_obj%mesh%setup_io(1,basepath=pathprefix)
ELSE
  CALL tw_obj%mesh%setup_io(1)
END IF
IF(save_debug)CALL tw_obj%save_debug()
END SUBROUTINE thincurr_setup_io
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_save_field(tw_ptr,vals,fieldname) BIND(C,NAME="thincurr_save_field")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: vals !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: fieldname(80) !< Needs docs
!
CHARACTER(LEN=80) :: name_tmp = ''
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
TYPE(tw_type), POINTER :: tw_obj
CALL c_f_pointer(tw_ptr, tw_obj)
CALL copy_string_rev(fieldname,name_tmp)
CALL c_f_pointer(vals, vals_tmp, [tw_obj%nelems])
!---Save plot fields
CALL tw_save_pfield(tw_obj,vals_tmp,TRIM(name_tmp))
END SUBROUTINE thincurr_save_field
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_save_scalar(tw_ptr,vals,fieldname) BIND(C,NAME="thincurr_save_scalar")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: vals !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: fieldname(80) !< Needs docs
!
CHARACTER(LEN=80) :: name_tmp = ''
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
TYPE(tw_type), POINTER :: tw_obj
CALL c_f_pointer(tw_ptr, tw_obj)
CALL copy_string_rev(fieldname,name_tmp)
CALL c_f_pointer(vals, vals_tmp, [tw_obj%mesh%np])
!---Save plot fields
CALL tw_obj%mesh%save_vertex_scalar(vals_tmp,TRIM(name_tmp))
END SUBROUTINE thincurr_save_scalar
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_cross_coupling(tw_ptr1,tw_ptr2,Mmat,cache_file,error_str) BIND(C,NAME="thincurr_cross_coupling")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr1 !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr2 !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: Mmat !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: cache_file(80) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(200) !< Needs docs
!
REAL(8), CONTIGUOUS, POINTER, DIMENSION(:,:) :: Mmat_tmp
CHARACTER(LEN=OFT_PATH_SLEN) :: filename = ''
TYPE(tw_type), POINTER :: tw_obj1,tw_obj2
CALL copy_string('',error_str)
CALL c_f_pointer(tw_ptr1, tw_obj1)
CALL c_f_pointer(tw_ptr2, tw_obj2)
CALL c_f_pointer(Mmat, Mmat_tmp, [tw_obj2%nelems,tw_obj1%nelems])
!
CALL copy_string_rev(cache_file,filename)
IF(TRIM(filename)=='')THEN
  CALL tw_compute_LmatDirect(tw_obj1,Mmat_tmp,col_model=tw_obj2)
ELSE
  CALL tw_compute_LmatDirect(tw_obj1,Mmat_tmp,col_model=tw_obj2,save_file=filename)
END IF
END SUBROUTINE thincurr_cross_coupling
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_cross_eval(tw_ptr1,tw_ptr2,nrhs,vec1,vec2,error_str) BIND(C,NAME="thincurr_cross_eval")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr1 !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr2 !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: nrhs
TYPE(c_ptr), VALUE, INTENT(in) :: vec1 !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: vec2 !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(200) !< Needs docs
!
REAL(8), CONTIGUOUS, POINTER, DIMENSION(:,:) :: vec1_tmp,vec2_tmp
TYPE(tw_type), POINTER :: tw_obj1,tw_obj2
CALL copy_string('',error_str)
CALL c_f_pointer(tw_ptr1, tw_obj1)
CALL c_f_pointer(tw_ptr2, tw_obj2)
CALL c_f_pointer(vec1, vec1_tmp, [tw_obj1%nelems,nrhs])
CALL c_f_pointer(vec2, vec2_tmp, [tw_obj2%nelems,nrhs])
!
CALL tw_compute_Lmat_MF(tw_obj1,tw_obj2,nrhs,vec1_tmp,vec2_tmp)
END SUBROUTINE thincurr_cross_eval
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_Lmat(tw_ptr,use_hodlr,Lmat_ptr,cache_file,error_str) BIND(C,NAME="thincurr_Lmat")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
LOGICAL(KIND=c_bool), VALUE, INTENT(in) :: use_hodlr
TYPE(c_ptr), INTENT(out) :: Lmat_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: cache_file(80) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(200) !< Needs docs
!
CHARACTER(LEN=OFT_PATH_SLEN) :: filename = ''
TYPE(tw_type), POINTER :: tw_obj
TYPE(oft_tw_hodlr_op), POINTER :: hodlr_op
CALL copy_string('',error_str)
CALL c_f_pointer(tw_ptr, tw_obj)
!
IF(use_hodlr)THEN
  ALLOCATE(hodlr_op)
  hodlr_op%tw_obj=>tw_obj
  CALL hodlr_op%setup(.TRUE.)
  CALL hodlr_op%compute_L()
  Lmat_ptr=C_LOC(hodlr_op)
ELSE
  CALL copy_string_rev(cache_file,filename)
  IF(TRIM(filename)=='')THEN
    CALL tw_compute_LmatDirect(tw_obj,tw_obj%Lmat)
  ELSE
    CALL tw_compute_LmatDirect(tw_obj,tw_obj%Lmat,save_file=filename)
  END IF
  Lmat_ptr=C_LOC(tw_obj%Lmat)
END IF
END SUBROUTINE thincurr_Lmat
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_Mcoil(tw_ptr,Mc_ptr,cache_file,error_str) BIND(C,NAME="thincurr_Mcoil")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
TYPE(c_ptr), INTENT(out) :: Mc_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: cache_file(80) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(200) !< Needs docs
!
CHARACTER(LEN=OFT_PATH_SLEN) :: filename = ''
TYPE(tw_type), POINTER :: tw_obj
CALL copy_string('',error_str)
CALL c_f_pointer(tw_ptr, tw_obj)
!
CALL copy_string_rev(cache_file,filename)
IF(TRIM(filename)=='')THEN
  CALL tw_compute_Ael2dr(tw_obj)
ELSE
  CALL tw_compute_Ael2dr(tw_obj,save_file=filename)
END IF
Mc_ptr=C_LOC(tw_obj%Ael2dr)
END SUBROUTINE thincurr_Mcoil
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_Msensor(tw_ptr,sensor_file,Ms_ptr,Msc_ptr,nsensors,cache_file,error_str) BIND(C,NAME="thincurr_Msensor")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: sensor_file(80) !< Needs docs
TYPE(c_ptr), INTENT(out) :: Ms_ptr !< Needs docs
TYPE(c_ptr), INTENT(out) :: Msc_ptr !< Needs docs
INTEGER(KIND=c_int), INTENT(out) :: nsensors
CHARACTER(KIND=c_char), INTENT(out) :: cache_file(80) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(200) !< Needs docs
!
CHARACTER(LEN=OFT_PATH_SLEN) :: sensor_filename = ''
CHARACTER(LEN=OFT_PATH_SLEN) :: filename = ''
TYPE(tw_type), POINTER :: tw_obj
TYPE(tw_sensors) :: sensors
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: jumper_nsets
CALL copy_string('',error_str)
CALL c_f_pointer(tw_ptr, tw_obj)
!
CALL copy_string_rev(sensor_file,sensor_filename)
NULLIFY(jumper_nsets)
CALL tw_load_sensors(TRIM(sensor_filename),tw_obj,sensors,jumper_nsets)
!
CALL copy_string_rev(cache_file,filename)
IF(TRIM(filename)=='')THEN
  CALL tw_compute_mutuals(tw_obj,sensors%nfloops,sensors%floops)
ELSE
  CALL tw_compute_mutuals(tw_obj,sensors%nfloops,sensors%floops,save_file=filename)
END IF
Ms_ptr=C_LOC(tw_obj%Ael2sen)
Msc_ptr=C_LOC(tw_obj%Adr2sen)
nsensors=sensors%nfloops
END SUBROUTINE thincurr_Msensor
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_Rmat(tw_ptr,copy_out,Rmat,error_str) BIND(C,NAME="thincurr_Rmat")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
LOGICAL(KIND=c_bool), VALUE,  INTENT(in) :: copy_out !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: Rmat !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(200) !< Needs docs
!
INTEGER(4) :: i,j
REAL(8), CONTIGUOUS, POINTER, DIMENSION(:,:) :: Rmat_tmp
TYPE(tw_type), POINTER :: tw_obj
CALL copy_string('',error_str)
CALL c_f_pointer(tw_ptr, tw_obj)
!
CALL tw_compute_Rmat(tw_obj,.TRUE.)
!
IF(copy_out)THEN
  CALL c_f_pointer(Rmat, Rmat_tmp, [tw_obj%nelems,tw_obj%nelems])
  Rmat_tmp=0.d0
  DO i=1,tw_obj%Rmat%nr
    DO j=tw_obj%Rmat%kr(i),tw_obj%Rmat%kr(i+1)-1
      Rmat_tmp(i,tw_obj%Rmat%lc(j))=tw_obj%Rmat%M(j)
    END DO
  END DO
END IF
END SUBROUTINE thincurr_Rmat
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_curr_regmat(tw_ptr,Rmat,error_str) BIND(C,NAME="thincurr_curr_regmat")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: Rmat !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(200) !< Needs docs
!
INTEGER(4) :: i,j,k,jj,pt,ih,ihp,ihc
REAL(8) :: rcurr(3),ftmp(3),gop(3,3),area,norm(3)
REAL(8), POINTER, DIMENSION(:,:) :: Rmat_tmp
TYPE(tw_type), POINTER :: tw_obj
CALL copy_string('',error_str)
CALL c_f_pointer(tw_ptr, tw_obj)
CALL c_f_pointer(Rmat, Rmat_tmp, [3*tw_obj%mesh%nc,tw_obj%nelems])
!
Rmat_tmp=0.d0
!---Avg to cells
ftmp=1.d0/3.d0
!$omp parallel do private(j,pt,rcurr,norm,gop,area)
DO i=1,tw_obj%mesh%nc
  CALL tw_obj%mesh%jacobian(i,ftmp,gop,area)
  CALL tw_obj%mesh%norm(i,ftmp,norm)
  DO j=1,3
    pt=tw_obj%pmap(tw_obj%mesh%lc(j,i))
    IF(pt==0)CYCLE
    rcurr = cross_product(gop(:,j),norm)
    Rmat_tmp((i-1)*3+1,pt) = Rmat_tmp((i-1)*3+1,pt) + rcurr(1)*SQRT(tw_obj%mesh%ca(i))
    Rmat_tmp((i-1)*3+2,pt) = Rmat_tmp((i-1)*3+2,pt) + rcurr(2)*SQRT(tw_obj%mesh%ca(i))
    Rmat_tmp((i-1)*3+3,pt) = Rmat_tmp((i-1)*3+3,pt) + rcurr(3)*SQRT(tw_obj%mesh%ca(i))
  END DO
END DO
!
DO ih=1,tw_obj%nholes
DO ihp=1,tw_obj%hmesh(ih)%n
DO ihc=tw_obj%hmesh(ih)%kpc(ihp),tw_obj%hmesh(ih)%kpc(ihp+1)-1
  i=ABS(tw_obj%hmesh(ih)%lpc(ihc))
  DO j=1,3
    IF(tw_obj%mesh%lc(j,i)==tw_obj%hmesh(ih)%lp(ihp))EXIT
  END DO
  CALL tw_obj%mesh%jacobian(i,ftmp,gop,area)
  CALL tw_obj%mesh%norm(i,ftmp,norm)
  rcurr = cross_product(gop(:,j),norm)*SIGN(1,tw_obj%hmesh(ih)%lpc(ihc))
  pt = tw_obj%np_active+ih
  Rmat_tmp((i-1)*3+1,pt) = Rmat_tmp((i-1)*3+1,pt) + rcurr(1)*SQRT(tw_obj%mesh%ca(i))
  Rmat_tmp((i-1)*3+2,pt) = Rmat_tmp((i-1)*3+2,pt) + rcurr(2)*SQRT(tw_obj%mesh%ca(i))
  Rmat_tmp((i-1)*3+3,pt) = Rmat_tmp((i-1)*3+3,pt) + rcurr(3)*SQRT(tw_obj%mesh%ca(i))
END DO
END DO
END DO
END SUBROUTINE thincurr_curr_regmat
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_eigenvalues(tw_ptr,direct,neigs,eig_vals,eig_vec,hodlr_ptr,error_str) BIND(C,NAME="thincurr_eigenvalues")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
LOGICAL(KIND=c_bool), VALUE, INTENT(in) :: direct !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: neigs !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: eig_vals !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: eig_vec !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: hodlr_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(200) !< Needs docs
!
REAL(8), CONTIGUOUS, POINTER :: vals_tmp(:),vec_tmp(:,:)
TYPE(tw_type), POINTER :: tw_obj
TYPE(oft_tw_hodlr_op), POINTER :: hodlr_op
CALL c_f_pointer(tw_ptr, tw_obj)
IF(tw_obj%nelems<=0)THEN
  CALL copy_string('Invalid ThinCurr model, may not be setup yet',error_str)
  RETURN
END IF
IF((.NOT.ASSOCIATED(tw_obj%Lmat)).AND.(.NOT.c_associated(hodlr_ptr)))THEN
  CALL copy_string('Inductance matrix required, but not computed',error_str)
  RETURN
END IF
IF(.NOT.ASSOCIATED(tw_obj%Rmat))THEN
  CALL copy_string('Resistance matrix required, but not computed',error_str)
  RETURN
END IF
CALL copy_string('',error_str)
CALL c_f_pointer(eig_vals, vals_tmp, [neigs])
CALL c_f_pointer(eig_vec, vec_tmp, [tw_obj%nelems,neigs])
!---Run eigenvalue analysis
IF(direct)THEN
  CALL lr_eigenmodes_direct(tw_obj,neigs,vals_tmp,vec_tmp)
ELSE
  IF(c_associated(hodlr_ptr))THEN
    CALL c_f_pointer(hodlr_ptr, hodlr_op)
    CALL lr_eigenmodes_arpack(tw_obj,neigs,vals_tmp,vec_tmp,hodlr_op=hodlr_op)
  ELSE
    CALL lr_eigenmodes_arpack(tw_obj,neigs,vals_tmp,vec_tmp)
  END IF
END IF
END SUBROUTINE thincurr_eigenvalues
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_freq_response(tw_ptr,direct,fr_limit,freq,fr_driver,hodlr_ptr,error_str) BIND(C,NAME="thincurr_freq_response")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
LOGICAL(KIND=c_bool), VALUE, INTENT(in) :: direct !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: fr_limit !< Needs docs
REAL(KIND=c_double), VALUE, INTENT(in) :: freq !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: fr_driver !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: hodlr_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(200) !< Needs docs
!
REAL(8), CONTIGUOUS, POINTER :: driver_tmp(:,:)
TYPE(tw_type), POINTER :: tw_obj
TYPE(oft_tw_hodlr_op), POINTER :: hodlr_op
CALL c_f_pointer(tw_ptr, tw_obj)
IF(tw_obj%nelems<=0)THEN
  CALL copy_string('Invalid ThinCurr model, may not be setup yet',error_str)
  RETURN
END IF
IF(direct.AND.(c_associated(hodlr_ptr)))THEN
  CALL copy_string('"direct=True" not supported with HODLR compression',error_str)
  RETURN
END IF
IF((.NOT.ASSOCIATED(tw_obj%Lmat)).AND.(.NOT.c_associated(hodlr_ptr)))THEN
  CALL copy_string('Inductance matrix required, but not computed',error_str)
  RETURN
END IF
IF(.NOT.ASSOCIATED(tw_obj%Rmat))THEN
  CALL copy_string('Resistance matrix required, but not computed',error_str)
  RETURN
END IF
CALL copy_string('',error_str)
CALL c_f_pointer(fr_driver, driver_tmp, [tw_obj%nelems,2])
!---Run eigenvalue analysis
IF(c_associated(hodlr_ptr))THEN
  CALL c_f_pointer(hodlr_ptr, hodlr_op)
  CALL frequency_response(tw_obj,LOGICAL(direct),fr_limit,freq,driver_tmp,hodlr_op=hodlr_op)
ELSE
  CALL frequency_response(tw_obj,LOGICAL(direct),fr_limit,freq,driver_tmp)
END IF
END SUBROUTINE thincurr_freq_response
END MODULE thincurr_f