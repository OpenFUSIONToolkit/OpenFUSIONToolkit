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
USE oft_io, ONLY: hdf5_create_file, hdf5_field_get_sizes, hdf5_read, hdf5_field_exist, &
  xdmf_plot_file
!--Grid
USE oft_trimesh_type, ONLY: oft_trimesh
USE oft_mesh_native, ONLY: r_mem, lc_mem, reg_mem, native_read_nodesets, native_read_sidesets
!---Linear Algebra
USE oft_la_base, ONLY: oft_vector
!---
USE fem_utils, ONLY: fem_interp
USE thin_wall, ONLY: tw_type, tw_save_pfield, tw_compute_LmatDirect, tw_compute_Rmat, &
  tw_compute_Ael2dr, tw_sensors, tw_compute_mutuals, tw_load_sensors, tw_compute_Lmat_MF, &
  tw_recon_curr, tw_compute_Bops
USE thin_wall_hodlr, ONLY: oft_tw_hodlr_op
USE thin_wall_solvers, ONLY: lr_eigenmodes_arpack, lr_eigenmodes_direct, frequency_response, &
  tw_reduce_model, run_td_sim, plot_td_sim
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
SUBROUTINE thincurr_setup(mesh_file,np,r_loc,nc,lc_loc,reg_loc,pmap_loc,jumper_start_in,tw_ptr,sizes,error_str,xml_ptr) BIND(C,NAME="thincurr_setup")
CHARACTER(KIND=c_char), INTENT(in) :: mesh_file(OFT_PATH_SLEN) !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: r_loc !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: lc_loc !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: reg_loc !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: pmap_loc !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: np !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: nc !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: jumper_start_in !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: sizes !< Needs docs
TYPE(c_ptr), INTENT(out) :: tw_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: xml_ptr !< Needs docs
!
LOGICAL :: success,is_2d
INTEGER(4) :: i,ndims,ierr,jumper_start
integer(i4), allocatable, dimension(:) :: dim_sizes
INTEGER(i4), POINTER, DIMENSION(:) :: sizes_tmp,pmap_tmp,reg_tmp
INTEGER(i4), POINTER, DIMENSION(:,:) :: lc_tmp
REAL(8), POINTER, DIMENSION(:,:) :: r_tmp
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
NULLIFY(mesh_nsets,hole_nsets,mesh_ssets)
IF(np>0)THEN
  ALLOCATE(tw_obj)
  ALLOCATE(oft_trimesh::tw_obj%mesh)
  CALL tw_obj%mesh%setup(-1,.FALSE.)
  !
  tw_obj%mesh%np=np
  ALLOCATE(tw_obj%mesh%r(3,tw_obj%mesh%np))
  CALL c_f_pointer(r_loc, r_tmp, [3,tw_obj%mesh%np])
  tw_obj%mesh%r=r_tmp
  !
  tw_obj%mesh%nc=nc
  ALLOCATE(tw_obj%mesh%lc(3,tw_obj%mesh%nc))
  CALL c_f_pointer(lc_loc, lc_tmp, [3,tw_obj%mesh%nc])
  tw_obj%mesh%lc=lc_tmp
  !
  ALLOCATE(tw_obj%mesh%reg(tw_obj%mesh%nc))
  IF(c_associated(reg_loc))THEN
    CALL c_f_pointer(reg_loc, reg_tmp, [tw_obj%mesh%nc])
    tw_obj%mesh%reg=reg_tmp
  ELSE
    tw_obj%mesh%reg=1
  END IF
ELSE
  !---------------------------------------------------------------------------
  ! Load model from file
  !---------------------------------------------------------------------------
  ALLOCATE(tw_obj)
  ALLOCATE(oft_trimesh::tw_obj%mesh)
  CALL tw_obj%mesh%setup(-1,.FALSE.)
  CALL hdf5_field_get_sizes(TRIM(filename),"mesh/R",ndims,dim_sizes)
  IF(ndims<0)THEN
    CALL copy_string('Point list ("mesh/R") not present in mesh file',error_str)
    RETURN
  END IF
  tw_obj%mesh%np=dim_sizes(2)
  is_2d=(dim_sizes(1)==2)
  DEALLOCATE(dim_sizes)
  CALL hdf5_field_get_sizes(TRIM(filename),"mesh/LC",ndims,dim_sizes)
  IF(ndims<0)THEN
    CALL copy_string('Cell list ("mesh/LC") not present in mesh file',error_str)
    RETURN
  END IF
  tw_obj%mesh%nc=dim_sizes(2)
  DEALLOCATE(dim_sizes)
  ALLOCATE(tw_obj%mesh%r(3,tw_obj%mesh%np))
  IF(is_2d)THEN
    tw_obj%mesh%r=0.d0
    ALLOCATE(r_tmp(2,tw_obj%mesh%np))
    CALL hdf5_read(r_tmp,TRIM(filename),"mesh/R",success)
    IF(.NOT.success)THEN
      CALL copy_string('Error reading point list from mesh file',error_str)
      RETURN
    END IF
    tw_obj%mesh%r(1:2,:)=r_tmp
    DEALLOCATE(r_tmp)
  ELSE
    CALL hdf5_read(tw_obj%mesh%r,TRIM(filename),"mesh/R",success)
    IF(.NOT.success)THEN
      CALL copy_string('Error reading point list from mesh file',error_str)
      RETURN
    END IF
  END IF
  ALLOCATE(tw_obj%mesh%lc(3,tw_obj%mesh%nc))
  CALL hdf5_read(tw_obj%mesh%lc,TRIM(filename),"mesh/LC",success)
  IF(.NOT.success)THEN
    CALL copy_string('Error reading cell list from mesh file',error_str)
    RETURN
  END IF
  ALLOCATE(tw_obj%mesh%reg(tw_obj%mesh%nc))
  IF(hdf5_field_exist(TRIM(filename),'mesh/REG'))THEN
    CALL hdf5_read(tw_obj%mesh%reg,TRIM(filename),"mesh/REG",success)
    IF(.NOT.success)THEN
      CALL copy_string('Error reading region ID from mesh file',error_str)
      RETURN
    END IF
  ELSE
    tw_obj%mesh%reg=1
  END IF
  !
  IF(hdf5_field_exist(TRIM(filename),'thincurr/periodicity/pmap'))THEN
    ALLOCATE(tw_obj%pmap(tw_obj%mesh%np))
    CALL hdf5_read(tw_obj%pmap,TRIM(filename),'thincurr/periodicity/pmap',success)
    IF(.NOT.success)THEN
      CALL copy_string('Error reading periodicity information from mesh file',error_str)
      RETURN
    END IF
  END IF
  !---Read nodesets and define holes and jumpers
  CALL native_read_nodesets(mesh_nsets,native_filename=TRIM(filename))
  jumper_start=jumper_start_in
  IF(jumper_start/=0)THEN
    ndims=SIZE(mesh_nsets)
    IF(ABS(jumper_start)>ndims)THEN
      CALL copy_string('"jumper_start" exceeds number of nodesets in file',error_str)
      RETURN
    END IF
    IF(jumper_start<0)jumper_start=ndims+1+jumper_start
    hole_nsets=>mesh_nsets(1:jumper_start-1)
    ALLOCATE(tw_obj%jumper_nsets(ndims-jumper_start+1))
    DO i=jumper_start,ndims
      tw_obj%jumper_nsets(i-jumper_start+1)%n=mesh_nsets(i)%n
      ALLOCATE(tw_obj%jumper_nsets(i-jumper_start+1)%v(tw_obj%jumper_nsets(i-jumper_start+1)%n))
      tw_obj%jumper_nsets(i-jumper_start+1)%v=mesh_nsets(i)%v
    END DO
  ELSE
    hole_nsets=>mesh_nsets
  END IF
  !---Read sidesets and copy to closures
  CALL native_read_sidesets(mesh_ssets,native_filename=TRIM(filename))
  IF(ASSOCIATED(mesh_ssets))THEN
    IF(mesh_ssets(1)%n>0)THEN
      tw_obj%nclosures=mesh_ssets(1)%n
      ALLOCATE(tw_obj%closures(tw_obj%nclosures))
      tw_obj%closures=mesh_ssets(1)%v
    END IF
    ndims=SIZE(mesh_ssets)
    DO i=1,ndims
      IF(ASSOCIATED(mesh_ssets(i)%v))DEALLOCATE(mesh_ssets(i)%v)
    END DO
    DEALLOCATE(mesh_ssets)
  END IF
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
!---Deallocate nodesets
IF(ASSOCIATED(mesh_nsets))THEN
  ndims=SIZE(mesh_nsets)
  DO i=1,ndims
    IF(ASSOCIATED(mesh_nsets(i)%v))DEALLOCATE(mesh_nsets(i)%v)
  END DO
  DEALLOCATE(mesh_nsets)
END IF
!
tw_ptr=C_LOC(tw_obj)
CALL c_f_pointer(sizes, sizes_tmp, [9])
sizes_tmp=[tw_obj%mesh%np,tw_obj%mesh%ne,tw_obj%mesh%nc,tw_obj%mesh%nreg,tw_obj%np_active, &
  tw_obj%nholes,tw_obj%n_vcoils,tw_obj%nelems,tw_obj%n_icoils]
END SUBROUTINE thincurr_setup
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_setup_io(tw_ptr,basepath,save_debug,error_str) BIND(C,NAME="thincurr_setup_io")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: basepath(OFT_PATH_SLEN) !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: save_debug !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
!
INTEGER(4) :: i,j,k,npts,nedges
INTEGER(4), ALLOCATABLE :: lctmp(:,:)
REAL(r8), ALLOCATABLE :: rtmp(:,:)
CHARACTER(LEN=OFT_PATH_SLEN) :: pathprefix = ''
TYPE(tw_type), POINTER :: tw_obj
CALL copy_string('',error_str)
CALL c_f_pointer(tw_ptr, tw_obj)
CALL copy_string_rev(basepath,pathprefix)
!---Setup I/0
IF(TRIM(pathprefix)/='')THEN
  CALL tw_obj%xdmf%setup('thincurr',pathprefix)
  CALL tw_obj%mesh%setup_io(tw_obj%xdmf,1)
ELSE
  CALL tw_obj%xdmf%setup('thincurr')
  CALL tw_obj%mesh%setup_io(tw_obj%xdmf,1)
END IF
IF(tw_obj%n_vcoils>0)THEN
  npts=0
  nedges=0
  DO i=1,tw_obj%n_vcoils
    DO j=1,tw_obj%vcoils(i)%ncoils
      npts=npts+tw_obj%vcoils(i)%coils(j)%npts
      nedges=nedges+tw_obj%vcoils(i)%coils(j)%npts-1
    END DO
  END DO
  ALLOCATE(rtmp(3,npts),lctmp(2,nedges))
  npts=0
  nedges=0
  DO i=1,tw_obj%n_vcoils
    DO j=1,tw_obj%vcoils(i)%ncoils
      DO k=1,tw_obj%vcoils(i)%coils(j)%npts
        npts=npts+1
        rtmp(:,npts)=tw_obj%vcoils(i)%coils(j)%pts(:,k)
        IF(k>1)THEN
          nedges=nedges+1
          lctmp(:,nedges)=[npts-1,npts]-1
        END IF
      END DO
    END DO
  END DO
  CALL tw_obj%xdmf%add_mesh(10,rtmp,lctmp,'vcoils')
  DEALLOCATE(rtmp,lctmp)
END IF
IF(tw_obj%n_icoils>0)THEN
  npts=0
  nedges=0
  DO i=1,tw_obj%n_icoils
    DO j=1,tw_obj%icoils(i)%ncoils
      npts=npts+tw_obj%icoils(i)%coils(j)%npts
      nedges=nedges+tw_obj%icoils(i)%coils(j)%npts-1
    END DO
  END DO
  ALLOCATE(rtmp(3,npts),lctmp(2,nedges))
  npts=0
  nedges=0
  DO i=1,tw_obj%n_icoils
    DO j=1,tw_obj%icoils(i)%ncoils
      DO k=1,tw_obj%icoils(i)%coils(j)%npts
        npts=npts+1
        rtmp(:,npts)=tw_obj%icoils(i)%coils(j)%pts(:,k)
        IF(k>1)THEN
          nedges=nedges+1
          lctmp(:,nedges)=[npts-1,npts]-1
        END IF
      END DO
    END DO
  END DO
  CALL tw_obj%xdmf%add_mesh(10,rtmp,lctmp,'icoils')
  DEALLOCATE(rtmp,lctmp)
END IF
IF(save_debug)CALL tw_obj%save_debug()
END SUBROUTINE thincurr_setup_io
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_save_field(tw_ptr,vals,fieldname) BIND(C,NAME="thincurr_save_field")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: vals !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: fieldname(OFT_SLEN) !< Needs docs
!
CHARACTER(LEN=OFT_SLEN) :: name_tmp = ''
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
SUBROUTINE thincurr_recon_curr(tw_ptr,vals,curr,format) BIND(C,NAME="thincurr_recon_curr")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: vals !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: curr !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: format !< Needs docs
INTEGER(4) :: i,j
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
REAL(8), POINTER, DIMENSION(:,:) :: ptvec,cellvec
TYPE(tw_type), POINTER :: tw_obj
DEBUG_STACK_PUSH
CALL c_f_pointer(tw_ptr, tw_obj)
CALL c_f_pointer(vals, vals_tmp, [tw_obj%nelems])
IF(format==1)THEN
  CALL c_f_pointer(curr, cellvec, [3,tw_obj%mesh%nc])
ELSE IF(format==2)THEN
  CALL c_f_pointer(curr, ptvec, [3,tw_obj%mesh%np])
  ALLOCATE(cellvec(3,tw_obj%mesh%nc))
END IF
!---Avg to cells
CALL tw_recon_curr(tw_obj,vals_tmp,cellvec)
!---Avg to points
IF(format==2)THEN
  DO i=1,tw_obj%mesh%np
    ptvec(:,i)=0.d0
    DO j=tw_obj%mesh%kpc(i),tw_obj%mesh%kpc(i+1)-1
      ptvec(:,i) = ptvec(:,i) + cellvec(:,tw_obj%mesh%lpc(j))*tw_obj%mesh%ca(tw_obj%mesh%lpc(j))/3.d0
    END DO
    ptvec(:,i) = ptvec(:,i)/tw_obj%mesh%va(i)
  END DO
  DEALLOCATE(cellvec)
END IF
DEBUG_STACK_POP
END SUBROUTINE thincurr_recon_curr
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_recon_field(tw_ptr,pot,coils,field,hodlr_ptr) BIND(C,NAME="thincurr_recon_field")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: pot !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: coils !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: field !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: hodlr_ptr !< Needs docs
INTEGER(4) :: k,j,jj
REAL(8) :: tmp
REAL(8), POINTER, DIMENSION(:) :: pot_tmp,coils_tmp,vtmp
REAL(8), POINTER, DIMENSION(:,:) :: field_tmp
TYPE(tw_type), POINTER :: tw_obj
CLASS(oft_vector), POINTER :: u,Bx,By,Bz
TYPE(oft_tw_hodlr_op), POINTER :: hodlr_op
DEBUG_STACK_PUSH
CALL c_f_pointer(tw_ptr, tw_obj)
CALL c_f_pointer(pot, pot_tmp, [tw_obj%nelems])
CALL c_f_pointer(coils, coils_tmp, [tw_obj%n_icoils])
CALL c_f_pointer(field, field_tmp, [3,tw_obj%mesh%np])
!
IF(c_associated(hodlr_ptr))THEN
  CALL c_f_pointer(hodlr_ptr, hodlr_op)
  IF(.NOT.ASSOCIATED(hodlr_op%aca_B_dense))CALL hodlr_op%compute_B()
  CALL tw_obj%Uloc%new(u)
  CALL tw_obj%Uloc_pts%new(Bx)
  CALL tw_obj%Uloc_pts%new(By)
  CALL tw_obj%Uloc_pts%new(Bz)
  CALL u%restore_local(pot_tmp)
  CALL hodlr_op%apply_bop(u,Bx,By,Bz)
  NULLIFY(vtmp)
  CALL Bx%get_local(vtmp)
  field_tmp(1,:)=vtmp
  CALL By%get_local(vtmp)
  field_tmp(2,:)=vtmp
  CALL Bz%get_local(vtmp)
  field_tmp(3,:)=vtmp
  DEALLOCATE(vtmp)
  CALL u%delete()
  CALL Bx%delete()
  CALL By%delete()
  CALL Bz%delete()
  DEALLOCATE(u,Bx,By,Bz)
ELSE
  !$omp parallel do private(j,jj,tmp)
  DO k=1,tw_obj%mesh%np
    DO jj=1,3
      tmp=0.d0
      !$omp simd reduction(+:tmp)
      DO j=1,tw_obj%nelems
        tmp=tmp+pot_tmp(j)*tw_obj%Bel(j,k,jj)
      END DO
      field_tmp(jj,k)=tmp
    END DO
  END DO
END IF
IF(tw_obj%n_icoils>0)THEN
  !$omp parallel do private(k,tmp) collapse(2)
  DO j=1,tw_obj%n_icoils
    DO jj=1,3
      !$omp simd
      DO k=1,tw_obj%mesh%np
        field_tmp(jj,k)=field_tmp(jj,k)+coils_tmp(j)*tw_obj%Bdr(k,j,jj)
      END DO
    END DO
  END DO
END IF
DEBUG_STACK_POP
END SUBROUTINE thincurr_recon_field
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_save_scalar(tw_ptr,vals,fieldname) BIND(C,NAME="thincurr_save_scalar")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: vals !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: fieldname(OFT_SLEN) !< Needs docs
!
CHARACTER(LEN=OFT_SLEN) :: name_tmp = ''
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
TYPE(tw_type), POINTER :: tw_obj
CALL c_f_pointer(tw_ptr, tw_obj)
CALL copy_string_rev(fieldname,name_tmp)
CALL c_f_pointer(vals, vals_tmp, [tw_obj%mesh%np])
!---Save plot fields
CALL tw_obj%mesh%save_vertex_scalar(vals_tmp,tw_obj%xdmf,TRIM(name_tmp))
END SUBROUTINE thincurr_save_scalar
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_scale_va(tw_ptr,vals,div_flag) BIND(C,NAME="thincurr_scale_va")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: vals !< Needs docs
LOGICAL, VALUE, INTENT(in) :: div_flag !< Needs docs
REAL(8), POINTER, DIMENSION(:) :: vals_tmp
TYPE(tw_type), POINTER :: tw_obj
CALL c_f_pointer(tw_ptr, tw_obj)
CALL c_f_pointer(vals, vals_tmp, [tw_obj%mesh%np])
IF(div_flag)THEN
  vals_tmp = vals_tmp/tw_obj%mesh%va
ELSE
  vals_tmp = vals_tmp*tw_obj%mesh%va
END IF
END SUBROUTINE thincurr_scale_va
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_cross_coupling(tw_ptr1,tw_ptr2,Mmat,cache_file,error_str) BIND(C,NAME="thincurr_cross_coupling")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr1 !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr2 !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: Mmat !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: cache_file(OFT_PATH_SLEN) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
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
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
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
CHARACTER(KIND=c_char), INTENT(in) :: cache_file(OFT_PATH_SLEN) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
!
CHARACTER(LEN=OFT_PATH_SLEN) :: filename = ''
TYPE(tw_type), POINTER :: tw_obj
TYPE(oft_tw_hodlr_op), POINTER :: hodlr_op
CALL c_f_pointer(tw_ptr, tw_obj)
IF((tw_obj%n_vcoils>0).AND.(.NOT.ASSOCIATED(tw_obj%Acoil2coil)))THEN
  CALL copy_string('Coil mutuals required if, # of Vcoils > 0',error_str)
  RETURN
END IF
CALL copy_string('',error_str)
!
CALL copy_string_rev(cache_file,filename)
IF(use_hodlr)THEN
  ALLOCATE(hodlr_op)
  hodlr_op%tw_obj=>tw_obj
  CALL hodlr_op%setup(.TRUE.)
  WRITE(*,*)
  IF(TRIM(filename)=='')THEN
    CALL hodlr_op%compute_L()
  ELSE
    CALL hodlr_op%compute_L(save_file=filename)
  END IF
  Lmat_ptr=C_LOC(hodlr_op)
ELSE
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
SUBROUTINE thincurr_Bmat(tw_ptr,hodlr_ptr,Bmat_ptr,Bdr_ptr,cache_file,error_str) BIND(C,NAME="thincurr_Bmat")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: hodlr_ptr !< Needs docs
TYPE(c_ptr), INTENT(out) :: Bmat_ptr !< Needs docs
TYPE(c_ptr), INTENT(out) :: Bdr_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: cache_file(OFT_PATH_SLEN) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
!
CHARACTER(LEN=OFT_PATH_SLEN) :: filename = ''
TYPE(tw_type), POINTER :: tw_obj
TYPE(oft_tw_hodlr_op), POINTER :: hodlr_op
CALL c_f_pointer(tw_ptr, tw_obj)
CALL copy_string('',error_str)
!
CALL copy_string_rev(cache_file,filename)
IF(c_associated(hodlr_ptr))THEN
  CALL c_f_pointer(hodlr_ptr,hodlr_op)
  IF(TRIM(filename)=='')THEN
    CALL hodlr_op%compute_B()
  ELSE
    CALL hodlr_op%compute_B(save_file=filename)
  END IF
  Bmat_ptr=c_null_ptr
  Bdr_ptr=C_LOC(hodlr_op%Icoil_Bmat)
ELSE
  IF(TRIM(filename)=='')THEN
    CALL tw_compute_Bops(tw_obj)
  ELSE
    CALL tw_compute_Bops(tw_obj,save_file=filename)
  END IF
  Bmat_ptr=C_LOC(tw_obj%Bel)
  Bdr_ptr=C_LOC(tw_obj%Bdr)
END IF
END SUBROUTINE thincurr_Bmat
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_Mcoil(tw_ptr,Mc_ptr,cache_file,error_str) BIND(C,NAME="thincurr_Mcoil")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
TYPE(c_ptr), INTENT(out) :: Mc_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: cache_file(OFT_PATH_SLEN) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
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
SUBROUTINE thincurr_Msensor(tw_ptr,sensor_file,Ms_ptr,Msc_ptr,nsensors,njumpers,sensor_ptr,cache_file,error_str) BIND(C,NAME="thincurr_Msensor")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: sensor_file(OFT_PATH_SLEN) !< Needs docs
TYPE(c_ptr), INTENT(out) :: Ms_ptr !< Needs docs
TYPE(c_ptr), INTENT(out) :: Msc_ptr !< Needs docs
INTEGER(KIND=c_int), INTENT(out) :: nsensors
INTEGER(KIND=c_int), INTENT(out) :: njumpers
TYPE(c_ptr), INTENT(inout) :: sensor_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: cache_file(OFT_PATH_SLEN) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
!
INTEGER(4) :: i
CHARACTER(LEN=OFT_PATH_SLEN) :: sensor_filename = ''
CHARACTER(LEN=OFT_PATH_SLEN) :: filename = ''
TYPE(tw_type), POINTER :: tw_obj
TYPE(tw_sensors), POINTER :: sensors
CALL copy_string('',error_str)
CALL c_f_pointer(tw_ptr, tw_obj)
!---Deallocate existing object
IF(c_associated(sensor_ptr))THEN
  CALL c_f_pointer(sensor_ptr, sensors)
  DO i=1,sensors%nfloops
    DEALLOCATE(sensors%floops(i)%r)
  END DO
  DEALLOCATE(sensors%floops)
  DO i=1,sensors%njumpers
    DEALLOCATE(sensors%jumpers(i)%hole_facs,sensors%jumpers(i)%points)
  END DO
  DEALLOCATE(sensors%jumpers)
END IF
ALLOCATE(sensors)
!
CALL copy_string_rev(sensor_file,sensor_filename)
CALL tw_load_sensors(TRIM(sensor_filename),tw_obj,sensors)
!
CALL copy_string_rev(cache_file,filename)
IF(TRIM(filename)=='')THEN
  CALL tw_compute_mutuals(tw_obj,sensors%nfloops,sensors%floops)
ELSE
  CALL tw_compute_mutuals(tw_obj,sensors%nfloops,sensors%floops,save_file=filename)
END IF
Ms_ptr=C_LOC(tw_obj%Ael2sen)
Msc_ptr=C_LOC(tw_obj%Adr2sen)
sensor_ptr=C_LOC(sensors)
nsensors=sensors%nfloops
njumpers=sensors%njumpers
END SUBROUTINE thincurr_Msensor
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_get_sensor_name(sensor_ptr,sensor_ind,sensor_name,error_str) BIND(C,NAME="thincurr_get_sensor_name")
TYPE(c_ptr), VALUE, INTENT(in) :: sensor_ptr !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: sensor_ind !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: sensor_name(40) !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(200) !< Needs docs
TYPE(tw_sensors), POINTER :: sensors
CALL copy_string('',error_str)
CALL c_f_pointer(sensor_ptr, sensors)
CALL copy_string(sensors%floops(sensor_ind)%name,sensor_name)
END SUBROUTINE thincurr_get_sensor_name
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_get_eta(tw_ptr,eta_ptr,error_str)BIND(C,NAME="thincurr_get_eta")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: eta_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
!
INTEGER(4) :: nreg_mesh
REAL(8), POINTER :: res_tmp(:)
TYPE(tw_type), POINTER :: tw_obj
CALL copy_string('',error_str)
CALL c_f_pointer(tw_ptr, tw_obj)
!
CALL c_f_pointer(eta_ptr, res_tmp, [tw_obj%mesh%nreg])
res_tmp=tw_obj%Eta_reg*mu0
END SUBROUTINE thincurr_get_eta
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_set_eta(tw_ptr,eta_ptr,error_str)BIND(C,NAME="thincurr_set_eta")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: eta_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
!
INTEGER(4) :: nreg_mesh
REAL(8), POINTER :: res_tmp(:)
TYPE(tw_type), POINTER :: tw_obj
CALL copy_string('',error_str)
CALL c_f_pointer(tw_ptr, tw_obj)
!
CALL c_f_pointer(eta_ptr, res_tmp, [tw_obj%mesh%nreg])
tw_obj%Eta_reg=res_tmp/mu0
END SUBROUTINE thincurr_set_eta
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_Rmat(tw_ptr,copy_out,Rmat,error_str) BIND(C,NAME="thincurr_Rmat")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
LOGICAL(KIND=c_bool), VALUE,  INTENT(in) :: copy_out !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: Rmat !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
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
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
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
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
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
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
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
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_time_domain(tw_ptr,direct,dt,nsteps,cg_tol,timestep_cn,nstatus,nplot, &
  vec_ic,sensor_ptr,ncurr,curr_ptr,nvolt,volt_ptr,volts_full,sensor_vals_ptr,hodlr_ptr,error_str) BIND(C,NAME="thincurr_time_domain")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
LOGICAL(KIND=c_bool), VALUE, INTENT(in) :: direct !< Needs docs
REAL(KIND=c_double), VALUE, INTENT(in) :: dt !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: nsteps !< Needs docs
REAL(KIND=c_double), VALUE, INTENT(in) :: cg_tol !< Needs docs
LOGICAL(KIND=c_bool), VALUE, INTENT(in) :: timestep_cn !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: nstatus !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: nplot !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: vec_ic !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: sensor_ptr !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: ncurr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: curr_ptr !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: nvolt !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: volt_ptr !< Needs docs
LOGICAL(KIND=c_bool), VALUE, INTENT(in) :: volts_full !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: sensor_vals_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: hodlr_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
!
LOGICAL :: pm_save
REAL(8), CONTIGUOUS, POINTER :: ic_tmp(:),curr_waveform(:,:),volt_waveform(:,:),sensor_waveform(:,:)
TYPE(tw_type), POINTER :: tw_obj
TYPE(tw_sensors), POINTER :: sensors
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
IF(c_associated(sensor_ptr))THEN
  CALL c_f_pointer(sensor_ptr, sensors)
ELSE
  ALLOCATE(sensors)
END IF
IF(ncurr>0)THEN
  CALL c_f_pointer(curr_ptr, curr_waveform, [ncurr,tw_obj%n_icoils+1])
ELSE
  NULLIFY(curr_waveform)
END IF
IF(nvolt>0)THEN
  IF(volts_full)THEN
    CALL c_f_pointer(volt_ptr, volt_waveform, [nvolt,tw_obj%nelems+1])
  ELSE
    CALL c_f_pointer(volt_ptr, volt_waveform, [nvolt,tw_obj%n_vcoils+1])
  END IF
ELSE
  NULLIFY(volt_waveform)
END IF
IF(c_associated(sensor_vals_ptr))THEN
  IF(.NOT.c_associated(sensor_ptr))THEN
    CALL copy_string('Sensor object required with sensor waveform',error_str)
    RETURN
  END IF
  IF(.NOT.ASSOCIATED(volt_waveform))THEN
    CALL copy_string('Voltage waveform required with sensor waveform',error_str)
    RETURN
  END IF
  CALL c_f_pointer(sensor_vals_ptr, sensor_waveform, [nvolt,sensors%nfloops+1])
ELSE
  NULLIFY(sensor_waveform)
END IF
CALL c_f_pointer(vec_ic, ic_tmp, [tw_obj%nelems])
!---Run eigenvalue analysis
pm_save=oft_env%pm; oft_env%pm=.FALSE.
IF(c_associated(hodlr_ptr))THEN
  CALL c_f_pointer(hodlr_ptr, hodlr_op)
  CALL run_td_sim(tw_obj,dt,nsteps,ic_tmp,LOGICAL(direct),cg_tol,LOGICAL(timestep_cn), &
    nstatus,nplot,sensors,curr_waveform,volt_waveform,sensor_waveform,hodlr_op=hodlr_op)
ELSE
  CALL run_td_sim(tw_obj,dt,nsteps,ic_tmp,LOGICAL(direct),cg_tol,LOGICAL(timestep_cn), &
    nstatus,nplot,sensors,curr_waveform,volt_waveform,sensor_waveform)
END IF
oft_env%pm=pm_save
END SUBROUTINE thincurr_time_domain
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_time_domain_plot(tw_ptr,compute_B,rebuild_sensors,nsteps,nplot, &
  sensor_ptr,sensor_vals_ptr,nsensor,hodlr_ptr,error_str) BIND(C,NAME="thincurr_time_domain_plot")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
LOGICAL(KIND=c_bool), VALUE, INTENT(in) :: compute_B !< Needs docs
LOGICAL(KIND=c_bool), VALUE, INTENT(in) :: rebuild_sensors !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: nsteps !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: nplot !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: sensor_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: sensor_vals_ptr !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: nsensor !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: hodlr_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
!
LOGICAL :: pm_save
REAL(8), CONTIGUOUS, POINTER :: sensor_waveform(:,:)
TYPE(tw_type), POINTER :: tw_obj
TYPE(tw_sensors), POINTER :: sensors
TYPE(oft_tw_hodlr_op), POINTER :: hodlr_op
CALL c_f_pointer(tw_ptr, tw_obj)
IF(tw_obj%nelems<=0)THEN
  CALL copy_string('Invalid ThinCurr model, may not be setup yet',error_str)
  RETURN
END IF
CALL copy_string('',error_str)
IF(c_associated(sensor_ptr))THEN
  CALL c_f_pointer(sensor_ptr, sensors)
ELSE
  ALLOCATE(sensors)
END IF
IF(nsensor>0)THEN
  IF(.NOT.c_associated(sensor_ptr))THEN
    CALL copy_string('Sensor object required with sensor waveform',error_str)
    RETURN
  END IF
  CALL c_f_pointer(sensor_vals_ptr, sensor_waveform, [nsensor,sensors%nfloops+1])
ELSE
  NULLIFY(sensor_waveform)
END IF
!---Run eigenvalue analysis
pm_save=oft_env%pm; oft_env%pm=.FALSE.
IF(c_associated(hodlr_ptr))THEN
  CALL c_f_pointer(hodlr_ptr, hodlr_op)
  CALL plot_td_sim(tw_obj,nsteps,nplot,sensors,LOGICAL(compute_B),LOGICAL(rebuild_sensors),sensor_waveform,hodlr_op=hodlr_op)
ELSE
  CALL plot_td_sim(tw_obj,nsteps,nplot,sensors,LOGICAL(compute_B),LOGICAL(rebuild_sensors),sensor_waveform)
END IF
oft_env%pm=pm_save
END SUBROUTINE thincurr_time_domain_plot
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE thincurr_reduce_model(tw_ptr,filename,neigs,eig_vec,compute_B,sensor_ptr,hodlr_ptr,error_str) BIND(C,NAME="thincurr_reduce_model")
TYPE(c_ptr), VALUE, INTENT(in) :: tw_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: filename(OFT_PATH_SLEN) !< Needs docs
INTEGER(KIND=c_int), VALUE, INTENT(in) :: neigs !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: eig_vec !< Needs docs
LOGICAL(KIND=c_bool), VALUE, INTENT(in) :: compute_B !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: sensor_ptr !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: hodlr_ptr !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Needs docs
!
REAL(8), CONTIGUOUS, POINTER :: vals_tmp(:),vec_tmp(:,:)
CHARACTER(LEN=OFT_PATH_SLEN) :: h5_path = 'none'
TYPE(tw_type), POINTER :: tw_obj
TYPE(tw_sensors), POINTER :: sensors
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
IF((tw_obj%n_icoils>0).AND.(.NOT.ASSOCIATED(tw_obj%Ael2dr)))THEN
  CALL copy_string('Coil mutuals required if, # of Icoils > 0',error_str)
  RETURN
END IF
CALL copy_string('',error_str)
IF(c_associated(sensor_ptr))THEN
  CALL c_f_pointer(sensor_ptr, sensors)
ELSE
  ALLOCATE(sensors)
END IF
CALL copy_string_rev(filename,h5_path)
CALL c_f_pointer(eig_vec, vec_tmp, [tw_obj%nelems,neigs])
IF(c_associated(hodlr_ptr))THEN
  CALL c_f_pointer(hodlr_ptr, hodlr_op)
  CALL tw_reduce_model(tw_obj,sensors,neigs,vec_tmp,h5_path,LOGICAL(compute_B),hodlr_op=hodlr_op)
ELSE
  CALL tw_reduce_model(tw_obj,sensors,neigs,vec_tmp,h5_path,LOGICAL(compute_B))
END IF
IF(.NOT.c_associated(sensor_ptr))DEALLOCATE(sensors)
END SUBROUTINE thincurr_reduce_model
END MODULE thincurr_f