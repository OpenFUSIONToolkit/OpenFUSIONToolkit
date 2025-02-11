!---------------------------------------------------------------------------
!> @file oft_base_f.F90
!
!> @defgroup doxy_oft_python Python
!! Python interface for the Open FUSION Toolkit
!
!> Fortran part of Python wrappers for OFT
!!
!! @authors Chris Hansen
!! @date May 2023
!! @ingroup doxy_oft_python
!---------------------------------------------------------------------------
MODULE oft_base_f
USE iso_c_binding, ONLY: c_int, c_double, c_char, c_loc, c_null_char, c_ptr, &
    c_f_pointer, c_bool, c_null_ptr
USE oft_base
USE oft_mesh_type, ONLY: smesh
USE oft_mesh_type, ONLY: mesh
USE oft_mesh_native, ONLY: r_mem, lc_mem, reg_mem
USE multigrid_build, ONLY: multigrid_construct, multigrid_construct_surf
IMPLICIT NONE
#include "local.h"
CONTAINS
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE copy_string(f_string,c_string)
CHARACTER(LEN=*), INTENT(in) :: f_string !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: c_string(*) !< Needs docs
INTEGER(4) :: i
DO i=1,LEN(f_string)
    c_string(i)=f_string(i:i)
END DO
c_string(LEN(f_string)+1)=c_null_char
END SUBROUTINE copy_string
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE copy_string_rev(c_string,f_string)
CHARACTER(KIND=c_char), INTENT(in) :: c_string(*) !< Needs docs
CHARACTER(LEN=*), INTENT(inout) :: f_string !< Needs docs
INTEGER(4) :: i
f_string=''
DO i=1,LEN(f_string)
    IF(c_string(i)==c_null_char)EXIT
    f_string(i:i)=c_string(i)
END DO
IF(i>LEN(f_string))CALL oft_warn("No termination character found when copying C string")
END SUBROUTINE copy_string_rev
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE oftpy_init(nthreads,input_file,slens) BIND(C,NAME="oftpy_init")
INTEGER(c_int), VALUE, INTENT(in) :: nthreads !< Needs docs
CHARACTER(KIND=c_char), INTENT(in) :: input_file(OFT_PATH_SLEN) !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: slens !< String lengths
INTEGER(4), POINTER, DIMENSION(:) :: slens_tmp
IF(oft_env%ifile/='none')RETURN
CALL copy_string_rev(input_file,oft_env%ifile)
CALL oft_init(nthreads)
CALL c_f_pointer(slens, slens_tmp, [4])
slens_tmp=[OFT_MPI_PLEN,OFT_SLEN,OFT_PATH_SLEN,OFT_ERROR_SLEN]
END SUBROUTINE oftpy_init
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE oftpy_load_xml(xml_file,oft_node_ptr) BIND(C,NAME="oftpy_load_xml")
CHARACTER(KIND=c_char), INTENT(in) :: xml_file(OFT_PATH_SLEN) !< Needs docs
TYPE(c_ptr), INTENT(out) :: oft_node_ptr !< Needs docs
#ifdef HAVE_XML
INTEGER(i4) :: ierr
LOGICAL :: rst
CHARACTER(LEN=OFT_PATH_SLEN) :: xml_filename = 'none'
TYPE(fox_node), POINTER :: doc,oft_node
!---Test for existence of XML file
CALL copy_string_rev(xml_file,xml_filename)
INQUIRE(FILE=TRIM(xml_filename),exist=rst)
IF(.NOT.rst)RETURN
doc=>fox_parseFile(TRIM(xml_filename),iostat=ierr)
IF(ierr/=0)RETURN
CALL xml_get_element(doc,"oft",oft_node,ierr)
IF(ierr/=0)RETURN
oft_node_ptr=C_LOC(oft_node)
#else
oft_node_ptr=C_NULL_PTR
#endif
END SUBROUTINE oftpy_load_xml
!------------------------------------------------------------------------------
!> Set debug verbosity level
!------------------------------------------------------------------------------
SUBROUTINE oftpy_set_debug(debug_level) BIND(C,NAME="oftpy_set_debug")
INTEGER(c_int), VALUE, INTENT(in) :: debug_level !< New value for debug level (must be in range [0,3])
oft_env%debug=debug_level
END SUBROUTINE oftpy_set_debug
!------------------------------------------------------------------------------
!> Set the number of OpenMP threads to use
!------------------------------------------------------------------------------
SUBROUTINE oftpy_set_nthreads(nthreads) BIND(C,NAME="oftpy_set_nthreads")
INTEGER(c_int), VALUE, INTENT(in) :: nthreads !< Number of threads to use for subsequent OpenMP parallel regions
CHARACTER(LEN=4) :: thrd_str,proc_str
IF(nthreads>omp_get_num_procs())THEN
    WRITE(thrd_str,'(I4.4)')nthreads
    WRITE(proc_str,'(I4.4)')omp_get_num_procs()
    CALL oft_warn("Number of requested threads ("//thrd_str//") exceeds number of available processors ("//proc_str//")")
END IF
CALL omp_set_num_threads(nthreads)
END SUBROUTINE oftpy_set_nthreads
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE oft_setup_smesh(ndim,np,r_loc,npc,nc,lc_loc,reg_loc,nregs) BIND(C,NAME="oft_setup_smesh")
TYPE(c_ptr), VALUE, INTENT(in) :: r_loc !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: lc_loc !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: reg_loc !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: np !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: ndim !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: npc !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: nc !< Needs docs
INTEGER(c_int), INTENT(out) :: nregs !< Needs docs
integer(i4), POINTER :: lc_tmp(:,:),reg_tmp(:)
real(r8), POINTER :: r_tmp(:,:)
IF(ndim>0)THEN
    ALLOCATE(r_mem(ndim,np))
    CALL c_f_pointer(r_loc, r_tmp, [ndim,np])
    r_mem=r_tmp
    ALLOCATE(lc_mem(npc,nc))
    CALL c_f_pointer(lc_loc, lc_tmp, [npc,nc])
    lc_mem=lc_tmp
    ALLOCATE(reg_mem(nc))
    CALL c_f_pointer(reg_loc, reg_tmp, [nc])
    reg_mem=reg_tmp
END IF
!---Setup Mesh
CALL multigrid_construct_surf
nregs=smesh%nreg
END SUBROUTINE oft_setup_smesh
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE oft_setup_vmesh(ndim,np,r_loc,npc,nc,lc_loc,reg_loc,nregs) BIND(C,NAME="oft_setup_vmesh")
TYPE(c_ptr), VALUE, INTENT(in) :: r_loc !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: lc_loc !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: reg_loc !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: np !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: ndim !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: npc !< Needs docs
INTEGER(c_int), VALUE, INTENT(in) :: nc !< Needs docs
INTEGER(c_int), INTENT(out) :: nregs !< Needs docs
integer(i4), POINTER :: lc_tmp(:,:),reg_tmp(:)
real(r8), POINTER :: r_tmp(:,:)
IF(ndim>0)THEN
    ALLOCATE(r_mem(ndim,np))
    CALL c_f_pointer(r_loc, r_tmp, [ndim,np])
    r_mem=r_tmp
    ALLOCATE(lc_mem(npc,nc))
    CALL c_f_pointer(lc_loc, lc_tmp, [npc,nc])
    lc_mem=lc_tmp
    ALLOCATE(reg_mem(nc))
    CALL c_f_pointer(reg_loc, reg_tmp, [nc])
    reg_mem=reg_tmp
END IF
!---Setup Mesh
CALL multigrid_construct
nregs=mesh%nreg
END SUBROUTINE oft_setup_vmesh
END MODULE oft_base_f