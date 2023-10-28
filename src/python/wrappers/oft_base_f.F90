!---------------------------------------------------------------------------
!> @file oft_base_f.F90
!
!> @defgroup doxy_oft_python Python
!! Python interface for the OpenFUSIONToolkit
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
CONTAINS
SUBROUTINE copy_string(f_string,c_string)
CHARACTER(LEN=*), INTENT(in) :: f_string
CHARACTER(KIND=c_char), INTENT(out) :: c_string(*)
INTEGER(4) :: i
DO i=1,LEN(f_string)
    c_string(i)=f_string(i:i)
END DO
c_string(LEN(f_string)+1)=c_null_char
END SUBROUTINE copy_string

SUBROUTINE copy_string_rev(c_string,f_string)
CHARACTER(KIND=c_char), INTENT(in) :: c_string(*)
CHARACTER(LEN=*), INTENT(inout) :: f_string
INTEGER(4) :: i
f_string=''
DO i=1,LEN(f_string)
    IF(c_string(i)==c_null_char)EXIT
    f_string(i:i)=c_string(i)
END DO
END SUBROUTINE copy_string_rev
!
SUBROUTINE oftpy_init(nthreads) BIND(C,NAME="oftpy_init")
INTEGER(c_int), VALUE, INTENT(in) :: nthreads
IF(oft_env%nbase>0)RETURN
oft_env%ifile='oftpyin'
CALL oft_init(nthreads)
END SUBROUTINE oftpy_init
!
SUBROUTINE oft_setup_smesh(ndim,np,r_loc,npc,nc,lc_loc,reg_loc,nregs) BIND(C,NAME="oft_setup_smesh")
TYPE(c_ptr), VALUE, INTENT(in) :: r_loc,lc_loc,reg_loc
INTEGER(c_int), VALUE, INTENT(in) :: np,ndim,npc,nc
INTEGER(c_int), INTENT(out) :: nregs
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
!---------------------------------------------------------------------------
! Setup Mesh
!---------------------------------------------------------------------------
CALL multigrid_construct_surf
nregs=smesh%nreg
END SUBROUTINE oft_setup_smesh
!
SUBROUTINE oft_setup_vmesh(ndim,np,r_loc,npc,nc,lc_loc,reg_loc,nregs) BIND(C,NAME="oft_setup_vmesh")
TYPE(c_ptr), VALUE, INTENT(in) :: r_loc,lc_loc,reg_loc
INTEGER(c_int), VALUE, INTENT(in) :: np,ndim,npc,nc
INTEGER(c_int), INTENT(out) :: nregs
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
!---------------------------------------------------------------------------
! Setup Mesh
!---------------------------------------------------------------------------
CALL multigrid_construct
nregs=mesh%nreg
END SUBROUTINE oft_setup_vmesh
END MODULE oft_base_f