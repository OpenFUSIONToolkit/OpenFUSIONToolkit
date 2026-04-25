!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------
MODULE oft_base_f
USE iso_c_binding, ONLY: c_int, c_double, c_char, c_loc, c_null_char, c_ptr, &
    c_f_pointer, c_bool, c_null_ptr, c_funptr, c_associated, c_f_procpointer
USE oft_base
USE oft_mesh_native, ONLY: r_mem, lc_mem, reg_mem
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct, multigrid_construct_surf
IMPLICIT NONE
#include "local.h"
CONTAINS
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE copy_string(f_string,c_string)
CHARACTER(LEN=*), INTENT(in) :: f_string !< Needs docs
CHARACTER(KIND=c_char), INTENT(out) :: c_string(*) !< Needs docs
INTEGER(4) :: i
DO i=1,LEN(f_string)
    c_string(i)=f_string(i:i)
END DO
c_string(LEN(f_string)+1)=c_null_char
END SUBROUTINE copy_string
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE oftpy_init(nthreads,quiet,input_file,slens,abort_fun) BIND(C,NAME="oftpy_init")
INTEGER(c_int), VALUE, INTENT(in) :: nthreads !< Needs docs
LOGICAL(c_bool), VALUE, INTENT(in) :: quiet !< If `True`, do not print OFT environment information on initialization
CHARACTER(KIND=c_char), INTENT(in) :: input_file(OFT_PATH_SLEN) !< Needs docs
TYPE(c_ptr), VALUE, INTENT(in) :: slens !< String lengths
TYPE(c_funptr), VALUE, INTENT(in) :: abort_fun !< Abort callback for Python
INTEGER(4), POINTER, DIMENSION(:) :: slens_tmp
LOGICAL :: quiet_f
IF(oft_env%ifile/='none')RETURN
CALL copy_string_rev(input_file,oft_env%ifile)
quiet_f=quiet
CALL oft_init(nthreads=nthreads,quiet=quiet_f)
CALL c_f_pointer(slens, slens_tmp, [4])
IF(c_associated(abort_fun))CALL c_f_procpointer(abort_fun,oft_abort_cb)
slens_tmp=[OFT_MPI_PLEN,OFT_SLEN,OFT_PATH_SLEN,OFT_ERROR_SLEN]
END SUBROUTINE oftpy_init
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE oftpy_load_xml(xml_file,oft_node_ptr) BIND(C,NAME="oftpy_load_xml")
CHARACTER(KIND=c_char), INTENT(in) :: xml_file(OFT_PATH_SLEN) !< Needs docs
TYPE(c_ptr), INTENT(out) :: oft_node_ptr !< Needs docs
INTEGER(i4) :: ierr
LOGICAL :: rst
CHARACTER(LEN=OFT_PATH_SLEN) :: xml_filename = 'none'
TYPE(xml_doc) :: doc
TYPE(xml_node), POINTER :: oft_node
!---Test for existence of XML file
CALL copy_string_rev(xml_file,xml_filename)
INQUIRE(FILE=TRIM(xml_filename),exist=rst)
IF(.NOT.rst)RETURN
CALL xml_parsefile(TRIM(xml_filename),doc,ierr)
IF(ierr/=0)RETURN
oft_node_ptr=C_LOC(doc%root)
END SUBROUTINE oftpy_load_xml
!---------------------------------------------------------------------------------
!> Set debug verbosity level
!---------------------------------------------------------------------------------
SUBROUTINE oftpy_set_debug(debug_level) BIND(C,NAME="oftpy_set_debug")
INTEGER(c_int), VALUE, INTENT(in) :: debug_level !< New value for debug level (must be in range [0,3])
oft_env%debug=debug_level
END SUBROUTINE oftpy_set_debug
!---------------------------------------------------------------------------------
!> Set the number of OpenMP threads to use
!---------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------
!> Setup surface mesh from specified coordinate, connectivity, and region arrays
!---------------------------------------------------------------------------------
SUBROUTINE oft_setup_smesh(ndim,np,r_loc,npc,nc,lc_loc,reg_loc,nregs,mesh_ptr) BIND(C,NAME="oft_setup_smesh")
TYPE(c_ptr), VALUE, INTENT(in) :: r_loc !< Pointer to point locations array `(3,np)`
TYPE(c_ptr), VALUE, INTENT(in) :: lc_loc !< Pointer to cell connectivity array `(npc,nc)`
TYPE(c_ptr), VALUE, INTENT(in) :: reg_loc !< Pointer to cell region array `(nc)`
INTEGER(c_int), VALUE, INTENT(in) :: np !< Number of points in surface mesh
INTEGER(c_int), VALUE, INTENT(in) :: ndim !< Spatial dimension of the mesh coordinates (Note: `r` array is always `(3,np)`)
INTEGER(c_int), VALUE, INTENT(in) :: npc !< Number of points per cell (3 for triangles, 4 for quads)
INTEGER(c_int), VALUE, INTENT(in) :: nc !< Number of cells in surface mesh
INTEGER(c_int), INTENT(out) :: nregs !< Number of regions in surface mesh
TYPE(c_ptr), INTENT(out) :: mesh_ptr !< Pointer to `multigrid_mesh` object containing the constructed surface mesh
integer(i4), POINTER :: lc_tmp(:,:),reg_tmp(:)
real(r8), POINTER :: r_tmp(:,:)
TYPE(multigrid_mesh), POINTER :: mg_mesh
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
ALLOCATE(mg_mesh)
CALL multigrid_construct_surf(mg_mesh)
nregs=mg_mesh%smesh%nreg
mesh_ptr=C_LOC(mg_mesh)
END SUBROUTINE oft_setup_smesh
!---------------------------------------------------------------------------------
!> Return references to arrays defining an existing surface mesh.
!!
!! This routine exposes pointers to the coordinate, connectivity, and region
!! arrays stored inside the Fortran `multigrid_mesh` object referenced by
!! `mesh_ptr`. No data are copied and callers must not free or reallocate the objects.
!---------------------------------------------------------------------------------
SUBROUTINE oft_smesh_get(mesh_ptr,ndim,np,r_loc,npc,nc,lc_loc,reg_loc,nreg,error_str) BIND(C,NAME="oft_smesh_get")
TYPE(c_ptr), VALUE, INTENT(in) :: mesh_ptr !< Pointer to `multigrid_mesh` object containing the desired surface mesh
INTEGER(c_int), INTENT(out) :: ndim !< Spatial dimension of the mesh coordinates (Note: `r` array is always `(3,np)`)
INTEGER(c_int), INTENT(out) :: np !< Number of points in surface mesh
TYPE(c_ptr), INTENT(out) :: r_loc !< Pointer to point locations array `(3,np)`
INTEGER(c_int), INTENT(out) :: npc !< Number of points per cell (3 for triangles, 4 for quads)
INTEGER(c_int), INTENT(out) :: nc !< Number of cells in surface mesh
TYPE(c_ptr), INTENT(out) :: lc_loc !< Pointer to cell connectivity array `(npc,nc)`
TYPE(c_ptr), INTENT(out) :: reg_loc !< Pointer to cell region array `(nc)`
INTEGER(c_int), INTENT(out) :: nreg !< Number of regions in surface mesh
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error information
TYPE(multigrid_mesh), POINTER :: mg_mesh
IF(.NOT.c_associated(mesh_ptr))THEN
  CALL copy_string('Mesh object not associated',error_str)
  RETURN
END IF
CALL copy_string('',error_str)
CALL c_f_pointer(mesh_ptr,mg_mesh)
!---Read values
ndim=mg_mesh%smesh%dim
np=mg_mesh%smesh%np
npc=mg_mesh%smesh%cell_np
nc=mg_mesh%smesh%nc
nreg=mg_mesh%smesh%nreg
!---Get pointers
r_loc=C_LOC(mg_mesh%smesh%r)
lc_loc=C_LOC(mg_mesh%smesh%lc)
reg_loc=C_LOC(mg_mesh%smesh%reg)
END SUBROUTINE oft_smesh_get
!---------------------------------------------------------------------------------
!> Setup volume mesh from specified coordinate, connectivity, and region arrays
!---------------------------------------------------------------------------------
SUBROUTINE oft_setup_vmesh(np,r_loc,npc,nc,lc_loc,reg_loc,nregs,mesh_ptr) BIND(C,NAME="oft_setup_vmesh")
TYPE(c_ptr), VALUE, INTENT(in) :: r_loc !< Pointer to point locations array `(3,np)`
TYPE(c_ptr), VALUE, INTENT(in) :: lc_loc !< Pointer to cell connectivity array `(npc,nc)`
TYPE(c_ptr), VALUE, INTENT(in) :: reg_loc !< Pointer to cell region array `(nc)`
INTEGER(c_int), VALUE, INTENT(in) :: np !< Number of points in volume mesh
INTEGER(c_int), VALUE, INTENT(in) :: npc !< Number of points per cell (4 for tetrahedra, 8 for hexahedra)
INTEGER(c_int), VALUE, INTENT(in) :: nc !< Number of cells in volume mesh
INTEGER(c_int), INTENT(out) :: nregs !< Number of regions in volume mesh
TYPE(c_ptr), INTENT(out) :: mesh_ptr !< Pointer to `multigrid_mesh` object containing the constructed volume mesh
integer(i4), POINTER :: lc_tmp(:,:),reg_tmp(:)
real(r8), POINTER :: r_tmp(:,:)
TYPE(multigrid_mesh), POINTER :: mg_mesh
IF(np>0)THEN
    ALLOCATE(r_mem(3,np))
    CALL c_f_pointer(r_loc, r_tmp, [3,np])
    r_mem=r_tmp
    ALLOCATE(lc_mem(npc,nc))
    CALL c_f_pointer(lc_loc, lc_tmp, [npc,nc])
    lc_mem=lc_tmp
    ALLOCATE(reg_mem(nc))
    CALL c_f_pointer(reg_loc, reg_tmp, [nc])
    reg_mem=reg_tmp
END IF
!---Setup Mesh
ALLOCATE(mg_mesh)
CALL multigrid_construct(mg_mesh)
nregs=mg_mesh%mesh%nreg
mesh_ptr=C_LOC(mg_mesh)
END SUBROUTINE oft_setup_vmesh
!---------------------------------------------------------------------------------
!> Return references to arrays defining an existing surface mesh.
!!
!! This routine exposes pointers to the coordinate, connectivity, and region
!! arrays stored inside the Fortran `multigrid_mesh` object referenced by
!! `mesh_ptr`. No data are copied and callers must not free or reallocate the objects.
!---------------------------------------------------------------------------------
SUBROUTINE oft_vmesh_get(mesh_ptr,np,r_loc,npc,nc,lc_loc,reg_loc,nreg,error_str) BIND(C,NAME="oft_vmesh_get")
TYPE(c_ptr), VALUE, INTENT(in) :: mesh_ptr !< Pointer to `multigrid_mesh` object containing the desired volume mesh
INTEGER(c_int), INTENT(out) :: np !< Number of points in volume mesh
TYPE(c_ptr), INTENT(out) :: r_loc !< Pointer to point locations array `(3,np)`
INTEGER(c_int), INTENT(out) :: npc !< Number of points per cell (4 for tetrahedra, 8 for hexahedra)
INTEGER(c_int), INTENT(out) :: nc !< Number of cells in volume mesh
TYPE(c_ptr), INTENT(out) :: lc_loc !< Pointer to cell connectivity array `(npc,nc)`
TYPE(c_ptr), INTENT(out) :: reg_loc !< Pointer to cell region array `(nc)`
INTEGER(c_int), INTENT(out) :: nreg !< Number of regions in volume mesh
CHARACTER(KIND=c_char), INTENT(out) :: error_str(OFT_ERROR_SLEN) !< Error information
TYPE(multigrid_mesh), POINTER :: mg_mesh
IF(.NOT.c_associated(mesh_ptr))THEN
  CALL copy_string('Mesh object not associated',error_str)
  RETURN
END IF
CALL copy_string('',error_str)
CALL c_f_pointer(mesh_ptr,mg_mesh)
!---Read values
np=mg_mesh%mesh%np
npc=mg_mesh%mesh%cell_np
nc=mg_mesh%mesh%nc
nreg=mg_mesh%mesh%nreg
!---Get pointers
r_loc=C_LOC(mg_mesh%mesh%r)
lc_loc=C_LOC(mg_mesh%mesh%lc)
reg_loc=C_LOC(mg_mesh%mesh%reg)
END SUBROUTINE oft_vmesh_get
END MODULE oft_base_f