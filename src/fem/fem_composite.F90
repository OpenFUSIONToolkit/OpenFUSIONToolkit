!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!> @file fem_composite.F90
!
!> Classes and infrastructure for composite FE representations
!!
!! @authors Chris Hansen
!! @date April 2014
!! @ingroup doxy_oft_fem
!------------------------------------------------------------------------------
MODULE fem_composite
USE oft_base
USE oft_stitching, ONLY: oft_seam, seam_list, destory_seam
USE oft_io, ONLY: hdf5_rst, hdf5_write, hdf5_read, hdf5_rst_destroy, hdf5_create_file, &
  hdf5_field_exist
!---
USE oft_la_base, ONLY: oft_vector, oft_map, map_list, oft_graph_ptr, &
  oft_matrix, oft_matrix_ptr, oft_local_mat, oft_ml_vecspace
USE oft_native_la, ONLY: oft_native_vector, native_vector_cast, &
  native_vector_slice_push, native_vector_slice_pop
USE oft_la_utils, ONLY: create_vector, combine_matrices, create_matrix, &
  create_identity_graph
!---
USE fem_base, ONLY: oft_fem_ptr, oft_ml_fem_ptr, fem_max_levels, &
  fem_common_linkage, fem_idx_ver, fem_idx_path
IMPLICIT NONE
#include "local.h"
!------------------------------------------------------------------------------
!> Composite FE type
!------------------------------------------------------------------------------
TYPE, PUBLIC :: oft_fem_comp_type
  INTEGER(i4) :: nfields = 0 !< Number of fields in composite representation
  CHARACTER(LEN=4), POINTER, DIMENSION(:) :: field_tags => NULL() !< Character tags for fields
  TYPE(oft_fem_ptr), POINTER, DIMENSION(:) :: fields => NULL() !< Individual representations
  TYPE(oft_seam), POINTER :: linkage => NULL() !< Global linkage information
  TYPE(oft_map), POINTER :: map(:) => NULL() !< Linear algebra mapping
  CLASS(oft_vector), POINTER :: cache_PETSc => NULL() !< PETSc vector cache
  CLASS(oft_vector), POINTER :: cache_native => NULL() !< Native vector cache
CONTAINS
  !> Create vector for FE representation
  PROCEDURE :: vec_create => fem_vec_create
  !> Save vector to HDF5 file
  PROCEDURE :: vec_save => fem_vec_save
  !> Load vector from HDF5 file
  PROCEDURE :: vec_load => fem_vec_load
  !> Create matrix for FE representation
  PROCEDURE :: mat_create => fem_mat_create
  !> Create matrix for FE representation
  PROCEDURE :: mat_setup_local => fem_mat_setup_local
  !> Create matrix for FE representation
  PROCEDURE :: mat_destroy_local => fem_mat_destroy_local
  !> Create matrix for FE representation
  PROCEDURE :: mat_zero_local => fem_mat_zero_local
  !> Create matrix for FE representation
  PROCEDURE :: mat_zero_local_rows => fem_mat_zero_local_rows
  !> Create matrix for FE representation
  PROCEDURE :: mat_add_local => fem_mat_add_local
  !> Destroy object
  PROCEDURE :: delete => fem_comp_delete
END TYPE oft_fem_comp_type
!------------------------------------------------------------------------------
!> Composite FE type pointer container
!------------------------------------------------------------------------------
TYPE, PUBLIC :: oft_fem_comp_ptr
  TYPE(oft_fem_comp_type), POINTER :: fe => NULL() !< Composite finite element object
END TYPE oft_fem_comp_ptr
!------------------------------------------------------------------------------
!> Multi-level composite FE type
!------------------------------------------------------------------------------
TYPE, PUBLIC :: oft_ml_fem_comp_type
  INTEGER(i4) :: nfields = 0 !< Number of fields in composite representation
  INTEGER(i4) :: nlevels = 0 !< Number of FE levels
  INTEGER(i4) :: level = 0 !< Current FE level
  INTEGER(i4) :: abs_level = 0 !< Asoblute FE refinement level
  INTEGER(i4) :: blevel = 0 !< FE base level
  INTEGER(i4) :: minlev = 1 !< Minimum level constructed
  TYPE(oft_fem_comp_ptr) :: levels(fem_max_levels) !< Composite FE objects for each level
  TYPE(oft_fem_comp_type), POINTER :: current_level => NULL() !< Composite FE object for current level
  TYPE(oft_ml_fem_ptr), POINTER, DIMENSION(:) :: ml_fields => NULL() !< Individual ML FE objects for each field
  CHARACTER(LEN=4), POINTER, DIMENSION(:) :: field_tags => NULL() !< Character tags for fields
  TYPE(oft_matrix_ptr) :: interp_matrices(fem_max_levels) !< Level-to-level interpolation matrices
CONTAINS
  !> Setup composite FE structure from sub-field FE definitions
  PROCEDURE :: setup => ml_fem_setup
  !> Destroy object
  PROCEDURE :: delete => ml_fem_comp_delete
  !> Create vector for FE representation
  PROCEDURE :: vec_create => ml_fem_vec_create
  !> Build interpolation operators from sub-fields
  PROCEDURE :: build_interp => ml_fem_build_interp
  !> Set level in ML framework if available
  PROCEDURE :: set_level => ml_fem_set_level
END TYPE oft_ml_fem_comp_type
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, PUBLIC, extends(oft_ml_vecspace) :: oft_ml_fe_comp_vecspace
  class(oft_ml_fem_comp_type), pointer :: ML_FE_rep => NULL() !< ML FE representation
  !> Needs docs
  PROCEDURE(ml_fe_base_pop), POINTER :: base_pop => NULL()
  !> Needs docs
  PROCEDURE(ml_fe_base_push), POINTER :: base_push => NULL()
contains
  !> Needs docs
  PROCEDURE :: vec_create => ml_fe_vecspace_create
  !> Needs docs
  PROCEDURE :: interp => ml_fe_vecspace_interp
  !> Needs docs
  PROCEDURE :: inject => ml_fe_vecspace_inject
end type oft_ml_fe_comp_vecspace
ABSTRACT INTERFACE
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
  SUBROUTINE ml_fe_base_push(self,afine,acors)
    IMPORT oft_ml_fe_comp_vecspace, oft_vector
    CLASS(oft_ml_fe_comp_vecspace), INTENT(inout) :: self
    CLASS(oft_vector), INTENT(inout) :: afine
    CLASS(oft_vector), INTENT(inout) :: acors
  END SUBROUTINE ml_fe_base_push
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
  SUBROUTINE ml_fe_base_pop(self,acors,afine)
    IMPORT oft_ml_fe_comp_vecspace, oft_vector
    CLASS(oft_ml_fe_comp_vecspace), INTENT(inout) :: self
    CLASS(oft_vector), INTENT(inout) :: acors
    CLASS(oft_vector), INTENT(inout) :: afine
  END SUBROUTINE ml_fe_base_pop
END INTERFACE
CONTAINS
!------------------------------------------------------------------------------
!> Create weight vector for FE representation
!------------------------------------------------------------------------------
subroutine fem_vec_create(self,new,cache,native)
class(oft_fem_comp_type), intent(inout) :: self
class(oft_vector), pointer, intent(out) :: new !< Vector to create
logical, optional, intent(in) :: cache !< Allow caching (optional)
logical, optional, intent(in) :: native !< Force native representation (optional)
INTEGER(i4) :: i
TYPE(map_list), POINTER :: maps(:)
TYPE(seam_list), POINTER :: stitches(:)
logical :: do_cache,force_native
DEBUG_STACK_PUSH
do_cache=.TRUE.
force_native=.FALSE.
IF(PRESENT(cache))do_cache=cache
IF(PRESENT(native))force_native=native
!---
IF(use_petsc.AND.(.NOT.force_native))THEN
  IF(ASSOCIATED(self%cache_PETSc))THEN
    CALL self%cache_PETSc%new(new)
  ELSE
    !---Create vector
    ALLOCATE(stitches(self%nfields),maps(self%nfields))
    DO i=1,self%nfields
      stitches(i)%s=>self%fields(i)%fe%linkage
      maps(i)%m=>self%fields(i)%fe%map
    END DO
    CALL create_vector(new,stitches,maps)
    DEALLOCATE(stitches,maps)
    !---
    IF(do_cache)THEN
      CALL new%new(self%cache_PETSc)
      IF(.NOT.ASSOCIATED(self%linkage))self%linkage=>self%cache_PETSc%stitch_info
      IF(.NOT.ASSOCIATED(self%map))self%map=>self%cache_PETSc%map
    END IF
  END IF
ELSE
  IF(ASSOCIATED(self%cache_native))THEN
    CALL self%cache_native%new(new)
  ELSE
    !---Create vector
    ALLOCATE(stitches(self%nfields),maps(self%nfields))
    DO i=1,self%nfields
      stitches(i)%s=>self%fields(i)%fe%linkage
      maps(i)%m=>self%fields(i)%fe%map
    END DO
    CALL create_vector(new,stitches,maps,native=native)
    DEALLOCATE(stitches,maps)
    !---
    IF(do_cache)THEN
      CALL new%new(self%cache_native)
      IF(.NOT.ASSOCIATED(self%linkage))self%linkage=>self%cache_native%stitch_info
      IF(.NOT.ASSOCIATED(self%map))self%map=>self%cache_native%map
    END IF
  END IF
END IF
DEBUG_STACK_POP
end subroutine fem_vec_create
!------------------------------------------------------------------------------
!> Save a Lagrange scalar field to a HDF5 restart file
!------------------------------------------------------------------------------
subroutine fem_vec_save(self,source,filename,path,append)
class(oft_fem_comp_type), intent(inout) :: self
class(oft_vector), target, intent(inout) :: source !< Source field
character(LEN=*), intent(in) :: filename !< Name of destination file
character(LEN=*), intent(in) :: path !< Field label in file
logical, optional, intent(in) :: append !< Append to file instead of creating?
INTEGER(i4) :: i
real(r8), pointer, dimension(:) :: valtmp
class(oft_vector), pointer :: outfield
class(oft_native_vector), pointer :: outvec
type(hdf5_rst) :: rst_info
logical :: do_append
DEBUG_STACK_PUSH
do_append=.FALSE.
IF(PRESENT(append))do_append=append
IF(oft_env%head_proc)THEN
  IF(do_append)THEN
    IF(.NOT.hdf5_field_exist(filename,fem_idx_path))CALL hdf5_write(fem_idx_ver,filename,fem_idx_path)
  ELSE
    CALL hdf5_create_file(filename)
    CALL hdf5_write(fem_idx_ver,filename,fem_idx_path)
  END IF
END IF
NULLIFY(valtmp)
!---
IF(oft_debug_print(1))WRITE(*,'(6A)')oft_indent,'Writing "',TRIM(path), &
  '" to file "',TRIM(filename),'"'
!---
DO i=1,self%nfields
  IF(oft_debug_print(2))WRITE(*,'(4X,3A)')'Field -> "', &
  TRIM(path)//'_'//TRIM(self%field_tags(i)),'"'
  CALL self%fields(i)%fe%vec_create(outfield,native=.TRUE.)
  IF(.NOT.native_vector_cast(outvec,outfield))CALL oft_abort('Failed to create "outfield".', &
    'fem_vec_save',__FILE__)
  !---
  CALL source%get_local(valtmp,i)
  CALL outfield%restore_local(valtmp)
  CALL native_vector_slice_push(outvec,self%fields(i)%fe%global%le,rst_info)
  CALL hdf5_write(rst_info,filename,TRIM(path)//'_'//TRIM(self%field_tags(i)))
  CALL hdf5_rst_destroy(rst_info)
  !---
  CALL outfield%delete
  DEALLOCATE(outfield,valtmp)
END DO
DEBUG_STACK_POP
end subroutine fem_vec_save
!------------------------------------------------------------------------------
!> Load a Lagrange scalar field from a HDF5 restart file
!------------------------------------------------------------------------------
subroutine fem_vec_load(self,source,filename,path,err_flag)
class(oft_fem_comp_type), intent(inout) :: self
class(oft_vector), target, intent(inout) :: source !< Destination field
character(LEN=*), intent(in) :: filename !< Name of source file
character(LEN=*), intent(in) :: path !< ield path in file
integer(i4), optional, intent(out) :: err_flag !< Error flag
INTEGER(i4) :: i,version
integer(i8), pointer, dimension(:) :: lge
real(r8), pointer, dimension(:) :: valtmp
class(oft_vector), pointer :: infield
class(oft_native_vector), pointer :: invec
type(hdf5_rst) :: rst_info
logical :: success
DEBUG_STACK_PUSH
CALL hdf5_read(version,filename,fem_idx_path,success)
IF(.NOT.success)THEN
  IF(PRESENT(err_flag))THEN
    err_flag=2
    DEBUG_STACK_POP
    RETURN
  ELSE
    CALL oft_abort("OFT_idx_Version not found, legacy files not supported","fem_vec_load",__FILE__)
  END IF
END IF
!---
NULLIFY(valtmp)
IF(oft_debug_print(1))WRITE(*,'(2X,5A)')'Reading "',TRIM(path), &
  '" from file "',TRIM(filename),'"'
DO i=1,self%nfields
  IF(oft_debug_print(2))WRITE(*,'(4X,3A)')'  Field -> "', &
    TRIM(path)//'_'//TRIM(self%field_tags(i)),'"'
  CALL self%fields(i)%fe%vec_create(infield,native=.TRUE.)
  IF(.NOT.native_vector_cast(invec,infield))CALL oft_abort('Failed to create "infield".', &
    'fem_vec_load',__FILE__)
  lge=>self%fields(i)%fe%global%le
  !---
  CALL native_vector_slice_push(invec,lge,rst_info,alloc_only=.TRUE.)
  IF(PRESENT(err_flag))THEN
    CALL hdf5_read(rst_info,filename,TRIM(path)//'_'//TRIM(self%field_tags(i)),success=success)
    IF(.NOT.success)THEN
      err_flag=-i
      CALL infield%delete
      NULLIFY(lge)
      DEALLOCATE(infield)
      CALL hdf5_rst_destroy(rst_info)
      EXIT
    END IF
  ELSE
    CALL hdf5_read(rst_info,filename,TRIM(path)//'_'//TRIM(self%field_tags(i)))
  END IF
  CALL native_vector_slice_pop(invec,lge,rst_info)
  CALL infield%get_local(valtmp)
  CALL source%restore_local(valtmp,i)
  !---
  CALL hdf5_rst_destroy(rst_info)
  CALL infield%delete
  DEALLOCATE(infield,valtmp)
END DO
NULLIFY(lge)
DEBUG_STACK_POP
end subroutine fem_vec_load
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine fem_mat_create(self,new,mask)
CLASS(oft_fem_comp_type), INTENT(inout) :: self
CLASS(oft_matrix), POINTER, INTENT(out) :: new
INTEGER(i4), OPTIONAL, INTENT(in) :: mask(:,:)
INTEGER(i4) :: i,j,k,nknown_graphs
INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: mat_mask,graph_ids
CLASS(oft_vector), POINTER :: tmp_vec
TYPE(oft_graph_ptr), ALLOCATABLE :: graphs(:,:),known_graphs(:)
DEBUG_STACK_PUSH
!---
IF(oft_debug_print(2))WRITE(*,'(2X,A)')'Building composite FE matrix'
ALLOCATE(mat_mask(self%nfields,self%nfields))
mat_mask=1
IF(PRESENT(mask))mat_mask=mask
ALLOCATE(graphs(self%nfields,self%nfields))
!---Populate known graphs
ALLOCATE(known_graphs(self%nfields*self%nfields))
ALLOCATE(graph_ids(2,self%nfields*self%nfields))
graph_ids=0
nknown_graphs=0
DO i=1,self%nfields
  DO k=1,nknown_graphs
    IF(ALL(graph_ids(:,k)==(/self%fields(i)%fe%type,self%fields(i)%fe%type/)))EXIT
  END DO
  IF(k<=nknown_graphs)CYCLE
  !---
  IF(oft_debug_print(3))WRITE(*,'(4X,A,2I4)')'Building graph',i,i
  nknown_graphs=nknown_graphs+1
  graph_ids(:,nknown_graphs)=(/self%fields(i)%fe%type,self%fields(i)%fe%type/)
  ALLOCATE(known_graphs(nknown_graphs)%g)
  known_graphs(nknown_graphs)%g%nr=self%fields(i)%fe%ne
  known_graphs(nknown_graphs)%g%nrg=self%fields(i)%fe%global%ne
  known_graphs(nknown_graphs)%g%nc=self%fields(i)%fe%ne
  known_graphs(nknown_graphs)%g%ncg=self%fields(i)%fe%global%ne
  known_graphs(nknown_graphs)%g%nnz=self%fields(i)%fe%nee
  known_graphs(nknown_graphs)%g%kr=>self%fields(i)%fe%kee
  known_graphs(nknown_graphs)%g%lc=>self%fields(i)%fe%lee
END DO
!---Set graphs
DO i=1,self%nfields
  DO j=1,self%nfields
    IF(mat_mask(i,j)==0)CYCLE
    IF(mat_mask(i,j)==2)THEN
      IF(i/=j)CALL oft_abort('Identity only valid on diagonal.', &
      'fem_mat_create',__FILE__)
      !---Setup identity graph
      CALL self%fields(i)%fe%vec_create(tmp_vec)
      CALL create_identity_graph(graphs(i,j)%g,tmp_vec)
      CALL tmp_vec%delete
      DEALLOCATE(tmp_vec)
      CYCLE
    END IF
    DO k=1,nknown_graphs
      IF(ALL(graph_ids(:,k)==(/self%fields(i)%fe%type, &
      self%fields(j)%fe%type/)))EXIT
    END DO
    IF(k<=nknown_graphs)THEN
      IF(oft_debug_print(3))WRITE(*,'(4X,A,2I4)')'Using known graph ',i,j
      graphs(i,j)%g=>known_graphs(k)%g
    ELSE
      IF(oft_debug_print(3))WRITE(*,'(4X,A,2I4)')'Building graph ',i,j
      nknown_graphs=nknown_graphs+1
      graph_ids(:,nknown_graphs)=(/self%fields(i)%fe%type, &
      self%fields(j)%fe%type/)
      ALLOCATE(known_graphs(nknown_graphs)%g)
      known_graphs(nknown_graphs)%g%nr=self%fields(i)%fe%ne
      known_graphs(nknown_graphs)%g%nrg=self%fields(i)%fe%global%ne
      known_graphs(nknown_graphs)%g%nc=self%fields(j)%fe%ne
      known_graphs(nknown_graphs)%g%ncg=self%fields(j)%fe%global%ne
      CALL fem_common_linkage(self%fields(i)%fe,self%fields(j)%fe, &
        known_graphs(nknown_graphs)%g%nnz,known_graphs(nknown_graphs)%g%kr, &
        known_graphs(nknown_graphs)%g%lc)
      graphs(i,j)%g=>known_graphs(nknown_graphs)%g
    END IF
  END DO
END DO
!---
CALL self%vec_create(tmp_vec)
CALL create_matrix(new,graphs,tmp_vec,tmp_vec)
CALL tmp_vec%delete
DO i=1,nknown_graphs
  DEALLOCATE(known_graphs(i)%g)
END DO
DEALLOCATE(graphs,known_graphs,mat_mask,graph_ids,tmp_vec)
DEBUG_STACK_POP
end subroutine fem_mat_create
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine fem_mat_setup_local(self,mloc,mask)
CLASS(oft_fem_comp_type), INTENT(inout) :: self
TYPE(oft_local_mat), INTENT(inout) :: mloc(:,:)
INTEGER(i4), OPTIONAL, INTENT(in) :: mask(:,:)
INTEGER(i4) :: i,j,k,nknown_graphs
INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: mat_mask,graph_ids,known_graphs
DEBUG_STACK_PUSH
ALLOCATE(mat_mask(self%nfields,self%nfields))
mat_mask=1
IF(PRESENT(mask))mat_mask=mask
!---Populate known graphs
ALLOCATE(known_graphs(2,self%nfields*self%nfields))
ALLOCATE(graph_ids(2,self%nfields*self%nfields))
graph_ids=0
nknown_graphs=0
DO i=1,self%nfields
  DO j=1,self%nfields
    DO k=1,nknown_graphs
      IF(ALL(graph_ids(:,k)==(/self%fields(i)%fe%type,self%fields(j)%fe%type/)))EXIT
    END DO
    IF(k<=nknown_graphs)CYCLE
    nknown_graphs=nknown_graphs+1
    graph_ids(:,nknown_graphs)=(/self%fields(i)%fe%type,self%fields(j)%fe%type/)
    known_graphs(:,nknown_graphs)=(/i,j/)
  END DO
END DO
!---
DO i=1,self%nfields
  DO j=1,self%nfields
    IF(mat_mask(i,j)==1)THEN
      ALLOCATE(mloc(i,j)%m(self%fields(i)%fe%nce,self%fields(j)%fe%nce))
      DO k=1,nknown_graphs
        IF(ALL(graph_ids(:,k)==(/self%fields(i)%fe%type,self%fields(j)%fe%type/)))EXIT
      END DO
      IF(k<=nknown_graphs)THEN
        mloc(i,j)%ind=>mloc(known_graphs(1,k),known_graphs(2,k))%ind
      ELSE
        ALLOCATE(mloc(i,j)%ind(self%fields(i)%fe%nce,self%fields(j)%fe%nce))
      END IF
    END IF
  END DO
END DO
DEALLOCATE(mat_mask,known_graphs,graph_ids)
DEBUG_STACK_POP
end subroutine fem_mat_setup_local
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine fem_mat_destroy_local(self,mloc)
CLASS(oft_fem_comp_type), INTENT(inout) :: self
TYPE(oft_local_mat), INTENT(inout) :: mloc(:,:)
INTEGER(i4) :: i,j
DEBUG_STACK_PUSH
!---
DO i=1,self%nfields
  DO j=1,self%nfields
    IF(ASSOCIATED(mloc(i,j)%m))deallocate(mloc(i,j)%m)
    IF(ASSOCIATED(mloc(i,j)%ind))deallocate(mloc(i,j)%ind)
  END DO
END DO
DEBUG_STACK_POP
end subroutine fem_mat_destroy_local
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine fem_mat_zero_local(self,mloc)
CLASS(oft_fem_comp_type), INTENT(inout) :: self
TYPE(oft_local_mat), INTENT(inout) :: mloc(:,:)
INTEGER(i4) :: i,j
DEBUG_STACK_PUSH
!---
DO i=1,self%nfields
  DO j=1,self%nfields
    IF(ASSOCIATED(mloc(i,j)%m))mloc(i,j)%m=0.d0
    IF(ASSOCIATED(mloc(i,j)%ind))mloc(i,j)%ind(1,1)=0
  END DO
END DO
DEBUG_STACK_POP
end subroutine fem_mat_zero_local
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine fem_mat_zero_local_rows(self,mloc,flag,irow)
CLASS(oft_fem_comp_type), INTENT(inout) :: self
TYPE(oft_local_mat), INTENT(inout) :: mloc(:,:)
LOGICAL, DIMENSION(:), INTENT(IN) :: flag
INTEGER(i4), INTENT(IN) :: irow
INTEGER(i4) :: i,j
DEBUG_STACK_PUSH
!---
DO i=1,self%fields(irow)%fe%nce
  IF(flag(i))THEN
    DO j=1,self%nfields
      IF(ASSOCIATED(mloc(irow,j)%m))mloc(irow,j)%m(i,:)=0.d0
    END DO
  END IF
END DO
DEBUG_STACK_POP
end subroutine fem_mat_zero_local_rows
!------------------------------------------------------------------------------
!> Add local contributions to full matrix in a thread-safe way
!------------------------------------------------------------------------------
subroutine fem_mat_add_local(self,mat,mloc,iloc,tlocks)
CLASS(oft_fem_comp_type), INTENT(inout) :: self
CLASS(oft_matrix), INTENT(inout) :: mat !< Full matrix
TYPE(oft_local_mat), INTENT(in) :: mloc(:,:) !< Local matrix
TYPE(oft_1d_int), INTENT(IN) :: iloc(:) !< Local FE entries
INTEGER(KIND=omp_lock_kind), INTENT(INOUT) :: tlocks(:) !< OpenMP row thread locks
INTEGER(i4) :: i,j
DEBUG_STACK_PUSH
!---Add matrix components
DO i=1,self%nfields
  CALL omp_set_lock(tlocks(i))
  DO j=1,self%nfields
    IF(ASSOCIATED(mloc(i,j)%m))THEN
      IF(ASSOCIATED(mloc(i,j)%ind))THEN
        CALL mat%add_values(iloc(i)%v,iloc(j)%v, &
        mloc(i,j)%m,self%fields(i)%fe%nce,self%fields(j)%fe%nce,i,j,mloc(i,j)%ind)
      ELSE
        CALL mat%add_values(iloc(i)%v,iloc(j)%v, &
        mloc(i,j)%m,self%fields(i)%fe%nce,self%fields(j)%fe%nce,i,j)
      END IF
    END IF
  END DO
  CALL omp_unset_lock(tlocks(i))
END DO
DEBUG_STACK_POP
end subroutine fem_mat_add_local
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine fem_comp_delete(self)
class(oft_fem_comp_type), intent(inout) :: self
! NULLIFY(self%linkage%be,self%linkage%lbe)
IF(ASSOCIATED(self%linkage))THEN
  CALL destory_seam(self%linkage)
  DEALLOCATE(self%linkage)
END IF
!---
IF(ASSOCIATED(self%map))THEN
  DEALLOCATE(self%map)
  ! IF(ASSOCIATED(self%map%gbe))DEALLOCATE(self%map%gbe)
  ! IF(ASSOCIATED(self%map%slice))DEALLOCATE(self%map%slice)
  ! IF(ASSOCIATED(self%map%lge))DEALLOCATE(self%map%lge)
END IF
IF(ASSOCIATED(self%cache_PETSc))THEN
  CALL self%cache_PETSc%delete()
  DEALLOCATE(self%cache_PETSc)
END IF
IF(ASSOCIATED(self%cache_native))THEN
  CALL self%cache_native%delete()
  DEALLOCATE(self%cache_native)
END IF
end subroutine fem_comp_delete
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine ml_fem_setup(self)
class(oft_ml_fem_comp_type), intent(inout) :: self
type(oft_fem_comp_type), pointer :: fe_tmp
INTEGER(i4) :: i,j
DEBUG_STACK_PUSH
DO j=1,self%nlevels
  ALLOCATE(self%levels(j)%fe)
  fe_tmp=>self%levels(j)%fe
  fe_tmp%nfields=self%nfields
  ALLOCATE(fe_tmp%fields(fe_tmp%nfields))
  ALLOCATE(fe_tmp%field_tags(fe_tmp%nfields))
  DO i=1,self%nfields
    fe_tmp%fields(i)%fe=>self%ml_fields(i)%ml%levels(j)%fe
    fe_tmp%field_tags(i)=self%field_tags(i)
  END DO
END DO
!---Set to highest level
self%level=self%nlevels
self%current_level=>self%levels(self%level)%fe
DEBUG_STACK_POP
end subroutine ml_fem_setup
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine ml_fem_comp_delete(self)
class(oft_ml_fem_comp_type), intent(inout) :: self
INTEGER(i4) :: i
DO i=self%minlev,self%nlevels
  IF(ASSOCIATED(self%levels(i)%fe))THEN
    CALL self%levels(i)%fe%delete()
    DEALLOCATE(self%levels(i)%fe)
  END IF
  IF(ASSOCIATED(self%interp_matrices(i)%m))THEN
    CALL self%interp_matrices(i)%m%delete()
    DEALLOCATE(self%interp_matrices(i)%m)
  END IF
END DO
IF(ASSOCIATED(self%field_tags))DEALLOCATE(self%field_tags)
IF(ASSOCIATED(self%ml_fields))DEALLOCATE(self%ml_fields)
NULLIFY(self%current_level)
self%nlevels=0
self%minlev=1
end subroutine ml_fem_comp_delete
!------------------------------------------------------------------------------
!> Create weight vector for FE representation
!------------------------------------------------------------------------------
subroutine ml_fem_vec_create(self,new,level,cache,native)
class(oft_ml_fem_comp_type), intent(inout) :: self
class(oft_vector), pointer, intent(out) :: new !< Vector to create
integer(i4), optional, intent(in) :: level !< FE level for init (optional)
logical, optional, intent(in) :: cache !< Allow caching (optional)
logical, optional, intent(in) :: native !< Force native representation (optional)
logical :: do_cache,force_native
DEBUG_STACK_PUSH
do_cache=.TRUE.
force_native=.FALSE.
IF(PRESENT(cache))do_cache=cache
IF(PRESENT(native))force_native=native
!---Create vector on current or new level
IF(PRESENT(level))THEN
  IF(level>self%nlevels.OR.level<=0)CALL oft_abort('Invalid FE level change requested', &
                                                   'ml_fem_vec_create',__FILE__)
  CALL self%levels(level)%fe%vec_create(new,cache=cache,native=native)
ELSE
  CALL self%current_level%vec_create(new,cache=cache,native=native)
END IF
DEBUG_STACK_POP
end subroutine ml_fem_vec_create
!------------------------------------------------------------------------------
!> Set the current level for a ML Compsite-FE structure
!------------------------------------------------------------------------------
subroutine ml_fem_set_level(self,level,propogate)
class(oft_ml_fem_comp_type), intent(inout) :: self
integer(i4), intent(in) :: level !< Desired level
logical, optional, intent(in) :: propogate
INTEGER(i4) :: i
DEBUG_STACK_PUSH
IF(level>self%nlevels.OR.level<=0)THEN
  CALL oft_abort('Invalid FE level change requested', &
    'ml_fem_set_level',__FILE__)
END IF
!---Update level
self%level=level
self%current_level=>self%levels(self%level)%fe
self%abs_level=self%level
IF(self%level>self%blevel.AND.self%blevel>0)self%abs_level=self%level-1
IF(PRESENT(propogate))THEN
  IF(propogate)THEN
    DO i=1,self%nfields
      CALL self%ml_fields(i)%ml%set_level(level)
    END DO
  END IF
END IF
DEBUG_STACK_POP
end subroutine ml_fem_set_level
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine ml_fem_build_interp(self,minlev)
class(oft_ml_fem_comp_type), intent(inout) :: self
integer(i4), optional, intent(in) :: minlev
class(oft_vector), pointer :: fvec,cvec
type(oft_graph_ptr), POINTER :: graphs(:,:)
type(oft_matrix_ptr), POINTER :: mats(:,:)
INTEGER(i4) :: i,j,levmin
DEBUG_STACK_PUSH
levmin=2
IF(PRESENT(minlev))levmin=minlev+1
OUTER: DO j=levmin,self%nlevels
!------------------------------------------------------------------------------
! Create composite matrix
!------------------------------------------------------------------------------
  CALL self%set_level(j)
  !---Specify child graphs
  ALLOCATE(graphs(self%nfields,self%nfields))
  DO i=1,self%nfields
    IF(.NOT.ASSOCIATED(self%ml_fields(i)%ml%interp_graphs(j)%g))THEN
      DEALLOCATE(graphs)
      CYCLE OUTER
    END IF
    graphs(i,i)%g=>self%ml_fields(i)%ml%interp_graphs(j)%g
  END DO
  !---Get coarse and fine vectors
  CALL self%vec_create(cvec,level=self%level-1)
  CALL self%vec_create(fvec)
  !---Construct matrix
  CALL create_matrix(self%interp_matrices(j)%m,graphs,fvec,cvec)
  DEALLOCATE(graphs)
!------------------------------------------------------------------------------
! Combine child matrices into composite matrix
!------------------------------------------------------------------------------
  !---Specify child matrices
  ALLOCATE(mats(self%nfields,self%nfields))
  DO i=1,self%nfields
    IF(.NOT.ASSOCIATED(self%ml_fields(i)%ml%interp_matrices(j)%m))THEN
      CALL oft_abort('Sub-field matrix not allocated.','ml_fem_build_interp',__FILE__)
    END IF
    mats(i,i)%m=>self%ml_fields(i)%ml%interp_matrices(j)%m
  END DO
  !---Combine matrices
  CALL combine_matrices(mats,self%nfields,self%nfields,self%interp_matrices(j)%m)
  DEALLOCATE(mats)
  CALL self%interp_matrices(j)%m%assemble
  !---Delete temporaries
  CALL cvec%delete
  CALL fvec%delete
  DEALLOCATE(cvec,fvec)
END DO OUTER
DEBUG_STACK_POP
end subroutine ml_fem_build_interp
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine ml_fe_vecspace_create(self,new,level,cache,native)
class(oft_ml_fe_comp_vecspace), intent(inout) :: self
class(oft_vector), pointer, intent(out) :: new
integer(i4), optional, intent(in) :: level
logical, optional, intent(in) :: cache
logical, optional, intent(in) :: native
DEBUG_STACK_PUSH
CALL self%ML_FE_rep%vec_create(new,level=level,cache=cache,native=native)
DEBUG_STACK_POP
end subroutine ml_fe_vecspace_create
!------------------------------------------------------------------------------
!> Interpolate a coarse level Lagrange scalar field to the next finest level
!!
!! @note The global Lagrange level in incremented by one in this subroutine
!------------------------------------------------------------------------------
subroutine ml_fe_vecspace_interp(self,acors,afine)
class(oft_ml_fe_comp_vecspace), intent(inout) :: self
class(oft_vector), intent(inout) :: acors !< Vector to interpolate
class(oft_vector), intent(inout) :: afine !< Fine vector from interpolation
DEBUG_STACK_PUSH
!---Step one level up
call self%ML_FE_rep%set_level(self%ML_FE_rep%level+1)
call afine%set(0.d0)
if(self%ML_FE_rep%level==self%ML_FE_rep%blevel+1)then
  IF(.NOT.ASSOCIATED(self%base_pop))CALL oft_abort("Base transfer not defined","ml_fe_vecspace_interp",__FILE__)
  call self%base_pop(acors,afine)
  DEBUG_STACK_POP
  return
end if
CALL self%ML_FE_rep%interp_matrices(self%ML_FE_rep%level)%m%apply(acors,afine)
DEBUG_STACK_POP
end subroutine ml_fe_vecspace_interp
!------------------------------------------------------------------------------
!> Interpolate a coarse level Lagrange scalar field to the next finest level
!!
!! @note The global Lagrange level in incremented by one in this subroutine
!------------------------------------------------------------------------------
subroutine ml_fe_vecspace_inject(self,afine,acors)
class(oft_ml_fe_comp_vecspace), intent(inout) :: self
class(oft_vector), intent(inout) :: afine !< Fine vector from interpolation
class(oft_vector), intent(inout) :: acors !< Vector to interpolate
DEBUG_STACK_PUSH
! Step down level down
call self%ML_FE_rep%set_level(self%ML_FE_rep%level-1)
call acors%set(0.d0)
if(self%ML_FE_rep%level==self%ML_FE_rep%blevel)then
  IF(.NOT.ASSOCIATED(self%base_push))CALL oft_abort("Base transfer not defined","ml_fe_vecspace_inject",__FILE__)
  call self%base_push(afine,acors)
  DEBUG_STACK_POP
  return
end if
CALL self%ML_FE_rep%interp_matrices(self%ML_FE_rep%level+1)%m%applyT(afine,acors)
DEBUG_STACK_POP
end subroutine ml_fe_vecspace_inject
END MODULE fem_composite
