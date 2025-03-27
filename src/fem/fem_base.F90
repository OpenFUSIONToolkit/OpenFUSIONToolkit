!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!> @file fem_base.F90
!
!> @defgroup doxy_oft_fem Finite Element
!! Finite element constructions for the Open FUSION Toolkit
!
!> Base FEM class and functions for construction of FE linkage
!! - FEM base class
!! - FE linkage construction
!! - FE cell dof mapping
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_fem
!------------------------------------------------------------------------------
module fem_base
USE oft_base
USE oft_sort, ONLY: sort_array, search_array, sort_matrix
USE oft_quadrature
USE oft_mesh_type, ONLY: oft_mesh, oft_bmesh, oft_init_seam
USE multigrid, ONLY: multigrid_mesh, multigrid_level
USE oft_stitching, ONLY: oft_seam, seam_list, oft_stitch_check, destory_seam
USE oft_io, ONLY: hdf5_rst, hdf5_write, hdf5_read, hdf5_rst_destroy, hdf5_create_file, &
  hdf5_field_exist
!---
USE oft_la_base, ONLY: oft_vector, oft_map, map_list, oft_graph, &
  oft_graph_ptr, oft_matrix, oft_matrix_ptr, oft_ml_vecspace
USE oft_native_la, ONLY: oft_native_vector, native_vector_cast, &
  native_vector_slice_push, native_vector_slice_pop
USE oft_la_utils, ONLY: create_vector, create_matrix
IMPLICIT NONE
#include "local.h"
PRIVATE
!---
integer(i4), public, parameter :: fem_max_levels=10 !< Maximum number of FE levels
integer(i4), public, parameter :: fem_idx_ver=1 !< File version for array indexing
character(LEN=16), public, parameter :: fem_idx_path="OFT_idx_Version" !< HDF5 field name
!------------------------------------------------------------------------------
!> Base FE type
!------------------------------------------------------------------------------
TYPE, ABSTRACT, PUBLIC :: oft_afem_type
  INTEGER(i4) :: order = -1 !< FE rep order
  INTEGER(i4) :: dim = -1 !< Dimension of FE rep (ex. 3 for 3-vector or 1 for scalar)
  INTEGER(i4) :: type = -1 !< FE type
  INTEGER(i4) :: ne = -1 !< Number of total elements
  INTEGER(i4) :: nce = -1 !< Number of elements per cell
  INTEGER(i4) :: nnodes = -1 !< Number of quadrature nodes on mesh
  INTEGER(i4) :: nee = -1 !< Number of element to element interactions
  INTEGER(i4) :: nec = -1 !< Number of element to cell interactions
  INTEGER(i4) :: nbe = -1 !< Number of boundary elements
  INTEGER(i4) :: necmax = -1 !< Maximum number of neighbors cells for one element
  LOGICAL, POINTER, DIMENSION(:) :: be => NULL() !< Boundary element flag
  LOGICAL, POINTER, DIMENSION(:) :: ce => NULL() !< Corner edge flag
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: kee => NULL() !< Pointer to element connectivity list
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:) :: lee => NULL() !< Element to element connectivity list
  INTEGER(i4), POINTER, DIMENSION(:) :: kec => NULL() !< Pointer to cell connectivity list
  INTEGER(i4), POINTER, DIMENSION(:) :: lec => NULL() !< Element to cell connectivity list
  INTEGER(i4), POINTER, DIMENSION(:) :: lbe => NULL() !< List of boundary elements
  INTEGER(i8), POINTER, DIMENSION(:) :: legacy_lge => NULL() !< Legacy global element list
  INTEGER(i4), POINTER, DIMENSION(:) :: bc => NULL() !< Boundary condition type
  TYPE(dof_map), POINTER, DIMENSION(:) :: cmap => NULL() !< Mapping from face index to dof type
  TYPE(oft_quad_type) :: quad !< Global quadrature structure
  TYPE(oft_seam), POINTER :: linkage => NULL() !< Global linkage information
  TYPE(oft_map), POINTER :: map => NULL() !< Linear algebra mapping
  TYPE(fem_mpi_global), POINTER :: global => NULL() !< Global index information
  CLASS(oft_vector), POINTER :: cache_PETSc => NULL() !< PETSc vector cache
  CLASS(oft_vector), POINTER :: cache_native => NULL() !< Native vector cache
CONTAINS
  !> Setup FE representation
  PROCEDURE(afem_setup), DEFERRED :: setup
  !> Get weight indices for a given cell
  PROCEDURE(afem_ncdofs), DEFERRED :: ncdofs
  !> Create vector for FE representation
  PROCEDURE :: vec_create => afem_vec_create
  !> Save vector to HDF5 file
  PROCEDURE :: vec_save => afem_vec_save
  !> Load vector from HDF5 file
  PROCEDURE :: vec_load => afem_vec_load
  !> Create matrix for FE representation
  PROCEDURE :: mat_create => afem_mat_create
  !> Destory FE type
  PROCEDURE :: delete => afem_delete
END TYPE oft_afem_type
!
ABSTRACT INTERFACE
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
  subroutine afem_setup(self,quad_order)
    import oft_afem_type, i4
    class(oft_afem_type), intent(inout) :: self
    integer(i4), intent(in) :: quad_order
  end subroutine afem_setup
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------  
  subroutine afem_ncdofs(self,cell,dofs)
    import oft_afem_type, i4
    class(oft_afem_type), intent(in) :: self
    integer(i4), intent(in) :: cell
    integer(i4), intent(inout) :: dofs(:)
  end subroutine afem_ncdofs
END INTERFACE
!------------------------------------------------------------------------------
!> Base FE type for boundary (triangle) meshes
!------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_afem_type) :: oft_bfem_type
  INTEGER(i4) :: gstruct(3) = -1 !< Geometric mapping array
  CLASS(oft_bmesh), POINTER :: mesh => NULL() !< Structure containing bound mesh
  TYPE(fem_parent), POINTER :: parent => NULL() !< Global index information
CONTAINS
  !> Setup FE representation
  PROCEDURE :: setup => bfem_setup
  !> Get weight indices for a given cell
  PROCEDURE :: ncdofs => bfem_ncdofs
  !> Destory FE type
  PROCEDURE :: delete => bfem_delete
END TYPE oft_bfem_type
!------------------------------------------------------------------------------
!> Base FE type
!------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_afem_type) :: oft_fem_type
  INTEGER(i4) :: gstruct(4) = -1 !< Geometric mapping array
  CLASS(oft_mesh), POINTER :: mesh => NULL() !< Structure containing bound mesh
CONTAINS
  !> Setup FE representation
  PROCEDURE :: setup => fem_setup
  !> Get weight indices for a given cell
  PROCEDURE :: ncdofs => fem_ncdofs
  !> Destory FE type
  PROCEDURE :: delete => fem_delete
END TYPE oft_fem_type
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, PUBLIC :: oft_fem_ptr
  CLASS(oft_afem_type), pointer :: fe => NULL() !< Finite element object
end type oft_fem_ptr
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
TYPE, PUBLIC :: oft_ml_fem_type
  INTEGER(i4) :: nlevels = 0 !< Number of FE levels
  INTEGER(i4) :: level = 0 !< Current FE level
  INTEGER(i4) :: abs_level = 0 !< Asoblute FE refinement level
  INTEGER(i4) :: blevel = 0 !< FE base level
  INTEGER(i4) :: minlev = 1 !< Lowest level
  TYPE(multigrid_mesh), POINTER :: ml_mesh => NULL() !< Structure containing bound ML mesh
  CLASS(oft_afem_type), POINTER :: current_level => NULL()
  TYPE(oft_fem_ptr) :: levels(fem_max_levels)
  TYPE(oft_graph_ptr) :: interp_graphs(fem_max_levels)
  TYPE(oft_matrix_ptr) :: interp_matrices(fem_max_levels)
CONTAINS
  !> Create vector for FE representation
  PROCEDURE :: vec_create => ml_fem_vec_create
  !> Set level in ML framework if available
  PROCEDURE :: set_level => ml_fem_set_level
  !> Destory FE type
  PROCEDURE :: delete => ml_fem_delete
END TYPE oft_ml_fem_type
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, PUBLIC :: oft_ml_fem_ptr
  TYPE(oft_ml_fem_type), pointer :: ml => NULL() !< ML object
end type oft_ml_fem_ptr
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, PUBLIC, extends(oft_ml_vecspace) :: oft_ml_fe_vecspace
  class(oft_ml_fem_type), pointer :: ML_FE_rep => NULL() !< ML FE representation
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
end type oft_ml_fe_vecspace
ABSTRACT INTERFACE
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
  SUBROUTINE ml_fe_base_push(self,afine,acors)
    IMPORT oft_ml_fe_vecspace, oft_vector
    CLASS(oft_ml_fe_vecspace), INTENT(inout) :: self
    CLASS(oft_vector), INTENT(inout) :: afine
    CLASS(oft_vector), INTENT(inout) :: acors
  END SUBROUTINE ml_fe_base_push
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
  SUBROUTINE ml_fe_base_pop(self,acors,afine)
    IMPORT oft_ml_fe_vecspace, oft_vector
    CLASS(oft_ml_fe_vecspace), INTENT(inout) :: self
    CLASS(oft_vector), INTENT(inout) :: acors
    CLASS(oft_vector), INTENT(inout) :: afine
  END SUBROUTINE ml_fe_base_pop
END INTERFACE
!------------------------------------------------------------------------------
!> Cell DOF type information
!------------------------------------------------------------------------------
type :: dof_map
  integer(i4) :: type = 0 !< Geometry type
  integer(i4) :: el = 0 !< Element id
  integer(i4) :: ind = 0 !< Sub-index
end type dof_map
!------------------------------------------------------------------------------
!> Global vector information and indicies
!------------------------------------------------------------------------------
type :: fem_mpi_global
  integer(8) :: ne = 0 !< Global element count
  integer(8), pointer, dimension(:) :: le => NULL() !< Global index of elements (ne)
  logical, pointer, dimension(:) :: gbe => NULL() !< Global boundary element flag (ne)
end type fem_mpi_global
!------------------------------------------------------------------------------
!> Parent FE information and indicies
!------------------------------------------------------------------------------
type :: fem_parent
  integer(i4) :: ne = 0 !< Parent element count
  integer(i4), pointer, dimension(:) :: le => NULL() !< Parent index of elements (ne)
end type fem_parent
!---
PUBLIC fem_common_linkage, fem_delete, bfem_delete
contains
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine afem_delete(self)
class(oft_afem_type), intent(inout) :: self
IF(ASSOCIATED(self%be))DEALLOCATE(self%be)
IF(ASSOCIATED(self%ce))DEALLOCATE(self%ce)
IF(ASSOCIATED(self%kee))DEALLOCATE(self%kee)
IF(ASSOCIATED(self%lee))DEALLOCATE(self%lee)
IF(ASSOCIATED(self%kec))DEALLOCATE(self%kec)
IF(ASSOCIATED(self%lec))DEALLOCATE(self%lec)
IF(ASSOCIATED(self%lbe))DEALLOCATE(self%lbe)
IF(ASSOCIATED(self%legacy_lge))DEALLOCATE(self%legacy_lge)
IF(ASSOCIATED(self%bc))DEALLOCATE(self%bc)
IF(ASSOCIATED(self%cmap))DEALLOCATE(self%cmap)
CALL self%quad%delete()
NULLIFY(self%linkage%be,self%linkage%lbe)
CALL destory_seam(self%linkage)
DEALLOCATE(self%linkage)
!---
IF(ASSOCIATED(self%global))THEN
  IF(ASSOCIATED(self%global%le))DEALLOCATE(self%global%le)
  IF(ASSOCIATED(self%global%gbe))DEALLOCATE(self%global%gbe)
END IF
!---
IF(ASSOCIATED(self%map))THEN
  !IF(ASSOCIATED(self%map%gbe))DEALLOCATE(self%map%gbe)
  IF(ASSOCIATED(self%map%slice))DEALLOCATE(self%map%slice)
  !IF(ASSOCIATED(self%map%lge))DEALLOCATE(self%map%lge)
  !self%map%lge=>self%global%le
  !self%map%gbe=>self%global%gbe
END IF
IF(ASSOCIATED(self%cache_PETSc))THEN
  CALL self%cache_PETSc%delete()
  DEALLOCATE(self%cache_PETSc)
END IF
IF(ASSOCIATED(self%cache_native))THEN
  CALL self%cache_native%delete()
  DEALLOCATE(self%cache_native)
END IF
end subroutine afem_delete
!------------------------------------------------------------------------------
!> Compute element to element linkage for a FE representation
!!
!! Creates a CSR graph representing the interaction between elements in a single
!! finite element representation. This interaction comes from performing a volume
!! integration of test functions against basis functions, which in Galerkin
!! finite element are the same sets. As a result interactions are present for
!! any two elements who share a common cell.
!!
!! @note The graph is constructed in the @ref fem_base::oft_fem_type::nee "nee",
!! @ref fem_base::oft_fem_type::kee "kee", and @ref fem_base::oft_fem_type::lee "lee"
!! fields of `self`
!------------------------------------------------------------------------------
subroutine afem_self_linkage(self)
class(oft_afem_type), intent(inout) :: self !< Finite element representation
integer(i4) :: i,j,k,offset,stack
integer(i4) :: jeg,jsg,jee,jse,js,je,jn
integer(i4) :: mycounts(4)
integer(i4), allocatable :: lcx(:,:),nr(:)
! class(oft_mesh), pointer :: mesh
DEBUG_STACK_PUSH
!---
! mesh=>self%mesh
allocate(nr(self%ne+1))
allocate(self%kee(self%ne+1))
nr=0 ! initialize interaction counter
!$omp parallel private(js,je,jn,j,k,lcx)
allocate(lcx(self%nce,self%necmax+1))
!$omp do
do i=1,self%ne ! loop over dofs
  js=self%kec(i)
  je=self%kec(i+1)-1
  IF(je<js)CYCLE
  do j=js,je  ! loop over neighbor cells
    k=self%lec(j)
    call self%ncdofs(k,lcx(:,j-js+1)) ! get DOFs
  end do
  jn=(je-js+1)
  call sort_array(lcx(:,1:jn),int(self%nce,4),jn)
  !
  nr(i)=1
  jn=lcx(1,1)
  do j=js,je ! loop over neighbor cells
    do k=1,self%nce ! loop over cell dofs
      if(lcx(k,j-js+1)/=jn)then
        nr(i)=nr(i)+1
        jn=lcx(k,j-js+1)
      end if
    end do
  end do
end do
!$omp master
self%nee=sum(nr)
self%kee(self%ne+1)=self%nee+1
do i=self%ne,1,-1 ! cumulative point to point count
  self%kee(i)=self%kee(i+1)-nr(i)
end do
if(self%kee(1)/=1)call oft_abort('Bad element to element count','fem_self_linkage',__FILE__)
allocate(self%lee(self%nee))
nr=0 ! reset interaction counter
!$omp end master
!$omp barrier
!$omp do
do i=1,self%ne ! loop over points
  js=self%kec(i)
  je=self%kec(i+1)-1
  IF(je<js)CYCLE
  do j=js,je  ! loop over neighbor cells
    k=self%lec(j)
    call self%ncdofs(k,lcx(:,j-js+1)) ! get DOFs
  end do
  jn=(je-js+1)
  call sort_array(lcx(:,1:jn),int(self%nce,4),jn)
  !
  nr(i)=1
  jn=lcx(1,1)
  self%lee(self%kee(i))=jn
  do j=js,je ! loop over neighbor cells
    do k=1,self%nce ! loop over corners
      if(lcx(k,j-js+1)/=jn)then
        self%lee(self%kee(i)+nr(i))=lcx(k,j-js+1)
        nr(i)=nr(i)+1
        jn=lcx(k,j-js+1)
      end if
    end do
  end do
end do
deallocate(lcx)
!$omp end parallel
deallocate(nr)
DEBUG_STACK_POP
end subroutine afem_self_linkage
!------------------------------------------------------------------------------
!> Compute element to element linkage between two FE representations
!!
!! Creates a CSR graph representing the interaction between elements of two different
!! finite element representations. This listed is constructed as in
!! @ref fem_base::fem_self_linkage "fem_self_linkage", however in this case
!! `self` is used as the test functions and `other` is used as the basis set
!! in the Galerkin intergral. These correspond to the row and columns of the CSR
!! graph respectively
!------------------------------------------------------------------------------
subroutine fem_common_linkage(self,other,nee,kee,lee)
class(oft_afem_type), intent(in) :: self !< Finite element representation for test set (rows)
class(oft_afem_type), intent(in) :: other !< Finite element representation for basis set (columns)
integer(i4), intent(out) :: nee !< Number of entries in graph
integer(i4), pointer, intent(out) :: kee(:) !< Row pointer into column list [self%ne+1]
integer(i4), pointer, intent(out) :: lee(:) !< Column list [nee]
integer(i4) :: i,j,k,offset,stack
integer(i4) :: jeg,jsg,jee,jse,necmax,js,je,jn
integer(i4) :: mycounts(4)
integer(i4), allocatable :: lcx(:,:),nr(:)
! class(oft_mesh), pointer :: mesh
DEBUG_STACK_PUSH
!---
! mesh=>self%mesh
allocate(nr(self%ne+1))
allocate(kee(self%ne+1))
nr=0 ! initialize interaction counter
!$omp parallel private(js,je,jn,j,k,lcx)
allocate(lcx(other%nce,self%necmax+1))
!$omp do
do i=1,self%ne ! loop over dofs
  js=self%kec(i)
  je=self%kec(i+1)-1
  IF(je<js)CYCLE
  do j=js,je  ! loop over neighbor cells
    k=self%lec(j)
    call other%ncdofs(k,lcx(:,j-js+1)) ! get DOFs
  end do
  jn=(je-js+1)
  call sort_array(lcx(:,1:jn),int(other%nce,4),jn)
  !
  nr(i)=1
  jn=lcx(1,1)
  do j=js,je ! loop over neighbor cells
    do k=1,other%nce ! loop over cell dofs
      if(lcx(k,j-js+1)/=jn)then
        nr(i)=nr(i)+1
        jn=lcx(k,j-js+1)
      end if
    end do
  end do
end do
!$omp master
nee=sum(nr)
kee(self%ne+1)=nee+1
do i=self%ne,1,-1 ! cumulative point to point count
  kee(i)=kee(i+1)-nr(i)
end do
if(kee(1)/=1)call oft_abort('Bad element to element count','fem_common_linkage',__FILE__)
allocate(lee(nee))
nr=0 ! reset interaction counter
!$omp end master
!$omp barrier
!$omp do
do i=1,self%ne ! loop over points
  js=self%kec(i)
  je=self%kec(i+1)-1
  IF(je<js)CYCLE
  do j=js,je  ! loop over neighbor cells
    k=self%lec(j)
    call other%ncdofs(k,lcx(:,j-js+1)) ! get DOFs
  end do
  jn=(je-js+1)
  call sort_array(lcx(:,1:jn),int(other%nce,4),jn)
  !
  nr(i)=1
  jn=lcx(1,1)
  lee(kee(i))=jn
  do j=js,je ! loop over neighbor cells
    do k=1,other%nce ! loop over cell dofs
      if(lcx(k,j-js+1)/=jn)then
        lee(kee(i)+nr(i))=lcx(k,j-js+1)
        nr(i)=nr(i)+1
        jn=lcx(k,j-js+1)
      end if
    end do
  end do
end do
deallocate(lcx)
!$omp end parallel
deallocate(nr)
CALL afem_fill_lgraph(self,other,nee,kee,lee)
DEBUG_STACK_POP
end subroutine fem_common_linkage
!---------------------------------------------------------------------------------
!> Supplement local graph with interactions from other processors
!!
!! Due to domain decomposition not all matrix elements for the local block
!! may be present in the graph constructed from the local mesh. This subroutine
!! expands the local matrix graph to include all entries in the local block
!! (interaction between elements owned by the local processor)
!---------------------------------------------------------------------------------
subroutine afem_fill_lgraph(row,col,nee,kee,lee)
class(oft_afem_type), intent(in) :: row !< FE representation for row space
class(oft_afem_type), intent(in) :: col !< FE representation for column space
integer(i4), intent(inout) :: nee !< Number of entries in graph
integer(i4), pointer, intent(inout) :: kee(:) !< Row pointer into column list [self%ne+1]
integer(i4), pointer, intent(inout) :: lee(:) !< Column list [nee]
type :: eout
  integer(i8), pointer, dimension(:,:) :: le => NULL()
end type eout
type(eout), allocatable, dimension(:) :: lerecv,lesend
TYPE(oft_graph) :: new_graph
integer(i4) :: neel,nnz_new,nproc_con,nbe,nbetmp,nbemax
integer(i4) :: i,j,k,m,mm,ierr,jp,jn,jp_new,jn_new,jk,neout,netmp
integer(i4), allocatable, dimension(:) :: isort,lsort,ncon,nrecv
integer(i4), allocatable, dimension(:) :: row_lg_ind,col_lg_ind,lle_tmp,lrtmp
integer(i8), allocatable, dimension(:) :: row_lg_sorted,col_lg_sorted
integer(i8), pointer, dimension(:,:) :: letmp,leout
logical :: periodic
logical, allocatable, dimension(:) :: glob_irow,glob_icol
logical, allocatable, dimension(:,:) :: row_skips,col_skips
DEBUG_STACK_PUSH
IF(oft_debug_print(2))WRITE(*,'(2A)')oft_indent,'Filling local graph'
CALL oft_increase_indent
periodic=(row%map%per.OR.col%map%per)
!---Create sorted global list
ALLOCATE(row_lg_ind(row%ne),row_lg_sorted(row%ne))
row_lg_ind=(/(i,i=1,row%ne)/)
row_lg_sorted=row%global%le
!---Flag redundant rows
IF(periodic)THEN
  DO i=row%linkage%kle(0),row%linkage%kle(1)-1
    j=row%linkage%lbe(row%linkage%lle(1,i))
    k=row%linkage%lbe(row%linkage%lle(2,i))
    IF(j>k)row_lg_sorted(ABS(j))=-ABS(row_lg_sorted(ABS(j)))
  END DO
  DO i=1,row%nbe
    IF(row%linkage%leo(i).AND.(row_lg_sorted(row%lbe(i))<0))THEN
      CALL oft_abort('Bad redundant row ownership','fem_fill_lgraph',__FILE__)
    END IF
  END DO
END IF
CALL sort_array(row_lg_sorted,row_lg_ind,row%ne)
!---
ALLOCATE(col_lg_ind(col%ne),col_lg_sorted(col%ne))
col_lg_ind=(/(i,i=1,col%ne)/)
col_lg_sorted=col%global%le
!---Flag redundant columns
IF(periodic)THEN
  DO i=col%linkage%kle(0),col%linkage%kle(1)-1
    j=col%linkage%lbe(col%linkage%lle(1,i))
    k=col%linkage%lbe(col%linkage%lle(2,i))
    IF(j>k)col_lg_sorted(ABS(j))=-ABS(col_lg_sorted(ABS(j)))
  END DO
  DO i=1,col%nbe
    IF(col%linkage%leo(i).AND.(col_lg_sorted(col%lbe(i))<0))THEN
      CALL oft_abort('Bad redundant column ownership','fem_fill_lgraph',__FILE__)
    END IF
  END DO
END IF
CALL sort_array(col_lg_sorted,col_lg_ind,col%ne)
!---
nbetmp=0
ALLOCATE(glob_irow(row%ne),glob_icol(col%ne))
glob_irow=.FALSE.
glob_icol=.FALSE.
IF(periodic)THEN
  DO i=row%linkage%kle(0),row%linkage%kle(1)-1
    k = row%linkage%lbe(row%linkage%lle(2,i))
    glob_irow(k)=.TRUE.
  END DO
  DO i=col%linkage%kle(0),col%linkage%kle(1)-1
    k = col%linkage%lbe(col%linkage%lle(2,i))
    glob_icol(k)=.TRUE.
  END DO
END IF
!$omp parallel do private(j,k) reduction(+:nbetmp)
DO i=1,row%nbe
  k=row%lbe(i)
  IF(glob_irow(k))THEN
    nbetmp=nbetmp+(kee(k+1)-kee(k))
  ELSE
    DO j=kee(k),kee(k+1)-1
      IF(col%be(lee(j)))nbetmp=nbetmp+1
    END DO
  END IF
END DO
IF(periodic)THEN
  !$omp parallel do private(j) reduction(+:nbetmp)
  DO i=1,row%ne
    IF(row%be(i))CYCLE
    DO j=kee(i),kee(i+1)-1
      IF(glob_icol(lee(j)))nbetmp=nbetmp+1
    END DO
  END DO
END IF
nbe=nbetmp
nbemax=nbetmp
IF(.NOT.row%linkage%full)nbemax=oft_mpi_max(nbetmp)
IF(oft_debug_print(3))WRITE(*,'(2A,I8)')oft_indent,'Max # of boundary elements',nbemax
!---Allocate temporary Send/Recv arrays
ALLOCATE(lerecv(0:row%linkage%nproc_con),lesend(0:row%linkage%nproc_con))
ALLOCATE(lesend(0)%le(nbemax,2))
!---Setup sorted linkage lists
IF(.NOT.row%linkage%full)THEN
  !---Row vector space
  ALLOCATE(lle_tmp(row%linkage%nle))
  lle_tmp=row%linkage%lle(1,:)
  ALLOCATE(row_skips(row%linkage%nproc_con,row%ne))
  row_skips=.TRUE.
  !$omp parallel do private(jp,jn,i,jk) schedule(dynamic,1)
  DO j=1,row%linkage%nproc_con
    jp=row%linkage%kle(j)
    jn=row%linkage%kle(j+1)-jp
    IF(jn>0)THEN
      CALL sort_array(lle_tmp(jp:jp+jn-1),jn)
      DO i=1,row%nbe
        jk=search_array(i,lle_tmp(jp:jp+jn-1),jn)
        IF(jk/=0)row_skips(j,row%lbe(i))=.FALSE.
      END DO
    END IF
  END DO
  DEALLOCATE(lle_tmp)
  !---Column vector space
  ALLOCATE(lle_tmp(col%linkage%nle))
  lle_tmp=col%linkage%lle(1,:)
  ALLOCATE(col_skips(row%linkage%nproc_con,col%ne))
  col_skips=.TRUE.
  !$omp parallel do private(jp,jn,i,jk) schedule(dynamic,1)
  DO j=1,row%linkage%nproc_con
    jp=col%linkage%kle(j)
    jn=col%linkage%kle(j+1)-jp
    IF(jn>0)THEN
      CALL sort_array(lle_tmp(jp:jp+jn-1),jn)
      DO i=1,col%nbe
        jk=search_array(i,lle_tmp(jp:jp+jn-1),jn)
        IF(jk/=0)col_skips(j,col%lbe(i))=.FALSE.
      END DO
    END IF
  END DO
  DEALLOCATE(lle_tmp)
END IF
!---
ALLOCATE(ncon(0:row%linkage%nproc_con),nrecv(0:row%linkage%nproc_con))
ncon=0
DO i=1,row%nbe
  k=row%lbe(i)
  IF(glob_irow(k))THEN
    DO j=kee(k),kee(k+1)-1
      ncon(0)=ncon(0)+1
      lesend(0)%le(ncon(0),:) = (/row%global%le(k),col%global%le(lee(j))/)
      IF(.NOT.row%linkage%full)THEN
        DO m=1,row%linkage%nproc_con
          IF(row_skips(m,k).AND.col_skips(m,lee(j)))CYCLE
          ncon(m)=ncon(m)+1
        END DO
      END IF
    END DO
  ELSE
    DO j=kee(k),kee(k+1)-1
      IF(col%be(lee(j)))THEN
        ncon(0)=ncon(0)+1
        lesend(0)%le(ncon(0),:) = (/row%global%le(k),col%global%le(lee(j))/)
        IF(.NOT.row%linkage%full)THEN
          DO m=1,row%linkage%nproc_con
            IF(row_skips(m,k).OR.col_skips(m,lee(j)))CYCLE
            ncon(m)=ncon(m)+1
          END DO
        END IF
      END IF
    END DO
  END IF
END DO
!---
IF(periodic)THEN
  DO i=1,row%ne
    IF(row%be(i))CYCLE
    DO j=kee(i),kee(i+1)-1
      IF(glob_icol(lee(j)))THEN
        ncon(0)=ncon(0)+1
        lesend(0)%le(ncon(0),:) = (/row%global%le(i),col%global%le(lee(j))/)
        IF(.NOT.row%linkage%full)THEN
          DO m=1,row%linkage%nproc_con
            IF(row_skips(m,i).AND.col_skips(m,lee(j)))CYCLE
            ncon(m)=ncon(m)+1
          END DO
        END IF
      END IF
    END DO
  END DO
END IF
nrecv(0)=nbemax
IF(.NOT.row%linkage%full)THEN
#ifdef HAVE_MPI
  DO j=1,row%linkage%nproc_con
    CALL MPI_ISEND(ncon(j),1,OFT_MPI_I4,row%linkage%proc_con(j),1,oft_env%COMM,row%linkage%send_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','fem_fill_lgraph',__FILE__)
    CALL MPI_IRECV(nrecv(j),1,OFT_MPI_I4,row%linkage%proc_con(j),1,oft_env%COMM,row%linkage%recv_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','fem_fill_lgraph',__FILE__)
  END DO
  CALL oft_mpi_waitall(row%linkage%nproc_con,row%linkage%recv_reqs,ierr)
  DO i=1,row%linkage%nproc_con
    ALLOCATE(lesend(i)%le(ncon(i),2))
    lesend(i)%le=0
    ALLOCATE(lerecv(i)%le(nrecv(i),2))
  END DO
  CALL oft_mpi_waitall(row%linkage%nproc_con,row%linkage%send_reqs,ierr)
  ncon=0
  DO i=1,row%nbe
    k=row%lbe(i)
    IF(glob_irow(k))THEN
      DO j=kee(k),kee(k+1)-1
        DO m=1,row%linkage%nproc_con
          IF(row_skips(m,k).AND.col_skips(m,lee(j)))CYCLE
          ncon(m)=ncon(m)+1
          lesend(m)%le(ncon(m),:) = (/row%global%le(k),col%global%le(lee(j))/)
        END DO
      END DO
    ELSE
      DO j=kee(k),kee(k+1)-1
        IF(col%be(lee(j)))THEN
          DO m=1,row%linkage%nproc_con
            IF(row_skips(m,k).OR.col_skips(m,lee(j)))CYCLE
            ncon(m)=ncon(m)+1
            lesend(m)%le(ncon(m),:) = (/row%global%le(k),col%global%le(lee(j))/)
          END DO
        END IF
      END DO
    END IF
  END DO
  !---
  IF(periodic)THEN
    DO i=1,row%ne
      IF(row%be(i))CYCLE
      DO j=kee(i),kee(i+1)-1
        IF(glob_icol(lee(j)))THEN
          DO m=1,row%linkage%nproc_con
            IF(row_skips(m,i).AND.col_skips(m,lee(j)))CYCLE
            ncon(m)=ncon(m)+1
            lesend(m)%le(ncon(m),:) = (/row%global%le(i),col%global%le(lee(j))/)
          END DO
        END IF
      END DO
    END DO
  END IF
#else
CALL oft_abort("Distributed linkage requires MPI","fem_fill_lgraph",__FILE__)
#endif
END IF
DEALLOCATE(glob_irow,glob_icol)
IF(.NOT.row%linkage%full)DEALLOCATE(row_skips,col_skips)
nnz_new=nee
!---
IF(.NOT.row%linkage%full)THEN
#ifdef HAVE_MPI
  CALL oft_mpi_barrier(ierr) ! Wait for all processes
  !---Create Send and Recv calls
  do j=1,row%linkage%nproc_con
    call MPI_ISEND(lesend(j)%le,2*ncon(j),OFT_MPI_I8,row%linkage%proc_con(j),1,oft_env%COMM,row%linkage%send_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','fem_fill_lgraph',__FILE__)
    call MPI_IRECV(lerecv(j)%le,2*nrecv(j),OFT_MPI_I8,row%linkage%proc_con(j),1,oft_env%COMM,row%linkage%recv_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','fem_fill_lgraph',__FILE__)
  end do
  !---Loop over each connected processor
  do while(.TRUE.)
    IF(oft_mpi_check_reqs(row%linkage%nproc_con,row%linkage%recv_reqs))EXIT ! All recieves have been processed
    CALL oft_mpi_waitany(row%linkage%nproc_con,row%linkage%recv_reqs,j,ierr) ! Wait for completed recieve
    IF(ierr/=0)CALL oft_abort('Error in MPI_WAITANY','fem_fill_lgraph',__FILE__)
    letmp=>lerecv(j)%le ! Point dummy input array to current Recv array
    netmp=nrecv(j)
    ncon(j)=0
    !---Convert to local indexing
    DO m=1,netmp
      jp=INT(search_array(letmp(m,1),row_lg_sorted,INT(row%ne,8)),4)
      IF(jp>0)THEN
        letmp(m,1)=row_lg_ind(jp)
        jn=INT(search_array(letmp(m,2),col_lg_sorted,INT(col%ne,8)),4)
        IF(jn>0)THEN
          letmp(m,2)=col_lg_ind(jn)
          nnz_new=nnz_new+1
        ELSE
          letmp(m,1)=HUGE(INT(1,8))
        END IF
      ELSE
        letmp(m,1)=HUGE(INT(1,8))
      END IF
    END DO
    !---
    CALL sort_array(letmp(:,1),letmp(:,2),INT(netmp,8))
  end do
#else
CALL oft_abort("Distributed linkage requires MPI","fem_fill_lgraph",__FILE__)
#endif
END IF
!---Local connectivity
ncon(0)=0
leout=>lesend(0)%le
lerecv(0)%le=>lesend(0)%le
jp=row%linkage%kle(0)
jn=row%linkage%kle(0+1)-jp
IF(jn>0)THEN
  DO m=1,nbe
    jp=INT(search_array(leout(m,1),row_lg_sorted,INT(row%ne,8)),4)
    IF(jp>0)THEN
      leout(m,1)=row_lg_ind(jp)
      jn=INT(search_array(leout(m,2),col_lg_sorted,INT(col%ne,8)),4)
      IF(jn>0)THEN
        leout(m,2)=col_lg_ind(jn)
        nnz_new=nnz_new+1
      ELSE
        leout(m,1)=HUGE(INT(1,8))
      END IF
    ELSE
      leout(m,1)=HUGE(INT(1,8))
    END IF
  END DO
  CALL sort_array(leout(:,1),leout(:,2),INT(nbe,8))
  DO i=nbe+1,nbemax
    leout(i,1)=HUGE(INT(1,8))
  END DO
ELSE
  DO i=1,nbemax
    leout(i,1)=HUGE(INT(1,8))
  END DO
END IF
DEALLOCATE(row_lg_sorted,col_lg_sorted)
!---Create updated graph
ALLOCATE(lrtmp(nnz_new),new_graph%kr(row%ne+1))
nnz_new=0
new_graph%kr=0
ncon=1
nproc_con=row%linkage%nproc_con
IF(row%linkage%full)nproc_con=0
DO i=1,row%ne
  DO j=kee(i),kee(i+1)-1
    new_graph%kr(i)=new_graph%kr(i)+1
    lrtmp(nnz_new+new_graph%kr(i))=lee(j)
  END DO
  !---
  DO m=0,nproc_con
    letmp=>lerecv(m)%le
    IF(ncon(m)>nrecv(m))CYCLE
    DO WHILE(letmp(ncon(m),1)<=i)
      IF(letmp(ncon(m),1)==i)THEN
        new_graph%kr(i)=new_graph%kr(i)+1
        lrtmp(nnz_new+new_graph%kr(i))=INT(letmp(ncon(m),2),4)
      END IF
      ncon(m)=ncon(m)+1
      IF(ncon(m)>nrecv(m))EXIT
    END DO
  END DO
  nnz_new=nnz_new+new_graph%kr(i)
END DO
DEALLOCATE(ncon,nrecv)
new_graph%kr(row%ne+1)=nnz_new+1
do i=row%ne,1,-1 ! cumulative unique point linkage count
  new_graph%kr(i)=new_graph%kr(i+1)-new_graph%kr(i)
end do
if(new_graph%kr(1)/=1)call oft_abort('Bad new graph setup','fem_fill_lgraph', &
__FILE__)
!---
DO i=1,row%ne
  jp=new_graph%kr(i)
  jn=new_graph%kr(i+1)-1
  CALL sort_array(lrtmp(jp:jn),jn-jp+1)
  m=lrtmp(jp)
  new_graph%kr(i)=1
  DO j=jp+1,jn
    IF(m==lrtmp(j))THEN
      lrtmp(j)=-1
    ELSE
      new_graph%kr(i)=new_graph%kr(i)+1
      m=lrtmp(j)
    END IF
  END DO
END DO
!---
new_graph%kr(row%ne+1)=0
new_graph%nnz=SUM(new_graph%kr)
new_graph%kr(row%ne+1)=new_graph%nnz+1
do i=row%ne,1,-1 ! cumulative unique point linkage count
  new_graph%kr(i)=new_graph%kr(i+1)-new_graph%kr(i)
end do
if(new_graph%kr(1)/=1)call oft_abort('Bad new graph setup','fem_fill_lgraph', &
__FILE__)
ALLOCATE(new_graph%lc(new_graph%nnz))
j=0
DO i=1,nnz_new
  IF(lrtmp(i)<0)CYCLE
  j=j+1
  new_graph%lc(j)=lrtmp(i)
END DO
if(new_graph%nnz/=j)call oft_abort('Bad new graph setup','fem_fill_lgraph', &
__FILE__)
DEALLOCATE(lesend(0)%le)
IF(.NOT.row%linkage%full)THEN
#ifdef HAVE_MPI
  CALL oft_mpi_waitall(row%linkage%nproc_con,row%linkage%send_reqs,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','fem_fill_lgraph',__FILE__)
  CALL oft_mpi_barrier(ierr)
  DO i=1,row%linkage%nproc_con
    DEALLOCATE(lesend(i)%le)
    DEALLOCATE(lerecv(i)%le)
  END DO
  DEALLOCATE(lesend,lerecv)
#else
CALL oft_abort("Distributed linkage requires MPI","fem_fill_lgraph",__FILE__)
#endif
END IF
!---Replace graph
nee=new_graph%nnz
kee=new_graph%kr
DEALLOCATE(new_graph%kr)
DEALLOCATE(lee)
lee=>new_graph%lc
NULLIFY(new_graph%lc)
! IF(oft_debug_print(3))WRITE(*,'(2A)')oft_indent,'Done'
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine afem_fill_lgraph
!------------------------------------------------------------------------------
!> Constructs a finite element representation on the associated volume mesh
!!
!! Constructs a finite element representation from a specified geometric mapping,
!! provided by @ref fem_base::oft_fem_type::gstruct "gstruct", and the assigned
!! @ref fem_base::oft_fem_type::mesh "mesh". The result is a fully defined finite
!! element type which my be used to define weight vectors, matrix graphs, etc.
!! required to employ the finite element method on the chosen tetrahedral grid.
!!
!! @param[in,out] self FE representation to construct
!! @param[in] quad_order Desired quadrature order
!------------------------------------------------------------------------------
subroutine fem_setup(self,quad_order)
class(oft_fem_type), intent(inout) :: self
integer(i4), intent(in) :: quad_order
integer(i4) :: i,j,k,offset,stack
integer(i4) :: jeg,jsg,jee,jse,necmax,js,je,jn
integer(i4) :: mycounts(4)
integer(i4), allocatable :: lcx(:,:),nr(:)
class(oft_mesh), pointer :: mesh
DEBUG_STACK_PUSH
if(oft_debug_print(2))write(*,'(2A)')oft_indent,'Starting FE setup'
CALL oft_increase_indent
!---
mesh=>self%mesh
!---Set quadrature order
CALL mesh%quad_rule(quad_order, self%quad)
!---Set number of node pts
self%nnodes=self%quad%np*mesh%nc
!---Compute DOFs from geometric linkage
mycounts=(/mesh%np,mesh%ne,mesh%nf,mesh%nc/)
self%ne=sum(self%gstruct*mycounts)
!---Compute boundary count from geometric linkage
mycounts=(/mesh%nbp,mesh%nbe,mesh%nbf,0/)
self%nbe=sum(self%gstruct*mycounts)
allocate(self%be(self%ne),self%lbe(self%nbe))
self%be=.FALSE.
!---
allocate(self%ce(self%ne))
self%ce=.FALSE.
!---Compute cell DOFs from geometric linkage
mycounts(1:3)=[mesh%cell_np,mesh%cell_ne,mesh%cell_nf]
mycounts(4)=1
self%nce=sum(self%gstruct*mycounts)
!---Count element to cell interactions
allocate(self%kec(self%ne+1))
stack=1
self%nbe=0
!---Point stacking
do i=1,mesh%np
  offset=(i-1)*self%gstruct(1)
  do j=1,self%gstruct(1)
    self%kec(j+offset)=stack
    stack=stack+mesh%kpc(i+1)-mesh%kpc(i)
    if(mesh%bp(i))then
      self%be(j+offset)=.TRUE.
      self%nbe=self%nbe+1
      self%lbe(self%nbe)=j+offset
    end if
    IF(mesh%cp(i))then
      self%ce(j+offset)=.TRUE.
    END IF
  end do
end do
!---Edge stacking
do i=1,mesh%ne
  offset=(i-1)*self%gstruct(2)+self%gstruct(1)*mesh%np
  do j=1,self%gstruct(2)
    self%kec(j+offset)=stack
    stack=stack+mesh%kec(i+1)-mesh%kec(i)
    if(mesh%be(i))then
      self%be(j+offset)=.TRUE.
      self%nbe=self%nbe+1
      self%lbe(self%nbe)=j+offset
    end if
    IF(mesh%ce(i))then
      self%ce(j+offset)=.TRUE.
    END IF
  end do
end do
!---Face stacking
do i=1,mesh%nf
  offset=(i-1)*self%gstruct(3)+self%gstruct(2)*mesh%ne+self%gstruct(1)*mesh%np
  do j=1,self%gstruct(3)
    self%kec(j+offset)=stack
    if(mesh%lfc(2,i)/=0)then
      stack=stack+2
    else
      stack=stack+1
    end if
    if(mesh%bf(i))then
      self%be(j+offset)=.TRUE.
      self%nbe=self%nbe+1
      self%lbe(self%nbe)=j+offset
    end if
  end do
end do
!---Cell stacking
do i=1,mesh%nc
  offset=(i-1)*self%gstruct(4)+self%gstruct(3)*mesh%nf+self%gstruct(2)*mesh%ne+self%gstruct(1)*mesh%np
  do j=1,self%gstruct(4)
    self%kec(j+offset)=stack
    stack=stack+1
  end do
end do
self%kec(self%ne+1)=stack
self%nec=stack-1
!---Assemble element to cell interaction list
allocate(self%lec(self%nec))
necmax=0
!$omp parallel private(offset,j,jee,jse,jeg,jsg) reduction(max:necmax)
!---Point stacking
!$omp do
do i=1,mesh%np
  offset=(i-1)*self%gstruct(1)
  do j=1,self%gstruct(1)
    jee=self%kec(j+offset+1)-1; jse=self%kec(j+offset)
    jeg=mesh%kpc(i+1)-1; jsg=mesh%kpc(i)
    necmax=max(necmax,jee-jse)
    self%lec(jse:jee)=mesh%lpc(jsg:jeg)
  end do
end do
!---Edge stacking
!$omp do
do i=1,mesh%ne
  offset=(i-1)*self%gstruct(2)+self%gstruct(1)*mesh%np
  do j=1,self%gstruct(2)
    jee=self%kec(j+offset+1)-1; jse=self%kec(j+offset)
    jeg=mesh%kec(i+1)-1; jsg=mesh%kec(i)
    necmax=max(necmax,jee-jse)
    self%lec(jse:jee)=mesh%lec(jsg:jeg)
  end do
end do
!---Face stacking
!$omp do
do i=1,mesh%nf
  offset=(i-1)*self%gstruct(3)+self%gstruct(2)*mesh%ne+self%gstruct(1)*mesh%np
  do j=1,self%gstruct(3)
    jee=self%kec(j+offset+1)-1; jse=self%kec(j+offset)
    if(self%mesh%lfc(2,i)/=0)then
      self%lec(jse:jee)=mesh%lfc(:,i)
    else
      self%lec(jse:jee)=mesh%lfc(1,i)
    end if
  end do
end do
!---Cell stacking
!$omp do
do i=1,mesh%nc
  offset=(i-1)*self%gstruct(4)+self%gstruct(3)*mesh%nf+self%gstruct(2)*mesh%ne &
  + self%gstruct(1)*mesh%np
  do j=1,self%gstruct(4)
    jse=self%kec(j+offset)
    self%lec(jse)=i
  end do
end do
!$omp end parallel
self%necmax=necmax
if(oft_debug_print(2))write(*,'(2A)')oft_indent,'Creating interaction graph'
! call fem_self_linkage(self)
call afem_self_linkage(self)
if(oft_debug_print(2))write(*,'(2A)')oft_indent,'Creating cell mapping'
call fem_ncdofs_map(self)
if(oft_debug_print(2))write(*,'(2A)')oft_indent,'Constructing MPI linkage'
call fem_global_linkage(self)
CALL afem_fill_lgraph(self,self,self%nee,self%kee,self%lee)
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine fem_setup
!------------------------------------------------------------------------------
!> Destroy FE object
!------------------------------------------------------------------------------
SUBROUTINE fem_delete(self)
CLASS(oft_fem_type), INTENT(inout) :: self
DEBUG_STACK_PUSH
CALL afem_delete(self)
NULLIFY(self%mesh)
DEBUG_STACK_POP
END SUBROUTINE fem_delete
!------------------------------------------------------------------------------
!> Compute FE global context and stitching information
!!
!! Sets up structures and information for distributed meshes. Primarily this
!! supports the creation of vector and matrix objects for the full mesh during
!! linear solves.
!!
!! @param[in,out] self Finite element representation
!------------------------------------------------------------------------------
subroutine fem_global_linkage(self)
class(oft_fem_type), intent(inout) :: self
type :: eout
  integer(i8), pointer, dimension(:) :: le => NULL()
end type eout
type(eout), allocatable, dimension(:) :: lerecv,lesend
integer(i8), pointer, dimension(:) :: letmp,leout
integer(i8), allocatable, dimension(:) :: lsort
integer(i4), allocatable :: linktmp(:,:,:),isort(:),ncon(:)
integer(i4) :: neel
integer(i4) :: kk,m,mm,l,request
integer(i4) :: nbetmp,nbemax
integer(i8) :: netmp,nemax
integer(i4) :: i,j,k,offset,stack1,stack2
integer(i4) :: je,js,jee,jse,ierr
integer(i8) :: offsetg
integer(i8) :: mycounts(4)
integer(i4), allocatable, dimension(:) :: bpi,bei,bfi
class(oft_mesh), pointer :: mesh
#ifdef HAVE_MPI
integer(i4) :: stat(MPI_STATUS_SIZE)
#endif
DEBUG_STACK_PUSH
!---
mesh=>self%mesh
IF(.NOT.mesh%fullmesh)CALL oft_mpi_barrier(ierr) ! Wait for all processes
CALL oft_increase_indent
!---Copy seam information from mesh
allocate(self%linkage,self%global)
CALL oft_init_seam(mesh,self%linkage)
!---Determine global element count
mycounts=(/mesh%global%np,mesh%global%ne,mesh%global%nf,mesh%global%nc/)
self%global%ne=sum(self%gstruct*mycounts)
allocate(self%global%le(self%ne))
!---
IF(mesh%fullmesh)THEN
  self%linkage%full=.TRUE.
ELSE
  self%linkage%full=.FALSE.
END IF
allocate(self%global%gbe(self%ne))
self%global%gbe=.FALSE.
!$omp parallel private(offset,offsetg,j)
!---Point stacking
!$omp do
do i=1,mesh%np
  offset=(i-1)*self%gstruct(1)
  offsetg=(mesh%global%lp(i)-1)*self%gstruct(1)
  do j=1,self%gstruct(1)
    self%global%le(j+offset)=j+offsetg
    self%global%gbe(j+offset)=mesh%global%gbp(i)
  end do
end do
!---Edge stacking
!$omp do
do i=1,mesh%ne
  offset = (i-1)*self%gstruct(2)+self%gstruct(1)*mesh%np
  offsetg = (ABS(mesh%global%le(i))-1)*self%gstruct(2)+self%gstruct(1)*mesh%global%np
  do j=1,self%gstruct(2)
    self%global%le(j+offset)=j+offsetg
    self%global%gbe(j+offset)=mesh%global%gbe(i)
  end do
end do
!---Face stacking
!$omp do
do i=1,mesh%nf
  offset = (i-1)*self%gstruct(3)+self%gstruct(2)*mesh%ne+self%gstruct(1)*mesh%np
  offsetg = (mesh%global%lf(i)-1)*self%gstruct(3)+self%gstruct(2)*mesh%global%ne+self%gstruct(1)*mesh%global%np
  do j=1,self%gstruct(3)
    self%global%le(j+offset)=j+offsetg
    self%global%gbe(j+offset)=mesh%global%gbf(i)
  end do
end do
!---Cell stacking
!$omp do
do i=1,mesh%nc
  offset = (i-1)*self%gstruct(4)+self%gstruct(3)*mesh%nf+self%gstruct(2)*mesh%ne+self%gstruct(1)*mesh%np
  offsetg = (mesh%global%lc(i)-1)*self%gstruct(4)+self%gstruct(3)*mesh%global%nf &
  + self%gstruct(2)*mesh%global%ne+self%gstruct(1)*mesh%global%np
  do j=1,self%gstruct(4)
    self%global%le(j+offset)=j+offsetg
  end do
end do
!$omp end parallel
!---
ALLOCATE(linktmp(2,5*self%nbe,0:self%linkage%nproc_con))
linktmp = 0
ALLOCATE(ncon(0:self%linkage%nproc_con))
ncon=0
ALLOCATE(self%linkage%kle(0:self%linkage%nproc_con+1)) ! Allocate point linkage arrays
self%linkage%kle=0
self%linkage%nbe=self%nbe
self%linkage%lbe=>self%lbe
self%linkage%be=>self%be
linktmp=0
!---Determine maximum boundary element count
nbetmp=self%nbe
nbemax=self%nbe
IF(.NOT.mesh%fullmesh)nbemax=oft_mpi_max(nbetmp)
if(oft_debug_print(3))write(*,'(2A,I8)')oft_indent,'Max # of boundary elements',nbemax
self%linkage%nbemax=nbemax
!---Allocate temporary Send/Recv arrays
ALLOCATE(lesend(self%linkage%nproc_con+1),lerecv(self%linkage%nproc_con+1))
allocate(lesend(1)%le(nbemax))
IF(.NOT.mesh%fullmesh)THEN
  do i=1,self%linkage%nproc_con
    allocate(lerecv(i)%le(nbemax))
  end do
END IF
!---Point dummy output array to Send array
leout=>lesend(1)%le
leout=0
!$omp parallel do
do i=1,self%nbe
  leout(i) = self%global%le(self%lbe(i)) ! Populate with global point index
end do
IF(.NOT.mesh%fullmesh)THEN
#ifdef HAVE_MPI
  CALL oft_mpi_barrier(ierr) ! Wait for all processes
  !---Point dummy Send arrays to main Send array
  do i=2,self%linkage%nproc_con
    lesend(i)%le=>lesend(1)%le
  end do
  !---Create Send and Recv calls
  do j=1,self%linkage%nproc_con
    CALL MPI_ISEND(lesend(j)%le,nbemax,OFT_MPI_I8,self%linkage%proc_con(j),1,oft_env%COMM,self%linkage%send_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','fem_global_linkage',__FILE__)
    CALL MPI_IRECV(lerecv(j)%le,nbemax,OFT_MPI_I8,self%linkage%proc_con(j),1,oft_env%COMM,self%linkage%recv_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','fem_global_linkage',__FILE__)
  end do
  !---Loop over each connected processor
  do while(.TRUE.)
    IF(oft_mpi_check_reqs(self%linkage%nproc_con,self%linkage%recv_reqs))EXIT ! All recieves have been processed
    CALL oft_mpi_waitany(self%linkage%nproc_con,self%linkage%recv_reqs,j,ierr) ! Wait for completed recieve
    IF(ierr/=0)CALL oft_abort('Error in MPI_WAITANY','fem_global_linkage',__FILE__)
    letmp=>lerecv(j)%le ! Point dummy input array to current Recv array
    !---Determine location of boundary points on other processors
    neel=0
    !$omp parallel do private(mm) reduction(+:neel)
    do m=1,self%nbe ! Loop over boundary points
      do mm=1,nbemax ! Loop over input points
        if(leout(m)==letmp(mm))then ! Found match
          !$omp critical
          ncon(j)=ncon(j)+1
          linktmp(:,ncon(j),j)=(/mm,m/)
          !$omp end critical
          neel=neel+1
        end if
        if(letmp(mm)==0)exit ! End of input points
      end do
    end do
    self%linkage%kle(j)=neel
  end do
#else
CALL oft_abort("Distributed linkage requires MPI","fem_global_linkage",__FILE__)
#endif
END IF
!---Local connectivity
!---Determine location of boundary points on other processors
neel=0
IF(self%mesh%periodic%nper>0)THEN
  !$omp parallel do private(mm) reduction(+:neel)
  do m=1,self%nbe ! Loop over boundary points
    do mm=m+1,self%nbe ! Loop over input points
      if(leout(m)==leout(mm))then ! Found match
        !$omp critical
        ncon(0)=ncon(0)+1
        linktmp(:,ncon(0),0)=(/mm,m/)
        ncon(0)=ncon(0)+1
        linktmp(:,ncon(0),0)=(/m,mm/)
        !$omp end critical
        neel=neel+2
      end if
      if(leout(mm)==0)exit ! End of input points
    end do
  end do
END IF
self%linkage%kle(0)=neel
!---Condense linkage to sparse rep
self%linkage%nle=sum(self%linkage%kle)
self%linkage%kle(self%linkage%nproc_con+1)=self%linkage%nle+1
do i=self%linkage%nproc_con,0,-1 ! cumulative unique point linkage count
  self%linkage%kle(i)=self%linkage%kle(i+1)-self%linkage%kle(i)
end do
!---
if(self%linkage%kle(0)/=1)call oft_abort('Bad element linkage count','fem_global_linkage',__FILE__)
!---
allocate(self%linkage%lle(2,self%linkage%nle))
allocate(self%linkage%send(0:self%linkage%nproc_con),self%linkage%recv(0:self%linkage%nproc_con))
!---
!$omp parallel private(j,m,lsort,isort)
ALLOCATE(lsort(MAXVAL(ncon)),isort(MAXVAL(ncon)))
!$omp do
do i=0,self%linkage%nproc_con
  !---
  DO j=1,ncon(i)
    m=linktmp(2,j,i)
    lsort(j)=m
    isort(j)=j
    !---
    IF(i>0)THEN
      IF(self%linkage%proc_con(i)<oft_env%rank)lsort(j)=linktmp(1,j,i)
    END IF
  END DO
  !---
  CALL sort_array(lsort,isort,ncon(i))
  DO m=0,ncon(i)-1
    self%linkage%lle(:,m+self%linkage%kle(i)) = &
    (/linktmp(2,isort(m+1),i),linktmp(1,isort(m+1),i)/)
  END DO
  !---Allocate permanent stitching arrays
  self%linkage%send(i)%n=j; self%linkage%recv(i)%n=j
  allocate(self%linkage%send(i)%v(self%linkage%send(i)%n))
  allocate(self%linkage%recv(i)%v(self%linkage%recv(i)%n))
end do
DEALLOCATE(lsort,isort)
!$omp end parallel
!---Propogate ownership
allocate(bpi(mesh%np),bei(mesh%ne),bfi(mesh%nf))
CALL get_inverse_map(mesh%lbp,mesh%nbp,bpi,mesh%np)
CALL get_inverse_map(mesh%lbe,mesh%nbe,bei,mesh%ne)
CALL get_inverse_map(mesh%lbf,mesh%nbf,bfi,mesh%nf)
allocate(self%linkage%leo(self%nbe))
mycounts(1)=self%gstruct(1)*mesh%np
mycounts(2)=mycounts(1)+self%gstruct(2)*mesh%ne
mycounts(3)=mycounts(2)+self%gstruct(3)*mesh%nf
!$omp parallel do private(k,j)
do i=1,self%nbe
  k=self%lbe(i)
  if(k<=mycounts(1))then
    j=FLOOR(REAL(k-1,8)/self%gstruct(1))+1
    self%linkage%leo(i)=mesh%pstitch%leo(bpi(j)) !mesh%linkage%lpo(bpi(j))
  else if(k>mycounts(1).AND.k<=mycounts(2))then
    j=FLOOR(REAL(k-mycounts(1)-1,8)/self%gstruct(2))+1
    self%linkage%leo(i)=mesh%estitch%leo(bei(j)) !mesh%linkage%leo(bei(j))
  else if(k>mycounts(2).AND.k<=mycounts(3))then
    j=FLOOR(REAL(k-mycounts(2)-1,8)/self%gstruct(3))+1
    self%linkage%leo(i)=mesh%fstitch%leo(bfi(j)) !mesh%linkage%lfo(bfi(j))
  else
    call oft_abort('Invalid boundary element','fem_global_linkage',__FILE__)
  end if
end do
IF(.NOT.mesh%fullmesh)THEN
#ifdef HAVE_MPI
  CALL oft_mpi_waitall(self%linkage%nproc_con,self%linkage%send_reqs,ierr) ! Wait for all sends to complete
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','fem_global_linkage',__FILE__)
  !---Deallocate temporary work arrays
  do i=1,self%linkage%nproc_con
    deallocate(lerecv(i)%le)
  end do
#else
CALL oft_abort("Distributed linkage requires MPI","fem_global_linkage",__FILE__)
#endif
END IF
CALL oft_stitch_check(self%linkage)
deallocate(linktmp,ncon)
deallocate(bpi,bei,bfi)
DEALLOCATE(lesend(1)%le)
DEALLOCATE(lesend,lerecv)
!---Construct map
ALLOCATE(self%map)
self%map%per=(mesh%periodic%nper>0)
self%map%offset=0
self%map%n=self%ne
self%map%ng=self%global%ne
self%map%nslice=0
DO i=1,self%ne
  IF(self%linkage%be(i))CYCLE
  self%map%nslice=self%map%nslice+1
END DO
DO i=1,self%linkage%nbe
  IF(.NOT.self%linkage%leo(i))CYCLE
  self%map%nslice=self%map%nslice+1
END DO
ALLOCATE(self%map%slice(self%map%nslice))
self%map%nslice=0
DO i=1,self%ne
  IF(self%linkage%be(i))CYCLE
  self%map%nslice=self%map%nslice+1
  self%map%slice(self%map%nslice)=i
END DO
DO i=1,self%linkage%nbe
  IF(.NOT.self%linkage%leo(i))CYCLE
  self%map%nslice=self%map%nslice+1
  self%map%slice(self%map%nslice)=self%linkage%lbe(i)
END DO
CALL sort_array(self%map%slice,self%map%nslice)
self%map%lge=>self%global%le
self%map%gbe=>self%global%gbe
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine fem_global_linkage
!------------------------------------------------------------------------------
!> Retrieve the indices of elements beloning to a given cell
!------------------------------------------------------------------------------
subroutine fem_ncdofs(self,cell,dofs)
class(oft_fem_type), intent(in) :: self !< Finite element representation
integer(i4), intent(in) :: cell !< Desired cell in mesh
integer(i4), intent(inout) :: dofs(:) !< Indices of cell elements [self%nce]
integer(i4) :: coffset,eoffset,boffset,i,j
DEBUG_STACK_PUSH
!---Point DOFs
do i=1,self%mesh%cell_np
  coffset=(i-1)*self%gstruct(1)
  do j=1,self%gstruct(1)
    dofs(j+coffset)=j+(self%mesh%lc(i,cell)-1)*self%gstruct(1)
  end do
end do
boffset=self%gstruct(1)*self%mesh%cell_np
!---Edge DOFs
do i=1,self%mesh%cell_ne
  coffset=(i-1)*self%gstruct(2)+boffset
  eoffset=self%gstruct(1)*self%mesh%np
  do j=1,self%gstruct(2)
    dofs(j+coffset)=j+(abs(self%mesh%lce(i,cell))-1)*self%gstruct(2)+eoffset
  end do
end do
boffset=boffset+self%gstruct(2)*self%mesh%cell_ne
!---Face DOFs
do i=1,self%mesh%cell_nf
  coffset=(i-1)*self%gstruct(3)+boffset
  eoffset=self%gstruct(2)*self%mesh%ne+self%gstruct(1)*self%mesh%np
  do j=1,self%gstruct(3)
    dofs(j+coffset)=j+(abs(self%mesh%lcf(i,cell))-1)*self%gstruct(3)+eoffset
  end do
end do
boffset=boffset+self%gstruct(3)*self%mesh%cell_nf
!---Cell DOFs
do i=1,self%gstruct(4)
  coffset=boffset
  eoffset=self%gstruct(3)*self%mesh%nf &
  + self%gstruct(2)*self%mesh%ne+self%gstruct(1)*self%mesh%np
  dofs(i+coffset)=i+(cell-1)*self%gstruct(4)+eoffset
end do
DEBUG_STACK_POP
end subroutine fem_ncdofs
!------------------------------------------------------------------------------
!> Construct cell element mapping
!!
!! Sets up the structure @ref fem_base::oft_fem_type::cmap "cmap", which defines
!! the local type, index, and geometric linkage of DOFs in a cell
!------------------------------------------------------------------------------
subroutine fem_ncdofs_map(self)
class(oft_fem_type), intent(inout) :: self !< Finite element representation
integer(i4) :: coffset,boffset,i,j
DEBUG_STACK_PUSH
allocate(self%cmap(self%nce))
!---Point DOFs
do i=1,self%mesh%cell_np
  coffset=(i-1)*self%gstruct(1)
  do j=1,self%gstruct(1)
    self%cmap(j+coffset)%type=1
    self%cmap(j+coffset)%el=i
    self%cmap(j+coffset)%ind=j
  end do
end do
boffset=self%gstruct(1)*self%mesh%cell_np
!---Edge DOFs
do i=1,self%mesh%cell_ne
  coffset=(i-1)*self%gstruct(2)+boffset
  do j=1,self%gstruct(2)
    self%cmap(j+coffset)%type=2
    self%cmap(j+coffset)%el=i
    self%cmap(j+coffset)%ind=j
  end do
end do
boffset=boffset+self%gstruct(2)*self%mesh%cell_ne
!---Face DOFs
do i=1,self%mesh%cell_nf
  coffset=(i-1)*self%gstruct(3)+boffset
  do j=1,self%gstruct(3)
    self%cmap(j+coffset)%type=3
    self%cmap(j+coffset)%el=i
    self%cmap(j+coffset)%ind=j
  end do
end do
boffset=boffset+self%gstruct(3)*self%mesh%cell_nf
!---Cell DOFs
do i=1,self%gstruct(4)
  coffset=i+boffset
  self%cmap(coffset)%type=4
  self%cmap(coffset)%el=1
  self%cmap(coffset)%ind=i
end do
DEBUG_STACK_POP
end subroutine fem_ncdofs_map
!------------------------------------------------------------------------------
!> Create weight vector for FE representation
!------------------------------------------------------------------------------
subroutine afem_vec_create(self,new,cache,native)
class(oft_afem_type), intent(inout) :: self
class(oft_vector), pointer, intent(out) :: new !< Vector to create
logical, optional, intent(in) :: cache !< Allow caching (optional)
logical, optional, intent(in) :: native !< Force native representation (optional)
TYPE(map_list) :: maps(1)
TYPE(seam_list) :: stitches(1)
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
    stitches(1)%s=>self%linkage
    maps(1)%m=>self%map
    CALL create_vector(new,stitches,maps)
    IF(do_cache)CALL new%new(self%cache_PETSc)
  END IF
ELSE
  IF(ASSOCIATED(self%cache_native))THEN
    CALL self%cache_native%new(new)
  ELSE
    !---Create vector
    stitches(1)%s=>self%linkage
    maps(1)%m=>self%map
    CALL create_vector(new,stitches,maps,native=native)
    IF(do_cache)CALL new%new(self%cache_native)
  END IF
END IF
DEBUG_STACK_POP
end subroutine afem_vec_create
!------------------------------------------------------------------------------
!> Save a Lagrange scalar field to a HDF5 restart file
!------------------------------------------------------------------------------
subroutine afem_vec_save(self,source,filename,path,append)
class(oft_afem_type), intent(inout) :: self
class(oft_vector), target, intent(inout) :: source !< Source field
character(*), intent(in) :: filename !< Name of destination file
character(*), intent(in) :: path !< Field label in file
logical, optional, intent(in) :: append !< Append to file instead of creating?
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
!---
NULLIFY(valtmp)
IF(oft_debug_print(1))WRITE(*,'(6A)')oft_indent,'Writing "',TRIM(path), &
  '" to file "',TRIM(filename),'"'
CALL self%vec_create(outfield,native=.TRUE.)
IF(.NOT.native_vector_cast(outvec,outfield))CALL oft_abort('Failed to create "outfield".', &
  'fem_vec_save',__FILE__)
!---
CALL source%get_local(valtmp,1)
CALL outfield%restore_local(valtmp)
CALL native_vector_slice_push(outvec,self%global%le,rst_info)
CALL hdf5_write(rst_info,filename,path)
!---
CALL outfield%delete
DEALLOCATE(outfield,valtmp)
CALL hdf5_rst_destroy(rst_info)
DEBUG_STACK_POP
end subroutine afem_vec_save
!------------------------------------------------------------------------------
!> Load a Lagrange scalar field from a HDF5 restart file
!------------------------------------------------------------------------------
subroutine afem_vec_load(self,source,filename,path,err_flag)
class(oft_afem_type), intent(inout) :: self
class(oft_vector), target, intent(inout) :: source !< Destination vector
character(*), intent(in) :: filename !< Name of source file
character(*), intent(in) :: path !< Field path in file
integer(i4), optional, intent(out) :: err_flag !< Error flag
integer(i4) :: version
integer(i8), pointer, dimension(:) :: lge
real(r8), pointer, dimension(:) :: valtmp
class(oft_vector), pointer :: infield
class(oft_native_vector), pointer :: invec
type(hdf5_rst) :: rst_info
logical :: success
DEBUG_STACK_PUSH
CALL hdf5_read(version,filename,fem_idx_path,success=success)
IF(.NOT.success)THEN
  IF(PRESENT(err_flag))THEN
    err_flag=2
    DEBUG_STACK_POP
    RETURN
  ELSE
    CALL oft_abort("OFT_idx_Version not found, legacy files not supported","afem_vec_load",__FILE__)
  END IF
END IF
!---
NULLIFY(valtmp)
IF(oft_env%head_proc.AND.oft_env%pm)WRITE(*,*)'Reading "',TRIM(path), &
  '" from file "',TRIM(filename),'"'
CALL self%vec_create(infield,native=.TRUE.)
IF(.NOT.native_vector_cast(invec,infield))CALL oft_abort('Failed to create "infield".', &
  'fem_vec_load',__FILE__)
!---
lge=>self%global%le
CALL native_vector_slice_push(invec,lge,rst_info,alloc_only=.TRUE.)
IF(PRESENT(err_flag))THEN
  CALL hdf5_read(rst_info,filename,path,success=success)
  IF(success)THEN
    err_flag=0
  ELSE
    err_flag=-1
    CALL infield%delete
    NULLIFY(lge)
    DEALLOCATE(infield)
    CALL hdf5_rst_destroy(rst_info)
    RETURN
  END IF
ELSE
  CALL hdf5_read(rst_info,filename,path)
END IF
CALL native_vector_slice_pop(invec,lge,rst_info)
CALL infield%get_local(valtmp)
CALL source%restore_local(valtmp)
!---
CALL infield%delete
NULLIFY(lge)
DEALLOCATE(infield,valtmp)
CALL hdf5_rst_destroy(rst_info)
DEBUG_STACK_POP
end subroutine afem_vec_load
!------------------------------------------------------------------------------
!> Create self-matrix for FE representation
!------------------------------------------------------------------------------
subroutine afem_mat_create(self,new)
CLASS(oft_afem_type), INTENT(inout) :: self
CLASS(oft_matrix), POINTER, INTENT(out) :: new !< Matrix to create
CLASS(oft_vector), POINTER :: tmp_vec
TYPE(oft_graph_ptr) :: graphs(1,1)
DEBUG_STACK_PUSH
CALL self%vec_create(tmp_vec)
ALLOCATE(graphs(1,1)%g)
graphs(1,1)%g%nr=self%ne
graphs(1,1)%g%nrg=self%global%ne
graphs(1,1)%g%nc=self%ne
graphs(1,1)%g%ncg=self%global%ne
graphs(1,1)%g%nnz=self%nee
graphs(1,1)%g%kr=>self%kee
graphs(1,1)%g%lc=>self%lee
!---
CALL create_matrix(new,graphs,tmp_vec,tmp_vec)
CALL tmp_vec%delete
DEALLOCATE(graphs(1,1)%g,tmp_vec)
DEBUG_STACK_POP
end subroutine afem_mat_create
!------------------------------------------------------------------------------
!> Create weight vector for FE representation
!------------------------------------------------------------------------------
subroutine ml_fem_vec_create(self,new,level,cache,native)
class(oft_ml_fem_type), intent(inout) :: self
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
!> Set the current level for a ML FE structure
!------------------------------------------------------------------------------
subroutine ml_fem_set_level(self,level)
class(oft_ml_fem_type), intent(inout) :: self
integer(i4), intent(in) :: level !< Desired level
DEBUG_STACK_PUSH
IF(level>self%nlevels.OR.level<=0)THEN
  WRITE(*,*)level,self%nlevels
  CALL oft_abort('Invalid FE level change requested', &
  'ml_fem_set_level',__FILE__)
END IF
!---Update level
self%level=level
self%current_level=>self%levels(self%level)%fe
self%abs_level=self%level
IF((self%level>self%blevel).AND.(self%blevel>0))self%abs_level=self%level-1
!---Set grid level
if(level<self%ml_mesh%mgdim)then
  call multigrid_level(self%ml_mesh,level)
else
  call multigrid_level(self%ml_mesh,self%ml_mesh%mgdim)
end if
DEBUG_STACK_POP
end subroutine ml_fem_set_level
!------------------------------------------------------------------------------
!> Destroy boundary multi-level FE object
!------------------------------------------------------------------------------
SUBROUTINE ml_fem_delete(self)
CLASS(oft_ml_fem_type), INTENT(inout) :: self
INTEGER(i4) :: i
DEBUG_STACK_PUSH
DO i=self%minlev,self%nlevels
  IF(ASSOCIATED(self%levels(i)%fe))THEN
    CALL self%levels(i)%fe%delete()
    DEALLOCATE(self%levels(i)%fe)
  END IF
  IF(ASSOCIATED(self%interp_matrices(i)%m))THEN
    CALL self%interp_matrices(i)%m%delete()
    DEALLOCATE(self%interp_matrices(i)%m)
  END IF
  IF(ASSOCIATED(self%interp_graphs(i)%g))THEN
    ! TODO: Perform actual deallocation
    DEALLOCATE(self%interp_graphs(i)%g)
  END IF
END DO
NULLIFY(self%ml_mesh,self%current_level)
self%nlevels=0
self%minlev=1
DEBUG_STACK_POP
END SUBROUTINE ml_fem_delete
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine ml_fe_vecspace_create(self,new,level,cache,native)
class(oft_ml_fe_vecspace), intent(inout) :: self
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
class(oft_ml_fe_vecspace), intent(inout) :: self
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
class(oft_ml_fe_vecspace), intent(inout) :: self
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
!------------------------------------------------------------------------------
!> Constructs a finite element representation on the associated surface mesh
!!
!! This subroutine is the equilivalent of @ref fem_base::fem_setup "fem_setup"
!! for trangular grids and @ref fem_base::oft_bfem_type. Generally this method
!! is used to construct a finite element representation for the boundary mesh,
!! however it may be used with arbitrary triangular grids
!------------------------------------------------------------------------------
subroutine bfem_setup(self,quad_order)
class(oft_bfem_type), intent(inout) :: self !< FE representation to construct
integer(i4), intent(in) :: quad_order !< Desired quadrature order
integer(i4) :: i,j,k,offset,stack
integer(i4) :: jeg,jsg,jee,jse,nefmax,js,je,jn
integer(i4) :: mycounts(3)
integer(i4), allocatable :: lcx(:,:),nr(:)
class(oft_bmesh), pointer :: mesh
DEBUG_STACK_PUSH
if(oft_debug_print(2))write(*,'(2A)')oft_indent,'Starting boundary FE setup'
CALL oft_increase_indent
!---
mesh=>self%mesh
IF(mesh%skip)THEN
  if(oft_debug_print(2))write(*,'(2A)')oft_indent,'Creating interaction graph: skip'
  if(oft_debug_print(2))write(*,'(2A)')oft_indent,'Creating cell mapping: skip'
  if(oft_debug_print(2))write(*,'(2A)')oft_indent,'Constructing MPI Linkage: skip'
  call bfem_global_linkage(self)
  CALL oft_decrease_indent
  DEBUG_STACK_POP
  RETURN
END IF
!---Set quadrature order
CALL mesh%quad_rule(quad_order, self%quad)
!---Set number of node pts
self%nnodes=self%quad%np*mesh%nc
!---Compute DOFs from geometric linkage
mycounts=(/mesh%np,mesh%ne,mesh%nc/)
self%ne=sum(self%gstruct*mycounts)
!---Compute boundary count from geometric linkage
mycounts=(/mesh%nbp,mesh%nbe,0/)
self%nbe=sum(self%gstruct*mycounts)
allocate(self%be(self%ne),self%lbe(self%nbe))
self%be=.FALSE.
!---
allocate(self%ce(self%ne))
self%ce=.FALSE.
!---Compute cell DOFs from geometric linkage
mycounts(1:2)=mesh%cell_np
mycounts(3)=1
self%nce=sum(self%gstruct*mycounts)
!---Count element to cell interactions
allocate(self%kec(self%ne+1))
stack=1
self%nbe=0
!---Point stacking
do i=1,mesh%np
  offset=(i-1)*self%gstruct(1)
  do j=1,self%gstruct(1)
    self%kec(j+offset)=stack
    stack=stack+mesh%kpc(i+1)-mesh%kpc(i)
    if(mesh%bp(i))then
      self%be(j+offset)=.TRUE.
      self%nbe=self%nbe+1
      self%lbe(self%nbe)=j+offset
    end if
  end do
end do
!---Edge stacking
do i=1,mesh%ne
  offset=(i-1)*self%gstruct(2)+self%gstruct(1)*mesh%np
  do j=1,self%gstruct(2)
    self%kec(j+offset)=stack
    stack=stack+(mesh%kec(i+1)-mesh%kec(i))
    ! if(mesh%lef(2,i)/=0)then
    !   stack=stack+2
    ! else
    !   stack=stack+1
    ! end if
    if(mesh%be(i))then
      self%be(j+offset)=.TRUE.
      self%nbe=self%nbe+1
      self%lbe(self%nbe)=j+offset
    end if
  end do
end do
!---Face stacking
do i=1,mesh%nc
  offset=(i-1)*self%gstruct(3)+self%gstruct(2)*mesh%ne+self%gstruct(1)*mesh%np
  do j=1,self%gstruct(3)
    self%kec(j+offset)=stack
    stack=stack+1
  end do
end do
IF(self%nbe/=SIZE(self%lbe))CALL oft_abort('Bad boundary count','bfem_setup',__FILE__)
self%kec(self%ne+1)=stack
self%nec=stack-1
!---Assemble element to cell interaction list
allocate(self%lec(self%nec))
nefmax=0
!$omp parallel private(offset,j,jee,jse,jeg,jsg) reduction(max:nefmax)
!---Point stacking
!$omp do
do i=1,mesh%np
  offset=(i-1)*self%gstruct(1)
  do j=1,self%gstruct(1)
    jee=self%kec(j+offset+1)-1; jse=self%kec(j+offset)
    jeg=mesh%kpc(i+1)-1; jsg=mesh%kpc(i)
    nefmax=MAX(nefmax,jee-jse)
    self%lec(jse:jee)=mesh%lpc(jsg:jeg)
  end do
end do
!---Edge stacking
!$omp do
do i=1,mesh%ne
  offset=(i-1)*self%gstruct(2)+self%gstruct(1)*mesh%np
  nefmax=MAX(nefmax,2)
  do j=1,self%gstruct(2)
    jee=self%kec(j+offset+1)-1; jse=self%kec(j+offset)
    jeg=mesh%kec(i+1)-1; jsg=mesh%kec(i)
    self%lec(jse:jee)=mesh%lec(jsg:jeg)
    ! if(mesh%lef(2,i)/=0)then
    !   self%lec(jse:jee)=mesh%lef(:,i)
    ! else
    !   self%lec(jse:jee)=mesh%lef(1,i)
    ! end if
  end do
end do
!---Face stacking
!$omp do
do i=1,mesh%nc
  offset=(i-1)*self%gstruct(3)+self%gstruct(2)*mesh%ne+self%gstruct(1)*mesh%np
  do j=1,self%gstruct(3)
    jse=self%kec(j+offset)
    self%lec(jse)=i
  end do
end do
!$omp end parallel
self%necmax=nefmax
if(oft_debug_print(2))write(*,'(2A)')oft_indent,'Creating interaction graph'
call afem_self_linkage(self)
if(oft_debug_print(2))write(*,'(2A)')oft_indent,'Creating cell mapping'
call bfem_nfdofs_map(self)
if(oft_debug_print(2))write(*,'(2A)')oft_indent,'Constructing MPI Linkage'
call bfem_global_linkage(self)
call afem_fill_lgraph(self,self,self%nee,self%kee,self%lee)
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine bfem_setup
!------------------------------------------------------------------------------
!> Destroy boundary FE object
!------------------------------------------------------------------------------
SUBROUTINE bfem_delete(self)
CLASS(oft_bfem_type), INTENT(inout) :: self
DEBUG_STACK_PUSH
CALL afem_delete(self)
IF(ASSOCIATED(self%parent))THEN
  self%parent%ne=0
  IF(ASSOCIATED(self%parent%le))DEALLOCATE(self%parent%le)
END IF
NULLIFY(self%mesh)
DEBUG_STACK_POP
END SUBROUTINE bfem_delete
!------------------------------------------------------------------------------
!> Compute FE global context and stitching information
!!
!! This subroutine is the equilivalent of @ref fem_base::fem_global_linkage
!! "fem_global_linkage" for triangular grids and @ref fem_base::oft_bfem_type
!------------------------------------------------------------------------------
subroutine bfem_global_linkage(self)
class(oft_bfem_type), intent(inout) :: self !< Finite element representation
type :: eout
  integer(i8), pointer, dimension(:) :: le => NULL()
end type eout
type(eout), allocatable, dimension(:) :: lerecv,lesend
LOGICAL :: has_parent
integer(i8), pointer, dimension(:) :: letmp,leout
integer(i8), allocatable, dimension(:) :: lsort
integer(i4), allocatable :: linktmp(:,:,:),isort(:),ncon(:)
integer(i4) :: neel,mycounts4(3)
integer(i4) :: kk,m,mm,l,request
integer(i4) :: nbetmp,nbemax
integer(i8) :: netmp,nemax
integer(i4) :: i,j,k,offset,offsetp,stack1,stack2
integer(i4) :: je,js,jee,jse,ierr
integer(i8) :: offsetg
integer(i8) :: mycounts(3)
integer(i4), allocatable, dimension(:) :: bpi,bei,bfi
class(oft_bmesh), pointer :: mesh
character(LEN=2) :: pltnum
#ifdef HAVE_MPI
integer(i4) :: stat(MPI_STATUS_SIZE)
#endif
DEBUG_STACK_PUSH
CALL oft_increase_indent
!---
mesh=>self%mesh
has_parent=ASSOCIATED(mesh%parent)
IF(mesh%skip)THEN
  allocate(self%linkage,self%global)
  CALL oft_init_seam(mesh,self%linkage)
  IF(.NOT.mesh%fullmesh)CALL oft_mpi_barrier(ierr) ! Wait for all processes
  self%linkage%skip=.TRUE.
  !---Determine global element count
  mycounts=(/mesh%global%np,mesh%global%ne,mesh%global%nc/)
  self%global%ne=sum(self%gstruct*mycounts)
  !---Determine parent element count
  IF(has_parent)THEN
    allocate(self%parent)
    mycounts4=(/mesh%parent%np,mesh%parent%ne,mesh%parent%nf/)
    self%parent%ne=sum(self%gstruct*mycounts4)
  END IF
  IF(mesh%fullmesh)THEN
    self%linkage%full=.TRUE.
  ELSE
    self%linkage%full=.FALSE.
  END IF
  !---
  allocate(self%linkage%kle(0:self%linkage%nproc_con+1))
  self%linkage%kle=0
  nbetmp=0_i4
  nbemax=oft_mpi_max(nbetmp)
  if(oft_debug_print(3))write(*,'(2A,I8)')oft_indent,'Max # of boundary elements',nbemax
  self%linkage%nbemax=nbemax
  !---Allocate temporary Send/Recv arrays
  ALLOCATE(lesend(self%linkage%nproc_con+1),lerecv(self%linkage%nproc_con+1))
  allocate(lesend(1)%le(nbemax))
  IF(.NOT.mesh%fullmesh)THEN
    do i=1,self%linkage%nproc_con
      allocate(lerecv(i)%le(nbemax))
    end do
  END IF
  !---Point dummy output array to Send array
  leout=>lesend(1)%le
  leout = -1
  IF(.NOT.mesh%fullmesh)THEN
#ifdef HAVE_MPI
    CALL oft_mpi_barrier(ierr) ! Wait for all processes
    !---Point dummy Send arrays to main Send array
    do i=2,self%linkage%nproc_con
      lesend(i)%le=>lesend(1)%le
    end do
    !---Create Send and Recv calls
    do j=1,self%linkage%nproc_con
      CALL MPI_ISEND(lesend(j)%le,nbemax,OFT_MPI_I8,self%linkage%proc_con(j),1,oft_env%COMM,self%linkage%send_reqs(j),ierr)
      IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','bfem_global_linkage',__FILE__)
      CALL MPI_IRECV(lerecv(j)%le,nbemax,OFT_MPI_I8,self%linkage%proc_con(j),1,oft_env%COMM,self%linkage%recv_reqs(j),ierr)
      IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','bfem_global_linkage',__FILE__)
    end do
    !---Loop over each connected processor
    do while(.TRUE.)
      IF(oft_mpi_check_reqs(self%linkage%nproc_con,self%linkage%recv_reqs))EXIT ! All recieves have been processed
      CALL oft_mpi_waitany(self%linkage%nproc_con,self%linkage%recv_reqs,j,ierr) ! Wait for completed recieve
      IF(ierr/=0)CALL oft_abort('Error in MPI_WAITANY','bfem_global_linkage',__FILE__)
    end do
#else
CALL oft_abort("Distributed linkage requires MPI","bfem_global_linkage",__FILE__)
#endif
  END IF
  !---Condense linkage to sparse rep
  self%linkage%nle=SUM(self%linkage%kle)
  self%linkage%kle(self%linkage%nproc_con+1)=self%linkage%nle+1
  do i=self%linkage%nproc_con,0,-1 ! cumulative unique point linkage count
    self%linkage%kle(i)=self%linkage%kle(i+1)-self%linkage%kle(i)
  end do
  !---
  if(self%linkage%kle(0)/=1)call oft_abort('Bad element linkage count(skip)','bfem_global_linkage',__FILE__)
  !---
  allocate(self%linkage%send(0:self%linkage%nproc_con),self%linkage%recv(0:self%linkage%nproc_con))
  do i=0,self%linkage%nproc_con
    !---Allocate permanent stitching arrays
    self%linkage%send(i)%n=0; self%linkage%recv(i)%n=0
  end do
  IF(.NOT.mesh%fullmesh)THEN
#ifdef HAVE_MPI
    CALL oft_mpi_waitall(self%linkage%nproc_con,self%linkage%send_reqs,ierr) ! Wait for all sends to complete
    IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','bfem_global_linkage',__FILE__)
    !---Deallocate temporary work arrays
    do i=1,self%linkage%nproc_con
      deallocate(lerecv(i)%le)
    end do
#else
CALL oft_abort("Distributed linkage requires MPI","bfem_global_linkage",__FILE__)
#endif
  END IF
  CALL oft_stitch_check(self%linkage)
  DEALLOCATE(lesend(1)%le)
  DEALLOCATE(lesend,lerecv)
  !---Construct map
  ALLOCATE(self%map)
  self%map%offset=0
  self%map%n=0
  self%map%ng=self%global%ne
  self%map%nslice=0
  ALLOCATE(self%map%slice(1))
  self%map%nslice=0
  self%map%slice=0
  self%map%lge=>self%global%le
  self%map%gbe=>self%global%gbe
  CALL oft_decrease_indent
  DEBUG_STACK_POP
  RETURN
END IF
!---Set MPI transfer arrays
allocate(self%linkage,self%global)
CALL oft_init_seam(mesh,self%linkage)
!
IF(has_parent)ALLOCATE(self%parent)
IF(.NOT.mesh%fullmesh)CALL oft_mpi_barrier(ierr) ! Wait for all processes
!---Determine global element count
mycounts=(/mesh%global%np,mesh%global%ne,mesh%global%nc/)
self%global%ne=sum(self%gstruct*mycounts)
allocate(self%global%le(self%ne))
!---Determine parent element count
IF(has_parent)THEN
  mycounts4=(/mesh%parent%np,mesh%parent%ne,mesh%parent%nf/)
  self%parent%ne=sum(self%gstruct*mycounts4)
  allocate(self%parent%le(self%ne))
END IF
!---
IF(mesh%fullmesh)THEN
  self%linkage%full=.TRUE.
ELSE
  self%linkage%full=.FALSE.
END IF
allocate(self%global%gbe(self%ne))
self%global%gbe=.FALSE.
!$omp parallel private(offset,offsetg,offsetp,j)
!---Point stacking
!$omp do
do i=1,mesh%np
  offset=(i-1)*self%gstruct(1)
  offsetg=(mesh%global%lp(i)-1)*self%gstruct(1)
  IF(has_parent)offsetp=(mesh%parent%lp(i)-1)*self%gstruct(1)
  do j=1,self%gstruct(1)
    self%global%le(j+offset)=j+offsetg
    IF(has_parent)THEN
      self%parent%le(j+offset)=j+offsetp
    ELSE
      self%global%gbe(j+offset)=mesh%global%gbp(i)
    END IF
  end do
end do
!$omp end do nowait
!---Edge stacking
!$omp do
do i=1,mesh%ne
  offset = (i-1)*self%gstruct(2)+self%gstruct(1)*mesh%np
  offsetg = (ABS(mesh%global%le(i))-1)*self%gstruct(2)+self%gstruct(1)*mesh%global%np
  IF(has_parent)offsetp = (ABS(mesh%parent%le(i))-1)*self%gstruct(2)+self%gstruct(1)*mesh%parent%np
  do j=1,self%gstruct(2)
    self%global%le(j+offset)=j+offsetg
    IF(has_parent)THEN
      self%parent%le(j+offset)=j+offsetp
    ELSE
      self%global%gbe(j+offset)=mesh%global%gbe(i)
    END IF
  end do
end do
!$omp end do nowait
!---Cell stacking
!$omp do
do i=1,mesh%nc
  offset = (i-1)*self%gstruct(3)+self%gstruct(2)*mesh%ne+self%gstruct(1)*mesh%np
  offsetg = (ABS(mesh%global%lc(i))-1)*self%gstruct(3)+self%gstruct(2)*mesh%global%ne &
  + self%gstruct(1)*mesh%global%np
  IF(has_parent)offsetp = (ABS(mesh%parent%lf(i))-1)*self%gstruct(3)+self%gstruct(2)*mesh%parent%ne &
    + self%gstruct(1)*mesh%parent%np
  do j=1,self%gstruct(3)
    self%global%le(j+offset)=j+offsetg
    IF(has_parent)THEN
      self%parent%le(j+offset)=j+offsetp
    ! ELSE
    !   self%global%gbe(j+offset)=mesh%global%gbc(i)
    END IF
  end do
end do
!$omp end parallel
ALLOCATE(linktmp(2,5*self%nbe,0:self%linkage%nproc_con))
linktmp = 0
ALLOCATE(ncon(0:self%linkage%nproc_con))
ncon=0
ALLOCATE(self%linkage%kle(0:self%linkage%nproc_con+1)) ! Allocate point linkage arrays
self%linkage%kle=0
self%linkage%nbe=self%nbe
self%linkage%lbe=>self%lbe
self%linkage%be=>self%be
linktmp=0
!---Determine maximum boundary element count
nbetmp=self%nbe
nbemax=self%nbe
IF(.NOT.mesh%fullmesh)nbemax=oft_mpi_max(nbetmp)
if(oft_debug_print(3))write(*,'(2A,I8)')oft_indent,'Max # of boundary elements',nbemax
self%linkage%nbemax=nbemax
!---Allocate temporary Send/Recv arrays
ALLOCATE(lesend(self%linkage%nproc_con+1),lerecv(self%linkage%nproc_con+1))
allocate(lesend(1)%le(nbemax))
IF(.NOT.mesh%fullmesh)THEN
  do i=1,self%linkage%nproc_con
    allocate(lerecv(i)%le(nbemax))
  end do
END IF
!---Point dummy output array to Send array
leout=>lesend(1)%le
leout=0
!$omp parallel do
do i=1,self%nbe
  leout(i) = self%global%le(self%lbe(i)) ! Populate with global point index
end do
IF(.NOT.mesh%fullmesh)THEN
#ifdef HAVE_MPI
  CALL oft_mpi_barrier(ierr) ! Wait for all processes
  !---Point dummy Send arrays to main Send array
  do i=2,self%linkage%nproc_con
    lesend(i)%le=>lesend(1)%le
  end do
  !---Create Send and Recv calls
  do j=1,self%linkage%nproc_con
    CALL MPI_ISEND(lesend(j)%le,nbemax,OFT_MPI_I8,self%linkage%proc_con(j),1,oft_env%COMM,self%linkage%send_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_ISEND','bfem_global_linkage',__FILE__)
    CALL MPI_IRECV(lerecv(j)%le,nbemax,OFT_MPI_I8,self%linkage%proc_con(j),1,oft_env%COMM,self%linkage%recv_reqs(j),ierr)
    IF(ierr/=0)CALL oft_abort('Error in MPI_IRECV','bfem_global_linkage',__FILE__)
  end do
  !---Loop over each connected processor
  do while(.TRUE.)
    IF(oft_mpi_check_reqs(self%linkage%nproc_con,self%linkage%recv_reqs))EXIT ! All recieves have been processed
    CALL oft_mpi_waitany(self%linkage%nproc_con,self%linkage%recv_reqs,j,ierr) ! Wait for completed recieve
    IF(ierr/=0)CALL oft_abort('Error in MPI_WAITANY','bfem_global_linkage',__FILE__)
    letmp=>lerecv(j)%le ! Point dummy input array to current Recv array
    !---Determine location of boundary points on other processors
    neel=0
    !$omp parallel do private(mm) reduction(+:neel)
    do m=1,self%nbe ! Loop over boundary points
      do mm=1,nbemax ! Loop over input points
        if(leout(m)==letmp(mm))then ! Found match
          !linktmp(j,m)=mm ! Populate linkage array
          !$omp critical
          ncon(j)=ncon(j)+1
          linktmp(:,ncon(j),j)=(/mm,m/)
          !$omp end critical
          neel=neel+1
          !exit
        end if
        if(letmp(mm)==0)exit ! End of input points
      end do
    end do
    self%linkage%kle(j)=neel
  end do
#else
CALL oft_abort("Distributed linkage requires MPI","bfem_global_linkage",__FILE__)
#endif
END IF
!---Local connectivity
!---Determine location of boundary points on other processors
neel=0
!$omp parallel do private(mm) reduction(+:neel)
do m=1,self%nbe ! Loop over boundary points
  do mm=1,self%nbe ! Loop over input points
    if(m==mm)CYCLE
    if(leout(m)==leout(mm))then ! Found match
      !$omp critical
      ncon(0)=ncon(0)+1
      linktmp(:,ncon(0),0)=(/mm,m/)
      !$omp end critical
      neel=neel+1
    end if
    if(leout(mm)==0)exit ! End of input points
  end do
end do
self%linkage%kle(0)=neel
!---Condense linkage to sparse rep
self%linkage%nle=sum(self%linkage%kle)
self%linkage%kle(self%linkage%nproc_con+1)=self%linkage%nle+1
do i=self%linkage%nproc_con,0,-1 ! cumulative unique point linkage count
  self%linkage%kle(i)=self%linkage%kle(i+1)-self%linkage%kle(i)
end do
!---
if(self%linkage%kle(0)/=1)call oft_abort('Bad element linkage count','bfem_global_linkage',__FILE__)
!---
allocate(self%linkage%lle(2,self%linkage%nle))
allocate(self%linkage%send(0:self%linkage%nproc_con),self%linkage%recv(0:self%linkage%nproc_con))
!---
!$omp parallel private(j,m,lsort,isort)
ALLOCATE(lsort(MAXVAL(ncon)),isort(MAXVAL(ncon)))
!$omp do
do i=0,self%linkage%nproc_con
  !---
  DO j=1,ncon(i)
    m=linktmp(2,j,i)
    lsort(j)=m
    isort(j)=j
    !---
    IF(i>0)THEN
      IF(self%linkage%proc_con(i)<oft_env%rank)lsort(j)=linktmp(1,j,i)
    END IF
  END DO
  !---
  CALL sort_array(lsort,isort,ncon(i))
  DO m=0,ncon(i)-1
    self%linkage%lle(:,m+self%linkage%kle(i)) = &
    (/linktmp(2,isort(m+1),i),linktmp(1,isort(m+1),i)/)
  END DO
  !---Allocate permanent stitching arrays
  self%linkage%send(i)%n=ncon(i); self%linkage%recv(i)%n=ncon(i)
  allocate(self%linkage%send(i)%v(self%linkage%send(i)%n))
  allocate(self%linkage%recv(i)%v(self%linkage%recv(i)%n))
end do
DEALLOCATE(lsort,isort)
!$omp end parallel
!---Propogate ownership
allocate(bpi(mesh%np),bei(mesh%ne))
CALL get_inverse_map(mesh%lbp,mesh%nbp,bpi,mesh%np)
CALL get_inverse_map(mesh%lbe,mesh%nbe,bei,mesh%ne)
allocate(self%linkage%leo(self%nbe))
mycounts(1)=self%gstruct(1)*mesh%np
mycounts(2)=mycounts(1)+self%gstruct(2)*mesh%ne
!$omp parallel do private(k,j)
do i=1,self%nbe
  k=self%lbe(i)
  if(k<=mycounts(1))then
    j=FLOOR(REAL(k-1,8)/self%gstruct(1))+1
    self%linkage%leo(i)=mesh%pstitch%leo(bpi(j))
  else if(k>mycounts(1).AND.k<=mycounts(2))then
    j=FLOOR(REAL(k-mycounts(1)-1,8)/self%gstruct(2))+1
    self%linkage%leo(i)=mesh%estitch%leo(bei(j))
  else
    call oft_abort('Invalid boundary element','bfem_global_linkage',__FILE__)
  end if
end do
!---Synchronize
IF(.NOT.mesh%fullmesh)THEN
  CALL oft_mpi_waitall(self%linkage%nproc_con,self%linkage%send_reqs,ierr) ! Wait for all sends to complete
  IF(ierr/=0)CALL oft_abort('Error in MPI_WAITALL','bfem_global_linkage',__FILE__)
  !---Deallocate temporary work arrays
  do i=1,self%linkage%nproc_con
    deallocate(lerecv(i)%le)
  end do
END IF
CALL oft_stitch_check(self%linkage)
deallocate(linktmp,ncon)
deallocate(bpi,bei)
DEALLOCATE(lesend(1)%le)
DEALLOCATE(lesend,lerecv)
!---Construct map
ALLOCATE(self%map)
self%map%offset=0
self%map%n=self%ne
self%map%ng=self%global%ne
self%map%nslice=0
DO i=1,self%ne
  IF(self%linkage%be(i))CYCLE
  self%map%nslice=self%map%nslice+1
END DO
DO i=1,self%linkage%nbe
  IF(.NOT.self%linkage%leo(i))CYCLE
  self%map%nslice=self%map%nslice+1
END DO
ALLOCATE(self%map%slice(self%map%nslice))
self%map%nslice=0
DO i=1,self%ne
  IF(self%linkage%be(i))CYCLE
  self%map%nslice=self%map%nslice+1
  self%map%slice(self%map%nslice)=i
END DO
DO i=1,self%linkage%nbe
  IF(.NOT.self%linkage%leo(i))CYCLE
  self%map%nslice=self%map%nslice+1
  self%map%slice(self%map%nslice)=self%linkage%lbe(i)
END DO
self%map%lge=>self%global%le
self%map%gbe=>self%global%gbe
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine bfem_global_linkage
!------------------------------------------------------------------------------
!> Retrieve the indices of elements beloning to a given face
!------------------------------------------------------------------------------
subroutine bfem_ncdofs(self,cell,dofs)
class(oft_bfem_type), intent(in) :: self !< Finite element representation
integer(i4), intent(in) :: cell !< Desired cell in mesh
integer(i4), intent(inout) :: dofs(:) !< Indices of face elements [self%nfe]
integer(i4) :: foffset,eoffset,boffset,i,j
DEBUG_STACK_PUSH
!---Point DOFs
do i=1,self%mesh%cell_np
  foffset=(i-1)*self%gstruct(1)
  eoffset=(self%mesh%lc(i,cell)-1)*self%gstruct(1)
  do j=1,self%gstruct(1)
    dofs(j+foffset)=j+eoffset
  end do
end do
boffset=self%gstruct(1)*self%mesh%cell_np
!---Edge DOFs
do i=1,self%mesh%cell_ne
  foffset=(i-1)*self%gstruct(2)+boffset
  eoffset=(abs(self%mesh%lce(i,cell))-1)*self%gstruct(2)+self%gstruct(1)*self%mesh%np
  do j=1,self%gstruct(2)
    dofs(j+foffset)=j+eoffset
  end do
end do
boffset=boffset+self%gstruct(2)*self%mesh%cell_ne
!---Face DOFs
do i=1,self%gstruct(3)
  foffset=boffset
  eoffset=self%gstruct(2)*self%mesh%ne+self%gstruct(1)*self%mesh%np
  dofs(i+foffset)=i+(cell-1)*self%gstruct(3)+eoffset
end do
DEBUG_STACK_POP
end subroutine bfem_ncdofs
!------------------------------------------------------------------------------
!> Construct face element mapping
!!
!! Sets up the structure @ref fem_base::oft_bfem_type::fmap "fmap", which defines
!! the local type, index, and geometric linkage of DOFs on a face
!------------------------------------------------------------------------------
subroutine bfem_nfdofs_map(self)
class(oft_bfem_type), intent(inout) :: self !< Finite element representation
integer(i4) :: foffset,boffset,i,j
DEBUG_STACK_PUSH
!---
allocate(self%cmap(self%nce))
!---Point DOFs
do i=1,self%mesh%cell_np
  foffset=(i-1)*self%gstruct(1)
  do j=1,self%gstruct(1)
    self%cmap(j+foffset)%type=1
    self%cmap(j+foffset)%el=i
    self%cmap(j+foffset)%ind=j
  end do
end do
boffset=self%gstruct(1)*self%mesh%cell_np
!---Edge DOFs
do i=1,self%mesh%cell_ne
  foffset=(i-1)*self%gstruct(2)+boffset
  do j=1,self%gstruct(2)
    self%cmap(j+foffset)%type=2
    self%cmap(j+foffset)%el=i
    self%cmap(j+foffset)%ind=j
  end do
end do
boffset=boffset+self%gstruct(2)*self%mesh%cell_ne
!---Face DOFs
do i=1,self%gstruct(3)
  foffset=i+boffset
  self%cmap(foffset)%type=3
  self%cmap(foffset)%el=1
  self%cmap(foffset)%ind=i
end do
DEBUG_STACK_POP
end subroutine bfem_nfdofs_map
end module fem_base
