!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_la_utils.F90
!
!> Matrix and vector management routines
!!
!! @authors Chris Hansen
!! @date December 2012
!! @ingroup doxy_oft_lin_alg
!---------------------------------------------------------------------------
MODULE oft_la_utils
USE oft_local
USE oft_base
USE oft_sort, ONLY: sort_array, search_array
USE oft_stitching, ONLY: oft_seam, seam_list, oft_global_stitch, &
  oft_stitch_check
USE oft_la_base, ONLY: oft_vector, oft_cvector, oft_map, map_list, &
  oft_matrix, oft_matrix_ptr, oft_cmatrix, oft_cmatrix_ptr, &
  oft_graph, oft_graph_ptr, oft_matrix_map
USE oft_native_la, ONLY: oft_native_vector, oft_native_cvector, &
  oft_native_matrix, oft_native_cmatrix, native_cvector_cast
#ifdef HAVE_PETSC
USE oft_petsc_la, ONLY: oft_petsc_vector, oft_petsc_vector_cast, oft_petsc_matrix, &
  oft_petsc_matrix_cast
#include "petsc/finclude/petscmat.h"
#undef IS
#undef Mat
use petscmat
#endif
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------
!> Create a new vector by combining a set of vectors.
!!
!! @param[in,out] vec Resulting vector
!! @param[in,out] stitch_info Array of seam structures
!! @param[in,out] maps Mapping from sub-vectors into full vector
!! @param[in] native Force native representation
!---------------------------------------------------------------------------
INTERFACE create_vector
  MODULE PROCEDURE create_vector_real
  MODULE PROCEDURE create_vector_comp
END INTERFACE create_vector
!---------------------------------------------------------------------------
!> Create a matrix using a set of non-overlapping graphs.
!!
!! Native and PETSc (real only)s matrices are supported.
!!
!! @param[in,out] mat Resulting matrix
!! @param[in] ingraphs Array of graphs representing submatrices
!! @param[in] row_vec Vector representing matrix rows
!! @param[in] col_vec Vector representing matrix columns
!! @param[in] native Force native representation [optional]
!---------------------------------------------------------------------------
INTERFACE create_matrix
  MODULE PROCEDURE create_matrix_real
  MODULE PROCEDURE create_matrix_comp
END INTERFACE create_matrix
!---------------------------------------------------------------------------
!> Combine a set of non-overlapping sub-matrices into a single matrix.
!!
!! Native and PETSc (real only) matrices are supported. The matrix should be created using
!! \ref create_matrix before this subroutine is called.
!!
!! @param[in] mats Array of sub-matrices (nr,nc)
!! @param[in] nr Number of row blocks
!! @param[in] nc Number of column blocks
!! @param[in,out] mat Resulting matrix
!---------------------------------------------------------------------------
INTERFACE combine_matrices
  MODULE PROCEDURE combine_matrices_real
  MODULE PROCEDURE combine_matrices_comp
END INTERFACE combine_matrices
CONTAINS
!------------------------------------------------------------------------------
!> Real implementation for \ref create_vector
!------------------------------------------------------------------------------
SUBROUTINE create_vector_real(vec,stitch_info,maps,native)
CLASS(oft_vector), POINTER, INTENT(inout) :: vec
TYPE(seam_list), INTENT(inout) :: stitch_info(:)
TYPE(map_list), INTENT(inout) :: maps(:)
LOGICAL, OPTIONAL, INTENT(in) :: native
!---
LOGICAL :: force_native
INTEGER(i4) :: i,offset,soffset,nblocks
INTEGER(i8) :: offsetg
DEBUG_STACK_PUSH
force_native=.FALSE.
IF(PRESENT(native))force_native=native
!---Determine vector type
IF(use_petsc.AND.(.NOT.force_native))THEN
#ifdef HAVE_PETSC
  ALLOCATE(oft_petsc_vector::vec)
#else
  CALL oft_abort("Not compiled with PETSc","create_vector_real",__FILE__)
#endif
ELSE
  ALLOCATE(oft_native_vector::vec)
END IF
!---Common setup
nblocks=SIZE(stitch_info)
vec%nblocks=nblocks
vec%n=0
vec%nslice=0
vec%ng=0
!---Allocate and setup map structure
ALLOCATE(vec%map(nblocks))
offset=0
soffset=0
offsetg=0
DO i=1,nblocks
  vec%map(i)%per=maps(i)%m%per
  vec%map(i)%offset=offset
  vec%map(i)%soffset=soffset
  vec%map(i)%offsetg=offsetg
  vec%map(i)%n=maps(i)%m%n
  vec%map(i)%ng=maps(i)%m%ng
  vec%map(i)%nslice=maps(i)%m%nslice
  vec%map(i)%slice=>maps(i)%m%slice
  vec%map(i)%lge=>maps(i)%m%lge
  offset=offset+vec%map(i)%n
  soffset=soffset+vec%map(i)%nslice
  offsetg=offsetg+vec%map(i)%ng
  !---
  vec%n=vec%n+vec%map(i)%n
  vec%nslice=vec%nslice+vec%map(i)%nslice
  vec%ng=vec%ng+vec%map(i)%ng
END DO
ALLOCATE(vec%stitch_info)
CALL condense_stitch(stitch_info,vec%map,vec%stitch_info)
!---Setup for given backend
SELECT TYPE(this=>vec)
  TYPE is(oft_native_vector)
    CALL setup_native_vec(this)
    IF(oft_debug_print(2))WRITE(*,'(2X,A)')'Native real vector created'
#ifdef HAVE_PETSC
  TYPE IS(oft_petsc_vector)
    CALL setup_petsc_vec(this)
    IF(oft_debug_print(2))WRITE(*,'(2X,A)')'PETSc real vector created'
#endif
  CLASS DEFAULT
    CALL oft_abort('Error in vector allocation.','create_vector_real',__FILE__)
END SELECT
DEBUG_STACK_POP
CONTAINS
!
SUBROUTINE setup_native_vec(this)
TYPE(oft_native_vector), INTENT(inout) :: this
INTEGER(i4) :: i
!---Create vector storage
ALLOCATE(this%v(this%n))
this%v=0.d0
!---Creating interior list
this%stitch_info%nie=this%n-this%stitch_info%nbe
ALLOCATE(this%stitch_info%lie(this%stitch_info%nie))
this%stitch_info%nie=0
DO i=1,this%n
  IF(this%stitch_info%be(i))CYCLE
  this%stitch_info%nie=this%stitch_info%nie+1
  this%stitch_info%lie(this%stitch_info%nie)=i
END DO
END SUBROUTINE setup_native_vec
!
#ifdef HAVE_PETSC
SUBROUTINE setup_petsc_vec(this)
TYPE(oft_petsc_vector), INTENT(inout) :: this
INTEGER(i4) :: i,j,offset,offset1,ierr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: is0
REAL(r8), ALLOCATABLE, DIMENSION(:) :: vtmp
TYPE(tis) :: istemp
IF(this%stitch_info%full)THEN
  CALL VecCreate(PETSC_COMM_SELF,this%v,ierr)
  CALL VecSetType(this%v,VECSEQ,ierr)
  CALL VecSetSizes(this%v,this%nslice,PETSC_DETERMINE,ierr)
ELSE
  CALL VecCreate(oft_env%COMM,this%v,ierr)
  CALL VecSetType(this%v,VECMPI,ierr)
  CALL VecSetSizes(this%v,this%nslice,PETSC_DETERMINE,ierr)
END IF
CALL VecCreate(PETSC_COMM_SELF,this%vloc,ierr)
CALL VecSetType(this%vloc,VECSEQ,ierr)
CALL VecSetSizes(this%vloc,this%n,PETSC_DETERMINE,ierr)
!---
CALL VecGetOwnershipRange(this%v,offset1,offset,ierr)
!---
ALLOCATE(is0(this%n))
offset=0
DO i=1,nblocks
  ALLOCATE(vtmp(this%map(i)%n))
  DO j=1,this%map(i)%nslice
    vtmp(this%map(i)%slice(j))=REAL(offset1+this%map(i)%soffset+j,8)
  END DO
  CALL oft_global_stitch(stitch_info(i)%s,vtmp,0)
  is0(this%map(i)%offset+1:this%map(i)%offset+this%map(i)%n)=INT(vtmp-1,4)
  DEALLOCATE(vtmp)
END DO
!---
CALL ISCreateGeneral(PETSC_COMM_SELF,this%n,is0,PETSC_COPY_VALUES,this%lis,ierr)
CALL ISCreateStride(PETSC_COMM_SELF,this%n,0,1,istemp,ierr)
CALL VecScatterCreate(this%v,this%lis,this%vloc,istemp,this%get,ierr)
CALL ISDestroy(istemp,ierr)
DEALLOCATE(is0)
END SUBROUTINE setup_petsc_vec
#endif
END SUBROUTINE create_vector_real
!------------------------------------------------------------------------------
!> Complex implementation for \ref create_vector
!------------------------------------------------------------------------------
SUBROUTINE create_vector_comp(vec,stitch_info,maps,native)
CLASS(oft_cvector), POINTER, INTENT(inout) :: vec
TYPE(seam_list), INTENT(inout) :: stitch_info(:)
TYPE(map_list), INTENT(inout) :: maps(:)
LOGICAL, OPTIONAL, INTENT(in) :: native
!---
LOGICAL :: force_native
INTEGER(i4) :: i,offset,soffset,nblocks
INTEGER(i8) :: offsetg
DEBUG_STACK_PUSH
force_native=.FALSE.
IF(PRESENT(native))force_native=native
!---Determine vector type
IF(use_petsc.AND.(.NOT.force_native))THEN
#ifdef HAVE_PETSC
  ! ALLOCATE(oft_petsc_cvector::vec)
  CALL oft_abort("Complex vectors not supported with PETSc","create_vector_comp",__FILE__)
#else
  CALL oft_abort("Not compiled with PETSc","create_vector_comp",__FILE__)
#endif
ELSE
  ALLOCATE(oft_native_cvector::vec)
END IF
!---Common setup
nblocks=SIZE(stitch_info)
vec%nblocks=nblocks
vec%n=0
vec%nslice=0
vec%ng=0
!---Allocate and setup map structure
ALLOCATE(vec%map(nblocks))
offset=0
soffset=0
offsetg=0
DO i=1,nblocks
  vec%map(i)%per=maps(i)%m%per
  vec%map(i)%offset=offset
  vec%map(i)%soffset=soffset
  vec%map(i)%offsetg=offsetg
  vec%map(i)%n=maps(i)%m%n
  vec%map(i)%ng=maps(i)%m%ng
  vec%map(i)%nslice=maps(i)%m%nslice
  vec%map(i)%slice=>maps(i)%m%slice
  vec%map(i)%lge=>maps(i)%m%lge
  offset=offset+vec%map(i)%n
  soffset=soffset+vec%map(i)%nslice
  offsetg=offsetg+vec%map(i)%ng
  !---
  vec%n=vec%n+vec%map(i)%n
  vec%nslice=vec%nslice+vec%map(i)%nslice
  vec%ng=vec%ng+vec%map(i)%ng
END DO
ALLOCATE(vec%stitch_info)
CALL condense_stitch(stitch_info,vec%map,vec%stitch_info)
!---Setup for given backend
SELECT TYPE(this=>vec)
  TYPE is(oft_native_cvector)
    CALL setup_native_vec(this)
    IF(oft_debug_print(2))WRITE(*,'(2X,A)')'Native complex vector created'
  CLASS DEFAULT
    CALL oft_abort('Error in vector allocation.','create_vector_comp',__FILE__)
END SELECT
DEBUG_STACK_POP
CONTAINS
SUBROUTINE setup_native_vec(this)
TYPE(oft_native_cvector), INTENT(inout) :: this
INTEGER(i4) :: i
!---Create vector storage
ALLOCATE(this%v(this%n))
this%v=(0.d0,0.d0)
!---Creating interior list
this%stitch_info%nie=this%n-this%stitch_info%nbe
ALLOCATE(this%stitch_info%lie(this%stitch_info%nie))
this%stitch_info%nie=0
DO i=1,this%n
  IF(this%stitch_info%be(i))CYCLE
  this%stitch_info%nie=this%stitch_info%nie+1
  this%stitch_info%lie(this%stitch_info%nie)=i
END DO
END SUBROUTINE setup_native_vec
END SUBROUTINE create_vector_comp
!------------------------------------------------------------------------------
!> Combine seam information for a set of vectors
!------------------------------------------------------------------------------
SUBROUTINE condense_stitch(stitch_info,map,stitcher)
TYPE(seam_list), INTENT(in) :: stitch_info(:) !< Array of seam structures
type(oft_map), INTENT(in) :: map(:) !< Mapping from sub-vectors into full vector
TYPE(oft_seam), INTENT(inout) :: stitcher !< Resulting seam structure for full vector
INTEGER(i4) :: i,j,k,m,nblocks,offset1,offset2
DEBUG_STACK_PUSH
!---
IF(oft_debug_print(3))WRITE(*,'(4X,A)')'Condensing seam structures'
nblocks=SIZE(stitch_info)
!---
stitcher%nproc_con=stitch_info(1)%s%nproc_con
stitcher%proc_split=stitch_info(1)%s%proc_split
stitcher%proc_con=>stitch_info(1)%s%proc_con
stitcher%send_reqs=>stitch_info(1)%s%send_reqs
stitcher%recv_reqs=>stitch_info(1)%s%recv_reqs
!---Condense stitching 
stitcher%full=.TRUE.
stitcher%skip=.TRUE.
stitcher%nbe=0
offset1=0
DO i=1,nblocks
  stitcher%full=stitcher%full.AND.stitch_info(i)%s%full
  stitcher%skip=stitcher%skip.AND.stitch_info(i)%s%skip
  stitcher%nbe=stitcher%nbe+stitch_info(i)%s%nbe
  IF(stitcher%nproc_con/=stitch_info(i)%s%nproc_con)CALL oft_abort("Inconsistent number of processor connections","condense_stitch",__FILE__)
  IF(stitcher%proc_split/=stitch_info(i)%s%proc_split)CALL oft_abort("Inconsistent processor connectivity split","condense_stitch",__FILE__)
  IF(stitcher%nproc_con>0)THEN
    IF(.NOT.ALL(stitcher%proc_con==stitch_info(i)%s%proc_con))CALL oft_abort("Inconsistent processor connectivity","condense_stitch",__FILE__)
  END IF
  !---
  offset1=offset1+map(i)%n
END DO
!---
ALLOCATE(stitcher%be(offset1))
stitcher%be=.FALSE.
!---
ALLOCATE(stitcher%lbe(stitcher%nbe))
ALLOCATE(stitcher%leo(stitcher%nbe))
offset1=0
offset2=0
DO i=1,nblocks
  DO j=1,map(i)%n
    stitcher%be(map(i)%offset+j)=stitch_info(i)%s%be(j)
  END DO
  DO j=1,stitch_info(i)%s%nbe
    stitcher%lbe(offset1+j)=map(i)%offset+stitch_info(i)%s%lbe(j)
    stitcher%leo(offset1+j)=stitch_info(i)%s%leo(j)
  END DO
  offset1=offset1+stitch_info(i)%s%nbe
END DO
!---
allocate(stitcher%kle(0:stitcher%nproc_con+1)) ! Allocate point linkage arrays
stitcher%kle=0
DO i=1,nblocks
  DO j=0,stitcher%nproc_con
    stitcher%kle(j)=stitcher%kle(j)+stitch_info(i)%s%kle(j+1)-stitch_info(i)%s%kle(j)
  END DO
END DO
!---Condense linkage to sparse rep
stitcher%nle=SUM(stitcher%kle)
stitcher%kle(stitcher%nproc_con+1)=stitcher%nle+1
do i=stitcher%nproc_con,0,-1 ! cumulative unique point linkage count
  stitcher%kle(i)=stitcher%kle(i+1)-stitcher%kle(i)
end do
if(stitcher%kle(0)/=1)call oft_abort('Bad element linkage count','fem_global_linkage',__FILE__)
!---Construct seam lists
allocate(stitcher%lle(2,stitcher%nle))
allocate(stitcher%send(0:stitcher%nproc_con),stitcher%recv(0:stitcher%nproc_con))
DO j=1,stitcher%nproc_con
  m=stitcher%kle(j)
  offset1=0
  stitcher%nbemax=0
  DO i=1,nblocks
    DO k=stitch_info(i)%s%kle(j),stitch_info(i)%s%kle(j+1)-1
      stitcher%lle(:,m)=(/stitch_info(i)%s%lle(1,k)+offset1,-1/)
      m=m+1
    END DO
    offset1=offset1+stitch_info(i)%s%nbe
    stitcher%nbemax=stitcher%nbemax+stitch_info(i)%s%nbemax
  END DO
  !---Allocate permanent stitching arrays
  stitcher%send(j)%n=m-stitcher%kle(j); stitcher%recv(j)%n=m-stitcher%kle(j)
  allocate(stitcher%send(j)%v(stitcher%send(j)%n))
  allocate(stitcher%recv(j)%v(stitcher%recv(j)%n))
END DO
!---Construct self-linkage
m=stitcher%kle(0)
offset1=0
stitcher%nbemax=0
DO i=1,nblocks
  DO k=stitch_info(i)%s%kle(0),stitch_info(i)%s%kle(1)-1
    stitcher%lle(:,m)=(/stitch_info(i)%s%lle(1,k)+offset1,stitch_info(i)%s%lle(2,k)+offset1/)
    m=m+1
  END DO
  offset1=offset1+stitch_info(i)%s%nbe
  stitcher%nbemax=stitcher%nbemax+stitch_info(i)%s%nbemax
END DO
!---Allocate permanent stitching arrays
stitcher%send(0)%n=m-stitcher%kle(0); stitcher%recv(0)%n=m-stitcher%kle(0)
allocate(stitcher%send(0)%v(stitcher%send(0)%n))
allocate(stitcher%recv(0)%v(stitcher%recv(0)%n))
!---Verify stitching structure
CALL oft_stitch_check(stitcher)
DEBUG_STACK_POP
END SUBROUTINE condense_stitch
!------------------------------------------------------------------------------
!> Real implementation for \ref create_matrix
!------------------------------------------------------------------------------
SUBROUTINE create_matrix_real(mat,ingraphs,row_vec,col_vec,native)
CLASS(oft_matrix), POINTER, INTENT(inout) :: mat
TYPE(oft_graph_ptr), INTENT(in) :: ingraphs(:,:)
CLASS(oft_vector), POINTER, INTENT(in) :: row_vec
CLASS(oft_vector), POINTER, INTENT(in) :: col_vec
LOGICAL, OPTIONAL, INTENT(in) :: native
LOGICAL :: force_native
INTEGER(i4) :: i,j,ni,nj,offset,soffset
INTEGER(i8) :: offsetg
DEBUG_STACK_PUSH
force_native=.FALSE.
IF(PRESENT(native))force_native=native
!---
IF(use_petsc.AND.(.NOT.force_native))THEN
#ifdef HAVE_PETSC
  ALLOCATE(oft_petsc_matrix::mat)
#else
  CALL oft_abort("Not compiled with PETSc","create_matrix_real",__FILE__)
#endif
ELSE
  ALLOCATE(oft_native_matrix::mat)
END IF
!---Commond setup tasks
ni=SIZE(ingraphs,DIM=1)
nj=SIZE(ingraphs,DIM=2)
mat%ni=ni
mat%nj=nj
!---Allocate and setup map structure
ALLOCATE(mat%i_map(mat%ni))
ALLOCATE(mat%j_map(mat%nj))
!---
offset=0
soffset=0
offsetg=0
DO i=1,mat%ni
  mat%i_map(i)%offset=offset
  mat%i_map(i)%soffset=soffset
  mat%i_map(i)%offsetg=offsetg
  mat%i_map(i)%per=row_vec%map(i)%per
  mat%i_map(i)%n=row_vec%map(i)%n
  mat%i_map(i)%ng=row_vec%map(i)%ng
  mat%i_map(i)%nslice=row_vec%map(i)%nslice
  mat%i_map(i)%slice=>row_vec%map(i)%slice
  mat%i_map(i)%lge=>row_vec%map(i)%lge
  offset=offset+mat%i_map(i)%n
  soffset=soffset+mat%i_map(i)%nslice
  offsetg=offsetg+mat%i_map(i)%ng
END DO
mat%nr=offset
mat%nrg=offsetg
mat%nrslice=soffset
!---
offset=0
soffset=0
offsetg=0
DO j=1,mat%nj
  mat%j_map(j)%offset=offset
  mat%j_map(j)%soffset=soffset
  mat%j_map(j)%offsetg=offsetg
  mat%j_map(j)%per=col_vec%map(j)%per
  mat%j_map(j)%n=col_vec%map(j)%n
  mat%j_map(j)%ng=col_vec%map(j)%ng
  mat%j_map(j)%nslice=col_vec%map(j)%nslice
  mat%j_map(j)%slice=>col_vec%map(j)%slice
  mat%j_map(j)%lge=>col_vec%map(j)%lge
  offset=offset+mat%j_map(j)%n
  soffset=soffset+mat%j_map(j)%nslice
  offsetg=offsetg+mat%j_map(j)%ng
END DO
mat%nc=offset
mat%ncg=offsetg
mat%ncslice=soffset
!---
SELECT TYPE(this=>mat)
  TYPE is(oft_native_matrix)
    CALL setup_native(this)
    IF(oft_debug_print(2))WRITE(*,'(2X,A)')'Native real CRS matrix created'
#ifdef HAVE_PETSC
  TYPE is(oft_petsc_matrix)
    CALL setup_petsc(this)
    IF(oft_debug_print(2))WRITE(*,'(2X,A)')'PETSc real matrix created'
#endif
  CLASS DEFAULT
    CALL oft_abort('Error in matrix allocation.','create_matrix_real',__FILE__)
END SELECT
DEBUG_STACK_POP
CONTAINS
!
SUBROUTINE setup_native(this)
TYPE(oft_native_matrix), INTENT(inout) :: this
TYPE(oft_graph), POINTER :: outgraph
CALL condense_graph(ingraphs,outgraph,this%map,row_vec%map,col_vec%map)
this%nnz=outgraph%nnz
this%kr=>outgraph%kr
this%lc=>outgraph%lc
DEALLOCATE(outgraph)
ALLOCATE(this%M(this%nnz))
this%M=0.d0
END SUBROUTINE setup_native
#ifdef HAVE_PETSC
SUBROUTINE setup_petsc(this)
TYPE(oft_petsc_matrix), INTENT(inout) :: this
INTEGER(i4) :: ii,jj,ierr,m,n,nnz_loc,nnz_max
INTEGER(i4), ALLOCATABLE :: nnz_dist(:)
REAL(r8), POINTER :: vals(:) => NULL()
CLASS(oft_petsc_vector), POINTER :: rowpv,colpv
CLASS(oft_vector), POINTER :: tmpv
IF(.NOT.oft_petsc_vector_cast(rowpv,row_vec))CALL oft_abort('"row_vec" is not a PETSc object.', &
  'create_matrix_real',__FILE__)
IF(.NOT.oft_petsc_vector_cast(colpv,col_vec))CALL oft_abort('"col_vec" is not a PETSc object.', &
  'create_matrix_real',__FILE__)
!---
ALLOCATE(nnz_dist(this%nrslice))
nnz_dist=0
CALL row_vec%new(tmpv)
NULLIFY(vals)
CALL tmpv%get_local(vals)
!---
DO i=1,this%ni
  nnz_loc=0
  m=0
  DO j=1,this%nj
    IF(ASSOCIATED(ingraphs(i,j)%g))THEN
      DO ii=1,this%i_map(i)%n
        vals(ii+this%i_map(i)%offset) = vals(ii+this%i_map(i)%offset) &
        + REAL(ingraphs(i,j)%g%kr(ii+1) - ingraphs(i,j)%g%kr(ii),8)
      END DO
    END IF
  END DO
END DO
!---
CALL tmpv%restore_local(vals,add=.TRUE.)
CALL tmpv%get_local(vals)
!---
DO i=1,this%ni
  DO ii=1,this%i_map(i)%nslice
    nnz_dist(ii+this%i_map(i)%soffset)=CEILING(vals(this%i_map(i)%slice(ii)+this%i_map(i)%offset),4)
  END DO
END DO
CALL tmpv%restore_local(vals)
DEALLOCATE(vals)
CALL tmpv%delete
!---
IF((this%nr==this%nrg).AND.(this%nc==this%ncg))THEN
  CALL MatCreate(PETSC_COMM_SELF,this%M,ierr)
  CALL MatSetType(this%M,MATSEQAIJ,ierr)
  CALL MatSetSizes(this%M,this%nr,this%nc,PETSC_DETERMINE,PETSC_DETERMINE,ierr)
  CALL MatSeqAIJSetPreallocation(this%M,PETSC_DEFAULT_INTEGER,INT(MIN(this%ncg,nnz_dist),4),ierr)
ELSE
  CALL MatCreate(oft_env%COMM,this%M,ierr)
  CALL MatSetType(this%M,MATMPIAIJ,ierr)
  CALL MatSetSizes(this%M,this%nrslice,this%ncslice,PETSC_DETERMINE,PETSC_DETERMINE,ierr)
  CALL MatMPIAIJSetPreallocation(this%M,PETSC_DEFAULT_INTEGER,MIN(this%ncslice,nnz_dist), &
    PETSC_DEFAULT_INTEGER,INT(MIN(this%ncg-this%ncslice,nnz_dist),4),ierr)
END IF
DEALLOCATE(nnz_dist)
!---
CALL ISLocalToGlobalMappingCreateIS(rowpv%lis,this%r_lis,ierr)
CALL ISLocalToGlobalMappingCreateIS(colpv%lis,this%c_lis,ierr)
CALL MatSetLocalToGlobalMapping(this%M,this%r_lis,this%c_lis,ierr)
!---
CALL MatGetSize(this%M,m,n,ierr)
END SUBROUTINE setup_petsc
#endif
END SUBROUTINE create_matrix_real
!------------------------------------------------------------------------------
!> Real implementation for \ref create_matrix
!------------------------------------------------------------------------------
SUBROUTINE create_matrix_comp(mat,ingraphs,row_vec,col_vec,native)
CLASS(oft_cmatrix), POINTER, INTENT(inout) :: mat
TYPE(oft_graph_ptr), INTENT(in) :: ingraphs(:,:)
CLASS(oft_cvector), POINTER, INTENT(in) :: row_vec
CLASS(oft_cvector), POINTER, INTENT(in) :: col_vec
LOGICAL, OPTIONAL, INTENT(in) :: native
LOGICAL :: force_native
INTEGER(i4) :: i,j,ni,nj,offset,soffset
INTEGER(i8) :: offsetg
DEBUG_STACK_PUSH
force_native=.FALSE.
IF(PRESENT(native))force_native=native
!---
IF(use_petsc.AND.(.NOT.force_native))THEN
#ifdef HAVE_PETSC
  ! ALLOCATE(oft_petsc_cmatrix::mat)
  CALL oft_abort("Complex matrices not supported with PETSc","create_matrix_comp",__FILE__)
#else
  CALL oft_abort("Not compiled with PETSc","create_matrix_comp",__FILE__)
#endif
ELSE
  ALLOCATE(oft_native_cmatrix::mat)
END IF
!---Commond setup tasks
ni=SIZE(ingraphs,DIM=1)
nj=SIZE(ingraphs,DIM=2)
mat%ni=ni
mat%nj=nj
!---Allocate and setup map structure
ALLOCATE(mat%i_map(mat%ni))
ALLOCATE(mat%j_map(mat%nj))
!---
offset=0
soffset=0
offsetg=0
mat%nrslice=0
mat%nr=0
mat%nrg=0
DO i=1,mat%ni
  mat%i_map(i)%offset=offset
  mat%i_map(i)%soffset=soffset
  mat%i_map(i)%offsetg=offsetg
  mat%i_map(i)%per=row_vec%map(i)%per
  mat%i_map(i)%n=row_vec%map(i)%n
  mat%i_map(i)%ng=row_vec%map(i)%ng
  mat%i_map(i)%nslice=row_vec%map(i)%nslice
  mat%i_map(i)%slice=>row_vec%map(i)%slice
  mat%i_map(i)%lge=>row_vec%map(i)%lge
  offset=offset+mat%i_map(i)%n
  soffset=soffset+mat%i_map(i)%nslice
  offsetg=offsetg+mat%i_map(i)%ng
  !---
  mat%nr=mat%nr+mat%i_map(i)%n
  mat%nrg=mat%nrg+mat%i_map(i)%ng
  mat%nrslice=mat%nrslice+mat%i_map(i)%nslice
END DO
!---
offset=0
soffset=0
offsetg=0
mat%ncslice=0
mat%nc=0
mat%ncg=0
DO j=1,mat%nj
  mat%j_map(j)%offset=offset
  mat%j_map(j)%soffset=soffset
  mat%j_map(j)%offsetg=offsetg
  mat%j_map(j)%per=col_vec%map(j)%per
  mat%j_map(j)%n=col_vec%map(j)%n
  mat%j_map(j)%ng=col_vec%map(j)%ng
  mat%j_map(j)%nslice=col_vec%map(j)%nslice
  mat%j_map(j)%slice=>col_vec%map(j)%slice
  mat%j_map(j)%lge=>col_vec%map(j)%lge
  offset=offset+mat%j_map(j)%n
  soffset=soffset+mat%j_map(j)%nslice
  offsetg=offsetg+mat%j_map(j)%ng
  !---
  mat%nc=mat%nc+mat%j_map(j)%n
  mat%ncg=mat%ncg+mat%j_map(j)%ng
  mat%ncslice=mat%ncslice+mat%j_map(j)%nslice
END DO
!---
SELECT TYPE(this=>mat)
  TYPE is(oft_native_cmatrix)
    CALL setup_native(this)
    IF(oft_debug_print(2))WRITE(*,'(2X,A)')'Native complex CRS matrix created'
  CLASS DEFAULT
    CALL oft_abort('Error in matrix allocation.','create_matrix_real',__FILE__)
END SELECT
DEBUG_STACK_POP
CONTAINS
!
SUBROUTINE setup_native(this)
TYPE(oft_native_cmatrix), INTENT(inout) :: this
CLASS(oft_native_cvector), POINTER :: rowv,colv
TYPE(oft_graph), POINTER :: outgraph
INTEGER(i4) :: ierr
INTEGER(i8) :: gsize,gtmp,gmax,gmin
!---Allocate and setup map structure
IF(.NOT.native_cvector_cast(rowv,row_vec))CALL oft_abort('"row_vec" is not a native vector.', &
  'create_matrix_real',__FILE__)
IF(.NOT.native_cvector_cast(colv,col_vec))CALL oft_abort('"col_vec" is not a native vector.', &
  'create_matrix_real',__FILE__)
!---
CALL condense_graph(ingraphs,outgraph,this%map,rowv%map,colv%map)
this%nr=outgraph%nr
this%nrg=outgraph%nrg
this%nc=outgraph%nc
this%ncg=outgraph%ncg
this%nnz=outgraph%nnz
this%kr=>outgraph%kr
this%lc=>outgraph%lc
DEALLOCATE(outgraph)
ALLOCATE(this%M(this%nnz))
this%M=0.d0
gtmp=this%nnz
IF(rowv%stitch_info%full)THEN
  gsize=gtmp
  gmax=gtmp
  gmin=gtmp
ELSE
#ifdef HAVE_MPI
  CALL MPI_ALLREDUCE(gtmp,gsize,1,OFT_MPI_I8,MPI_SUM,oft_env%COMM,ierr)
  CALL MPI_ALLREDUCE(gtmp,gmin,1,OFT_MPI_I8,MPI_MIN,oft_env%COMM,ierr)
  CALL MPI_ALLREDUCE(gtmp,gmax,1,OFT_MPI_I8,MPI_MAX,oft_env%COMM,ierr)
#else
  CALL oft_abort("Not built with MPI","create_matrix_real",__FILE__)
#endif
END IF
END SUBROUTINE setup_native
END SUBROUTINE create_matrix_comp
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE csr_remove_redundant(nr,kr,nnz,lc)
INTEGER(4), INTENT(in) :: nr
INTEGER(4), INTENT(inout) :: kr(nr+1)
INTEGER(4), INTENT(inout) :: nnz
INTEGER(4), POINTER, INTENT(inout) :: lc(:)
INTEGER(4) :: i,j,nremove,js
INTEGER(4), POINTER :: lctmp(:)
nremove=0
DO i=1,nr
  js=kr(i)+1
  IF(kr(i+1)-kr(i)>0)lc(kr(i)-nremove)=lc(kr(i))
  kr(i)=kr(i)-nremove
  DO j=js,kr(i+1)-1
    IF(lc(j)==lc(j-1))THEN
      nremove=nremove+1
    ELSE
      lc(j-nremove)=lc(j)
    END IF
  END DO
END DO
kr(nr+1)=kr(nr+1)-nremove
nnz=kr(nr+1)-1
lctmp=>lc
ALLOCATE(lc(kr(nr+1)-1))
lc=lctmp(1:kr(nr+1)-1)
DEALLOCATE(lctmp)
END SUBROUTINE csr_remove_redundant
!------------------------------------------------------------------------------
!> Combine a set of non-overlapping CRS-graphs into a graph
!------------------------------------------------------------------------------
SUBROUTINE condense_graph(ingraphs,outgraph,maps,row_map,col_map)
TYPE(oft_graph_ptr), INTENT(in) :: ingraphs(:,:) !< Array of graphs representing submatrices
TYPE(oft_graph), POINTER, INTENT(inout) :: outgraph !< Resulting graph
TYPE(oft_matrix_map), POINTER, INTENT(out) :: maps(:,:) !< Mapping from sub-graphs into full graph
TYPE(oft_map), DIMENSION(:), INTENT(in) :: row_map !< Vector representing matrix rows
TYPE(oft_map), DIMENSION(:), INTENT(in) :: col_map !< Vector representing matrix rows
!---
INTEGER(i4) :: ni,nj
INTEGER(i4) :: i,j,jj,l,k
INTEGER(i4) :: offset1,offset2
INTEGER(i4), ALLOCATABLE :: krtmp(:)
DEBUG_STACK_PUSH
!---
IF(oft_debug_print(3))WRITE(*,'(4X,A)')'Condensing matrix graphs'
ni=SIZE(ingraphs,DIM=1)
nj=SIZE(ingraphs,DIM=2)
!---Condense stitching info
ALLOCATE(maps(ni,nj))
ALLOCATE(outgraph)
!---
outgraph%nnz=0
outgraph%nr=0
outgraph%nrg=0
outgraph%nc=0
outgraph%ncg=0
DO i=1,ni
  DO j=1,nj
    IF(ASSOCIATED(ingraphs(i,j)%g))THEN
      !---
      outgraph%nnz=outgraph%nnz+ingraphs(i,j)%g%nnz
      !---
      maps(i,j)%nnz=ingraphs(i,j)%g%nnz
      ALLOCATE(maps(i,j)%lc_map(maps(i,j)%nnz))
    ELSE
      maps(i,j)%nnz=0
    END IF
  END DO
END DO
DO i=1,ni
  !---
  outgraph%nr=outgraph%nr+row_map(i)%n
  outgraph%nrg=outgraph%nrg+row_map(i)%ng
END DO
DO j=1,nj
  !---
  outgraph%nc=outgraph%nc+col_map(j)%n
  outgraph%ncg=outgraph%ncg+col_map(j)%ng
END DO
!---
ALLOCATE(outgraph%kr(outgraph%nr+1))
ALLOCATE(krtmp(outgraph%nr))
ALLOCATE(outgraph%lc(outgraph%nnz))
outgraph%kr=0
offset2=0
DO i=1,ni
  DO j=1,nj
    !---
    IF(ASSOCIATED(ingraphs(i,j)%g))THEN
      ALLOCATE(maps(i,j)%ext(2,ingraphs(i,j)%g%nr))
      maps(i,j)%ext=0
      maps(i,j)%kr=>ingraphs(i,j)%g%kr
      maps(i,j)%lc=>ingraphs(i,j)%g%lc
      DO k=1,ingraphs(i,j)%g%nr
        outgraph%kr(offset2+k)=outgraph%kr(offset2+k) + &
        ingraphs(i,j)%g%kr(k+1) - ingraphs(i,j)%g%kr(k)
      END DO
    END IF
  END DO
  offset2=offset2+row_map(i)%n
END DO
!---Generate row pointer
outgraph%kr(outgraph%nr+1)=outgraph%nnz+1
do i=outgraph%nr,1,-1 ! cumulative unique point linkage count
  outgraph%kr(i)=outgraph%kr(i+1)-outgraph%kr(i)
end do
IF((outgraph%nnz>0).AND.(outgraph%kr(1)/=1))call oft_abort('Bad element linkage count', &
  'condense_graph',__FILE__)
!---
krtmp=0
offset1=0
offset2=0
DO i=1,ni
  offset1=0
  DO j=1,nj
    !---
    IF(ASSOCIATED(ingraphs(i,j)%g))THEN
      DO k=1,ingraphs(i,j)%g%nr
        DO l=ingraphs(i,j)%g%kr(k),ingraphs(i,j)%g%kr(k+1)-1
          jj = outgraph%kr(offset2+k)+krtmp(offset2+k)+l-ingraphs(i,j)%g%kr(k)
          outgraph%lc(jj) = offset1+ingraphs(i,j)%g%lc(l)
          maps(i,j)%lc_map(l) = jj
        END DO
        maps(i,j)%ext(1,k) = outgraph%kr(offset2+k)+krtmp(offset2+k)
        maps(i,j)%ext(2,k) = outgraph%kr(offset2+k)+krtmp(offset2+k) &
          + ingraphs(i,j)%g%kr(k+1)-1-ingraphs(i,j)%g%kr(k)
        krtmp(offset2+k)=krtmp(offset2+k) + &
        ingraphs(i,j)%g%kr(k+1) - ingraphs(i,j)%g%kr(k)
      END DO
    END IF
    offset1=offset1+col_map(j)%n
  END DO
  offset2=offset2+row_map(i)%n
END DO
DEBUG_STACK_POP
END SUBROUTINE condense_graph
!------------------------------------------------------------------------------
!> Modify a CSR graph by adding dense blocks
!------------------------------------------------------------------------------
subroutine graph_add_dense_blocks(graph_in,graph_out,dense_flag,dense_nodes)
TYPE(oft_graph), INTENT(inout) :: graph_in !< Input graph to augment
TYPE(oft_graph), INTENT(inout) :: graph_out !< Output graph
INTEGER(4), INTENT(in) :: dense_flag(:) !< Index of dense blocks (0 for elements outside dense blocks)
TYPE(oft_1d_int), INTENT(in) :: dense_nodes(:) !< Indices for each dense block
INTEGER(4) :: i,j,k,nr,offset
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: ltmp
DEBUG_STACK_PUSH
graph_out%nr=graph_in%nr
graph_out%nrg=graph_in%nrg
graph_out%nc=graph_in%nc
graph_out%ncg=graph_in%ncg
graph_out%nnz=0
ALLOCATE(graph_out%kr(graph_out%nr+1))
graph_out%kr=0
ALLOCATE(ltmp(graph_in%nc))
DO i=1,graph_in%nr
  IF(dense_flag(i)/=0)THEN
    ltmp=2*graph_in%nc
    offset=dense_nodes(dense_flag(i))%n
    ltmp(1:offset)=dense_nodes(dense_flag(i))%v
    DO j=graph_in%kr(i),graph_in%kr(i+1)-1
      ltmp(j-graph_in%kr(i)+offset+1)=graph_in%lc(j)
    END DO
    CALL sort_array(ltmp,offset+graph_in%kr(i+1)-graph_in%kr(i))
    graph_out%kr(i)=1
    DO j=2,offset+graph_in%kr(i+1)-graph_in%kr(i)
      IF(ltmp(j)>ltmp(j-1))graph_out%kr(i)=graph_out%kr(i)+1
    END DO
  ELSE
    graph_out%kr(i)=graph_in%kr(i+1)-graph_in%kr(i)
  END IF
END DO
graph_out%nnz=SUM(graph_out%kr)
graph_out%kr(graph_out%nr+1)=graph_out%nnz+1
DO i=graph_out%nr,1,-1
  graph_out%kr(i)=graph_out%kr(i+1)-graph_out%kr(i)
END DO
IF(graph_out%kr(1)/=1)CALL oft_abort('Bad new graph setup','graph_add_dense_blocks', &
  __FILE__)
ALLOCATE(graph_out%lc(graph_out%nnz))
DO i=1,graph_in%nr
  IF(dense_flag(i)/=0)THEN
    ltmp=2*graph_in%nc
    offset=dense_nodes(dense_flag(i))%n
    ltmp(1:offset)=dense_nodes(dense_flag(i))%v
    DO j=graph_in%kr(i),graph_in%kr(i+1)-1
      ltmp(j-graph_in%kr(i)+offset+1)=graph_in%lc(j)
    END DO
    CALL sort_array(ltmp,offset+graph_in%kr(i+1)-graph_in%kr(i))
    nr=0
    graph_out%lc(graph_out%kr(i))=ltmp(1)
    DO j=2,offset+graph_in%kr(i+1)-graph_in%kr(i)
      IF(ltmp(j)>ltmp(j-1))THEN
      nr=nr+1
      graph_out%lc(graph_out%kr(i)+nr)=ltmp(j)
      END IF
    END DO
  ELSE
    DO j=graph_in%kr(i),graph_in%kr(i+1)-1
      graph_out%lc(graph_out%kr(i)+j-graph_in%kr(i))=graph_in%lc(j)
    END DO
  END IF
END DO
!---
DEALLOCATE(ltmp)
DEBUG_STACK_POP
end subroutine graph_add_dense_blocks
!------------------------------------------------------------------------------
!> Modify a CSR graph by adding dense columns at specified indices
!------------------------------------------------------------------------------
subroutine graph_add_full_col(graph_in,graph_out,nadd,nodes_add)
TYPE(oft_graph), INTENT(inout) :: graph_in !< Input graph to augment
TYPE(oft_graph), INTENT(inout) :: graph_out !< Output graph
INTEGER(4), INTENT(in) :: nadd !< Number of columns to add
INTEGER(4), INTENT(in) :: nodes_add(nadd) !< Indices of columns to add
INTEGER(4) :: i,j,k,nr,iadd
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: ltmp
DEBUG_STACK_PUSH
graph_out%nr=graph_in%nr
graph_out%nrg=graph_in%nrg
graph_out%nc=graph_in%nc
graph_out%ncg=graph_in%ncg
graph_out%nnz=0
ALLOCATE(graph_out%kr(graph_in%nr+1))
graph_out%kr=0
ALLOCATE(ltmp(graph_in%nc))
DO i=1,graph_in%nr
  graph_out%kr(i)=graph_in%kr(i+1)-graph_in%kr(i)
  DO iadd=1,nadd
    j=search_array(nodes_add(iadd),graph_in%lc(graph_in%kr(i):graph_in%kr(i+1)-1), &
      graph_in%kr(i+1)-graph_in%kr(i))
    IF(j==0)graph_out%kr(i)=graph_out%kr(i)+1
  END DO
END DO
graph_out%nnz=SUM(graph_out%kr)
graph_out%kr(graph_in%nr+1)=graph_out%nnz+1
DO i=graph_in%nr,1,-1
  graph_out%kr(i)=graph_out%kr(i+1)-graph_out%kr(i)
END DO
IF(graph_out%kr(1)/=1)CALL oft_abort('Bad new graph setup','graph_add_full_col', &
  __FILE__)
ALLOCATE(graph_out%lc(graph_out%nnz))
DO i=1,graph_in%nr
  DO j=graph_in%kr(i),graph_in%kr(i+1)-1
    graph_out%lc(graph_out%kr(i)+j-graph_in%kr(i))=graph_in%lc(j)
  END DO
  k=0
  DO iadd=1,nadd
    j=search_array(nodes_add(iadd),graph_in%lc(graph_in%kr(i):graph_in%kr(i+1)-1), &
      graph_in%kr(i+1)-graph_in%kr(i))
    IF(j==0)THEN
      k=k+1
      graph_out%lc(graph_out%kr(i+1)-k)=nodes_add(iadd)
    END IF
  END DO
  IF((graph_out%kr(i+1)-graph_out%kr(i))>(graph_in%kr(i+1)-graph_in%kr(i)))THEN
    CALL sort_array(graph_out%lc(graph_out%kr(i):graph_out%kr(i+1)-1), &
      graph_out%kr(i+1)-graph_out%kr(i))
  END IF
END DO
!---
DEALLOCATE(ltmp)
DEBUG_STACK_POP
end subroutine graph_add_full_col
!------------------------------------------------------------------------------
!> Real implementation for \ref combine_matrices
!------------------------------------------------------------------------------
SUBROUTINE combine_matrices_real(mats,nr,nc,mat)
TYPE(oft_matrix_ptr), INTENT(in) :: mats(:,:)
INTEGER(i4), INTENT(in) :: nr
INTEGER(i4), INTENT(in) :: nc
CLASS(oft_matrix), POINTER, INTENT(inout) :: mat
DEBUG_STACK_PUSH
SELECT TYPE(this=>mat)
  TYPE is(oft_native_matrix)
    CALL combine_native(this)
    IF(oft_debug_print(2))WRITE(*,'(2X,A)')'Native real CRS matrices combined'
#ifdef HAVE_PETSC
  TYPE IS(oft_petsc_matrix)
    CALL combine_petsc(this)
    IF(oft_debug_print(2))WRITE(*,'(2X,A)')'PETSc real matrices combined'
#endif
  CLASS DEFAULT
    CALL oft_abort('Error in matrix allocation.','combine_matrices_real',__FILE__)
END SELECT
DEBUG_STACK_POP
CONTAINS
!
SUBROUTINE combine_native(this)
TYPE(oft_native_matrix), INTENT(inout) :: this
INTEGER(i4) :: i,j,ii,jj,jp,jn
this%M=0.d0
DO i=1,nr
  DO j=1,nc
    IF(.NOT.ASSOCIATED(mats(i,j)%m))CYCLE
    SELECT TYPE(cmat=>mats(i,j)%m)
      TYPE IS(oft_native_matrix)
        DO ii=1,this%i_map(i)%n
          jp=this%map(i,j)%kr(ii)
          jn=this%map(i,j)%kr(ii+1)-1
          DO jj=jp,jn
            this%M(this%map(i,j)%lc_map(jj)) = cmat%M(jj)
          END DO
        END DO
    CLASS DEFAULT
      CALL oft_abort("Source matrix of different type","combine_matrices_real",__FILE__)
    END SELECT
  END DO
END DO
END SUBROUTINE combine_native
!
#ifdef HAVE_PETSC
SUBROUTINE combine_petsc(this)
TYPE(oft_petsc_matrix), INTENT(inout) :: this
TYPE(tmat), ALLOCATABLE :: pmats(:,:)
INTEGER(i4) :: i,j,ii,jj,ierr,m,n,offset,proc
INTEGER(i4) :: mcomm,roffset,coffset,rgoffset,cgoffset,ncols,maxcols,rows(1),ncslice
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: rstart,cstart,cols,coltmp
INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: goffsets,offsets,tmpoff
REAL(r8), POINTER :: vals(:),varray(:,:)
INTEGER(i4), POINTER :: int_tmp(:)
REAL(r8), POINTER :: real_tmp(:)
mcomm=oft_env%COMM
IF(this%nr==this%nrg)mcomm=PETSC_COMM_SELF
CALL MatGetOwnershipRange(this%m,m,n,ierr)
rgoffset=m
!---
ALLOCATE(pmats(nr,nc))
ALLOCATE(rstart(nr),cstart(nc))
rstart=-1
cstart=-1
!---
DO i=1,nr
  DO j=1,nc
    pmats(i,j)=PETSC_NULL_MAT
    IF(.NOT.ASSOCIATED(mats(i,j)%m))CYCLE
    SELECT TYPE(pmat=>mats(i,j)%m)
      TYPE IS(oft_petsc_matrix)
        IF(rstart(i)<0)rstart(i)=pmat%i_map(1)%nslice
        IF(cstart(j)<0)cstart(j)=pmat%j_map(1)%nslice
      CLASS DEFAULT
        CALL oft_abort("Incorrect matrix type","combine_matrices_real::combine_petsc",__FILE__)
    END SELECT
  END DO
END DO
!---
ALLOCATE(offsets(nc,oft_env%nprocs+1),tmpoff(nc,oft_env%nprocs+1))
tmpoff=0
DO i=1,nr
  DO j=1,nc
    IF(ASSOCIATED(mats(i,j)%m))THEN
      SELECT TYPE(pmat=>mats(i,j)%m)
        TYPE IS(oft_petsc_matrix)
          CALL MatGetOwnershipRangeColumn(pmat%m,m,n,ierr)
          tmpoff(j,oft_env%rank+1)=m
        CLASS DEFAULT
          CALL oft_abort("Incorrect matrix type","combine_matrices_real::combine_petsc",__FILE__)
      END SELECT
    END IF
  END DO
END DO
IF(this%nr/=this%nrg)THEN
  CALL MPI_ALLREDUCE(tmpoff,offsets,nc*(oft_env%nprocs+1),OFT_MPI_I4,MPI_SUM,oft_env%COMM,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','combine_matrices_real',__FILE__)
  offsets(:,oft_env%nprocs+1)=offsets(:,oft_env%nprocs)+INT(this%ncg,4)
  DEALLOCATE(tmpoff)
  ALLOCATE(goffsets(nc,oft_env%nprocs),tmpoff(nc,oft_env%nprocs))
  tmpoff=0
  CALL MatGetOwnershipRangeColumn(this%m,m,n,ierr)
  tmpoff(1,oft_env%rank+1)=m
  DO j=2,nc
    tmpoff(j,oft_env%rank+1)=tmpoff(j-1,oft_env%rank+1)+this%j_map(j-1)%nslice
  END DO
  CALL MPI_ALLREDUCE(tmpoff,goffsets,nc*(oft_env%nprocs),OFT_MPI_I4,MPI_SUM,oft_env%COMM,ierr)
  IF(ierr/=0)CALL oft_abort('Error in MPI_ALLREDUCE','combine_matrices_real',__FILE__)
  DEALLOCATE(tmpoff)
ELSE
  offsets=-1
  offsets(:,oft_env%rank+1)=tmpoff(:,oft_env%rank+1)
  offsets(:,oft_env%rank+2)=INT(this%ncg,4)
  ALLOCATE(goffsets(nc,oft_env%nprocs))
  goffsets(1,:)=0
  DO j=2,nc
    goffsets(j,:)=goffsets(j-1,:)+this%j_map(j-1)%nslice
  END DO
  DEALLOCATE(tmpoff)
END IF
!---
DO i=1,nr
  CALL MatGetOwnershipRangeColumn(this%m,m,n,ierr)
  cgoffset=m
  ncslice=0
  DO j=1,nc
    IF(ASSOCIATED(mats(i,j)%m))THEN
    SELECT TYPE(pmat=>mats(i,j)%m)
      TYPE IS(oft_petsc_matrix)
        CALL MatGetOwnershipRange(pmat%m,m,n,ierr)
        roffset=m
        CALL MatGetOwnershipRangeColumn(pmat%m,m,n,ierr)
        coffset=m
        pmats(i,j)=pmat%m
        maxcols=0
        NULLIFY(int_tmp,real_tmp)
        DO ii=1,pmat%i_map(1)%nslice
          CALL MatGetRow(pmat%m,roffset+ii-1,ncols,int_tmp,real_tmp,ierr)
          maxcols=MAX(ncols,maxcols)
          CALL MatRestoreRow(pmat%m,roffset+ii-1,ncols,int_tmp,real_tmp,ierr)
        END DO
        ALLOCATE(varray(maxcols,1),cols(maxcols),coltmp(maxcols))
        cols=0
        vals=>varray(:,1)
        m=0
        DO ii=1,pmat%i_map(1)%nslice
          rows=ii-1
          CALL MatGetRow(pmat%m,roffset+ii-1,ncols,cols,vals,ierr)
          proc=1
          DO jj=1,ncols
            DO WHILE(cols(jj)>=offsets(j,proc+1))
              proc=proc+1
              IF(proc>oft_env%nprocs)WRITE(*,*)proc,cols(jj),offsets(j,proc)
            END DO
            coltmp(jj)=cols(jj)-offsets(j,proc)+goffsets(j,proc)
          END DO
          CALL MatSetValues(this%M,1,rows+rgoffset,ncols, &
            coltmp,varray,ADD_VALUES,ierr)
          CALL MatRestoreRow(pmat%m,roffset+ii-1,ncols,cols,vals,ierr)
        END DO
        DEALLOCATE(varray,cols,coltmp)
    END SELECT
    END IF
    cgoffset=cgoffset+cstart(j)
  END DO
  rgoffset=rgoffset+rstart(i)
END DO
DEALLOCATE(rstart,cstart,goffsets,offsets)
CALL this%assemble()
END SUBROUTINE combine_petsc
#endif
END SUBROUTINE combine_matrices_real
!------------------------------------------------------------------------------
!> Real implementation for \ref combine_matrices
!------------------------------------------------------------------------------
SUBROUTINE combine_matrices_comp(mats,nr,nc,mat)
TYPE(oft_cmatrix_ptr), INTENT(in) :: mats(:,:)
INTEGER(i4), INTENT(in) :: nr
INTEGER(i4), INTENT(in) :: nc
CLASS(oft_cmatrix), POINTER, INTENT(inout) :: mat
DEBUG_STACK_PUSH
SELECT TYPE(this=>mat)
  TYPE is(oft_native_cmatrix)
    CALL combine_native(this)
    IF(oft_debug_print(2))WRITE(*,'(2X,A)')'Native complex CRS matrices combined'
  CLASS DEFAULT
    CALL oft_abort('Error in matrix allocation.','combine_matrices_comp',__FILE__)
END SELECT
DEBUG_STACK_POP
CONTAINS
!
SUBROUTINE combine_native(this)
TYPE(oft_native_cmatrix), INTENT(inout) :: this
INTEGER(i4) :: i,j,ii,jj,jp,jn
this%M=(0.d0,0.d0)
DO i=1,nr
  DO j=1,nc
    IF(.NOT.ASSOCIATED(mats(i,j)%m))CYCLE
    SELECT TYPE(cmat=>mats(i,j)%m)
      TYPE IS(oft_native_cmatrix)
        DO ii=1,this%i_map(i)%n
          jp=this%map(i,j)%kr(ii)
          jn=this%map(i,j)%kr(ii+1)-1
          DO jj=jp,jn
            this%M(this%map(i,j)%lc_map(jj)) = cmat%M(jj)
          END DO
        END DO
    END SELECT
  END DO
END DO
END SUBROUTINE combine_native
END SUBROUTINE combine_matrices_comp
!------------------------------------------------------------------------------
!> Create an identity graph for a given vector
!------------------------------------------------------------------------------
SUBROUTINE create_identity_graph(outgraph,vec)
TYPE(oft_graph), POINTER, INTENT(inout) :: outgraph !< Resulting graph
CLASS(oft_vector), POINTER, INTENT(in) :: vec !< Vector representing matrix rows/columns
INTEGER(i4) :: i
DEBUG_STACK_PUSH
!---Setup graph
ALLOCATE(outgraph)
outgraph%nnz=vec%n
outgraph%nr=vec%n
outgraph%nrg=vec%ng
outgraph%nc=vec%n
outgraph%ncg=vec%ng
!---Create indexing
ALLOCATE(outgraph%kr(outgraph%nr+1))
ALLOCATE(outgraph%lc(outgraph%nnz))
outgraph%kr=[(i,i=1,outgraph%nr+1)]
outgraph%lc=[(i,i=1,outgraph%nr)]
DEBUG_STACK_POP
END SUBROUTINE create_identity_graph
END MODULE oft_la_utils
