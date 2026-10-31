!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file thin_wall_hodlr.F90
!
!> Hierarchical Off-Diagonal Low Rank matrix approximation functionality for ThinCurr
!!  - SVD compression
!!  - Adaptive Cross Approximation
!!
!! @authors Chris Hansen
!! @date Feb 2024
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
MODULE thin_wall_hodlr
USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_float, c_sizeof
USE thin_wall
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
type :: oft_tw_block
  INTEGER(4) :: np = 0
  INTEGER(4) :: nelems = 0
  INTEGER(4) :: ncells = 0
  INTEGER(4) :: parent = -1
  REAL(8) :: center(3) = 0.d0
  REAL(8) :: extent = 0.d0
  INTEGER(4), POINTER :: ielem(:) => NULL()
  INTEGER(4), POINTER :: inv_map(:) => NULL()
  INTEGER(4), POINTER :: ipts(:) => NULL()
  INTEGER(4), POINTER :: icell(:) => NULL()
end type oft_tw_block
INTEGER(4), PARAMETER :: tw_hole_aca_rmin = 4 !< Minimum ACA rank for groups of cells far from a hole
REAL(8), PARAMETER :: tw_hole_aca_sep = 1.d0 !< Separation ratio (distance/(sum of extents)) for ACA of hole columns
REAL(8), PARAMETER :: tw_hole_aca_tol = 1.d-6 !< Relative ACA tolerance for hole columns
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
type oft_tw_level
  INTEGER(4) :: nblocks = 0
  INTEGER(4), POINTER, DIMENSION(:,:) :: mat_mask => NULL()
  TYPE(oft_tw_block), POINTER, DIMENSION(:) :: blocks => NULL()
end type oft_tw_level
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
type, extends(oft_noop_matrix) :: oft_tw_hodlr_op
  INTEGER(4) :: nlevels = 0 !< Number of levels in block heirarchy
  INTEGER(4) :: nblocks = 0 !< Number of blocks on the lowest level
  INTEGER(4) :: ndense = 0 !< Number of diagonal interactions
  INTEGER(4) :: nsparse = 0 !< Number of compressed off-diagonal interactions
  INTEGER(4) :: leaf_target = 1500 !< Target size for leaves on lowest level
  INTEGER(4) :: aca_min_its = 20 !< Minimum number of ACA+ iterations
  INTEGER(4) :: min_rank = 10 !< Minimum rank of compressed off-diagonal matrices (>= 1)
  REAL(r8) :: times(4) = 0 !< Timing information for diagonal and sparse matrix products
  REAL(8) :: aca_sep_thresh(2) = [1.1d0,1.1d0] !< Separation threshold for blocks ACA+ vs SVD [lowest level, all other levels]
  REAL(8) :: L_svd_tol = -1.d0 !< SVD tolerance for inductance matrix
  REAL(8) :: L_aca_tol = -1.d0 !< ACA+ tolerance for inductance matrix
  REAL(8) :: B_svd_tol = -1.d0 !< SVD tolerance for B-field operators
  REAL(8) :: B_aca_tol = -1.d0 !< ACA+ tolerance for B-field operators
  LOGICAL :: pre_single = .FALSE. !< Compute, store, and apply block-Jacobi preconditioners in single precision?
  LOGICAL :: pre_coarse = .FALSE. !< Add a global coarse correction to block-Jacobi preconditioners (two-level preconditioner)?
  INTEGER(4) :: pre_neigs = 0 !< Number of slowest L/R eigenmodes added to the coarse space (time-dependent runs only; requires pre_coarse and ARPACK)
  INTEGER(4) :: pre_evecs_hash = 0 !< Hash of the resistance matrix used for `pre_evecs`
  REAL(r8), POINTER, DIMENSION(:,:) :: pre_evecs => NULL() !< Slowest L/R eigenmodes for the coarse space, kept between runs [nelems,:]
  INTEGER(4), POINTER, DIMENSION(:,:) :: dense_blocks => NULL() !< Indices of diagonal interactions
  INTEGER(4), POINTER, DIMENSION(:,:) :: sparse_blocks => NULL() !< Indices of compressed off-diagonal interactions
  INTEGER(4), POINTER, DIMENSION(:,:) :: L_thread_ptr => NULL()
  INTEGER(4), POINTER, DIMENSION(:,:) :: B_thread_ptr => NULL()
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:) :: aca_U_mats => NULL() !< U matrices for compressed off-diagonal interactions (L)
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:) :: aca_V_mats => NULL() !< V matrices for compressed off-diagonal interactions (L)
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:) :: aca_dense => NULL() !< Fallback dense matrices for compressed off-diagonal interactions (L)
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:) :: dense_mats => NULL() !< Dense matrices for diagonal interactions (L)
  TYPE(oft_native_dense_matrix) :: hole_Vcoil_mat !< Dense coupling matrix to holes and Vcoils (L)
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:,:) :: aca_BU_mats => NULL() !< U matrices for compressed off-diagonal interactions (B-field)
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:,:) :: aca_BV_mats => NULL() !< V matrices for compressed off-diagonal interactions (B-field)
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:,:) :: aca_B_dense => NULL() !< Fallback dense matrices for compressed off-diagonal interactions (B-field)
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:,:) :: dense_B_mats => NULL() !< Dense matrices for diagonal interactions (B-field)
  REAL(8), POINTER, DIMENSION(:,:,:) :: hole_Vcoil_Bmat => NULL() !< Dense coupling matrix to holes and Vcoils (B-field)
  REAL(8), POINTER, DIMENSION(:,:,:) :: Icoil_Bmat => NULL() !< Dense coupling matrix to Icoils (B-field) [np,n_icoils,3]
  TYPE(oft_tw_block), POINTER, DIMENSION(:) :: blocks => NULL() !< REMOVE
  TYPE(oft_tw_level), POINTER, DIMENSION(:) :: levels => NULL() !< Block heirarchy
  type(tw_type), pointer :: tw_obj => NULL()
contains
  !> Setup HODLR by performing partitioning and tagging block-block interactions
  procedure :: setup => tw_hodlr_setup
  !> Compute compressed inductance operator
  procedure :: compute_L => tw_hodlr_Lcompute
  !> Compute compressed B-field operators
  procedure :: compute_B => tw_hodlr_Bcompute
  !> Apply the matrix (L)
  procedure :: apply_real => tw_hodlr_Lapply
  !> Apply B-field operator
  procedure :: apply_bop => tw_hodlr_Bapply
  !> Assemble matrix (L)
  procedure :: assemble => tw_hodlr_Lassemble
end type oft_tw_hodlr_op
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
TYPE :: cmat_container
  COMPLEX(c8), POINTER, DIMENSION(:,:) :: M => NULL() !< Double precision matrix
  COMPLEX(c4), POINTER, DIMENSION(:,:) :: M4 => NULL() !< Single precision matrix
END TYPE
!------------------------------------------------------------------------------
!> Complex block-Jacobi preconditioner for ThinCurr HODLR matrices
!------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_csolver) :: oft_tw_hodlr_bjpre
  LOGICAL :: refactor = .TRUE. !< Refactor solution on next application
  INTEGER(4) :: max_block_size = 0
  COMPLEX(c8) :: alpha = (1.d0,0.d0)
  COMPLEX(c8) :: beta =  (1.d0,0.d0)
  TYPE(oft_tw_hodlr_op), POINTER :: mf_obj => NULL()
  TYPE(oft_native_matrix), POINTER :: Rmat => NULL()
  TYPE(cmat_container), POINTER, DIMENSION(:) :: inverse_mats => NULL()
CONTAINS
  !> Solve system
  PROCEDURE :: apply => bjprecond_apply
  !> Clean-up internal storage
  PROCEDURE :: delete => bjprecond_delete
END TYPE oft_tw_hodlr_bjpre
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
TYPE :: rmat_container
  REAL(r8), POINTER, DIMENSION(:,:) :: M => NULL() !< Double precision matrix
  REAL(r4), POINTER, DIMENSION(:,:) :: M4 => NULL() !< Single precision matrix
END TYPE
!------------------------------------------------------------------------------
!> Real block-Jacobi preconditioner for ThinCurr HODLR matrices
!------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_solver) :: oft_tw_hodlr_rbjpre
  LOGICAL :: refactor = .TRUE. !< Refactor solution on next application
  INTEGER(4) :: coarse_order = 1 !< Polynomial order of the coarse space on each block (1 or 2)
  INTEGER(4) :: max_block_size = 0
  INTEGER(4) :: nz = 0 !< Number of coarse vectors
  INTEGER(4), POINTER, DIMENSION(:) :: z_block => NULL() !< Dense block of each coarse vector (0: unit vector, -k: global vector k in z_gvals)
  INTEGER(4), POINTER, DIMENSION(:) :: z_unit => NULL() !< Element of each unit coarse vector (0 for block and global vectors)
  REAL(r8), POINTER, DIMENSION(:,:) :: z_vals => NULL() !< Values of each block coarse vector on its block
  REAL(r8), POINTER, DIMENSION(:,:) :: z_glob => NULL() !< Global vectors to add to the coarse space [nelems,:] (optional, not owned)
  REAL(r8), POINTER, DIMENSION(:,:) :: z_gvals => NULL() !< Values of each global coarse vector [nelems_exp,:]
  INTEGER(4) :: nzg = 0 !< Number of global coarse vectors (the last `nzg` of the `nz`, columns of `z_gvals`)
  REAL(r8), POINTER, DIMENSION(:,:) :: zlz => NULL() !< Z^T L Z
  REAL(r8), POINTER, DIMENSION(:,:) :: zrz => NULL() !< Z^T R Z
  REAL(r8), POINTER, DIMENSION(:,:) :: einv => NULL() !< (Z^T (alpha*L + beta*R) Z)^-1
  REAL(r8) :: alpha = 1.d0
  REAL(r8) :: beta = 1.d0
  REAL(r8) :: times = 0.d0
  TYPE(oft_tw_hodlr_op), POINTER :: mf_obj => NULL()
  TYPE(oft_native_matrix), POINTER :: Rmat => NULL()
  TYPE(rmat_container), POINTER, DIMENSION(:) :: inverse_mats => NULL()
CONTAINS
  !> Solve system
  PROCEDURE :: apply => rbjprecond_apply
  !> Clean-up internal storage
  PROCEDURE :: delete => rbjprecond_delete
END TYPE oft_tw_hodlr_rbjpre
! Helper struct for HODLR matrices settings hash
TYPE, BIND(c) :: HODLR_hash
  INTEGER(c_int) :: nsparse
  REAL(c_float) :: svd_tol
  REAL(c_float) :: aca_tol
END TYPE HODLR_hash
!------------------------------------------------------------------------------
!> Inductance block between two blocks of elements (rows: elements of `iblock`, columns: elements of `jblock`)
!------------------------------------------------------------------------------
TYPE, EXTENDS(tw_aca_source) :: tw_aca_Lblock
  CLASS(oft_tw_hodlr_op), POINTER :: op => NULL() !< HODLR operator
  INTEGER(4) :: ilevel = 0 !< Level of the blocks
  INTEGER(4) :: iblock = 0 !< Row block
  INTEGER(4) :: jblock = 0 !< Column block
  TYPE(oft_tw_block) :: one !< Single element block
  REAL(8), ALLOCATABLE :: buf(:,:) !< Work array
CONTAINS
  PROCEDURE :: row => tw_aca_Lblock_row
  PROCEDURE :: col => tw_aca_Lblock_col
END TYPE tw_aca_Lblock
!------------------------------------------------------------------------------
!> B-field block stacked over the three directions (rows: elements of `iblock`, columns: points of `jblock` times
!! direction). A column (one point and direction) needs all three directions at that point, so they are kept for
!! the other directions of the same point.
!------------------------------------------------------------------------------
TYPE, EXTENDS(tw_aca_source) :: tw_aca_B3block
  CLASS(oft_tw_hodlr_op), POINTER :: op => NULL() !< HODLR operator
  INTEGER(4) :: ilevel = 0 !< Level of the blocks
  INTEGER(4) :: iblock = 0 !< Element block
  INTEGER(4) :: jblock = 0 !< Point block
  INTEGER(4) :: nslot = 0 !< Number of points with stored columns
  INTEGER(4), ALLOCATABLE :: slot(:) !< Storage slot of each block point (0: not computed)
  REAL(8), ALLOCATABLE :: cols(:,:,:) !< Stored columns of computed points [nelems,3,:]
  TYPE(oft_tw_block) :: one_elem !< Single element block
  TYPE(oft_tw_block) :: one_pt !< Single point block
CONTAINS
  PROCEDURE :: row => tw_aca_B3block_row
  PROCEDURE :: col => tw_aca_B3block_col
END TYPE tw_aca_B3block
!------------------------------------------------------------------------------
!> Cell-pair kernel between a group of cells (rows) and the cells around a hole (columns)
!------------------------------------------------------------------------------
TYPE, EXTENDS(tw_aca_source) :: tw_aca_cellblock
  CLASS(oft_bmesh), POINTER :: mesh => NULL() !< Mesh
  INTEGER(4), POINTER :: rcells(:) => NULL() !< Row cells
  INTEGER(4), POINTER :: ccells(:) => NULL() !< Column cells
CONTAINS
  PROCEDURE :: row => tw_aca_cellblock_row
  PROCEDURE :: col => tw_aca_cellblock_col
END TYPE tw_aca_cellblock
CONTAINS
!------------------------------------------------------------------------------
!> Hole and V-coil columns of the inductance matrix, with ACA for groups of cells far from each hole
!!
!! A hole column is \f$ K C_h \f$: the cell-pair kernel K between all cells and the cells around hole h, times the
!! current vector of each of those cells. Cells are split into disjoint groups by the leaf block of their first
!! vertex in a leaf; cells with no vertex in any leaf form one more group, which is always summed directly. Groups
!! well separated from a hole's cells use an ACA approximation of their block of K; the others, and
!! groups where ACA does not converge at low rank, are summed directly. V-coil columns are copied from the coil
!! couplings.
!------------------------------------------------------------------------------
SUBROUTINE tw_hodlr_LmatHole(self,Lmat)
CLASS(oft_tw_hodlr_op), INTENT(inout) :: self !< HODLR operator
REAL(8), CONTIGUOUS, POINTER, INTENT(inout) :: Lmat(:,:) !< Hole and V-coil columns [nelems_exp,nholes+n_vcoils]
INTEGER(4) :: i,ii,j,jj,h,ib,ig,nleaf,ngroup,ncell,nstrip,jk
INTEGER(4), ALLOCATABLE :: cell_owner(:),vert_leaf(:),kgroup(:),lgroup(:),kstrip(:),lstrip(:),cnt(:)
REAL(8), ALLOCATABLE :: cstrip(:,:),gcen(:,:),gext(:),scen(:,:),sext(:),V(:,:)
TYPE(tw_type), POINTER :: tw
CLASS(oft_bmesh), POINTER :: mesh
tw=>self%tw_obj
mesh=>tw%mesh
Lmat=0.d0
!---Disjoint cell groups from the leaf blocks
nleaf=self%levels(self%nlevels)%nblocks
!---Leaf block of each vertex (the first leaf containing it), then each cell goes to the leaf of its first vertex
!   in a leaf; cells with no active vertex form their own group, always summed directly
ALLOCATE(vert_leaf(mesh%np),cell_owner(mesh%nc))
vert_leaf=0
DO ib=nleaf,1,-1
  DO j=1,self%levels(self%nlevels)%blocks(ib)%nelems
    vert_leaf(tw%lpmap_inv(self%levels(self%nlevels)%blocks(ib)%ielem(j)))=ib
  END DO
END DO
cell_owner=nleaf+1
DO i=1,mesh%nc
  DO ii=1,3
    IF(vert_leaf(mesh%lc(ii,i))>0)THEN
      cell_owner(i)=vert_leaf(mesh%lc(ii,i))
      EXIT
    END IF
  END DO
END DO
DEALLOCATE(vert_leaf)
ngroup=nleaf+1
ALLOCATE(kgroup(ngroup+1),lgroup(mesh%nc),cnt(ngroup))
cnt=0
DO i=1,mesh%nc
  cnt(cell_owner(i))=cnt(cell_owner(i))+1
END DO
kgroup(1)=1
DO ig=1,ngroup
  kgroup(ig+1)=kgroup(ig)+cnt(ig)
END DO
cnt=0
DO i=1,mesh%nc
  ig=cell_owner(i)
  lgroup(kgroup(ig)+cnt(ig))=i
  cnt(ig)=cnt(ig)+1
END DO
ALLOCATE(gcen(3,ngroup),gext(ngroup))
DO ig=1,ngroup
  CALL cell_cluster(lgroup(kgroup(ig):kgroup(ig+1)-1),gcen(:,ig),gext(ig))
END DO
!---Cells around each hole with their current vectors (sign * basis vector), combined per cell
DEALLOCATE(cnt)
ALLOCATE(cnt(tw%nholes))
cnt=0
DO i=1,mesh%nc
  DO ii=tw%kfh(i),tw%kfh(i+1)-1
    cnt(ABS(tw%lfh(1,ii)))=cnt(ABS(tw%lfh(1,ii)))+1
  END DO
END DO
ALLOCATE(kstrip(tw%nholes+1))
kstrip(1)=1
DO h=1,tw%nholes
  kstrip(h+1)=kstrip(h)+cnt(h)
END DO
ALLOCATE(lstrip(kstrip(tw%nholes+1)-1),cstrip(3,kstrip(tw%nholes+1)-1))
cnt=0
DO i=1,mesh%nc
  DO ii=tw%kfh(i),tw%kfh(i+1)-1
    h=ABS(tw%lfh(1,ii))
    IF(cnt(h)>0)THEN
      IF(lstrip(kstrip(h)+cnt(h)-1)==i)THEN
        cstrip(:,kstrip(h)+cnt(h)-1)=cstrip(:,kstrip(h)+cnt(h)-1)+SIGN(1,tw%lfh(1,ii))*tw%qbasis(:,tw%lfh(2,ii),i)
        CYCLE
      END IF
    END IF
    lstrip(kstrip(h)+cnt(h))=i
    cstrip(:,kstrip(h)+cnt(h))=SIGN(1,tw%lfh(1,ii))*tw%qbasis(:,tw%lfh(2,ii),i)
    cnt(h)=cnt(h)+1
  END DO
END DO
ALLOCATE(scen(3,tw%nholes),sext(tw%nholes))
DO h=1,tw%nholes
  CALL cell_cluster(lstrip(kstrip(h):kstrip(h)+cnt(h)-1),scen(:,h),sext(h))
END DO
!---Columns
ALLOCATE(V(3,mesh%nc))
DO h=1,tw%nholes
  nstrip=cnt(h)
  IF(nstrip==0)CYCLE
  !$omp parallel do schedule(dynamic,1) private(ncell)
  DO ig=1,ngroup
    ncell=kgroup(ig+1)-kgroup(ig)
    IF(ncell==0)CYCLE
    IF((ig<=nleaf).AND.(MIN(ncell,nstrip)/3>=tw_hole_aca_rmin).AND. &
      (magnitude(gcen(:,ig)-scen(:,h))/(gext(ig)+sext(h))>tw_hole_aca_sep))THEN
      IF(group_aca(lgroup(kgroup(ig):kgroup(ig+1)-1),lstrip(kstrip(h):kstrip(h)+nstrip-1), &
        cstrip(:,kstrip(h):kstrip(h)+nstrip-1),V))CYCLE
    END IF
    CALL group_direct(lgroup(kgroup(ig):kgroup(ig+1)-1),lstrip(kstrip(h):kstrip(h)+nstrip-1), &
      cstrip(:,kstrip(h):kstrip(h)+nstrip-1),V)
  END DO
  !---Scatter to element and hole rows
  DO j=1,mesh%nc
    DO jj=1,3
      jk=tw%pmap(mesh%lc(jj,j))
      IF(jk<=0)CYCLE
      Lmat(jk,h)=Lmat(jk,h)+DOT_PRODUCT(tw%qbasis(:,jj,j),V(:,j))
    END DO
    DO jj=tw%kfh(j),tw%kfh(j+1)-1
      jk=ABS(tw%lfh(1,jj))+tw%np_active
      Lmat(jk,h)=Lmat(jk,h)+SIGN(1,tw%lfh(1,jj))*DOT_PRODUCT(tw%qbasis(:,tw%lfh(2,jj),j),V(:,j))
    END DO
  END DO
END DO
!---Passive coils
DO i=1,tw%np_active+tw%nholes
  DO j=1,tw%n_vcoils
    Lmat(i,tw%nholes+j)=tw%Ael2coil(i,j)
  END DO
END DO
DO i=1,tw%nholes
  DO j=1,tw%n_vcoils
    Lmat(tw%np_active+tw%nholes+j,i)=tw%Ael2coil(tw%np_active+i,j)
  END DO
END DO
DO i=1,tw%n_vcoils
  DO j=1,tw%n_vcoils
    Lmat(tw%np_active+tw%nholes+j,tw%nholes+i)=tw%Acoil2coil(i,j)
  END DO
END DO
DEALLOCATE(cell_owner,kgroup,lgroup,cnt,gcen,gext,kstrip,lstrip,cstrip,scen,sext,V)
Lmat=Lmat/(4.d0*pi)
CONTAINS
!---Center and extent of a set of cells
SUBROUTINE cell_cluster(cells,c,e)
INTEGER(4), INTENT(in) :: cells(:)
REAL(8), INTENT(out) :: c(3),e
INTEGER(4) :: m,n
c=0.d0; e=0.d0
IF(SIZE(cells)==0)RETURN
DO m=1,SIZE(cells)
  c=c+SUM(mesh%r(:,mesh%lc(:,cells(m))),DIM=2)/3.d0
END DO
c=c/REAL(SIZE(cells),8)
DO m=1,SIZE(cells)
  DO n=1,3
    e=MAX(e,magnitude(mesh%r(:,mesh%lc(n,cells(m)))-c))
  END DO
END DO
END SUBROUTINE cell_cluster
!---Direct sum V(:,j) = sum_i c_i K(j,i) for the cells of a group
SUBROUTINE group_direct(gcells,scells,c,Vout)
INTEGER(4), INTENT(in) :: gcells(:),scells(:)
REAL(8), INTENT(in) :: c(:,:)
REAL(8), INTENT(inout) :: Vout(:,:)
INTEGER(4) :: m,n
DO m=1,SIZE(gcells)
  Vout(:,gcells(m))=0.d0
  DO n=1,SIZE(scells)
    Vout(:,gcells(m))=Vout(:,gcells(m))+c(:,n)*tw_aca_cell_kernel(mesh,gcells(m),scells(n))
  END DO
END DO
END SUBROUTINE group_direct
!---ACA+ (tw_aca_plus) of K(group,hole cells), then V = U (W c^T); .FALSE. if not converged at low rank
FUNCTION group_aca(gcells,scells,c,Vout) RESULT(success)
INTEGER(4), TARGET, INTENT(in) :: gcells(:),scells(:)
REAL(8), INTENT(in) :: c(:,:)
REAL(8), INTENT(inout) :: Vout(:,:)
LOGICAL :: success
INTEGER(4) :: m,n,kk,ir
REAL(8), ALLOCATABLE :: U(:,:),W(:,:),y(:,:)
TYPE(tw_aca_cellblock) :: src
m=SIZE(gcells); n=SIZE(scells)
src%mesh=>mesh
src%rcells=>gcells
src%ccells=>scells
!---Relative tolerance, at most MIN(m,n)/3 terms
CALL tw_aca_plus(src,m,n,tw_hole_aca_tol,tw_hole_aca_rmin,MIN(m,n)/3+1,1,1,U,W,kk,rel_tol=.TRUE.)
success=(kk>0)
IF(success)THEN
  ALLOCATE(y(kk,3))
  y=MATMUL(W,TRANSPOSE(c))
  DO ir=1,m
    Vout(:,gcells(ir))=MATMUL(U(ir,:),y)
  END DO
  DEALLOCATE(U,W,y)
END IF
END FUNCTION group_aca
END SUBROUTINE tw_hodlr_LmatHole
!---------------------------------------------------------------------------------
!> Compute mutual inductance matrix between two thin-wall models
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_Lmatblock(row_obj,col_obj,Lmat,row_block,col_block,sym)
TYPE(tw_type), INTENT(in) :: row_obj !< Thin-wall model object for rows
TYPE(tw_type), INTENT(in) :: col_obj !< Thin-wall model object for columns
TYPE(oft_tw_block), INTENT(in) :: row_block !< Block of rows to compute
TYPE(oft_tw_block), INTENT(in) :: col_block !< Block of columns to compute
REAL(8), INTENT(out) :: Lmat(:,:) !< Mutual inductance matrix
LOGICAL, OPTIONAL, INTENT(in) :: sym !< Row and column blocks are the same (diagonal block): compute each cell pair once
INTEGER(4) :: i,ii,j,jj,ik,jk,irow,jcol,jstart
LOGICAL :: is_sym
INTEGER(8) :: counts(4)
REAL(8) :: tmp,dl_min,f(3)
REAL(8) :: area_i,evec_i(3,3),pts_i(3,3)
REAL(8) :: area_j,evec_j(3,3),pts_j(3,3)
CLASS(oft_bmesh), POINTER :: rmesh,cmesh
DEBUG_STACK_PUSH
!
! WRITE(*,*)'Building element<->element inductance matrix block'
Lmat=0.d0
! CALL mytimer%tick
!---Setup quadrature
!
rmesh=>row_obj%mesh
cmesh=>col_obj%mesh
!---
counts=0
f=1.d0/3.d0
is_sym=.FALSE.
IF(PRESENT(sym))is_sym=sym
jstart=1
DO irow=1,row_block%ncells
    i = row_block%icell(irow)
    ! CALL rmesh%jacobian(i,f,rgop,area_i)
    ! CALL rmesh%norm(i,f,norm_i)
    area_i=rmesh%ca(i)
    DO ii=1,3
    pts_i(:,ii)=rmesh%r(:,rmesh%lc(ii,i))
    ! evec_i(:,ii)=cross_product(rgop(:,ii),norm_i)
    evec_i(:,ii)=row_obj%qbasis(:,ii,i)
    END DO
    IF(is_sym)jstart=irow
    !---Compute inter-edge inductances
    DO jcol=jstart,col_block%ncells
    j = col_block%icell(jcol)
    ! CALL cmesh%jacobian(j,f,cgop,area_j)
    ! CALL cmesh%norm(j,f,norm_j)
    area_j=cmesh%ca(j)
    DO jj=1,3
        pts_j(:,jj)=cmesh%r(:,cmesh%lc(jj,j))
        ! evec_j(:,jj)=cross_product(cgop(:,jj),norm_j)
        evec_j(:,jj)=col_obj%qbasis(:,jj,j)
    END DO
    !---Cell-pair integral of 1/|r-r'|
    tmp=tw_cell_pair_kernel(pts_i,area_i,pts_j,area_j)
    !
    DO ii=1,3
        ik=row_block%inv_map(rmesh%lc(ii,i))
        IF(ik==0)CYCLE
        DO jj=1,3
        jk=col_block%inv_map(cmesh%lc(jj,j))
        IF(jk==0)CYCLE
        dl_min = DOT_PRODUCT(evec_i(:,ii),evec_j(:,jj))*tmp
        ! !$omp atomic
        Lmat(jk,ik) = Lmat(jk,ik) + dl_min
        IF(is_sym.AND.(jcol/=irow))Lmat(ik,jk) = Lmat(ik,jk) + dl_min ! Mirrored pair (j,i)
        END DO
    END DO
    END DO
END DO
Lmat = Lmat/(4.d0*pi)
! elapsed_time=mytimer%tock()
! WRITE(*,*)'  Time = ',elapsed_time
DEBUG_STACK_POP
END SUBROUTINE tw_compute_Lmatblock
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_Bops_hole(self,Bop,Bop_dr,B_dx)
TYPE(tw_type), INTENT(inout) :: self
REAL(8), INTENT(out) :: Bop(:,:,:) !< Magnetic field evaluation matrix
REAL(8), INTENT(out) :: Bop_dr(:,:,:) !< Magnetic field evaluation matrix
REAL(8), INTENT(in) :: B_dx !< Spatial step size for finite difference evaluation of B-field
REAL(r8) :: evec_i(3,3),pts_i(3,3),pt_j(3),diffvec(3),ecc(3)
REAL(r8) :: cvec(3),cpt(3),tmp,area_i,norm_j(3),f(3)
REAL(r8), ALLOCATABLE :: atmp(:,:,:)
INTEGER(4) :: i,ii,j,jj,ik,k,kk
CLASS(oft_bmesh), POINTER :: bmesh
bmesh=>self%mesh
Bop=0.d0
f=1.d0/3.d0
! WRITE(*,*)'Building element->element magnetic reconstruction operator'
!$omp parallel private(ii,j,jj,ik,pts_i,tmp,pt_j,evec_i, &
!$omp atmp,i,area_i,norm_j,diffvec)
ALLOCATE(atmp(3,3,bmesh%np))
!$omp do schedule(dynamic,100)
DO i=1,bmesh%nc
  IF(self%kfh(i+1)-self%kfh(i)==0)CYCLE
  ! CALL bmesh%jacobian(i,f,rgop,area_i)
  ! CALL bmesh%norm(i,f,norm_i)
  area_i=bmesh%ca(i)
  DO ii=1,3
    pts_i(:,ii)=bmesh%r(:,bmesh%lc(ii,i))
    ! evec_i(:,ii)=cross_product(rgop(:,ii),norm_i)
    evec_i(:,ii)=self%qbasis(:,ii,i)
  END DO
  CALL bmesh%norm(i,f,norm_j)
  !---Compute sensor couplings
  atmp=0.d0
  DO j=1,bmesh%np
    !---Field from the cell's three basis functions
    atmp(:,:,j)=tw_cell_point_Bkernel(pts_i,evec_i,area_i,norm_j,bmesh%r(:,j),bmesh%va(j),B_dx)
  END DO
  DO ii=self%kfh(i),self%kfh(i+1)-1
    ik=ABS(self%lfh(1,ii))!+self%np_active
    DO j=1,bmesh%np
      DO jj=1,3
        tmp=SIGN(1,self%lfh(1,ii))*atmp(jj,self%lfh(2,ii),j)
        !$omp atomic
        Bop(j,ik,jj) = Bop(j,ik,jj) + tmp
      END DO
    END DO
  END DO
END DO
DEALLOCATE(atmp)
!$omp end parallel
!
! WRITE(*,*)'Building vcoil->element magnetic reconstruction operator'
IF(self%n_vcoils>0)THEN
  !$omp parallel do private(ii,j,k,kk,pt_j,ecc,diffvec,cvec,cpt)
  DO i=1,bmesh%np
    pt_j=bmesh%r(:,i)
    !---Compute driver contributions
    DO j=1,self%n_vcoils
      ecc = 0.d0
      IF(self%vcoils(j)%sens_mask)CYCLE
      DO k=1,self%vcoils(j)%ncoils
        IF(self%vcoils(j)%coils(k)%is_ring)THEN
          diffvec=tw_ring_bvec(self%vcoils(j)%coils(k)%ring_rz,pt_j)
        ELSE
          diffvec=0.d0
          DO kk=2,self%vcoils(j)%coils(k)%npts
            cvec = self%vcoils(j)%coils(k)%pts(:,kk)-self%vcoils(j)%coils(k)%pts(:,kk-1)
            cpt = (self%vcoils(j)%coils(k)%pts(:,kk)+self%vcoils(j)%coils(k)%pts(:,kk-1))/2.d0
            diffvec = diffvec + cross_product(cvec,pt_j-cpt)/SUM((pt_j-cpt)**2)**1.5d0
          END DO
        END IF
        ecc=ecc+self%vcoils(j)%scales(k)*diffvec
      END DO
      DO jj=1,3
        Bop(i,self%nholes+j,jj) = &
          Bop(i,self%nholes+j,jj) + ecc(jj)
      END DO
    END DO
  END DO
END IF
Bop=Bop/(4.d0*pi)
!
Bop_dr=0.d0
! WRITE(*,*)'Building icoil->element magnetic reconstruction operator'
IF(self%n_icoils>0)THEN
  !$omp parallel do private(ii,j,k,kk,pt_j,ecc,diffvec,cvec,cpt)
  DO i=1,bmesh%np
    pt_j=bmesh%r(:,i)
    !---Compute driver contributions
    DO j=1,self%n_icoils
      ecc = 0.d0
      IF(self%icoils(j)%sens_mask)CYCLE
      DO k=1,self%icoils(j)%ncoils
        IF(self%icoils(j)%coils(k)%is_ring)THEN
          diffvec=tw_ring_bvec(self%icoils(j)%coils(k)%ring_rz,pt_j)
        ELSE
          diffvec=0.d0
          DO kk=2,self%icoils(j)%coils(k)%npts
            cvec = self%icoils(j)%coils(k)%pts(:,kk)-self%icoils(j)%coils(k)%pts(:,kk-1)
            cpt = (self%icoils(j)%coils(k)%pts(:,kk)+self%icoils(j)%coils(k)%pts(:,kk-1))/2.d0
            diffvec = diffvec + cross_product(cvec,pt_j-cpt)/SUM((pt_j-cpt)**2)**1.5d0
          END DO
        END IF
        ecc=ecc+self%icoils(j)%scales(k)*diffvec
      END DO
      DO jj=1,3
        Bop_dr(i,j,jj) = Bop_dr(i,j,jj) + ecc(jj)
      END DO
    END DO
  END DO
END IF
Bop_dr=Bop_dr*mu0/(4.d0*pi)
END SUBROUTINE tw_compute_Bops_hole
!---------------------------------------------------------------------------------
!> Compute all three directions of a block of the B-field operator in one pass
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_Bops_block3(self,Bop,row_block,col_block,B_dx)
TYPE(tw_type), INTENT(inout) :: self
REAL(8), INTENT(out) :: Bop(:,:,:) !< Magnetic field evaluation matrix [col_block%np,row_block%nelems,3]
TYPE(oft_tw_block), INTENT(in) :: row_block !< Block of rows to compute
TYPE(oft_tw_block), INTENT(in) :: col_block !< Block of columns to compute
REAL(8), INTENT(in) :: B_dx !< Spatial step size for finite difference evaluation of B-field
REAL(r8) :: evec_i(3,3),pts_i(3,3)
REAL(r8) :: area_i,norm_j(3),f(3)
REAL(r8), ALLOCATABLE :: atmp(:,:,:)
INTEGER(4) :: i,ii,j,jj,ik,irow,jcol
CLASS(oft_bmesh), POINTER :: bmesh
bmesh=>self%mesh
Bop=0.d0
f=1.d0/3.d0
ALLOCATE(atmp(3,3,col_block%np))
DO irow=1,row_block%ncells
    i = row_block%icell(irow)
    area_i=bmesh%ca(i)
    DO ii=1,3
    pts_i(:,ii)=bmesh%r(:,bmesh%lc(ii,i))
    evec_i(:,ii)=self%qbasis(:,ii,i)
    END DO
    CALL bmesh%norm(i,f,norm_j)
    !---Compute sensor couplings
    atmp=0.d0
    DO jcol=1,col_block%np
    j=col_block%ipts(jcol)
    !---Field from the cell's three basis functions
    atmp(:,:,jcol)=tw_cell_point_Bkernel(pts_i,evec_i,area_i,norm_j,bmesh%r(:,j),bmesh%va(j),B_dx)
    END DO
    DO ii=1,3
        ik=row_block%inv_map(bmesh%lc(ii,i))
        IF(ik==0)CYCLE
        DO jj=1,3
            DO j=1,col_block%np
                Bop(j,ik,jj) = Bop(j,ik,jj) + atmp(jj,ii,j)
            END DO
        END DO
    END DO
END DO
Bop=Bop/(4.d0*pi)
DEALLOCATE(atmp)
END SUBROUTINE tw_compute_Bops_block3
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE tw_hodlr_setup(self,required)
class(oft_tw_hodlr_op), intent(inout) :: self
LOGICAL, INTENT(in) :: required
INTEGER(4) :: i,j,k,l,nlast,branch_count,pind_i,pind_j,max_pts,top_level,level,io_unit,ierr
INTEGER(4) :: size_out,nblocks_progress(9),nblocks_complete,counts(4),iblock,jblock,sparse_count(2)
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: kc_tmp,kr_tmp
INTEGER(4), ALLOCATABLE, DIMENSION(:,:) :: block_interaction,row_count
INTEGER(8) :: mat_size(3),compressed_size
LOGICAL :: is_masked,is_neighbor
LOGICAL, ALLOCATABLE, DIMENSION(:) :: near_closure,cell_mark
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: dense_blocks
REAL(8) :: corners_i(3,8),corners_j(3,8),rmin(3),elapsed_time,avg_size,L_approx_scale,B_approx_scale
REAL(8), ALLOCATABLE, DIMENSION(:) :: point_block
TYPE(oft_native_dense_matrix), POINTER :: mat_tmp => NULL()
type(oft_timer) :: mytimer
!---Tuning options
INTEGER(4) :: target_size = 1500
INTEGER(4) :: aca_min_its = 20
INTEGER(4) :: min_rank = 10
REAL(8) :: aca_sep_thresh(2) = [1.1d0,1.1d0]
REAL(8) :: L_svd_tol = -1.d99
REAL(8) :: L_aca_rel_tol = -1.d0
REAL(8) :: B_svd_tol = -1.d99
REAL(8) :: B_aca_rel_tol = -1.d0
LOGICAL :: pre_single = .FALSE.
LOGICAL :: pre_coarse = .FALSE.
INTEGER(4) :: pre_neigs = 0
NAMELIST/thincurr_hodlr_options/target_size,min_rank,aca_min_its,aca_sep_thresh,L_svd_tol, &
  L_aca_rel_tol,B_svd_tol,B_aca_rel_tol,pre_single,pre_coarse,pre_neigs
!
IF(self%tw_obj%kpmap_inv(self%tw_obj%np_active+1)/=self%tw_obj%np_active+1)THEN
  WRITE(*,*)self%tw_obj%kpmap_inv(self%tw_obj%np_active+1),self%tw_obj%np_active+1
  CALL oft_abort( &
    "HODLR not supported with periodic grids","tw_hodlr_setup",__FILE__)
END IF
!---
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit,thincurr_hodlr_options, IOSTAT=ierr)
CLOSE(io_unit)
if(ierr<0)THEN
  IF(required)THEN
    CALL oft_abort('No "thincurr_hodlr_options" found in input file.', &
      'tw_hodlr_setup',__FILE__)
  ELSE
    WRITE(*,'(2A)')oft_indent,'No "thincurr_hodlr_options" found in input file, using full representation'
    RETURN
  END IF
END IF
if(ierr>0)call oft_abort('Error parsing "thincurr_hodlr_options" in input file.', &
  'tw_hodlr_setup',__FILE__)
IF((L_svd_tol>0.d0).AND.(min_rank>aca_min_its))CALL oft_abort('"min_rank" must be <= "aca_min_its"', &
  "tw_hodlr_setup",__FILE__)
IF(min_rank<1)CALL oft_abort('"min_rank" must be >= 1',"tw_hodlr_setup",__FILE__)
!
self%leaf_target=target_size
self%min_rank=min_rank
self%aca_min_its=aca_min_its
self%aca_sep_thresh=aca_sep_thresh
self%pre_single=pre_single
self%pre_coarse=pre_coarse
self%pre_neigs=pre_neigs
!
L_approx_scale = 0.5d0*self%tw_obj%mesh%hrms
IF(L_svd_tol<0.d0)THEN
  IF(L_svd_tol<-1.d98)CALL oft_abort('L_svd_tol must be set in thincurr_hodlr_options', &
    'tw_hodlr_setup',__FILE__)
  L_svd_tol=ABS(L_svd_tol)*L_approx_scale
END IF
self%L_svd_tol=L_svd_tol
self%L_aca_tol=L_aca_rel_tol*self%L_svd_tol
!
B_approx_scale = 0.5d0/self%tw_obj%mesh%hrms
IF(B_svd_tol<0.d0)THEN
  IF(B_svd_tol<-1.d98)CALL oft_abort('B_svd_tol must be set in thincurr_hodlr_options', &
    'tw_hodlr_setup',__FILE__)
  B_svd_tol=ABS(B_svd_tol)*B_approx_scale
END IF
self%B_svd_tol=B_svd_tol
self%B_aca_tol=B_aca_rel_tol*self%B_svd_tol
!---Partition mesh
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Partitioning grid for block low rank compressed operators'
CALL oft_increase_indent
CALL tw_part_mesh(self%tw_obj,self%leaf_target,self%nlevels,self%levels)
sparse_count=0
DO i=1,self%nlevels
  ! WRITE(*,*)i,self%levels(i)%nblocks
  ALLOCATE(self%levels(i)%mat_mask(self%levels(i)%nblocks,self%levels(i)%nblocks))
  self%levels(i)%mat_mask = 0
  IF(i>1)THEN
    DO j=1,self%levels(i)%nblocks
      DO k=1,self%levels(i)%nblocks
        IF(self%levels(i-1)%mat_mask(self%levels(i)%blocks(j)%parent,self%levels(i)%blocks(k)%parent)/=0)THEN
          ! WRITE(*,*)'Masking',i,j,k
          self%levels(i)%mat_mask(j,k)=-1
          self%levels(i)%mat_mask(k,j)=-1
        END IF
      END DO
    END DO
  END IF
  DO j=1,self%levels(i)%nblocks
    DO k=j+1,self%levels(i)%nblocks
      ! WRITE(*,*)i,j,k,magnitude(self%levels(i)%blocks(j)%center-self%levels(i)%blocks(k)%center), &
      ! (self%levels(i)%blocks(j)%extent+self%levels(i)%blocks(k)%extent)
      IF(self%levels(i)%mat_mask(j,k)/=0)CYCLE
      IF(magnitude(self%levels(i)%blocks(j)%center-self%levels(i)%blocks(k)%center)/ &
        (self%levels(i)%blocks(j)%extent+self%levels(i)%blocks(k)%extent) > self%aca_sep_thresh(2))THEN
        self%levels(i)%mat_mask(j,k)=2
        self%levels(i)%mat_mask(k,j)=-2
      ELSE IF(i==self%nlevels)THEN
        IF(magnitude(self%levels(i)%blocks(j)%center-self%levels(i)%blocks(k)%center)/ &
          (self%levels(i)%blocks(j)%extent+self%levels(i)%blocks(k)%extent) > self%aca_sep_thresh(1))THEN
          self%levels(i)%mat_mask(j,k)=2
          self%levels(i)%mat_mask(k,j)=-2
        ELSE
          self%levels(i)%mat_mask(j,k)=3
          self%levels(i)%mat_mask(k,j)=-3
        END IF
      END IF
    END DO
    IF(i==self%nlevels)self%levels(i)%mat_mask(j,j)=1
  END DO
  counts(1) = SUM(COUNT(self%levels(i)%mat_mask==0,DIM=1))
  counts(2) = SUM(COUNT(self%levels(i)%mat_mask==1,DIM=1))
  counts(3) = SUM(COUNT(self%levels(i)%mat_mask==2,DIM=1))
  counts(4) = SUM(COUNT(self%levels(i)%mat_mask==3,DIM=1))
  ! WRITE(*,*)'  ',counts
  self%ndense=self%ndense+counts(2)
  self%nsparse=self%nsparse+counts(3)+counts(4)
  sparse_count=sparse_count+counts(3:4)
END DO
! WRITE(*,*)self%ndense,self%nsparse
ALLOCATE(self%dense_blocks(3,self%ndense),self%sparse_blocks(3,self%nsparse))
self%ndense=0
self%nsparse=0
DO i=1,self%nlevels
  DO j=1,self%levels(i)%nblocks
    DO k=1,self%levels(i)%nblocks
      IF(self%levels(i)%mat_mask(j,k)==1)THEN
        self%ndense=self%ndense+1
        self%dense_blocks(:,self%ndense)=[i,j,k]
      ELSE IF(self%levels(i)%mat_mask(j,k)>=2)THEN
        self%nsparse=self%nsparse+1
        self%sparse_blocks(:,self%nsparse)=[i,j,k]
      END IF
    END DO
  END DO
END DO
! WRITE(*,*)self%ndense,self%nsparse
! CALL oft_abort("","",__FILE__)
self%nblocks=self%levels(self%nlevels)%nblocks
self%blocks=>self%levels(self%nlevels)%blocks
WRITE(*,'(2A,I6)')oft_indent,'nBlocks =        ',self%nblocks
ALLOCATE(point_block(self%tw_obj%mesh%np))
point_block=0.d0
avg_size=0.d0
DO i=1,self%nblocks
  avg_size=avg_size+self%blocks(i)%nelems
  ! DO j=1,self%blocks(i)%nelems
  !     point_block(self%tw_obj%lpmap_inv(self%blocks(i)%ielem))=i*1.d0
  ! END DO
  DO j=1,self%blocks(i)%np
    point_block(self%blocks(i)%ipts)=i*1.d0
  END DO
END DO
WRITE(*,'(2A,I6)')oft_indent,'Avg block size = ',INT(avg_size/REAL(self%nblocks,8),4)
IF(self%L_aca_tol<=0.d0)THEN
  sparse_count=[0,sparse_count(1)+sparse_count(2)]
END IF
WRITE(*,'(2A,I6)')oft_indent,'# of SVD  =      ',sparse_count(2)
WRITE(*,'(2A,I6)')oft_indent,'# of ACA+ =      ',sparse_count(1)
IF(self%pre_single)WRITE(*,'(2A)')oft_indent,'Using single precision block preconditioners'
IF(self%pre_coarse)WRITE(*,'(2A)')oft_indent,'Using two-level block preconditioners'
!---Diagnostic block field (only when plot output has been set up)
IF(self%tw_obj%xdmf%n_grids>0)CALL self%tw_obj%mesh%save_vertex_scalar(point_block,self%tw_obj%xdmf,'ACA_Block')
ALLOCATE(cell_mark(self%tw_obj%mesh%nc))
DO j=1,self%nlevels
  DO i=1,self%levels(j)%nblocks
    IF(oft_debug_print(1))WRITE(*,'(A,I4,I6)')oft_indent,i,self%levels(j)%blocks(i)%nelems
    !
    IF(ASSOCIATED(self%levels(j)%blocks(i)%inv_map))CYCLE
    ALLOCATE(self%levels(j)%blocks(i)%inv_map(self%tw_obj%mesh%np))
    self%levels(j)%blocks(i)%inv_map = 0
    DO l=1,self%levels(j)%blocks(i)%nelems
      k=self%tw_obj%lpmap_inv(self%levels(j)%blocks(i)%ielem(l))
      self%levels(j)%blocks(i)%inv_map(k)=l
    END DO
    !
    cell_mark=.FALSE.
    DO l=1,self%levels(j)%blocks(i)%nelems
      pind_j = self%tw_obj%lpmap_inv(self%levels(j)%blocks(i)%ielem(l))
      DO k=self%tw_obj%mesh%kpc(pind_j),self%tw_obj%mesh%kpc(pind_j+1)-1
        cell_mark(self%tw_obj%mesh%lpc(k))=.TRUE.
      END DO
    END DO
    self%levels(j)%blocks(i)%ncells=COUNT(cell_mark)
    ALLOCATE(self%levels(j)%blocks(i)%icell(self%levels(j)%blocks(i)%ncells))
    self%levels(j)%blocks(i)%ncells=0
    DO l=1,self%tw_obj%mesh%nc
      IF(cell_mark(l))THEN
        self%levels(j)%blocks(i)%ncells=self%levels(j)%blocks(i)%ncells+1
        self%levels(j)%blocks(i)%icell(self%levels(j)%blocks(i)%ncells)=l
      END IF
    END DO
  END DO
END DO
DEALLOCATE(cell_mark)
CALL oft_decrease_indent
END SUBROUTINE tw_hodlr_setup
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE tw_hodlr_Lcompute(self,save_file)
class(oft_tw_hodlr_op), target, intent(inout) :: self
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: save_file
INTEGER(4) :: i,j,k,l,nlast,branch_count,pind_i,pind_j,max_pts,top_level,level
INTEGER(4) :: nblocks_progress(9),nblocks_complete,counts(4),iblock,jblock
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: kc_tmp,kr_tmp
INTEGER(4), ALLOCATABLE, DIMENSION(:,:) :: block_interaction,row_count
INTEGER(8) :: mat_size(3),compressed_size,size_out
INTEGER(8), ALLOCATABLE, DIMENSION(:) :: dense_sizes,sparse_sizes
LOGICAL :: is_masked,is_neighbor,mat_updated
LOGICAL, ALLOCATABLE, DIMENSION(:) :: near_closure,cell_mark
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: dense_blocks
REAL(8) :: corners_i(3,8),corners_j(3,8),rmin(3),elapsed_time,avg_size,curr_size,target_size
REAL(8), ALLOCATABLE, DIMENSION(:) :: point_block
TYPE(oft_native_dense_matrix), POINTER :: mat_tmp => NULL()
type(oft_timer) :: mytimer
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Building block low rank inductance operator'
CALL oft_increase_indent
WRITE(*,'(2A,ES9.2)')oft_indent,'SVD tolerance:  ',self%L_svd_tol
WRITE(*,'(2A,ES9.2)')oft_indent,'ACA+ tolerance: ',self%L_aca_tol
IF(self%L_svd_tol<0.d0)CALL oft_abort('"L_svd_tol" must be > 0','tw_Lmat_MF_Lcompute',__FILE__)
!---Build matrix with ACA+ or SVD compression of off-diagonal blocks
CALL mytimer%tick()
ALLOCATE(self%dense_mats(self%ndense))
ALLOCATE(self%aca_U_mats(self%nsparse))
ALLOCATE(self%aca_V_mats(self%nsparse))
ALLOCATE(self%aca_dense(self%nsparse))
WRITE(*,'(2A)')oft_indent,'Building hole and Vcoil columns'
IF(self%tw_obj%nholes+self%tw_obj%n_vcoils>0)THEN
  ALLOCATE(self%hole_Vcoil_mat%M(tw_nelems_exp(self%tw_obj),self%tw_obj%nholes+self%tw_obj%n_vcoils))
  CALL tw_hodlr_LmatHole(self,self%hole_Vcoil_mat%M)
END IF
IF(PRESENT(save_file))CALL load_from_file()
compressed_size=0
mat_updated=.FALSE.
WRITE(*,'(2A)')oft_indent,'Building diagonal blocks'
ALLOCATE(dense_sizes(self%ndense),sparse_sizes(self%nsparse))
dense_sizes=0
sparse_sizes=0
!$omp parallel private(i,level,j,k,size_out,avg_size) reduction(+:compressed_size) &
!$omp reduction(.OR.:mat_updated)
!$omp single
nblocks_progress = INT(self%ndense*[0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0],4)
nblocks_complete = 0
!$omp end single
!$omp do schedule(dynamic,1)
DO i=1,self%ndense
  level = self%dense_blocks(1,i)
  j = self%dense_blocks(2,i)
  k = self%dense_blocks(3,i)
  !---Build full representation
  IF(.NOT.ASSOCIATED(self%dense_mats(i)%M))THEN
    mat_updated=.TRUE.
    ALLOCATE(self%dense_mats(i)%M(self%levels(level)%blocks(j)%nelems,self%levels(level)%blocks(k)%nelems))
    CALL tw_compute_Lmatblock(self%tw_obj,self%tw_obj,self%dense_mats(i)%M, &
      self%levels(level)%blocks(k),self%levels(level)%blocks(j),sym=(j==k))
  END IF
  dense_sizes(i) = self%levels(level)%blocks(j)%nelems*self%levels(level)%blocks(k)%nelems
  compressed_size = compressed_size + dense_sizes(i)
  !$omp critical
  nblocks_complete = nblocks_complete + 1
  DO k=1,9
    IF(nblocks_complete==nblocks_progress(k))WRITE(*,'(A,2X,I2,A1)')oft_indent,k*10,'%'
  END DO
  !$omp end critical
END DO
!$omp single
IF(self%L_aca_tol>0.d0)THEN
  WRITE(*,'(2A)')oft_indent,'Building off-diagonal blocks using ACA+ and SVD compression'
ELSE
  WRITE(*,'(2A)')oft_indent,'Building off-diagonal blocks with SVD compression'
END IF
nblocks_progress = INT(self%nsparse*[0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0],4)
nblocks_complete = 0
!$omp end single
!$omp do schedule(dynamic,1)
DO i=1,self%nsparse
  level = self%sparse_blocks(1,i)
  j = self%sparse_blocks(2,i)
  k = self%sparse_blocks(3,i)
  !---Build full representation
  IF(ASSOCIATED(self%aca_dense(i)%M))THEN
    size_out=self%levels(level)%blocks(k)%nelems*self%levels(level)%blocks(j)%nelems
  ELSE IF(ASSOCIATED(self%aca_U_mats(i)%M))THEN
    size_out=SIZE(self%aca_V_mats(i)%M,DIM=1)*(self%levels(level)%blocks(k)%nelems+self%levels(level)%blocks(j)%nelems)
  ELSE
    mat_updated=.TRUE.
    size_out=-1
    IF((self%L_aca_tol>0.d0).AND.self%levels(level)%mat_mask(j,k)==2)THEN
      CALL aca_approx(i,self%L_aca_tol,size_out)
      IF(size_out>0)CALL tw_lowrank_recompress(self%aca_U_mats(i),self%aca_V_mats(i),self%L_svd_tol,self%min_rank, &
        self%aca_min_its,size_out)
      IF(size_out==-1)WRITE(*,'(2A,2I6)')oft_indent,'ACA+ failure',self%levels(level)%blocks(j)%nelems, &
        self%levels(level)%blocks(k)%nelems
    END IF
    !---IF ACA+ failed or diagonal compute dense matrix
    IF(size_out<0)THEN
      ALLOCATE(self%aca_dense(i)%M(self%levels(level)%blocks(k)%nelems,self%levels(level)%blocks(j)%nelems))
      CALL tw_compute_Lmatblock(self%tw_obj,self%tw_obj,self%aca_dense(i)%M, &
        self%levels(level)%blocks(j),self%levels(level)%blocks(k))
      CALL compress_block(i,self%L_svd_tol,size_out)
    END IF
  END IF
  sparse_sizes(i) = size_out
  compressed_size = compressed_size + sparse_sizes(i)
  !$omp critical
  nblocks_complete = nblocks_complete + 1
  DO k=1,9
    IF(nblocks_complete==nblocks_progress(k))WRITE(*,'(A,2X,I2,A1)')oft_indent,k*10,'%'
  END DO
  !$omp end critical
END DO
!$omp end parallel
elapsed_time=mytimer%tock()
WRITE(*,'(2A,F6.1,A,ES9.2,A,ES9.2,A)')oft_indent,'Compressed size:',compressed_size*1.d2/(REAL(self%tw_obj%np_active,8)**2), &
    "%  (",REAL(compressed_size,8),"/",REAL(self%tw_obj%np_active,8)**2,")"
WRITE(*,'(3A)')oft_indent,'Time = ',time_to_string(elapsed_time)
!---Building thread blocking
ALLOCATE(self%L_thread_ptr(2,oft_env%nthreads+1))
self%L_thread_ptr=0
self%L_thread_ptr(:,1)=1
target_size=SUM(sparse_sizes)/REAL(oft_env%nthreads,8)
curr_size=0.d0
i=1
DO j=1,self%nsparse
  curr_size=curr_size+sparse_sizes(j)
  IF(curr_size>REAL(i,8)*target_size)THEN
    self%L_thread_ptr(1,i+1)=j+1
    i=i+1
  END IF
END DO
target_size=SUM(dense_sizes)/REAL(oft_env%nthreads,8)
curr_size=0.d0
i=1
DO j=1,self%ndense
  curr_size=curr_size+dense_sizes(j)
  IF(curr_size>REAL(i,8)*target_size)THEN
    self%L_thread_ptr(2,i+1)=j+1
    i=i+1
  END IF
END DO
DO j=1,oft_env%nthreads
  IF(self%L_thread_ptr(1,j+1)==0)self%L_thread_ptr(1,j+1)=self%nsparse+1
  IF(self%L_thread_ptr(2,j+1)==0)self%L_thread_ptr(2,j+1)=self%ndense+1
END DO
self%L_thread_ptr(:,oft_env%nthreads+1)=[self%nsparse,self%ndense]+1
DEALLOCATE(sparse_sizes,dense_sizes)
!---Save to file if requested
IF(mat_updated.AND.PRESENT(save_file))CALL save_to_file()
CALL oft_decrease_indent
CONTAINS
SUBROUTINE save_to_file()
INTEGER(4) :: ierr,io_unit,hash_tmp(6),file_counts(6)
CHARACTER(LEN=5) :: matrix_id
TYPE(HODLR_hash), TARGET :: sparse_hash
!---Save computed matrix to file
IF(TRIM(save_file)/='none')THEN
  hash_tmp(1) = tw_nelems_exp(self%tw_obj)
  hash_tmp(2) = self%tw_obj%mesh%nc
  hash_tmp(3) = self%ndense
  sparse_hash%nsparse = self%nsparse
  sparse_hash%aca_tol = REAL(self%L_aca_tol,c_float)
  sparse_hash%svd_tol = REAL(self%L_svd_tol,c_float)
  hash_tmp(4) = oft_simple_hash(C_LOC(sparse_hash),INT(c_sizeof(sparse_hash),8),starting_value=tw_cache_ver)
  DO i=1,self%nlevels
    hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(i)%mat_mask),INT(4*self%levels(i)%nblocks*self%levels(i)%nblocks,8),starting_value=hash_tmp(4))
  END DO
  hash_tmp(5) = oft_simple_hash(C_LOC(self%tw_obj%mesh%lc),INT(4*3*self%tw_obj%mesh%nc,8))
  hash_tmp(6) = oft_simple_hash(C_LOC(self%tw_obj%mesh%r),INT(8*3*self%tw_obj%mesh%np,8))
  WRITE(*,'(3A)')oft_indent,'Saving HODLR matrix to file: ',TRIM(save_file)
  CALL hdf5_create_file(TRIM(save_file))
  CALL hdf5_write(hash_tmp,TRIM(save_file),'MODEL_hash')
  !
  hash_tmp(6)=0
  DO i=1,self%ndense
    level = self%dense_blocks(1,i)
    j = self%dense_blocks(2,i)
    k = self%dense_blocks(3,i)
    hash_tmp(1:3) = [level,j,k]
    hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ielem),INT(4*self%levels(level)%blocks(j)%nelems,8))
    hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ielem),INT(4*self%levels(level)%blocks(k)%nelems,8))
    WRITE(matrix_id,'(I5)')i
    CALL hdf5_write(hash_tmp,TRIM(save_file),'DENSE_hash_'//TRIM(ADJUSTL(matrix_id)))
    CALL hdf5_write(self%dense_mats(i)%M,TRIM(save_file),'DENSE_'//TRIM(ADJUSTL(matrix_id)))
  END DO
  !
  DO i=1,self%nsparse
    level = self%sparse_blocks(1,i)
    j = self%sparse_blocks(2,i)
    k = self%sparse_blocks(3,i)
    hash_tmp(1:3) = [level,j,k]
    hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ielem),INT(4*self%levels(level)%blocks(j)%nelems,8))
    hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ielem),INT(4*self%levels(level)%blocks(k)%nelems,8))
    WRITE(matrix_id,'(I5)')i
    IF(ASSOCIATED(self%aca_dense(i)%M))THEN
      hash_tmp(6)=-1
      CALL hdf5_write(hash_tmp,TRIM(save_file),'SPARSE_hash_'//TRIM(ADJUSTL(matrix_id)))
      CALL hdf5_write(self%aca_dense(i)%M,TRIM(save_file),'SPARSE_'//TRIM(ADJUSTL(matrix_id)))
    ELSE
      hash_tmp(6)=SIZE(self%aca_V_mats(i)%M,DIM=1)
      CALL hdf5_write(hash_tmp,TRIM(save_file),'SPARSE_hash_'//TRIM(ADJUSTL(matrix_id)))
      CALL hdf5_write(self%aca_U_mats(i)%M,TRIM(save_file),'SPARSE_U_'//TRIM(ADJUSTL(matrix_id)))
      CALL hdf5_write(self%aca_V_mats(i)%M,TRIM(save_file),'SPARSE_V_'//TRIM(ADJUSTL(matrix_id)))
    END IF
  END DO
  ! TODO: Save create appropriate hash to enable saving Vcoil matrix
  ! WRITE(io_unit)hash_tmp
  ! WRITE(io_unit)self%hole_Vcoil_mat%M
  ! CLOSE(io_unit)
END IF
END SUBROUTINE save_to_file
!
SUBROUTINE load_from_file()
LOGICAL :: exists
INTEGER(4) :: ierr,io_unit,hash_tmp(6),file_counts(6)
CHARACTER(LEN=5) :: matrix_id
TYPE(HODLR_hash), TARGET :: sparse_hash
!---Try to load from file
IF(TRIM(save_file)/='none')THEN
  INQUIRE(FILE=TRIM(save_file),EXIST=exists)
  IF(exists)THEN
    hash_tmp(1) = tw_nelems_exp(self%tw_obj)
    hash_tmp(2) = self%tw_obj%mesh%nc
    hash_tmp(3) = self%ndense
    sparse_hash%nsparse = self%nsparse
    sparse_hash%aca_tol = REAL(self%L_aca_tol,c_float)
    sparse_hash%svd_tol = REAL(self%L_svd_tol,c_float)
    hash_tmp(4) = oft_simple_hash(C_LOC(sparse_hash),INT(c_sizeof(sparse_hash),8),starting_value=tw_cache_ver)
    DO i=1,self%nlevels
      hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(i)%mat_mask),INT(4*self%levels(i)%nblocks*self%levels(i)%nblocks,8),starting_value=hash_tmp(4))
    END DO
    hash_tmp(5) = oft_simple_hash(C_LOC(self%tw_obj%mesh%lc),INT(4*3*self%tw_obj%mesh%nc,8))
    hash_tmp(6) = oft_simple_hash(C_LOC(self%tw_obj%mesh%r),INT(8*3*self%tw_obj%mesh%np,8))
    WRITE(*,'(3A)')oft_indent,'Reading HODLR matrix from file: ',TRIM(save_file)
    CALL oft_increase_indent
    CALL hdf5_read(file_counts,TRIM(save_file),'MODEL_hash',success=exists)
    IF(exists.AND.ALL(file_counts==hash_tmp))THEN
      hash_tmp(6)=0
      DO i=1,self%ndense
        level = self%dense_blocks(1,i)
        j = self%dense_blocks(2,i)
        k = self%dense_blocks(3,i)
        hash_tmp(1:3) = [level,j,k]
        hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ielem),INT(4*self%levels(level)%blocks(j)%nelems,8))
        hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ielem),INT(4*self%levels(level)%blocks(k)%nelems,8))
        WRITE(matrix_id,'(I5)')i
        CALL hdf5_read(file_counts,TRIM(save_file),'DENSE_hash_'//TRIM(ADJUSTL(matrix_id)),success=exists)
        IF(exists.AND.ALL(file_counts==hash_tmp))THEN
          ALLOCATE(self%dense_mats(i)%M(self%levels(level)%blocks(j)%nelems,self%levels(level)%blocks(k)%nelems))
          CALL hdf5_read(self%dense_mats(i)%M,TRIM(save_file),'DENSE_'//TRIM(ADJUSTL(matrix_id)),success=exists)
          IF(.NOT.exists)DEALLOCATE(self%dense_mats(i)%M)
        END IF
      END DO
      !
      DO i=1,self%nsparse
        level = self%sparse_blocks(1,i)
        j = self%sparse_blocks(2,i)
        k = self%sparse_blocks(3,i)
        hash_tmp(1:3) = [level,j,k]
        hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ielem),INT(4*self%levels(level)%blocks(j)%nelems,8))
        hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ielem),INT(4*self%levels(level)%blocks(k)%nelems,8))
        WRITE(matrix_id,'(I5)')i
        CALL hdf5_read(file_counts,TRIM(save_file),'SPARSE_hash_'//TRIM(ADJUSTL(matrix_id)),success=exists)
        IF(exists.AND.ALL(file_counts(1:5)==hash_tmp(1:5)))THEN
          IF(file_counts(6)<0)THEN
            ALLOCATE(self%aca_dense(i)%M(self%levels(level)%blocks(k)%nelems,self%levels(level)%blocks(j)%nelems))
            CALL hdf5_read(self%aca_dense(i)%M,TRIM(save_file),'SPARSE_'//TRIM(ADJUSTL(matrix_id)),success=exists)
            IF(.NOT.exists)DEALLOCATE(self%aca_dense(i)%M)
          ELSE
            ALLOCATE(self%aca_U_mats(i)%M(self%levels(level)%blocks(k)%nelems,file_counts(6)))
            CALL hdf5_read(self%aca_U_mats(i)%M,TRIM(save_file),'SPARSE_U_'//TRIM(ADJUSTL(matrix_id)),success=exists)
            IF(.NOT.exists)THEN
              DEALLOCATE(self%aca_U_mats(i)%M)
            ELSE
              ALLOCATE(self%aca_V_mats(i)%M(file_counts(6),self%levels(level)%blocks(j)%nelems))
              CALL hdf5_read(self%aca_V_mats(i)%M,TRIM(save_file),'SPARSE_V_'//TRIM(ADJUSTL(matrix_id)),success=exists)
              IF(.NOT.exists)DEALLOCATE(self%aca_U_mats(i)%M,self%aca_V_mats(i)%M)
            END IF
          END IF
        END IF
      END DO
    ELSE
      WRITE(*,'(3A)')oft_indent,'Ignoring stored matrix, Model hashes do not match'
      IF(oft_debug_print(1))THEN
        WRITE(*,*)'File: ',file_counts
        WRITE(*,*)'Model:',hash_tmp
      END IF
    END IF
    CALL oft_decrease_indent
  END IF
END IF
END SUBROUTINE load_from_file
!
SUBROUTINE compress_block(isparse,tol,size_out)
INTEGER(4), INTENT(in) :: isparse
REAL(8), INTENT(in) :: tol
INTEGER(8), INTENT(out) :: size_out
INTEGER(4) :: M,N,rank
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: U,V
!
M = SIZE(self%aca_dense(isparse)%M,DIM=1)
N = SIZE(self%aca_dense(isparse)%M,DIM=2)
CALL tw_lowrank_qrcp(self%aca_dense(isparse)%M,tol,self%min_rank,U,V,rank)
IF(rank>0)THEN
  ALLOCATE(self%aca_U_mats(isparse)%M(M,rank),self%aca_V_mats(isparse)%M(rank,N))
  self%aca_U_mats(isparse)%M=U
  self%aca_V_mats(isparse)%M=V
  DEALLOCATE(self%aca_dense(isparse)%M)
  size_out=INT(rank,8)*INT(M+N,8)
ELSE
  size_out=INT(M,8)*INT(N,8)
END IF
END SUBROUTINE compress_block
!
SUBROUTINE aca_approx(isparse,tol,size_out)
INTEGER(4), INTENT(in) :: isparse
REAL(8), INTENT(in) :: tol
INTEGER(8), INTENT(out) :: size_out
!
INTEGER(4) :: M,N,MIN_DIM,max_iter,Iref,Jref,nterms
INTEGER(8) :: full_size
REAL(8) :: tmp
REAL(8), ALLOCATABLE :: U(:,:),V(:,:)
TYPE(tw_aca_Lblock) :: src
!
src%op=>self
src%ilevel=self%sparse_blocks(1,isparse)
src%iblock=self%sparse_blocks(2,isparse)
src%jblock=self%sparse_blocks(3,isparse)
size_out=-2
M = self%levels(src%ilevel)%blocks(src%jblock)%nelems
N = self%levels(src%ilevel)%blocks(src%iblock)%nelems
MIN_DIM = min(M,N)
IF(MIN_DIM<=self%aca_min_its)RETURN
full_size = INT(M,8)*INT(N,8)
! If we haven't converged before running for max_iter, we'll stop anyway.
! We want the sum of the two matrices to always be smaller than the dense one
max_iter = MIN(self%leaf_target*2,INT(full_size/INT(M+N,8),4)-1)
size_out=-1
ALLOCATE(src%buf(MAX(M,N),1))
! Choose our starting random reference row and column.
CALL random_number(tmp)
Iref = INT(tmp*(N-1),4)+1
CALL random_number(tmp)
Jref = INT(tmp*(M-1),4)+1
!---Rows: elements of iblock, columns: elements of jblock; stored as M ~ U_ACA V_ACA (columns: iblock)
CALL tw_aca_plus(src,N,M,tol,self%aca_min_its,max_iter,Iref,Jref,U,V,nterms)
IF(nterms>0)THEN
  ALLOCATE(self%aca_U_mats(isparse)%M(M,nterms))
  self%aca_U_mats(isparse)%M=TRANSPOSE(V)
  ALLOCATE(self%aca_V_mats(isparse)%M(nterms,N))
  self%aca_V_mats(isparse)%M=TRANSPOSE(U)
  size_out=INT(nterms,8)*INT(M+N,8)
  DEALLOCATE(U,V)
END IF
DEALLOCATE(src%buf)
IF(ASSOCIATED(src%one%ielem))DEALLOCATE(src%one%ielem,src%one%icell,src%one%inv_map)
END SUBROUTINE aca_approx
END SUBROUTINE tw_hodlr_Lcompute
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE tw_hodlr_Bcompute(self,save_file)
class(oft_tw_hodlr_op), target, intent(inout) :: self
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: save_file
INTEGER(4) :: i,ii,j,k,l,nlast,branch_count,pind_i,pind_j,max_pts,top_level,level,dim,flip
INTEGER(4) :: nblocks_progress(9),nblocks_complete,counts(4),iblock,jblock
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: kc_tmp,kr_tmp
INTEGER(4), ALLOCATABLE, DIMENSION(:,:) :: block_interaction,row_count
INTEGER(8) :: mat_size(3),compressed_size,size_out
INTEGER(8), ALLOCATABLE, DIMENSION(:) :: dense_sizes,sparse_sizes
LOGICAL :: is_masked,is_neighbor,mat_updated
LOGICAL, ALLOCATABLE, DIMENSION(:) :: near_closure,cell_mark
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: dense_blocks
REAL(8) :: corners_i(3,8),corners_j(3,8),rmin(3),elapsed_time,avg_size,curr_size,target_size
REAL(8), ALLOCATABLE, DIMENSION(:) :: point_block
REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: Btmp
TYPE(oft_native_dense_matrix), POINTER :: mat_tmp => NULL()
type(oft_timer) :: mytimer
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Building block low rank magnetic field operator'
CALL oft_increase_indent
WRITE(*,'(2A,ES9.2)')oft_indent,'SVD tolerance:  ',self%B_svd_tol
WRITE(*,'(2A,ES9.2)')oft_indent,'ACA+ tolerance: ',self%B_aca_tol
IF(self%B_svd_tol<0.d0)CALL oft_abort('"B_svd_tol" must be > 0','tw_Lmat_MF_Bcompute',__FILE__)
!---Build matrix with SVD compression of off-diagonal blocks
! approximation with ACA+ planned
CALL mytimer%tick()
ALLOCATE(self%dense_B_mats(self%ndense,3))
ALLOCATE(self%aca_BU_mats(2*self%nsparse,3))
ALLOCATE(self%aca_BV_mats(2*self%nsparse,3))
ALLOCATE(self%aca_B_dense(2*self%nsparse,3))
!
WRITE(*,'(2A)')oft_indent,'Building hole, Vcoil and Icoil columns'
IF(MAX(self%tw_obj%n_icoils,self%tw_obj%nholes+self%tw_obj%n_vcoils)>0)THEN
  ALLOCATE(self%hole_Vcoil_Bmat(self%tw_obj%mesh%np,MAX(1,self%tw_obj%nholes+self%tw_obj%n_vcoils),3))
  ALLOCATE(self%Icoil_Bmat(self%tw_obj%mesh%np,MAX(1,self%tw_obj%n_icoils),3))
  CALL tw_compute_Bops_hole(self%tw_obj,self%hole_Vcoil_Bmat,self%Icoil_Bmat,self%tw_obj%B_dx)
  self%tw_obj%Bdr=>self%Icoil_Bmat
END IF
!
IF(PRESENT(save_file))CALL read_from_file()
WRITE(*,'(2A)')oft_indent,'Building diagonal blocks'
compressed_size=0
mat_updated=.FALSE.
ALLOCATE(dense_sizes(self%ndense),sparse_sizes(self%nsparse))
dense_sizes=0
sparse_sizes=0
!$omp parallel private(i,ii,level,j,k,size_out,avg_size,dim,flip,Btmp) reduction(+:compressed_size) &
!$omp reduction(.OR.:mat_updated)
!$omp single
nblocks_progress = INT(self%ndense*[0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0],4)
nblocks_complete = 0
!$omp end single
!$omp do schedule(dynamic,1)
DO i=1,self%ndense
  level = self%dense_blocks(1,i)
  j = self%dense_blocks(2,i)
  k = self%dense_blocks(3,i)
  !---Build full representation
  IF(.NOT.(ASSOCIATED(self%dense_B_mats(i,1)%M).AND.ASSOCIATED(self%dense_B_mats(i,2)%M) &
    .AND.ASSOCIATED(self%dense_B_mats(i,3)%M)))THEN
    mat_updated=.TRUE.
    ALLOCATE(Btmp(self%levels(level)%blocks(j)%np,self%levels(level)%blocks(k)%nelems,3))
    CALL tw_compute_Bops_block3(self%tw_obj,Btmp, &
      self%levels(level)%blocks(k),self%levels(level)%blocks(j),self%tw_obj%B_dx)
    DO dim=1,3
      IF(.NOT.ASSOCIATED(self%dense_B_mats(i,dim)%M))THEN
        ALLOCATE(self%dense_B_mats(i,dim)%M(self%levels(level)%blocks(j)%np,self%levels(level)%blocks(k)%nelems))
        self%dense_B_mats(i,dim)%M=Btmp(:,:,dim)
      END IF
    END DO
    DEALLOCATE(Btmp)
  END IF
  dense_sizes(i) = 3*self%levels(level)%blocks(j)%np*self%levels(level)%blocks(k)%nelems
  compressed_size = compressed_size + dense_sizes(i)
  !$omp critical
  nblocks_complete = nblocks_complete + 1
  DO k=1,9
    IF(nblocks_complete==nblocks_progress(k))WRITE(*,'(A,2X,I2,A1)')oft_indent,k*10,'%'
  END DO
  !$omp end critical
END DO
!$omp single
IF(self%B_aca_tol>0.d0)THEN
  WRITE(*,'(2A)')oft_indent,'Building off-diagonal blocks using ACA+ and SVD compression'
ELSE
  WRITE(*,'(2A)')oft_indent,'Building off-diagonal blocks with SVD compression'
END IF
nblocks_progress = INT(self%nsparse*[0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0],4)
nblocks_complete = 0
!$omp end single
!$omp do schedule(dynamic,1)
DO i=1,self%nsparse
  level = self%sparse_blocks(1,i)
  j = self%sparse_blocks(2,i)
  k = self%sparse_blocks(3,i)
  sparse_sizes(i) = 0
  DO flip=1,2
    ii=(i-1)*2+flip
    IF(flip==2)THEN
      k = self%sparse_blocks(2,i)
      j = self%sparse_blocks(3,i)
    END IF
    IF(ALL([(ASSOCIATED(self%aca_B_dense(ii,dim)%M).OR.ASSOCIATED(self%aca_BU_mats(ii,dim)%M),dim=1,3)]))THEN
      !---Loaded from file
      DO dim=1,3
        IF(ASSOCIATED(self%aca_B_dense(ii,dim)%M))THEN
          size_out=self%levels(level)%blocks(k)%nelems*self%levels(level)%blocks(j)%nelems
        ELSE
          size_out=SIZE(self%aca_BV_mats(ii,dim)%M,DIM=1)*(self%levels(level)%blocks(k)%nelems+self%levels(level)%blocks(j)%nelems)
        END IF
        sparse_sizes(i) = sparse_sizes(i) + size_out
      END DO
    ELSE
      mat_updated=.TRUE.
      DO dim=1,3
        IF(ASSOCIATED(self%aca_B_dense(ii,dim)%M))DEALLOCATE(self%aca_B_dense(ii,dim)%M)
        IF(ASSOCIATED(self%aca_BU_mats(ii,dim)%M))DEALLOCATE(self%aca_BU_mats(ii,dim)%M,self%aca_BV_mats(ii,dim)%M)
      END DO
      !---All three directions with one ACA+ (shared pivots), then per-direction SVD recompression
      size_out=-1
      IF((self%B_aca_tol>0.d0).AND.ABS(self%levels(level)%mat_mask(j,k))==2)THEN
        CALL aca_approx3(ii,level,j,k,self%B_aca_tol,size_out)
        IF(size_out==-1)WRITE(*,'(2A,2I4)')oft_indent,'ACA+ failure',self%levels(level)%blocks(j)%np, &
          self%levels(level)%blocks(k)%nelems
      END IF
      IF(size_out>0)THEN
        DO dim=1,3
          CALL tw_lowrank_recompress(self%aca_BU_mats(ii,dim),self%aca_BV_mats(ii,dim),self%B_svd_tol,self%min_rank, &
            self%min_rank,size_out)
          sparse_sizes(i) = sparse_sizes(i) + size_out
        END DO
      ELSE
        !---ACA+ failed or not used: dense block (all three directions in one pass) and SVD compression
        ALLOCATE(Btmp(self%levels(level)%blocks(k)%np,self%levels(level)%blocks(j)%nelems,3))
        CALL tw_compute_Bops_block3(self%tw_obj,Btmp, &
          self%levels(level)%blocks(j),self%levels(level)%blocks(k),self%tw_obj%B_dx)
        DO dim=1,3
          ALLOCATE(self%aca_B_dense(ii,dim)%M(self%levels(level)%blocks(k)%np,self%levels(level)%blocks(j)%nelems))
          self%aca_B_dense(ii,dim)%M=Btmp(:,:,dim)
          CALL compress_block(ii,j,k,dim,self%B_svd_tol,size_out)
          sparse_sizes(i) = sparse_sizes(i) + size_out
        END DO
      END IF
    END IF
    IF(ALLOCATED(Btmp))DEALLOCATE(Btmp)
  END DO
  compressed_size = compressed_size + sparse_sizes(i)
  !$omp critical
  nblocks_complete = nblocks_complete + 1
  DO k=1,9
    IF(nblocks_complete==nblocks_progress(k))WRITE(*,'(A,2X,I2,A1)')oft_indent,k*10,'%'
  END DO
  !$omp end critical
END DO
!$omp end parallel
elapsed_time=mytimer%tock()
WRITE(*,'(2A,F6.1,A,ES9.2,A,ES9.2,A)')oft_indent,'Compressed size:', &
  compressed_size*1.d2/(3*REAL(self%tw_obj%np_active,8)*REAL(self%tw_obj%mesh%np,8)), &
  "%  (",REAL(compressed_size,8),"/",3*REAL(self%tw_obj%np_active,8)*REAL(self%tw_obj%mesh%np,8),")"
WRITE(*,'(3A)')oft_indent,'Time = ',time_to_string(elapsed_time)
!---Building thread blocking
ALLOCATE(self%B_thread_ptr(2,oft_env%nthreads+1))
self%B_thread_ptr=0
self%B_thread_ptr(:,1)=1
target_size=SUM(sparse_sizes)/REAL(oft_env%nthreads,8)
curr_size=0.d0
i=1
DO j=1,self%nsparse
  curr_size=curr_size+sparse_sizes(j)
  IF(curr_size>REAL(i,8)*target_size)THEN
    self%B_thread_ptr(1,i+1)=j+1
    i=i+1
  END IF
END DO
target_size=SUM(dense_sizes)/REAL(oft_env%nthreads,8)
curr_size=0.d0
i=1
DO j=1,self%ndense
  curr_size=curr_size+dense_sizes(j)
  IF(curr_size>REAL(i,8)*target_size)THEN
    self%B_thread_ptr(2,i+1)=j+1
    i=i+1
  END IF
END DO
DO j=1,oft_env%nthreads
  IF(self%B_thread_ptr(1,j+1)==0)self%B_thread_ptr(1,j+1)=self%nsparse+1
  IF(self%B_thread_ptr(2,j+1)==0)self%B_thread_ptr(2,j+1)=self%ndense+1
END DO
self%B_thread_ptr(:,oft_env%nthreads+1)=[self%nsparse,self%ndense]+1
DEALLOCATE(sparse_sizes,dense_sizes)
!
IF(mat_updated.AND.PRESENT(save_file))CALL save_to_file()
CALL oft_decrease_indent
CONTAINS
SUBROUTINE save_to_file()
INTEGER(4) :: ierr,io_unit,hash_tmp(6),file_counts(6)
CHARACTER(LEN=1) :: dim_id
CHARACTER(LEN=5) :: matrix_id
TYPE(HODLR_hash), TARGET :: sparse_hash
!---Save computed matrix to file
IF(TRIM(save_file)/='none')THEN
  hash_tmp(1) = tw_nelems_exp(self%tw_obj)
  hash_tmp(2) = self%tw_obj%mesh%nc
  hash_tmp(3) = self%ndense
  sparse_hash%nsparse = self%nsparse
  sparse_hash%aca_tol = REAL(self%B_aca_tol,c_float)
  sparse_hash%svd_tol = REAL(self%B_svd_tol,c_float)
  hash_tmp(4) = oft_simple_hash(C_LOC(sparse_hash),INT(c_sizeof(sparse_hash),8),starting_value=tw_cache_ver)
  DO i=1,self%nlevels
    hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(i)%mat_mask),INT(4*self%levels(i)%nblocks*self%levels(i)%nblocks,8),starting_value=hash_tmp(4))
  END DO
  hash_tmp(5) = oft_simple_hash(C_LOC(self%tw_obj%mesh%lc),INT(4*3*self%tw_obj%mesh%nc,8))
  hash_tmp(6) = oft_simple_hash(C_LOC(self%tw_obj%mesh%r),INT(8*3*self%tw_obj%mesh%np,8))
  WRITE(*,'(3A)')oft_indent,'Saving HODLR matrix from file: ',TRIM(save_file)
  CALL hdf5_create_file(TRIM(save_file))
  CALL hdf5_write(hash_tmp,TRIM(save_file),'MODEL_hash')
  !
  hash_tmp(6)=0
  DO i=1,self%ndense
    level = self%dense_blocks(1,i)
    j = self%dense_blocks(2,i)
    k = self%dense_blocks(3,i)
    hash_tmp(1:3) = [level,j,k]
    hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ipts),INT(4*self%levels(level)%blocks(j)%np,8))
    hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ielem),INT(4*self%levels(level)%blocks(k)%nelems,8))
    WRITE(matrix_id,'(I5)')i
    CALL hdf5_write(hash_tmp,TRIM(save_file),'DENSE_hash_'//TRIM(ADJUSTL(matrix_id)))
    DO dim=1,3
      WRITE(dim_id,'(I1)')dim
      CALL hdf5_write(self%dense_B_mats(i,dim)%M,TRIM(save_file),'DENSE_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id)
    END DO
  END DO
  !
  DO i=1,self%nsparse
    level = self%sparse_blocks(1,i)
    j = self%sparse_blocks(2,i)
    k = self%sparse_blocks(3,i)
    DO flip=1,2
      ii=(i-1)*2+flip
      IF(flip==2)THEN
        k = self%sparse_blocks(2,i)
        j = self%sparse_blocks(3,i)
      END IF
      WRITE(matrix_id,'(I5)')ii
      DO dim=1,3
        hash_tmp(1:3) = [level,j,k]
        hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ielem),INT(4*self%levels(level)%blocks(j)%nelems,8))
        hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ipts),INT(4*self%levels(level)%blocks(k)%np,8))
        WRITE(dim_id,'(I1)')dim
        IF(ASSOCIATED(self%aca_B_dense(ii,dim)%M))THEN
          hash_tmp(6)=-1
          CALL hdf5_write(hash_tmp,TRIM(save_file),'SPARSE_hash_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id)
          CALL hdf5_write(self%aca_B_dense(ii,dim)%M,TRIM(save_file),'SPARSE_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id)
        ELSE
          hash_tmp(6)=SIZE(self%aca_BV_mats(ii,dim)%M,DIM=1)
          CALL hdf5_write(hash_tmp,TRIM(save_file),'SPARSE_hash_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id)
          CALL hdf5_write(self%aca_BU_mats(ii,dim)%M,TRIM(save_file),'SPARSE_U_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id)
          CALL hdf5_write(self%aca_BV_mats(ii,dim)%M,TRIM(save_file),'SPARSE_V_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id)
        END IF
      END DO
    END DO
  END DO
  ! TODO: Save create appropriate hash to enable saving coil matrices
  ! WRITE(io_unit)hash_tmp
  ! WRITE(io_unit)self%hole_Vcoil_mat%M
  ! CLOSE(io_unit)
END IF
END SUBROUTINE save_to_file
SUBROUTINE read_from_file()
INTEGER(4) :: ierr,io_unit,hash_tmp(6),file_counts(6)
LOGICAL :: exists
CHARACTER(LEN=1) :: dim_id
CHARACTER(LEN=5) :: matrix_id
TYPE(HODLR_hash), TARGET :: sparse_hash
!---Try to load from file
IF(TRIM(save_file)/='none')THEN
  INQUIRE(FILE=TRIM(save_file),EXIST=exists)
  IF(exists)THEN
    hash_tmp(1) = tw_nelems_exp(self%tw_obj)
    hash_tmp(2) = self%tw_obj%mesh%nc
    hash_tmp(3) = self%ndense
    sparse_hash%nsparse = self%nsparse
    sparse_hash%aca_tol = REAL(self%B_aca_tol,c_float)
    sparse_hash%svd_tol = REAL(self%B_svd_tol,c_float)
    hash_tmp(4) = oft_simple_hash(C_LOC(sparse_hash),INT(c_sizeof(sparse_hash),8),starting_value=tw_cache_ver)
    DO i=1,self%nlevels
      hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(i)%mat_mask),INT(4*self%levels(i)%nblocks*self%levels(i)%nblocks,8),starting_value=hash_tmp(4))
    END DO
    hash_tmp(5) = oft_simple_hash(C_LOC(self%tw_obj%mesh%lc),INT(4*3*self%tw_obj%mesh%nc,8))
    hash_tmp(6) = oft_simple_hash(C_LOC(self%tw_obj%mesh%r),INT(8*3*self%tw_obj%mesh%np,8))
    WRITE(*,'(3A)')oft_indent,'Reading HODLR matrix from file: ',TRIM(save_file)
    CALL oft_increase_indent
    CALL hdf5_read(file_counts,TRIM(save_file),'MODEL_hash',success=exists)
    ! OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
    ! READ(io_unit, IOSTAT=ierr)file_counts
    IF(exists.AND.ALL(file_counts==hash_tmp))THEN
      hash_tmp(6)=0
      DO i=1,self%ndense
        level = self%dense_blocks(1,i)
        j = self%dense_blocks(2,i)
        k = self%dense_blocks(3,i)
        hash_tmp(1:3) = [level,j,k]
        hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ipts),INT(4*self%levels(level)%blocks(j)%np,8))
        hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ielem),INT(4*self%levels(level)%blocks(k)%nelems,8))
        WRITE(matrix_id,'(I5)')i
        CALL hdf5_read(file_counts,TRIM(save_file),'DENSE_hash_'//TRIM(ADJUSTL(matrix_id)),success=exists)
        ! READ(io_unit, IOSTAT=ierr)file_counts
        IF(exists.AND.ALL(file_counts==hash_tmp))THEN
          DO dim=1,3
            ALLOCATE(self%dense_B_mats(i,dim)%M(self%levels(level)%blocks(j)%np,self%levels(level)%blocks(k)%nelems))
            WRITE(dim_id,'(I1)')dim
            CALL hdf5_read(self%dense_B_mats(i,dim)%M,TRIM(save_file),'DENSE_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id,success=exists)
            ! READ(io_unit, IOSTAT=ierr)self%dense_B_mats(i,dim)%M
            IF(.NOT.exists)DEALLOCATE(self%dense_B_mats(i,dim)%M)
          END DO
        END IF
      END DO
      !
      DO i=1,self%nsparse
        level = self%sparse_blocks(1,i)
        j = self%sparse_blocks(2,i)
        k = self%sparse_blocks(3,i)
        DO flip=1,2
          ii=(i-1)*2+flip
          IF(flip==2)THEN
            k = self%sparse_blocks(2,i)
            j = self%sparse_blocks(3,i)
          END IF
          WRITE(matrix_id,'(I5)')ii
          DO dim=1,3
            WRITE(dim_id,'(I1)')dim
            hash_tmp(1:3) = [level,j,k]
            hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ielem),INT(4*self%levels(level)%blocks(j)%nelems,8))
            hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ipts),INT(4*self%levels(level)%blocks(k)%np,8))
            CALL hdf5_read(file_counts,TRIM(save_file),'SPARSE_hash_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id,success=exists)
            ! READ(io_unit, IOSTAT=ierr)file_counts
            IF(exists.AND.ALL(file_counts(1:5)==hash_tmp(1:5)))THEN
              IF(file_counts(6)<0)THEN
                ALLOCATE(self%aca_B_dense(ii,dim)%M(self%levels(level)%blocks(k)%np,self%levels(level)%blocks(j)%nelems))
                CALL hdf5_read(self%aca_B_dense(ii,dim)%M,TRIM(save_file),'SPARSE_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id,success=exists)
                ! READ(io_unit, IOSTAT=ierr)self%aca_B_dense(ii,dim)%M
                IF(.NOT.exists)DEALLOCATE(self%aca_B_dense(ii,dim)%M)
              ELSE
                ALLOCATE(self%aca_BU_mats(ii,dim)%M(self%levels(level)%blocks(k)%np,file_counts(6)))
                CALL hdf5_read(self%aca_BU_mats(ii,dim)%M,TRIM(save_file),'SPARSE_U_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id,success=exists)
                ! READ(io_unit, IOSTAT=ierr)self%aca_BU_mats(ii,dim)%M
                IF(.NOT.exists)THEN
                  DEALLOCATE(self%aca_BU_mats(ii,dim)%M)
                ELSE
                  ALLOCATE(self%aca_BV_mats(ii,dim)%M(file_counts(6),self%levels(level)%blocks(j)%nelems))
                  CALL hdf5_read(self%aca_BV_mats(ii,dim)%M,TRIM(save_file),'SPARSE_V_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id,success=exists)
                  ! READ(io_unit, IOSTAT=ierr)self%aca_BV_mats(ii,dim)%M
                  IF(.NOT.exists)DEALLOCATE(self%aca_BU_mats(ii,dim)%M,self%aca_BV_mats(ii,dim)%M)
                END IF
              END IF
            END IF
          END DO
        END DO
      END DO
      ! CLOSE(io_unit)
    ELSE
      WRITE(*,'(2A)')oft_indent,'Ignoring stored matrix: Model hashes do not match'
      IF(oft_debug_print(1))THEN
        WRITE(*,*)'File: ',file_counts
        WRITE(*,*)'Model:',hash_tmp
      END IF
      ! CLOSE(io_unit)
    END IF
    CALL oft_decrease_indent
  END IF
END IF
END SUBROUTINE read_from_file
!
SUBROUTINE compress_block(isparse,iblock,jblock,dim,tol,size_out)
INTEGER(4), INTENT(in) :: isparse,iblock,jblock,dim
REAL(8), INTENT(in) :: tol
INTEGER(8), INTENT(out) :: size_out
INTEGER(4) :: M,N,rank
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: U,V
!
M = SIZE(self%aca_B_dense(isparse,dim)%M,DIM=1)
N = SIZE(self%aca_B_dense(isparse,dim)%M,DIM=2)
CALL tw_lowrank_qrcp(self%aca_B_dense(isparse,dim)%M,tol,self%min_rank,U,V,rank)
IF(rank>0)THEN
  ALLOCATE(self%aca_BU_mats(isparse,dim)%M(M,rank),self%aca_BV_mats(isparse,dim)%M(rank,N))
  self%aca_BU_mats(isparse,dim)%M=U
  self%aca_BV_mats(isparse,dim)%M=V
  DEALLOCATE(self%aca_B_dense(isparse,dim)%M)
  size_out=INT(rank,8)*INT(M+N,8)
ELSE
  size_out=INT(M,8)*INT(N,8)
END IF
END SUBROUTINE compress_block
!
!---ACA+ approximation of all three directions of a B-field block at once
!   Runs ACA+ (tw_aca_plus) on the stacked block [Bx; By; Bz] (block points, one set per direction) with shared
!   pivots. One element at all points gives all three directions in one kernel evaluation, and the three directions
!   at a point are computed together and kept (tw_aca_B3block_col). The stacked factors are split per direction
!   (shared right factor) for recompression with tw_lowrank_recompress.
SUBROUTINE aca_approx3(isparse,ilevel,iblock,jblock,tol,size_out)
INTEGER(4), INTENT(in) :: isparse,ilevel,iblock,jblock
REAL(8), INTENT(in) :: tol
INTEGER(8), INTENT(out) :: size_out
INTEGER(4) :: M,M3,N,MIN_DIM,max_iter,Iref,Jref,nterms,d
INTEGER(8) :: full_size
REAL(8) :: tmp
REAL(8), ALLOCATABLE :: U(:,:),V(:,:)
TYPE(tw_aca_B3block) :: src
!
size_out=-2
M = self%levels(ilevel)%blocks(jblock)%np
M3 = 3*M
N = self%levels(ilevel)%blocks(iblock)%nelems
MIN_DIM = min(M,N)
IF(MIN_DIM<=self%aca_min_its)RETURN
full_size = INT(M3,8)*INT(N,8)
!---Stop if the stacked factors would not be smaller than the dense stacked block, or if a direction's factor
!   (M x nterms) would have more columns than rows
max_iter = MIN(3*self%leaf_target*2,INT(full_size/INT(M3+N,8),4)-1,MIN_DIM-1)
size_out=-1
src%op=>self
src%ilevel=ilevel
src%iblock=iblock
src%jblock=jblock
!---Random reference element (row) and stacked point column
CALL random_number(tmp)
Iref = INT(tmp*(N-1),4)+1
CALL random_number(tmp)
Jref = INT(tmp*(M3-1),4)+1
CALL tw_aca_plus(src,N,M3,tol,self%aca_min_its,max_iter,Iref,Jref,U,V,nterms)
IF(nterms>0)THEN
  !---Split per direction: B_d ~ BU_d BV, with BU_d = V(:,direction d)^T and the shared BV = U^T
  DO d=1,3
    ALLOCATE(self%aca_BU_mats(isparse,d)%M(M,nterms))
    self%aca_BU_mats(isparse,d)%M=TRANSPOSE(V(:,(d-1)*M+1:d*M))
    ALLOCATE(self%aca_BV_mats(isparse,d)%M(nterms,N))
    self%aca_BV_mats(isparse,d)%M=TRANSPOSE(U)
  END DO
  size_out=INT(nterms,8)*INT(M3+N,8)
  DEALLOCATE(U,V)
END IF
IF(ASSOCIATED(src%one_elem%ielem))DEALLOCATE(src%one_elem%ielem,src%one_elem%icell,src%one_elem%inv_map)
IF(ASSOCIATED(src%one_pt%ipts))DEALLOCATE(src%one_pt%ipts)
END SUBROUTINE aca_approx3
END SUBROUTINE tw_hodlr_Bcompute
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
subroutine tw_hodlr_Lapply(self,a,b)
class(oft_tw_hodlr_op), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a
class(oft_vector), intent(inout) :: b
integer(4) :: i,j,k,iblock,jblock,level
real(8), allocatable :: row_tmp(:),col_tmp(:),int_tmp1(:),int_tmp2(:)
REAL(8), POINTER :: avals(:),bvals(:),ared(:),bred(:)
NULLIFY(avals,bvals,ared,bred)
CALL a%get_local(ared)
CALL b%set(0.d0)
CALL b%get_local(bred)
IF(ASSOCIATED(self%tw_obj%junc))THEN
  ALLOCATE(avals(tw_nelems_exp(self%tw_obj)),bvals(tw_nelems_exp(self%tw_obj)))
  CALL tw_junction_expand(self%tw_obj,ared,avals)
  bvals=0.d0
  DEALLOCATE(ared)
ELSE
  avals=>ared
  bvals=>bred
END IF
!$omp parallel private(j,k,iblock,jblock,level,row_tmp,col_tmp,int_tmp1,int_tmp2)
ALLOCATE(col_tmp(self%tw_obj%np_active),row_tmp(self%tw_obj%np_active))
ALLOCATE(int_tmp1(self%tw_obj%np_active),int_tmp2(self%tw_obj%np_active))
int_tmp2=0.d0
! !$omp do schedule(dynamic,4)
DO j=self%L_thread_ptr(1,oft_tid+1),self%L_thread_ptr(1,oft_tid+2)-1
  level = self%sparse_blocks(1,j)
  iblock = self%sparse_blocks(2,j)
  jblock = self%sparse_blocks(3,j)
  IF(ASSOCIATED(self%aca_dense(j)%M))THEN
    !$omp simd
    DO k=1,self%levels(level)%blocks(iblock)%nelems
      col_tmp(k)=avals(self%levels(level)%blocks(iblock)%ielem(k))
    END DO
    CALL dgemv('N',self%levels(level)%blocks(jblock)%nelems,self%levels(level)%blocks(iblock)%nelems,1.d0,self%aca_dense(j)%M, &
      self%levels(level)%blocks(jblock)%nelems,col_tmp,1,0.d0,row_tmp,1)
    !$omp simd
    DO k=1,self%levels(level)%blocks(jblock)%nelems
      int_tmp2(self%levels(level)%blocks(jblock)%ielem(k))=int_tmp2(self%levels(level)%blocks(jblock)%ielem(k))  + row_tmp(k)
    END DO
    !
    !$omp simd
    DO k=1,self%levels(level)%blocks(jblock)%nelems
      col_tmp(k)=avals(self%levels(level)%blocks(jblock)%ielem(k))
    END DO
    CALL dgemv('T',self%levels(level)%blocks(jblock)%nelems,self%levels(level)%blocks(iblock)%nelems,1.d0,self%aca_dense(j)%M, &
      self%levels(level)%blocks(jblock)%nelems,col_tmp,1,0.d0,row_tmp,1)
    !$omp simd
    DO k=1,self%levels(level)%blocks(iblock)%nelems
      int_tmp2(self%levels(level)%blocks(iblock)%ielem(k))=int_tmp2(self%levels(level)%blocks(iblock)%ielem(k)) + row_tmp(k)
    END DO
  ELSE
    !---Upper side
    !$omp simd
    DO k=1,self%levels(level)%blocks(iblock)%nelems
      col_tmp(k)=avals(self%levels(level)%blocks(iblock)%ielem(k))
    END DO
    k=SIZE(self%aca_V_mats(j)%M,DIM=1)
    CALL dgemv('N',k,self%levels(level)%blocks(iblock)%nelems,1.d0,self%aca_V_mats(j)%M, &
      k,col_tmp,1,0.d0,int_tmp1,1)
    CALL dgemv('N',self%levels(level)%blocks(jblock)%nelems,k,1.d0,self%aca_U_mats(j)%M, &
      self%levels(level)%blocks(jblock)%nelems,int_tmp1,1,0.d0,row_tmp,1)
    !$omp simd
    DO k=1,self%levels(level)%blocks(jblock)%nelems
      int_tmp2(self%levels(level)%blocks(jblock)%ielem(k))=int_tmp2(self%levels(level)%blocks(jblock)%ielem(k))+row_tmp(k)
    END DO
    !---Apply other side of diagonal
    !$omp simd
    DO k=1,self%levels(level)%blocks(jblock)%nelems
      col_tmp(k)=avals(self%levels(level)%blocks(jblock)%ielem(k))
    END DO
    k=SIZE(self%aca_V_mats(j)%M,DIM=1)
    CALL dgemv('T',self%levels(level)%blocks(jblock)%nelems,k,1.d0,self%aca_U_mats(j)%M, &
      self%levels(level)%blocks(jblock)%nelems,col_tmp,1,0.d0,int_tmp1,1)
    CALL dgemv('T',k,self%levels(level)%blocks(iblock)%nelems,1.d0,self%aca_V_mats(j)%M, &
      k,int_tmp1,1,0.d0,row_tmp,1)
    !$omp simd
    DO k=1,self%levels(level)%blocks(iblock)%nelems
      int_tmp2(self%levels(level)%blocks(iblock)%ielem(k))=int_tmp2(self%levels(level)%blocks(iblock)%ielem(k))+row_tmp(k)
    END DO
  END IF
END DO
!!$omp end do nowait
! !$omp do schedule(dynamic,1)
DO j=self%L_thread_ptr(2,oft_tid+1),self%L_thread_ptr(2,oft_tid+2)-1
  level = self%dense_blocks(1,j)
  iblock = self%dense_blocks(2,j)
  jblock = self%dense_blocks(3,j)
  !$omp simd
  DO k=1,self%levels(level)%blocks(iblock)%nelems
    col_tmp(k)=avals(self%levels(level)%blocks(iblock)%ielem(k))
  END DO
  CALL dgemv('N',self%levels(level)%blocks(jblock)%nelems,self%levels(level)%blocks(iblock)%nelems,1.d0,self%dense_mats(j)%M, &
    self%levels(level)%blocks(jblock)%nelems,col_tmp,1,0.d0,row_tmp,1)
  !$omp simd
  DO k=1,self%levels(level)%blocks(jblock)%nelems
    int_tmp2(self%levels(level)%blocks(jblock)%ielem(k))=int_tmp2(self%levels(level)%blocks(jblock)%ielem(k))+row_tmp(k)
  END DO
END DO
!!$omp end do nowait
!$omp critical
bvals(1:self%tw_obj%np_active)=bvals(1:self%tw_obj%np_active)+int_tmp2
!$omp end critical
DEALLOCATE(col_tmp,row_tmp,int_tmp1,int_tmp2)
!$omp end parallel
IF(self%tw_obj%nholes+self%tw_obj%n_vcoils>0)THEN
  CALL dgemv('N',tw_nelems_exp(self%tw_obj),self%tw_obj%nholes+self%tw_obj%n_vcoils,1.d0,self%hole_Vcoil_mat%M, &
    tw_nelems_exp(self%tw_obj),avals(self%tw_obj%np_active+1:tw_nelems_exp(self%tw_obj)),1,1.d0,bvals,1)
  avals(self%tw_obj%np_active+1:tw_nelems_exp(self%tw_obj))=0.d0 ! Prevent double diagonal contributions
  CALL dgemv('T',tw_nelems_exp(self%tw_obj),self%tw_obj%nholes+self%tw_obj%n_vcoils,1.d0,self%hole_Vcoil_mat%M, &
    tw_nelems_exp(self%tw_obj),avals,1,1.d0,bvals(self%tw_obj%np_active+1:tw_nelems_exp(self%tw_obj)),1)
END IF
IF(ASSOCIATED(self%tw_obj%junc))THEN
  CALL tw_junction_reduce_vec(self%tw_obj,bvals,bred)
  DEALLOCATE(bvals)
  bvals=>bred
END IF
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
end subroutine tw_hodlr_Lapply
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
subroutine tw_hodlr_Bapply(self,a,bx,by,bz)
class(oft_tw_hodlr_op), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a
class(oft_vector), intent(inout) :: bx,by,bz
integer(4) :: i,j,jj,k,dim,iblock,jblock,level,flip
real(8), allocatable :: row_tmp(:),col_tmp(:),int_tmp1(:),int_tmp2(:,:)
REAL(8), POINTER :: avals(:),bvals(:,:),btmp(:),ared(:)
NULLIFY(avals,bvals,ared)
CALL a%get_local(ared)
IF(ASSOCIATED(self%tw_obj%junc))THEN
  ALLOCATE(avals(tw_nelems_exp(self%tw_obj)))
  CALL tw_junction_expand(self%tw_obj,ared,avals)
  DEALLOCATE(ared)
ELSE
  avals=>ared
END IF
CALL bx%set(0.d0)
CALL by%set(0.d0)
CALL bz%set(0.d0)
ALLOCATE(bvals(bx%n,3))
btmp=>bvals(:,1)
CALL bx%get_local(btmp)
btmp=>bvals(:,2)
CALL by%get_local(btmp)
btmp=>bvals(:,3)
CALL bz%get_local(btmp)
!$omp parallel private(j,k,dim,jj,flip,iblock,jblock,level,row_tmp,col_tmp,int_tmp1,int_tmp2)
ALLOCATE(col_tmp(self%tw_obj%np_active),row_tmp(self%tw_obj%mesh%np))
ALLOCATE(int_tmp1(self%tw_obj%mesh%np),int_tmp2(self%tw_obj%mesh%np,3))
int_tmp2=0.d0
!!$omp do schedule(dynamic,1)
DO j=self%B_thread_ptr(1,oft_tid+1),self%B_thread_ptr(1,oft_tid+2)-1
  level = self%sparse_blocks(1,j)
  iblock = self%sparse_blocks(2,j)
  jblock = self%sparse_blocks(3,j)
  DO flip=1,2
    jj=(j-1)*2+flip
    IF(flip==2)THEN
      jblock = self%sparse_blocks(2,j)
      iblock = self%sparse_blocks(3,j)
    END IF
    DO dim=1,3
      IF(ASSOCIATED(self%aca_B_dense(jj,dim)%M))THEN
        !$omp simd
        DO k=1,self%levels(level)%blocks(iblock)%nelems
          col_tmp(k)=avals(self%levels(level)%blocks(iblock)%ielem(k))
        END DO
        CALL dgemv('N',self%levels(level)%blocks(jblock)%np,self%levels(level)%blocks(iblock)%nelems,1.d0,self%aca_B_dense(jj,dim)%M, &
          self%levels(level)%blocks(jblock)%np,col_tmp,1,0.d0,row_tmp,1)
        !$omp simd
        DO k=1,self%levels(level)%blocks(jblock)%np
          int_tmp2(self%levels(level)%blocks(jblock)%ipts(k),dim)=int_tmp2(self%levels(level)%blocks(jblock)%ipts(k),dim)  + row_tmp(k)
        END DO
      ELSE
        !$omp simd
        DO k=1,self%levels(level)%blocks(iblock)%nelems
          col_tmp(k)=avals(self%levels(level)%blocks(iblock)%ielem(k))
        END DO
        k=SIZE(self%aca_BV_mats(jj,dim)%M,DIM=1)
        CALL dgemv('N',k,self%levels(level)%blocks(iblock)%nelems,1.d0,self%aca_BV_mats(jj,dim)%M, &
          k,col_tmp,1,0.d0,int_tmp1,1)
        CALL dgemv('N',self%levels(level)%blocks(jblock)%np,k,1.d0,self%aca_BU_mats(jj,dim)%M, &
          self%levels(level)%blocks(jblock)%np,int_tmp1,1,0.d0,row_tmp,1)
        !$omp simd
        DO k=1,self%levels(level)%blocks(jblock)%np
          int_tmp2(self%levels(level)%blocks(jblock)%ipts(k),dim)=int_tmp2(self%levels(level)%blocks(jblock)%ipts(k),dim)+row_tmp(k)
        END DO
      END IF
    END DO
  END DO
END DO
! !$omp end do nowait
!!$omp do schedule(static)
DO j=self%B_thread_ptr(2,oft_tid+1),self%B_thread_ptr(2,oft_tid+2)-1
  level = self%dense_blocks(1,j)
  iblock = self%dense_blocks(2,j)
  jblock = self%dense_blocks(3,j)
  !$omp simd
  DO k=1,self%levels(level)%blocks(iblock)%nelems
    col_tmp(k)=avals(self%levels(level)%blocks(iblock)%ielem(k))
  END DO
  DO dim=1,3
    CALL dgemv('N',self%levels(level)%blocks(jblock)%np,self%levels(level)%blocks(iblock)%nelems,1.d0,self%dense_B_mats(j,dim)%M, &
      self%levels(level)%blocks(jblock)%np,col_tmp,1,0.d0,row_tmp,1)
    !$omp simd
    DO k=1,self%levels(level)%blocks(jblock)%np
      int_tmp2(self%levels(level)%blocks(jblock)%ipts(k),dim)=int_tmp2(self%levels(level)%blocks(jblock)%ipts(k),dim)+row_tmp(k)
    END DO
  END DO
END DO
!$omp critical
DO dim=1,3
  bvals(1:self%tw_obj%mesh%np,dim)=bvals(1:self%tw_obj%mesh%np,dim)+int_tmp2(:,dim)
END DO
!$omp end critical
DEALLOCATE(col_tmp,row_tmp,int_tmp1,int_tmp2)
!$omp end parallel
IF((self%tw_obj%nholes+self%tw_obj%n_vcoils)>0)THEN
  DO dim=1,3
    CALL dgemv('N',self%tw_obj%mesh%np,self%tw_obj%nholes+self%tw_obj%n_vcoils,1.d0,self%hole_Vcoil_Bmat(:,:,dim), &
      self%tw_obj%mesh%np,avals(self%tw_obj%np_active+1:tw_nelems_exp(self%tw_obj)),1,1.d0,bvals(:,dim),1)
  END DO
END IF
CALL bx%restore_local(bvals(:,1))
CALL by%restore_local(bvals(:,2))
CALL bz%restore_local(bvals(:,3))
DEALLOCATE(avals,bvals)
end subroutine tw_hodlr_Bapply
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
subroutine tw_hodlr_Lassemble(self,diag)
class(oft_tw_hodlr_op), intent(inout) :: self
class(oft_vector), optional, target, intent(inout) :: diag
integer(4) :: i,j,level,iblock
REAL(8), POINTER, DIMENSION(:) :: vals,vred
IF(present(diag))THEN
  IF(associated(self%D))CALL self%D%delete
  CALL diag%new(self%D)
  NULLIFY(vals)
  CALL self%D%get_local(vals)
  IF(ASSOCIATED(self%tw_obj%junc))THEN
    vred=>vals
    ALLOCATE(vals(tw_nelems_exp(self%tw_obj)))
    vals=0.d0
  END IF
  !---Get hole diagonal values
  IF(self%tw_obj%nholes+self%tw_obj%n_vcoils>0)THEN
    DO i=1,self%tw_obj%nholes+self%tw_obj%n_vcoils
      vals(i+self%tw_obj%np_active)=self%hole_Vcoil_mat%M(i+self%tw_obj%np_active,i)
    END DO
  END IF
  DO i=1,self%ndense
    level=self%dense_blocks(1,i)
    iblock=self%dense_blocks(2,i)
    DO j=1,self%levels(level)%blocks(iblock)%nelems
      vals(self%levels(level)%blocks(iblock)%ielem(j))=self%dense_mats(i)%M(j,j)
    END DO
  END DO
  IF(ASSOCIATED(self%tw_obj%junc))THEN
    CALL tw_junction_reduce_diag(self%tw_obj,vals,vred)
    CALL tw_hodlr_junction_cross_diag(self,vred)
    DEALLOCATE(vals)
    vals=>vred
  END IF
  CALL self%D%restore_local(vals)
  DEALLOCATE(vals)
  CALL diag%add(0.d0,1.d0,self%D)
END IF
self%nr=self%tw_obj%nelems; self%nrg=self%tw_obj%nelems
self%nc=self%tw_obj%nelems; self%ncg=self%tw_obj%nelems
end subroutine tw_hodlr_Lassemble
!---------------------------------------------------------------------------------
!> Add the cross terms of the reduced diagonal of a junction model, sum over e/=f of P(e,a) L(e,f) P(f,a),
!! for the columns of P with more than one entry (entries of L are taken from the HODLR blocks)
!---------------------------------------------------------------------------------
SUBROUTINE tw_hodlr_junction_cross_diag(self,d)
CLASS(oft_tw_hodlr_op), INTENT(inout) :: self !< HODLR operator
REAL(8), INTENT(inout) :: d(:) !< Reduced diagonal [nelems]
TYPE :: level_index
  INTEGER(4), ALLOCATABLE :: ind(:,:)
END TYPE level_index
INTEGER(4) :: a,i,j,k,l,np
INTEGER(4), ALLOCATABLE :: pt_kr(:),pt_lc(:),blk(:,:),pos(:,:),dense_ind(:)
REAL(8), ALLOCATABLE :: pt_val(:)
TYPE(level_index), ALLOCATABLE :: sparse_ind(:)
CALL tw_junction_transpose(self%tw_obj,pt_kr,pt_lc,pt_val)
IF(ALL(pt_kr(2:)-pt_kr(:SIZE(pt_kr)-1)<2))THEN
  DEALLOCATE(pt_kr,pt_lc,pt_val)
  RETURN
END IF
np=self%tw_obj%np_active
!---Block and position of each vertex on each level, and stored interaction of each block pair
ALLOCATE(blk(self%nlevels,np),pos(self%nlevels,np),sparse_ind(self%nlevels))
DO l=1,self%nlevels
  DO i=1,self%levels(l)%nblocks
    DO j=1,self%levels(l)%blocks(i)%nelems
      blk(l,self%levels(l)%blocks(i)%ielem(j))=i
      pos(l,self%levels(l)%blocks(i)%ielem(j))=j
    END DO
  END DO
  ALLOCATE(sparse_ind(l)%ind(self%levels(l)%nblocks,self%levels(l)%nblocks))
  sparse_ind(l)%ind=0
END DO
DO i=1,self%nsparse
  sparse_ind(self%sparse_blocks(1,i))%ind(self%sparse_blocks(2,i),self%sparse_blocks(3,i))=i
END DO
ALLOCATE(dense_ind(self%levels(self%nlevels)%nblocks))
DO i=1,self%ndense
  dense_ind(self%dense_blocks(2,i))=i
END DO
!---Cross terms
DO a=1,SIZE(d)
  DO k=pt_kr(a),pt_kr(a+1)-1
    DO j=pt_kr(a),pt_kr(a+1)-1
      IF(j==k)CYCLE
      d(a)=d(a)+pt_val(k)*pt_val(j)*entry(pt_lc(k),pt_lc(j))
    END DO
  END DO
END DO
DO l=1,self%nlevels
  DEALLOCATE(sparse_ind(l)%ind)
END DO
DEALLOCATE(pt_kr,pt_lc,pt_val,blk,pos,dense_ind,sparse_ind)
CONTAINS
!> Entry (e,f) of the expanded inductance matrix
FUNCTION entry(e,f) RESULT(val)
INTEGER(4), INTENT(in) :: e !< Row
INTEGER(4), INTENT(in) :: f !< Column
REAL(8) :: val
INTEGER(4) :: lev,ib,jb,isp,r,c
IF(f>np)THEN
  val=self%hole_Vcoil_mat%M(e,f-np)
  RETURN
ELSE IF(e>np)THEN
  val=self%hole_Vcoil_mat%M(f,e-np)
  RETURN
END IF
DO lev=1,self%nlevels
  ib=blk(lev,e)
  jb=blk(lev,f)
  IF(ib==jb)CYCLE
  isp=MAX(sparse_ind(lev)%ind(ib,jb),sparse_ind(lev)%ind(jb,ib))
  IF(isp==0)CYCLE
  !---Stored as rows of the second block, columns of the first
  IF(blk(lev,e)==self%sparse_blocks(3,isp))THEN
    r=pos(lev,e); c=pos(lev,f)
  ELSE
    r=pos(lev,f); c=pos(lev,e)
  END IF
  IF(ASSOCIATED(self%aca_dense(isp)%M))THEN
    val=self%aca_dense(isp)%M(r,c)
  ELSE
    val=DOT_PRODUCT(self%aca_U_mats(isp)%M(r,:),self%aca_V_mats(isp)%M(:,c))
  END IF
  RETURN
END DO
val=self%dense_mats(dense_ind(blk(self%nlevels,e)))%M(pos(self%nlevels,e),pos(self%nlevels,f))
END FUNCTION entry
END SUBROUTINE tw_hodlr_junction_cross_diag
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE tw_part_mesh(self,leaf_target,nlevels,levels)
TYPE(tw_type), INTENT(inout) :: self
INTEGER(4), INTENT(in) :: leaf_target
INTEGER(4), INTENT(out) :: nlevels
TYPE(oft_tw_level), POINTER, INTENT(out) :: levels(:)
INTEGER(4) :: i,j
REAL(8) :: xs
! INTEGER(4), ALLOCATABLE, DIMENSION(:) :: child_inds
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: active_pts
!---------------------------------------------------------------------------------
!> Need docs
!---------------------------------------------------------------------------------
TYPE :: tw_oct_tree
  INTEGER(4) :: npts = 0
  INTEGER(4) :: nchildren = 0
  INTEGER(4) :: depth = -1
  INTEGER(4) :: dind = -1
  INTEGER(4) :: pind = -1
  REAL(8) :: bounds(2,3) = 0.d0
  REAL(8) :: extent = 1.d99
  REAL(8) :: center(3) = 0.d0
  INTEGER(4), POINTER, CONTIGUOUS, DIMENSION(:) :: pts => NULL()
  TYPE(tw_oct_tree), POINTER, DIMENSION(:) :: children => NULL()
  TYPE(tw_oct_tree), POINTER :: parent => NULL()
END TYPE tw_oct_tree

TYPE(tw_oct_tree) :: mesh_tree

!---Initialize tree
mesh_tree%npts=self%mesh%np !self%np_active
ALLOCATE(mesh_tree%pts(mesh_tree%npts))
ALLOCATE(active_pts(3,mesh_tree%npts))
DO i=1,self%mesh%np
  ! IF(self%pmap(i)==0)CYCLE
  ! mesh_tree%pts(self%pmap(i))=self%pmap(i)
  ! active_pts(:,self%pmap(i))=self%mesh%r(:,i)
  mesh_tree%pts(i)=i
  active_pts(:,i)=self%mesh%r(:,i)
END DO
!---find bounding box around mesh
DO i=1,3
  mesh_tree%bounds(1,i)=MINVAL(active_pts(i,:))
  mesh_tree%bounds(2,i)=MAXVAL(active_pts(i,:))
  IF(mesh_tree%bounds(2,i)-mesh_tree%bounds(1,i)<1.d-3)THEN
  mesh_tree%bounds(1,i)=mesh_tree%bounds(1,i)-5.d-4
  mesh_tree%bounds(2,i)=mesh_tree%bounds(2,i)+5.d-4
  END IF
  xs = mesh_tree%bounds(2,i) - mesh_tree%bounds(1,i)
  mesh_tree%bounds(1,i) = mesh_tree%bounds(1,i) - xs*1.d-2
  mesh_tree%bounds(2,i) = mesh_tree%bounds(2,i) + xs*1.d-2
END DO
! mesh_tree%extents=mesh_tree%bounds
!mesh_tree%center=SUM(mesh_tree%bounds,DIM=1)/2.d0
mesh_tree%center=SUM(active_pts,DIM=2)/REAL(mesh_tree%npts,8)
mesh_tree%depth=1
!---Recursively subdivide
nlevels=1
IF(mesh_tree%npts>leaf_target)THEN
  CALL subdivide_leaf(mesh_tree,2)
ELSE
  nlevels=1
END IF
DEALLOCATE(active_pts)
ALLOCATE(levels(nlevels))
CALL count_levels(mesh_tree,1)
DO i=1,nlevels
  ALLOCATE(levels(i)%blocks(levels(i)%nblocks))
  levels(i)%nblocks=0
END DO
CALL fill_levels(mesh_tree,1,1)
CONTAINS
!
RECURSIVE SUBROUTINE count_levels(leaf,depth)
TYPE(tw_oct_tree), INTENT(inout) :: leaf
INTEGER(4), INTENT(in) :: depth
INTEGER(4) :: j
!---Ensure we are at the bottom of each tree
IF(leaf%nchildren>0)THEN
  DO j=1,leaf%nchildren
    CALL count_levels(leaf%children(j),depth+1)
  END DO
ELSE
  DO j=depth+1,nlevels
    levels(j)%nblocks=levels(j)%nblocks+1
  END DO
END IF
!---Add leaf if non-empty
IF(leaf%npts>0)levels(depth)%nblocks=levels(depth)%nblocks+1
END SUBROUTINE count_levels
!
RECURSIVE SUBROUTINE fill_levels(leaf,depth,parent)
TYPE(tw_oct_tree), INTENT(inout) :: leaf
INTEGER(4), INTENT(in) :: depth,parent
INTEGER(4) :: i,j,block_id
!---Add leaf if non-empty
IF(leaf%npts==0)RETURN
levels(depth)%nblocks=levels(depth)%nblocks+1
!
block_id=levels(depth)%nblocks
levels(depth)%blocks(block_id)%np=leaf%npts
ALLOCATE(levels(depth)%blocks(block_id)%ipts(leaf%npts))
levels(depth)%blocks(block_id)%ipts=leaf%pts
levels(depth)%blocks(block_id)%nelems=0
DO i=1,leaf%npts
  IF(self%pmap(leaf%pts(i))==0)CYCLE
  levels(depth)%blocks(block_id)%nelems=levels(depth)%blocks(block_id)%nelems+1
END DO
ALLOCATE(levels(depth)%blocks(block_id)%ielem(levels(depth)%blocks(block_id)%nelems))
levels(depth)%blocks(block_id)%nelems=0
DO i=1,leaf%npts
  IF(self%pmap(leaf%pts(i))==0)CYCLE
  levels(depth)%blocks(block_id)%nelems=levels(depth)%blocks(block_id)%nelems+1
  levels(depth)%blocks(block_id)%ielem(levels(depth)%blocks(block_id)%nelems)=self%pmap(leaf%pts(i))
END DO
! levels(depth)%blocks(block_id)%ielem=leaf%pts
levels(depth)%blocks(block_id)%center=leaf%center
levels(depth)%blocks(block_id)%extent=leaf%extent
levels(depth)%blocks(block_id)%parent=parent
!---Ensure we are at the bottom of each tree
IF(leaf%nchildren>0)THEN
  DO j=1,leaf%nchildren
    CALL fill_levels(leaf%children(j),depth+1,block_id)
  END DO
  DEALLOCATE(leaf%children)
  NULLIFY(leaf%parent)
ELSE
  DO j=depth+1,nlevels
    levels(j)%nblocks=levels(j)%nblocks+1
    !
    block_id=levels(j)%nblocks
    levels(j)%blocks(block_id)%np=leaf%npts
    ALLOCATE(levels(j)%blocks(block_id)%ipts(leaf%npts))
    levels(j)%blocks(block_id)%ipts=leaf%pts
    levels(j)%blocks(block_id)%nelems=0
    DO i=1,leaf%npts
      IF(self%pmap(leaf%pts(i))==0)CYCLE
      levels(j)%blocks(block_id)%nelems=levels(j)%blocks(block_id)%nelems+1
    END DO
    ALLOCATE(levels(j)%blocks(block_id)%ielem(levels(j)%blocks(block_id)%nelems))
    levels(j)%blocks(block_id)%nelems=0
    DO i=1,leaf%npts
      IF(self%pmap(leaf%pts(i))==0)CYCLE
      levels(j)%blocks(block_id)%nelems=levels(j)%blocks(block_id)%nelems+1
      levels(j)%blocks(block_id)%ielem(levels(j)%blocks(block_id)%nelems)=self%pmap(leaf%pts(i))
    END DO
    ! ALLOCATE(levels(j)%blocks(block_id)%ielem(leaf%npts))
    ! levels(j)%blocks(block_id)%ielem=leaf%pts
    levels(j)%blocks(block_id)%center=leaf%center
    levels(j)%blocks(block_id)%extent=leaf%extent
    levels(j)%blocks(block_id)%parent=levels(j-1)%nblocks
  END DO
END IF
!---Destroy leaf on way back up
leaf%npts=0
DEALLOCATE(leaf%pts)
END SUBROUTINE fill_levels
!---Helper function
RECURSIVE SUBROUTINE subdivide_leaf(leaf,depth)
TYPE(tw_oct_tree), TARGET, INTENT(inout) :: leaf
INTEGER(4), INTENT(in) :: depth
INTEGER(4) :: j,k,kk,i,ii,npts,sep_max(1)
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: child_mark
REAL(8) :: xm,ym,zm,rleft(3),rdenom(3),locs(3),xs,extents(2,3),split_center(3)
REAL(8), ALLOCATABLE, DIMENSION(:) :: vals
REAL(8), PARAMETER :: pt_tol=0.d-6
TYPE(tw_oct_tree) :: child_tmp(2)
! WRITE(*,*)'Subdividing block',leaf%ncells,nblocks
IF(nlevels>20)CALL oft_abort("Maximum recursion", "subdivide_leaf", __FILE__)
nlevels=MAX(nlevels,depth)
ALLOCATE(vals(leaf%npts),child_mark(leaf%npts))
child_mark=0
locs=0.d0
DO j=1,3
  vals=active_pts(j,leaf%pts)
  CALL sort_array(vals,child_mark,leaf%npts)
  split_center(j)=vals(INT(leaf%npts/2,4))
  locs(j)=SUM((vals-split_center(j))**2)
END DO
DEALLOCATE(vals,child_mark)
!---Split cell in half in each direction
xm = split_center(1) !(leaf%bounds(1,1)+leaf%bounds(2,1))/2.d0
ym = split_center(2) !(leaf%bounds(1,2)+leaf%bounds(2,2))/2.d0
zm = split_center(3) !(leaf%bounds(1,3)+leaf%bounds(2,3))/2.d0
DO k=1,leaf%npts
  IF(ABS(active_pts(1,leaf%pts(k))-xm)<1.d-8)xm=xm+(leaf%bounds(2,1)-leaf%bounds(1,1))*1.d-4
  IF(ABS(active_pts(2,leaf%pts(k))-ym)<1.d-8)ym=ym+(leaf%bounds(2,2)-leaf%bounds(1,2))*1.d-4
  IF(ABS(active_pts(3,leaf%pts(k))-zm)<1.d-8)zm=zm+(leaf%bounds(2,3)-leaf%bounds(1,3))*1.d-4
END DO
!---Split cell in half
child_tmp(1)%bounds=leaf%bounds
child_tmp(2)%bounds=leaf%bounds
sep_max=MAXLOC(locs)
j=sep_max(1)
!WRITE(*,*)'Split',j,locs
locs=[xm,ym,zm]
child_tmp(1)%bounds(1,j)=locs(j)
child_tmp(2)%bounds(2,j)=locs(j)
! WRITE(*,*)0
! WRITE(*,*)leaf%bounds(:,1)
! WRITE(*,*)leaf%bounds(:,2)
! WRITE(*,*)leaf%bounds(:,3)
! DO j=1,8
!   WRITE(*,*)j
!   WRITE(*,*)child_tmp(j)%bounds(:,1)
!   WRITE(*,*)child_tmp(j)%bounds(:,2)
!   WRITE(*,*)child_tmp(j)%bounds(:,3)
! END DO
!---
ALLOCATE(child_mark(leaf%npts))
child_mark=0
leaf%nchildren=0
DO j=1,2
  child_tmp(j)%depth=depth
  child_tmp(j)%parent=>leaf
  ! child_mark=0
  rleft = child_tmp(j)%bounds(1,:)
  rdenom = child_tmp(j)%bounds(2,:)-child_tmp(j)%bounds(1,:)
  DO k=1,leaf%npts
    locs = (active_pts(:,leaf%pts(k))-rleft)/rdenom
    IF(ALL(locs>=-pt_tol).AND.ALL(locs<1.d0+pt_tol).AND.child_mark(k)==0)child_mark(k)=j
  END DO
  child_tmp(j)%npts=COUNT(child_mark==j)
  IF(child_tmp(j)%npts==0)CYCLE
  ! IF(child_tmp(j)%npts<200)WRITE(*,*)'Small leaf',child_tmp(j)%npts
  leaf%nchildren=leaf%nchildren+1
  ALLOCATE(child_tmp(j)%pts(child_tmp(j)%npts))
  extents(1,:)=1.d99
  extents(2,:)=-1.d99
  child_tmp(j)%npts=0
  child_tmp(j)%center=0.d0
  DO k=1,leaf%npts
    IF(child_mark(k)/=j)CYCLE
    child_tmp(j)%npts=child_tmp(j)%npts+1
    child_tmp(j)%pts(child_tmp(j)%npts)=leaf%pts(k)
    child_tmp(j)%center=child_tmp(j)%center+active_pts(:,leaf%pts(k))
    DO kk=1,3
      extents(1,kk)=MIN(extents(1,kk),active_pts(kk,leaf%pts(k)))
      extents(2,kk)=MAX(extents(2,kk),active_pts(kk,leaf%pts(k)))
    END DO
  END DO
  DO kk=1,3
    IF(extents(2,kk)-extents(1,kk)<1.d-8)THEN
      extents(1,kk)=extents(1,kk)-5.d-9
      extents(2,kk)=extents(2,kk)+5.d-9
    END IF
    xs = extents(2,kk) - extents(1,kk)
    extents(1,kk) = MAX(child_tmp(j)%bounds(1,kk),extents(1,kk) - xs*1.d-2)
    extents(2,kk) = MIN(child_tmp(j)%bounds(2,kk),extents(2,kk) + xs*1.d-2)
    !---Shrink bounds to match extent
    child_tmp(j)%bounds(1,kk) = MAX(child_tmp(j)%bounds(1,kk),extents(1,kk))
    child_tmp(j)%bounds(2,kk) = MIN(child_tmp(j)%bounds(2,kk),extents(2,kk))
  END DO
  child_tmp(j)%center=child_tmp(j)%center/REAL(child_tmp(j)%npts,8)
  child_tmp(j)%extent=0.d0
  DO k=1,child_tmp(j)%npts
    child_tmp(j)%extent=MAX(child_tmp(j)%extent,magnitude(active_pts(:,child_tmp(j)%pts(k))-child_tmp(j)%center))
  END DO
END DO
IF(ANY(child_mark==0))THEN
  WRITE(*,*)leaf%depth
  WRITE(*,*)leaf%bounds(1,:)
  WRITE(*,*)leaf%bounds(2,:)
  WRITE(*,*)"---"
  DO k=1,leaf%npts
    IF(child_mark(k)==0)THEN
      WRITE(*,*)active_pts(:,leaf%pts(k))
      DO j=1,2!8
        rleft = child_tmp(j)%bounds(1,:)
        rdenom = child_tmp(j)%bounds(2,:)-child_tmp(j)%bounds(1,:)
        locs = (active_pts(:,leaf%pts(k))-rleft)/rdenom
        WRITE(*,*)j,locs
      END DO
    END IF
  END DO
  CALL oft_abort("Lost element in partitioning","tw_part_mesh::subdivide_leaf",__FILE__)
END IF
DEALLOCATE(child_mark)
!---Subdivide further
IF(leaf%nchildren==0)RETURN
! WRITE(*,*)'Parent',depth,leaf%nchildren,leaf%ncells,branch_target==0
ALLOCATE(leaf%children(leaf%nchildren))
leaf%nchildren=0
DO j=1,2
  IF(child_tmp(j)%npts>0)THEN
    leaf%nchildren=leaf%nchildren+1
    CALL clone_tree(child_tmp(j),leaf%children(leaf%nchildren))
    ! WRITE(*,*)'  Child block',child_tmp(j)%ncells,child_tmp(j)%ncells<=branch_target
    ! DO kk=1,3
    !   WRITE(*,*)'    ',child_tmp(j)%bounds(:,kk)
    ! END DO
    IF(leaf%children(leaf%nchildren)%npts>leaf_target)THEN
      CALL subdivide_leaf(leaf%children(leaf%nchildren),depth+1)
    ! ELSE
    !   nleaves=nleaves+1
    END IF
  END IF
END DO
! WRITE(*,*)'Exit',depth
END SUBROUTINE subdivide_leaf
SUBROUTINE clone_tree(in_tree,out_tree)
TYPE(tw_oct_tree), INTENT(inout) :: in_tree,out_tree
out_tree%dind=in_tree%dind
out_tree%pind=in_tree%pind
out_tree%depth=in_tree%depth
out_tree%bounds=in_tree%bounds
out_tree%extent=in_tree%extent
out_tree%center=in_tree%center
out_tree%npts=in_tree%npts
out_tree%nchildren=in_tree%nchildren
out_tree%pts=>in_tree%pts
out_tree%children=>in_tree%children
out_tree%parent=>in_tree%parent
END SUBROUTINE clone_tree
END SUBROUTINE tw_part_mesh
!------------------------------------------------------------------------------
!> Low-rank compression of a dense block by truncated pivoted QR and a small SVD
!!
!! Householder QR with column pivoting is stopped once the norm of the remaining block \f$ \|R_{22}\|_F \f$ (known
!! exactly from the column norms) is below a tenth of `tol`. The small factor is then compressed with an SVD and
!! truncated with the same rule as a full SVD: the smallest rank whose discarded part, including \f$ R_{22} \f$, has
!! a Frobenius norm below `tol`. Costs \f$ O(M N r) \f$ instead of \f$ O(M N \min(M,N)) \f$ for a full SVD.
!!
!! Pivoted QR can need a few more columns than the SVD rank to reach the tolerance, so it may run past the
!! break-even rank (by `MAX(8,kbe/8)` columns); the block is compressed only if the final rank is below it.
!------------------------------------------------------------------------------
SUBROUTINE tw_lowrank_qrcp(A,tol,min_rank,U,V,rank)
REAL(8), INTENT(in) :: A(:,:) !< Dense block [M,N]
REAL(8), INTENT(in) :: tol !< Absolute tolerance on the Frobenius norm of the error
INTEGER(4), INTENT(in) :: min_rank !< Minimum rank to keep
REAL(8), ALLOCATABLE, INTENT(out) :: U(:,:) !< Left factor [M,rank]
REAL(8), ALLOCATABLE, INTENT(out) :: V(:,:) !< Right factor, singular values included [rank,N]
INTEGER(4), INTENT(out) :: rank !< Rank of the factors (-1 if the block is not worth compressing)
INTEGER(4) :: M,N,MIN_DIM,kmin,kbe,kqr,r,j,p,i,LWORK,INFO
INTEGER(4), ALLOCATABLE :: perm(:),IWORK(:)
REAL(8) :: aii,tmp,etrail,tstop
REAL(8), ALLOCATABLE :: W(:,:),cn(:),cn0(:),tau(:),work(:),Q(:,:),Rs(:,:),S(:),Ur(:,:),Vt(:,:)
M=SIZE(A,1); N=SIZE(A,2)
MIN_DIM=MIN(M,N)
kmin=MAX(1,MIN(min_rank,MIN_DIM))
kbe=INT(INT(M,8)*INT(N,8)/INT(M+N,8),4) ! Break-even rank: no savings at or above
kqr=MIN(MIN_DIM,kbe+MAX(8,kbe/8)) ! Most QR columns, past kbe since QR can need a few more columns than the SVD
rank=-1
ALLOCATE(W(M,N),perm(N),cn(N),cn0(N),tau(MIN_DIM),work(MAX(M,N)))
W=A
DO j=1,N
  perm(j)=j
  cn(j)=SUM(W(:,j)**2)
END DO
cn0=cn
tstop=(0.1d0*tol)**2
etrail=SUM(cn)
!---Pivoted Householder QR, stopped once the remaining block is small enough
r=0
DO WHILE(r<MIN_DIM)
  IF((etrail<=tstop).AND.(r>=kmin))EXIT
  IF(r>=kqr)EXIT
  r=r+1
  p=r-1+MAXLOC(cn(r:N),DIM=1)
  IF(p/=r)THEN
    work(1:M)=W(:,r); W(:,r)=W(:,p); W(:,p)=work(1:M)
    j=perm(r); perm(r)=perm(p); perm(p)=j
    tmp=cn(r); cn(r)=cn(p); cn(p)=tmp
    tmp=cn0(r); cn0(r)=cn0(p); cn0(p)=tmp
  END IF
  CALL DLARFG(M-r+1,W(r,r),W(MIN(r+1,M),r),1,tau(r))
  IF(r<N)THEN
    aii=W(r,r)
    W(r,r)=1.d0
    CALL DLARF('L',M-r+1,N-r,W(r,r),1,tau(r),W(r,r+1),M,work)
    W(r,r)=aii
  END IF
  !---Downdate the remaining column norms (recomputed where cancellation sets in)
  DO j=r+1,N
    IF(r>=M)THEN
      cn(j)=0.d0
    ELSE
      tmp=cn(j)-W(r,j)**2
      IF(tmp<=SQRT(EPSILON(1.d0))*cn0(j))THEN
        cn(j)=SUM(W(r+1:M,j)**2)
        cn0(j)=cn(j)
      ELSE
        cn(j)=tmp
      END IF
    END IF
  END DO
  etrail=0.d0
  IF(r<N)etrail=SUM(cn(r+1:N))
END DO
IF((etrail>tol**2).OR.(r==0))THEN ! Remaining block above tol after kqr columns: not compressible
  DEALLOCATE(W,perm,cn,cn0,tau,work)
  RETURN
END IF
!---Q factor and R factor (columns in the original order)
ALLOCATE(Q(M,r),Rs(r,N))
Q=W(:,1:r)
LWORK=-1
CALL DORGQR(M,r,r,Q,M,tau,work,LWORK,INFO)
LWORK=INT(work(1))
DEALLOCATE(work); ALLOCATE(work(LWORK))
CALL DORGQR(M,r,r,Q,M,tau,work,LWORK,INFO)
IF(INFO/=0)CALL oft_abort("Forming Q failed","tw_lowrank_qrcp",__FILE__)
Rs=0.d0
DO j=1,N
  DO i=1,MIN(j,r)
    Rs(i,perm(j))=W(i,j)
  END DO
END DO
DEALLOCATE(W)
!---SVD of the small factor and truncation, including the remaining block
ALLOCATE(S(r),Ur(r,r),Vt(r,N),IWORK(8*r))
LWORK=-1
CALL DGESDD('S',r,N,Rs,r,S,Ur,r,Vt,r,work,LWORK,IWORK,INFO)
LWORK=INT(work(1))
DEALLOCATE(work); ALLOCATE(work(LWORK))
CALL DGESDD('S',r,N,Rs,r,S,Ur,r,Vt,r,work,LWORK,IWORK,INFO)
IF(INFO/=0)CALL oft_abort("The algorithm computing SVD failed to converge.","tw_lowrank_qrcp",__FILE__)
tmp=etrail
DO i=r,1,-1
  tmp=tmp+S(i)**2
  IF(SQRT(tmp)>tol)EXIT
END DO
i=MIN(r,MAX(kmin,i))
IF(i<kbe)THEN
  rank=i
  ALLOCATE(U(M,i),V(i,N))
  CALL DGEMM('N','N',M,i,r,1.d0,Q,M,Ur,r,0.d0,U,M)
  DO j=1,i
    V(j,:)=S(j)*Vt(j,:)
  END DO
END IF
DEALLOCATE(perm,cn,cn0,tau,work,Q,Rs,S,Ur,Vt,IWORK)
END SUBROUTINE tw_lowrank_qrcp
!------------------------------------------------------------------------------
!> Point a single element block at element `ie` of `blk`
!------------------------------------------------------------------------------
SUBROUTINE tw_aca_set_elem(tw,one,blk,ie)
TYPE(tw_type), INTENT(in) :: tw !< Thin-wall model
TYPE(oft_tw_block), INTENT(inout) :: one !< Single element block
TYPE(oft_tw_block), INTENT(in) :: blk !< Block
INTEGER(4), INTENT(in) :: ie !< Element of the block
INTEGER(4) :: k1,k2
IF(.NOT.ASSOCIATED(one%ielem))THEN
  one%nelems=1
  ALLOCATE(one%ielem(1),one%icell(200),one%inv_map(tw%mesh%np))
  one%ielem=1
  one%icell=0
  one%inv_map=0
END IF
one%inv_map(tw%lpmap_inv(one%ielem(1)))=0
one%ielem=blk%ielem(ie)
one%inv_map(tw%lpmap_inv(one%ielem(1)))=1
k1=tw%mesh%kpc(tw%lpmap_inv(one%ielem(1)))
k2=tw%mesh%kpc(tw%lpmap_inv(one%ielem(1))+1)-1
one%ncells=k2-k1+1
one%icell(1:one%ncells)=tw%mesh%lpc(k1:k2)
END SUBROUTINE tw_aca_set_elem
!------------------------------------------------------------------------------
!> Row: element `i` of `iblock` against the elements of `jblock`
!------------------------------------------------------------------------------
SUBROUTINE tw_aca_Lblock_row(self,i,vec)
CLASS(tw_aca_Lblock), INTENT(inout) :: self !< Block source
INTEGER(4), INTENT(in) :: i !< Row (element of `iblock`)
REAL(8), INTENT(out) :: vec(:) !< Row of the block [ncols]
INTEGER(4) :: n
n=SIZE(vec)
CALL tw_aca_set_elem(self%op%tw_obj,self%one,self%op%levels(self%ilevel)%blocks(self%iblock),i)
CALL tw_compute_Lmatblock(self%op%tw_obj,self%op%tw_obj,self%buf(1:n,1:1),self%one, &
  self%op%levels(self%ilevel)%blocks(self%jblock))
vec=self%buf(1:n,1)
END SUBROUTINE tw_aca_Lblock_row
!------------------------------------------------------------------------------
!> Column: element `j` of `jblock` against the elements of `iblock`
!------------------------------------------------------------------------------
SUBROUTINE tw_aca_Lblock_col(self,i,vec)
CLASS(tw_aca_Lblock), INTENT(inout) :: self !< Block source
INTEGER(4), INTENT(in) :: i !< Column (element of `jblock`)
REAL(8), INTENT(out) :: vec(:) !< Column of the block [nrows]
INTEGER(4) :: n
n=SIZE(vec)
CALL tw_aca_set_elem(self%op%tw_obj,self%one,self%op%levels(self%ilevel)%blocks(self%jblock),i)
CALL tw_compute_Lmatblock(self%op%tw_obj,self%op%tw_obj,self%buf(1:n,1:1),self%one, &
  self%op%levels(self%ilevel)%blocks(self%iblock))
vec=self%buf(1:n,1)
END SUBROUTINE tw_aca_Lblock_col
!------------------------------------------------------------------------------
!> Row: field of element `i` of `iblock` at all points of `jblock`, stacked over the three directions
!------------------------------------------------------------------------------
SUBROUTINE tw_aca_B3block_row(self,i,vec)
CLASS(tw_aca_B3block), INTENT(inout) :: self !< Block source
INTEGER(4), INTENT(in) :: i !< Row (element of `iblock`)
REAL(8), INTENT(out) :: vec(:) !< Row of the block [ncols]
INTEGER(4) :: d,M
REAL(8), ALLOCATABLE :: B3(:,:,:)
M=self%op%levels(self%ilevel)%blocks(self%jblock)%np
CALL tw_aca_set_elem(self%op%tw_obj,self%one_elem,self%op%levels(self%ilevel)%blocks(self%iblock),i)
ALLOCATE(B3(M,1,3))
CALL tw_compute_Bops_block3(self%op%tw_obj,B3,self%one_elem,self%op%levels(self%ilevel)%blocks(self%jblock), &
  self%op%tw_obj%B_dx)
DO d=1,3
  vec((d-1)*M+1:d*M)=B3(:,1,d)
END DO
DEALLOCATE(B3)
END SUBROUTINE tw_aca_B3block_row
!------------------------------------------------------------------------------
!> Column: field direction `(i-1)/M+1` at point `MOD(i-1,M)+1` of `jblock` for all elements of `iblock`
!!
!! All three directions are computed together and kept, so the other directions of the same point are free.
!------------------------------------------------------------------------------
SUBROUTINE tw_aca_B3block_col(self,i,vec)
CLASS(tw_aca_B3block), INTENT(inout) :: self !< Block source
INTEGER(4), INTENT(in) :: i !< Stacked column (point and direction)
REAL(8), INTENT(out) :: vec(:) !< Column of the block [nrows]
INTEGER(4) :: M,N,p
REAL(8), ALLOCATABLE :: B3(:,:,:),tmp(:,:,:)
M=self%op%levels(self%ilevel)%blocks(self%jblock)%np
N=self%op%levels(self%ilevel)%blocks(self%iblock)%nelems
p=MOD(i-1,M)+1
IF(.NOT.ALLOCATED(self%slot))THEN
  ALLOCATE(self%slot(M),self%cols(N,3,16),self%one_pt%ipts(1))
  self%slot=0
  self%one_pt%np=1
END IF
IF(self%slot(p)==0)THEN
  IF(self%nslot==SIZE(self%cols,3))THEN
    ALLOCATE(tmp(N,3,2*self%nslot))
    tmp(:,:,1:self%nslot)=self%cols
    CALL MOVE_ALLOC(tmp,self%cols)
  END IF
  self%nslot=self%nslot+1
  self%slot(p)=self%nslot
  self%one_pt%ipts=self%op%levels(self%ilevel)%blocks(self%jblock)%ipts(p)
  ALLOCATE(B3(1,N,3))
  CALL tw_compute_Bops_block3(self%op%tw_obj,B3,self%op%levels(self%ilevel)%blocks(self%iblock),self%one_pt, &
    self%op%tw_obj%B_dx)
  self%cols(:,:,self%nslot)=B3(1,:,:)
  DEALLOCATE(B3)
END IF
vec=self%cols(:,(i-1)/M+1,self%slot(p))
END SUBROUTINE tw_aca_B3block_col
!------------------------------------------------------------------------------
!> Kernel between cells `jc` (target) and `ic` (source)
!------------------------------------------------------------------------------
FUNCTION tw_aca_cell_kernel(mesh,jc,ic) RESULT(val)
CLASS(oft_bmesh), INTENT(in) :: mesh !< Mesh
INTEGER(4), INTENT(in) :: jc !< Target cell
INTEGER(4), INTENT(in) :: ic !< Source cell
REAL(8) :: val !< Cell-pair kernel (@ref thin_wall::tw_cell_pair_kernel "tw_cell_pair_kernel")
REAL(8) :: pi_(3,3),pj_(3,3)
INTEGER(4) :: m
DO m=1,3
  pi_(:,m)=mesh%r(:,mesh%lc(m,ic))
  pj_(:,m)=mesh%r(:,mesh%lc(m,jc))
END DO
val=tw_cell_pair_kernel(pi_,mesh%ca(ic),pj_,mesh%ca(jc))
END FUNCTION tw_aca_cell_kernel
!------------------------------------------------------------------------------
!> Row: row cell `i` against all column cells
!------------------------------------------------------------------------------
SUBROUTINE tw_aca_cellblock_row(self,i,vec)
CLASS(tw_aca_cellblock), INTENT(inout) :: self !< Block source
INTEGER(4), INTENT(in) :: i !< Row (index in `rcells`)
REAL(8), INTENT(out) :: vec(:) !< Row of the block [ncols]
INTEGER(4) :: j
DO j=1,SIZE(self%ccells)
  vec(j)=tw_aca_cell_kernel(self%mesh,self%rcells(i),self%ccells(j))
END DO
END SUBROUTINE tw_aca_cellblock_row
!------------------------------------------------------------------------------
!> Column: all row cells against column cell `i`
!------------------------------------------------------------------------------
SUBROUTINE tw_aca_cellblock_col(self,i,vec)
CLASS(tw_aca_cellblock), INTENT(inout) :: self !< Block source
INTEGER(4), INTENT(in) :: i !< Column (index in `ccells`)
REAL(8), INTENT(out) :: vec(:) !< Column of the block [nrows]
INTEGER(4) :: j
DO j=1,SIZE(self%rcells)
  vec(j)=tw_aca_cell_kernel(self%mesh,self%rcells(j),self%ccells(i))
END DO
END SUBROUTINE tw_aca_cellblock_col
!
SUBROUTINE get_qr(Mat,Q,R)
REAL(8), INTENT(in) :: Mat(:,:)
REAL(8), ALLOCATABLE, INTENT(out) :: Q(:,:),R(:,:)
!
INTEGER(4) :: i,j,M,N,K,LWORK,INFO
REAL(8), ALLOCATABLE, DIMENSION(:) :: WORK,TAU
M = SIZE(Mat,DIM=1)
N = SIZE(Mat,DIM=2)
K = MIN(M,N)
ALLOCATE(Q(M,N),TAU(MIN(M,N)))
Q = Mat
LWORK = -1
ALLOCATE(WORK(1))
CALL dgeqrf( M, N, Q, M, TAU, WORK, LWORK, INFO )
IF(INFO/=0)CALL oft_abort("QR factorization failed","tw_Lmat_MF_setup::get_qr",__FILE__)
LWORK = INT(WORK(1))
DEALLOCATE(WORK)
ALLOCATE(WORK(LWORK))
CALL dgeqrf( M, N, Q, M, TAU, WORK, LWORK, INFO )
IF(INFO/=0)CALL oft_abort("QR factorization failed","tw_Lmat_MF_setup::get_qr",__FILE__)
DEALLOCATE(WORK)
ALLOCATE(R(K,K))
R=0.d0
DO i=1,M
  DO j=i,N
    R(i,j)=Q(i,j)
  END DO
END DO
ALLOCATE(WORK(1))
LWORK=-1
CALL DORGQR( M, N, N, Q, M, TAU, WORK, LWORK, INFO )
IF(INFO/=0)CALL oft_abort("QR factorization failed","tw_Lmat_MF_setup::get_qr",__FILE__)
LWORK = INT(WORK(1))
DEALLOCATE(WORK)
ALLOCATE(WORK(LWORK))
CALL DORGQR( M, N, N, Q, M, TAU, WORK, LWORK, INFO )
IF(INFO/=0)CALL oft_abort("QR factorization failed","tw_Lmat_MF_setup::get_qr",__FILE__)
DEALLOCATE(WORK,TAU)
END SUBROUTINE get_qr
!------------------------------------------------------------------------------
!> Recompress a low-rank factorization \f$ U V \f$ (from ACA+) with the SVD truncation rule
!!
!! QR factorizations \f$ U = Q_U R_U \f$, \f$ V^T = Q_V R_V \f$ and an SVD of the small core \f$ R_U R_V^T \f$; the
!! factors are replaced by those of the smallest rank whose discarded singular values have a Frobenius norm below
!! `tol` (at least `min_rank`), with the singular values included in U.
!------------------------------------------------------------------------------
SUBROUTINE tw_lowrank_recompress(Umat,Vmat,tol,min_rank,warn_rank,size_out)
TYPE(oft_native_dense_matrix), INTENT(inout) :: Umat !< Left factor [M,k], replaced by [M,rank]
TYPE(oft_native_dense_matrix), INTENT(inout) :: Vmat !< Right factor [k,N], replaced by [rank,N]
REAL(8), INTENT(in) :: tol !< Absolute tolerance on the Frobenius norm of the discarded part
INTEGER(4), INTENT(in) :: min_rank !< Minimum rank to keep
INTEGER(4), INTENT(in) :: warn_rank !< Warn if the rank is not reduced and k is above this
INTEGER(8), INTENT(out) :: size_out !< Size of the new factors, rank*(M+N)
INTEGER(4) :: i,j,M,N,K,LWORK,INFO
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: IWORK
REAL(8) :: frob_K
REAL(8), ALLOCATABLE, DIMENSION(:) :: S,WORK
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: U,VT,Atmp,QU,QV,RU,RV
!---QR factorizations of U and V^T
CALL get_qr(Umat%M,QU,RU)
ALLOCATE(Atmp(SIZE(Vmat%M,2),SIZE(Vmat%M,1)))
Atmp=TRANSPOSE(Vmat%M)
CALL get_qr(Atmp,QV,RV)
DEALLOCATE(Atmp)
M=SIZE(QU,DIM=1)
N=SIZE(QV,DIM=1)
K=SIZE(RU,DIM=1)
ALLOCATE(Atmp(K,K))
CALL DGEMM('N','T',K,K,K,1.d0,RU,K,RV,K,0.d0,Atmp,K)
!---SVD of the core
ALLOCATE(S(K),U(K,K),VT(K,K),WORK(1),IWORK(8*K))
LWORK=-1
CALL DGESDD('S',K,K,Atmp,K,S,U,K,VT,K,WORK,LWORK,IWORK,INFO)
IF(INFO/=0)CALL oft_abort("The algorithm computing SVD failed to converge.","tw_lowrank_recompress",__FILE__)
LWORK=INT(WORK(1))
DEALLOCATE(WORK)
ALLOCATE(WORK(LWORK))
CALL DGESDD('S',K,K,Atmp,K,S,U,K,VT,K,WORK,LWORK,IWORK,INFO)
IF(INFO/=0)CALL oft_abort("The algorithm computing SVD failed to converge.","tw_lowrank_recompress",__FILE__)
!---Truncate at desired accuracy
frob_K=0.d0
DO i=K,1,-1
  frob_K=frob_K+S(i)**2
  IF(SQRT(frob_K)>tol)EXIT
END DO
i=MAX(MIN(min_rank,K),i)
IF((i==K).AND.(K>warn_rank))WRITE(*,'(2A,3ES14.6)')oft_indent,'SVD recompression failed',S(1),S(K),tol
!---Replace U and V: U = Q_U U_r S_r, V = VT_r Q_V^T
DEALLOCATE(Umat%M,Vmat%M)
ALLOCATE(Umat%M(M,i),Vmat%M(i,N))
DO j=1,i
  U(:,j)=U(:,j)*S(j)
END DO
CALL DGEMM('N','N',M,i,K,1.d0,QU,M,U,K,0.d0,Umat%M,M)
CALL DGEMM('N','T',i,N,K,1.d0,VT,K,QV,N,0.d0,Vmat%M,i)
size_out=INT(i,8)*INT(M+N,8)
DEALLOCATE(WORK,IWORK,S,U,VT,Atmp,QU,QV,RU,RV)
END SUBROUTINE tw_lowrank_recompress
!------------------------------------------------------------------------------
!> Precondition a linear system using a Block-Jacobi method (complex)
!!
!! Each diagonal block is inverted explicitly using its LU factorization. If `pre_single`
!! is set on the HODLR operator each block is rounded to single precision before
!! factorization, and the inverse is computed, stored, and applied in single precision.
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE bjprecond_apply(self,u,g)
CLASS(oft_tw_hodlr_bjpre), INTENT(inout) :: self !< Solver object
CLASS(oft_cvector), INTENT(inout) :: u !< Guess/Solution field
CLASS(oft_cvector), INTENT(inout) :: g !< RHS/Residual field
!---
INTEGER(4) :: i,j,k,l,n,level,iblock,max_size
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: imap,hole_inds
COMPLEX(8), POINTER, DIMENSION(:) :: utmp,gtmp,ue,ge
COMPLEX(8), ALLOCATABLE, DIMENSION(:,:) :: vloc
COMPLEX(4), ALLOCATABLE, DIMENSION(:,:) :: vloc4
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%inverse_mats))ALLOCATE(self%inverse_mats(self%mf_obj%ndense+1))
IF(self%refactor)THEN
  max_size = 0
  !$omp parallel private(imap,i,j,k,l,n,level,iblock) reduction(MAX:max_size)
  ALLOCATE(imap(tw_nelems_exp(self%mf_obj%tw_obj)))
  !$omp do schedule(static)
  DO i=1,self%mf_obj%ndense
    level=self%mf_obj%dense_blocks(1,i)
    iblock=self%mf_obj%dense_blocks(2,i)
    n=self%mf_obj%levels(level)%blocks(iblock)%nelems
    max_size = MAX(max_size,n)
    IF(.NOT.ASSOCIATED(self%inverse_mats(i)%M))ALLOCATE(self%inverse_mats(i)%M(n,n))
    ! Explicit loop avoids array temporaries on (small) OpenMP thread stacks
    !$omp simd collapse(2)
    DO k=1,n
      DO j=1,n
        self%inverse_mats(i)%M(j,k)=self%alpha*self%mf_obj%dense_mats(i)%M(j,k)
      END DO
    END DO
    imap=0
    DO j=1,n
      imap(self%mf_obj%levels(level)%blocks(iblock)%ielem(j))=j
    END DO
    DO j=1,n
      k=self%mf_obj%levels(level)%blocks(iblock)%ielem(j)
      DO l=self%Rmat%kr(k),self%Rmat%kr(k+1)-1
        IF(imap(self%Rmat%lc(l))>0)THEN
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) = &
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) + self%beta*self%Rmat%M(l)
        END IF
      END DO
    END DO
    CALL invert_block(self%inverse_mats(i),n)
  END DO
  !$omp end do
  DEALLOCATE(imap)
  !$omp end parallel
  self%max_block_size = max_size
  !---Get holes
  ALLOCATE(imap(tw_nelems_exp(self%mf_obj%tw_obj)))
  n=self%mf_obj%tw_obj%nholes+self%mf_obj%tw_obj%n_vcoils
  IF(n>0)THEN
    i = self%mf_obj%ndense+1
    self%max_block_size = MAX(self%max_block_size,n)
    IF(.NOT.ASSOCIATED(self%inverse_mats(i)%M))ALLOCATE(self%inverse_mats(i)%M(n,n))
    imap=0
    DO j=1,n
      DO k=1,n
        self%inverse_mats(i)%M(j,k) = &
          self%alpha*self%mf_obj%hole_Vcoil_mat%M(j+self%mf_obj%tw_obj%np_active,k)
      END DO
      imap(j+self%mf_obj%tw_obj%np_active)=j
    END DO
    DO j=1,n
      k=j+self%mf_obj%tw_obj%np_active
      DO l=self%Rmat%kr(k),self%Rmat%kr(k+1)-1
        IF(imap(self%Rmat%lc(l))>0)THEN
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) = &
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) + self%beta*self%Rmat%M(l)
        END IF
      END DO
    END DO
    CALL invert_block(self%inverse_mats(i),n)
  END IF
  DEALLOCATE(imap)
  self%refactor=.FALSE.
END IF
!
NULLIFY(utmp,gtmp)
CALL u%get_local(utmp)
CALL g%get_local(gtmp)
utmp=(0.d0,0.d0)
IF(ASSOCIATED(self%mf_obj%tw_obj%junc))THEN
  ALLOCATE(ue(tw_nelems_exp(self%mf_obj%tw_obj)),ge(tw_nelems_exp(self%mf_obj%tw_obj)))
  CALL tw_junction_expand_c(self%mf_obj%tw_obj,gtmp,ge)
  ue=(0.d0,0.d0)
ELSE
  ue=>utmp
  ge=>gtmp
END IF
n = self%mf_obj%tw_obj%nholes+self%mf_obj%tw_obj%n_vcoils
ALLOCATE(hole_inds(n))
DO j=1,n
  hole_inds(j)=self%mf_obj%tw_obj%np_active+j
END DO
!$omp parallel private(vloc,vloc4,i,level,iblock,n)
IF(self%mf_obj%pre_single)THEN
  ALLOCATE(vloc4(self%max_block_size,2))
ELSE
  ALLOCATE(vloc(self%max_block_size,2))
END IF
DO i=self%mf_obj%L_thread_ptr(2,oft_tid+1),self%mf_obj%L_thread_ptr(2,oft_tid+2)-1
  level=self%mf_obj%dense_blocks(1,i)
  iblock=self%mf_obj%dense_blocks(2,i)
  n = self%mf_obj%levels(level)%blocks(iblock)%nelems
  CALL apply_block(self%inverse_mats(i),n,self%mf_obj%levels(level)%blocks(iblock)%ielem,vloc,vloc4)
END DO
!$omp single
n = self%mf_obj%tw_obj%nholes+self%mf_obj%tw_obj%n_vcoils
IF(n>0)CALL apply_block(self%inverse_mats(self%mf_obj%ndense+1),n,hole_inds,vloc,vloc4)
!$omp end single
IF(self%mf_obj%pre_single)THEN
  DEALLOCATE(vloc4)
ELSE
  DEALLOCATE(vloc)
END IF
!$omp end parallel
IF(ASSOCIATED(self%mf_obj%tw_obj%junc))THEN
  CALL tw_junction_reduce_vec_c(self%mf_obj%tw_obj,ue,utmp)
  DEALLOCATE(ue,ge)
END IF
CALL u%restore_local(utmp)
CALL g%set((0.d0,0.d0))
DEALLOCATE(utmp,gtmp,hole_inds)
DEBUG_STACK_POP
CONTAINS
!------------------------------------------------------------------------------
!> Invert assembled block, converting to single precision if requested
!------------------------------------------------------------------------------
SUBROUTINE invert_block(mat,nrows)
TYPE(cmat_container), INTENT(inout) :: mat !< Block to invert (`M` must be populated)
INTEGER(4), INTENT(in) :: nrows !< Size of block
INTEGER(4) :: error,j,k
IF(self%mf_obj%pre_single)THEN
  IF(.NOT.ASSOCIATED(mat%M4))ALLOCATE(mat%M4(nrows,nrows))
  ! Explicit loop avoids array temporaries on (small) OpenMP thread stacks
  !$omp simd collapse(2)
  DO k=1,nrows
    DO j=1,nrows
      mat%M4(j,k)=CMPLX(mat%M(j,k),KIND=c4)
    END DO
  END DO
  DEALLOCATE(mat%M)
  CALL lapack_matinv(nrows,mat%M4,error)
  IF(error/=0)CALL oft_abort("Single precision block inversion failed, block may be ill-conditioned (try pre_single=F)", &
    "bjprecond_apply",__FILE__)
ELSE
  CALL lapack_matinv(nrows,mat%M,error)
  IF(error/=0)CALL oft_abort("Block inversion failed","bjprecond_apply",__FILE__)
END IF
END SUBROUTINE invert_block
!------------------------------------------------------------------------------
!> Apply preconditioner for a single block
!------------------------------------------------------------------------------
SUBROUTINE apply_block(mat,nrows,inds,vloc,vloc4)
TYPE(cmat_container), INTENT(in) :: mat !< Inverted block
INTEGER(4), INTENT(in) :: nrows !< Size of block
INTEGER(4), INTENT(in) :: inds(nrows) !< Global indices of block rows
COMPLEX(8), ALLOCATABLE, INTENT(inout) :: vloc(:,:) !< Double precision work array [max_block_size,2]
COMPLEX(4), ALLOCATABLE, INTENT(inout) :: vloc4(:,:) !< Single precision work array [max_block_size,2]
INTEGER(4) :: j
REAL(8) :: scale
IF(self%mf_obj%pre_single)THEN
  !---Scale RHS to avoid underflow/overflow in single precision
  scale=0.d0
  DO j=1,nrows
    scale=MAX(scale,ABS(REAL(ge(inds(j)),8)),ABS(AIMAG(ge(inds(j)))))
  END DO
  IF(scale<=0.d0)RETURN ! Zero RHS for this block (corresponding entries of ue are already zero)
  DO j=1,nrows
    vloc4(j,1) = CMPLX(ge(inds(j))/scale,KIND=c4)
  END DO
  CALL cgemv('N',nrows,nrows,(1.0_r4,0.0_r4),mat%M4,nrows,vloc4(:,1),1,(0.0_r4,0.0_r4),vloc4(:,2),1)
  DO j=1,nrows
    ue(inds(j)) = CMPLX(vloc4(j,2),KIND=c8)*scale
  END DO
ELSE
  DO j=1,nrows
    vloc(j,1) = ge(inds(j))
  END DO
  CALL zgemv('N',nrows,nrows,(1.d0,0.d0),mat%M,nrows,vloc(:,1),1,(0.d0,0.d0),vloc(:,2),1)
  DO j=1,nrows
    ue(inds(j)) = vloc(j,2)
  END DO
END IF
END SUBROUTINE apply_block
END SUBROUTINE bjprecond_apply
!------------------------------------------------------------------------------
!> Destroy Block-Jacobi preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine bjprecond_delete(self,propogate)
class(oft_tw_hodlr_bjpre), intent(inout) :: self
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
INTEGER(i4) :: i
DEBUG_STACK_PUSH
!---Destroy local solvers
IF(ASSOCIATED(self%inverse_mats))THEN
  DO i=1,SIZE(self%inverse_mats)
    IF(ASSOCIATED(self%inverse_mats(i)%M))DEALLOCATE(self%inverse_mats(i)%M)
    IF(ASSOCIATED(self%inverse_mats(i)%M4))DEALLOCATE(self%inverse_mats(i)%M4)
  END DO
  DEALLOCATE(self%inverse_mats)
END IF
!---Reset
self%max_block_size=0
self%refactor=.TRUE.
NULLIFY(self%mf_obj)
DEBUG_STACK_POP
end subroutine bjprecond_delete
!------------------------------------------------------------------------------
!> Precondition a linear system using a Block-Jacobi method (real)
!!
!! Each diagonal block is inverted explicitly using its Cholesky factorization. If `pre_single`
!! is set on the HODLR operator each block is rounded to single precision before
!! factorization, and the inverse is computed, stored, and applied in single precision.
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE rbjprecond_apply(self,u,g)
CLASS(oft_tw_hodlr_rbjpre), INTENT(inout) :: self !< Solver object
CLASS(oft_vector), INTENT(inout) :: u !< Guess/Solution field
CLASS(oft_vector), INTENT(inout) :: g !< RHS/Residual field
!---
INTEGER(4) :: i,j,k,l,n,level,iblock,max_size,info
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: imap,hole_inds
REAL(8), POINTER, DIMENSION(:) :: utmp,gtmp,ue,ge
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: vloc
REAL(4), ALLOCATABLE, DIMENSION(:,:) :: vloc4
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%inverse_mats))ALLOCATE(self%inverse_mats(self%mf_obj%ndense+1))
IF(self%refactor)THEN
  max_size = 0
  !$omp parallel private(imap,i,j,k,l,n,level,iblock) reduction(MAX:max_size)
  ALLOCATE(imap(tw_nelems_exp(self%mf_obj%tw_obj)))
  !$omp do schedule(static)
  DO i=1,self%mf_obj%ndense
    level=self%mf_obj%dense_blocks(1,i)
    iblock=self%mf_obj%dense_blocks(2,i)
    n=self%mf_obj%levels(level)%blocks(iblock)%nelems
    max_size = MAX(max_size,n)
    IF(.NOT.ASSOCIATED(self%inverse_mats(i)%M))ALLOCATE(self%inverse_mats(i)%M(n,n))
    ! Explicit loop avoids array temporaries on (small) OpenMP thread stacks
    !$omp simd collapse(2)
    DO k=1,n
      DO j=1,n
        self%inverse_mats(i)%M(j,k)=self%alpha*self%mf_obj%dense_mats(i)%M(j,k)
      END DO
    END DO
    imap=0
    DO j=1,n
      imap(self%mf_obj%levels(level)%blocks(iblock)%ielem(j))=j
    END DO
    DO j=1,n
      k=self%mf_obj%levels(level)%blocks(iblock)%ielem(j)
      DO l=self%Rmat%kr(k),self%Rmat%kr(k+1)-1
        IF(imap(self%Rmat%lc(l))>0)THEN
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) = &
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) + self%beta*self%Rmat%M(l)
        END IF
      END DO
    END DO
    CALL invert_block(self%inverse_mats(i),n)
  END DO
  !$omp end do
  DEALLOCATE(imap)
  !$omp end parallel
  self%max_block_size = max_size
  !---Get holes
  ALLOCATE(imap(tw_nelems_exp(self%mf_obj%tw_obj)))
  n=self%mf_obj%tw_obj%nholes+self%mf_obj%tw_obj%n_vcoils
  IF(n>0)THEN
    i = self%mf_obj%ndense+1
    self%max_block_size = MAX(self%max_block_size,n)
    IF(.NOT.ASSOCIATED(self%inverse_mats(i)%M))ALLOCATE(self%inverse_mats(i)%M(n,n))
    imap=0
    DO j=1,n
      DO k=1,n
        self%inverse_mats(i)%M(j,k) = &
          self%alpha*self%mf_obj%hole_Vcoil_mat%M(j+self%mf_obj%tw_obj%np_active,k)
      END DO
      imap(j+self%mf_obj%tw_obj%np_active)=j
    END DO
    DO j=1,n
      k=j+self%mf_obj%tw_obj%np_active
      DO l=self%Rmat%kr(k),self%Rmat%kr(k+1)-1
        IF(imap(self%Rmat%lc(l))>0)THEN
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) = &
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) + self%beta*self%Rmat%M(l)
        END IF
      END DO
    END DO
    CALL invert_block(self%inverse_mats(i),n)
  END IF
  DEALLOCATE(imap)
  !---Coarse correction: Z^T (alpha*L + beta*R) Z, from the stored Z^T L Z and Z^T R Z
  IF(self%mf_obj%pre_coarse)THEN
    IF(.NOT.ASSOCIATED(self%zlz))CALL rbjprecond_coarse_setup(self)
    IF(self%nz>0)THEN
      IF(.NOT.ASSOCIATED(self%einv))ALLOCATE(self%einv(self%nz,self%nz))
      self%einv=self%alpha*self%zlz+self%beta*self%zrz
      CALL lapack_cholesky(self%nz,self%einv,info)
      IF(info/=0)CALL oft_abort("coarse factorization failed","rbjprecond_apply",__FILE__)
    END IF
  END IF
  self%refactor=.FALSE.
END IF
!
NULLIFY(utmp,gtmp)
CALL u%get_local(utmp)
CALL g%get_local(gtmp)
utmp=0.d0
IF(ASSOCIATED(self%mf_obj%tw_obj%junc))THEN
  ALLOCATE(ue(tw_nelems_exp(self%mf_obj%tw_obj)),ge(tw_nelems_exp(self%mf_obj%tw_obj)))
  CALL tw_junction_expand(self%mf_obj%tw_obj,gtmp,ge)
  ue=0.d0
ELSE
  ue=>utmp
  ge=>gtmp
END IF
n = self%mf_obj%tw_obj%nholes+self%mf_obj%tw_obj%n_vcoils
ALLOCATE(hole_inds(n))
DO j=1,n
  hole_inds(j)=self%mf_obj%tw_obj%np_active+j
END DO
!$omp parallel private(vloc,vloc4,i,level,iblock,n)
IF(self%mf_obj%pre_single)THEN
  ALLOCATE(vloc4(self%max_block_size,2))
ELSE
  ALLOCATE(vloc(self%max_block_size,2))
END IF
DO i=self%mf_obj%L_thread_ptr(2,oft_tid+1),self%mf_obj%L_thread_ptr(2,oft_tid+2)-1
  level=self%mf_obj%dense_blocks(1,i)
  iblock=self%mf_obj%dense_blocks(2,i)
  n = self%mf_obj%levels(level)%blocks(iblock)%nelems
  CALL apply_block(self%inverse_mats(i),n,self%mf_obj%levels(level)%blocks(iblock)%ielem,vloc,vloc4)
END DO
!$omp single
n = self%mf_obj%tw_obj%nholes+self%mf_obj%tw_obj%n_vcoils
IF(n>0)CALL apply_block(self%inverse_mats(self%mf_obj%ndense+1),n,hole_inds,vloc,vloc4)
!$omp end single
IF(self%mf_obj%pre_single)THEN
  DEALLOCATE(vloc4)
ELSE
  DEALLOCATE(vloc)
END IF
!$omp end parallel
!---Coarse correction on the expanded residual/solution (same as gtmp/utmp without junctions)
IF(self%mf_obj%pre_coarse.AND.(self%nz>0))CALL rbjprecond_coarse_apply(self,ge,ue)
IF(ASSOCIATED(self%mf_obj%tw_obj%junc))THEN
  CALL tw_junction_reduce_vec(self%mf_obj%tw_obj,ue,utmp)
  DEALLOCATE(ue,ge)
END IF
CALL u%restore_local(utmp)
CALL g%set(0.d0)
DEALLOCATE(utmp,gtmp,hole_inds)
DEBUG_STACK_POP
CONTAINS
!------------------------------------------------------------------------------
!> Invert assembled block, converting to single precision if requested
!------------------------------------------------------------------------------
SUBROUTINE invert_block(mat,nrows)
TYPE(rmat_container), INTENT(inout) :: mat !< Block to invert (`M` must be populated)
INTEGER(4), INTENT(in) :: nrows !< Size of block
INTEGER(4) :: error,j,k
IF(self%mf_obj%pre_single)THEN
  IF(.NOT.ASSOCIATED(mat%M4))ALLOCATE(mat%M4(nrows,nrows))
  ! Explicit loop avoids array temporaries on (small) OpenMP thread stacks
  !$omp simd collapse(2)
  DO k=1,nrows
    DO j=1,nrows
      mat%M4(j,k)=REAL(mat%M(j,k),r4)
    END DO
  END DO
  DEALLOCATE(mat%M)
  CALL lapack_cholesky(nrows,mat%M4,error)
  IF(error/=0)CALL oft_abort("Single precision block factorization failed, block may be ill-conditioned (try pre_single=F)", &
    "rbjprecond_apply",__FILE__)
ELSE
  CALL lapack_cholesky(nrows,mat%M,error)
  IF(error/=0)CALL oft_abort("Block factorization failed","rbjprecond_apply",__FILE__)
END IF
END SUBROUTINE invert_block
!------------------------------------------------------------------------------
!> Apply preconditioner for a single block
!------------------------------------------------------------------------------
SUBROUTINE apply_block(mat,nrows,inds,vloc,vloc4)
TYPE(rmat_container), INTENT(in) :: mat !< Inverted block
INTEGER(4), INTENT(in) :: nrows !< Size of block
INTEGER(4), INTENT(in) :: inds(nrows) !< Global indices of block rows
REAL(8), ALLOCATABLE, INTENT(inout) :: vloc(:,:) !< Double precision work array [max_block_size,2]
REAL(4), ALLOCATABLE, INTENT(inout) :: vloc4(:,:) !< Single precision work array [max_block_size,2]
INTEGER(4) :: j
REAL(8) :: scale
IF(self%mf_obj%pre_single)THEN
  !---Scale RHS to avoid underflow/overflow in single precision
  scale=0.d0
  DO j=1,nrows
    scale=MAX(scale,ABS(ge(inds(j))))
  END DO
  IF(scale<=0.d0)RETURN ! Zero RHS for this block (corresponding entries of ue are already zero)
  DO j=1,nrows
    vloc4(j,1) = REAL(ge(inds(j))/scale,r4)
  END DO
  CALL sgemv('N',nrows,nrows,1.0_r4,mat%M4,nrows,vloc4(:,1),1,0.0_r4,vloc4(:,2),1)
  DO j=1,nrows
    ue(inds(j)) = REAL(vloc4(j,2),r8)*scale
  END DO
ELSE
  DO j=1,nrows
    vloc(j,1) = ge(inds(j))
  END DO
  CALL dgemv('N',nrows,nrows,1.d0,mat%M,nrows,vloc(:,1),1,0.d0,vloc(:,2),1)
  DO j=1,nrows
    ue(inds(j)) = vloc(j,2)
  END DO
END IF
END SUBROUTINE apply_block
END SUBROUTINE rbjprecond_apply
!------------------------------------------------------------------------------
!> Build the coarse space of the two-level preconditioner and its projections of L and R
!!
!! On each diagonal (leaf) block: a constant and the three coordinates of the element positions, orthonormalised
!! on the block (directions the block does not span are dropped); and a unit vector for each hole and V-coil.
!! These are the smooth, long-range current patterns block-Jacobi cannot represent. The global vectors in `z_glob`
!! (e.g. the slowest L/R eigenmodes), if present, are added after orthonormalisation against the rest of the space.
!!
!! The vectors are defined on the element indices of the blocks. For meshes split at junctions these are the
!! expanded elements, and the coarse space is \f$ Z = P^T Z_{exp} \f$: each vector is reduced before the L and R
!! products, and the products are expanded before the projection, so that \f$ Z^T A Z \f$ is exact. The global
!! vectors are the exception: they are given on the model unknowns and expanded with the right inverse of \f$ P^T \f$.
!------------------------------------------------------------------------------
SUBROUTINE rbjprecond_coarse_setup(self)
CLASS(oft_tw_hodlr_rbjpre), INTENT(inout) :: self !< Solver object
INTEGER(4) :: i,j,k,m,c,n,level,iblock,nmax,nunit,nfun,nglob
REAL(r8) :: nrm,nrm0,cen(3),ext
REAL(r8), ALLOCATABLE :: zb(:,:)
REAL(8), POINTER, DIMENSION(:) :: vals,wvals
REAL(8), ALLOCATABLE, DIMENSION(:) :: rvals,wexp
CLASS(oft_vector), POINTER :: zvec,azvec
TYPE(tw_type), POINTER :: tw
tw=>self%mf_obj%tw_obj
nunit=tw%nholes+tw%n_vcoils
nmax=0
DO i=1,self%mf_obj%ndense
  nmax=MAX(nmax,self%mf_obj%levels(self%mf_obj%dense_blocks(1,i))%blocks(self%mf_obj%dense_blocks(2,i))%nelems)
END DO
nfun=4
IF(self%coarse_order>1)nfun=10
nglob=0
IF(ASSOCIATED(self%z_glob))nglob=SIZE(self%z_glob,DIM=2)
ALLOCATE(self%z_block(nfun*self%mf_obj%ndense+nunit+nglob),self%z_unit(nfun*self%mf_obj%ndense+nunit+nglob))
ALLOCATE(self%z_vals(nmax,nfun*self%mf_obj%ndense+nunit),zb(nmax,nfun))
self%z_vals=0.d0
c=0
DO i=1,self%mf_obj%ndense
  level=self%mf_obj%dense_blocks(1,i)
  iblock=self%mf_obj%dense_blocks(2,i)
  n=self%mf_obj%levels(level)%blocks(iblock)%nelems
  cen=0.d0
  DO j=1,n
    cen=cen+tw%mesh%r(:,tw%lpmap_inv(self%mf_obj%levels(level)%blocks(iblock)%ielem(j)))
  END DO
  cen=cen/REAL(n,8)
  ext=1.d-30
  DO j=1,n
    zb(j,1)=1.d0
    zb(j,2:4)=tw%mesh%r(:,tw%lpmap_inv(self%mf_obj%levels(level)%blocks(iblock)%ielem(j)))-cen
    ext=MAX(ext,MAXVAL(ABS(zb(j,2:4))))
  END DO
  zb(1:n,2:4)=zb(1:n,2:4)/ext
  IF(nfun==10)THEN
    zb(1:n,5)=zb(1:n,2)**2; zb(1:n,6)=zb(1:n,3)**2; zb(1:n,7)=zb(1:n,4)**2
    zb(1:n,8)=zb(1:n,2)*zb(1:n,3); zb(1:n,9)=zb(1:n,2)*zb(1:n,4); zb(1:n,10)=zb(1:n,3)*zb(1:n,4)
  END IF
  !---Modified Gram-Schmidt on the block, dropping directions the block does not span
  nrm0=SQRT(REAL(n,8))
  DO m=1,nfun
    DO k=1,m-1
      IF(SUM(zb(1:n,k)**2)>0.d0)zb(1:n,m)=zb(1:n,m)-DOT_PRODUCT(zb(1:n,k),zb(1:n,m))*zb(1:n,k)
    END DO
    nrm=SQRT(SUM(zb(1:n,m)**2))
    IF(nrm>1.d-6*nrm0)THEN
      zb(1:n,m)=zb(1:n,m)/nrm
      c=c+1
      self%z_block(c)=i
      self%z_unit(c)=0
      self%z_vals(1:n,c)=zb(1:n,m)
    ELSE
      zb(1:n,m)=0.d0
    END IF
  END DO
END DO
DO j=1,nunit
  c=c+1
  self%z_block(c)=0
  self%z_unit(c)=tw%np_active+j
END DO
self%nz=c
self%nzg=0
DEALLOCATE(zb)
!---Global vectors, orthonormalised against the block/unit vectors (orthonormal) and each other
IF(nglob>0)THEN
  ALLOCATE(self%z_gvals(tw_nelems_exp(tw),nglob),zb(self%nz+nglob,1))
  k=0
  DO m=1,nglob
    !---Global vectors are given on the model unknowns; on junction meshes they are expanded with the right inverse
    !   of P^T (z = P (P^T P)^-1 x), so that P^T z = x before orthonormalisation
    IF(ASSOCIATED(tw%junc))THEN
      CALL junction_right_inverse(self%z_glob(:,m),self%z_gvals(:,k+1))
    ELSE
      self%z_gvals(:,k+1)=self%z_glob(:,m)
    END IF
    nrm0=SQRT(SUM(self%z_gvals(:,k+1)**2))
    IF(nrm0<=0.d0)CYCLE
    DO j=1,2
      CALL rbjprecond_coarse_restrict(self,self%z_gvals(:,k+1),zb(1:self%nz,1))
      DO i=1,self%nz
        CALL rbjprecond_coarse_expand(self,i,-zb(i,1),self%z_gvals(:,k+1))
      END DO
    END DO
    nrm=SQRT(SUM(self%z_gvals(:,k+1)**2))
    IF(nrm>1.d-6*nrm0)THEN
      self%z_gvals(:,k+1)=self%z_gvals(:,k+1)/nrm
      k=k+1
      self%nz=self%nz+1
      self%nzg=k
      self%z_block(self%nz)=-k
      self%z_unit(self%nz)=0
    END IF
  END DO
  DEALLOCATE(zb)
END IF
!---Z^T L Z and Z^T R Z, one L and R product per coarse vector
ALLOCATE(self%zlz(self%nz,self%nz),self%zrz(self%nz,self%nz))
CALL tw%Uloc%new(zvec)
CALL tw%Uloc%new(azvec)
NULLIFY(vals,wvals)
ALLOCATE(vals(tw_nelems_exp(tw)))
IF(ASSOCIATED(tw%junc))ALLOCATE(rvals(tw%nelems),wexp(tw_nelems_exp(tw)))
DO c=1,self%nz
  vals=0.d0
  CALL rbjprecond_coarse_expand(self,c,1.d0,vals)
  IF(ASSOCIATED(tw%junc))THEN
    CALL tw_junction_reduce_vec(tw,vals,rvals)
    CALL zvec%restore_local(rvals)
  ELSE
    CALL zvec%restore_local(vals)
  END IF
  !---L and R of the model (reduced for junction meshes; the preconditioner's own R is the expanded one)
  CALL self%mf_obj%apply(zvec,azvec)
  CALL azvec%get_local(wvals)
  CALL coarse_project(self%zlz(:,c))
  CALL tw%Rmat%apply(zvec,azvec)
  CALL azvec%get_local(wvals)
  CALL coarse_project(self%zrz(:,c))
END DO
DEALLOCATE(vals,wvals)
IF(ASSOCIATED(tw%junc))DEALLOCATE(rvals,wexp)
CALL zvec%delete(); CALL azvec%delete()
DEALLOCATE(zvec,azvec)
WRITE(*,'(2A,I6,A)')oft_indent,'Two-level preconditioner: ',self%nz,' coarse vectors'
CONTAINS
!> Z^T w for a product w of the model operator (expanded first for junction meshes)
SUBROUTINE coarse_project(zw)
REAL(r8), INTENT(out) :: zw(:)
IF(ASSOCIATED(tw%junc))THEN
  CALL tw_junction_expand(tw,wvals,wexp)
  CALL rbjprecond_coarse_restrict(self,wexp,zw)
ELSE
  CALL rbjprecond_coarse_restrict(self,wvals,zw)
END IF
END SUBROUTINE coarse_project
!---Expanded vector z = P (P^T P)^-1 x, so that P^T z = x (Jacobi-preconditioned CG on P^T P, which is sparse and SPD;
!   it is not diagonal, as hole and junction unknowns span many elements)
SUBROUTINE junction_right_inverse(x,z)
REAL(r8), INTENT(in) :: x(:) !< Vector on the model unknowns [nelems]
REAL(r8), INTENT(out) :: z(:) !< Expanded vector [nelems_exp]
INTEGER(4) :: it,kk
REAL(r8) :: rz,rzold,alf,xnrm
REAL(r8), ALLOCATABLE :: d(:),y(:),r(:),pv(:),ap(:),zr(:)
ALLOCATE(d(tw%nelems),y(tw%nelems),r(tw%nelems),pv(tw%nelems),ap(tw%nelems),zr(tw%nelems))
!---Jacobi preconditioner, diag(P^T P)
d=0.d0
DO kk=1,tw%junc%p_kr(tw%junc%nelems_exp+1)-1
  d(tw%junc%p_lc(kk))=d(tw%junc%p_lc(kk))+tw%junc%p_val(kk)**2
END DO
WHERE(d<=0.d0)d=1.d0
xnrm=SQRT(SUM(x**2))
y=0.d0; r=x
zr=r/d; pv=zr; rz=DOT_PRODUCT(r,zr)
DO it=1,1000
  IF(SQRT(SUM(r**2))<=1.d-12*xnrm)EXIT
  CALL tw_junction_expand(tw,pv,z)
  CALL tw_junction_reduce_vec(tw,z,ap)
  alf=rz/DOT_PRODUCT(pv,ap)
  y=y+alf*pv; r=r-alf*ap
  zr=r/d; rzold=rz; rz=DOT_PRODUCT(r,zr)
  pv=zr+(rz/rzold)*pv
END DO
IF(SQRT(SUM(r**2))>1.d-12*xnrm)CALL oft_warn('Junction expansion of a global coarse vector did not converge')
CALL tw_junction_expand(tw,y,z)
DEALLOCATE(d,y,r,pv,ap,zr)
END SUBROUTINE junction_right_inverse
END SUBROUTINE rbjprecond_coarse_setup
!------------------------------------------------------------------------------
!> Add a*z_c to an element vector
!------------------------------------------------------------------------------
SUBROUTINE rbjprecond_coarse_expand(self,c,a,vals)
CLASS(oft_tw_hodlr_rbjpre), INTENT(in) :: self !< Solver object
INTEGER(4), INTENT(in) :: c !< Index of the coarse vector
REAL(r8), INTENT(in) :: a !< Coefficient
REAL(r8), INTENT(inout) :: vals(:) !< Element vector [nelems_exp]
INTEGER(4) :: j,level,iblock
IF(self%z_block(c)<0)THEN
  vals=vals+a*self%z_gvals(:,-self%z_block(c))
ELSE IF(self%z_block(c)==0)THEN
  vals(self%z_unit(c))=vals(self%z_unit(c))+a
ELSE
  level=self%mf_obj%dense_blocks(1,self%z_block(c))
  iblock=self%mf_obj%dense_blocks(2,self%z_block(c))
  DO j=1,self%mf_obj%levels(level)%blocks(iblock)%nelems
    vals(self%mf_obj%levels(level)%blocks(iblock)%ielem(j))=vals(self%mf_obj%levels(level)%blocks(iblock)%ielem(j)) &
      + a*self%z_vals(j,c)
  END DO
END IF
END SUBROUTINE rbjprecond_coarse_expand
!------------------------------------------------------------------------------
!> Coarse projection w = Z^T v of an element vector
!------------------------------------------------------------------------------
SUBROUTINE rbjprecond_coarse_restrict(self,v,w)
CLASS(oft_tw_hodlr_rbjpre), INTENT(in) :: self !< Solver object
REAL(r8), CONTIGUOUS, INTENT(in) :: v(:) !< Element vector [nelems_exp]
REAL(r8), CONTIGUOUS, INTENT(out) :: w(:) !< Coarse vector [nz]
INTEGER(4) :: c,j,level,iblock
DO c=1,self%nz-self%nzg
  IF(self%z_block(c)==0)THEN
    w(c)=v(self%z_unit(c))
  ELSE
    level=self%mf_obj%dense_blocks(1,self%z_block(c))
    iblock=self%mf_obj%dense_blocks(2,self%z_block(c))
    w(c)=0.d0
    DO j=1,self%mf_obj%levels(level)%blocks(iblock)%nelems
      w(c)=w(c)+self%z_vals(j,c)*v(self%mf_obj%levels(level)%blocks(iblock)%ielem(j))
    END DO
  END IF
END DO
!---Global vectors together: w_g = Z_g^T v
IF(self%nzg>0)CALL dgemv('T',SIZE(self%z_gvals,1),self%nzg,1.d0,self%z_gvals,SIZE(self%z_gvals,1), &
  v,1,0.d0,w(self%nz-self%nzg+1),1)
END SUBROUTINE rbjprecond_coarse_restrict
!------------------------------------------------------------------------------
!> Add the coarse correction u += Z (Z^T A Z)^-1 Z^T g
!------------------------------------------------------------------------------
SUBROUTINE rbjprecond_coarse_apply(self,g,u)
CLASS(oft_tw_hodlr_rbjpre), INTENT(in) :: self !< Solver object
REAL(r8), CONTIGUOUS, INTENT(in) :: g(:) !< Residual [nelems_exp]
REAL(r8), CONTIGUOUS, INTENT(inout) :: u(:) !< Preconditioned residual [nelems_exp]
INTEGER(4) :: c
REAL(r8), ALLOCATABLE :: w(:),y(:)
ALLOCATE(w(self%nz),y(self%nz))
CALL rbjprecond_coarse_restrict(self,g,w)
CALL dgemv('N',self%nz,self%nz,1.d0,self%einv,self%nz,w,1,0.d0,y,1)
DO c=1,self%nz-self%nzg
  CALL rbjprecond_coarse_expand(self,c,y(c),u)
END DO
!---Global vectors together: u += Z_g y_g
IF(self%nzg>0)CALL dgemv('N',SIZE(self%z_gvals,1),self%nzg,1.d0,self%z_gvals,SIZE(self%z_gvals,1), &
  y(self%nz-self%nzg+1),1,1.d0,u,1)
DEALLOCATE(w,y)
END SUBROUTINE rbjprecond_coarse_apply
!------------------------------------------------------------------------------
!> Destroy Block-Jacobi preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine rbjprecond_delete(self,propogate)
class(oft_tw_hodlr_rbjpre), intent(inout) :: self
LOGICAL, optional, intent(in) :: propogate !< Update matrix non-zero pattern? (optional)
INTEGER(i4) :: i
DEBUG_STACK_PUSH
!---Destroy local solvers
IF(ASSOCIATED(self%inverse_mats))THEN
  DO i=1,SIZE(self%inverse_mats)
    IF(ASSOCIATED(self%inverse_mats(i)%M))DEALLOCATE(self%inverse_mats(i)%M)
    IF(ASSOCIATED(self%inverse_mats(i)%M4))DEALLOCATE(self%inverse_mats(i)%M4)
  END DO
  DEALLOCATE(self%inverse_mats)
END IF
IF(ASSOCIATED(self%z_block))DEALLOCATE(self%z_block,self%z_unit,self%z_vals)
NULLIFY(self%z_glob) ! Not owned
IF(ASSOCIATED(self%z_gvals))DEALLOCATE(self%z_gvals)
IF(ASSOCIATED(self%zlz))DEALLOCATE(self%zlz,self%zrz)
IF(ASSOCIATED(self%einv))DEALLOCATE(self%einv)
self%nz=0
!---Reset
self%max_block_size=0
self%refactor=.TRUE.
NULLIFY(self%mf_obj)
DEBUG_STACK_POP
end subroutine rbjprecond_delete
END MODULE thin_wall_hodlr
