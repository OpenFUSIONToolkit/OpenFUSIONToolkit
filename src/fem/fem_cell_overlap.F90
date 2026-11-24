!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file fem_cell_overlap.F90
!
!> Cell-adjacency overlapping subdomains for FE matrices (Schwarz methods)
!!
!! Implements @ref oft_la_base::oft_cell_overlap_builder "oft_cell_overlap_builder" for
!! a single FE representation or a composite one on a volume mesh. Each process owns
!! `nparts` blocks. A block is made of:
!! - **Core cells:** the local cells with cell part `cpart==b`, plus every local cell
!!   containing a local DoF of color `dcolor==b` (so that the DoFs assigned to a block
!!   are always in the block). With `nparts==1` the core is every local cell.
!! - **Block cells:** the core grown by `layers` layers of adjacent cells, where two
!!   cells are adjacent if they share at least one point, two points (an edge) or
!!   `face_np` points (a face) depending on `adj`. Layers are grown over the cells of
!!   all processes ("extended cell graph" made of the local cells plus halo cells
!!   received from other processes).
!! - **Block DoFs:** every DoF of every field in the closure of the block cells.
!! - **Assigned DoFs:** DoFs owned by this process with color `dcolor==b` (all owned
!!   DoFs if `nparts==1`). Every owned DoF is assigned to exactly one block.
!!
!! The cell part and DoF color come from the same call to
!! @ref oft_mesh_local::mesh_local_partition "mesh_local_partition" (METIS on the local
!! cells, then a vote for points, edges and faces), see
!! @ref fem_utils::fem_partition_cells "fem_partition_cells".
!!
!! Setup algorithm of @ref fem_cell_overlap_build (collective):
!! - A: foreign base cells within `base_layers` point layers of this process' base
!!   cells are found on the replicated base mesh (no communication,
!!   @ref oft_mesh_overlap::mesh_base_candidates "mesh_base_candidates").
!! - B: their fine cells (global cell id and global point ids) are requested from the
!!   owners (one request/reply `MPI_ALLTOALL(V)` exchange).
!! - C: block cells are grown over the extended cell graph (OpenMP over blocks).
!! - D: the DoFs (global row id and owner) of the halo cells used by any block are
!!   requested from the cells' owners (one request/reply exchange).
!! - E: the union DoF list and block/assigned positions are formed and validated.
!!
!! Steps B and D are skipped on shared meshes (`mesh%fullmesh`), where every cell is
!! local. Global row ids always come from the matrix row maps
!! (`i_map(f)%offsetg+i_map(f)%lge`), and owners from synchronizing the rank id by
!! ownership on the row space seam. Periodic row maps are not supported.
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup doxy_oft_fem
!---------------------------------------------------------------------------------
MODULE fem_cell_overlap
USE oft_base
USE oft_sort, ONLY: sort_array, search_array
USE oft_stitching, ONLY: oft_global_stitch
USE oft_scatter_base, ONLY: oft_exchange_i8, oft_exchange_report, oft_exchange_nsent
USE oft_mesh_type, ONLY: oft_mesh
USE multigrid, ONLY: multigrid_mesh
USE oft_mesh_overlap, ONLY: mesh_overlap_nchild, mesh_base_gid_inverse, mesh_base_candidates, &
  mesh_cell_base_index
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_cell_overlap_builder, oft_overlap_dofs, &
  overlap_dofs_delete, OFT_CELL_ADJ_POINT, OFT_CELL_ADJ_EDGE, OFT_CELL_ADJ_FACE
USE oft_native_la, ONLY: oft_native_vector, native_vector_cast
USE fem_base, ONLY: oft_afem_type, oft_fem_type, oft_fem_ptr
USE fem_composite, ONLY: oft_fem_comp_type
USE fem_utils, ONLY: fem_partition_cells, comp_fem_partition_cells
IMPLICIT NONE
#include "local.h"
PRIVATE
!---------------------------------------------------------------------------------
!> Sorted list of global cell ids
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: oft_cell_gid_set
  INTEGER(i4) :: n = 0 !< Number of cells
  INTEGER(i8), ALLOCATABLE, DIMENSION(:) :: v !< Global cell ids (sorted) [n]
END TYPE oft_cell_gid_set
!---------------------------------------------------------------------------------
!> Integer list (setup work storage)
!---------------------------------------------------------------------------------
TYPE :: ilist
  INTEGER(i4) :: n = 0 !< Number of values
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: v !< Values
END TYPE ilist
!---------------------------------------------------------------------------------
!> Cell-adjacency overlap builder for a single or composite FE representation
!!
!! Created with @ref fem_cell_overlap_create. The FE representation, mesh and
!! multi-level mesh are referenced, not owned. `core_cells` and `blk_cells` hold the
!! result of the last @ref fem_cell_overlap_build "build" for diagnostics.
!---------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_cell_overlap_builder) :: oft_fem_cell_overlap
  INTEGER(i4) :: nparts = 0 !< Number of blocks on this process
  INTEGER(i4) :: nfields = 0 !< Number of FE fields
  INTEGER(i4) :: ndofs = 0 !< Total number of local DoFs (all fields)
  INTEGER(i4) :: level = 0 !< Level of `mesh` in `ml_mesh`
  INTEGER(i4) :: nhalo = 0 !< Number of halo cells received in the last build
  INTEGER(i4) :: nhalo_used = 0 !< Number of halo cells used by a block in the last build
  TYPE(multigrid_mesh), POINTER :: ml_mesh => NULL() !< Multi-level mesh
  CLASS(oft_mesh), POINTER :: mesh => NULL() !< Mesh of the FE representation (one level of `ml_mesh`)
  CLASS(oft_afem_type), POINTER :: fe => NULL() !< Single FE representation (if not composite)
  CLASS(oft_fem_comp_type), POINTER :: fe_comp => NULL() !< Composite FE representation (if composite)
  TYPE(oft_fem_ptr), ALLOCATABLE, DIMENSION(:) :: fields !< FE representation of each field [nfields]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: cpart !< Block (part) of each local cell [mesh%nc]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: dcolor !< Block (color) of each local DoF [ndofs]
  TYPE(oft_cell_gid_set), ALLOCATABLE, DIMENSION(:) :: core_cells !< Global ids of each block's core cells [nparts]
  TYPE(oft_cell_gid_set), ALLOCATABLE, DIMENSION(:) :: blk_cells !< Global ids of each block's cells [nparts]
CONTAINS
  !> Build the overlapping DoF sets (collective)
  PROCEDURE :: build => fem_cell_overlap_build
  !> Number of blocks on this process
  PROCEDURE :: nblocks => fem_cell_overlap_nblocks
  !> Release builder resources
  PROCEDURE :: delete => fem_cell_overlap_delete
END TYPE oft_fem_cell_overlap
!---------------------------------------------------------------------------------
!> Create a cell-adjacency overlap builder
!!
!! Specific versions for a single FE representation (`CLASS(oft_afem_type)`, which
!! must be a volume `oft_fem_type`) and a composite one (`CLASS(oft_fem_comp_type)`).
!---------------------------------------------------------------------------------
INTERFACE fem_cell_overlap_create
  MODULE PROCEDURE fem_cell_overlap_create_single
  MODULE PROCEDURE fem_cell_overlap_create_comp
END INTERFACE fem_cell_overlap_create
PUBLIC fem_cell_overlap_create
CONTAINS
!---------------------------------------------------------------------------------
!> Create a cell-adjacency overlap builder for a single FE representation
!!
!! See @ref fem_cell_overlap_init for the arguments.
!---------------------------------------------------------------------------------
SUBROUTINE fem_cell_overlap_create_single(self,ml_mesh,fe,nparts,dcolor)
CLASS(oft_fem_cell_overlap), INTENT(inout) :: self !< Builder object
TYPE(multigrid_mesh), TARGET, INTENT(in) :: ml_mesh !< Multi-level mesh containing the FE mesh
CLASS(oft_afem_type), TARGET, INTENT(in) :: fe !< FE representation (volume `oft_fem_type`)
INTEGER(i4), INTENT(in) :: nparts !< Number of blocks on this process
INTEGER(i4), OPTIONAL, INTENT(out) :: dcolor(:) !< Block (color) of each local DoF [fe%ne] (optional)
DEBUG_STACK_PUSH
CALL fem_cell_overlap_delete(self)
self%fe=>fe
self%nfields=1
ALLOCATE(self%fields(1))
self%fields(1)%fe=>fe
CALL fem_cell_overlap_init(self,ml_mesh,nparts,dcolor)
DEBUG_STACK_POP
END SUBROUTINE fem_cell_overlap_create_single
!---------------------------------------------------------------------------------
!> Create a cell-adjacency overlap builder for a composite FE representation
!!
!! See @ref fem_cell_overlap_init for the arguments.
!---------------------------------------------------------------------------------
SUBROUTINE fem_cell_overlap_create_comp(self,ml_mesh,fe,nparts,dcolor)
CLASS(oft_fem_cell_overlap), INTENT(inout) :: self !< Builder object
TYPE(multigrid_mesh), TARGET, INTENT(in) :: ml_mesh !< Multi-level mesh containing the FE mesh
CLASS(oft_fem_comp_type), TARGET, INTENT(in) :: fe !< Composite FE representation (volume fields)
INTEGER(i4), INTENT(in) :: nparts !< Number of blocks on this process
INTEGER(i4), OPTIONAL, INTENT(out) :: dcolor(:) !< Block (color) of each local DoF [sum of field ne] (optional)
INTEGER(i4) :: i
DEBUG_STACK_PUSH
CALL fem_cell_overlap_delete(self)
self%fe_comp=>fe
self%nfields=fe%nfields
IF(self%nfields<1)CALL oft_abort('Composite FE without fields','fem_cell_overlap_create',__FILE__)
ALLOCATE(self%fields(self%nfields))
DO i=1,self%nfields
  self%fields(i)%fe=>fe%fields(i)%fe
END DO
CALL fem_cell_overlap_init(self,ml_mesh,nparts,dcolor)
DEBUG_STACK_POP
END SUBROUTINE fem_cell_overlap_create_comp
!---------------------------------------------------------------------------------
!> Common part of @ref fem_cell_overlap_create
!!
!! Checks that every field is a volume FE representation on the same mesh, which is
!! one of the levels of `ml_mesh`, and computes the cell part `cpart` and DoF color
!! `dcolor`. With `nparts>1` both come from one call to
!! @ref fem_utils::fem_partition_cells "fem_partition_cells" (or the composite
!! version). With `nparts==1` all cells and DoFs are in block 1 (no partitioning).
!! The DoF color is returned in `dcolor` (if present) so callers can use it as the
!! matrix coloring.
!---------------------------------------------------------------------------------
SUBROUTINE fem_cell_overlap_init(self,ml_mesh,nparts,dcolor)
CLASS(oft_fem_cell_overlap), INTENT(inout) :: self !< Builder object
TYPE(multigrid_mesh), TARGET, INTENT(in) :: ml_mesh !< Multi-level mesh containing the FE mesh
INTEGER(i4), INTENT(in) :: nparts !< Number of blocks on this process
INTEGER(i4), OPTIONAL, INTENT(out) :: dcolor(:) !< Block (color) of each local DoF [ndofs] (optional)
INTEGER(i4) :: i
DEBUG_STACK_PUSH
IF(nparts<1)CALL oft_abort('Number of blocks must be positive','fem_cell_overlap_create',__FILE__)
self%ml_mesh=>ml_mesh
self%nparts=nparts
self%ndofs=0
DO i=1,self%nfields
  SELECT TYPE(this=>self%fields(i)%fe)
  CLASS IS(oft_fem_type)
    IF(i==1)THEN
      self%mesh=>this%mesh
    ELSE
      IF(.NOT.ASSOCIATED(self%mesh,this%mesh))CALL oft_abort('Fields must share one mesh', &
        'fem_cell_overlap_create',__FILE__)
    END IF
    self%ndofs=self%ndofs+this%ne
  CLASS DEFAULT
    CALL oft_abort('Volume FE representation required','fem_cell_overlap_create',__FILE__)
  END SELECT
END DO
!---Level of the mesh (identified through its cell list, as polymorphic array elements
!   are not reliably comparable with ASSOCIATED)
self%level=0
DO i=1,ml_mesh%mgmax
  IF(ml_mesh%meshes(i)%nc/=self%mesh%nc)CYCLE
  IF(ml_mesh%meshes(i)%global%nc/=self%mesh%global%nc)CYCLE
  IF(self%mesh%nc>0)THEN
    IF(.NOT.ASSOCIATED(self%mesh%lc,ml_mesh%meshes(i)%lc))CYCLE
  END IF
  self%level=i
  EXIT
END DO
IF(self%level==0)CALL oft_abort('FE mesh is not a level of the multi-level mesh', &
  'fem_cell_overlap_create',__FILE__)
!---Cell parts and DoF colors
ALLOCATE(self%cpart(self%mesh%nc),self%dcolor(self%ndofs))
IF(nparts>1)THEN
  IF(ASSOCIATED(self%fe_comp))THEN
    CALL comp_fem_partition_cells(self%fe_comp,self%dcolor,nparts,self%cpart)
  ELSE
    SELECT TYPE(this=>self%fields(1)%fe)
    CLASS IS(oft_fem_type)
      CALL fem_partition_cells(this,self%dcolor,nparts,self%cpart)
    END SELECT
  END IF
  IF(self%mesh%nc>0)THEN
    IF(MINVAL(self%cpart)<1.OR.MAXVAL(self%cpart)>nparts)CALL oft_abort('Invalid cell partition', &
      'fem_cell_overlap_create',__FILE__)
  END IF
  IF(self%ndofs>0)THEN
    IF(MINVAL(self%dcolor)<1.OR.MAXVAL(self%dcolor)>nparts)CALL oft_abort('Invalid DoF coloring', &
      'fem_cell_overlap_create',__FILE__)
  END IF
ELSE
  self%cpart=1
  self%dcolor=1
END IF
IF(PRESENT(dcolor))THEN
  IF(SIZE(dcolor)<self%ndofs)CALL oft_abort('DoF color output too small','fem_cell_overlap_create',__FILE__)
  dcolor(1:self%ndofs)=self%dcolor
END IF
DEBUG_STACK_POP
END SUBROUTINE fem_cell_overlap_init
!---------------------------------------------------------------------------------
!> Number of blocks on this process
!---------------------------------------------------------------------------------
FUNCTION fem_cell_overlap_nblocks(self) RESULT(nb)
CLASS(oft_fem_cell_overlap), INTENT(inout) :: self !< Builder object
INTEGER(i4) :: nb !< Number of blocks
nb=self%nparts
END FUNCTION fem_cell_overlap_nblocks
!---------------------------------------------------------------------------------
!> Release builder resources (FE representation and meshes are only dereferenced)
!---------------------------------------------------------------------------------
SUBROUTINE fem_cell_overlap_delete(self)
CLASS(oft_fem_cell_overlap), INTENT(inout) :: self !< Builder object
DEBUG_STACK_PUSH
CALL clear_cell_sets(self)
IF(ALLOCATED(self%fields))DEALLOCATE(self%fields)
IF(ALLOCATED(self%cpart))DEALLOCATE(self%cpart)
IF(ALLOCATED(self%dcolor))DEALLOCATE(self%dcolor)
NULLIFY(self%ml_mesh,self%mesh,self%fe,self%fe_comp)
self%nparts=0
self%nfields=0
self%ndofs=0
self%level=0
self%nhalo=0
self%nhalo_used=0
DEBUG_STACK_POP
END SUBROUTINE fem_cell_overlap_delete
!---------------------------------------------------------------------------------
!> Release the diagnostic cell sets of the last build
!---------------------------------------------------------------------------------
SUBROUTINE clear_cell_sets(self)
CLASS(oft_fem_cell_overlap), INTENT(inout) :: self !< Builder object
IF(ALLOCATED(self%core_cells))DEALLOCATE(self%core_cells)
IF(ALLOCATED(self%blk_cells))DEALLOCATE(self%blk_cells)
END SUBROUTINE clear_cell_sets
!---------------------------------------------------------------------------------
!> Build the overlapping DoF sets of the blocks of this process (collective)
!!
!! Implements steps A to E described in @ref fem_cell_overlap. `mat` must be a square
!! matrix on the FE representation (one row block per field, in field order, with
!! local rows numbered as the FE DoFs) and `u` a native vector of its row space.
!!
!! Requirements checked (abort otherwise): no periodic row maps; row maps consistent
!! with the FE fields; unique global ids of local rows; uniform refinement above the
!! transfer level; every received reply has the expected size; consistent DoF owners;
!! every owned DoF assigned to exactly one block and assigned DoFs within their block.
!---------------------------------------------------------------------------------
SUBROUTINE fem_cell_overlap_build(self,mat,u,layers,adj,base_layers,dofs)
CLASS(oft_fem_cell_overlap), INTENT(inout) :: self !< Builder object
CLASS(oft_matrix), INTENT(inout) :: mat !< Matrix defining the row maps
CLASS(oft_vector), INTENT(inout) :: u !< Row-space vector (seam information)
INTEGER(i4), INTENT(in) :: layers !< Number of layers of cells to grow by
INTEGER(i4), INTENT(in) :: adj !< Cell adjacency (`OFT_CELL_ADJ_*`)
INTEGER(i4), INTENT(in) :: base_layers !< Number of base mesh point layers bounding the halo (<0: `layers`)
TYPE(oft_overlap_dofs), INTENT(inout) :: dofs !< Resulting DoF sets
CLASS(oft_native_vector), POINTER :: uv
CLASS(oft_mesh), POINTER :: mesh,bmesh
INTEGER(i4) :: i,j,k,l,m,b,e,e2,f,q,np,me,nc,ncp,nce_tot,thr,kb,ncand,nref,nfine,next,nhalo, &
  npair,nup,pass,ntl,nused,nu,nown,istart,iend,nb,ioff,row
INTEGER(i8) :: n8,nfine8
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: owner,cand,cown,cb,inv,scnt,rcnt,kbc,lbc,ibc,eown, &
  iperm,upos,kp,lpx,knb,lnb,nnb,cnt,tl,list,hit,nsel,lperm
INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: lrow,hdown
INTEGER(i8), ALLOCATABLE, DIMENSION(:) :: grow,gsorted,sbuf,rbuf,egid,pk,dk,lsorted
INTEGER(i8), ALLOCATABLE, DIMENSION(:,:) :: epts,hdof
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: dtmp
REAL(r8), ALLOCATABLE, DIMENSION(:) :: rown_r
LOGICAL, ALLOCATABLE, DIMENSION(:) :: needed,mark,inb
INTEGER(i4), POINTER, DIMENSION(:) :: lcpart
TYPE(ilist), ALLOCATABLE, DIMENSION(:) :: core,bx
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ldofs
TYPE(oft_timer) :: stimer
INTEGER(i8) :: nsent0
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%mesh))CALL oft_abort('Builder not created','fem_cell_overlap_build',__FILE__)
IF(.NOT.native_vector_cast(uv,u))CALL oft_abort('Native vector required','fem_cell_overlap_build',__FILE__)
IF(layers<0)CALL oft_abort('Negative number of layers','fem_cell_overlap_build',__FILE__)
SELECT CASE(adj)
CASE(OFT_CELL_ADJ_POINT)
  thr=1
CASE(OFT_CELL_ADJ_EDGE)
  thr=2
CASE(OFT_CELL_ADJ_FACE)
  thr=self%mesh%face_np
CASE DEFAULT
  CALL oft_abort('Invalid cell adjacency type','fem_cell_overlap_build',__FILE__)
END SELECT
CALL overlap_dofs_delete(dofs)
CALL clear_cell_sets(self)
CALL stimer%tick
nsent0=oft_exchange_nsent
mesh=>self%mesh
np=oft_env%nprocs
me=oft_env%rank
nc=mesh%nc
ncp=mesh%cell_np
nb=self%nparts
!---------------------------------------------------------------------------------
! Row maps: global id and owner of each local row
!---------------------------------------------------------------------------------
IF(mat%ni/=self%nfields)CALL oft_abort('Matrix row blocks do not match FE fields','fem_cell_overlap_build',__FILE__)
IF(mat%nr/=self%ndofs)CALL oft_abort('Matrix rows do not match FE DoFs','fem_cell_overlap_build',__FILE__)
ioff=0
nce_tot=0
DO f=1,self%nfields
  IF(mat%i_map(f)%per)CALL oft_abort('Periodic row maps are not supported','fem_cell_overlap_build',__FILE__)
  IF(mat%i_map(f)%n/=self%fields(f)%fe%ne.OR.mat%i_map(f)%offset/=ioff) &
    CALL oft_abort('Matrix row map does not match FE field','fem_cell_overlap_build',__FILE__)
  ioff=ioff+self%fields(f)%fe%ne
  nce_tot=nce_tot+self%fields(f)%fe%nce
END DO
ALLOCATE(grow(mat%nr))
DO f=1,self%nfields
  DO j=1,mat%i_map(f)%n
    grow(mat%i_map(f)%offset+j)=mat%i_map(f)%offsetg+mat%i_map(f)%lge(j)
  END DO
END DO
ALLOCATE(gsorted(mat%nr))
gsorted=grow
n8=INT(mat%nr,i8)
CALL sort_array(gsorted,n8)
DO i=2,mat%nr
  IF(gsorted(i)==gsorted(i-1))CALL oft_abort('Duplicate global row ids','fem_cell_overlap_build',__FILE__)
END DO
DEALLOCATE(gsorted)
ALLOCATE(rown_r(mat%nr),owner(mat%nr))
rown_r=REAL(me,r8)
CALL oft_global_stitch(uv%stitch_info,rown_r,0)
owner=NINT(rown_r)
DEALLOCATE(rown_r)
!---Local rows of each local cell (all fields, ncdofs order)
ALLOCATE(lrow(nce_tot,nc))
DO f=1,self%nfields
  ALLOCATE(ldofs(self%fields(f)%fe%nce))
  ioff=0
  DO i=1,f-1
    ioff=ioff+self%fields(i)%fe%nce
  END DO
  DO i=1,nc
    CALL self%fields(f)%fe%ncdofs(i,ldofs)
    lrow(ioff+1:ioff+self%fields(f)%fe%nce,i)=mat%i_map(f)%offset+ldofs
  END DO
  DEALLOCATE(ldofs)
END DO
!---------------------------------------------------------------------------------
! Steps A and B: halo cell topology from the owners of candidate base cells
!---------------------------------------------------------------------------------
nhalo=0
ALLOCATE(scnt(0:np-1),rcnt(0:np-1))
IF(.NOT.mesh%fullmesh)THEN
  bmesh=>self%ml_mesh%meshes(self%ml_mesh%nbase)
  lcpart=>self%ml_mesh%meshes(self%ml_mesh%nbase+1)%base%lcpart
  kb=layers
  IF(base_layers>=0)kb=base_layers
  CALL mesh_base_candidates(self%ml_mesh,kb,ncand,cand,cown)
  CALL mesh_cell_base_index(self%ml_mesh,mesh,cb,nref)
  CALL mesh_base_gid_inverse(self%ml_mesh,inv)
  nfine=mesh_overlap_nchild**nref
  !---Request candidate base cells (base global ids) from their owners
  scnt=0
  ALLOCATE(sbuf(MAX(1,ncand)))
  DO i=1,ncand
    scnt(cown(i))=scnt(cown(i))+1
    sbuf(i)=bmesh%global%lc(cand(i))
  END DO
  CALL oft_exchange_i8(scnt,sbuf,rcnt,rbuf)
  DEALLOCATE(sbuf)
  !---Local cells of each base cell
  ALLOCATE(kbc(bmesh%nc+1),lbc(MAX(1,nc)))
  kbc=0
  DO i=1,nc
    kbc(cb(i)+1)=kbc(cb(i)+1)+1
  END DO
  kbc(1)=1
  DO i=1,bmesh%nc
    kbc(i+1)=kbc(i+1)+kbc(i)
  END DO
  ALLOCATE(ibc(bmesh%nc))
  ibc=kbc(1:bmesh%nc)
  DO i=1,nc
    lbc(ibc(cb(i)))=i
    ibc(cb(i))=ibc(cb(i))+1
  END DO
  DEALLOCATE(ibc)
  !---Reply with the fine cells of each requested base cell: [gid, cell_np point gids]
  m=SUM(rcnt)
  ALLOCATE(sbuf(MAX(1,m*nfine*(1+ncp))))
  k=0
  l=0
  DO q=0,np-1
    scnt(q)=0
    DO i=1,rcnt(q)
      l=l+1
      IF(rbuf(l)<1_i8.OR.rbuf(l)>INT(bmesh%nc,i8))CALL oft_abort('Invalid base cell request', &
        'fem_cell_overlap_build',__FILE__)
      j=inv(rbuf(l))
      IF(lcpart(j)/=me+1)CALL oft_abort('Base cell requested from wrong owner','fem_cell_overlap_build',__FILE__)
      IF(kbc(j+1)-kbc(j)/=nfine)CALL oft_abort('Unexpected number of fine cells in base cell', &
        'fem_cell_overlap_build',__FILE__)
      DO m=kbc(j),kbc(j+1)-1
        e=lbc(m)
        k=k+1
        sbuf(k)=mesh%global%lc(e)
        DO ioff=1,ncp
          k=k+1
          sbuf(k)=mesh%global%lp(mesh%lc(ioff,e))
        END DO
      END DO
      scnt(q)=scnt(q)+nfine*(1+ncp)
    END DO
  END DO
  DEALLOCATE(rbuf,kbc,lbc)
  ALLOCATE(nsel(0:np-1))
  nsel=0
  DO i=1,ncand
    nsel(cown(i))=nsel(cown(i))+1
  END DO
  CALL oft_exchange_i8(scnt,sbuf,rcnt,rbuf)
  DEALLOCATE(sbuf)
  DO q=0,np-1
    IF(rcnt(q)/=nsel(q)*nfine*(1+ncp))CALL oft_abort('Unexpected halo topology reply size', &
      'fem_cell_overlap_build',__FILE__)
  END DO
  nhalo=SUM(rcnt)/(1+ncp)
END IF
CALL oft_exchange_report('cell overlap A-B (halo topology)',stimer,nsent0)
!---------------------------------------------------------------------------------
! Step C: extended cell graph and block growth
!---------------------------------------------------------------------------------
next=nc+nhalo
ALLOCATE(egid(next),epts(ncp,next),eown(next))
DO i=1,nc
  egid(i)=mesh%global%lc(i)
  DO j=1,ncp
    epts(j,i)=mesh%global%lp(mesh%lc(j,i))
  END DO
  eown(i)=me
END DO
IF(nhalo>0)THEN
  k=0
  e=nc
  DO q=0,np-1
    DO i=1,rcnt(q)/(1+ncp)
      e=e+1
      k=k+1
      egid(e)=rbuf(k)
      epts(:,e)=rbuf(k+1:k+ncp)
      k=k+ncp
      eown(e)=q
    END DO
  END DO
END IF
IF(ALLOCATED(rbuf))DEALLOCATE(rbuf)
!---Cell global ids must be unique (halo cells are never local)
ALLOCATE(lsorted(next),lperm(next))
lsorted=egid
lperm=[(i,i=1,next)]
IF(next>0)CALL sort_array(lsorted,lperm,next)
DO i=2,next
  IF(lsorted(i)==lsorted(i-1))CALL oft_abort('Duplicate cell global ids','fem_cell_overlap_build',__FILE__)
END DO
!---Point to cell CSR keyed by global point id
npair=ncp*next
ALLOCATE(pk(MAX(1,npair)),iperm(MAX(1,npair)),upos(MAX(1,npair)))
DO e=1,next
  DO j=1,ncp
    pk((e-1)*ncp+j)=epts(j,e)
    iperm((e-1)*ncp+j)=(e-1)*ncp+j
  END DO
END DO
IF(npair>0)CALL sort_array(pk,iperm,npair)
nup=0
DO i=1,npair
  IF(i==1)THEN
    nup=1
  ELSE IF(pk(i)/=pk(i-1))THEN
    nup=nup+1
  END IF
  upos(iperm(i))=nup
END DO
ALLOCATE(kp(nup+1),lpx(MAX(1,npair)))
kp=0
DO i=1,npair
  kp(upos(i)+1)=kp(upos(i)+1)+1
END DO
kp(1)=1
DO i=1,nup
  kp(i+1)=kp(i+1)+kp(i)
END DO
ALLOCATE(cnt(nup))
cnt=kp(1:nup)
DO i=1,npair
  lpx(cnt(upos(i)))=(i-1)/ncp+1
  cnt(upos(i))=cnt(upos(i))+1
END DO
DEALLOCATE(cnt,pk,iperm)
!---Cell adjacency: cells sharing at least thr points (two passes: count, fill)
ALLOCATE(nnb(next),knb(next+1))
DO pass=1,2
  !$omp parallel private(e,e2,j,k,i,m,ntl,cnt,tl)
  ALLOCATE(cnt(MAX(1,next)),tl(MAX(1,next)))
  cnt=0
  !$omp do schedule(static)
  DO e=1,next
    ntl=0
    DO j=1,ncp
      k=upos((e-1)*ncp+j)
      DO i=kp(k),kp(k+1)-1
        e2=lpx(i)
        IF(e2==e)CYCLE
        IF(cnt(e2)==0)THEN
          ntl=ntl+1
          tl(ntl)=e2
        END IF
        cnt(e2)=cnt(e2)+1
      END DO
    END DO
    m=0
    DO i=1,ntl
      IF(cnt(tl(i))>=thr)THEN
        m=m+1
        IF(pass==2)lnb(knb(e)+m-1)=tl(i)
      END IF
      cnt(tl(i))=0
    END DO
    IF(pass==1)nnb(e)=m
  END DO
  !$omp end do
  DEALLOCATE(cnt,tl)
  !$omp end parallel
  IF(pass==1)THEN
    knb(1)=1
    DO e=1,next
      knb(e+1)=knb(e)+nnb(e)
    END DO
    ALLOCATE(lnb(MAX(1,knb(next+1)-1)))
  END IF
END DO
DEALLOCATE(nnb,kp,lpx,upos)
!---Block cores: cells of the part plus cells holding a DoF of the block's color
ALLOCATE(core(nb),inb(nb))
DO pass=1,2
  DO b=1,nb
    IF(pass==2)ALLOCATE(core(b)%v(core(b)%n))
    core(b)%n=0
  END DO
  DO i=1,nc
    inb=.FALSE.
    inb(self%cpart(i))=.TRUE.
    DO j=1,nce_tot
      inb(self%dcolor(lrow(j,i)))=.TRUE.
    END DO
    DO b=1,nb
      IF(inb(b))THEN
        core(b)%n=core(b)%n+1
        IF(pass==2)core(b)%v(core(b)%n)=i
      END IF
    END DO
  END DO
END DO
DEALLOCATE(inb)
!---Grow each block by layers of adjacent cells
ALLOCATE(bx(nb))
!$omp parallel private(b,i,j,k,l,e,e2,m,istart,iend,mark,list)
ALLOCATE(mark(MAX(1,next)),list(MAX(1,next)))
mark=.FALSE.
!$omp do schedule(dynamic,1)
DO b=1,nb
  m=0
  DO i=1,core(b)%n
    e=core(b)%v(i)
    mark(e)=.TRUE.
    m=m+1
    list(m)=e
  END DO
  istart=1
  iend=m
  DO l=1,layers
    DO i=istart,iend
      e=list(i)
      DO k=knb(e),knb(e+1)-1
        e2=lnb(k)
        IF(.NOT.mark(e2))THEN
          mark(e2)=.TRUE.
          m=m+1
          list(m)=e2
        END IF
      END DO
    END DO
    IF(m==iend)EXIT
    istart=iend+1
    iend=m
  END DO
  bx(b)%n=m
  ALLOCATE(bx(b)%v(m))
  bx(b)%v=list(1:m)
  DO i=1,m
    mark(list(i))=.FALSE.
  END DO
  IF(m>0)CALL sort_array(bx(b)%v,m)
END DO
!$omp end do
DEALLOCATE(mark,list)
!$omp end parallel
DEALLOCATE(knb,lnb)
!---Halo cells used by any block
ALLOCATE(needed(next))
needed=.FALSE.
DO b=1,nb
  DO i=1,bx(b)%n
    needed(bx(b)%v(i))=.TRUE.
  END DO
END DO
!---Diagnostic cell sets (global cell ids)
ALLOCATE(self%core_cells(nb),self%blk_cells(nb))
DO b=1,nb
  self%core_cells(b)%n=core(b)%n
  ALLOCATE(self%core_cells(b)%v(core(b)%n))
  DO i=1,core(b)%n
    self%core_cells(b)%v(i)=egid(core(b)%v(i))
  END DO
  n8=INT(core(b)%n,i8)
  CALL sort_array(self%core_cells(b)%v,n8)
  self%blk_cells(b)%n=bx(b)%n
  ALLOCATE(self%blk_cells(b)%v(bx(b)%n))
  DO i=1,bx(b)%n
    self%blk_cells(b)%v(i)=egid(bx(b)%v(i))
  END DO
  n8=INT(bx(b)%n,i8)
  CALL sort_array(self%blk_cells(b)%v,n8)
END DO
CALL oft_exchange_report('cell overlap C (layers)',stimer,oft_exchange_nsent)
nsent0=oft_exchange_nsent
!---------------------------------------------------------------------------------
! Step D: DoFs (global id and owner) of the used halo cells from their owners
!---------------------------------------------------------------------------------
self%nhalo=nhalo
self%nhalo_used=COUNT(needed(nc+1:next))
ALLOCATE(hdof(nce_tot,MAX(1,nhalo)),hdown(nce_tot,MAX(1,nhalo)))
IF(.NOT.mesh%fullmesh)THEN
  !---Requests in increasing halo index (grouped by owner)
  scnt=0
  nsel=0
  ALLOCATE(sbuf(MAX(1,self%nhalo_used)))
  k=0
  DO e=nc+1,next
    IF(.NOT.needed(e))CYCLE
    k=k+1
    sbuf(k)=egid(e)
    scnt(eown(e))=scnt(eown(e))+1
    nsel(eown(e))=nsel(eown(e))+1
  END DO
  CALL oft_exchange_i8(scnt,sbuf,rcnt,rbuf)
  DEALLOCATE(sbuf)
  !---Reply with [global row id, owner] for every DoF of each requested local cell
  m=SUM(rcnt)
  ALLOCATE(sbuf(MAX(1,m*2*nce_tot)))
  k=0
  l=0
  n8=INT(next,i8)
  DO q=0,np-1
    scnt(q)=0
    DO i=1,rcnt(q)
      l=l+1
      j=INT(search_array(rbuf(l),lsorted,n8),i4)
      IF(j==0)CALL oft_abort('Requested cell not found','fem_cell_overlap_build',__FILE__)
      e=lperm(j)
      IF(e>nc)CALL oft_abort('Requested cell is not local','fem_cell_overlap_build',__FILE__)
      DO j=1,nce_tot
        row=lrow(j,e)
        sbuf(k+1)=grow(row)
        sbuf(k+2)=INT(owner(row),i8)
        k=k+2
      END DO
      scnt(q)=scnt(q)+2*nce_tot
    END DO
  END DO
  DEALLOCATE(rbuf)
  CALL oft_exchange_i8(scnt,sbuf,rcnt,rbuf)
  DEALLOCATE(sbuf)
  DO q=0,np-1
    IF(rcnt(q)/=nsel(q)*2*nce_tot)CALL oft_abort('Unexpected halo DoF reply size', &
      'fem_cell_overlap_build',__FILE__)
  END DO
  !---Replies arrive grouped by owner, in request order
  k=0
  DO e=nc+1,next
    IF(.NOT.needed(e))CYCLE
    DO j=1,nce_tot
      hdof(j,e-nc)=rbuf(k+1)
      hdown(j,e-nc)=INT(rbuf(k+2),i4)
      IF(hdown(j,e-nc)<0.OR.hdown(j,e-nc)>=np)CALL oft_abort('Invalid halo DoF owner', &
        'fem_cell_overlap_build',__FILE__)
      k=k+2
    END DO
  END DO
  DEALLOCATE(rbuf,nsel)
END IF
DEALLOCATE(lsorted,lperm)
!---------------------------------------------------------------------------------
! Step E: union, block and assigned DoF positions
!---------------------------------------------------------------------------------
nused=COUNT(needed)
ALLOCATE(dk(MAX(1,nused*nce_tot)),dtmp(MAX(1,nused*nce_tot)),iperm(MAX(1,nused*nce_tot)))
k=0
DO e=1,next
  IF(.NOT.needed(e))CYCLE
  DO j=1,nce_tot
    k=k+1
    IF(e<=nc)THEN
      dk(k)=grow(lrow(j,e))
      dtmp(k)=owner(lrow(j,e))
    ELSE
      dk(k)=hdof(j,e-nc)
      dtmp(k)=hdown(j,e-nc)
    END IF
    iperm(k)=k
  END DO
END DO
IF(k>0)CALL sort_array(dk,iperm,k)
ALLOCATE(lperm(MAX(1,k)))
nu=0
DO i=1,k
  IF(nu>0)THEN
    IF(dk(i)==dk(nu))THEN
      IF(dtmp(iperm(i))/=lperm(nu))CALL oft_abort('Inconsistent DoF owners', &
        'fem_cell_overlap_build',__FILE__)
      CYCLE
    END IF
  END IF
  nu=nu+1
  dk(nu)=dk(i)
  lperm(nu)=dtmp(iperm(i))
END DO
dofs%n=nu
dofs%nblocks=nb
ALLOCATE(dofs%gdof(nu),dofs%odof(nu))
dofs%gdof=dk(1:nu)
dofs%odof=lperm(1:nu)
DEALLOCATE(lperm)
DEALLOCATE(dk,dtmp,iperm)
!---Block DoF positions
ALLOCATE(dofs%blk(nb),dofs%asg(nb))
n8=INT(nu,i8)
!$omp parallel do private(b,i,j,e,m,l,dtmp) schedule(dynamic,1)
DO b=1,nb
  ALLOCATE(dtmp(MAX(1,bx(b)%n*nce_tot)))
  m=0
  DO i=1,bx(b)%n
    e=bx(b)%v(i)
    DO j=1,nce_tot
      m=m+1
      IF(e<=nc)THEN
        dtmp(m)=INT(search_array(grow(lrow(j,e)),dofs%gdof,n8),i4)
      ELSE
        dtmp(m)=INT(search_array(hdof(j,e-nc),dofs%gdof,n8),i4)
      END IF
    END DO
  END DO
  IF(m>0)CALL sort_array(dtmp,m)
  l=0
  DO i=1,m
    IF(l>0)THEN
      IF(dtmp(i)==dtmp(l))CYCLE
    END IF
    l=l+1
    dtmp(l)=dtmp(i)
  END DO
  dofs%blk(b)%n=l
  ALLOCATE(dofs%blk(b)%v(l))
  dofs%blk(b)%v=dtmp(1:l)
  DEALLOCATE(dtmp)
END DO
IF(nu>0)THEN
  DO b=1,nb
    IF(dofs%blk(b)%n>0)THEN
      IF(dofs%blk(b)%v(1)<1)CALL oft_abort('Block DoF missing from union','fem_cell_overlap_build',__FILE__)
    END IF
  END DO
END IF
!---Assigned DoFs: owned DoFs by color (every owned row exactly once)
ALLOCATE(hit(MAX(1,nu)),dtmp(MAX(1,mat%nr)))
hit=0
DO b=1,nb
  dofs%asg(b)%n=0
END DO
nown=0
DO i=1,mat%nr
  dtmp(i)=0
  IF(owner(i)/=me)CYCLE
  nown=nown+1
  j=INT(search_array(grow(i),dofs%gdof,n8),i4)
  IF(j==0)CALL oft_abort('Owned DoF missing from union','fem_cell_overlap_build',__FILE__)
  IF(dofs%odof(j)/=me)CALL oft_abort('Inconsistent owner of owned DoF','fem_cell_overlap_build',__FILE__)
  hit(j)=hit(j)+1
  dtmp(i)=j
  b=self%dcolor(i)
  dofs%asg(b)%n=dofs%asg(b)%n+1
END DO
DO j=1,nu
  IF((dofs%odof(j)==me).NEQV.(hit(j)==1))CALL oft_abort('Owned DoF not assigned exactly once', &
    'fem_cell_overlap_build',__FILE__)
END DO
DO b=1,nb
  ALLOCATE(dofs%asg(b)%v(dofs%asg(b)%n))
  dofs%asg(b)%n=0
END DO
DO i=1,mat%nr
  IF(dtmp(i)==0)CYCLE
  b=self%dcolor(i)
  dofs%asg(b)%n=dofs%asg(b)%n+1
  dofs%asg(b)%v(dofs%asg(b)%n)=dtmp(i)
END DO
DO b=1,nb
  IF(dofs%asg(b)%n>0)CALL sort_array(dofs%asg(b)%v,dofs%asg(b)%n)
  DO i=1,dofs%asg(b)%n
    IF(search_array(dofs%asg(b)%v(i),dofs%blk(b)%v,dofs%blk(b)%n)==0) &
      CALL oft_abort('Assigned DoF outside its block','fem_cell_overlap_build',__FILE__)
  END DO
END DO
IF(SUM(dofs%asg(:)%n)/=nown)CALL oft_abort('Owned DoF count mismatch','fem_cell_overlap_build',__FILE__)
DEALLOCATE(hit,dtmp)
CALL oft_exchange_report('cell overlap D-E (halo DoFs)',stimer,nsent0)
!---Clean up
DO b=1,nb
  IF(ALLOCATED(core(b)%v))DEALLOCATE(core(b)%v)
  IF(ALLOCATED(bx(b)%v))DEALLOCATE(bx(b)%v)
END DO
DEALLOCATE(core,bx,needed,egid,epts,eown,hdof,hdown,grow,owner,lrow,scnt,rcnt)
IF(ALLOCATED(cand))DEALLOCATE(cand,cown)
IF(ALLOCATED(cb))DEALLOCATE(cb)
IF(ALLOCATED(inv))DEALLOCATE(inv)
DEBUG_STACK_POP
END SUBROUTINE fem_cell_overlap_build
END MODULE fem_cell_overlap
