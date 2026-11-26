!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_cell_overlap.F90
!
!> Regression tests for the cell-adjacency overlap builder (`fem_cell_overlap`) on a
!! cube, with a second order Lagrange representation (single FE) and its vector
!! (composite, 3 fields) version.
!!
!! Every build is checked exactly against a brute force reference computed from all
!! cells of the mesh gathered on every process:
!! - Block cells equal the core cells grown by k layers on the global cell graph, and
!!   the cores follow their definition (cell part, plus cells holding a DoF of the
!!   block's color)
!! - Block DoFs are the DoFs of the block cells, the union is the union of the blocks
!!   and the owners match the owners of the row space
!! - Every owned DoF is assigned to exactly one block, within the block, by color
!! - Face <= edge <= point and layer k <= k+1 (as DoF sets), and with no overlap and
!!   one block the union is exactly the local rows
!! - A larger base mesh bound (`base_layers`) gives the same result
!!
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_cell_overlap
USE oft_base
USE oft_sort, ONLY: sort_array, search_array
USE oft_mesh_cube, ONLY: mesh_cube_id
USE oft_mesh_type, ONLY: oft_mesh
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE oft_stitching, ONLY: oft_global_stitch
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_overlap_dofs, overlap_dofs_delete, &
  OFT_CELL_ADJ_POINT, OFT_CELL_ADJ_EDGE, OFT_CELL_ADJ_FACE
USE oft_native_la, ONLY: oft_native_vector, native_vector_cast
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE oft_lag_basis, ONLY: oft_lag_setup
USE fem_cell_overlap, ONLY: oft_fem_cell_overlap, fem_cell_overlap_create
IMPLICIT NONE
TYPE :: set8
  INTEGER(i4) :: n = 0
  INTEGER(i8), ALLOCATABLE, DIMENSION(:) :: v
END TYPE set8
INTEGER(i4) :: io_unit,nfail,nchecks,ncp,ncg,nce_tot,me
TYPE(multigrid_mesh), TARGET :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_lag
TYPE(oft_ml_fem_comp_type), TARGET :: ML_vlag
CLASS(oft_mesh), POINTER :: mesh
CLASS(oft_vector), POINTER :: u => NULL()
CLASS(oft_matrix), POINTER :: mat => NULL()
!---Global reference data (all cells, sorted by global cell id)
INTEGER(i8), ALLOCATABLE, DIMENSION(:) :: gcell
INTEGER(i8), ALLOCATABLE, DIMENSION(:,:) :: gdofs
INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: gdown
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: knbr,lnbr,cshare
!---Local row data
INTEGER(i8), ALLOCATABLE, DIMENSION(:) :: grow,gsort
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: owner,gperm
INTEGER(i4), ALLOCATABLE, DIMENSION(:,:) :: lrow
CALL oft_init
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
CALL oft_lag_setup(mg_mesh,2,ML_lag,ML_vlag_obj=ML_vlag,minlev=-1)
oft_env%pm=.FALSE.
me=oft_env%rank
mesh=>mg_mesh%mesh
ncp=mesh%cell_np
nfail=0
nchecks=0
CALL build_cell_graph()
!---Single FE (scalar Lagrange)
CALL ML_lag%current_level%mat_create(mat)
CALL ML_lag%current_level%vec_create(u)
CALL setup_rows()
CALL run_cases(.FALSE.,1)
CALL run_cases(.FALSE.,3)
CALL mat%delete
DEALLOCATE(mat)
CALL u%delete
DEALLOCATE(u)
!---Composite FE (vector Lagrange, 3 fields)
CALL ML_vlag%current_level%mat_create(mat)
CALL ML_vlag%current_level%vec_create(u)
CALL setup_rows()
CALL run_cases(.TRUE.,1)
CALL run_cases(.TRUE.,3)
IF(oft_env%head_proc)THEN
  WRITE(*,'(A,I8,A,I8)')'Checks: ',nchecks,' failed: ',nfail
  OPEN(NEWUNIT=io_unit,FILE='cell_overlap.results')
  WRITE(io_unit,'(I0,1X,I0)')nfail,nchecks
  CLOSE(io_unit)
END IF
CALL oft_finalize
CONTAINS
! Global cell graph: all cells with their global point ids, neighbours sharing >=1 point
! (with the number of shared points)
SUBROUTINE build_cell_graph()
INTEGER(i4) :: i,j,k,l,c,c2,ierr,npair,nup,nn
INTEGER(i4), ALLOCATABLE :: counts(:),displs(:),iperm(:),upos(:),kp(:),lp(:),cnt(:),tl(:),cperm(:)
INTEGER(i8), ALLOCATABLE :: lbuf(:),gbuf(:),pk(:),gpts(:,:),gc(:)
ALLOCATE(lbuf((1+ncp)*mesh%nc))
DO i=1,mesh%nc
  lbuf((i-1)*(1+ncp)+1)=mesh%global%lc(i)
  DO j=1,ncp
    lbuf((i-1)*(1+ncp)+1+j)=mesh%global%lp(mesh%lc(j,i))
  END DO
END DO
CALL allgather8(lbuf,gbuf)
ncg=SIZE(gbuf)/(1+ncp)
IF(INT(ncg,i8)/=mesh%global%nc)CALL oft_abort('Cell count mismatch','build_cell_graph',__FILE__)
ALLOCATE(gc(ncg),cperm(ncg),gpts(ncp,ncg),gcell(ncg))
DO i=1,ncg
  gc(i)=gbuf((i-1)*(1+ncp)+1)
  cperm(i)=i
END DO
CALL sort_array(gc,cperm,ncg)
DO i=1,ncg
  gcell(i)=gc(i)
  gpts(:,i)=gbuf((cperm(i)-1)*(1+ncp)+2:cperm(i)*(1+ncp))
  IF(i>1)THEN
    IF(gcell(i)==gcell(i-1))CALL oft_abort('Duplicate cell','build_cell_graph',__FILE__)
  END IF
END DO
!---Point to cell map
npair=ncp*ncg
ALLOCATE(pk(npair),iperm(npair),upos(npair))
DO i=1,npair
  pk(i)=gpts(MOD(i-1,ncp)+1,(i-1)/ncp+1)
  iperm(i)=i
END DO
CALL sort_array(pk,iperm,npair)
nup=0
DO i=1,npair
  IF(i==1)THEN
    nup=1
  ELSE IF(pk(i)/=pk(i-1))THEN
    nup=nup+1
  END IF
  upos(iperm(i))=nup
END DO
ALLOCATE(kp(nup+1),lp(npair),cnt(MAX(nup,ncg)),tl(ncg))
kp=0
DO i=1,npair
  kp(upos(i)+1)=kp(upos(i)+1)+1
END DO
kp(1)=1
DO i=1,nup
  kp(i+1)=kp(i+1)+kp(i)
END DO
cnt(1:nup)=kp(1:nup)
DO i=1,npair
  lp(cnt(upos(i)))=(i-1)/ncp+1
  cnt(upos(i))=cnt(upos(i))+1
END DO
!---Neighbours (brute force with counts of shared points)
ALLOCATE(knbr(ncg+1))
cnt=0
knbr(1)=1
DO l=1,2
  DO c=1,ncg
    nn=0
    DO j=1,ncp
      k=upos((c-1)*ncp+j)
      DO i=kp(k),kp(k+1)-1
        c2=lp(i)
        IF(c2==c)CYCLE
        IF(cnt(c2)==0)THEN
          nn=nn+1
          tl(nn)=c2
        END IF
        cnt(c2)=cnt(c2)+1
      END DO
    END DO
    DO i=1,nn
      IF(l==2)THEN
        lnbr(knbr(c)+i-1)=tl(i)
        cshare(knbr(c)+i-1)=cnt(tl(i))
      END IF
      cnt(tl(i))=0
    END DO
    IF(l==1)knbr(c+1)=knbr(c)+nn
  END DO
  IF(l==1)ALLOCATE(lnbr(knbr(ncg+1)-1),cshare(knbr(ncg+1)-1))
END DO
DEALLOCATE(lbuf,gbuf,gc,cperm,gpts,pk,iperm,upos,kp,lp,cnt,tl)
IF(oft_env%head_proc)WRITE(*,'(A,I8,A,I8)')'Global cells: ',ncg,' neighbour pairs: ',knbr(ncg+1)-1
END SUBROUTINE build_cell_graph
! Local row data (global id, owner) and global DoFs of every cell
SUBROUTINE setup_rows()
INTEGER(i4) :: i,j,f,ioff,nce
INTEGER(i4), ALLOCATABLE :: ldofs(:),cperm(:)
INTEGER(i8), ALLOCATABLE :: lbuf(:),gbuf(:),gc(:)
REAL(r8), ALLOCATABLE :: rown(:)
CLASS(oft_native_vector), POINTER :: uv
IF(ALLOCATED(grow))DEALLOCATE(grow,gsort,gperm,owner,lrow,gdofs,gdown)
ALLOCATE(grow(mat%nr),gsort(mat%nr),gperm(mat%nr),owner(mat%nr),rown(mat%nr))
DO f=1,mat%ni
  DO j=1,mat%i_map(f)%n
    grow(mat%i_map(f)%offset+j)=mat%i_map(f)%offsetg+mat%i_map(f)%lge(j)
  END DO
END DO
gsort=grow
gperm=[(i,i=1,mat%nr)]
CALL sort_array(gsort,gperm,mat%nr)
IF(.NOT.native_vector_cast(uv,u))CALL oft_abort('Native vector required','setup_rows',__FILE__)
rown=REAL(me,r8)
CALL oft_global_stitch(uv%stitch_info,rown,0)
owner=NINT(rown)
!---Local rows of each cell (ncdofs order, all fields)
nce_tot=0
DO f=1,mat%ni
  IF(mat%ni==1)THEN
    nce_tot=nce_tot+ML_lag%current_level%nce
  ELSE
    nce_tot=nce_tot+ML_vlag%current_level%fields(f)%fe%nce
  END IF
END DO
ALLOCATE(lrow(nce_tot,mesh%nc))
ioff=0
DO f=1,mat%ni
  IF(mat%ni==1)THEN
    nce=ML_lag%current_level%nce
  ELSE
    nce=ML_vlag%current_level%fields(f)%fe%nce
  END IF
  ALLOCATE(ldofs(nce))
  DO i=1,mesh%nc
    IF(mat%ni==1)THEN
      CALL ML_lag%current_level%ncdofs(i,ldofs)
    ELSE
      CALL ML_vlag%current_level%fields(f)%fe%ncdofs(i,ldofs)
    END IF
    lrow(ioff+1:ioff+nce,i)=mat%i_map(f)%offset+ldofs
  END DO
  ioff=ioff+nce
  DEALLOCATE(ldofs)
END DO
!---Gather [cell gid, (dof gid, owner) x nce_tot] of all cells
ALLOCATE(lbuf((1+2*nce_tot)*mesh%nc))
DO i=1,mesh%nc
  ioff=(i-1)*(1+2*nce_tot)
  lbuf(ioff+1)=mesh%global%lc(i)
  DO j=1,nce_tot
    lbuf(ioff+2*j)=grow(lrow(j,i))
    lbuf(ioff+2*j+1)=INT(owner(lrow(j,i)),i8)
  END DO
END DO
CALL allgather8(lbuf,gbuf)
ALLOCATE(gdofs(nce_tot,ncg),gdown(nce_tot,ncg),gc(ncg),cperm(ncg))
DO i=1,ncg
  gc(i)=gbuf((i-1)*(1+2*nce_tot)+1)
  cperm(i)=i
END DO
CALL sort_array(gc,cperm,ncg)
DO i=1,ncg
  IF(gc(i)/=gcell(i))CALL oft_abort('Cell list mismatch','setup_rows',__FILE__)
  ioff=(cperm(i)-1)*(1+2*nce_tot)
  DO j=1,nce_tot
    gdofs(j,i)=gbuf(ioff+2*j)
    gdown(j,i)=INT(gbuf(ioff+2*j+1),i4)
  END DO
END DO
DEALLOCATE(lbuf,gbuf,gc,cperm,rown)
END SUBROUTINE setup_rows
! All builds for one FE representation and number of blocks
SUBROUTINE run_cases(comp,nparts)
LOGICAL, INTENT(in) :: comp
INTEGER(i4), INTENT(in) :: nparts
TYPE(oft_fem_cell_overlap) :: ovl
TYPE(oft_overlap_dofs) :: dofs
TYPE(set8), ALLOCATABLE :: res(:,:,:)
INTEGER(i4) :: lay,ia,b,i,nb,nh,nhu,nbmax
INTEGER(i4), ALLOCATABLE :: dcolor(:)
CHARACTER(LEN=64) :: tag
ALLOCATE(dcolor(mat%nr))
IF(comp)THEN
  CALL fem_cell_overlap_create(ovl,mg_mesh,ML_vlag%current_level,nparts,dcolor)
ELSE
  CALL fem_cell_overlap_create(ovl,mg_mesh,ML_lag%current_level,nparts,dcolor)
END IF
nb=ovl%nblocks()
CALL check('nblocks',nb==nparts)
CALL check('dcolor output',ALL(dcolor==ovl%dcolor))
ALLOCATE(res(nb,0:2,3))
DO lay=0,2
  DO ia=1,3
    WRITE(tag,'(A,L1,A,I1,A,I1,A,I1)')'comp=',comp,' nparts=',nparts,' layers=',lay,' adj=',ia
    CALL ovl%build(mat,u,lay,ia,-1,dofs)
    CALL check_build(ovl,dofs,lay,ia,TRIM(tag))
    nh=oft_mpi_sum(ovl%nhalo)
    nhu=oft_mpi_sum(ovl%nhalo_used)
    nbmax=oft_mpi_max(MAXVAL(ovl%blk_cells(:)%n))
    IF(oft_env%head_proc)WRITE(*,'(2A,4(A,I8))')TRIM(tag),':',' halo=',nh,' used=',nhu, &
      ' max block cells=',nbmax,' of',ncg
    DO b=1,nb
      res(b,lay,ia)%n=dofs%blk(b)%n
      ALLOCATE(res(b,lay,ia)%v(dofs%blk(b)%n))
      DO i=1,dofs%blk(b)%n
        res(b,lay,ia)%v(i)=dofs%gdof(dofs%blk(b)%v(i))
      END DO
    END DO
    IF(lay==0.AND.nparts==1)THEN
      CALL check('Union = local rows (layers=0, nparts=1) '//TRIM(tag), &
        same8(dofs%gdof,dofs%n,gsort,mat%nr))
    END IF
  END DO
END DO
!---Monotonicity: face <= edge <= point, layer k <= k+1
DO b=1,nb
  DO lay=0,2
    CALL check('face <= edge',subset8(res(b,lay,3),res(b,lay,2)))
    CALL check('edge <= point',subset8(res(b,lay,2),res(b,lay,1)))
    IF(lay<2)THEN
      DO ia=1,3
        CALL check('layer k <= k+1',subset8(res(b,lay,ia),res(b,lay+1,ia)))
      END DO
    END IF
  END DO
END DO
!---Larger base mesh bound gives the same blocks
WRITE(tag,'(A,L1,A,I1,A)')'comp=',comp,' nparts=',nparts,' layers=1 adj=1 base_layers=3'
CALL ovl%build(mat,u,1,OFT_CELL_ADJ_POINT,3,dofs)
CALL check_build(ovl,dofs,1,OFT_CELL_ADJ_POINT,TRIM(tag))
DO b=1,nb
  CALL check('base_layers>layers unchanged '//TRIM(tag), &
    same8([(dofs%gdof(dofs%blk(b)%v(i)),i=1,dofs%blk(b)%n)],dofs%blk(b)%n,res(b,1,1)%v,res(b,1,1)%n))
END DO
CALL overlap_dofs_delete(dofs)
CALL ovl%delete
CALL check('delete',(.NOT.ALLOCATED(ovl%cpart)).AND.(.NOT.ALLOCATED(ovl%blk_cells)))
DEALLOCATE(res,dcolor)
END SUBROUTINE run_cases
! Exact checks of one build against the brute force reference
SUBROUTINE check_build(ovl,dofs,layers,adj,tag)
TYPE(oft_fem_cell_overlap), INTENT(in) :: ovl
TYPE(oft_overlap_dofs), INTENT(in) :: dofs
INTEGER(i4), INTENT(in) :: layers,adj
CHARACTER(LEN=*), INTENT(in) :: tag
INTEGER(i4) :: b,i,j,k,l,c,thr,m,nf,nown,nn,p,r
INTEGER(i4), ALLOCATABLE :: lvl(:),list(:),hit(:)
INTEGER(i8), ALLOCATABLE :: cref(:),dref(:),un(:),dk(:)
LOGICAL :: ok,inb
SELECT CASE(adj)
CASE(OFT_CELL_ADJ_POINT)
  thr=1
CASE(OFT_CELL_ADJ_EDGE)
  thr=2
CASE DEFAULT
  thr=mesh%face_np
END SELECT
CALL check('nblocks '//tag,dofs%nblocks==ovl%nparts.AND.SIZE(dofs%blk)==ovl%nparts.AND.SIZE(dofs%asg)==ovl%nparts)
!---Union is sorted unique
ok=.TRUE.
DO i=2,dofs%n
  IF(dofs%gdof(i)<=dofs%gdof(i-1))ok=.FALSE.
END DO
CALL check('union sorted unique '//tag,ok)
ALLOCATE(lvl(ncg),list(ncg))
lvl=-1
ALLOCATE(un(0))
DO b=1,ovl%nparts
  !---Core definition
  nn=0
  DO c=1,mesh%nc
    inb=(ovl%cpart(c)==b)
    DO j=1,nce_tot
      IF(ovl%dcolor(lrow(j,c))==b)inb=.TRUE.
    END DO
    IF(inb)nn=nn+1
  END DO
  ALLOCATE(cref(nn))
  nn=0
  DO c=1,mesh%nc
    inb=(ovl%cpart(c)==b)
    DO j=1,nce_tot
      IF(ovl%dcolor(lrow(j,c))==b)inb=.TRUE.
    END DO
    IF(inb)THEN
      nn=nn+1
      cref(nn)=mesh%global%lc(c)
    END IF
  END DO
  IF(nn>0)CALL sort_array(cref,INT(nn,i8))
  CALL check('core cells '//tag,same8(ovl%core_cells(b)%v,ovl%core_cells(b)%n,cref,nn))
  IF(ovl%nparts==1)CALL check('core = all cells '//tag,nn==mesh%nc)
  !---Reference BFS on the global cell graph
  m=0
  DO i=1,nn
    c=INT(search_array(cref(i),gcell,INT(ncg,i8)),i4)
    lvl(c)=0
    m=m+1
    list(m)=c
  END DO
  k=1
  DO l=1,layers
    nf=m
    DO i=k,nf
      c=list(i)
      DO j=knbr(c),knbr(c+1)-1
        IF(cshare(j)<thr)CYCLE
        IF(lvl(lnbr(j))<0)THEN
          lvl(lnbr(j))=l
          m=m+1
          list(m)=lnbr(j)
        END IF
      END DO
    END DO
    k=nf+1
  END DO
  DEALLOCATE(cref)
  ALLOCATE(cref(m),dref(m*nce_tot))
  DO i=1,m
    cref(i)=gcell(list(i))
    lvl(list(i))=-1
    DO j=1,nce_tot
      dref((i-1)*nce_tot+j)=gdofs(j,list(i))
    END DO
  END DO
  CALL sort_array(cref,INT(m,i8))
  CALL check('block cells '//tag,same8(ovl%blk_cells(b)%v,ovl%blk_cells(b)%n,cref,m))
  !---Block DoFs = DoFs of the block cells
  nn=unique8(dref,m*nce_tot)
  ALLOCATE(dk(dofs%blk(b)%n))
  DO i=1,dofs%blk(b)%n
    dk(i)=dofs%gdof(dofs%blk(b)%v(i))
  END DO
  ok=.TRUE.
  DO i=2,dofs%blk(b)%n
    IF(dofs%blk(b)%v(i)<=dofs%blk(b)%v(i-1))ok=.FALSE.
  END DO
  CALL check('block positions sorted '//tag,ok)
  CALL check('block DoFs '//tag,same8(dk,dofs%blk(b)%n,dref,nn))
  un=[un,dref(1:nn)]
  DEALLOCATE(cref,dref,dk)
END DO
!---Union = union of blocks
nn=unique8(un,SIZE(un))
CALL check('union of blocks '//tag,same8(dofs%gdof,dofs%n,un,nn))
DEALLOCATE(un)
!---Owners against the reference owners of all cells' DoFs
ok=.TRUE.
ALLOCATE(hit(dofs%n))
hit=-1
DO c=1,ncg
  DO j=1,nce_tot
    p=INT(search_array(gdofs(j,c),dofs%gdof,INT(dofs%n,i8)),i4)
    IF(p==0)CYCLE
    IF(hit(p)>=0.AND.hit(p)/=gdown(j,c))ok=.FALSE.
    hit(p)=gdown(j,c)
  END DO
END DO
DO p=1,dofs%n
  IF(hit(p)/=dofs%odof(p))ok=.FALSE.
END DO
CALL check('owners '//tag,ok)
!---Assigned DoFs: every owned row exactly once, within its block, by color
hit=0
ok=.TRUE.
DO b=1,ovl%nparts
  DO i=1,dofs%asg(b)%n
    p=dofs%asg(b)%v(i)
    IF(i>1)THEN
      IF(p<=dofs%asg(b)%v(i-1))ok=.FALSE.
    END IF
    hit(p)=hit(p)+1
    IF(search_array(p,dofs%blk(b)%v,dofs%blk(b)%n)==0)ok=.FALSE.
    IF(dofs%odof(p)/=me)ok=.FALSE.
    r=INT(search_array(dofs%gdof(p),gsort,INT(mat%nr,i8)),i4)
    IF(r==0)THEN
      ok=.FALSE.
    ELSE
      IF(ovl%dcolor(gperm(r))/=b)ok=.FALSE.
    END IF
  END DO
END DO
CALL check('assigned within block, owned, by color '//tag,ok)
nown=COUNT(owner==me)
ok=(SUM(hit)==nown)
DO i=1,mat%nr
  IF(owner(i)/=me)CYCLE
  p=INT(search_array(grow(i),dofs%gdof,INT(dofs%n,i8)),i4)
  IF(p==0)THEN
    ok=.FALSE.
  ELSE
    IF(hit(p)/=1)ok=.FALSE.
  END IF
END DO
CALL check('owned rows assigned once '//tag,ok)
DEALLOCATE(hit,lvl,list)
END SUBROUTINE check_build
! Sort and remove duplicates in place, returning the number of unique values
FUNCTION unique8(a,n) RESULT(nu)
INTEGER(i8), INTENT(inout) :: a(:)
INTEGER(i4), INTENT(in) :: n
INTEGER(i4) :: nu,i
nu=0
IF(n==0)RETURN
CALL sort_array(a,INT(n,i8))
nu=1
DO i=2,n
  IF(a(i)/=a(nu))THEN
    nu=nu+1
    a(nu)=a(i)
  END IF
END DO
END FUNCTION unique8
! Equality of two int8 lists
FUNCTION same8(a,na,b,nb) RESULT(eq)
INTEGER(i8), INTENT(in) :: a(:),b(:)
INTEGER(i4), INTENT(in) :: na,nb
LOGICAL :: eq
eq=(na==nb)
IF(eq.AND.na>0)eq=ALL(a(1:na)==b(1:nb))
END FUNCTION same8
! Is sorted set a a subset of sorted set b?
FUNCTION subset8(a,b) RESULT(sub)
TYPE(set8), INTENT(in) :: a,b
LOGICAL :: sub
INTEGER(i4) :: i
sub=.TRUE.
DO i=1,a%n
  IF(search_array(a%v(i),b%v,INT(b%n,i8))==0)THEN
    sub=.FALSE.
    RETURN
  END IF
END DO
END FUNCTION subset8
! Gather int8 lists of all processes (in rank order) on every process
SUBROUTINE allgather8(lbuf,gbuf)
INTEGER(i8), INTENT(in) :: lbuf(:)
INTEGER(i8), ALLOCATABLE, INTENT(out) :: gbuf(:)
INTEGER(i4) :: i,n,ierr
INTEGER(i4), ALLOCATABLE :: counts(:),displs(:)
n=SIZE(lbuf)
ALLOCATE(counts(oft_env%nprocs),displs(oft_env%nprocs))
#ifdef HAVE_MPI
CALL MPI_ALLGATHER(n,1,OFT_MPI_I4,counts,1,OFT_MPI_I4,oft_env%comm,ierr)
#else
counts=n
#endif
displs(1)=0
DO i=2,oft_env%nprocs
  displs(i)=displs(i-1)+counts(i-1)
END DO
ALLOCATE(gbuf(SUM(counts)))
#ifdef HAVE_MPI
CALL MPI_ALLGATHERV(lbuf,n,OFT_MPI_I8,gbuf,counts,displs,OFT_MPI_I8,oft_env%comm,ierr)
#else
gbuf=lbuf
#endif
DEALLOCATE(counts,displs)
END SUBROUTINE allgather8
! Record a check (failure on any process fails it)
SUBROUTINE check(name,passed)
CHARACTER(LEN=*), INTENT(in) :: name
LOGICAL, INTENT(in) :: passed
INTEGER(i4) :: nbad
nbad=0
IF(.NOT.passed)nbad=1
nbad=oft_mpi_sum(nbad)
nchecks=nchecks+1
IF(nbad>0)THEN
  nfail=nfail+1
  IF(oft_env%head_proc)WRITE(*,'(3A,I4,A)')'FAILED: ',name,' (on ',nbad,' processes)'
END IF
END SUBROUTINE check
END PROGRAM test_cell_overlap
