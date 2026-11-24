!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file native_overlap.F90
!
!> Overlapping subdomains across processes for native CRS matrices (Schwarz methods)
!!
!! Builds, on each process, the matrix of a subdomain made of the process' owned
!! rows plus `levels` layers of overlap grown along the global matrix graph, which
!! may include rows held only by other processes ("ghost" rows), together with scatter
!! plans moving vector values between the distributed vector and the subdomain.
!!
!! Matrices are stored unassembled per process, so rows are first made complete:
!! - The full local matrix (`Mfull`, see @ref oft_native_la::native_matrix_setup_full
!!   "native_matrix_setup_full") holds the assembled value of every entry whose row
!!   and column are both local, and these values are bitwise identical on every
!!   process holding the entry. Rows of DOFs that are not shared (interior rows) are
!!   therefore complete locally.
!! - A row shared with other processes (interface row) is incomplete everywhere: each
!!   process only has the columns it holds, and may lack entries that couple two of
!!   its DOFs through an element it does not hold. The complete row is the union of
!!   the local rows of all processes holding it (common entries agree). The owner of
!!   each DOF builds its complete row from the rows sent by the other holders.
!! - Every owned row is then complete at its owner, and the rows of all other
!!   subdomain DOFs (local copies owned elsewhere and ghosts) are requested from
!!   their owners, keyed by global index (`lge` of the row maps).
!!
!! Invariants relied on: every DOF has exactly one owner (`leo`); global indices are
!! unique per DOF (meshes with periodic images on one process are not supported);
!! the matrix pattern is fixed after setup (values are refreshed by
!! @ref native_overlap_update with fixed scatter plans); vectors are consistent across
!! copies when gathered.
!!
!! Setup exchanges index data with `MPI_ALLTOALL`/`MPI_ALLTOALLV` (collective);
!! numeric updates and vector transfers use point-to-point scatter plans.
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup doxy_oft_native_la
!---------------------------------------------------------------------------------
MODULE oft_native_overlap_mod
USE oft_local
USE oft_base
USE oft_sort, ONLY: sort_array, search_array
USE oft_stitching, ONLY: oft_seam, oft_global_stitch
USE oft_scatter
USE oft_la_base, ONLY: oft_vector
USE oft_native_la, ONLY: oft_native_matrix, oft_native_vector, oft_native_pattern, native_pattern_create, &
  native_matrix_set_pattern, native_matrix_setup_full, native_vector_cast
IMPLICIT NONE
#include "local.h"
PRIVATE
!---------------------------------------------------------------------------------
!> Growable integer list (setup work storage)
!---------------------------------------------------------------------------------
TYPE :: ilist
  INTEGER(i4) :: n = 0 !< Number of values
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: v !< Values
END TYPE ilist
!---------------------------------------------------------------------------------
!> Overlapping subdomain of a distributed native matrix
!!
!! Subdomain DOFs are numbered with the DOFs held locally first (in increasing local
!! index) followed by ghost DOFs (in increasing global index).
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: oft_native_overlap
  INTEGER(i4) :: levels = 0 !< Number of overlap layers
  INTEGER(i4) :: n = 0 !< Number of subdomain DOFs
  INTEGER(i4) :: nloc = 0 !< Number of subdomain DOFs held locally (numbered first)
  INTEGER(i4) :: nown = 0 !< Number of subdomain DOFs owned by this process
  INTEGER(i4) :: nr = 0 !< Number of local rows of the parent matrix
  INTEGER(i4) :: ncr = 0 !< Number of entries of the complete owned rows (`cvals(1:ncr)`)
  INTEGER(i8), ALLOCATABLE, DIMENSION(:) :: gsub !< Global index of each subdomain DOF [n]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: lsub !< Local index of each locally held subdomain DOF [nloc]
  INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: osub !< Owner rank of each subdomain DOF [n]
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: cvals !< Complete owned rows, then rows received from owners
  REAL(r8), ALLOCATABLE, DIMENSION(:) :: acc !< Accumulation storage for additive write-back [nr]
  TYPE(oft_scatter_plan) :: plan_rows !< Local `Mfull` values to complete owned rows (`cvals(1:ncr)`)
  TYPE(oft_scatter_plan) :: plan_ghost !< Complete owned rows to rows of non-owned subdomain DOFs
  TYPE(oft_scatter_plan) :: plan_vals !< Row values (`cvals`) to subdomain matrix values
  TYPE(oft_scatter_plan) :: plan_vec !< Distributed vector to subdomain vector (reverse: add back)
  TYPE(oft_scatter_plan) :: plan_own !< Owned subdomain DOFs (restricted write-back)
  TYPE(oft_native_matrix), POINTER :: amat => NULL() !< Subdomain matrix [n,n]
  CLASS(oft_native_matrix), POINTER :: parent => NULL() !< Parent (distributed) matrix
  TYPE(oft_seam), POINTER :: seam => NULL() !< Seam of the parent row space
END TYPE oft_native_overlap
PUBLIC native_overlap_setup, native_overlap_update, native_overlap_gather, native_overlap_scatter
PUBLIC native_overlap_delete
CONTAINS
!---------------------------------------------------------------------------------
!> Set up an overlapping subdomain of a distributed native matrix (collective)
!!
!! The subdomain is the set of DOFs owned by this process grown by `levels` layers
!! along the global graph of `mat`. Builds the subdomain matrix structure, the scatter
!! plans and the subdomain matrix values (@ref native_overlap_update).
!---------------------------------------------------------------------------------
SUBROUTINE native_overlap_setup(self,mat,u,levels)
TYPE(oft_native_overlap), INTENT(inout) :: self !< Overlap object
CLASS(oft_native_matrix), TARGET, INTENT(inout) :: mat !< Parent matrix (square)
CLASS(oft_vector), INTENT(inout) :: u !< Vector of the parent row space
INTEGER(i4), INTENT(in) :: levels !< Number of overlap layers
CLASS(oft_native_vector), POINTER :: uv
INTEGER(i4) :: i,j,k,m,q,r,np,me,nown_rows,ntmp,lev,nnz,jp,jn
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: owner,orow,olist,scnt,rcnt,ckr,co,gkr,go,ipos,itmp, &
  rrow,rown,krec,lrec,isort,dlist_own,dsrc,kr_sub,lpos,sub_of_d,newo
INTEGER(i4), POINTER, CONTIGUOUS, DIMENSION(:) :: kr,lc
INTEGER(i8), ALLOCATABLE, DIMENSION(:) :: gidx,gsorted,sbuf,rbuf,cg,gg,rg,tmp8,dlist,front,newg,dsorted
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: gperm,dperm
REAL(r8), ALLOCATABLE, DIMENSION(:) :: rown_r
TYPE(ilist), ALLOCATABLE, DIMENSION(:) :: rows_send,rows_recv,ghost_send,ghost_recv
TYPE(oft_native_pattern), POINTER :: pattern
INTEGER(i4) :: nd,nfront,nnew,ng_rows,ngv,nsub_nnz
DEBUG_STACK_PUSH
IF(.NOT.native_vector_cast(uv,u))CALL oft_abort('Native vector required','native_overlap_setup',__FILE__)
IF(mat%nr/=mat%nc)CALL oft_abort('Square matrix required','native_overlap_setup',__FILE__)
CALL native_overlap_delete(self)
self%parent=>mat
self%seam=>uv%stitch_info
self%levels=levels
self%nr=mat%nr
np=oft_env%nprocs
me=oft_env%rank
!---Global index of each local row
ALLOCATE(gidx(mat%nr))
DO i=1,mat%ni
  IF(mat%i_map(i)%per)CALL oft_abort('Periodic row maps are not supported','native_overlap_setup',__FILE__)
  DO j=1,mat%i_map(i)%n
    gidx(mat%i_map(i)%offset+j)=mat%i_map(i)%offsetg+mat%i_map(i)%lge(j)
  END DO
END DO
ALLOCATE(gsorted(mat%nr),gperm(mat%nr))
gsorted=gidx
gperm=[(i,i=1,mat%nr)]
IF(mat%nr>0)CALL sort_array(gsorted,gperm,mat%nr)
!---Owner of each local row (owner rank synchronized from the owner)
ALLOCATE(rown_r(mat%nr),owner(mat%nr))
rown_r=REAL(me,r8)
CALL oft_global_stitch(self%seam,rown_r,0)
owner=NINT(rown_r)
DEALLOCATE(rown_r)
!---Full local matrix
CALL native_matrix_setup_full(mat,u)
kr=>mat%kr
lc=>mat%lc
!---Owned rows
ALLOCATE(orow(mat%nr))
orow=0
nown_rows=0
DO i=1,mat%nr
  IF(owner(i)==me)THEN
    nown_rows=nown_rows+1
    orow(i)=nown_rows
  END IF
END DO
ALLOCATE(olist(nown_rows))
DO i=1,mat%nr
  IF(orow(i)>0)olist(orow(i))=i
END DO
!------------------------------------------------------------------------------
! Complete owned rows: other holders send their local rows to the owner
!------------------------------------------------------------------------------
ALLOCATE(scnt(0:np-1),rcnt(0:np-1),rows_send(0:np-1),rows_recv(0:np-1))
scnt=0
DO i=1,mat%nr
  IF(owner(i)==me)CYCLE
  scnt(owner(i))=scnt(owner(i))+2+2*(kr(i+1)-kr(i))
END DO
ALLOCATE(sbuf(MAX(1,SUM(scnt))),ipos(0:np-1))
CALL csr_offsets(scnt,ipos)
DO i=1,mat%nr
  q=owner(i)
  IF(q==me)CYCLE
  sbuf(ipos(q))=gidx(i)
  sbuf(ipos(q)+1)=kr(i+1)-kr(i)
  ipos(q)=ipos(q)+2
  DO j=kr(i),kr(i+1)-1
    sbuf(ipos(q))=gidx(lc(j))
    sbuf(ipos(q)+1)=owner(lc(j))
    ipos(q)=ipos(q)+2
    CALL ilist_push(rows_send(q),j) ! Value source (Mfull location) for the plan
  END DO
END DO
CALL oft_exchange_i8(scnt,sbuf,rcnt,rbuf)
DEALLOCATE(sbuf)
!---Received entries: owned row, global column and owner, in arrival order
m=0
k=1
DO WHILE(k<=SUM(rcnt))
  m=m+INT(rbuf(k+1),4)
  k=k+2+2*INT(rbuf(k+1),4)
END DO
ALLOCATE(rrow(MAX(1,m)),rg(MAX(1,m)),rown(MAX(1,m)))
m=0
k=1
DO WHILE(k<=SUM(rcnt))
  r=local_index(rbuf(k))
  IF(r==0)CALL oft_abort('Received row not held locally','native_overlap_setup',__FILE__)
  IF(orow(r)==0)CALL oft_abort('Received row not owned','native_overlap_setup',__FILE__)
  ntmp=INT(rbuf(k+1),4)
  DO j=1,ntmp
    m=m+1
    rrow(m)=orow(r)
    rg(m)=rbuf(k+2*j)
    rown(m)=INT(rbuf(k+2*j+1),4)
  END DO
  k=k+2+2*ntmp
END DO
!---Union of local and received columns of each owned row (sorted by global index)
ALLOCATE(krec(nown_rows+1))
krec=0
DO j=1,m
  krec(rrow(j)+1)=krec(rrow(j)+1)+1
END DO
krec(1)=1
DO r=1,nown_rows
  krec(r+1)=krec(r+1)+krec(r)
END DO
ALLOCATE(lrec(MAX(1,m)),itmp(nown_rows))
itmp=krec(1:nown_rows)
DO j=1,m
  lrec(itmp(rrow(j)))=j
  itmp(rrow(j))=itmp(rrow(j))+1
END DO
DEALLOCATE(itmp)
ALLOCATE(ckr(nown_rows+1))
ckr(1)=1
ALLOCATE(cg(MAX(1,kr(mat%nr+1)-1+m)),co(MAX(1,kr(mat%nr+1)-1+m)))
DO r=1,nown_rows
  i=olist(r)
  ntmp=(kr(i+1)-kr(i))+(krec(r+1)-krec(r))
  ALLOCATE(tmp8(ntmp),isort(ntmp),itmp(ntmp))
  ntmp=0
  DO j=kr(i),kr(i+1)-1
    ntmp=ntmp+1
    tmp8(ntmp)=gidx(lc(j))
    itmp(ntmp)=owner(lc(j))
  END DO
  DO k=krec(r),krec(r+1)-1
    ntmp=ntmp+1
    tmp8(ntmp)=rg(lrec(k))
    itmp(ntmp)=rown(lrec(k))
  END DO
  isort=[(k,k=1,ntmp)]
  CALL sort_array(tmp8,isort,ntmp)
  k=ckr(r)
  DO j=1,ntmp
    IF(j>1)THEN
      IF(tmp8(j)==tmp8(j-1))CYCLE
    END IF
    cg(k)=tmp8(j)
    co(k)=itmp(isort(j))
    k=k+1
  END DO
  ckr(r+1)=k
  DEALLOCATE(tmp8,isort,itmp)
END DO
self%ncr=ckr(nown_rows+1)-1
!---Plan from Mfull to complete owned rows: received entries (grouped by source
!---process in increasing rank order) and local entries
m=0
k=1
DO q=0,np-1
  ntmp=0
  DO WHILE(ntmp<rcnt(q))
    jn=INT(rbuf(k+1),4)
    DO j=1,jn
      m=m+1
      CALL ilist_push(rows_recv(q),complete_pos(rrow(m),rg(m)))
    END DO
    ntmp=ntmp+2+2*jn
    k=k+2+2*jn
  END DO
END DO
DEALLOCATE(rbuf)
ALLOCATE(lpos(MAX(1,kr(mat%nr+1)-1)),itmp(MAX(1,kr(mat%nr+1)-1)))
ntmp=0
DO r=1,nown_rows
  i=olist(r)
  DO j=kr(i),kr(i+1)-1
    ntmp=ntmp+1
    itmp(ntmp)=j
    lpos(ntmp)=complete_pos(r,gidx(lc(j)))
  END DO
END DO
CALL build_plan(self%plan_rows,rows_send,rows_recv,ntmp,itmp,lpos)
DEALLOCATE(lpos,itmp,rrow,rg,rown,krec,lrec)
!------------------------------------------------------------------------------
! Grow the subdomain from the owned rows, requesting rows of new DOFs from owners
!------------------------------------------------------------------------------
ALLOCATE(ghost_send(0:np-1),ghost_recv(0:np-1))
nd=nown_rows
ALLOCATE(dlist(MAX(1,nd)),newo(MAX(1,nd)),dsrc(MAX(1,nd)))
DO r=1,nown_rows
  dlist(r)=gidx(olist(r))
  newo(r)=me
  dsrc(r)=r
END DO
ng_rows=0
ngv=0
ALLOCATE(gkr(1),gg(1),go(1))
gkr(1)=1
nfront=nd
ALLOCATE(front(MAX(1,nfront)))
IF(nd>0)front(1:nd)=dlist(1:nd)
DO lev=1,levels
  !---Candidate DOFs: columns of the rows of the front not yet in the subdomain
  CALL sorted_copy(dlist,nd,dsorted)
  ALLOCATE(newg(MAX(1,count_front_cols())),itmp(MAX(1,count_front_cols())))
  nnew=0
  DO k=1,nfront
    j=dlist_index(front(k))
    CALL row_entries(dsrc(j))
  END DO
  !---Unique new DOFs
  IF(nnew>0)THEN
    ALLOCATE(isort(nnew),rown(nnew))
    isort=[(k,k=1,nnew)]
    rown=itmp(1:nnew)
    CALL sort_array(newg,isort,nnew)
    m=0
    DO k=1,nnew
      IF(m>0)THEN
        IF(newg(k)==newg(m))CYCLE
      END IF
      m=m+1
      newg(m)=newg(k)
      itmp(m)=rown(isort(k))
    END DO
    nnew=m
    DEALLOCATE(isort,rown)
  END IF
  !---Append to the subdomain and request their rows from their owners
  CALL grow_arrays(nd+nnew)
  DO k=1,nnew
    dlist(nd+k)=newg(k)
    newo(nd+k)=itmp(k)
    IF(itmp(k)==me)CALL oft_abort('Owned DOF outside the owned set','native_overlap_setup',__FILE__)
  END DO
  CALL request_rows(nnew,newg,itmp,nd)
  nd=nd+nnew
  DEALLOCATE(front)
  ALLOCATE(front(MAX(1,nnew)))
  IF(nnew>0)front(1:nnew)=newg(1:nnew)
  nfront=nnew
  DEALLOCATE(newg,itmp,dsorted)
END DO
!---Rows received from owners follow the complete owned rows in cvals
CALL build_plan(self%plan_ghost,ghost_send,ghost_recv,0,[0],[0])
ALLOCATE(self%cvals(MAX(1,self%ncr+ngv)))
self%cvals=0.d0
!------------------------------------------------------------------------------
! Subdomain numbering: local DOFs (increasing local index), then ghosts (by global index)
!------------------------------------------------------------------------------
self%n=nd
ALLOCATE(self%gsub(nd),self%osub(nd),sub_of_d(MAX(1,nd)),itmp(MAX(1,nd)))
self%nloc=0
DO k=1,nd
  itmp(k)=local_index(dlist(k))
  IF(itmp(k)>0)self%nloc=self%nloc+1
END DO
ALLOCATE(self%lsub(self%nloc))
!---Local DOFs by local index
ALLOCATE(isort(MAX(1,nd)),tmp8(MAX(1,nd)))
m=0
DO k=1,nd
  IF(itmp(k)==0)CYCLE
  m=m+1
  tmp8(m)=INT(itmp(k),i8)
  isort(m)=k
END DO
IF(m>0)CALL sort_array(tmp8,isort,m)
DO k=1,m
  sub_of_d(isort(k))=k
  self%lsub(k)=INT(tmp8(k),4)
END DO
!---Ghosts by global index
m=0
DO k=1,nd
  IF(itmp(k)/=0)CYCLE
  m=m+1
  tmp8(m)=dlist(k)
  isort(m)=k
END DO
IF(m>0)CALL sort_array(tmp8,isort,m)
DO k=1,m
  sub_of_d(isort(k))=self%nloc+k
END DO
DEALLOCATE(isort,tmp8,itmp)
DO k=1,nd
  self%gsub(sub_of_d(k))=dlist(k)
  self%osub(sub_of_d(k))=newo(k)
END DO
self%nown=COUNT(self%osub==me)
!---Sorted subdomain global indices for lookups
IF(ALLOCATED(dsorted))DEALLOCATE(dsorted)
IF(ALLOCATED(dperm))DEALLOCATE(dperm)
ALLOCATE(dsorted(MAX(1,nd)),dperm(MAX(1,nd)))
IF(nd>0)THEN
  dsorted(1:nd)=self%gsub
  dperm(1:nd)=[(k,k=1,nd)]
  CALL sort_array(dsorted,dperm,nd)
END IF
!------------------------------------------------------------------------------
! Subdomain matrix: entries of each subdomain row with columns in the subdomain
!------------------------------------------------------------------------------
ALLOCATE(kr_sub(nd+1),dlist_own(nd))
DO k=1,nd
  dlist_own(sub_of_d(k))=dsrc(k)
END DO
kr_sub(1)=1
DO r=1,nd
  kr_sub(r+1)=kr_sub(r)+count_sub_row(dlist_own(r))
END DO
nsub_nnz=kr_sub(nd+1)-1
ALLOCATE(lc(MAX(1,nsub_nnz)),lpos(MAX(1,nsub_nnz)))
DO r=1,nd
  CALL fill_sub_row(dlist_own(r),kr_sub(r))
  jp=kr_sub(r)
  jn=kr_sub(r+1)-jp
  IF(jn>1)CALL sort_array(lc(jp:jp+jn-1),lpos(jp:jp+jn-1),jn)
END DO
ALLOCATE(kr(nd+1))
kr=kr_sub
IF(nsub_nnz==0)THEN
  DEALLOCATE(lc)
  ALLOCATE(lc(0))
END IF
ALLOCATE(self%amat)
self%amat%force_local=.TRUE.
self%amat%nr=nd
self%amat%nc=nd
self%amat%nrg=nd
self%amat%ncg=nd
self%amat%nrslice=nd
self%amat%ncslice=nd
self%amat%ni=1
self%amat%nj=1
ALLOCATE(self%amat%i_map(1),self%amat%j_map(1))
self%amat%i_map(1)%n=nd
self%amat%j_map(1)%n=nd
CALL native_pattern_create(pattern,nd,nd,kr,lc)
CALL native_matrix_set_pattern(self%amat,pattern)
ALLOCATE(self%amat%M(MAX(1,nsub_nnz)))
CALL scatter_plan_setup_local(self%plan_vals,nsub_nnz,lpos)
DEALLOCATE(kr_sub,lpos,dlist_own)
!------------------------------------------------------------------------------
! Vector plans: local DOFs locally, ghosts from their owners
!------------------------------------------------------------------------------
DO q=0,np-1
  ghost_send(q)%n=0
  ghost_recv(q)%n=0
END DO
scnt=0
DO k=self%nloc+1,nd
  scnt(self%osub(k))=scnt(self%osub(k))+1
END DO
ALLOCATE(sbuf(MAX(1,SUM(scnt))))
CALL csr_offsets(scnt,ipos)
DO k=self%nloc+1,nd
  q=self%osub(k)
  sbuf(ipos(q))=self%gsub(k)
  ipos(q)=ipos(q)+1
  CALL ilist_push(ghost_recv(q),k)
END DO
CALL oft_exchange_i8(scnt,sbuf,rcnt,rbuf)
k=0
DO q=0,np-1
  DO j=1,rcnt(q)
    k=k+1
    r=local_index(rbuf(k))
    IF(r==0)CALL oft_abort('Requested value not held locally','native_overlap_setup',__FILE__)
    CALL ilist_push(ghost_send(q),r)
  END DO
END DO
DEALLOCATE(sbuf,rbuf)
CALL build_plan(self%plan_vec,ghost_send,ghost_recv,self%nloc,self%lsub,[(k,k=1,self%nloc)])
!---Owned subdomain DOFs (all held locally)
ALLOCATE(isort(MAX(1,self%nown)),itmp(MAX(1,self%nown)))
m=0
DO k=1,self%nloc
  IF(self%osub(k)/=me)CYCLE
  m=m+1
  isort(m)=self%lsub(k)
  itmp(m)=k
END DO
CALL scatter_plan_setup_local(self%plan_own,m,isort,itmp)
DEALLOCATE(isort,itmp)
ALLOCATE(self%acc(self%nr))
!---Cleanup
DEALLOCATE(gidx,gsorted,gperm,owner,orow,olist,scnt,rcnt,ipos,ckr,cg,co,gkr,gg,go)
DEALLOCATE(dlist,newo,dsrc,front,dsorted,dperm,sub_of_d)
IF(ALLOCATED(rbuf))DEALLOCATE(rbuf)
DEALLOCATE(rows_send,rows_recv,ghost_send,ghost_recv)
CALL native_overlap_update(self)
DEBUG_STACK_POP
CONTAINS
!> Local index of a global index (0 if not held locally)
FUNCTION local_index(g) RESULT(il)
INTEGER(i8), INTENT(in) :: g
INTEGER(i4) :: il
INTEGER(i8) :: ipos8
il=0
IF(mat%nr==0)RETURN
ipos8=search_array(g,gsorted,INT(mat%nr,i8))
IF(ipos8>0)il=gperm(ipos8)
END FUNCTION local_index
!> Location in cvals of column `g` of complete owned row `r`
FUNCTION complete_pos(r,g) RESULT(pos)
INTEGER(i4), INTENT(in) :: r
INTEGER(i8), INTENT(in) :: g
INTEGER(i4) :: pos
pos=INT(search_array(g,cg(ckr(r):ckr(r+1)-1),INT(ckr(r+1)-ckr(r),i8)),4)
IF(pos==0)CALL oft_abort('Entry missing from complete row','native_overlap_setup',__FILE__)
pos=pos+ckr(r)-1
END FUNCTION complete_pos
!> Position in dlist of a subdomain global index (0 if absent)
FUNCTION dlist_index(g) RESULT(pos)
INTEGER(i8), INTENT(in) :: g
INTEGER(i4) :: pos
INTEGER(i4) :: ilow,ihigh,imid
pos=0
ilow=1
ihigh=nd
DO WHILE(ilow<=ihigh)
  imid=(ilow+ihigh)/2
  IF(dsorted(imid)==g)THEN
    pos=dperm_tmp(imid)
    RETURN
  ELSE IF(dsorted(imid)<g)THEN
    ilow=imid+1
  ELSE
    ihigh=imid-1
  END IF
END DO
END FUNCTION dlist_index
!> Permutation of the sorted copy of dlist (kept in itmp-independent storage)
FUNCTION dperm_tmp(k) RESULT(j)
INTEGER(i4), INTENT(in) :: k
INTEGER(i4) :: j
j=dperm(k)
END FUNCTION dperm_tmp
!> Sorted copy of the first n entries of a list, with permutation in dperm
SUBROUTINE sorted_copy(list,n,sorted)
INTEGER(i8), INTENT(in) :: list(:)
INTEGER(i4), INTENT(in) :: n
INTEGER(i8), ALLOCATABLE, INTENT(out) :: sorted(:)
ALLOCATE(sorted(MAX(1,n)))
IF(ALLOCATED(dperm))DEALLOCATE(dperm)
ALLOCATE(dperm(MAX(1,n)))
IF(n>0)THEN
  sorted(1:n)=list(1:n)
  dperm(1:n)=[(k,k=1,n)]
  CALL sort_array(sorted,dperm,n)
END IF
END SUBROUTINE sorted_copy
!> Number of entries in the rows of the front (upper bound on new DOFs)
FUNCTION count_front_cols() RESULT(nc)
INTEGER(i4) :: nc,kk,jj,src
nc=0
DO kk=1,nfront
  jj=dlist_index(front(kk))
  src=dsrc(jj)
  IF(src>0)THEN
    nc=nc+ckr(src+1)-ckr(src)
  ELSE
    nc=nc+gkr(-src+1)-gkr(-src)
  END IF
END DO
END FUNCTION count_front_cols
!> Add columns of a row (owned row src>0, received row -src) that are not in the subdomain
SUBROUTINE row_entries(src)
INTEGER(i4), INTENT(in) :: src
INTEGER(i4) :: jj
IF(src>0)THEN
  DO jj=ckr(src),ckr(src+1)-1
    IF(dlist_index(cg(jj))>0)CYCLE
    nnew=nnew+1
    newg(nnew)=cg(jj)
    itmp(nnew)=co(jj)
  END DO
ELSE
  DO jj=gkr(-src),gkr(-src+1)-1
    IF(dlist_index(gg(jj))>0)CYCLE
    nnew=nnew+1
    newg(nnew)=gg(jj)
    itmp(nnew)=go(jj)
  END DO
END IF
END SUBROUTINE row_entries
!> Grow the subdomain lists to hold `nmax` DOFs
SUBROUTINE grow_arrays(nmax)
INTEGER(i4), INTENT(in) :: nmax
INTEGER(i8), ALLOCATABLE :: t8(:)
INTEGER(i4), ALLOCATABLE :: t4(:)
IF(SIZE(dlist)>=nmax)RETURN
ALLOCATE(t8(nmax))
t8(1:nd)=dlist(1:nd)
CALL MOVE_ALLOC(t8,dlist)
ALLOCATE(t4(nmax))
t4(1:nd)=newo(1:nd)
CALL MOVE_ALLOC(t4,newo)
ALLOCATE(t4(nmax))
t4(1:nd)=dsrc(1:nd)
CALL MOVE_ALLOC(t4,dsrc)
END SUBROUTINE grow_arrays
!> Request the complete rows of `nreq` DOFs (appended to dlist at `noff`) from their owners
SUBROUTINE request_rows(nreq,greq,oreq,noff)
INTEGER(i4), INTENT(in) :: nreq
INTEGER(i8), INTENT(in) :: greq(:)
INTEGER(i4), INTENT(in) :: oreq(:)
INTEGER(i4), INTENT(in) :: noff
INTEGER(i4) :: kk,jj,qq,rr,nn,ntot,mm
INTEGER(i4), ALLOCATABLE :: sc(:),rc(:),ip(:),ordr(:)
INTEGER(i8), ALLOCATABLE :: sb(:),rb(:),sb2(:),rb2(:)
INTEGER(i4), ALLOCATABLE :: t4(:)
INTEGER(i8), ALLOCATABLE :: t8(:)
ALLOCATE(sc(0:np-1),rc(0:np-1),ip(0:np-1))
!---Requested global indices, grouped by owner
sc=0
DO kk=1,nreq
  sc(oreq(kk))=sc(oreq(kk))+1
END DO
ALLOCATE(sb(MAX(1,nreq)),ordr(MAX(1,nreq)))
CALL csr_offsets(sc,ip)
DO kk=1,nreq
  sb(ip(oreq(kk)))=greq(kk)
  ordr(ip(oreq(kk)))=kk
  ip(oreq(kk))=ip(oreq(kk))+1
END DO
CALL oft_exchange_i8(sc,sb,rc,rb)
!---Owner: reply with the structure of each requested complete row and record value sources
ntot=0
DO kk=1,SUM(rc)
  rr=local_index(rb(kk))
  IF(rr==0)CALL oft_abort('Requested row not held','native_overlap_setup',__FILE__)
  IF(orow(rr)==0)CALL oft_abort('Requested row not owned','native_overlap_setup',__FILE__)
  ntot=ntot+1+2*(ckr(orow(rr)+1)-ckr(orow(rr)))
END DO
ALLOCATE(sb2(MAX(1,ntot)))
nn=0
kk=0
DO qq=0,np-1
  sc(qq)=0
  DO jj=1,rc(qq)
    kk=kk+1
    rr=orow(local_index(rb(kk)))
    nn=nn+1
    sb2(nn)=ckr(rr+1)-ckr(rr)
    sc(qq)=sc(qq)+1
    DO mm=ckr(rr),ckr(rr+1)-1
      sb2(nn+1)=cg(mm)
      sb2(nn+2)=co(mm)
      nn=nn+2
      sc(qq)=sc(qq)+2
      CALL ilist_push(ghost_send(qq),mm)
    END DO
  END DO
END DO
CALL oft_exchange_i8(sc,sb2,rc,rb2)
!---Requester: store received rows (in request order per owner) and their value locations
kk=0
nn=0
DO qq=0,np-1
  DO jj=1,COUNT(oreq(1:nreq)==qq)
    kk=kk+1
    ng_rows=ng_rows+1
    rr=INT(rb2(nn+1),4)
    !---Grow received row storage
    ALLOCATE(t4(ng_rows+1))
    t4(1:ng_rows)=gkr(1:ng_rows)
    t4(ng_rows+1)=gkr(ng_rows)+rr
    CALL MOVE_ALLOC(t4,gkr)
    ALLOCATE(t8(MAX(1,gkr(ng_rows+1)-1)))
    IF(gkr(ng_rows)>1)t8(1:gkr(ng_rows)-1)=gg(1:gkr(ng_rows)-1)
    CALL MOVE_ALLOC(t8,gg)
    ALLOCATE(t4(MAX(1,gkr(ng_rows+1)-1)))
    IF(gkr(ng_rows)>1)t4(1:gkr(ng_rows)-1)=go(1:gkr(ng_rows)-1)
    CALL MOVE_ALLOC(t4,go)
    DO mm=1,rr
      gg(gkr(ng_rows)+mm-1)=rb2(nn+2*mm)
      go(gkr(ng_rows)+mm-1)=INT(rb2(nn+2*mm+1),4)
      ngv=ngv+1
      CALL ilist_push(ghost_recv(qq),self%ncr+ngv)
    END DO
    nn=nn+1+2*rr
    !---Row source of the requested DOF
    dsrc(noff+ordr(kk))=-ng_rows
  END DO
END DO
DEALLOCATE(sc,rc,ip,sb,rb,sb2,rb2,ordr)
END SUBROUTINE request_rows
!> Number of entries of a row (source src) with columns in the subdomain
FUNCTION count_sub_row(src) RESULT(nc)
INTEGER(i4), INTENT(in) :: src
INTEGER(i4) :: nc,jj
nc=0
IF(src>0)THEN
  DO jj=ckr(src),ckr(src+1)-1
    IF(sub_index(cg(jj))>0)nc=nc+1
  END DO
ELSE
  DO jj=gkr(-src),gkr(-src+1)-1
    IF(sub_index(gg(jj))>0)nc=nc+1
  END DO
END IF
END FUNCTION count_sub_row
!> Fill the entries of a subdomain row starting at kfirst (columns and value locations)
SUBROUTINE fill_sub_row(src,kfirst)
INTEGER(i4), INTENT(in) :: src
INTEGER(i4), INTENT(in) :: kfirst
INTEGER(i4) :: jj,kk,sc_
kk=kfirst
IF(src>0)THEN
  DO jj=ckr(src),ckr(src+1)-1
    sc_=sub_index(cg(jj))
    IF(sc_==0)CYCLE
    lc(kk)=sc_
    lpos(kk)=jj
    kk=kk+1
  END DO
ELSE
  DO jj=gkr(-src),gkr(-src+1)-1
    sc_=sub_index(gg(jj))
    IF(sc_==0)CYCLE
    lc(kk)=sc_
    lpos(kk)=self%ncr+jj ! Received values are stored in the order of their entries
    kk=kk+1
  END DO
END IF
END SUBROUTINE fill_sub_row
!> Subdomain index of a global index (0 if not in the subdomain)
FUNCTION sub_index(g) RESULT(si)
INTEGER(i8), INTENT(in) :: g
INTEGER(i4) :: si
INTEGER(i8) :: p8
si=0
IF(nd==0)RETURN
p8=search_array(g,dsorted(1:nd),INT(nd,i8))
IF(p8>0)si=dperm(p8)
END FUNCTION sub_index
END SUBROUTINE native_overlap_setup
!---------------------------------------------------------------------------------
!> Refresh the subdomain matrix values from the parent matrix (collective)
!!
!! Updates the parent's full local matrix, completes the owned rows, transfers the
!! rows of non-owned subdomain DOFs from their owners and gathers the subdomain matrix
!! values, all with the scatter plans built by @ref native_overlap_setup.
!---------------------------------------------------------------------------------
SUBROUTINE native_overlap_update(self)
TYPE(oft_native_overlap), INTENT(inout) :: self !< Overlap object
DEBUG_STACK_PUSH
CALL self%parent%update_slice
CALL scatter_begin(self%plan_rows,self%parent%Mfull)
CALL scatter_end(self%plan_rows,self%cvals,OFT_SCATTER_INSERT)
CALL scatter_begin(self%plan_ghost,self%cvals)
CALL scatter_end(self%plan_ghost,self%cvals,OFT_SCATTER_INSERT)
CALL scatter_local(self%plan_vals,self%cvals,self%amat%M,OFT_SCATTER_INSERT)
DEBUG_STACK_POP
END SUBROUTINE native_overlap_update
!---------------------------------------------------------------------------------
!> Gather subdomain values from local values of a distributed vector (collective)
!---------------------------------------------------------------------------------
SUBROUTINE native_overlap_gather(self,x,xs)
TYPE(oft_native_overlap), INTENT(inout) :: self !< Overlap object
REAL(r8), INTENT(in) :: x(:) !< Local values of the distributed vector [nr]
REAL(r8), INTENT(inout) :: xs(:) !< Subdomain values [n]
DEBUG_STACK_PUSH
CALL scatter_begin(self%plan_vec,x)
CALL scatter_end(self%plan_vec,xs,OFT_SCATTER_INSERT)
DEBUG_STACK_POP
END SUBROUTINE native_overlap_gather
!---------------------------------------------------------------------------------
!> Write subdomain values back to local values of a distributed vector (collective)
!!
!! - `additive=.FALSE.` (restricted): owned entries are set from the subdomain; other
!!   local entries are left unchanged and must be synchronized by the caller
!!   (`stitch(0)`).
!! - `additive=.TRUE.`: `x` is set to the sum over all subdomains containing each DOF
!!   (contributions of ghosts are added on their owner, local contributions are
!!   summed across copies with `stitch(1)`), consistent on all copies.
!---------------------------------------------------------------------------------
SUBROUTINE native_overlap_scatter(self,xs,x,additive)
TYPE(oft_native_overlap), INTENT(inout) :: self !< Overlap object
REAL(r8), INTENT(in) :: xs(:) !< Subdomain values [n]
REAL(r8), INTENT(inout) :: x(:) !< Local values of the distributed vector [nr]
LOGICAL, INTENT(in) :: additive !< Sum contributions of all subdomains?
INTEGER(i4) :: i
DEBUG_STACK_PUSH
IF(additive)THEN
  !$omp parallel do simd if(self%nr>OFT_OMP_VTHRESH)
  DO i=1,self%nr
    self%acc(i)=0.d0
  END DO
  CALL scatter_begin(self%plan_vec,xs,reverse=.TRUE.)
  CALL scatter_end(self%plan_vec,self%acc,OFT_SCATTER_ADD)
  CALL oft_global_stitch(self%seam,self%acc,1)
  !$omp parallel do simd if(self%nr>OFT_OMP_VTHRESH)
  DO i=1,self%nr
    x(i)=self%acc(i)
  END DO
ELSE
  CALL scatter_local(self%plan_own,xs,x,OFT_SCATTER_INSERT,reverse=.TRUE.)
END IF
DEBUG_STACK_POP
END SUBROUTINE native_overlap_scatter
!---------------------------------------------------------------------------------
!> Release all storage of an overlap object
!---------------------------------------------------------------------------------
SUBROUTINE native_overlap_delete(self)
TYPE(oft_native_overlap), INTENT(inout) :: self !< Overlap object
IF(ALLOCATED(self%gsub))DEALLOCATE(self%gsub)
IF(ALLOCATED(self%lsub))DEALLOCATE(self%lsub)
IF(ALLOCATED(self%osub))DEALLOCATE(self%osub)
IF(ALLOCATED(self%cvals))DEALLOCATE(self%cvals)
IF(ALLOCATED(self%acc))DEALLOCATE(self%acc)
CALL scatter_plan_delete(self%plan_rows)
CALL scatter_plan_delete(self%plan_ghost)
CALL scatter_plan_delete(self%plan_vals)
CALL scatter_plan_delete(self%plan_vec)
CALL scatter_plan_delete(self%plan_own)
IF(ASSOCIATED(self%amat))THEN
  CALL self%amat%delete
  DEALLOCATE(self%amat)
END IF
NULLIFY(self%parent,self%seam)
self%n=0
self%nloc=0
self%nown=0
self%ncr=0
END SUBROUTINE native_overlap_delete
!---------------------------------------------------------------------------------
!> Starting offset of each process' part of a buffer grouped by process
!---------------------------------------------------------------------------------
SUBROUTINE csr_offsets(cnt,ipos)
INTEGER(i4), INTENT(in) :: cnt(0:) !< Number of entries for each process
INTEGER(i4), INTENT(out) :: ipos(0:) !< First entry for each process
INTEGER(i4) :: q
ipos(0)=1
DO q=1,SIZE(cnt)-1
  ipos(q)=ipos(q-1)+cnt(q-1)
END DO
END SUBROUTINE csr_offsets
!---------------------------------------------------------------------------------
!> Append a value to a growable list
!---------------------------------------------------------------------------------
SUBROUTINE ilist_push(list,val)
TYPE(ilist), INTENT(inout) :: list !< List
INTEGER(i4), INTENT(in) :: val !< Value to append
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: tmp
IF(.NOT.ALLOCATED(list%v))ALLOCATE(list%v(16))
IF(list%n>=SIZE(list%v))THEN
  ALLOCATE(tmp(2*SIZE(list%v)))
  tmp(1:list%n)=list%v(1:list%n)
  CALL MOVE_ALLOC(tmp,list%v)
END IF
list%n=list%n+1
list%v(list%n)=val
END SUBROUTINE ilist_push
!---------------------------------------------------------------------------------
!> Build a scatter plan from per-process send and receive lists and a local part
!!
!! Neighbours are the processes with entries in either list (the calling process
!! must not appear).
!---------------------------------------------------------------------------------
SUBROUTINE build_plan(plan,slist,rlist,nloc,lsrc,ldst)
TYPE(oft_scatter_plan), INTENT(inout) :: plan !< Plan
TYPE(ilist), INTENT(in) :: slist(0:) !< Source indices sent to each process
TYPE(ilist), INTENT(in) :: rlist(0:) !< Destination indices of values received from each process
INTEGER(i4), INTENT(in) :: nloc !< Number of local entries
INTEGER(i4), INTENT(in) :: lsrc(:) !< Local source indices [nloc]
INTEGER(i4), INTENT(in) :: ldst(:) !< Local destination indices [nloc]
INTEGER(i4) :: q,nn,j,ns,nr
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: nbr,ksend,krecv,ls,lr
nn=0
ns=0
nr=0
DO q=0,SIZE(slist)-1
  IF(slist(q)%n+rlist(q)%n==0)CYCLE
  IF(q==oft_env%rank)CALL oft_abort('Self exchange','build_plan',__FILE__)
  nn=nn+1
  ns=ns+slist(q)%n
  nr=nr+rlist(q)%n
END DO
ALLOCATE(nbr(MAX(1,nn)),ksend(nn+1),krecv(nn+1),ls(MAX(1,ns)),lr(MAX(1,nr)))
nn=0
ksend(1)=1
krecv(1)=1
DO q=0,SIZE(slist)-1
  IF(slist(q)%n+rlist(q)%n==0)CYCLE
  nn=nn+1
  nbr(nn)=q
  DO j=1,slist(q)%n
    ls(ksend(nn)+j-1)=slist(q)%v(j)
  END DO
  ksend(nn+1)=ksend(nn)+slist(q)%n
  DO j=1,rlist(q)%n
    lr(krecv(nn)+j-1)=rlist(q)%v(j)
  END DO
  krecv(nn+1)=krecv(nn)+rlist(q)%n
END DO
CALL scatter_plan_setup(plan,nn,nbr(1:nn),ksend,ls,krecv,lr,nloc,lsrc,ldst)
CALL scatter_plan_check(plan)
DEALLOCATE(nbr,ksend,krecv,ls,lr)
END SUBROUTINE build_plan
END MODULE oft_native_overlap_mod
