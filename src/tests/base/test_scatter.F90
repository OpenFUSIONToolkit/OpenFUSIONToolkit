!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_scatter.F90
!
!> Regression tests for scatter plans (@ref oft_scatter_base::oft_scatter_plan)
!!
!! Each process exchanges values with its neighbours on a ring (none on one
!! process). Every send, receive and local index list is a function of the ranks,
!! the neighbour and the entry, so each process computes the values it should
!! receive directly from the source values of its neighbours. Plans with unique and
!! with repeated indices (insert: the last entry wins; add: every entry counts), with
!! and without local weights, are applied forward and reverse, in insert and add
!! mode, to real and complex values, with and without a mask and the local part,
!! through `scatter_end`, the split `scatter_wait_recv`/`scatter_unpack`/
!! `scatter_wait_send` calls and `scatter_local`.
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup testing
!--------------------------------------------------------------------------------
MODULE test_scatter_helpers
USE oft_base
USE oft_scatter
IMPLICIT NONE
INTEGER(i4), PARAMETER :: m = 24 !< Length of the local arrays
INTEGER(i4), PARAMETER :: ns = 3 !< Entries exchanged with each neighbour
INTEGER(i4), PARAMETER :: nl = 4 !< Local entries
REAL(r8), PARAMETER :: wts(nl) = [0.5d0, 0.d0, 2.d0, -1.d0] !< Local weights (one zero)
INTEGER(i4) :: nfail = 0 !< Number of failed checks
INTEGER(i4) :: ncheck = 0 !< Number of checks
CONTAINS
! Neighbours of rank r on a ring of np processes (sorted, without r)
FUNCTION ring_nbrs(r,np) RESULT(nbr)
INTEGER(i4), INTENT(in) :: r,np
INTEGER(i4), ALLOCATABLE :: nbr(:)
INTEGER(i4) :: a,b
IF(np==1)THEN
  ALLOCATE(nbr(0))
ELSE IF(np==2)THEN
  nbr=[1-r]
ELSE
  a=MOD(r+np-1,np)
  b=MOD(r+1,np)
  nbr=[MIN(a,b),MAX(a,b)]
END IF
END FUNCTION ring_nbrs
! Position of rank b among the neighbours of rank a
FUNCTION nbr_pos(a,b,np) RESULT(j)
INTEGER(i4), INTENT(in) :: a,b,np
INTEGER(i4) :: j
INTEGER(i4), ALLOCATABLE :: nbr(:)
nbr=ring_nbrs(a,np)
DO j=1,SIZE(nbr)
  IF(nbr(j)==b)RETURN
END DO
CALL oft_abort('Not a neighbour','nbr_pos',__FILE__)
END FUNCTION nbr_pos
! Source index of entry k sent to the j-th neighbour (forward)
FUNCTION isend(rep,j,k) RESULT(i)
LOGICAL, INTENT(in) :: rep
INTEGER(i4), INTENT(in) :: j,k
INTEGER(i4) :: i
IF(rep)THEN
  i=1+MOD(j+2*k,4)
ELSE
  i=(j-1)*ns+k
END IF
END FUNCTION isend
! Destination index of entry k received from the j-th neighbour (forward)
FUNCTION irecv(rep,j,k) RESULT(i)
LOGICAL, INTENT(in) :: rep
INTEGER(i4), INTENT(in) :: j,k
INTEGER(i4) :: i
IF(rep)THEN
  i=1+MOD(j+k,3)
ELSE
  i=10+(j-1)*ns+k
END IF
END FUNCTION irecv
! Source index of local entry k (forward)
FUNCTION ilsrc(rep,k) RESULT(i)
LOGICAL, INTENT(in) :: rep
INTEGER(i4), INTENT(in) :: k
INTEGER(i4) :: i
IF(rep)THEN
  i=1+MOD(k,3)
ELSE
  i=2*ns+k
END IF
END FUNCTION ilsrc
! Destination index of local entry k (forward)
FUNCTION ildst(rep,k) RESULT(i)
LOGICAL, INTENT(in) :: rep
INTEGER(i4), INTENT(in) :: k
INTEGER(i4) :: i
IF(rep)THEN
  i=1+MOD(k,3)
ELSE
  i=16+k
END IF
END FUNCTION ildst
! Source value at index i on rank p
FUNCTION xval(p,i) RESULT(v)
INTEGER(i4), INTENT(in) :: p,i
COMPLEX(c8) :: v
v=CMPLX(100*p+i,-(p+1)*i,c8)
END FUNCTION xval
! Initial destination value at index i
FUNCTION yval(i) RESULT(v)
INTEGER(i4), INTENT(in) :: i
COMPLEX(c8) :: v
v=CMPLX(-i,0.25d0*i,c8)
END FUNCTION yval
! Mask of packed entry i
FUNCTION mval(i) RESULT(lm)
INTEGER(i4), INTENT(in) :: i
LOGICAL :: lm
lm=(MOD(i,3)/=0)
END FUNCTION mval
! Build the plan of this process
SUBROUTINE build_plan(plan,rep,weighted)
TYPE(oft_scatter_plan), INTENT(inout) :: plan
LOGICAL, INTENT(in) :: rep,weighted
INTEGER(i4) :: j,k,nn
INTEGER(i4), ALLOCATABLE :: nbr(:),kp(:),ls(:),lr(:),lsrc(:),ldst(:)
nbr=ring_nbrs(oft_env%rank,oft_env%nprocs)
nn=SIZE(nbr)
ALLOCATE(kp(nn+1),ls(MAX(1,nn*ns)),lr(MAX(1,nn*ns)),lsrc(nl),ldst(nl))
kp(1)=1
DO j=1,nn
  kp(j+1)=kp(j)+ns
  DO k=1,ns
    ls((j-1)*ns+k)=isend(rep,j,k)
    lr((j-1)*ns+k)=irecv(rep,j,k)
  END DO
END DO
DO k=1,nl
  lsrc(k)=ilsrc(rep,k)
  ldst(k)=ildst(rep,k)
END DO
IF(weighted)THEN
  CALL scatter_plan_setup(plan,nn,nbr,kp,ls,kp,lr,nl,lsrc,ldst,wts)
ELSE
  CALL scatter_plan_setup(plan,nn,nbr,kp,ls,kp,lr,nl,lsrc,ldst)
END IF
CALL scatter_plan_check(plan)
END SUBROUTINE build_plan
! Expected result of a scatter (communicating, or local part only with `lonly`)
SUBROUTINE reference(rep,weighted,rev,mode,masked,noloc,lonly,yref)
LOGICAL, INTENT(in) :: rep,weighted,rev,masked,noloc,lonly
INTEGER(i4), INTENT(in) :: mode
COMPLEX(c8), INTENT(out) :: yref(m)
INTEGER(i4) :: i,j,k,a,jb,r,np,nn,isrc,idst,ipack
INTEGER(i4), ALLOCATABLE :: nbr(:)
COMPLEX(c8) :: val
r=oft_env%rank
np=oft_env%nprocs
DO i=1,m
  yref(i)=yval(i)
END DO
nbr=ring_nbrs(r,np)
nn=SIZE(nbr)
!---Remote entries, in neighbour then entry order
IF(.NOT.lonly)THEN
  DO j=1,nn
    a=nbr(j)
    jb=nbr_pos(a,r,np)
    DO k=1,ns
      IF(rev)THEN
        isrc=irecv(rep,jb,k)
        idst=isend(rep,j,k)
      ELSE
        isrc=isend(rep,jb,k)
        idst=irecv(rep,j,k)
      END IF
      ipack=(jb-1)*ns+k
      val=xval(a,isrc)
      IF(masked.AND.(.NOT.mval(ipack)))val=(0.d0,0.d0)
      CALL put(idst,val)
    END DO
  END DO
END IF
!---Local entries
IF(noloc)RETURN
DO k=1,nl
  IF(rev)THEN
    isrc=ildst(rep,k)
    idst=ilsrc(rep,k)
  ELSE
    isrc=ilsrc(rep,k)
    idst=ildst(rep,k)
  END IF
  val=xval(r,isrc)
  IF(masked.AND.(.NOT.mval(nn*ns+k)))val=(0.d0,0.d0)
  IF(rev.AND.weighted)THEN
    IF(wts(k)==0.d0)CYCLE
    val=wts(k)*val
  END IF
  CALL put(idst,val)
END DO
CONTAINS
SUBROUTINE put(ii,vv)
INTEGER(i4), INTENT(in) :: ii
COMPLEX(c8), INTENT(in) :: vv
IF(mode==OFT_SCATTER_ADD)THEN
  yref(ii)=yref(ii)+vv
ELSE
  yref(ii)=vv
END IF
END SUBROUTINE put
END SUBROUTINE reference
! Compare a result with the reference
SUBROUTINE check(label,y,yref)
CHARACTER(LEN=*), INTENT(in) :: label
COMPLEX(c8), INTENT(in) :: y(m),yref(m)
REAL(r8) :: err
ncheck=ncheck+1
err=MAXVAL(ABS(y-yref))
IF(err>1.d-12*MAX(1.d0,MAXVAL(ABS(yref))))THEN
  nfail=nfail+1
  WRITE(*,'(A,I4,2A,ES12.4)')'Rank ',oft_env%rank,': FAILED ',label,err
END IF
END SUBROUTINE check
! Run one case on real and complex values
SUBROUTINE run_case(plan,rep,weighted,rev,mode,masked,noloc,split)
TYPE(oft_scatter_plan), INTENT(inout) :: plan
LOGICAL, INTENT(in) :: rep,weighted,rev,masked,noloc,split
INTEGER(i4), INTENT(in) :: mode
INTEGER(i4) :: i,np
LOGICAL :: mask(m)
REAL(r8) :: xr(m),yr(m)
COMPLEX(c8) :: xc(m),yc(m),yref(m)
CHARACTER(LEN=80) :: label
WRITE(label,'(A,6L2,I2)')'rep/wt/rev/mask/noloc/split, mode',rep,weighted,rev,masked,noloc,split,mode
DO i=1,m
  xc(i)=xval(oft_env%rank,i)
  yc(i)=yval(i)
  mask(i)=mval(i)
END DO
xr=REAL(xc,r8)
yr=REAL(yc,r8)
CALL reference(rep,weighted,rev,mode,masked,noloc,.FALSE.,yref)
!---Real values
IF(masked)THEN
  CALL scatter_begin(plan,xr,reverse=rev,mask=mask,nolocal=noloc)
ELSE
  CALL scatter_begin(plan,xr,reverse=rev,nolocal=noloc)
END IF
IF(split)THEN
  CALL scatter_wait_recv(plan,yr)
  CALL scatter_unpack(plan,yr,mode,nolocal=noloc)
  CALL scatter_wait_send(plan,yr)
ELSE
  CALL scatter_end(plan,yr,mode,nolocal=noloc)
END IF
CALL check('real '//TRIM(label),CMPLX(yr,0.d0,c8),CMPLX(REAL(yref,r8),0.d0,c8))
!---Complex values
IF(masked)THEN
  CALL scatter_begin(plan,xc,reverse=rev,mask=mask,nolocal=noloc)
ELSE
  CALL scatter_begin(plan,xc,reverse=rev,nolocal=noloc)
END IF
IF(split)THEN
  CALL scatter_wait_recv(plan,yc)
  CALL scatter_unpack(plan,yc,mode,nolocal=noloc)
  CALL scatter_wait_send(plan,yc)
ELSE
  CALL scatter_end(plan,yc,mode,nolocal=noloc)
END IF
CALL check('complex '//TRIM(label),yc,yref)
END SUBROUTINE run_case
! Local part only (scatter_local) on real and complex values
SUBROUTINE run_local(plan,rep,weighted,rev,mode)
TYPE(oft_scatter_plan), INTENT(inout) :: plan
LOGICAL, INTENT(in) :: rep,weighted,rev
INTEGER(i4), INTENT(in) :: mode
INTEGER(i4) :: i
REAL(r8) :: xr(m),yr(m)
COMPLEX(c8) :: xc(m),yc(m),yref(m)
CHARACTER(LEN=80) :: label
WRITE(label,'(A,3L2,I2)')'local rep/wt/rev, mode',rep,weighted,rev,mode
DO i=1,m
  xc(i)=xval(oft_env%rank,i)
  yc(i)=yval(i)
END DO
xr=REAL(xc,r8)
yr=REAL(yc,r8)
CALL reference(rep,weighted,rev,mode,.FALSE.,.FALSE.,.TRUE.,yref)
CALL scatter_local(plan,xr,yr,mode,reverse=rev)
CALL check('real '//TRIM(label),CMPLX(yr,0.d0,c8),CMPLX(REAL(yref,r8),0.d0,c8))
CALL scatter_local(plan,xc,yc,mode,reverse=rev)
CALL check('complex '//TRIM(label),yc,yref)
END SUBROUTINE run_local
END MODULE test_scatter_helpers

PROGRAM test_scatter
USE oft_base
USE oft_scatter
USE test_scatter_helpers
IMPLICIT NONE
INTEGER(i4) :: irep,iwt,idir,imode,imask,iloc,io_unit
LOGICAL :: rep,weighted,rev,split
INTEGER(i4) :: modes(2)
TYPE(oft_scatter_plan) :: plan
CALL oft_init
modes=[OFT_SCATTER_INSERT,OFT_SCATTER_ADD]
DO irep=0,1
  rep=(irep==1)
  DO iwt=0,1
    weighted=(iwt==1)
    CALL build_plan(plan,rep,weighted)
    !---Index uniqueness flags (unique lists: threaded loops; repeated: serial)
    ncheck=ncheck+1
    IF((plan%recv_unique.EQV.rep).OR.(plan%send_unique.EQV.rep))THEN
      nfail=nfail+1
      WRITE(*,'(A,I4,A,3L2)')'Rank ',oft_env%rank,': FAILED uniqueness flags',rep, &
        plan%recv_unique,plan%send_unique
    END IF
    DO idir=0,1
      rev=(idir==1)
      DO imode=1,2
        split=(imode==2)
        DO imask=0,1
          DO iloc=0,1
            CALL run_case(plan,rep,weighted,rev,modes(imode),imask==1,iloc==1,split)
          END DO
        END DO
        CALL run_local(plan,rep,weighted,rev,modes(imode))
      END DO
    END DO
    CALL scatter_plan_delete(plan)
  END DO
END DO
!---Report
nfail=oft_mpi_sum(nfail)
IF(oft_env%head_proc)THEN
  WRITE(*,'(A,I6,A,I4)')'Scatter plan checks per process: ',ncheck,', failed (all processes): ',nfail
  OPEN(NEWUNIT=io_unit,FILE='scatter.results')
  WRITE(io_unit,'(I0)')nfail
  WRITE(io_unit,'(I0)')ncheck
  CLOSE(io_unit)
END IF
CALL oft_finalize
END PROGRAM test_scatter
