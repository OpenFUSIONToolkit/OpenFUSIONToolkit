!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_fem_qcache.F90
!
!> Regression tests for the quadrature-point basis cache (\ref fem_quad_cache).
!!
!! For Lagrange (order p, values + gradients), H^1 (order p+1, values + gradients) and
!! H(Curl) (order p, values + curls) spaces on the finest mesh level, and for their p-1 levels
!! evaluated at the finest quadrature, both cache modes (FULL and REFERENCE) are compared with
!! the direct basis evaluators at every (cell, quadrature point), using the Jacobian evaluated
!! at every point as the reference geometry. REFERENCE mode is checked both with gathered table
!! rows ('ref') and in the canonical frame with per-cell orientation transforms ('refc'). Also
!! checked: field reconstruction from local coefficients, test-function accumulation, a local
!! matrix, det*w and physical quadrature points.
!!
!! Results are written to 'qcache.results' as "<name> <max normalized error>" lines.
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_fem_qcache
USE oft_base
USE oft_mesh_type, ONLY: oft_mesh
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE fem_base, ONLY: oft_ml_fem_type, oft_afem_type
USE oft_lag_basis, ONLY: oft_lag_setup, oft_scalar_fem, oft_lag_eval_all, oft_lag_geval_all
USE oft_h1_basis, ONLY: oft_h1_setup, oft_h1_fem, oft_h1_eval_all, oft_h1_geval_all
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_fem, oft_hcurl_eval_all, &
  oft_hcurl_ceval_all, oft_hcurl_get_cgops
USE fem_quad_cache
IMPLICIT NONE
INTEGER(i4) :: order,ierr,io_unit,i,j,est_seen(2),fst_seen(8)
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_lag,ML_h1,ML_hcurl
TYPE(oft_qcache_geom), TARGET :: geom
CLASS(oft_mesh), POINTER :: mesh
REAL(r8) :: err_detw,err_pt,gop(3,4),cgop(3,6),detw,det,pt(3),scale
LOGICAL :: test_coarse
NAMELIST/test_qcache_options/order,test_coarse
!---Initialize enviroment
CALL oft_init
order=1
test_coarse=.TRUE.
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_qcache_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
mesh=>mg_mesh%mesh
!---FE spaces (all p-levels on the finest mesh)
CALL oft_lag_setup(mg_mesh,order,ML_lag,minlev=mg_mesh%mgdim)
CALL oft_h1_setup(mg_mesh,order+1,ML_h1,minlev=mg_mesh%mgdim)
CALL oft_hcurl_setup(mg_mesh,order,ML_hcurl,minlev=mg_mesh%mgdim)
!---Geometry at the finest H(Curl) quadrature (as xMHD)
CALL geom%setup(mesh,ML_hcurl%levels(ML_hcurl%nlevels)%fe%quad,store_pts=.TRUE.,store_shared=.TRUE.)
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(A,I8,A,I4,A,I8,A,I8)')'Geometry: nc = ',geom%nc,', nqp = ',geom%nqp, &
    ', per-qp cells = ',geom%nvar,', shared cells = ',COUNT(geom%cell_shared)
END IF
!---Orientation states present on this mesh
est_seen=0
fst_seen=0
DO i=1,mesh%nc
  DO j=1,geom%ne
    est_seen(geom%orient(j,i))=est_seen(geom%orient(j,i))+1
  END DO
  DO j=geom%ne+1,geom%ne+geom%nf
    fst_seen(geom%orient(j,i))=fst_seen(geom%orient(j,i))+1
  END DO
END DO
IF(oft_env%head_proc)THEN
  WRITE(*,'(A,2I8)')'Edge states seen (local): ',est_seen
  WRITE(*,'(A,8I8)')'Face states seen (local): ',fst_seen(1:geom%nfs)
END IF
!---Geometry checks
err_detw=0.d0; err_pt=0.d0
!$omp parallel do private(j,gop,cgop,detw,det,pt) reduction(max:err_detw,err_pt)
DO i=1,mesh%nc
  DO j=1,geom%nqp
    CALL fem_qcache_geom_get(geom,i,j,gop,detw,cgop)
    CALL mesh%jacobian(i,geom%pts(:,j),gop,det)
    err_detw=MAX(err_detw,ABS(detw-det*geom%wts(j))/ABS(det*geom%wts(j)))
    pt=mesh%log2phys(i,geom%pts(:,j))
    err_pt=MAX(err_pt,MAXVAL(ABS(pt-geom%pt(:,j,i))))
  END DO
END DO
scale=MAXVAL(ABS(mesh%r))
err_detw=oft_mpi_max(err_detw)
err_pt=oft_mpi_max(err_pt/scale)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='qcache.results')
  WRITE(io_unit,'(A,ES12.4)')'geom_detw ',err_detw
  WRITE(io_unit,'(A,ES12.4)')'geom_pt ',err_pt
  CLOSE(io_unit)
END IF
!---Basis checks
CALL check_space(ML_lag%levels(ML_lag%nlevels)%fe,'lag')
CALL check_space(ML_h1%levels(ML_h1%nlevels)%fe,'h1')
CALL check_space(ML_hcurl%levels(ML_hcurl%nlevels)%fe,'hcurl')
IF(test_coarse.AND.order>1)THEN
  CALL check_space(ML_lag%levels(ML_lag%nlevels-1)%fe,'lag_pm1')
  CALL check_space(ML_h1%levels(ML_h1%nlevels-1)%fe,'h1_pm1')
  CALL check_space(ML_hcurl%levels(ML_hcurl%nlevels-1)%fe,'hcurl_pm1')
END IF
IF(order>1)CALL hcurl_span_check(ML_hcurl%levels(ML_hcurl%nlevels)%fe,ML_h1%levels(ML_h1%nlevels)%fe)
CALL geom%delete
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
! Evaluate quantity iq of all DOFs directly (out(nce,npc))
!------------------------------------------------------------------------------
SUBROUTINE direct_eval(fe,iq,cell,f,gop,cgop,out)
CLASS(oft_afem_type), INTENT(in) :: fe
INTEGER(i4), INTENT(in) :: iq,cell
REAL(r8), INTENT(in) :: f(:),gop(3,4),cgop(3,6)
REAL(r8), INTENT(out) :: out(:,:)
REAL(r8) :: rop1(fe%nce),rop3(3,fe%nce)
INTEGER(i4) :: c
SELECT TYPE(fe)
  CLASS IS(oft_scalar_fem)
    IF(iq==qcache_val)THEN
      CALL oft_lag_eval_all(fe,cell,f,rop1)
      out(:,1)=rop1
      RETURN
    END IF
    CALL oft_lag_geval_all(fe,cell,f,rop3,gop)
  CLASS IS(oft_h1_fem)
    IF(iq==qcache_val)THEN
      CALL oft_h1_eval_all(fe,cell,f,rop1)
      out(:,1)=rop1
      RETURN
    END IF
    CALL oft_h1_geval_all(fe,cell,f,rop3,gop)
  CLASS IS(oft_hcurl_fem)
    IF(iq==qcache_val)THEN
      CALL oft_hcurl_eval_all(fe,cell,f,rop3,gop)
    ELSE
      CALL oft_hcurl_ceval_all(fe,cell,f,rop3,cgop)
    END IF
END SELECT
DO c=1,3
  out(:,c)=rop3(c,:)
END DO
END SUBROUTINE direct_eval
!------------------------------------------------------------------------------
! Compare both cache modes with direct evaluation for one FE space
!------------------------------------------------------------------------------
SUBROUTINE check_space(fe_in,name)
CLASS(oft_afem_type), TARGET, INTENT(in) :: fe_in
CHARACTER(LEN=*), INTENT(in) :: name
TYPE(oft_fem_qcache) :: cache
INTEGER(i4), PARAMETER :: nmat_cells=64
INTEGER(i4) :: ipath,mode,iq,cell,m,j,c,npc,nce,nfld,io_unit,imat
INTEGER(i4), ALLOCATABLE :: rows(:),dofs(:)
REAL(r8) :: gop(3,4),cgop(3,6),det,emax(3),vmax(3),rmax(3),rscale(3),amax(3),ascale(3)
REAL(r8) :: flux(3),tmp,kmax,kscale,icmax(3),acmax(3)
REAL(r8), ALLOCATABLE :: coefc(:,:),dI(:,:,:,:),cI(:,:,:),flux3(:,:,:),res3(:,:),work(:)
REAL(r8), ALLOCATABLE :: dval(:,:),cval(:,:),coefs(:,:),dinterp(:,:),cinterp(:,:),dres(:,:),cres(:,:)
REAL(r8), ALLOCATABLE :: coef_glob(:,:),dmat(:,:),cmat(:,:)
CHARACTER(LEN=5) :: mode_name
CHARACTER(LEN=5), PARAMETER :: qnames(3)=(/'val  ','grad ','curl '/)
CHARACTER(LEN=5), PARAMETER :: pnames(3)=(/'full ','ref  ','refc '/)
LOGICAL :: has(3),canon
nce=fe_in%nce
nfld=2
!---Deterministic local coefficients
ALLOCATE(coef_glob(fe_in%ne,nfld))
DO j=1,fe_in%ne
  coef_glob(j,1)=SIN(1.37d0*REAL(j,8)+0.3d0)
  coef_glob(j,2)=COS(0.71d0*REAL(j,8))
END DO
! Paths: 1 = FULL, 2 = REFERENCE (gather rows), 3 = REFERENCE (canonical frame)
DO ipath=1,3
  mode_name=pnames(ipath)
  mode=qcache_mode_ref
  IF(ipath==1)mode=qcache_mode_full
  !---Lagrange/H1: values+gradients, H(Curl): values+curls
  IF(ipath/=3)CALL cache%setup(fe_in,geom,mode,need_val=.TRUE.,need_grad=.TRUE.,need_curl=.TRUE.)
  ! Canonical frame where available (gather path otherwise, see fem_qcache_canonical_ok)
  canon=(ipath==3).AND.fem_qcache_canonical_ok(cache)
  DO iq=1,3
    has(iq)=(cache%q(iq)%xform>=0)
  END DO
  imat=1
  IF(.NOT.has(1))imat=2
  IF(oft_env%head_proc)THEN
    WRITE(*,'(2X,A,1X,A,A,I4,A,I6,A,F10.2,A,F9.4,A)')TRIM(name),TRIM(mode_name), &
      ': nce = ',nce,', rows = ',cache%nrows,', mem = ',REAL(cache%mem(),8)/1.d6,' MB, setup = ', &
      cache%setup_time,' s'
    IF(ipath==2)WRITE(*,'(4X,A,4I6,A,ES10.2,A,2I4)')'orientation transforms (identity/signed perm/dense/none): ', &
      cache%nxf,', fit residual = ',cache%xf_resid,', max regressors/face DOFs = ',MAXVAL(cache%xf_nreg), &
      cache%ent_ndof(cache%nent-1)
  END IF
  emax=0.d0; vmax=0.d0; rmax=0.d0; rscale=0.d0; amax=0.d0; ascale=0.d0; kmax=0.d0; kscale=0.d0
  icmax=0.d0; acmax=0.d0
  !$omp parallel private(iq,cell,m,j,c,npc,rows,dofs,gop,cgop,det,flux,tmp,dval,cval,coefs, &
  !$omp dinterp,cinterp,dres,cres,dmat,cmat,coefc,dI,cI,flux3,res3,work) &
  !$omp reduction(max:emax,vmax,rmax,rscale,amax,ascale,kmax,kscale,icmax,acmax)
  ALLOCATE(rows(nce),dofs(nce),dval(nce,3),cval(nce,3),coefs(nce,nfld))
  ALLOCATE(coefc(nce,nfld),dI(3,geom%nqp,nfld,3),cI(3,geom%nqp,nfld),flux3(3,geom%nqp,1),res3(nce,1))
  ALLOCATE(work(fem_qcache_work_size(cache,nfld)))
  ALLOCATE(dinterp(3,nfld),cinterp(3,nfld),dres(nce,3),cres(nce,3),dmat(nce,nce),cmat(nce,nce))
  !$omp do schedule(dynamic,16)
  DO cell=1,mesh%nc
    CALL fem_qcache_cell_rows(cache,cell,rows)
    CALL fe_in%ncdofs(cell,dofs)
    DO j=1,nfld
      coefs(:,j)=coef_glob(dofs,j)
    END DO
    coefc=coefs
    IF(canon)CALL fem_qcache_coefs_to_ref(cache,cell,coefs,nfld)
    dres=0.d0; cres=0.d0; dmat=0.d0; cmat=0.d0
    DO m=1,geom%nqp
      gop=0.d0
      CALL mesh%jacobian(cell,geom%pts(:,m),gop,det)
      CALL oft_hcurl_get_cgops(gop,cgop)
      DO iq=1,3
        IF(.NOT.has(iq))CYCLE
        npc=cache%q(iq)%npc
        CALL direct_eval(fe_in,iq,cell,geom%pts(:,m),gop,cgop,dval)
        !---Field reconstruction
        IF(canon)THEN
          CALL fem_qcache_interp(cache,iq,cell,m,coefs,nfld,cinterp)
        ELSE
          CALL fem_qcache_interp(cache,iq,cell,m,coefs,nfld,cinterp,rows)
        END IF
        DO j=1,nfld
          DO c=1,npc
            dinterp(c,j)=SUM(coef_glob(dofs,j)*dval(:,c))
            dI(c,m,j,iq)=dinterp(c,j)
            rscale(iq)=MAX(rscale(iq),SUM(ABS(coef_glob(dofs,j)*dval(:,c))))
          END DO
        END DO
        rmax(iq)=MAX(rmax(iq),MAXVAL(ABS(cinterp(1:npc,:)-dinterp(1:npc,:))))
        !---Test-function accumulation (accumulated over qp, checked per cell)
        flux=(/0.3d0,-0.5d0,0.8d0/)*REAL(m,8)
        DO c=1,npc
          dres(:,iq)=dres(:,iq)+dval(:,c)*flux(c)
        END DO
        IF(canon)THEN
          CALL fem_qcache_accum(cache,iq,cell,m,flux,cres(:,iq))
        ELSE
          CALL fem_qcache_accum(cache,iq,cell,m,flux,cres(:,iq),rows)
        END IF
        tmp=0.d0
        DO c=1,npc
          tmp=MAX(tmp,MAXVAL(ABS(dval(:,c)*flux(c))))
        END DO
        ascale(iq)=MAX(ascale(iq),tmp*REAL(geom%nqp,8))
        !---Basis evaluation
        IF(canon)THEN
          CALL fem_qcache_eval(cache,iq,cell,m,cval)
        ELSE
          CALL fem_qcache_eval(cache,iq,cell,m,cval,rows)
        END IF
        !---Local matrix (first quantity, first cells only)
        IF(iq==imat.AND.cell<=nmat_cells)THEN
          DO c=1,npc
            DO j=1,nce
              dmat(:,j)=dmat(:,j)+dval(:,c)*dval(j,c)
              cmat(:,j)=cmat(:,j)+cval(:,c)*cval(j,c)
            END DO
          END DO
        END IF
        IF(canon)CALL fem_qcache_res_from_ref(cache,cell,cval,npc)
        emax(iq)=MAX(emax(iq),MAXVAL(ABS(cval(:,1:npc)-dval(:,1:npc))))
        vmax(iq)=MAX(vmax(iq),MAXVAL(ABS(dval(:,1:npc))))
      END DO
    END DO
    IF(canon)CALL fem_qcache_res_from_ref(cache,cell,cres,3)
    !---Cell-batched kernels (FULL and REFERENCE)
    IF(ipath<3)THEN
      DO iq=1,3
        IF(.NOT.has(iq))CYCLE
        npc=cache%q(iq)%npc
        CALL fem_qcache_interp_cell(cache,iq,cell,coefc,nfld,cI,work)
        icmax(iq)=MAX(icmax(iq),MAXVAL(ABS(cI(1:npc,:,:)-dI(1:npc,:,:,iq))))
        DO m=1,geom%nqp
          flux3(:,m,1)=(/0.3d0,-0.5d0,0.8d0/)*REAL(m,8)
        END DO
        res3=0.d0
        CALL fem_qcache_accum_cell(cache,iq,cell,flux3,1,res3,work)
        acmax(iq)=MAX(acmax(iq),MAXVAL(ABS(res3(:,1)-dres(:,iq))))
      END DO
    END IF
    DO iq=1,3
      IF(has(iq))amax(iq)=MAX(amax(iq),MAXVAL(ABS(cres(:,iq)-dres(:,iq))))
    END DO
    IF(cell<=nmat_cells)THEN
      IF(canon)CALL fem_qcache_mat_from_ref(cache,cache,cell,cmat)
      kmax=MAX(kmax,MAXVAL(ABS(cmat-dmat)))
      kscale=MAX(kscale,MAXVAL(ABS(dmat)))
    END IF
  END DO
  DEALLOCATE(rows,dofs,dval,cval,coefs,dinterp,cinterp,dres,cres,dmat,cmat,coefc,dI,cI,flux3,res3,work)
  !$omp end parallel
  !---Global reduction and output
  kmax=oft_mpi_max(kmax/kscale)
  IF(oft_env%head_proc)THEN
    WRITE(*,'(4X,A,ES12.4)')'local matrix err = ',kmax
    OPEN(NEWUNIT=io_unit,FILE='qcache.results',POSITION='APPEND')
    WRITE(io_unit,'(A,ES12.4)')TRIM(name)//'_'//TRIM(mode_name)//'_mat ',kmax
    IF(ipath==2)WRITE(io_unit,'(A,ES12.4)')TRIM(name)//'_xform_resid ',cache%xf_resid
    IF(ipath==2)WRITE(io_unit,'(A,I6)')TRIM(name)//'_xform_none ',cache%nxf(3)
    IF(ipath==3)WRITE(io_unit,'(A,I2)')TRIM(name)//'_canonical ',MERGE(1,0,canon)
    CLOSE(io_unit)
  END IF
  DO iq=1,3
    IF(.NOT.has(iq))CYCLE
    emax(iq)=oft_mpi_max(emax(iq)/vmax(iq))
    rmax(iq)=oft_mpi_max(rmax(iq)/rscale(iq))
    icmax(iq)=oft_mpi_max(icmax(iq)/rscale(iq))
    acmax(iq)=oft_mpi_max(acmax(iq)/ascale(iq))
    amax(iq)=oft_mpi_max(amax(iq)/ascale(iq))
    IF(oft_env%head_proc)THEN
      WRITE(*,'(4X,A,A,3ES12.4)')qnames(iq),' eval/interp/accum err = ',emax(iq),rmax(iq),amax(iq)
      OPEN(NEWUNIT=io_unit,FILE='qcache.results',POSITION='APPEND')
      WRITE(io_unit,'(A,ES12.4)')TRIM(name)//'_'//TRIM(mode_name)//'_'//TRIM(qnames(iq))//'_eval ',emax(iq)
      WRITE(io_unit,'(A,ES12.4)')TRIM(name)//'_'//TRIM(mode_name)//'_'//TRIM(qnames(iq))//'_interp ',rmax(iq)
      WRITE(io_unit,'(A,ES12.4)')TRIM(name)//'_'//TRIM(mode_name)//'_'//TRIM(qnames(iq))//'_accum ',amax(iq)
      IF(ipath<3)THEN
        WRITE(io_unit,'(A,ES12.4)')TRIM(name)//'_'//TRIM(mode_name)//'_'//TRIM(qnames(iq))//'_icell ',icmax(iq)
        WRITE(io_unit,'(A,ES12.4)')TRIM(name)//'_'//TRIM(mode_name)//'_'//TRIM(qnames(iq))//'_acell ',acmax(iq)
      END IF
      CLOSE(io_unit)
    END IF
  END DO
  IF(ipath/=2)CALL cache%delete
END DO
DEALLOCATE(coef_glob)
END SUBROUTINE check_space
!------------------------------------------------------------------------------
! Diagnostic: are rotated H(Curl) face functions combinations of canonical H(Curl) functions
! plus gradients of canonical H^1 functions of the same cell? Reports the maximum relative
! least-squares residual for H(Curl) values (and, with the same H(Curl) coefficients, curls).
!------------------------------------------------------------------------------
SUBROUTINE hcurl_span_check(fe_hcurl,fe_h1)
CLASS(oft_afem_type), TARGET, INTENT(in) :: fe_hcurl,fe_h1
TYPE(oft_fem_qcache) :: cc,ch
INTEGER(i4) :: i,s,j,k,m,mm,nd,nr,nm,info,lwork,nc_c,nc_h,nmc,io_unit
REAL(r8) :: rv,rc,rv_own
INTEGER(i4) :: rank
INTEGER(i4), ALLOCATABLE :: iwork(:)
REAL(r8), ALLOCATABLE :: amat(:,:),bmat(:,:),a0(:,:),b0(:,:),work(:),ac(:,:),bc(:,:),sv(:)
EXTERNAL :: dgelsd
CALL cc%setup(fe_hcurl,geom,qcache_mode_ref)
CALL ch%setup(fe_h1,geom,qcache_mode_ref,need_val=.FALSE.,need_grad=.TRUE.)
nc_c=cc%nce; nc_h=ch%nce
nm=geom%nlog*geom%nqp
nmc=geom%ncurl*geom%nqp
rv=0.d0; rc=0.d0
DO i=1,cc%nent
  nd=cc%ent_ndof(i)
  IF(nd==0.OR.cc%ent_nst(i)<=2)CYCLE ! faces only
  DO s=2,cc%ent_nst(i)
    nr=nc_c+nc_h
    ALLOCATE(amat(nm,nr),bmat(nm,nd),a0(nm,nr),b0(nm,nd),ac(nmc,nc_c),bc(nmc,nd))
    mm=0
    DO m=1,geom%nqp
      DO k=1,geom%nlog
        mm=mm+1
        amat(mm,1:nc_c)=cc%q(qcache_val)%ref(1:nc_c,k,m)
        amat(mm,nc_c+1:nr)=ch%q(qcache_grad)%ref(1:nc_h,k,m)
        DO j=1,nd
          bmat(mm,j)=cc%q(qcache_val)%ref(cc%ent_row0(s,i)+j,k,m)
        END DO
      END DO
    END DO
    mm=0
    DO m=1,geom%nqp
      DO k=1,geom%ncurl
        mm=mm+1
        ac(mm,:)=cc%q(qcache_curl)%ref(1:nc_c,k,m)
        DO j=1,nd
          bc(mm,j)=cc%q(qcache_curl)%ref(cc%ent_row0(s,i)+j,k,m)
        END DO
      END DO
    END DO
    a0=amat; b0=bmat
    ! SVD-based solve: the combined regressor set is rank deficient (eg. sum of H^1 vertex gradients)
    ALLOCATE(work(1),iwork(1),sv(nr))
    CALL dgelsd(nm,nr,nd,amat,nm,bmat,nm,sv,1.d-12,rank,work,-1,iwork,info)
    lwork=INT(work(1),4)
    k=iwork(1)
    DEALLOCATE(work,iwork); ALLOCATE(work(lwork),iwork(MAX(1,k)))
    CALL dgelsd(nm,nr,nd,amat,nm,bmat,nm,sv,1.d-12,rank,work,lwork,iwork,info)
    DEALLOCATE(work,iwork,sv)
    rv=MAX(rv,MAXVAL(ABS(b0-MATMUL(a0,bmat(1:nr,:))))/MAXVAL(ABS(a0)))
    rc=MAX(rc,MAXVAL(ABS(bc-MATMUL(ac,bmat(1:nc_c,:))))/MAXVAL(ABS(ac)))
    DEALLOCATE(amat,bmat,a0,b0,ac,bc)
  END DO
END DO
rv_own=REAL(cc%nxf(3),8)
IF(oft_env%head_proc)THEN
  WRITE(*,'(A,I6,2ES12.4)')'H(Curl) faces: untransformable states / fit residual with grad(H^1) (values, curls): ', &
    INT(rv_own),rv,rc
  OPEN(NEWUNIT=io_unit,FILE='qcache.results',POSITION='APPEND')
  WRITE(io_unit,'(A,I6)')'hcurl_span_none ',INT(rv_own)
  WRITE(io_unit,'(A,ES12.4)')'hcurl_span_with_grad_val ',rv
  WRITE(io_unit,'(A,ES12.4)')'hcurl_span_with_grad_curl ',rc
  CLOSE(io_unit)
END IF
CALL cc%delete
CALL ch%delete
END SUBROUTINE hcurl_span_check
END PROGRAM test_fem_qcache
