!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file fem_qcache_bench.F90
!
!> Kernel benchmark for the quadrature-point basis cache (\ref fem_quad_cache), using the
!! FE spaces of xMHD (Lagrange order p, H^1 order p+1, H(Curl) order p) on the finest level
!! with the finest H(Curl) quadrature rule, as \ref xmhd::xmhd_build_ops does.
!!
!! Three evaluation variants are timed over all cells and quadrature points:
!! - direct: `mesh%%jacobian` (once per straight cell, per qp otherwise) + `oft_*_eval_all` calls
!! - full: FULL mode cache
!! - ref: REFERENCE mode cache, table rows gathered per DOF
!! - refc: REFERENCE mode cache, canonical frame (per-cell orientation transform of coefficients
!!   and results, contiguous table access) for spaces where it is available, gather otherwise
!! - fullc / refcel: FULL / REFERENCE mode with the cell-batched kernels
!!   (`fem_qcache_interp_cell`, `fem_qcache_accum_cell`; basis via `fem_qcache_eval_cell`)
!!
!! each for three kernels:
!! - basis: physical basis values/gradients/curls of all three spaces into local arrays
!! - recon: basis + reconstruction of the xMHD fields from local coefficients (6 Lagrange
!!   fields with gradients, 1 grad(H^1) field, 1 H(Curl) field with curl), as in
!!   \ref xmhd::xmhd_interp_apply
!! - resid: recon + test-function accumulation of all fields (residual-like)
!!
!! Repetitions are interleaved (rep -> variant -> kernel). Output lines start with `QCACHE_BENCH`.
!! Not run by the regression suite; see `benchmarks/xmhd_qcache/run_qcache_bench.py`.
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM fem_qcache_bench
USE oft_base
USE oft_mesh_type, ONLY: oft_mesh, cell_is_curved
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
USE fem_base, ONLY: oft_ml_fem_type
USE oft_lag_basis, ONLY: oft_lag_setup, oft_scalar_fem, oft_lag_eval_all, oft_lag_geval_all
USE oft_h1_basis, ONLY: oft_h1_setup, oft_h1_fem, oft_h1_geval_all
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_fem, oft_hcurl_eval_all, &
  oft_hcurl_ceval_all, oft_hcurl_get_cgops
USE fem_quad_cache
IMPLICIT NONE
INTEGER(i4), PARAMETER :: nlf=6 ! Lagrange fields (V, N, Ti, Te)
INTEGER(i4) :: order,nreps,ierr,io_unit,r,v,k,i
INTEGER(i4) :: nce_l,nce_h,nce_c
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_lag,ML_h1,ML_hcurl
CLASS(oft_mesh), POINTER :: mesh
CLASS(oft_scalar_fem), POINTER :: lag
CLASS(oft_h1_fem), POINTER :: h1
CLASS(oft_hcurl_fem), POINTER :: hcurl
TYPE(oft_qcache_geom), TARGET :: geom
TYPE(oft_fem_qcache), TARGET :: cache(3,2) ! (lag/h1/hcurl, full/ref)
REAL(r8), ALLOCATABLE :: lcoef(:,:),hcoef(:),ccoef(:)
REAL(r8), ALLOCATABLE :: times(:,:,:),sums(:,:,:)
REAL(r8) :: elapsed,setup_t(3,2),mem(3,2)
TYPE(oft_timer) :: mytimer
CHARACTER(LEN=6), PARAMETER :: vnames(6)=(/'direct','full  ','ref   ','refc  ','fullc ','refcel'/)
CHARACTER(LEN=5), PARAMETER :: knames(3)=(/'basis','recon','resid'/)
NAMELIST/qcache_bench_options/order,nreps
!---Initialize enviroment
CALL oft_init
order=2
nreps=5
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,qcache_bench_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid and FE spaces (finest level only)
CALL multigrid_construct(mg_mesh)
mesh=>mg_mesh%mesh
CALL oft_lag_setup(mg_mesh,order,ML_lag,minlev=-1)
CALL oft_h1_setup(mg_mesh,order+1,ML_h1,minlev=-1)
CALL oft_hcurl_setup(mg_mesh,order,ML_hcurl,minlev=-1)
SELECT TYPE(this=>ML_lag%levels(ML_lag%nlevels)%fe)
  CLASS IS(oft_scalar_fem)
    lag=>this
END SELECT
SELECT TYPE(this=>ML_h1%levels(ML_h1%nlevels)%fe)
  CLASS IS(oft_h1_fem)
    h1=>this
END SELECT
SELECT TYPE(this=>ML_hcurl%levels(ML_hcurl%nlevels)%fe)
  CLASS IS(oft_hcurl_fem)
    hcurl=>this
END SELECT
nce_l=lag%nce; nce_h=h1%nce; nce_c=hcurl%nce
!---Coefficient vectors (local DOF numbering)
ALLOCATE(lcoef(lag%ne,nlf),hcoef(h1%ne),ccoef(hcurl%ne))
DO i=1,lag%ne
  DO k=1,nlf
    lcoef(i,k)=SIN(0.37d0*REAL(i,8)+REAL(k,8))
  END DO
END DO
DO i=1,h1%ne
  hcoef(i)=COS(0.53d0*REAL(i,8))
END DO
DO i=1,hcurl%ne
  ccoef(i)=SIN(0.29d0*REAL(i,8)+0.1d0)
END DO
!---Build caches
CALL geom%setup(mesh,hcurl%quad)
DO v=1,2
  CALL cache(1,v)%setup(lag,geom,v)
  CALL cache(2,v)%setup(h1,geom,v)
  CALL cache(3,v)%setup(hcurl,geom,v)
  DO k=1,3
    setup_t(k,v)=cache(k,v)%setup_time
    mem(k,v)=REAL(cache(k,v)%mem(),8)
  END DO
END DO
IF(oft_env%head_proc)THEN
  WRITE(*,'(A,2I8,I4,3I6,I8)')'QCACHE_BENCH info ncells,nqp,nthreads,nce_lag,nce_h1,nce_hcurl,nvar ', &
    mesh%nc,geom%nqp,oft_env%nthreads,nce_l,nce_h,nce_c,geom%nvar
  WRITE(*,'(A,ES12.4,I14)')'QCACHE_BENCH geom setup_s,bytes ',geom%setup_time,geom%mem()
  DO v=1,2
    WRITE(*,'(A,A,3ES12.4,3ES12.4)')'QCACHE_BENCH cache ',TRIM(vnames(v+1)),setup_t(:,v),mem(:,v)
  END DO
  DO k=1,3
    WRITE(*,'(A,I2,4I6,ES12.4)')'QCACHE_BENCH xform space,nident,nperm,ndense,nnone,resid ',k,cache(k,2)%nxf, &
      cache(k,2)%xf_resid
  END DO
END IF
!---Timed kernels (interleaved)
ALLOCATE(times(3,6,nreps),sums(3,6,nreps))
DO r=1,nreps
  DO v=1,6
    DO k=1,3
      CALL mytimer%tick
      IF(v==1)THEN
        CALL run_direct(k,sums(k,v,r))
      ELSE IF(v<=4)THEN
        CALL run_cache(MIN(v-1,2),v==4,k,sums(k,v,r))
      ELSE
        CALL run_cell(v-4,k,sums(k,v,r))
      END IF
      times(k,v,r)=mytimer%tock()
    END DO
  END DO
END DO
IF(oft_env%head_proc)THEN
  DO v=1,6
    DO k=1,3
      WRITE(*,'(A,A,1X,A,ES24.15)',ADVANCE='NO')'QCACHE_BENCH time ',TRIM(vnames(v)),knames(k),sums(k,v,1)
      DO r=1,nreps
        WRITE(*,'(ES12.4)',ADVANCE='NO')times(k,v,r)
      END DO
      WRITE(*,*)
    END DO
  END DO
END IF
DO v=1,2
  DO k=1,3
    CALL cache(k,v)%delete
  END DO
END DO
CALL geom%delete
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
! Direct evaluation (as xmhd_build_ops / xmhd_interp_apply)
!------------------------------------------------------------------------------
SUBROUTINE run_direct(kern,csum)
INTEGER(i4), INTENT(in) :: kern
REAL(r8), INTENT(out) :: csum
INTEGER(i4) :: cell,m,j,f
LOGICAL :: curved
REAL(r8) :: gop(3,4),cgop(3,6),det,detw,lval(nlf),lgrad(3,nlf),b0(3),j0(3),tmp,vtmp(3),sink
INTEGER(i4), ALLOCATABLE :: dofs(:)
REAL(r8), ALLOCATABLE :: lag_rop(:),lag_gop(:,:),h1_gop(:,:),hc_rop(:,:),hc_cop(:,:)
REAL(r8), ALLOCATABLE :: lc(:,:),hc(:),cc(:),lres(:,:),hres(:),cres(:)
sink=0.d0
!$omp parallel private(cell,m,j,f,curved,gop,cgop,det,detw,lval,lgrad,b0,j0,tmp,vtmp,dofs, &
!$omp lag_rop,lag_gop,h1_gop,hc_rop,hc_cop,lc,hc,cc,lres,hres,cres) reduction(+:sink)
ALLOCATE(dofs(MAX(nce_l,nce_h,nce_c)))
ALLOCATE(lag_rop(nce_l),lag_gop(3,nce_l),h1_gop(3,nce_h),hc_rop(3,nce_c),hc_cop(3,nce_c))
ALLOCATE(lc(nce_l,nlf),hc(nce_h),cc(nce_c),lres(nce_l,nlf),hres(nce_h),cres(nce_c))
gop=0.d0
!$omp do schedule(static)
DO cell=1,mesh%nc
  IF(kern>1)THEN
    CALL lag%ncdofs(cell,dofs)
    lc=lcoef(dofs(1:nce_l),:)
    CALL h1%ncdofs(cell,dofs)
    hc=hcoef(dofs(1:nce_h))
    CALL hcurl%ncdofs(cell,dofs)
    cc=ccoef(dofs(1:nce_c))
    lres=0.d0; hres=0.d0; cres=0.d0
  END IF
  curved=cell_is_curved(mesh,cell)
  DO m=1,geom%nqp
    IF(curved.OR.m==1)THEN
      CALL mesh%jacobian(cell,geom%pts(:,m),gop,det)
      CALL oft_hcurl_get_cgops(gop,cgop)
    END IF
    detw=det*geom%wts(m)
    CALL oft_lag_eval_all(lag,cell,geom%pts(:,m),lag_rop)
    CALL oft_lag_geval_all(lag,cell,geom%pts(:,m),lag_gop,gop)
    CALL oft_h1_geval_all(h1,cell,geom%pts(:,m),h1_gop,gop)
    CALL oft_hcurl_eval_all(hcurl,cell,geom%pts(:,m),hc_rop,gop)
    CALL oft_hcurl_ceval_all(hcurl,cell,geom%pts(:,m),hc_cop,cgop)
    IF(kern==1)THEN
      sink=sink+(lag_rop(1)+lag_gop(1,nce_l)+h1_gop(2,nce_h)+hc_rop(3,nce_c)+hc_cop(1,1))*detw
      CYCLE
    END IF
    !---Reconstruct fields
    DO f=1,nlf
      tmp=0.d0; vtmp=0.d0
      DO j=1,nce_l
        tmp=tmp+lc(j,f)*lag_rop(j)
        vtmp=vtmp+lc(j,f)*lag_gop(:,j)
      END DO
      lval(f)=tmp; lgrad(:,f)=vtmp
    END DO
    b0=0.d0
    DO j=1,nce_h
      b0=b0+hc(j)*h1_gop(:,j)
    END DO
    j0=0.d0
    DO j=1,nce_c
      b0=b0+cc(j)*hc_rop(:,j)
      j0=j0+cc(j)*hc_cop(:,j)
    END DO
    IF(kern==2)THEN
      sink=sink+(SUM(lval)+SUM(lgrad)+SUM(b0)+SUM(j0))*detw
      CYCLE
    END IF
    !---Test-function accumulation
    DO f=1,nlf
      tmp=lval(f)*detw
      vtmp=lgrad(:,f)*detw
      DO j=1,nce_l
        lres(j,f)=lres(j,f)+lag_rop(j)*tmp+DOT_PRODUCT(lag_gop(:,j),vtmp)
      END DO
    END DO
    vtmp=b0*detw
    DO j=1,nce_h
      hres(j)=hres(j)+DOT_PRODUCT(h1_gop(:,j),vtmp)
    END DO
    DO j=1,nce_c
      cres(j)=cres(j)+DOT_PRODUCT(hc_rop(:,j),vtmp)+DOT_PRODUCT(hc_cop(:,j),j0*detw)
    END DO
  END DO
  IF(kern==3)sink=sink+SUM(lres)+SUM(hres)+SUM(cres)
END DO
DEALLOCATE(dofs,lag_rop,lag_gop,h1_gop,hc_rop,hc_cop,lc,hc,cc,lres,hres,cres)
!$omp end parallel
csum=sink
END SUBROUTINE run_direct
!------------------------------------------------------------------------------
! Cached evaluation (v=1: FULL, v=2: REFERENCE; canon: canonical frame instead of row gather)
!------------------------------------------------------------------------------
SUBROUTINE run_cache(v,canon,kern,csum)
INTEGER(i4), INTENT(in) :: v,kern
LOGICAL, INTENT(in) :: canon
LOGICAL :: gl,gh,gc
REAL(r8), INTENT(out) :: csum
INTEGER(i4) :: cell,m,f
REAL(r8) :: gop(3,4),detw,lval(1,nlf),lgrad(3,nlf),bh(3,1),bc(3,1),jc(3,1),b0(3),sink,flux(3)
INTEGER(i4), ALLOCATABLE :: dofs(:),rl(:),rh(:),rc(:)
REAL(r8), ALLOCATABLE :: lag_rop(:,:),lag_gop(:,:),h1_gop(:,:),hc_rop(:,:),hc_cop(:,:)
REAL(r8), ALLOCATABLE :: lc(:,:),hc(:,:),cc(:,:),lres(:,:),hres(:,:),cres(:,:)
! Gather rows (gl/gh/gc) unless the canonical frame is requested and available for the space
gl=.NOT.(canon.AND.fem_qcache_canonical_ok(cache(1,v)))
gh=.NOT.(canon.AND.fem_qcache_canonical_ok(cache(2,v)))
gc=.NOT.(canon.AND.fem_qcache_canonical_ok(cache(3,v)))
sink=0.d0
!$omp parallel private(cell,m,f,gop,detw,lval,lgrad,bh,bc,jc,b0,flux,dofs,rl,rh,rc, &
!$omp lag_rop,lag_gop,h1_gop,hc_rop,hc_cop,lc,hc,cc,lres,hres,cres) reduction(+:sink)
ALLOCATE(dofs(MAX(nce_l,nce_h,nce_c)),rl(nce_l),rh(nce_h),rc(nce_c))
ALLOCATE(lag_rop(nce_l,1),lag_gop(nce_l,3),h1_gop(nce_h,3),hc_rop(nce_c,3),hc_cop(nce_c,3))
ALLOCATE(lc(nce_l,nlf),hc(nce_h,1),cc(nce_c,1),lres(nce_l,nlf),hres(nce_h,1),cres(nce_c,1))
!$omp do schedule(static)
DO cell=1,mesh%nc
  IF(gl)CALL fem_qcache_cell_rows(cache(1,v),cell,rl)
  IF(gh)CALL fem_qcache_cell_rows(cache(2,v),cell,rh)
  IF(gc)CALL fem_qcache_cell_rows(cache(3,v),cell,rc)
  IF(kern>1)THEN
    CALL lag%ncdofs(cell,dofs)
    lc=lcoef(dofs(1:nce_l),:)
    CALL h1%ncdofs(cell,dofs)
    hc(:,1)=hcoef(dofs(1:nce_h))
    CALL hcurl%ncdofs(cell,dofs)
    cc(:,1)=ccoef(dofs(1:nce_c))
    IF(.NOT.gl)CALL fem_qcache_coefs_to_ref(cache(1,v),cell,lc,nlf)
    IF(.NOT.gh)CALL fem_qcache_coefs_to_ref(cache(2,v),cell,hc,1)
    IF(.NOT.gc)CALL fem_qcache_coefs_to_ref(cache(3,v),cell,cc,1)
    lres=0.d0; hres=0.d0; cres=0.d0
  END IF
  DO m=1,geom%nqp
    CALL fem_qcache_geom_get(geom,cell,m,gop,detw)
    IF(kern==1)THEN
      CALL k_eval(cache(1,v),qcache_val,cell,m,lag_rop,rl,gl)
      CALL k_eval(cache(1,v),qcache_grad,cell,m,lag_gop,rl,gl)
      CALL k_eval(cache(2,v),qcache_grad,cell,m,h1_gop,rh,gh)
      CALL k_eval(cache(3,v),qcache_val,cell,m,hc_rop,rc,gc)
      CALL k_eval(cache(3,v),qcache_curl,cell,m,hc_cop,rc,gc)
      sink=sink+(lag_rop(1,1)+lag_gop(nce_l,1)+h1_gop(nce_h,2)+hc_rop(nce_c,3)+hc_cop(1,1))*detw
      CYCLE
    END IF
    !---Reconstruct fields
    CALL k_interp(cache(1,v),qcache_val,cell,m,lc,nlf,lval,rl,gl)
    CALL k_interp(cache(1,v),qcache_grad,cell,m,lc,nlf,lgrad,rl,gl)
    CALL k_interp(cache(2,v),qcache_grad,cell,m,hc,1,bh,rh,gh)
    CALL k_interp(cache(3,v),qcache_val,cell,m,cc,1,bc,rc,gc)
    CALL k_interp(cache(3,v),qcache_curl,cell,m,cc,1,jc,rc,gc)
    b0=bh(:,1)+bc(:,1)
    IF(kern==2)THEN
      sink=sink+(SUM(lval)+SUM(lgrad)+SUM(b0)+SUM(jc))*detw
      CYCLE
    END IF
    !---Test-function accumulation
    DO f=1,nlf
      flux(1)=lval(1,f)*detw
      CALL k_accum(cache(1,v),qcache_val,cell,m,flux,lres(:,f),rl,gl)
      flux=lgrad(:,f)*detw
      CALL k_accum(cache(1,v),qcache_grad,cell,m,flux,lres(:,f),rl,gl)
    END DO
    flux=b0*detw
    CALL k_accum(cache(2,v),qcache_grad,cell,m,flux,hres(:,1),rh,gh)
    CALL k_accum(cache(3,v),qcache_val,cell,m,flux,cres(:,1),rc,gc)
    flux=jc(:,1)*detw
    CALL k_accum(cache(3,v),qcache_curl,cell,m,flux,cres(:,1),rc,gc)
  END DO
  IF(kern==3)THEN
    IF(.NOT.gl)CALL fem_qcache_res_from_ref(cache(1,v),cell,lres,nlf)
    IF(.NOT.gh)CALL fem_qcache_res_from_ref(cache(2,v),cell,hres,1)
    IF(.NOT.gc)CALL fem_qcache_res_from_ref(cache(3,v),cell,cres,1)
    sink=sink+SUM(lres)+SUM(hres)+SUM(cres)
  END IF
END DO
DEALLOCATE(dofs,rl,rh,rc,lag_rop,lag_gop,h1_gop,hc_rop,hc_cop,lc,hc,cc,lres,hres,cres)
!$omp end parallel
csum=sink
END SUBROUTINE run_cache
!------------------------------------------------------------------------------
! Cell-batched cached evaluation (v=1: FULL, v=2: REFERENCE)
!------------------------------------------------------------------------------
SUBROUTINE run_cell(v,kern,csum)
INTEGER(i4), INTENT(in) :: v,kern
REAL(r8), INTENT(out) :: csum
INTEGER(i4) :: cell,m,f,nq,nw
REAL(r8) :: gop(3,4),detw,sink
INTEGER(i4), ALLOCATABLE :: dofs(:),rl(:),rh(:),rc(:)
REAL(r8), ALLOCATABLE :: lag_rop(:,:,:),lag_gop(:,:,:),h1_gop(:,:,:),hc_rop(:,:,:),hc_cop(:,:,:)
REAL(r8), ALLOCATABLE :: lc(:,:),hc(:,:),cc(:,:),lres(:,:),hres(:,:),cres(:,:),work(:)
REAL(r8), ALLOCATABLE :: lval(:,:,:),lgrad(:,:,:),bh(:,:,:),bc(:,:,:),jc(:,:,:),dw(:)
REAL(r8), ALLOCATABLE :: fl(:,:,:),fg(:,:,:),fb(:,:,:),fj(:,:,:)
nq=geom%nqp
nw=MAX(fem_qcache_work_size(cache(1,v),nlf),fem_qcache_work_size(cache(2,v),1),fem_qcache_work_size(cache(3,v),1))
sink=0.d0
!$omp parallel private(cell,m,f,gop,detw,dofs,rl,rh,rc,lag_rop,lag_gop,h1_gop,hc_rop,hc_cop, &
!$omp lc,hc,cc,lres,hres,cres,work,lval,lgrad,bh,bc,jc,dw,fl,fg,fb,fj) reduction(+:sink)
ALLOCATE(dofs(MAX(nce_l,nce_h,nce_c)),rl(nce_l),rh(nce_h),rc(nce_c),work(nw),dw(nq))
ALLOCATE(lc(nce_l,nlf),hc(nce_h,1),cc(nce_c,1),lres(nce_l,nlf),hres(nce_h,1),cres(nce_c,1))
IF(kern==1)THEN
  ALLOCATE(lag_rop(nce_l,1,nq),lag_gop(nce_l,3,nq),h1_gop(nce_h,3,nq),hc_rop(nce_c,3,nq),hc_cop(nce_c,3,nq))
ELSE
  ALLOCATE(lval(1,nq,nlf),lgrad(3,nq,nlf),bh(3,nq,1),bc(3,nq,1),jc(3,nq,1))
  ALLOCATE(fl(1,nq,nlf),fg(3,nq,nlf),fb(3,nq,1),fj(3,nq,1))
END IF
!$omp do schedule(static)
DO cell=1,mesh%nc
  DO m=1,nq
    CALL fem_qcache_geom_get(geom,cell,m,gop,dw(m))
  END DO
  IF(kern==1)THEN
    CALL fem_qcache_cell_rows(cache(1,v),cell,rl)
    CALL fem_qcache_cell_rows(cache(2,v),cell,rh)
    CALL fem_qcache_cell_rows(cache(3,v),cell,rc)
    CALL fem_qcache_eval_cell(cache(1,v),qcache_val,cell,lag_rop,rl)
    CALL fem_qcache_eval_cell(cache(1,v),qcache_grad,cell,lag_gop,rl)
    CALL fem_qcache_eval_cell(cache(2,v),qcache_grad,cell,h1_gop,rh)
    CALL fem_qcache_eval_cell(cache(3,v),qcache_val,cell,hc_rop,rc)
    CALL fem_qcache_eval_cell(cache(3,v),qcache_curl,cell,hc_cop,rc)
    DO m=1,nq
      sink=sink+(lag_rop(1,1,m)+lag_gop(nce_l,1,m)+h1_gop(nce_h,2,m)+hc_rop(nce_c,3,m)+hc_cop(1,1,m))*dw(m)
    END DO
    CYCLE
  END IF
  CALL lag%ncdofs(cell,dofs)
  lc=lcoef(dofs(1:nce_l),:)
  CALL h1%ncdofs(cell,dofs)
  hc(:,1)=hcoef(dofs(1:nce_h))
  CALL hcurl%ncdofs(cell,dofs)
  cc(:,1)=ccoef(dofs(1:nce_c))
  !---Reconstruct fields at all qp
  CALL fem_qcache_interp_cell(cache(1,v),qcache_val,cell,lc,nlf,lval,work)
  CALL fem_qcache_interp_cell(cache(1,v),qcache_grad,cell,lc,nlf,lgrad,work)
  CALL fem_qcache_interp_cell(cache(2,v),qcache_grad,cell,hc,1,bh,work)
  CALL fem_qcache_interp_cell(cache(3,v),qcache_val,cell,cc,1,bc,work)
  CALL fem_qcache_interp_cell(cache(3,v),qcache_curl,cell,cc,1,jc,work)
  IF(kern==2)THEN
    DO m=1,nq
      sink=sink+(SUM(lval(:,m,:))+SUM(lgrad(:,m,:))+SUM(bh(:,m,1)+bc(:,m,1))+SUM(jc(:,m,1)))*dw(m)
    END DO
    CYCLE
  END IF
  !---Pointwise fluxes, then test-function accumulation for all qp
  DO m=1,nq
    fl(:,m,:)=lval(:,m,:)*dw(m)
    fg(:,m,:)=lgrad(:,m,:)*dw(m)
    fb(:,m,1)=(bh(:,m,1)+bc(:,m,1))*dw(m)
    fj(:,m,1)=jc(:,m,1)*dw(m)
  END DO
  lres=0.d0; hres=0.d0; cres=0.d0
  CALL fem_qcache_accum_cell(cache(1,v),qcache_val,cell,fl,nlf,lres,work)
  CALL fem_qcache_accum_cell(cache(1,v),qcache_grad,cell,fg,nlf,lres,work)
  CALL fem_qcache_accum_cell(cache(2,v),qcache_grad,cell,fb,1,hres,work)
  CALL fem_qcache_accum_cell(cache(3,v),qcache_val,cell,fb,1,cres,work)
  CALL fem_qcache_accum_cell(cache(3,v),qcache_curl,cell,fj,1,cres,work)
  sink=sink+SUM(lres)+SUM(hres)+SUM(cres)
END DO
DEALLOCATE(dofs,rl,rh,rc,work,dw,lc,hc,cc,lres,hres,cres)
IF(kern==1)THEN
  DEALLOCATE(lag_rop,lag_gop,h1_gop,hc_rop,hc_cop)
ELSE
  DEALLOCATE(lval,lgrad,bh,bc,jc,fl,fg,fb,fj)
END IF
!$omp end parallel
csum=sink
END SUBROUTINE run_cell
!------------------------------------------------------------------------------
! Kernel wrappers: with (gather) or without (canonical frame) table rows
!------------------------------------------------------------------------------
! Evaluate with or without gathered rows
SUBROUTINE k_eval(c,iq,cell,m,out,rows,use_rows)
TYPE(oft_fem_qcache), INTENT(in) :: c
INTEGER(i4), INTENT(in) :: iq,cell,m
REAL(r8), CONTIGUOUS, INTENT(out) :: out(:,:)
INTEGER(i4), INTENT(in) :: rows(:)
LOGICAL, INTENT(in) :: use_rows
IF(use_rows)THEN
  CALL fem_qcache_eval(c,iq,cell,m,out,rows)
ELSE
  CALL fem_qcache_eval(c,iq,cell,m,out)
END IF
END SUBROUTINE k_eval
! Reconstruct with or without gathered rows
SUBROUTINE k_interp(c,iq,cell,m,coefs,nfld,out,rows,use_rows)
TYPE(oft_fem_qcache), INTENT(in) :: c
INTEGER(i4), INTENT(in) :: iq,nfld,cell,m
REAL(r8), CONTIGUOUS, INTENT(in) :: coefs(:,:)
REAL(r8), INTENT(out) :: out(:,:)
INTEGER(i4), INTENT(in) :: rows(:)
LOGICAL, INTENT(in) :: use_rows
IF(use_rows)THEN
  CALL fem_qcache_interp(c,iq,cell,m,coefs,nfld,out,rows)
ELSE
  CALL fem_qcache_interp(c,iq,cell,m,coefs,nfld,out)
END IF
END SUBROUTINE k_interp
! Accumulate with or without gathered rows
SUBROUTINE k_accum(c,iq,cell,m,flux,res,rows,use_rows)
TYPE(oft_fem_qcache), INTENT(in) :: c
INTEGER(i4), INTENT(in) :: iq,cell,m
REAL(r8), INTENT(in) :: flux(:)
REAL(r8), CONTIGUOUS, INTENT(inout) :: res(:)
INTEGER(i4), INTENT(in) :: rows(:)
LOGICAL, INTENT(in) :: use_rows
IF(use_rows)THEN
  CALL fem_qcache_accum(c,iq,cell,m,flux,res,rows)
ELSE
  CALL fem_qcache_accum(c,iq,cell,m,flux,res)
END IF
END SUBROUTINE k_accum
END PROGRAM fem_qcache_bench
