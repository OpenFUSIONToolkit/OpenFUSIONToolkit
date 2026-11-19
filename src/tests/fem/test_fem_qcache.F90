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
!! at every point as the reference geometry. Also checked: field reconstruction from local
!! coefficients, test-function accumulation, det*w and physical quadrature points.
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
INTEGER(i4) :: mode,iq,cell,m,j,c,npc,nce,nfld,io_unit
INTEGER(i4), ALLOCATABLE :: rows(:),dofs(:)
REAL(r8) :: gop(3,4),cgop(3,6),det,emax(3),vmax(3),rmax(3),rscale(3),amax(3),ascale(3)
REAL(r8) :: flux(3),tmp
REAL(r8), ALLOCATABLE :: dval(:,:),cval(:,:),coefs(:,:),dinterp(:,:),cinterp(:,:),dres(:),cres(:)
REAL(r8), ALLOCATABLE :: coef_glob(:,:)
CHARACTER(LEN=5) :: mode_name
CHARACTER(LEN=5), PARAMETER :: qnames(3)=(/'val  ','grad ','curl '/)
LOGICAL :: has(3)
nce=fe_in%nce
nfld=2
!---Deterministic local coefficients
ALLOCATE(coef_glob(fe_in%ne,nfld))
DO j=1,fe_in%ne
  coef_glob(j,1)=SIN(1.37d0*REAL(j,8)+0.3d0)
  coef_glob(j,2)=COS(0.71d0*REAL(j,8))
END DO
DO mode=qcache_mode_full,qcache_mode_ref
  IF(mode==qcache_mode_full)THEN
    mode_name='full'
  ELSE
    mode_name='ref'
  END IF
  !---Lagrange/H1: values+gradients, H(Curl): values+curls
  CALL cache%setup(fe_in,geom,mode,need_val=.TRUE.,need_grad=.TRUE.,need_curl=.TRUE.)
  DO iq=1,3
    has(iq)=(cache%q(iq)%xform>=0)
  END DO
  IF(oft_env%head_proc)WRITE(*,'(2X,A,1X,A,A,I4,A,I6,A,F10.2,A,F9.4,A)')TRIM(name),TRIM(mode_name), &
    ': nce = ',nce,', rows = ',cache%nrows,', mem = ',REAL(cache%mem(),8)/1.d6,' MB, setup = ', &
    cache%setup_time,' s'
  emax=0.d0; vmax=0.d0; rmax=0.d0; rscale=0.d0; amax=0.d0; ascale=0.d0
  !$omp parallel private(iq,cell,m,j,c,npc,rows,dofs,gop,cgop,det,flux,tmp,dval,cval,coefs, &
  !$omp dinterp,cinterp,dres,cres) reduction(max:emax,vmax,rmax,rscale,amax,ascale)
  ALLOCATE(rows(nce),dofs(nce),dval(nce,3),cval(nce,3),coefs(nce,nfld))
  ALLOCATE(dinterp(3,nfld),cinterp(3,nfld),dres(nce),cres(nce))
  !$omp do schedule(dynamic,16)
  DO cell=1,mesh%nc
    CALL fem_qcache_cell_rows(cache,cell,rows)
    CALL fe_in%ncdofs(cell,dofs)
    DO j=1,nfld
      coefs(:,j)=coef_glob(dofs,j)
    END DO
    DO m=1,geom%nqp
      gop=0.d0
      CALL mesh%jacobian(cell,geom%pts(:,m),gop,det)
      CALL oft_hcurl_get_cgops(gop,cgop)
      DO iq=1,3
        IF(.NOT.has(iq))CYCLE
        npc=cache%q(iq)%npc
        CALL direct_eval(fe_in,iq,cell,geom%pts(:,m),gop,cgop,dval)
        CALL fem_qcache_eval(cache,iq,cell,m,rows,cval)
        emax(iq)=MAX(emax(iq),MAXVAL(ABS(cval(:,1:npc)-dval(:,1:npc))))
        vmax(iq)=MAX(vmax(iq),MAXVAL(ABS(dval(:,1:npc))))
        !---Field reconstruction
        CALL fem_qcache_interp(cache,iq,cell,m,rows,coefs,nfld,cinterp)
        DO j=1,nfld
          DO c=1,npc
            dinterp(c,j)=SUM(coefs(:,j)*dval(:,c))
            rscale(iq)=MAX(rscale(iq),SUM(ABS(coefs(:,j)*dval(:,c))))
          END DO
        END DO
        rmax(iq)=MAX(rmax(iq),MAXVAL(ABS(cinterp(1:npc,:)-dinterp(1:npc,:))))
        !---Test-function accumulation
        flux=(/0.3d0,-0.5d0,0.8d0/)*REAL(m,8)
        dres=0.d0
        DO c=1,npc
          dres=dres+dval(:,c)*flux(c)
        END DO
        cres=0.d0
        CALL fem_qcache_accum(cache,iq,cell,m,rows,flux,cres)
        amax(iq)=MAX(amax(iq),MAXVAL(ABS(cres-dres)))
        tmp=0.d0
        DO c=1,npc
          tmp=MAX(tmp,MAXVAL(ABS(dval(:,c)*flux(c))))
        END DO
        ascale(iq)=MAX(ascale(iq),tmp)
      END DO
    END DO
  END DO
  DEALLOCATE(rows,dofs,dval,cval,coefs,dinterp,cinterp,dres,cres)
  !$omp end parallel
  !---Global reduction and output
  DO iq=1,3
    IF(.NOT.has(iq))CYCLE
    emax(iq)=oft_mpi_max(emax(iq)/vmax(iq))
    rmax(iq)=oft_mpi_max(rmax(iq)/rscale(iq))
    amax(iq)=oft_mpi_max(amax(iq)/ascale(iq))
    IF(oft_env%head_proc)THEN
      WRITE(*,'(4X,A,A,3ES12.4)')qnames(iq),' eval/interp/accum err = ',emax(iq),rmax(iq),amax(iq)
      OPEN(NEWUNIT=io_unit,FILE='qcache.results',POSITION='APPEND')
      WRITE(io_unit,'(A,ES12.4)')TRIM(name)//'_'//TRIM(mode_name)//'_'//TRIM(qnames(iq))//'_eval ',emax(iq)
      WRITE(io_unit,'(A,ES12.4)')TRIM(name)//'_'//TRIM(mode_name)//'_'//TRIM(qnames(iq))//'_interp ',rmax(iq)
      WRITE(io_unit,'(A,ES12.4)')TRIM(name)//'_'//TRIM(mode_name)//'_'//TRIM(qnames(iq))//'_accum ',amax(iq)
      CLOSE(io_unit)
    END IF
  END DO
  CALL cache%delete
END DO
DEALLOCATE(coef_glob)
END SUBROUTINE check_space
END PROGRAM test_fem_qcache
