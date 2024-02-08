!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file thincurr_fr.F90
!
!> Run thin wall frequency response simulations using ThinCurr
!!
!! **Option group:** `thincurr_fr_options`
!! |  Option                 |  Description  | Type [dim] |
!! |-------------------------|---------------|------------|
!! |  `mesh_file="none"`     |  Surface mesh filename (Cubit) | str |
!! |  `plot_run=F`           |  Produce plot files from stored restart files | bool |
!! |  `direct=T`             |  Use direct solver | bool |
!! |  `force_f0=0.`          |  Toroidal field (R0*B0) for force calculation | float |
!! |  `freq=1.d3`            |  Frequency for FR run | float |
!! |  `mode_file="none"`     |  DCON mode surface file from "mode_to_tw" | str |
!! |  `fr_limit=0`           |  Frequency limit for FR run (1->Inf, 2->Zero) | str |
!!
!! @authors Chris Hansen
!! @date Feb 2022
!! @ingroup doxy_thincurr
!---------------------------------------------------------------------------
PROGRAM thincurr_fr
USE oft_base
USE oft_io, ONLY: hdf5_create_timestep
USE oft_mesh_type, ONLY: smesh
USE oft_mesh_native, ONLY: native_read_nodesets, native_read_sidesets
#ifdef HAVE_NCDF
USE oft_mesh_cubit, ONLY: smesh_cubit_load, cubit_read_nodesets, cubit_read_sidesets
#endif
USE multigrid_build, ONLY: multigrid_construct_surf
!
USE oft_la_base, ONLY: oft_vector, oft_graph
USE oft_lu, ONLY: oft_lusolver, lapack_matinv
USE oft_native_la, ONLY: oft_native_vector, oft_native_matrix, partition_graph
USE oft_deriv_matrices, ONLY: oft_sum_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_gmres_solver, create_diag_pre, &
  create_native_solver
#ifdef HAVE_ARPACK
USE oft_arpack, ONLY: oft_iram_eigsolver
#endif
USE axi_green, ONLY: green
USE mhd_utils, ONLY: mu0
USE thin_wall
IMPLICIT NONE
INTEGER(4) :: nsensors = 0
TYPE(tw_type) :: tw_sim,mode_source
TYPE(tw_sensors) :: sensors
!
INTEGER(4) :: i,n,ierr,io_unit
REAL(8), POINTER, contiguous, DIMENSION(:,:) :: mode_driver,fr_driver,fr_sensor
REAL(8), ALLOCATABLE :: eig_rval(:),eig_ival(:),eig_vec(:,:)
CHARACTER(LEN=2) :: eig_tag
TYPE(oft_timer) :: mytimer
CLASS(oft_vector), POINTER :: uio
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_nsets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_ssets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: hole_nsets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: jumper_nsets => NULL()
!
INTEGER(4) :: fr_limit = 0
INTEGER(4) :: jumper_start = -1
LOGICAL :: direct = .TRUE.
LOGICAL :: save_L = .FALSE.
LOGICAL :: save_Mcoil = .FALSE.
LOGICAL :: save_Msen = .FALSE.
REAL(8) :: freq = 1.d3
REAL(8) :: force_f0 = 0.d0
CHARACTER(LEN=80) :: mode_file = 'none'
NAMELIST/thincurr_fr_options/direct,save_L,save_Mcoil,save_Msen,freq,mode_file,fr_limit,jumper_start
!---
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit,thincurr_fr_options, IOSTAT=ierr)
CLOSE(io_unit)
if(ierr<0)call oft_abort('No thin-wall options found in input file.', &
  'thincurr_fr',__FILE__)
if(ierr>0)call oft_abort('Error parsing thin-wall options in input file.', &
  'thincurr_fr',__FILE__)
!---Setup mesh
CALL multigrid_construct_surf
! ALLOCATE(mg_mesh)
! mg_mesh%mgmax=1
! mg_mesh%nbase=1
! oft_env%nbase=1
! mg_mesh%mgdim=mg_mesh%mgmax
! CALL smesh_cubit_load
SELECT CASE(smesh%cad_type)
CASE(0)
  CALL native_read_nodesets(mesh_nsets)
  CALL native_read_sidesets(mesh_ssets)
CASE(2)
#ifdef HAVE_NCDF
  CALL cubit_read_nodesets(mesh_nsets)
  CALL cubit_read_sidesets(mesh_ssets)
#endif
CASE DEFAULT
  CALL oft_abort("Unsupported mesh type","thincurr_fr",__FILE__)
END SELECT
IF(ASSOCIATED(mesh_ssets))THEN
  IF(mesh_ssets(1)%n>0)THEN
    tw_sim%nclosures=mesh_ssets(1)%n
    ALLOCATE(tw_sim%closures(tw_sim%nclosures))
    tw_sim%closures=mesh_ssets(1)%v
  END IF
END IF
tw_sim%mesh=>smesh
IF(jumper_start>0)THEN
  n=SIZE(mesh_nsets)
  hole_nsets=>mesh_nsets(1:jumper_start-1)
  jumper_nsets=>mesh_nsets(jumper_start:n)
ELSE
  hole_nsets=>mesh_nsets
END IF
CALL tw_sim%setup(hole_nsets)
!---Setup I/0
CALL smesh%setup_io(1)
IF(oft_debug_print(1))CALL tw_sim%save_debug()
!---------------------------------------------------------------------------
! Frequency-response run
!---------------------------------------------------------------------------
!---Load sensors
CALL tw_load_sensors(tw_sim,sensors,jumper_nsets)
!---Setup voltage source
ALLOCATE(fr_driver(tw_sim%nelems,2),fr_sensor(sensors%nfloops,2))
fr_driver=0.d0
fr_sensor=0.d0
IF(TRIM(mode_file)/="none")THEN ! DCON mode source
  !---Load DCON mesh
  CALL tw_load_mode(mode_file,mode_source,mode_driver)
  !---Compute voltage source term from DCON mode
  CALL tw_compute_Lmat_MF(mode_source,tw_sim,2,mode_driver,fr_driver)
  !---Compute mode-sensor coupling
  CALL tw_compute_mutuals(mode_source,sensors%nfloops,sensors%floops)
  CALL dgemv('N',sensors%nfloops,mode_source%nelems,1.d0,mode_source%Ael2sen, &
    sensors%nfloops,mode_driver(:,1),1,0.d0,fr_sensor(:,1),1)
  CALL dgemv('N',sensors%nfloops,mode_source%nelems,1.d0,mode_source%Ael2sen, &
    sensors%nfloops,mode_driver(:,2),1,0.d0,fr_sensor(:,2),1)
  DEALLOCATE(mode_source%Adr2sen,mode_source%Ael2sen,mode_driver)
  !---Compute sensor mutuals for main model
  CALL tw_compute_mutuals(tw_sim,sensors%nfloops,sensors%floops)
ELSE ! Driver coil as source
  IF(tw_sim%n_icoils==0)CALL oft_abort('No coils specified for FR run', &
    'thincurr_fr', __FILE__)
  IF(tw_sim%n_icoils>2)CALL oft_abort('More than two coils specified for FR run', &
    'thincurr_fr', __FILE__)
  !---------------------------------------------------------------------------
  ! Load or build coil to element mutual matrix
  !---------------------------------------------------------------------------
  IF(save_Mcoil)THEN
    CALL tw_compute_Ael2dr(tw_sim,'Mcoil.save')
  ELSE
    CALL tw_compute_Ael2dr(tw_sim)
  END IF
  !---------------------------------------------------------------------------
  ! Load or build sensor mutual matrices
  !---------------------------------------------------------------------------
  IF(save_Msen)THEN
    CALL tw_compute_mutuals(tw_sim,sensors%nfloops,sensors%floops,'Msen.save')
  ELSE
    CALL tw_compute_mutuals(tw_sim,sensors%nfloops,sensors%floops)
  END IF
  !---FR with coil 1
  fr_driver(:,1)=tw_sim%Ael2dr(:,1)
  fr_sensor(:,1)=tw_sim%Adr2sen(:,1)
  IF(tw_sim%n_icoils>1)THEN ! Second coil if present
    fr_driver(:,2)=tw_sim%Ael2dr(:,2)
    fr_sensor(:,2)=tw_sim%Adr2sen(:,2)
  END IF
END IF
!---------------------------------------------------------------------------
! Load or build element to element mutual matrix
!---------------------------------------------------------------------------
IF(save_L)THEN
  CALL tw_compute_LmatDirect(tw_sim,tw_sim,tw_sim%Lmat,'Lmat.save')
ELSE
  CALL tw_compute_LmatDirect(tw_sim,tw_sim,tw_sim%Lmat)
END IF
!---Setup resistivity matrix
CALL tw_compute_Rmat(tw_sim,.TRUE.)
!---Compute Frequency-response
CALL run_fr(tw_sim, freq, fr_driver, fr_sensor, fr_limit)
!---
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE run_fr
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE run_fr(self,freq,driver,senin,fr_limit)
TYPE(tw_type), INTENT(in) :: self
REAL(r8), INTENT(in) :: freq
REAL(r8), INTENT(inout) :: driver(:,:)
REAL(r8), INTENT(inout) :: senin(:,:)
INTEGER(i4), INTENT(in) :: fr_limit
INTEGER(i4) :: i,j,k,io_unit,info
REAL(r8) :: lam0
REAL(r8), ALLOCATABLE :: senout(:,:)
DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:) :: x,b
DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: Mmat
CLASS(oft_vector), POINTER :: uloc
TYPE(oft_graph) :: graph
LOGICAL :: pm_save
!---
WRITE(*,*)
WRITE(*,*)'Starting Frequency-response run'
!---Setup matrix
ALLOCATE(Mmat(self%nelems,self%nelems),b(self%nelems))
b=-(0.d0,1.d0)*driver(:,1) + (1.d0,0.d0)*driver(:,2) ! -i*L_e*I_e
SELECT CASE(fr_limit)
  CASE(0)
    WRITE(*,'(X,A,ES13.5)')'  Frequency [Hz] = ',freq
    b=b*freq*2.d0*pi ! Scale RHS by forcing frequency
    Mmat=(0.d0,1.d0)*freq*2.d0*pi*self%Lmat
    DO i=1,self%Rmat%nr
      DO j=self%Rmat%kr(i),self%Rmat%kr(i+1)-1
        Mmat(i,self%Rmat%lc(j))=Mmat(i,self%Rmat%lc(j)) + self%Rmat%M(j)
      END DO
    END DO
  CASE(1)
    WRITE(*,'(X,A)')'  Frequency -> Inf (L limit)'
    Mmat=(0.d0,1.d0)*self%Lmat
  CASE(2)
    WRITE(*,'(X,A)')'  Frequency -> 0   (R limit)'
    Mmat=(0.d0,0.d0)
    DO i=1,self%Rmat%nr
      DO j=self%Rmat%kr(i),self%Rmat%kr(i+1)-1
        Mmat(i,self%Rmat%lc(j))=Mmat(i,self%Rmat%lc(j)) + self%Rmat%M(j)
      END DO
    END DO
  CASE DEFAULT
    CALL oft_abort('Invalid "fr_limit" value (0,1,2)', 'run_fr', __FILE__)
END SELECT
!---Solve system
ALLOCATE(x(self%nelems))
x=(0.d0,0.d0)
IF(direct)THEN
  CALL lapack_matinv(self%nelems,Mmat,info)
  CALL zgemv('N',self%nelems,self%nelems,1.d0,Mmat,self%nelems,b,1,0.d0,x,1)
ELSE
  graph%nr=self%Rmat%nr
  graph%nnz=self%Rmat%nnz
  graph%kr=>self%Rmat%kr
  graph%lc=>self%Rmat%lc
  CALL gmres_comp(60,-1,self%nelems,Mmat,x,b,graph,MAX(1,INT(self%nelems/2000,4)),.TRUE.)
END IF
DEALLOCATE(Mmat,b)

!---Setup local vectors for I/O
CALL self%Uloc%new(uloc)
!---Save real part
driver(:,1)=REAL(x,8)
CALL uloc%restore_local(driver(:,1))
CALL tw_save_pfield(self,driver(:,1),'JRe')
CALL tw_rst_save(self,uloc,'pThinCurr_frRe.rst','U')
!---Save imaginary part
driver(:,2)=AIMAG(x)
CALL uloc%restore_local(driver(:,2))
CALL tw_save_pfield(self,driver(:,2),'JIm')
CALL tw_rst_save(self,uloc,'pThinCurr_frIm.rst','U')
!---Save probe signals
IF(sensors%nfloops>0)THEN
  ALLOCATE(senout(sensors%nfloops,2))
  senout=senin
  CALL dgemv('N',sensors%nfloops,self%nelems,1.d0,self%Ael2sen,sensors%nfloops,driver(:,1),1,1.d0,senout(:,1),1)
  CALL dgemv('N',sensors%nfloops,self%nelems,1.d0,self%Ael2sen,sensors%nfloops,driver(:,2),1,1.d0,senout(:,2),1)
  OPEN(NEWUNIT=io_unit,FILE='thincurr_fr.dat')
  DO i=1,sensors%nfloops
    WRITE(io_unit,'(A,4Es24.15)')sensors%floops(i)%name,senout(i,:),senin(i,:)
  END DO
  CLOSE(io_unit)
  DEALLOCATE(senout)
END IF

!---Cleanup
CALL uloc%delete
DEALLOCATE(x,uloc)
END SUBROUTINE run_fr
!------------------------------------------------------------------------------
! FUNCTION: vec_comp_dot
!------------------------------------------------------------------------------
!> Dot product with a second vector
!!
!! @param[in] a Second vector for dot product
!! @result \f$ \sum_i self_i a_i \f$
!------------------------------------------------------------------------------
function vec_comp_dot(n,a,b) result(dot)
integer(i4), intent(in) :: n
DOUBLE COMPLEX, intent(in) :: a(n),b(n)
INTEGER(i4) :: i
DOUBLE COMPLEX :: dot
dot=(0.d0,0.d0)
!$omp parallel do reduction(+:dot)
DO i=1,n
  dot=dot+CONJG(a(i))*b(i)
END DO
end function vec_comp_dot
!------------------------------------------------------------------------------
! SUBROUTINE: mat_comp
!------------------------------------------------------------------------------
!> Dot product with a second vector
!!
!! @param[in] a Second vector for dot product
!! @result \f$ \sum_i self_i a_i \f$
!------------------------------------------------------------------------------
subroutine mat_comp(n,A,x,b)
integer(i4), intent(in) :: n
DOUBLE COMPLEX, intent(in) :: A(n,n),x(n)
DOUBLE COMPLEX, intent(out) :: b(n)
integer(i4) :: i,j,i1,i2
b=(0.d0,0.d0)
IF(omp_in_parallel())THEN
  DO j=1,n
    !$omp simd
    DO i=1,n
      b(i)=b(i)+A(i,j)*x(j)
    END DO
  END DO
ELSE
  !$omp parallel private(i1,i2,i,j)
  CALL oft_thread_slice(oft_tid,oft_env%nthreads,n,i1,i2)
  DO j=1,n
    !$omp simd
    DO i=i1,i2
      b(i)=b(i)+A(i,j)*x(j)
    END DO
  END DO
  !$omp end parallel
END IF
end subroutine mat_comp
!------------------------------------------------------------------------------
! SUBROUTINE gmres_comp
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE gmres_comp(nrits,its,ncoils,A,u,g,graph,nlocal,loverlap)
TYPE :: mat_parts
  INTEGER(i4) :: n = 0 !< Number of values in set
  INTEGER(i4), POINTER, DIMENSION(:) :: ind => NULL() !< Indices
  DOUBLE COMPLEX, POINTER, CONTIGUOUS, DIMENSION(:) :: x => NULL() !< Indices
  DOUBLE COMPLEX, POINTER, CONTIGUOUS, DIMENSION(:) :: b => NULL() !< Indices
  DOUBLE COMPLEX, POINTER, CONTIGUOUS, DIMENSION(:,:) :: mat => NULL() !< Indices
END TYPE mat_parts
INTEGER(i4), INTENT(in) :: nrits,its,ncoils,nlocal
DOUBLE COMPLEX, INTENT(inout) :: A(ncoils,ncoils)
DOUBLE COMPLEX, CONTIGUOUS, intent(inout) :: u(:),g(:)
TYPE(oft_graph), INTENT(inout) :: graph
LOGICAL, INTENT(in) :: loverlap
DOUBLE COMPLEX, allocatable, dimension(:) :: r,w
DOUBLE COMPLEX, allocatable, dimension(:,:) :: v,z
DOUBLE COMPLEX, allocatable :: h(:,:),c(:),s(:),res(:)
DOUBLE COMPLEX :: delta,hkmi,ggin
REAL(r8) :: uu,uuold,gg,ggold,elapsed_time
integer(i4) :: i,j,jr,k,kk,nits,info
integer(i4), allocatable, dimension(:) :: part
LOGICAL, allocatable, dimension(:) :: eflag
TYPE(mat_parts), pointer, dimension(:) :: parts,parts_tmp
LOGICAL :: pm
TYPE(oft_timer) :: stimer
pm=oft_env%pm
oft_env%pm=.FALSE.
!---Preconditioner setup
IF(nlocal>1)THEN
  WRITE(*,"(A,I6)")" nLocal = ",nlocal
  ALLOCATE(part(ncoils))
  CALL partition_graph(graph,nlocal,part)
  !---Create partition mappings
  ALLOCATE(parts(nlocal))
  DO i=1,ncoils
    IF(part(i)<0)CYCLE
    parts(part(i))%n=parts(part(i))%n+1
  END DO
  DO i=1,nlocal
    ALLOCATE(parts(i)%ind(parts(i)%n))
    ! ALLOCATE(parts(i)%x(parts(i)%n))
    ! ALLOCATE(parts(i)%b(parts(i)%n))
    ! ALLOCATE(parts(i)%mat(parts(i)%n,parts(i)%n))
    parts(i)%n=0
  END DO
  DO i=1,ncoils
    IF(part(i)<0)CYCLE
    parts(part(i))%n=parts(part(i))%n+1
    parts(part(i))%ind(parts(part(i))%n) = i
  END DO
  IF(loverlap)THEN
    parts_tmp=>parts
    ALLOCATE(parts(nlocal))
    !$omp parallel private(j,jr,k,eflag)
    ALLOCATE(eflag(ncoils))
    !$omp do
    DO i=1,nlocal
      eflag=(part==i)
      DO j=1,parts_tmp(i)%n
        jr=parts_tmp(i)%ind(j)
        DO k=graph%kr(jr),graph%kr(jr+1)-1
          eflag(graph%lc(k))=.TRUE.
        END DO
      END DO
      DO j=1,ncoils
        IF(part(j)<0)eflag(j)=.FALSE.
      END DO
      parts(i)%n=COUNT(eflag)
      ALLOCATE(parts(i)%ind(parts(i)%n))
      k=0
      DO j=1,graph%nr
        IF(.NOT.eflag(j))CYCLE
        k=k+1
        parts(i)%ind(k)=j
      END DO
      DEALLOCATE(parts_tmp(i)%ind)
    END DO
    DEALLOCATE(eflag)
    !$omp end parallel
    DEALLOCATE(parts_tmp)
  END IF
  DO i=1,nlocal
    ! ALLOCATE(parts(i)%ind(parts(i)%n))
    ALLOCATE(parts(i)%x(parts(i)%n))
    ALLOCATE(parts(i)%b(parts(i)%n))
    ALLOCATE(parts(i)%mat(parts(i)%n,parts(i)%n))
    ! parts(i)%n=0
  END DO
  !---Get matrix slice
  !$omp parallel do private(j,k) schedule(static,1)
  DO i=1,nlocal
    DO j=1,parts(i)%n
      DO k=1,parts(i)%n
        parts(i)%mat(j,k)=A(parts(i)%ind(j),parts(i)%ind(k))
      END DO
    END DO
    CALL lapack_matinv(parts(i)%n,parts(i)%mat,info)
  END DO
  IF(loverlap)THEN
    !$omp parallel do private(j)
    DO i=1,nlocal
      DO j=1,parts(i)%n
        IF(part(parts(i)%ind(j))/=i)parts(i)%ind(j)=-parts(i)%ind(j)
      END DO
    END DO
  END IF
  DEALLOCATE(part)
END IF
!---
IF(pm)THEN
  WRITE(*,*)'Starting GMRES solver'
  CALL stimer%tick
END IF
!---
ALLOCATE(v(ncoils,nrits+1),z(ncoils,nrits+1))
ALLOCATE(r(ncoils),w(ncoils))
w=(0.d0,0.d0)
allocate(c(nrits+1),s(nrits+1),res(nrits+1),h(nrits+1,nrits))
c=(0.d0,0.d0); s=(0.d0,0.d0); res=(0.d0,0.d0); h=(0.d0,0.d0)
!
uu = REAL(vec_comp_dot(ncoils,u,u),8)
IF(uu>0.d0)CALL mat_comp(ncoils,A,u,w)
r=g-w
gg = REAL(vec_comp_dot(ncoils,r,r),8)
ggin=gg
100 FORMAT (I6,2ES14.6)
110 FORMAT (I6,3ES14.6)
IF(pm)WRITE(*,100)0,SQRT(uu),SQRT(gg)
nits=0
do j=1,ncoils
  IF(gg==0.d0)EXIT
  res=(0.d0,0.d0)
  res(1)=SQRT(gg)
  v(:,1)=r/res(1)
  do i=1,nrits
    !---Precondition search direction
    IF(nlocal>1)THEN
      !$omp parallel do private(kk) schedule(static,1)
      DO k=1,nlocal
        DO kk=1,parts(k)%n
          parts(k)%b(kk)=v(ABS(parts(k)%ind(kk)),i)
          parts(k)%x(kk)=(0.d0,0.d0)
        END DO
        CALL mat_comp(parts(k)%n,parts(k)%mat,parts(k)%b,parts(k)%x)
        DO kk=1,parts(k)%n
          IF(parts(k)%ind(kk)>0)z(parts(k)%ind(kk),i)=parts(k)%x(kk)
        END DO
      END DO
      ! CALL gmres_comp(4,4,ncoils,A,z(:,i),w,1)
    ELSE
      do k=1,ncoils
        z(k,i)=v(k,i)/A(k,k)
      end do
    END IF
    CALL mat_comp(ncoils,A,z(:,i),w)
    !---Arnoldi iteration
    do k=1,i
      h(k,i)=vec_comp_dot(ncoils,w,v(:,k))
      w=w-h(k,i)*v(:,k)
    end do
    h(i+1,i)=SQRT(vec_comp_dot(ncoils,w,w))
    v(:,i+1)=w/h(i+1,i)
    !---Apply Givens rotation
    do k=2,i
      hkmi = h(k-1,i)
      h(k-1,i) = c(k-1)*hkmi + s(k-1)*h(k,i)
      h(k,i) = -s(k-1)*hkmi + c(k-1)*h(k,i)
    enddo
    delta = SQRT( ABS(h(i,i))**2 + ABS(h(i+1,i))**2 )
    c(i) = h(i,i) / delta
    s(i) = h(i+1,i) / delta
    h(i,i) = c(i)*h(i,i) + s(i)*h(i+1,i)
    !---Update residual
    res(i+1) = -s(i)*res(i)
    res(i) = c(i)*res(i)
    nits=nits+1
    IF(i==nrits)EXIT
    IF(nits==its)EXIT
    IF(pm)WRITE(*,'(I6,14X,ES14.6)')nits,ABS(res(i+1))
  end do
  k=i
  do i=k,2,-1
    res(i) = res(i) / h(i,i)
    res(1:i-1) = res(1:i-1) - res(i)*h(1:i-1,i)
  end do
  res(1) = res(1) / h(1,1)
  ! Update iterate
  do i=1,k
    u=u+res(i)*z(:,i)
  enddo
  CALL mat_comp(ncoils,A,u,w)
  r=g-w
  IF(nits==its)EXIT
  ggold=gg; uuold=uu
  gg = REAL(vec_comp_dot(ncoils,r,r),8)
  uu = REAL(vec_comp_dot(ncoils,u,u),8)
  IF(pm)WRITE(*,110)nits,SQRT(uu),SQRT(gg),SQRT(gg/uu)
  IF(its==-1.AND.sqrt(gg/uu)<1.d-10)EXIT
  IF(its==-2.AND.sqrt(gg/uu)<1.d-15)EXIT
  IF(nits==its)EXIT
  IF(uu==uuold)EXIT
end do
g=r
DEALLOCATE(c,s,h,res)
DEALLOCATE(v,z,r,w)
!
IF(pm)THEN
  elapsed_time=stimer%tock()
  WRITE(*,*)'  Time = ',elapsed_time
END IF
IF(nlocal>1)THEN
  DO i=1,nlocal
    DEALLOCATE(parts(i)%ind)
    DEALLOCATE(parts(i)%x)
    DEALLOCATE(parts(i)%b)
    DEALLOCATE(parts(i)%mat)
  END DO
  DEALLOCATE(parts)
END IF
oft_env%pm=pm
END SUBROUTINE gmres_comp
END PROGRAM thincurr_fr
