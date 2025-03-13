!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file thin_wall_solvers.F90
!
!> Module for thin-wall modeling on 3D triangular meshes
!!
!!
!! @authors Chris Hansen
!! @date May 2017
!! @ingroup doxy_oft_physics
!---------------------------------------------------------------------------
MODULE thin_wall_solvers
USE oft_base
USE oft_sort, ONLY: sort_array
USE oft_io, ONLY: oft_bin_file, hdf5_add_string_attribute
!
USE oft_la_base, ONLY: oft_vector, oft_cvector, oft_matrix, oft_graph
USE oft_lu, ONLY: oft_lusolver, lapack_matinv, lapack_cholesky
USE oft_native_la, ONLY: oft_native_dense_matrix, partition_graph
USE oft_deriv_matrices, ONLY: oft_sum_matrix, oft_sum_cmatrix
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_native_gmres_csolver
USE oft_solver_utils, ONLY: create_diag_pre, create_native_solver, create_cg_solver
#ifdef HAVE_ARPACK
USE oft_arpack, ONLY: oft_irlm_eigsolver
#endif
USE mhd_utils, ONLY: mu0
USE thin_wall
USE thin_wall_hodlr, ONLY: oft_tw_hodlr_op, oft_tw_hodlr_bjpre, oft_tw_hodlr_rbjpre
IMPLICIT NONE
#include "local.h"
CONTAINS
!------------------------------------------------------------------------------
!> Compute L/R eigenmodes of ThinCurr model using a direct approach via LAPACK
!------------------------------------------------------------------------------
SUBROUTINE lr_eigenmodes_direct(self,neigs,eig_rval,eig_vec,eig_ival)
TYPE(tw_type), INTENT(in) :: self !< ThinCurr object
INTEGER(4), INTENT(in) :: neigs !< Number of eigenvalues to compute
REAL(8), INTENT(out) :: eig_rval(:) !< Real part of eigenvalues
REAL(8), INTENT(out) :: eig_vec(:,:) !< Eigenvectors [self%nelems,neigs]
REAL(8), OPTIONAL, INTENT(out) :: eig_ival(:) !< Imaginary part of eigenvalues (should be zero)
!---
INTEGER(i4) :: i,j,info,N,LDVL,LDVR,LDA,LWORK
REAL(r8) :: elapsed_time,chk_var
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: ipiv,sort_tmp
REAL(r8), ALLOCATABLE, DIMENSION(:) :: etatmp,WR,WI,WORK,W_tmp
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: Mmat,Acomp,VL,VR
CHARACTER(LEN=1) :: JOBVL,JOBVR
TYPE(oft_timer) :: loctimer
DEBUG_STACK_PUSH
!---Invert Mmat
ALLOCATE(Mmat(self%nelems,self%nelems))
Mmat=0.d0
DO i=1,self%Rmat%nr
  DO j=self%Rmat%kr(i),self%Rmat%kr(i+1)-1
    Mmat(i,self%Rmat%lc(j))=self%Rmat%M(j)
  END DO
END DO
! CALL lapack_matinv(self%nelems,Mmat,info)
CALL lapack_cholesky(self%nelems,Mmat,info)
!---
ALLOCATE(Acomp(self%nelems,self%nelems))
N = self%nelems
CALL dgemm('N','N',N,N,N,1.d0,Mmat,N,self%Lmat,N,0.d0,Acomp,N)
DEALLOCATE(Mmat)
!---Compute eigenvalues
JOBVL = 'N'
JOBVR = 'V'
N = self%nelems
LDA = N
LDVL = N
LDVR = N
WRITE(*,*)'Starting eigenvalue solve'
CALL loctimer%tick
ALLOCATE(WR(N),WI(N),VL(LDVL,N),VR(LDVR,N),WORK(1))
LWORK = -1
CALL DGEEV(JOBVL, JOBVR, N, Acomp, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
LWORK = INT(WORK(1),4)
IF(oft_debug_print(1))WRITE(*,*)'  Block size = ',lwork/N
DEALLOCATE(work)
ALLOCATE(work(lwork))
CALL DGEEV(JOBVL, JOBVR, N, Acomp, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
DEALLOCATE(VL,WORK,Acomp)
elapsed_time=loctimer%tock()
WRITE(*,*)'  Time = ',elapsed_time,INFO
!---Sort eigenvalues
ALLOCATE(sort_tmp(N))
sort_tmp=[(N+1-i,i=1,N)]
ALLOCATE(W_tmp(N))
W_tmp=ABS(WR)
sort_tmp=[(i,i=1,N)]
CALL sort_array(W_tmp,sort_tmp,N)
DEALLOCATE(W_tmp)
!---Copy to output
DO i=1,neigs
  j = sort_tmp(N+1-i)
  eig_rval(i)=WR(j)
  IF(PRESENT(eig_ival))eig_ival(i)=WI(j)
   eig_vec(:,i)=REAL(VR(:,j),8)
END DO
DEALLOCATE(WR,WI,VR,sort_tmp)
WRITE(*,*)'Eigenvalues'
DO i=1,MIN(neigs,5)
  WRITE(*,*)'  ',eig_rval(i)
END DO
DEBUG_STACK_POP
END SUBROUTINE lr_eigenmodes_direct
!------------------------------------------------------------------------------
!> Compute L/R eigenmodes of ThinCurr model using an iterative Lanczos method via ARPACK
!------------------------------------------------------------------------------
SUBROUTINE lr_eigenmodes_arpack(self,neigs,eig_rval,eig_vec,eig_ival,hodlr_op)
TYPE(tw_type), INTENT(in) :: self !< ThinCurr object
INTEGER(4), INTENT(in) :: neigs !< Number of eigenvalues to compute
REAL(8), INTENT(out) :: eig_rval(:) !< Real part of eigenvalues
REAL(8), INTENT(out) :: eig_vec(:,:) !< Eigenvectors [self%nelems,neigs]
REAL(8), OPTIONAL, INTENT(out) :: eig_ival(:) !< Imaginary part of eigenvalues (should be zero)
TYPE(oft_tw_hodlr_op), TARGET, OPTIONAL, INTENT(inout) :: hodlr_op !< HODLR L matrix
#ifdef HAVE_ARPACK
!---
INTEGER(4) :: i,j,k
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: sort_tmp
REAL(8) :: lam0,elapsed_time
REAL(8), ALLOCATABLE, DIMENSION(:) :: W_tmp
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: Mmat
LOGICAL :: pm_save
CLASS(oft_vector), POINTER :: uloc
TYPE(oft_native_dense_matrix), TARGET :: fmat,bmat
CLASS(oft_solver), POINTER :: linv
TYPE(oft_irlm_eigsolver) :: arsolver
TYPE(oft_timer) :: loctimer
DEBUG_STACK_PUSH
WRITE(*,*)'Starting eigenvalue solve'
!---Setup local vectors
CALL self%Uloc%new(uloc)
CALL self%Rmat%assemble(uloc)

!---Setup R^-1 solver
NULLIFY(linv)
CALL create_native_solver(linv,'lu')
SELECT TYPE(this=>linv)
CLASS IS(oft_lusolver)
  IF(this%package(1:4)=='none')THEN
    WRITE(*,*)'  Note: LU solver not available, falling back to CG'
    CALL linv%delete
    DEALLOCATE(linv)
  END IF
CLASS DEFAULT
  CALL oft_abort('Unable to allocate LU solver', 'run_eig', __FILE__)
END SELECT
IF(.NOT.ASSOCIATED(linv))THEN
  CALL create_native_solver(linv,'cg')
  linv%its=-2
  CALL create_diag_pre(linv%pre)
END IF
linv%A=>self%Rmat

!---Setup Arnoldi eig value solver
IF(PRESENT(hodlr_op))THEN
  arsolver%A=>hodlr_op
ELSE
  fmat%nr=self%nelems; fmat%nc=self%nelems
  fmat%nrg=self%nelems; fmat%ncg=self%nelems
  fmat%M => self%Lmat
  arsolver%A=>fmat
END IF
arsolver%M=>self%Rmat
arsolver%tol=1.E-5_r8
arsolver%Minv=>linv
arsolver%nev=neigs
arsolver%mode=2

!---Compute eigenvalues/vectors
CALL loctimer%tick
pm_save=oft_env%pm; oft_env%pm=.false.
CALL arsolver%apply(uloc,lam0)
oft_env%pm=pm_save
elapsed_time=loctimer%tock()
WRITE(*,*)'  Time = ',elapsed_time
IF(arsolver%info>=0)THEN
  !---Sort eigenvalues
  ALLOCATE(sort_tmp(arsolver%nev),W_tmp(arsolver%nev))
  W_tmp=ABS(arsolver%eig_val(1,1:arsolver%nev))
  sort_tmp=[(i,i=1,arsolver%nev)]
  CALL sort_array(W_tmp,sort_tmp,arsolver%nev)
  DEALLOCATE(W_tmp)
  !---Copy output
  DO i=1,neigs
    j = sort_tmp(arsolver%nev+1-i)
    eig_rval(i)=arsolver%eig_val(1,j)
    IF(PRESENT(eig_ival))eig_ival(i)=arsolver%eig_val(2,j)
    eig_vec(:,i)=arsolver%eig_vec(:,j)
  END DO
  DEALLOCATE(sort_tmp)
  ! !---Handle closures
  ! DO i=1,self%nclosures
  !   eig_vec(self%closures(i),1:neigs)=0.d0
  ! END DO
  WRITE(*,*)'Eigenvalues'
  DO i=1,MIN(neigs,5)
    WRITE(*,*)'  ',eig_rval(i)
  END DO
END IF
!---Cleanup
CALL uloc%delete()
CALL arsolver%delete()
DEALLOCATE(uloc)
#else
CALL oft_abort("Iterative eigenvalue solve requires ARPACK", "lr_eigenmodes_arpack", __FILE__)
#endif
DEBUG_STACK_POP
END SUBROUTINE lr_eigenmodes_arpack
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE frequency_response(self,direct,fr_limit,freq,driver,hodlr_op)
TYPE(tw_type), INTENT(in) :: self !< ThinCurr object
LOGICAL, INTENT(in) :: direct !< Use direct solver?
REAL(r8), INTENT(in) :: freq !< Frequency for calculation (if `fr_limit=0`)
INTEGER(i4), INTENT(in) :: fr_limit !< Limit for frequency response (0: Use `freq`, 1: Inductive, 2: Resistive)
REAL(r8), INTENT(inout) :: driver(:,:) !< Driver voltages (real, imaginary) [self%nelems,2]
TYPE(oft_tw_hodlr_op), TARGET, OPTIONAL, INTENT(inout) :: hodlr_op !< HODLR L matrix
INTEGER(i4) :: i,j,k,io_unit,info
REAL(r8) :: lam0
REAL(r8), ALLOCATABLE :: senout(:,:)
DOUBLE COMPLEX, POINTER, DIMENSION(:) :: x,b
DOUBLE COMPLEX, POINTER, DIMENSION(:,:) :: Mmat
TYPE(oft_sum_cmatrix), TARGET :: aca_Fmat
CLASS(oft_vector), POINTER :: uloc
CLASS(oft_cvector), POINTER :: aloc,bloc
TYPE(oft_graph) :: graph
TYPE(oft_native_gmres_csolver), TARGET :: frinv
TYPE(oft_tw_hodlr_bjpre), TARGET :: frinv_pre
TYPE(oft_bin_file) :: floop_hist
LOGICAL :: pm_save
!---
WRITE(*,*)
WRITE(*,*)'Starting Frequency-response run'
!---Setup matrix
ALLOCATE(b(self%nelems))
b=(1.d0,0.d0)*driver(:,1) + (0.d0,1.d0)*driver(:,2)
IF(PRESENT(hodlr_op))THEN
  SELECT CASE(fr_limit)
  CASE(0)
    WRITE(*,'(X,A,ES13.5)')'  Frequency [Hz] = ',freq
    aca_Fmat%rJ=>hodlr_op
    aca_Fmat%rK=>self%Rmat
    aca_Fmat%beta=(0.d0,1.d0)*freq*2.d0*pi
    aca_Fmat%alam=(1.d0,0.d0)
  CASE(1)
    WRITE(*,'(X,A)')'  Frequency -> Inf (L limit)'
    aca_Fmat%rJ=>hodlr_op
    aca_Fmat%rK=>self%Rmat
    aca_Fmat%beta=(0.d0,1.d0)
    aca_Fmat%alam=(0.d0,0.d0)
  CASE(2)
    WRITE(*,'(X,A)')'  Frequency -> 0   (R limit)'
    aca_Fmat%rJ=>self%Rmat
    aca_Fmat%rK=>self%Rmat
    aca_Fmat%beta=(0.d0,0.d0)
    aca_Fmat%alam=(1.d0,0.d0)
  CASE DEFAULT
    CALL oft_abort('Invalid "fr_limit" value (0,1,2)', 'run_fr', __FILE__)
  END SELECT
ELSE
  ALLOCATE(Mmat(self%nelems,self%nelems))
  SELECT CASE(fr_limit)
  CASE(0)
    WRITE(*,'(X,A,ES13.5)')'  Frequency [Hz] = ',freq
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
END IF
!---Solve system
ALLOCATE(x(self%nelems))
x=(0.d0,0.d0)
IF(direct)THEN
  CALL lapack_matinv(self%nelems,Mmat,info)
  CALL zgemv('N',self%nelems,self%nelems,(1.d0,0.d0),Mmat,self%nelems,b,1,(0.d0,0.d0),x,1)
ELSE
  IF(PRESENT(hodlr_op))THEN
    frinv%pre=>frinv_pre
    frinv_pre%mf_obj=>hodlr_op
    frinv_pre%Rmat=>self%Rmat
    frinv_pre%alpha=aca_Fmat%beta
    frinv_pre%beta=aca_Fmat%alam
    ! ALLOCATE(oft_diag_cscale::frinv%pre)
    frinv%A=>aca_Fmat
    frinv%nrits=60
    frinv%its=-1
    NULLIFY(aloc,bloc)
    CALL self%Uloc%new(aloc)
    CALL self%Uloc%new(bloc)
    CALL aca_Fmat%assemble(aloc)
    !---Solve system
    CALL aloc%restore_local(x)
    CALL bloc%restore_local(b)
    CALL frinv%apply(aloc,bloc)
    CALL aloc%get_local(x)
    !---Cleanup temporaries
    CALL aloc%delete()
    CALL bloc%delete()
    DEALLOCATE(aloc,bloc)
  ELSE
    graph%nr=self%Rmat%nr
    graph%nnz=self%Rmat%nnz
    graph%kr=>self%Rmat%kr
    graph%lc=>self%Rmat%lc
    CALL gmres_comp(60,-1,self%nelems,Mmat,x,b,graph,MAX(1,INT(self%nelems/2000,4)),.TRUE.)
    DEALLOCATE(Mmat)
  END IF
END IF
driver(:,1)=REAL(x,8)
driver(:,2)=AIMAG(x)
DEALLOCATE(b,x)
CONTAINS
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
!> Needs Docs
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE gmres_comp(nrits,its,ncoils,A,u,g,graph,nlocal,loverlap)
INTEGER(i4), INTENT(in) :: nrits,its,ncoils,nlocal
DOUBLE COMPLEX, INTENT(inout) :: A(ncoils,ncoils)
DOUBLE COMPLEX, CONTIGUOUS, intent(inout) :: u(:),g(:)
TYPE(oft_graph), INTENT(inout) :: graph
LOGICAL, INTENT(in) :: loverlap
!
TYPE :: mat_parts
  INTEGER(i4) :: n = 0 !< Number of values in set
  INTEGER(i4), POINTER, DIMENSION(:) :: ind => NULL() !< Indices
  DOUBLE COMPLEX, POINTER, CONTIGUOUS, DIMENSION(:) :: x => NULL() !< Indices
  DOUBLE COMPLEX, POINTER, CONTIGUOUS, DIMENSION(:) :: b => NULL() !< Indices
  DOUBLE COMPLEX, POINTER, CONTIGUOUS, DIMENSION(:,:) :: mat => NULL() !< Indices
END TYPE mat_parts
!
DOUBLE COMPLEX, allocatable, dimension(:) :: r,w
DOUBLE COMPLEX, allocatable, dimension(:,:) :: v,z
DOUBLE COMPLEX, allocatable :: h(:,:),c(:),s(:),res(:)
DOUBLE COMPLEX :: delta,hkmi
REAL(r8) :: uu,uuold,gg,ggold,ggin,elapsed_time
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
DO j=1,ncoils
  IF(gg==0.d0)EXIT
  res=(0.d0,0.d0)
  res(1)=SQRT(gg)
  v(:,1)=r/res(1)
  DO i=1,nrits
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
      DO k=1,ncoils
        z(k,i)=v(k,i)/A(k,k)
      END DO
    END IF
    CALL mat_comp(ncoils,A,z(:,i),w)
    !---Arnoldi iteration
    DO k=1,i
      h(k,i)=vec_comp_dot(ncoils,w,v(:,k))
      w=w-h(k,i)*v(:,k)
    END DO
    h(i+1,i)=SQRT(vec_comp_dot(ncoils,w,w))
    v(:,i+1)=w/h(i+1,i)
    !---Apply Givens rotation
    DO k=2,i
      hkmi = h(k-1,i)
      h(k-1,i) = c(k-1)*hkmi + s(k-1)*h(k,i)
      h(k,i) = -s(k-1)*hkmi + c(k-1)*h(k,i)
    END DO
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
    ! IF(pm)WRITE(*,'(I6,14X,ES14.6)')nits,ABS(res(i+1))
  END DO
  k=i
  DO i=k,2,-1
    res(i) = res(i) / h(i,i)
    res(1:i-1) = res(1:i-1) - res(i)*h(1:i-1,i)
  END DO
  res(1) = res(1) / h(1,1)
  ! Update iterate
  DO i=1,k
    u=u+res(i)*z(:,i)
  END DO
  CALL mat_comp(ncoils,A,u,w)
  r=g-w
  IF(nits==its)EXIT
  ggold=gg; uuold=uu
  gg = REAL(vec_comp_dot(ncoils,r,r),8)
  uu = REAL(vec_comp_dot(ncoils,u,u),8)
  IF(pm)WRITE(*,110)nits,SQRT(uu),SQRT(gg),SQRT(gg/uu)
  IF(its==-1.AND.sqrt(gg/uu)<1.d-7)EXIT
  IF(its==-2.AND.sqrt(gg/uu)<1.d-14)EXIT
  IF(nits==its)EXIT
  IF(uu==uuold)EXIT
END DO
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
END SUBROUTINE frequency_response
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE run_td_sim(self,dt,nsteps,vec,direct,lin_tol,use_cn,nstatus,nplot,sensors,curr_waveform,volt_waveform,sensor_vals,hodlr_op)
TYPE(tw_type), INTENT(in) :: self
REAL(8), INTENT(in) :: dt
INTEGER(4), INTENT(in) :: nsteps
REAL(8), INTENT(inout) :: vec(:)
LOGICAL, INTENT(in) :: direct
REAL(8), INTENT(in) :: lin_tol
LOGICAL, INTENT(in) :: use_cn
INTEGER(4), INTENT(in) :: nstatus
INTEGER(4), INTENT(in) :: nplot
TYPE(tw_sensors), INTENT(in) :: sensors
REAL(8), POINTER, INTENT(in) :: curr_waveform(:,:)
REAL(8), POINTER, INTENT(in) :: volt_waveform(:,:)
REAL(8), POINTER, INTENT(in) :: sensor_vals(:,:)
TYPE(oft_tw_hodlr_op), TARGET, OPTIONAL, INTENT(inout) :: hodlr_op !< HODLR L matrix
!---
INTEGER(4) :: i,j,k,ntimes_curr,ntimes_volt,ncols,itime,io_unit,neta,face,info,ind1,nits,int_inds(2)
REAL(8) :: uu,t,tmp,area,p2,p1,val_prev,dt_op,int_facs(2)
REAL(8), ALLOCATABLE, DIMENSION(:) :: icoil_curr,icoil_dcurr,pcoil_volt,senout,jumpout,eta_check
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: cc_vals
REAL(8), POINTER, DIMENSION(:) :: vals
CLASS(oft_vector), POINTER :: u,g,up
CLASS(oft_matrix), POINTER :: Lmat
TYPE(oft_native_dense_matrix), TARGET :: Lmat_dense,Minv
TYPE(oft_sum_matrix), TARGET :: fmat,bmat
CLASS(oft_solver), POINTER :: linv
TYPE(oft_tw_hodlr_rbjpre), TARGET :: linv_pre
TYPE(oft_bin_file) :: floop_hist,jumper_hist
LOGICAL :: exists,volt_full
CHARACTER(LEN=4) :: pltnum
CHARACTER(LEN=15) :: fmt_str
WRITE(*,*)
WRITE(*,*)'Starting simulation'
!---Setup coil waveform
IF(ASSOCIATED(curr_waveform))THEN
  ncols=SIZE(curr_waveform,DIM=2)
  IF(ncols-1/=self%n_icoils)CALL oft_abort('# of currents in waveform does not match # of icoils', &
    'run_td_sim',__FILE__)
  ntimes_curr=SIZE(curr_waveform,DIM=1)
  ALLOCATE(icoil_curr(ncols-1),icoil_dcurr(ncols-1))
ELSE
  ALLOCATE(icoil_curr(self%n_icoils),icoil_dcurr(self%n_icoils))
  DO j=1,self%n_icoils
    icoil_curr=0.d0
    icoil_dcurr=0.d0
  END DO
  ntimes_curr=0
END IF
!---Setup voltage waveform
volt_full=.FALSE.
IF(ASSOCIATED(volt_waveform))THEN
  ncols=SIZE(volt_waveform,DIM=2)
  IF((ncols-1/=self%n_vcoils).AND.(ncols-1/=self%nelems))THEN
    CALL oft_abort('# of voltages in waveform does not match # of vcoils or model size','run_td_sim',__FILE__)
  END IF
  IF(ncols-1==self%nelems)THEN
    ALLOCATE(pcoil_volt(self%nelems))
    volt_full=.TRUE.
  ELSE
    ALLOCATE(pcoil_volt(self%n_vcoils))
  END IF
  ntimes_volt=SIZE(volt_waveform,DIM=1)
ELSE
  ALLOCATE(pcoil_volt(self%n_vcoils))
  DO j=1,self%n_vcoils
    pcoil_volt=0.d0
  END DO
  ntimes_volt=0
END IF
!---
CALL self%Uloc%new(u)
CALL self%Uloc%new(up)
CALL self%Uloc%new(g)
!---Setup inductance matrix wrapper
IF(PRESENT(hodlr_op))THEN
  Lmat=>hodlr_op
ELSE
  Lmat_dense%nr=self%nelems; Lmat_dense%nc=self%nelems
  Lmat_dense%nrg=self%nelems; Lmat_dense%ncg=self%nelems
  Lmat_dense%M=>self%Lmat
  Lmat=>Lmat_dense
END IF
!---Setup backward matrix wrapper
IF(use_cn)THEN
  bmat%nr=self%nelems; bmat%nc=self%nelems
  bmat%nrg=self%nelems; bmat%ncg=self%nelems
  bmat%J=>Lmat
  bmat%K=>self%Rmat
  bmat%alam = -dt/2.d0 ! L - (dt/2)*R
  dt_op = dt/2.d0
ELSE
  dt_op = dt
END IF
!---Setup timestep solvers
IF(.NOT.direct)THEN
  !---Setup forward matrix wrapper
  fmat%nr=self%nelems; fmat%nc=self%nelems
  fmat%nrg=self%nelems; fmat%ncg=self%nelems
  fmat%J=>Lmat
  fmat%K=>self%Rmat
  fmat%alam = dt_op ! L + dt_op*R
  CALL fmat%assemble(u)
  !---Create CG solver for forward matrix
  CALL create_cg_solver(linv)
  linv%A=>fmat
  linv%its=-2
  linv%atol=lin_tol
  IF(PRESENT(hodlr_op))THEN
    linv_pre%mf_obj=>hodlr_op
    linv_pre%Rmat=>self%Rmat
    linv_pre%alpha=1.d0
    linv_pre%beta=fmat%alam
    linv%pre=>linv_pre
  ELSE
    CALL create_diag_pre(linv%pre)
  END IF
ELSE
  !---Setup dense matrix for inverse
  Minv%nr=self%nelems; Minv%nc=self%nelems
  Minv%nrg=self%nelems; Minv%ncg=self%nelems
  ALLOCATE(Minv%M(Minv%nr,Minv%nr))
  !---Build forward matrix [L + dt_op*R] (overwritten with inverse)
  Minv%M=self%Lmat
  DO i=1,self%Rmat%nr
    DO j=self%Rmat%kr(i),self%Rmat%kr(i+1)-1
      Minv%M(i,self%Rmat%lc(j))=Minv%M(i,self%Rmat%lc(j)) + dt_op*self%Rmat%M(j)
    END DO
  END DO
  WRITE(*,*)'Starting factorization'
  oft_env%pm=.TRUE.
  ! CALL lapack_matinv(Minv%nr,Minv%M,info)
  CALL lapack_cholesky(Minv%nr,Minv%M,info)
  oft_env%pm=.FALSE.
END IF
!---
ALLOCATE(vals(self%nelems))
vals=vec
CALL u%restore_local(vals)
t=0.d0
IF(ntimes_curr>0)THEN
  DO j=1,self%n_icoils
    icoil_curr(j)=linterp(curr_waveform(:,1),curr_waveform(:,j+1),ntimes_curr,t,1)
  END DO
END IF
!---
WRITE(pltnum,'(I4.4)')0
CALL tw_rst_save(self,u,'pThinCurr_'//pltnum//'.rst','U')
CALL hdf5_write(t,'pThinCurr_'//pltnum//'.rst','time')
CALL hdf5_write(icoil_curr,'pThinCurr_'//pltnum//'.rst','coil_currents')
!---Save sensor data for t=0
CALL u%get_local(vals)
IF(sensors%nfloops>0)THEN
  WRITE(fmt_str,'(I6)')sensors%nfloops+1
  fmt_str='('//TRIM(fmt_str)//'E24.15)'
  ALLOCATE(senout(sensors%nfloops+1))
  senout=0.d0
  IF(ntimes_curr>0)CALL dgemv('N',sensors%nfloops,self%n_icoils,1.d0,self%Adr2sen, &
    sensors%nfloops,icoil_curr,1,0.d0,senout(2),1)
  IF((ntimes_volt>0).AND.ASSOCIATED(sensor_vals))THEN
    CALL linterp_facs(sensor_vals(:,1),ntimes_volt,t,int_inds,int_facs,1)
    DO i=1,sensors%nfloops
      senout(i+1) = senout(i+1) + sensor_vals(int_inds(1),i+1)*int_facs(1) &
        + sensor_vals(int_inds(2),i+1)*int_facs(2)
    END DO
  END IF
  !---Setup history file
  IF(oft_env%head_proc)THEN
    floop_hist%filedesc = 'ThinCurr flux loop history file'
    CALL floop_hist%setup('floops.hist')
    CALL floop_hist%add_field('time', 'r8', desc="Simulation time [s]")
    DO i=1,sensors%nfloops
      CALL floop_hist%add_field(sensors%floops(i)%name, 'r8')
    END DO
    CALL floop_hist%write_header
    CALL floop_hist%open
    senout(1)=t
    CALL floop_hist%write(data_r8=senout)
  END IF
END IF
IF(sensors%njumpers>0)THEN
  ALLOCATE(jumpout(sensors%njumpers+1))
  DO j=1,sensors%njumpers
    tmp=0.d0
    val_prev=0.d0
    ind1=self%pmap(sensors%jumpers(j)%points(1))
    IF(ind1>0)val_prev=vals(ind1)
    DO k=1,sensors%jumpers(j)%np-1
      ind1=self%pmap(sensors%jumpers(j)%points(k+1))
      IF(ind1>0)THEN
        tmp=tmp+vals(ind1)-val_prev
        val_prev=vals(ind1)
      ELSE
        tmp=tmp-val_prev
        val_prev=0.d0
      END IF
    END DO
    DO k=1,self%nholes
      tmp=tmp+vals(self%np_active+k)*sensors%jumpers(j)%hole_facs(k)
    END DO
    jumpout(j+1)=tmp/mu0
  END DO
  !---Setup history file
  IF(oft_env%head_proc)THEN
    jumper_hist%filedesc = 'ThinCurr current jumper history file'
    CALL jumper_hist%setup('jumpers.hist')
    CALL jumper_hist%add_field('time', 'r8', desc="Simulation time [s]")
    DO i=1,sensors%njumpers
      CALL jumper_hist%add_field(sensors%jumpers(i)%name, 'r8')
    END DO
    CALL jumper_hist%write_header
    CALL jumper_hist%open
    jumpout(1)=t
    CALL jumper_hist%write(data_r8=jumpout)
  END IF
END IF
!---Advance system in time
CALL up%add(0.d0,1.d0,u)
DO i=1,nsteps
  !---Update driven coil dI/dt waveforms
  IF(use_cn)THEN
    CALL bmat%apply(u,g)
    IF(ntimes_curr>0)THEN
      DO j=1,self%n_icoils
        ! Start of step
        icoil_dcurr(j)=linterp(curr_waveform(:,1),curr_waveform(:,j+1),ntimes_curr,t+dt/4.d0,1)
        icoil_dcurr(j)=icoil_dcurr(j)-linterp(curr_waveform(:,1),curr_waveform(:,j+1),ntimes_curr,t-dt/4.d0,1)
        ! End of step
        icoil_dcurr(j)=icoil_dcurr(j)+linterp(curr_waveform(:,1),curr_waveform(:,j+1),ntimes_curr,t+dt*5.d0/4.d0,1)
        icoil_dcurr(j)=icoil_dcurr(j)-linterp(curr_waveform(:,1),curr_waveform(:,j+1),ntimes_curr,t+dt*3.d0/4.d0,1)
      END DO
    END IF
    IF(ntimes_volt>0)THEN
      IF(volt_full)THEN
        CALL linterp_facs(volt_waveform(:,1),ntimes_volt,t,int_inds,int_facs,1)
        pcoil_volt = (volt_waveform(int_inds(1),2:)*int_facs(1) &
          + volt_waveform(int_inds(2),2:)*int_facs(2))/2.d0
        CALL linterp_facs(volt_waveform(:,1),ntimes_volt,t+dt,int_inds,int_facs,1)
        pcoil_volt = pcoil_volt + (volt_waveform(int_inds(1),2:)*int_facs(1) &
          + volt_waveform(int_inds(2),2:)*int_facs(2))/2.d0
      ELSE
        DO j=1,self%n_vcoils
          ! Start of step
          pcoil_volt(j)=linterp(volt_waveform(:,1),volt_waveform(:,j+1),ntimes_volt,t,1)/2.d0
          pcoil_volt(j)=pcoil_volt(j)+linterp(volt_waveform(:,1),volt_waveform(:,j+1),ntimes_volt,t+dt,1)/2.d0
        END DO
      END IF
      pcoil_volt=pcoil_volt*dt
    END IF
  ELSE
    CALL Lmat%apply(u,g)
    IF(ntimes_curr>0)THEN
      DO j=1,self%n_icoils
        icoil_dcurr(j)=linterp(curr_waveform(:,1),curr_waveform(:,j+1),ntimes_curr,t+dt*5.d0/4.d0,1)
        icoil_dcurr(j)=icoil_dcurr(j)-linterp(curr_waveform(:,1),curr_waveform(:,j+1),ntimes_curr,t+dt*3.d0/4.d0,1)
      END DO
      icoil_dcurr=icoil_dcurr*2.d0
    END IF
    IF(ntimes_volt>0)THEN
      IF(volt_full)THEN
        CALL linterp_facs(volt_waveform(:,1),ntimes_volt,t+dt,int_inds,int_facs,1)
        pcoil_volt = volt_waveform(int_inds(1),2:)*int_facs(1) &
          + volt_waveform(int_inds(2),2:)*int_facs(2)
      ELSE
        DO j=1,self%n_vcoils
          pcoil_volt(j)=linterp(volt_waveform(:,1),volt_waveform(:,j+1),ntimes_volt,t+dt,1)
        END DO
      END IF
      pcoil_volt=pcoil_volt*dt
    END IF
  END IF
  uu=g%dot(g)
  CALL g%get_local(vals)
  IF(ntimes_curr>0)CALL dgemv('N',self%nelems,self%n_icoils,-1.d0,self%Ael2dr, &
    self%nelems,icoil_dcurr,1,1.d0,vals,1)
  IF(volt_full)THEN
    vals=vals+pcoil_volt
  ELSE
    DO j=1,self%n_vcoils
      vals(self%np_active+self%nholes+j)=vals(self%np_active+self%nholes+j)+pcoil_volt(j)
    END DO
  END IF
  CALL g%restore_local(vals)
  IF(direct)THEN
    CALL Minv%apply(g,u)
    nits=1
  ELSE
    CALL up%add(-1.d0,1.d0,u)
    CALL u%add(1.d0,1.d0,up)
    CALL up%add(0.d0,1.d0,u)
    CALL linv%apply(u,g)
    nits=linv%cits
  END IF
  uu=SQRT(u%dot(u))
  t=t+dt
  IF(MOD(i,nstatus)==0)WRITE(*,*)'Timestep ',i,REAL(t,4),REAL(uu,4),nits
  IF(MOD(i,nplot)==0)THEN
    WRITE(pltnum,'(I4.4)')i
    CALL tw_rst_save(self,u,'pThinCurr_'//pltnum//'.rst','U')
    CALL hdf5_write(t,'pThinCurr_'//pltnum//'.rst','time')
    CALL hdf5_write(icoil_curr,'pThinCurr_'//pltnum//'.rst','coil_currents')
  END IF
  IF(ntimes_curr>0)THEN
    DO j=1,self%n_icoils
      icoil_curr(j)=linterp(curr_waveform(:,1),curr_waveform(:,j+1),ntimes_curr,t,1)
    END DO
  END IF
  CALL u%get_local(vals)
  IF(sensors%nfloops>0)THEN
    CALL dgemv('N',sensors%nfloops,self%nelems,1.d0,self%Ael2sen,sensors%nfloops, &
      vals,1,0.d0,senout(2),1)
    IF(ntimes_curr>0)CALL dgemv('N',sensors%nfloops,self%n_icoils,1.d0,self%Adr2sen, &
      sensors%nfloops,icoil_curr,1,1.d0,senout(2),1)
    IF((ntimes_volt>0).AND.ASSOCIATED(sensor_vals))THEN
      CALL linterp_facs(sensor_vals(:,1),ntimes_volt,t,int_inds,int_facs,1)
      DO j=1,sensors%nfloops
        senout(j+1) = senout(j+1) + sensor_vals(int_inds(1),j+1)*int_facs(1) &
          + sensor_vals(int_inds(2),j+1)*int_facs(2)
      END DO
    END IF
    senout(1)=t
    CALL floop_hist%write(data_r8=senout)
  END IF
  IF(sensors%njumpers>0)THEN
    DO j=1,sensors%njumpers
      tmp=0.d0
      val_prev=0.d0
      ind1=self%pmap(sensors%jumpers(j)%points(1))
      IF(ind1>0)val_prev=vals(ind1)
      DO k=1,sensors%jumpers(j)%np-1
        ind1=self%pmap(sensors%jumpers(j)%points(k+1))
        IF(ind1>0)THEN
          tmp=tmp+vals(ind1)-val_prev
          val_prev=vals(ind1)
        ELSE
          tmp=tmp-val_prev
          val_prev=0.d0
        END IF
      END DO
      DO k=1,self%nholes
        tmp=tmp+vals(self%np_active+k)*sensors%jumpers(j)%hole_facs(k)
      END DO
      jumpout(j+1)=tmp/mu0
    END DO
    jumpout(1)=t
    CALL jumper_hist%write(data_r8=jumpout)
  END IF
END DO
!---Copy solution to input
CALL u%get_local(vals)
vec=vals
!---Cleanup
IF(sensors%nfloops>0)THEN
  CALL floop_hist%close
  DEALLOCATE(senout)
END IF
IF(sensors%njumpers>0)THEN
  CALL jumper_hist%close
  DEALLOCATE(jumpout)
END IF
CALL u%delete()
CALL up%delete()
CALL g%delete()
DEALLOCATE(vals,icoil_curr,icoil_dcurr,pcoil_volt,u,up,g)
IF(direct)THEN
  DEALLOCATE(minv%m)
ELSE
  CALL linv%pre%delete()
  IF(.NOT.PRESENT(hodlr_op))DEALLOCATE(linv%pre)
  CALL linv%delete()
  DEALLOCATE(linv)
END IF
END SUBROUTINE run_td_sim
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE plot_td_sim(self,nsteps,nplot,sensors,compute_B,rebuild_sensors,sensor_vals,hodlr_op)
TYPE(tw_type), INTENT(inout) :: self !< Needs Docs
INTEGER(4), INTENT(in) :: nsteps !< Needs Docs
INTEGER(4), INTENT(in) :: nplot !< Needs Docs
TYPE(tw_sensors), INTENT(in) :: sensors !< Needs Docs
LOGICAL, INTENT(in) :: compute_B !< Needs Docs
LOGICAL, INTENT(in) :: rebuild_sensors !< Needs Docs
REAL(8), POINTER, INTENT(in) :: sensor_vals(:,:)
TYPE(oft_tw_hodlr_op), TARGET, OPTIONAL, INTENT(inout) :: hodlr_op !< HODLR L matrix
!---
INTEGER(4) :: i,j,jj,k,ntimes_curr,ncols,itime,io_unit,face,ind1,ind2,int_inds(2)
REAL(8) :: uu,t,tmp,area,tmp2,val_prev,int_facs(2)
REAL(8), ALLOCATABLE, DIMENSION(:) :: coil_vec,senout,jumpout
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: cc_vals
REAL(8), POINTER, DIMENSION(:) :: vals,vtmp
CLASS(oft_vector), POINTER :: u,Bx,By,Bz
TYPE(oft_bin_file) :: floop_hist,jumper_hist
LOGICAL :: exists
CHARACTER(LEN=4) :: pltnum
WRITE(*,*)'Post-processing simulation'
ALLOCATE(coil_vec(MAX(1,self%n_icoils)))
CALL self%Uloc%new(u)
ALLOCATE(vals(self%nelems))
!---
IF(rebuild_sensors)THEN
  IF(sensors%nfloops>0)THEN
    ALLOCATE(senout(sensors%nfloops+1))
    senout = 0.d0
    !---Setup history file
    IF(oft_env%head_proc)THEN
      floop_hist%filedesc = 'ThinCurr flux loop history file'
      CALL floop_hist%setup('floops.hist')
      CALL floop_hist%add_field('time', 'r8', desc="Simulation time [s]")
      DO i=1,sensors%nfloops
        CALL floop_hist%add_field(sensors%floops(i)%name, 'r8')
      END DO
      CALL floop_hist%write_header
      CALL floop_hist%open
    END IF
  END IF
  IF(sensors%njumpers>0)THEN
    ALLOCATE(jumpout(sensors%njumpers+1))
    jumpout = 0.d0
    !---Setup history file
    IF(oft_env%head_proc)THEN
      jumper_hist%filedesc = 'ThinCurr current jumper history file'
      CALL jumper_hist%setup('jumpers.hist')
      CALL jumper_hist%add_field('time', 'r8', desc="Simulation time [s]")
      DO i=1,sensors%njumpers
        CALL jumper_hist%add_field(sensors%jumpers(i)%name, 'r8')
      END DO
      CALL jumper_hist%write_header
      CALL jumper_hist%open
    END IF
  END IF
END IF
!
IF(compute_B)THEN
  IF(.NOT.ALLOCATED(cc_vals))ALLOCATE(cc_vals(3,self%mesh%np))
  IF(PRESENT(hodlr_op))THEN
    IF(.NOT.ASSOCIATED(hodlr_op%aca_B_dense))CALL hodlr_op%compute_B()
    CALL self%Uloc_pts%new(Bx)
    CALL self%Uloc_pts%new(By)
    CALL self%Uloc_pts%new(Bz)
  ELSE
    IF(.NOT.ASSOCIATED(self%Bel))CALL tw_compute_Bops(self)
  END IF
END IF
!
CALL self%xdmf%clear_timesteps()
DO i=0,nsteps
  IF(MOD(i,nplot)/=0)CYCLE
  !
  WRITE(pltnum,'(I4.4)')i
  CALL tw_rst_load(u,'pThinCurr_'//pltnum//'.rst','U')
  CALL hdf5_read(t,'pThinCurr_'//pltnum//'.rst','time')
  IF(self%n_icoils>0)CALL hdf5_read(coil_vec,'pThinCurr_'//pltnum//'.rst','coil_currents')
  !
  CALL self%xdmf%add_timestep(t)
  CALL u%get_local(vals)
  CALL tw_save_pfield(self,vals,'J')
  !
  IF(compute_B)THEN
    IF(PRESENT(hodlr_op))THEN
      CALL hodlr_op%apply_bop(u,Bx,By,Bz)
      NULLIFY(vtmp)
      CALL Bx%get_local(vtmp)
      cc_vals(1,:)=vtmp
      CALL By%get_local(vtmp)
      cc_vals(2,:)=vtmp
      CALL Bz%get_local(vtmp)
      cc_vals(3,:)=vtmp
      ! IF(self%n_icoils>0)THEN
      !   !$omp parallel do private(k,tmp) collapse(2)
      !     DO j=1,self%n_icoils
      !     DO jj=1,3
      !         !$omp simd
      !         DO k=1,self%mesh%np
      !             cc_vals(jj,k)=cc_vals(jj,k)+coil_vec(j)*hodlr_op%Icoil_Bmat(k,j,jj)
      !         END DO
      !     END DO
      !     END DO
      ! END IF
    ELSE
      !$omp parallel do private(j,jj,tmp)
      DO k=1,self%mesh%np
        DO jj=1,3
          tmp=0.d0
          !$omp simd reduction(+:tmp)
          DO j=1,self%nelems
            tmp=tmp+vals(j)*self%Bel(j,k,jj)
          END DO
          cc_vals(jj,k)=tmp
        END DO
      END DO
    END IF
    IF(self%n_icoils>0)THEN
      !$omp parallel do private(j,tmp) collapse(2)
      DO k=1,self%mesh%np
        DO jj=1,3
          tmp = 0.d0
          !$omp simd reduction(+:tmp)
          DO j=1,self%n_icoils
            tmp=tmp+coil_vec(j)*self%Bdr(k,j,jj)
          END DO
          cc_vals(jj,k)=cc_vals(jj,k)+tmp
        END DO
      END DO
    ! END IF
    END IF
    CALL self%mesh%save_vertex_vector(cc_vals,self%xdmf,'B_v')
  END IF
  IF(rebuild_sensors)THEN
    !
    IF(sensors%nfloops>0)THEN
      CALL dgemv('N',sensors%nfloops,self%nelems,1.d0,self%Ael2sen, &
        sensors%nfloops,vals,1,0.d0,senout(2),1)
      IF(self%n_icoils>0)CALL dgemv('N',sensors%nfloops,self%n_icoils,1.d0,self%Adr2sen,sensors%nfloops, &
        coil_vec,1,1.d0,senout(2),1)
      IF(ASSOCIATED(sensor_vals))THEN
        CALL linterp_facs(sensor_vals(:,1),SIZE(sensor_vals,DIM=1),t,int_inds,int_facs,1)
        DO j=1,sensors%nfloops
          senout(j+1) = senout(j+1) + sensor_vals(int_inds(1),j+1)*int_facs(1) &
            + sensor_vals(int_inds(2),j+1)*int_facs(2)
        END DO
      END IF
      senout(1)=t
      CALL floop_hist%write(data_r8=senout)
    END IF
    IF(sensors%njumpers>0)THEN
      DO j=1,sensors%njumpers
        tmp=0.d0
        val_prev=0.d0
        ind1=self%pmap(sensors%jumpers(j)%points(1))
        IF(ind1>0)val_prev=vals(ind1)
        DO k=1,sensors%jumpers(j)%np-1
          ind1=self%pmap(sensors%jumpers(j)%points(k+1))
          IF(ind1>0)THEN
            tmp=tmp+vals(ind1)-val_prev
            val_prev=vals(ind1)
          ELSE
            tmp=tmp-val_prev
            val_prev=0.d0
          END IF
        END DO
        DO k=1,self%nholes
          tmp=tmp+vals(self%np_active+k)*sensors%jumpers(j)%hole_facs(k)
        END DO
        jumpout(j+1)=tmp/mu0
      END DO
      jumpout(1)=t
      CALL jumper_hist%write(data_r8=jumpout)
    END IF
  END IF
END DO
!---Cleanup
IF(rebuild_sensors)THEN
IF(sensors%nfloops>0)THEN
  CALL floop_hist%close()
  DEALLOCATE(senout)
END IF
IF(sensors%njumpers>0)THEN
  CALL jumper_hist%close()
  DEALLOCATE(jumpout)
END IF
END IF
CALL u%delete
DEALLOCATE(u)
DEALLOCATE(vals)
IF(compute_B)THEN
  CALL Bx%delete()
  CALL By%delete()
  CALL Bz%delete()
  DEALLOCATE(vtmp)
END IF
IF(self%n_icoils>0)DEALLOCATE(coil_vec)
END SUBROUTINE plot_td_sim
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE tw_reduce_model(self,sensors,neigs,eig_vec,filename,compute_B,hodlr_op)
TYPE(tw_type), INTENT(inout) :: self !< Needs docs
TYPE(tw_sensors), INTENT(in) :: sensors !< Sensor information
INTEGER(4), INTENT(in) :: neigs !< Needs docs
REAL(8), INTENT(in) :: eig_vec(:,:) !< Needs docs
CHARACTER(LEN=*), INTENT(in) :: filename !< Needs docs
LOGICAL, INTENT(in) :: compute_B !< Needs docs
TYPE(oft_tw_hodlr_op), TARGET, OPTIONAL, INTENT(inout) :: hodlr_op !< HODLR L matrix
!---
INTEGER(4) :: i,j
REAL(8), POINTER :: vals(:),pvals(:)
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: Mattmp,Mat_red,Bxtmp,Bytmp,Bztmp
CHARACTER(LEN=4) :: sub_coil_id
CLASS(oft_vector), POINTER :: atmp,btmp,Bx,By,Bz
character(LEN=80), dimension(1) :: description
CALL hdf5_create_file(TRIM(filename))
CALL hdf5_write(tw_idx_ver,TRIM(filename),tw_idx_path)
!---Setup temporaries
NULLIFY(vals)
CALL self%Uloc%new(atmp)
CALL self%Uloc%new(btmp)
CALL atmp%get_local(vals)
ALLOCATE(Mat_red(neigs,neigs),Mattmp(self%nelems,neigs))
!---Save reduction vectors
CALL hdf5_write(eig_vec,TRIM(filename),'Basis')
description=["Basis set used for model reduction"]
CALL hdf5_add_string_attribute(TRIM(filename),'Basis','description',description)
!---Reduce L matrix
IF(PRESENT(hodlr_op))THEN
  DO i=1,neigs
    vals=eig_vec(:,i)
    CALL atmp%restore_local(vals)
    CALL btmp%set(0.d0)
    CALL hodlr_op%apply(atmp,btmp)
    CALL btmp%get_local(vals)
    Mattmp(:,i) = vals
  END DO
ELSE
  CALL dgemm('N','N',self%nelems,neigs,self%nelems,1.d0, &
    self%Lmat,self%nelems,eig_vec,self%nelems,0.d0,Mattmp,self%nelems)
END IF
CALL dgemm('T','N',neigs,neigs,self%nelems,1.d0,eig_vec, &
  self%nelems,Mattmp,self%nelems,0.d0,Mat_red,neigs)
CALL hdf5_write(Mat_red,TRIM(filename),'L')
description=["Self-inductance matrix for reduced model"]
CALL hdf5_add_string_attribute(TRIM(filename),'L','description',description)
!---Reduce R matrix
DO i=1,neigs
  vals=eig_vec(:,i)
  CALL atmp%restore_local(vals)
  CALL btmp%set(0.d0)
  CALL self%Rmat%apply(atmp,btmp)
  CALL btmp%get_local(vals)
  Mattmp(:,i) = vals
END DO
CALL dgemm('T','N',neigs,neigs,self%nelems,1.d0, &
  eig_vec,self%nelems,Mattmp,self%nelems,0.d0,Mat_red,neigs)
CALL hdf5_write(Mat_red,TRIM(filename),'R')
description=["Resistance matrix for reduced model"]
CALL hdf5_add_string_attribute(TRIM(filename),'R','description',description)
CALL btmp%delete()
DEALLOCATE(Mat_red,Mattmp,btmp)
!---Reduce sensor coupling matrix
IF(sensors%nfloops>0)THEN
  ! Save sensor names and indexes
  CALL hdf5_create_group(TRIM(filename),'SENSORS')
  DO i=1,sensors%nfloops
    CALL hdf5_create_group(TRIM(filename), &
      'SENSORS/'//TRIM(sensors%floops(i)%name))
    CALL hdf5_write(i,TRIM(filename), &
      'SENSORS/'//TRIM(sensors%floops(i)%name)//'/index')
    CALL hdf5_write(sensors%floops(i)%scale_fac,TRIM(filename), &
      'SENSORS/'//TRIM(sensors%floops(i)%name)//'/scale')
    CALL hdf5_write(sensors%floops(i)%r,TRIM(filename), &
      'SENSORS/'//TRIM(sensors%floops(i)%name)//'/pts')
  END DO
  ! Save mutual coupling matrix
  ALLOCATE(Mat_red(sensors%nfloops,neigs))
  CALL dgemm('N','N',sensors%nfloops,neigs,self%nelems,1.d0, &
    self%Ael2sen,sensors%nfloops,eig_vec,self%nelems,0.d0,Mat_red,sensors%nfloops)
  CALL hdf5_write(Mat_red,TRIM(filename),'Ms')
  description=["Model to sensor mutual inductance matrix"]
  CALL hdf5_add_string_attribute(TRIM(filename),'Ms','description',description)
  DEALLOCATE(Mat_red)
END IF
!---Reduce coil coupling matrix
IF(self%n_icoils>0)THEN
  ! Save coil names and indexes
  CALL hdf5_create_group(TRIM(filename),'COILS')
  DO i=1,self%n_icoils
    CALL hdf5_create_group(TRIM(filename), &
      'COILS/'//TRIM(self%icoils(i)%name))
    CALL hdf5_write(i,TRIM(filename), &
      'COILS/'//TRIM(self%icoils(i)%name)//'/index')
    CALL hdf5_create_group(TRIM(filename), &
      'COILS/'//TRIM(self%icoils(i)%name)//'/SUB_COILS')
    DO j=1,self%icoils(i)%ncoils
      WRITE(sub_coil_id,'(I4.4)')j
      CALL hdf5_create_group(TRIM(filename), &
        'COILS/'//TRIM(self%icoils(i)%name)//'/SUB_COILS/'//sub_coil_id)
      CALL hdf5_write(self%icoils(i)%scales(j),TRIM(filename), &
        'COILS/'//TRIM(self%icoils(i)%name)//'/SUB_COILS/'//sub_coil_id//'/scale')
      CALL hdf5_write(self%icoils(i)%coils(j)%pts,TRIM(filename), &
        'COILS/'//TRIM(self%icoils(i)%name)//'/SUB_COILS/'//sub_coil_id//'/pts')
    END DO
  END DO
  ! Save mutual coupling matrix
  ALLOCATE(Mat_red(neigs,self%n_icoils))
  CALL dgemm('T','N',neigs,self%n_icoils,self%nelems,1.d0, &
    eig_vec,self%nelems,self%Ael2dr,self%nelems,0.d0,Mat_red,neigs)
  CALL hdf5_write(Mat_red,TRIM(filename),'Mc')
  description=["Model to coil mutual inductance matrix"]
  CALL hdf5_add_string_attribute(TRIM(filename),'Mc','description',description)
  DEALLOCATE(Mat_red)
  IF(sensors%nfloops>0)THEN
    CALL hdf5_write(self%Adr2sen,TRIM(filename),'Msc')
    description=["Coil to sensor mutual inductance matrix"]
    CALL hdf5_add_string_attribute(TRIM(filename),'Msc','description',description)
  END IF
END IF
!---Reduce B matrix
IF(compute_B)THEN
  CALL self%Uloc_pts%new(Bx)
  CALL self%Uloc_pts%new(By)
  CALL self%Uloc_pts%new(Bz)
  ALLOCATE(Bxtmp(self%mesh%np,neigs),Bytmp(self%mesh%np,neigs),Bztmp(self%mesh%np,neigs))
  IF(PRESENT(hodlr_op))THEN
    IF(.NOT.ASSOCIATED(hodlr_op%aca_B_dense))CALL hodlr_op%compute_B()
    NULLIFY(pvals)
    DO i=1,neigs
      vals=eig_vec(:,i)
      CALL atmp%restore_local(vals)
      CALL Bx%set(0.d0)
      CALL hodlr_op%apply_bop(atmp,Bx,By,Bz)
      CALL Bx%get_local(pvals)
      Bxtmp(:,i) = pvals
      CALL By%get_local(pvals)
      Bytmp(:,i) = pvals
      CALL Bz%get_local(pvals)
      Bztmp(:,i) = pvals
    END DO
    DEALLOCATE(pvals)
  ELSE
    IF(.NOT.ASSOCIATED(self%Bel))CALL tw_compute_Bops(self)
    CALL dgemm('N','N',self%mesh%np,neigs,self%nelems,1.d0, &
      self%Bel(:,:,1),self%mesh%np,eig_vec,self%nelems,0.d0,Bxtmp,self%mesh%np)
    CALL dgemm('N','N',self%mesh%np,neigs,self%nelems,1.d0, &
      self%Bel(:,:,2),self%mesh%np,eig_vec,self%nelems,0.d0,Bytmp,self%mesh%np)
    CALL dgemm('N','N',self%mesh%np,neigs,self%nelems,1.d0, &
      self%Bel(:,:,3),self%mesh%np,eig_vec,self%nelems,0.d0,Bztmp,self%mesh%np)
  END IF
  CALL hdf5_write(Bxtmp,TRIM(filename),'Bx')
  description=["X-component of B reconstruction matrix for model"]
  CALL hdf5_add_string_attribute(TRIM(filename),'Bx','description',description)
  CALL hdf5_write(Bytmp,TRIM(filename),'By')
  description=["Y-component of B reconstruction matrix for model"]
  CALL hdf5_add_string_attribute(TRIM(filename),'By','description',description)
  CALL hdf5_write(Bztmp,TRIM(filename),'Bz')
  description=["Z-component of B reconstruction matrix for model"]
  CALL hdf5_add_string_attribute(TRIM(filename),'Bz','description',description)
  DEALLOCATE(Bxtmp,Bytmp,Bztmp)
  CALL Bx%delete()
  CALL By%delete()
  CALL Bz%delete()
  DEALLOCATE(Bx,By,Bz)
  IF(self%n_icoils>0)THEN
    CALL hdf5_write(self%Bdr(:,:,1),TRIM(filename),'Bx_c')
    description=["X-component of B reconstruction matrix for coils"]
    CALL hdf5_add_string_attribute(TRIM(filename),'Bx_c','description',description)
    CALL hdf5_write(self%Bdr(:,:,2),TRIM(filename),'By_c')
    description=["Y-component of B reconstruction matrix for coils"]
    CALL hdf5_add_string_attribute(TRIM(filename),'By_c','description',description)
    CALL hdf5_write(self%Bdr(:,:,3),TRIM(filename),'Bz_c')
    description=["Z-component of B reconstruction matrix for coils"]
    CALL hdf5_add_string_attribute(TRIM(filename),'Bz_c','description',description)
  END IF
END IF
CALL atmp%delete()
DEALLOCATE(vals,atmp)
END SUBROUTINE tw_reduce_model
END MODULE thin_wall_solvers