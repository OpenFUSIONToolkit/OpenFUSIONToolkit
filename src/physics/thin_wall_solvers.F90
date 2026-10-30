!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file thin_wall_solvers.F90
!
!> Module for thin-wall modeling on 3D triangular meshes
!!
!!
!! @authors Chris Hansen
!! @date May 2017
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
MODULE thin_wall_solvers
USE oft_base
USE oft_sort, ONLY: sort_array
USE oft_io, ONLY: oft_bin_file, hdf5_add_string_attribute
!
USE oft_la_base, ONLY: oft_vector, oft_vector_ptr, oft_cvector, oft_matrix, oft_graph, &
  vector_extrapolate
USE oft_lu, ONLY: oft_lusolver, lapack_matinv, lapack_cholesky
USE oft_native_la, ONLY: oft_native_dense_matrix, oft_native_dense_cmatrix
USE oft_deriv_matrices, ONLY: oft_sum_matrix, oft_sum_cmatrix
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_native_gmres_csolver, oft_native_dense_cbjprecond, oft_diag_cscale
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
!---------------------------------------------------------------------------------
!> Compute L/R eigenmodes of ThinCurr model using a direct approach via LAPACK
!---------------------------------------------------------------------------------
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
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Starting eigenvalue solve (direct)'
CALL oft_increase_indent
CALL loctimer%tick
ALLOCATE(WR(N),WI(N),VL(LDVL,N),VR(LDVR,N),WORK(1))
LWORK = -1
CALL DGEEV(JOBVL, JOBVR, N, Acomp, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
LWORK = INT(WORK(1),4)
IF(oft_debug_print(1))WRITE(*,'(2A,I8)')oft_indent,'LU Block size = ',lwork/N
DEALLOCATE(work)
ALLOCATE(work(lwork))
CALL DGEEV(JOBVL, JOBVR, N, Acomp, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
DEALLOCATE(VL,WORK,Acomp)
elapsed_time=loctimer%tock()
WRITE(*,'(2A,ES12.4,I4)')oft_indent,'Time = ',elapsed_time,INFO
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
WRITE(*,'(2A)')oft_indent,'Eigenvalues'
CALL oft_increase_indent
DO i=1,MIN(neigs,5)
  WRITE(*,'(A,ES14.6)')oft_indent,eig_rval(i)
END DO
CALL oft_decrease_indent
CALL oft_decrease_indent
DEBUG_STACK_POP
END SUBROUTINE lr_eigenmodes_direct
!---------------------------------------------------------------------------------
!> Compute L/R eigenmodes of ThinCurr model using an iterative Lanczos method via ARPACK
!---------------------------------------------------------------------------------
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
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Starting eigenvalue solve (ARPACK)'
CALL oft_increase_indent
!---Setup local vectors
CALL self%Uloc%new(uloc)
CALL self%Rmat%assemble(uloc)

!---Setup R^-1 solver
NULLIFY(linv)
CALL create_native_solver(linv,'lu')
SELECT TYPE(this=>linv)
CLASS IS(oft_lusolver)
  IF(this%package(1:4)=='none')THEN
    WRITE(*,'(2A)')oft_indent,'Note: LU solver not available, falling back to CG'
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
WRITE(*,'(2A,ES12.4)')oft_indent,'Time = ',elapsed_time
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
  WRITE(*,'(2A)')oft_indent,'Eigenvalues'
  CALL oft_increase_indent
  DO i=1,MIN(neigs,5)
    WRITE(*,'(A,ES14.6)')oft_indent,eig_rval(i)
  END DO
  CALL oft_decrease_indent
END IF
!---Cleanup
CALL uloc%delete()
CALL arsolver%delete()
DEALLOCATE(uloc)
CALL oft_decrease_indent
#else
CALL oft_abort("Iterative eigenvalue solve requires ARPACK", "lr_eigenmodes_arpack", __FILE__)
#endif
DEBUG_STACK_POP
END SUBROUTINE lr_eigenmodes_arpack
!---------------------------------------------------------------------------------
!> Compute frequency response of a ThinCurr model to given driver voltages
!!
!! Solves \f$ (i \omega L + R) I = V \f$ for the complex currents \f$ I \f$, where
!! \f$ \omega = 2 \pi f \f$, or the inductive (\f$ i L I = V \f$) or resistive
!! (\f$ R I = V \f$) limit if `fr_limit/=0`. The system is solved either directly
!! using LAPACK (dense \f$ L \f$ only) or iteratively using GMRES. The iterative
!! solver is preconditioned using block-Jacobi on the HODLR diagonal blocks (if
!! `hodlr_op` is present) or on blocks from partitioning the graph of \f$ R \f$
!! (diagonal scaling for models with <=2000 elements).
!---------------------------------------------------------------------------------
SUBROUTINE frequency_response(self,direct,fr_limit,freq,driver,hodlr_op)
TYPE(tw_type), INTENT(in) :: self !< ThinCurr object
LOGICAL, INTENT(in) :: direct !< Use direct solver? (not supported with `hodlr_op`)
REAL(r8), INTENT(in) :: freq !< Frequency for response calculation [Hz] (unused if `fr_limit/=0`)
INTEGER(i4), INTENT(in) :: fr_limit !< Frequency limit for calculation (0: none, 1: inductive, 2: resistive)
REAL(r8), INTENT(inout) :: driver(:,:) !< Driver voltages (real, imaginary), overwritten with resulting currents [self%nelems,2]
TYPE(oft_tw_hodlr_op), TARGET, OPTIONAL, INTENT(inout) :: hodlr_op !< HODLR L matrix (optional)
INTEGER(i4) :: i,j,info
DOUBLE COMPLEX, POINTER, DIMENSION(:) :: x,b
TYPE(oft_sum_cmatrix), TARGET :: aca_Fmat
TYPE(oft_native_dense_cmatrix), TARGET :: dense_Fmat
CLASS(oft_cvector), POINTER :: aloc,bloc
TYPE(oft_graph) :: graph
TYPE(oft_native_gmres_csolver), TARGET :: frinv
TYPE(oft_tw_hodlr_bjpre), TARGET :: frinv_pre
TYPE(oft_native_dense_cbjprecond), TARGET :: frinv_dense_pre
TYPE(oft_diag_cscale), TARGET :: frinv_diag_pre
!---
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Starting frequency-domain run'
CALL oft_increase_indent
!---Setup matrix
ALLOCATE(b(self%nelems))
b=(1.d0,0.d0)*driver(:,1) + (0.d0,1.d0)*driver(:,2)
IF(PRESENT(hodlr_op))THEN
  SELECT CASE(fr_limit)
  CASE(0)
    WRITE(*,'(2A,ES13.5)')oft_indent,'Frequency [Hz] = ',freq
    aca_Fmat%rJ=>hodlr_op
    aca_Fmat%rK=>self%Rmat
    aca_Fmat%beta=(0.d0,1.d0)*freq*2.d0*pi
    aca_Fmat%alam=(1.d0,0.d0)
  CASE(1)
    WRITE(*,'(2A)')oft_indent,'Frequency -> Inf (L limit)'
    aca_Fmat%rJ=>hodlr_op
    aca_Fmat%rK=>self%Rmat
    aca_Fmat%beta=(0.d0,1.d0)
    aca_Fmat%alam=(0.d0,0.d0)
  CASE(2)
    WRITE(*,'(2A)')oft_indent,'Frequency -> 0   (R limit)'
    aca_Fmat%rJ=>self%Rmat
    aca_Fmat%rK=>self%Rmat
    aca_Fmat%beta=(0.d0,0.d0)
    aca_Fmat%alam=(1.d0,0.d0)
  CASE DEFAULT
    CALL oft_abort('Invalid "fr_limit" value (0,1,2)', 'run_fr', __FILE__)
  END SELECT
ELSE
  dense_Fmat%nr=self%nelems; dense_Fmat%nc=self%nelems
  dense_Fmat%nrg=self%nelems; dense_Fmat%ncg=self%nelems
  ALLOCATE(dense_Fmat%M(self%nelems,self%nelems))
  SELECT CASE(fr_limit)
  CASE(0)
    WRITE(*,'(2A,ES13.5)')oft_indent,'Frequency [Hz] = ',freq
    dense_Fmat%M=(0.d0,1.d0)*freq*2.d0*pi*self%Lmat
    DO i=1,self%Rmat%nr
      DO j=self%Rmat%kr(i),self%Rmat%kr(i+1)-1
        dense_Fmat%M(i,self%Rmat%lc(j))=dense_Fmat%M(i,self%Rmat%lc(j)) + self%Rmat%M(j)
      END DO
    END DO
  CASE(1)
    WRITE(*,'(2A)')oft_indent,'Frequency -> Inf (L limit)'
    dense_Fmat%M=(0.d0,1.d0)*self%Lmat
  CASE(2)
    WRITE(*,'(2A)')oft_indent,'Frequency -> 0   (R limit)'
    dense_Fmat%M=(0.d0,0.d0)
    DO i=1,self%Rmat%nr
      DO j=self%Rmat%kr(i),self%Rmat%kr(i+1)-1
        dense_Fmat%M(i,self%Rmat%lc(j))=dense_Fmat%M(i,self%Rmat%lc(j)) + self%Rmat%M(j)
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
  CALL lapack_matinv(self%nelems,dense_Fmat%M,info)
  IF(info/=0)CALL oft_abort('Matrix inversion failed','frequency_response',__FILE__)
  CALL zgemv('N',self%nelems,self%nelems,(1.d0,0.d0),dense_Fmat%M,self%nelems,b,1,(0.d0,0.d0),x,1)
  CALL dense_Fmat%delete()
ELSE
  NULLIFY(aloc,bloc)
  CALL self%Uloc%new(aloc)
  CALL self%Uloc%new(bloc)
  frinv%nrits=60
  frinv%its=-1
  IF(PRESENT(hodlr_op))THEN
    CALL aca_Fmat%assemble(aloc)
    frinv%pre=>frinv_pre
    frinv_pre%mf_obj=>hodlr_op
    frinv_pre%Rmat=>self%Rmat
    IF(ASSOCIATED(self%junc))frinv_pre%Rmat=>self%junc%Rmat_exp
    frinv_pre%alpha=aca_Fmat%beta
    frinv_pre%beta=aca_Fmat%alam
    frinv%A=>aca_Fmat
  ELSE
    frinv%A=>dense_Fmat
    !---Converge on |r|/|u| only (no absolute or relative tolerance)
    frinv%min_nits=10
    frinv%atol=0.d0
    frinv%rtol=0.d0
    !---Block-Jacobi preconditioner partitioned using the graph of R
    !   (about 2000 owned rows per block, more with overlap)
    frinv_dense_pre%nlocal=CEILING(REAL(self%nelems,r8)/2000.d0)
    IF(frinv_dense_pre%nlocal>1)THEN
      graph%nr=self%Rmat%nr
      graph%nnz=self%Rmat%nnz
      graph%kr=>self%Rmat%kr
      graph%lc=>self%Rmat%lc
      frinv_dense_pre%loverlap=.TRUE.
      CALL frinv_dense_pre%setup(dense_Fmat,graph)
      frinv%pre=>frinv_dense_pre
    ELSE
      CALL dense_Fmat%assemble(aloc)
      frinv%pre=>frinv_diag_pre
    END IF
  END IF
  !---Solve system
  CALL aloc%restore_local(x)
  CALL bloc%restore_local(b)
  CALL frinv%apply(aloc,bloc)
  CALL aloc%get_local(x)
  !---Cleanup temporaries
  CALL aloc%delete()
  CALL bloc%delete()
  DEALLOCATE(aloc,bloc)
  CALL frinv%delete()
  IF(PRESENT(hodlr_op))THEN
    CALL frinv_pre%delete()
  ELSE
    CALL frinv_dense_pre%delete()
    CALL frinv_diag_pre%delete()
    CALL dense_Fmat%delete()
  END IF
END IF
driver(:,1)=REAL(x,8)
driver(:,2)=AIMAG(x)
DEALLOCATE(b,x)
CALL oft_decrease_indent
END SUBROUTINE frequency_response
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE run_td_sim(self,times,nsteps,vec,direct,lin_tols,use_cn,nstatus,nplot,sensors,curr_waveform,volt_waveform,sensor_vals,hodlr_op)
TYPE(tw_type), INTENT(in) :: self
INTEGER(4), INTENT(in) :: nsteps
REAL(8), INTENT(in) :: times(nsteps+1)
REAL(8), INTENT(inout) :: vec(:)
LOGICAL, INTENT(in) :: direct
REAL(8), INTENT(in) :: lin_tols(2)
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
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: its_hist
REAL(8) :: uu,t,tmp,area,p2,p1,val_prev,dt,dt_old,dt_op,int_facs(2),elapsed_time
REAL(8), ALLOCATABLE, DIMENSION(:) :: icoil_curr,icoil_dcurr,pcoil_volt,senout,jumpout,eta_check,time_hist
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: cc_vals
REAL(8), POINTER, DIMENSION(:) :: vals
CLASS(oft_vector), POINTER :: u,g,up
CLASS(oft_matrix), POINTER :: Lmat
TYPE(oft_native_dense_matrix), TARGET :: Lmat_dense,Minv
TYPE(oft_sum_matrix), TARGET :: fmat,bmat
CLASS(oft_solver), POINTER :: linv
TYPE(oft_tw_hodlr_rbjpre), TARGET :: linv_pre
TYPE(oft_bin_file) :: floop_hist,jumper_hist
TYPE(oft_timer) :: solve_timer
LOGICAL :: exists,volt_full,pm_save
CHARACTER(LEN=4) :: pltnum
CHARACTER(LEN=15) :: fmt_str
CHARACTER(LEN=OFT_SLEN) :: hole_jumper_name,rst_file
!---Extrapolation fields
integer(i4) :: nextrap,maxextrap
real(r8), allocatable, dimension(:) :: extrapt
type(oft_vector_ptr), allocatable, dimension(:) :: extrap_fields
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Starting time-domain simulation'
CALL oft_increase_indent
WRITE(*,'(2A)')oft_indent,'timestep    time           sol_norm     nits    solver time'
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
!---Get first timestep size
dt=times(2)-times(1)
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
  linv%atol=lin_tols(1)
  linv%rtol=lin_tols(2)
  IF(PRESENT(hodlr_op))THEN
    linv_pre%mf_obj=>hodlr_op
    linv_pre%Rmat=>self%Rmat
    IF(ASSOCIATED(self%junc))linv_pre%Rmat=>self%junc%Rmat_exp
    linv_pre%alpha=1.d0
    linv_pre%beta=fmat%alam
    linv%pre=>linv_pre
  ELSE
    CALL create_diag_pre(linv%pre)
  END IF
  !---Setup extrapolation fields
  maxextrap=3 ! Number of points for extrapolation (0: No extrapolation, 2: Linear, 3: Quadratic)
  IF(maxextrap>0)THEN
    ALLOCATE(extrap_fields(maxextrap),extrapt(maxextrap))
    DO i=1,maxextrap
      CALL self%Uloc%new(extrap_fields(i)%f)
      extrapt(i)=0.d0
    END DO
    nextrap=0
  END IF
ELSE
  maxextrap=-1
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
  WRITE(*,'(2A)')oft_indent,'Starting factorization'
  pm_save=oft_env%pm; oft_env%pm=.TRUE.
  ! CALL lapack_matinv(Minv%nr,Minv%M,info)
  CALL lapack_cholesky(Minv%nr,Minv%M,info)
  oft_env%pm=pm_save
END IF
!---
ALLOCATE(vals(self%nelems))
vals=vec
CALL u%restore_local(vals)
t=times(1)
IF(ntimes_curr>0)THEN
  DO j=1,self%n_icoils
    icoil_curr(j)=linterp(curr_waveform(:,1),curr_waveform(:,j+1),ntimes_curr,t,1)
  END DO
END IF
!---
WRITE(pltnum,'(I4.4)')0
rst_file=TRIM(self%rst_prefix)//'pThinCurr_'//pltnum//'.rst'
CALL tw_rst_save(self,u,TRIM(rst_file),'U')
CALL hdf5_write(t,TRIM(rst_file),'time')
CALL hdf5_write(icoil_curr,TRIM(rst_file),'coil_currents')
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
    CALL floop_hist%setup(TRIM(self%rst_prefix)//'floops.hist')
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
IF(sensors%njumpers+self%nholes+self%n_vcoils>0)THEN
  ALLOCATE(jumpout(sensors%njumpers+self%nholes+self%n_vcoils+1))
  DO j=1,sensors%njumpers
    jumpout(j+1)=sensors%jumpers(j)%compute_current(self,vals)
  END DO
  DO j=1,self%nholes+self%n_vcoils
    jumpout(sensors%njumpers+j+1)=tw_hole_value(self,vals,j)/mu0
  END DO
  !---Setup history file
  IF(oft_env%head_proc)THEN
    jumper_hist%filedesc = 'ThinCurr current jumper history file'
    CALL jumper_hist%setup(TRIM(self%rst_prefix)//'jumpers.hist')
    CALL jumper_hist%add_field('time', 'r8', desc="Simulation time [s]")
    DO i=1,sensors%njumpers
      CALL jumper_hist%add_field(sensors%jumpers(i)%name, 'r8')
    END DO
    DO i=1,self%nholes
      WRITE(hole_jumper_name,'(A,I4.4)')'HOLE_',i
      CALL jumper_hist%add_field(hole_jumper_name, 'r8')
    END DO
    DO i=1,self%n_vcoils
      CALL jumper_hist%add_field(self%vcoils(i)%name, 'r8')
    END DO
    CALL jumper_hist%write_header
    CALL jumper_hist%open
    jumpout(1)=t
    CALL jumper_hist%write(data_r8=jumpout)
  END IF
END IF
!---Advance system in time
CALL up%add(0.d0,1.d0,u)
ALLOCATE(its_hist(nstatus),time_hist(nstatus))
its_hist=0
time_hist=0.d0
DO i=1,nsteps
  dt_old = dt
  dt = times(i+1)-times(i)
  !---Update matrices if timestep has changed
  IF(ABS((dt_old-dt)/dt_old)>1.d-6)THEN
    WRITE(*,'(2X,A,1ES10.2)')'Updating matrices with new timestep = ',dt
    IF(use_cn)THEN
      bmat%alam = -dt/2.d0  ! Update backward matrix timestep
      dt_op = dt/2.d0
    ELSE
      dt_op = dt
    END IF
    IF(direct)THEN
      !---Re-build forward matrix [L + dt_op*R] (overwritten with inverse)
      Minv%M=self%Lmat
      DO k=1,self%Rmat%nr
        DO j=self%Rmat%kr(k),self%Rmat%kr(k+1)-1
          Minv%M(k,self%Rmat%lc(j))=Minv%M(k,self%Rmat%lc(j)) + dt_op*self%Rmat%M(j)
        END DO
      END DO
      WRITE(*,'(2A)')oft_indent,'Starting factorization'
      pm_save=oft_env%pm; oft_env%pm=.TRUE.
      ! CALL lapack_matinv(Minv%nr,Minv%M,info)
      CALL lapack_cholesky(Minv%nr,Minv%M,info)
      oft_env%pm=pm_save
    ELSE
      fmat%alam = dt_op ! Update forward matrix timestep
      CALL fmat%assemble(g)
      IF(PRESENT(hodlr_op))THEN
        linv_pre%beta=fmat%alam ! Update preconditioner timestep
        linv_pre%refactor=.TRUE.
      END IF
    END IF
  END IF

  !---Update extrapolation fields
  IF(maxextrap>0)THEN
    DO j=maxextrap,2,-1
      CALL extrap_fields(j)%f%add(0.d0,1.d0,extrap_fields(j-1)%f)
      extrapt(j)=extrapt(j-1)
    END DO
    IF(nextrap<maxextrap)nextrap=nextrap+1
    CALL extrap_fields(1)%f%add(0.d0,1.d0,u)
    extrapt(1)=t
  END IF
  !---Update driven coil dI/dt waveforms
  IF(ntimes_curr>0)THEN
    DO j=1,self%n_icoils
      icoil_dcurr(j)=linterp(curr_waveform(:,1),curr_waveform(:,j+1),ntimes_curr,t+dt,1)
      icoil_dcurr(j)=icoil_dcurr(j)-linterp(curr_waveform(:,1),curr_waveform(:,j+1),ntimes_curr,t,1)
    END DO
  END IF
  IF(use_cn)THEN
    CALL bmat%apply(u,g)
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
      vals(self%nelems-self%n_vcoils+j)=vals(self%nelems-self%n_vcoils+j)+pcoil_volt(j)
    END DO
  END IF
  CALL g%restore_local(vals)
  CALL solve_timer%tick()
  IF(direct)THEN
    CALL Minv%apply(g,u)
    nits=1
  ELSE
    !---Extrapolate solution
    IF(maxextrap>0)CALL vector_extrapolate(extrapt,extrap_fields,nextrap,t+dt,u)
    CALL linv%apply(u,g)
    nits=linv%cits
  END IF
  elapsed_time=solve_timer%tock()
  IF(PRESENT(hodlr_op))THEN
    ! WRITE(*,'(A,5ES12.3)')'Timing: ',hodlr_op%times,linv_pre%times
    hodlr_op%times=0.d0
    linv_pre%times=0.d0
  END IF
  uu=SQRT(u%dot(u))
  t=times(i+1)
  its_hist(MOD(i,nstatus)+1)=nits
  time_hist(MOD(i,nstatus)+1)=elapsed_time
  IF(MOD(i,nstatus)==0)THEN
    nits=INT(SUM(its_hist)/REAL(nstatus,8),4)
    elapsed_time=SUM(time_hist)/REAL(nstatus,8)
    WRITE(*,'(2X,I6,ES16.6,ES14.4,2X,I6,F12.2)')i,t,uu,nits,elapsed_time
  END IF
  IF(MOD(i,nplot)==0)THEN
    WRITE(pltnum,'(I4.4)')i
    rst_file=TRIM(self%rst_prefix)//'pThinCurr_'//pltnum//'.rst'
    CALL tw_rst_save(self,u,TRIM(rst_file),'U')
    CALL hdf5_write(t,TRIM(rst_file),'time')
    CALL hdf5_write(icoil_curr,TRIM(rst_file),'coil_currents')
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
  IF(sensors%njumpers+self%nholes+self%n_vcoils>0)THEN
    DO j=1,sensors%njumpers
      jumpout(j+1)=sensors%jumpers(j)%compute_current(self,vals)
    END DO
    DO j=1,self%nholes+self%n_vcoils
      jumpout(sensors%njumpers+j+1)=tw_hole_value(self,vals,j)/mu0
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
IF(sensors%njumpers+self%nholes+self%n_vcoils>0)THEN
  CALL jumper_hist%close
  DEALLOCATE(jumpout)
END IF
CALL u%delete()
CALL up%delete()
CALL g%delete()
DEALLOCATE(vals,icoil_curr,icoil_dcurr,pcoil_volt,u,up,g)
IF(maxextrap>0)THEN
  DO i=1,maxextrap
    CALL extrap_fields(i)%f%delete()
  END DO
  DEALLOCATE(extrap_fields,extrapt)
END IF
IF(direct)THEN
  DEALLOCATE(minv%m)
ELSE
  CALL linv%pre%delete()
  IF(.NOT.PRESENT(hodlr_op))DEALLOCATE(linv%pre)
  CALL linv%delete()
  DEALLOCATE(linv)
END IF
CALL oft_decrease_indent
END SUBROUTINE run_td_sim
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
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
CHARACTER(LEN=OFT_SLEN) :: hole_jumper_name,rst_file
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Post-processing time-domain simulation'
CALL oft_increase_indent
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
      CALL floop_hist%setup(TRIM(self%rst_prefix)//'floops.hist')
      CALL floop_hist%add_field('time', 'r8', desc="Simulation time [s]")
      DO i=1,sensors%nfloops
        CALL floop_hist%add_field(sensors%floops(i)%name, 'r8')
      END DO
      CALL floop_hist%write_header
      CALL floop_hist%open
    END IF
  END IF
  IF(sensors%njumpers+self%nholes+self%n_vcoils>0)THEN
    ALLOCATE(jumpout(sensors%njumpers+self%nholes+self%n_vcoils+1))
    jumpout = 0.d0
    !---Setup history file
    IF(oft_env%head_proc)THEN
      jumper_hist%filedesc = 'ThinCurr current jumper history file'
      CALL jumper_hist%setup(TRIM(self%rst_prefix)//'jumpers.hist')
      CALL jumper_hist%add_field('time', 'r8', desc="Simulation time [s]")
      DO i=1,sensors%njumpers
        CALL jumper_hist%add_field(sensors%jumpers(i)%name, 'r8')
      END DO
      DO i=1,self%nholes
        WRITE(hole_jumper_name,'(A,I4.4)')'HOLE_',i
        CALL jumper_hist%add_field(hole_jumper_name, 'r8')
      END DO
      DO i=1,self%n_vcoils
        CALL jumper_hist%add_field(self%vcoils(i)%name, 'r8')
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
    IF(.NOT.ASSOCIATED(self%Bel))CALL tw_compute_Bops(self,self%B_dx)
  END IF
END IF
!
CALL self%xdmf%clear_timesteps()
DO i=0,nsteps
  IF(MOD(i,nplot)/=0)CYCLE
  !
  WRITE(pltnum,'(I4.4)')i
  rst_file=TRIM(self%rst_prefix)//'pThinCurr_'//pltnum//'.rst'
  CALL tw_rst_load(u,TRIM(rst_file),'U')
  CALL hdf5_read(t,TRIM(rst_file),'time')
  IF(self%n_icoils>0)CALL hdf5_read(coil_vec,TRIM(rst_file),'coil_currents')
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
      DEALLOCATE(vtmp)
    ELSE
      DO jj=1,3
        CALL dgemv('N',self%mesh%np,self%nelems,1.d0,self%Bel(:,:,jj),self%mesh%np, &
          vals,1,0.d0,cc_vals(jj,1),3)
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
    IF(sensors%njumpers+self%nholes+self%n_vcoils>0)THEN
      DO j=1,sensors%njumpers
        jumpout(j+1)=sensors%jumpers(j)%compute_current(self,vals)
      END DO
      DO j=1,self%nholes+self%n_vcoils
        jumpout(sensors%njumpers+j+1)=tw_hole_value(self,vals,j)/mu0
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
IF(sensors%njumpers+self%nholes+self%n_vcoils>0)THEN
  CALL jumper_hist%close()
  DEALLOCATE(jumpout)
END IF
END IF
CALL u%delete
DEALLOCATE(u)
DEALLOCATE(vals)
IF(compute_B.AND.PRESENT(hodlr_op))THEN
  CALL Bx%delete()
  CALL By%delete()
  CALL Bz%delete()
  DEALLOCATE(Bx,By,Bz)
END IF
IF(self%n_icoils>0)DEALLOCATE(coil_vec)
CALL oft_decrease_indent
END SUBROUTINE plot_td_sim
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
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
    IF(.NOT.ASSOCIATED(self%Bel))CALL tw_compute_Bops(self,self%B_dx)
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