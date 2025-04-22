!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file oft_arpack.F90
!
!> Interface to parallel ARPACK for eigenvalue problems
!! - IR-Arnoldi (General)
!! - IR-Lanczos (Symmetric)
!!
!! @authors Chris Hansen
!! @date June 2012
!! @ingroup doxy_oft_lin_alg
!------------------------------------------------------------------------------
MODULE oft_arpack
#ifdef HAVE_ARPACK
USE oft_base
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver, oft_eigsolver, oft_solver_bc
IMPLICIT NONE
#include "local.h"
PRIVATE
!------------------------------------------------------------------------------
!> Implicity Restarted Arnoldi Method
!------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_eigsolver) :: oft_iram_eigsolver
  INTEGER(i4) :: mode = 1 !< Operational mode
  INTEGER(i4) :: nev = 1 !< Number of eigenvalues to compute
  INTEGER(i4) :: ncv = 0 !< Size of Ritz space (determined in call to `apply`)
  INTEGER(i4) :: info = 0 !< Solver status/return code
  REAL(r8) :: tol = 1.E-10_r8 !< Solver tolerance
  REAL(r8), POINTER, DIMENSION(:,:) :: eig_val => NULL() !< Eigenvalues
  REAL(r8), POINTER, DIMENSION(:,:) :: eig_vec => NULL() !< Eigenvectors
  CHARACTER(LEN=2) :: which = 'LM' !< Spectrum search flag
  CLASS(oft_solver), POINTER :: Minv => NULL() !< RHS inversion operator
  class(oft_solver_bc), pointer :: bc => NULL() !< Boundary condition
  CLASS(oft_solver_bc), POINTER :: orthog => NULL() !< Orthogonalization
CONTAINS
  !> Solve system
  PROCEDURE :: apply => iram_eig_apply
  !> Solve system
  PROCEDURE :: max => iram_eig_max
  !> Clean-up internal storage
  PROCEDURE :: delete => iram_delete
END TYPE oft_iram_eigsolver
!------------------------------------------------------------------------------
!> Implicity Restarted Lanczos Method
!------------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_eigsolver) :: oft_irlm_eigsolver
  INTEGER(i4) :: mode = 1 !< Operational mode
  INTEGER(i4) :: nev = 1 !< Number of eigenvalues to compute
  INTEGER(i4) :: ncv = 0 !< Size of Ritz space (determined in call to `apply`)
  INTEGER(i4) :: info = 0 !< Solver status/return code
  REAL(r8) :: tol = 1.E-10_r8 !< Solver tolerance
  REAL(r8), POINTER, DIMENSION(:,:) :: eig_val => NULL() !< Eigenvalues
  REAL(r8), POINTER, DIMENSION(:,:) :: eig_vec => NULL() !< Eigenvectors
  CHARACTER(LEN=2) :: which = 'LM' !< Spectrum search flag
  CLASS(oft_solver), POINTER :: Minv => NULL() !< RHS inversion operator
  class(oft_solver_bc), pointer :: bc => NULL() !< Boundary condition
  CLASS(oft_solver_bc), POINTER :: orthog => NULL() !< Orthogonalization
CONTAINS
  !> Solve system
  PROCEDURE :: apply => irlm_eig_apply
  !> Solve system
  PROCEDURE :: max => irlm_eig_max
  !> Clean-up internal storage
  PROCEDURE :: delete => irlm_delete
END TYPE oft_irlm_eigsolver
INTERFACE
!------------------------------------------------------------------------------
!> Interface to pdnaupd from ARPACK
!!
!! Driver subroutine for non-symmetric eigenvalue problems using ARPACK's
!! Implicitly Restarted Arnoldi Method.
!------------------------------------------------------------------------------
  SUBROUTINE dnaupd(ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr, &
                    workd,workl,lworkl,info)
  INTEGER :: ido
  CHARACTER(LEN=1) :: bmat
  INTEGER :: n
  CHARACTER(LEN=2) :: which
  INTEGER :: nev
  DOUBLE PRECISION :: tol
  DOUBLE PRECISION, DIMENSION(n) :: resid
  INTEGER :: ncv
  DOUBLE PRECISION, DIMENSION(ldv,ncv) :: v
  INTEGER :: ldv
  INTEGER, DIMENSION(11) :: iparam
  INTEGER, DIMENSION(14) :: ipntr
  DOUBLE PRECISION, DIMENSION(3*n) :: workd
  DOUBLE PRECISION, DIMENSION(lworkl) :: workl
  INTEGER :: lworkl
  INTEGER :: info
  END SUBROUTINE dnaupd
!------------------------------------------------------------------------------
!> Interface to pdsaupd from ARPACK
!------------------------------------------------------------------------------
  SUBROUTINE dsaupd(ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr, &
                    workd,workl,lworkl,info)
  INTEGER :: ido
  CHARACTER(LEN=1) :: bmat
  INTEGER :: n
  CHARACTER(LEN=2) :: which
  INTEGER :: nev
  DOUBLE PRECISION :: tol
  DOUBLE PRECISION, DIMENSION(n) :: resid
  INTEGER :: ncv
  DOUBLE PRECISION, DIMENSION(ldv,ncv) :: v
  INTEGER :: ldv
  INTEGER, DIMENSION(11) :: iparam
  INTEGER, DIMENSION(11) :: ipntr
  DOUBLE PRECISION, DIMENSION(3*n) :: workd
  DOUBLE PRECISION, DIMENSION(lworkl) :: workl
  INTEGER :: lworkl
  INTEGER :: info
  END SUBROUTINE dsaupd
!------------------------------------------------------------------------------
!> Interface to dneupd from ARPACK
!------------------------------------------------------------------------------
  SUBROUTINE dneupd(rvec,howmny,select,dr,di,z,ldz,sigmar,sigmai,workev, &
                    bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl, &
                    lworkl,info)
  LOGICAL :: rvec
  CHARACTER(LEN=1) :: howmny
  LOGICAL, DIMENSION(ncv) :: select
  DOUBLE PRECISION, DIMENSION(nev+1) :: dr
  DOUBLE PRECISION, DIMENSION(nev+1) :: di
  DOUBLE PRECISION, DIMENSION(ldz,*) :: z
  INTEGER :: ldz
  DOUBLE PRECISION :: sigmar
  DOUBLE PRECISION :: sigmai
  DOUBLE PRECISION, DIMENSION(3*ncv) :: workev
  CHARACTER(LEN=1) :: bmat
  INTEGER :: n
  CHARACTER(LEN=2) :: which
  INTEGER :: nev
  DOUBLE PRECISION :: tol
  DOUBLE PRECISION, DIMENSION(n) :: resid
  INTEGER :: ncv
  DOUBLE PRECISION, DIMENSION(ldv,ncv) :: v
  INTEGER :: ldv
  INTEGER, DIMENSION(11) :: iparam
  INTEGER, DIMENSION(14) :: ipntr
  DOUBLE PRECISION, DIMENSION(3*n) :: workd
  DOUBLE PRECISION, DIMENSION(lworkl) :: workl
  INTEGER :: lworkl
  INTEGER :: info
  END SUBROUTINE dneupd
!------------------------------------------------------------------------------
!> Interface to dseupd from ARPACK
!------------------------------------------------------------------------------
  SUBROUTINE dseupd(rvec,howmny,select,d,z,ldz,sigma,bmat,n,which,nev,tol, &
                    resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
  LOGICAL :: rvec
  CHARACTER(LEN=1) :: howmny
  LOGICAL, DIMENSION(ncv) :: select
  DOUBLE PRECISION, DIMENSION(nev) :: d
  DOUBLE PRECISION, DIMENSION(ldz,*) :: z
  INTEGER :: ldz
  DOUBLE PRECISION :: sigma
  CHARACTER(LEN=1) :: bmat
  INTEGER :: n
  CHARACTER(LEN=2) :: which
  INTEGER :: nev
  DOUBLE PRECISION :: tol
  DOUBLE PRECISION, DIMENSION(n) :: resid
  INTEGER :: ncv
  DOUBLE PRECISION, DIMENSION(ldv,ncv) :: v
  INTEGER :: ldv
  INTEGER, DIMENSION(11) :: iparam
  INTEGER, DIMENSION(11) :: ipntr
  DOUBLE PRECISION, DIMENSION(3*n) :: workd
  DOUBLE PRECISION, DIMENSION(lworkl) :: workl
  INTEGER :: lworkl
  INTEGER :: info
  END SUBROUTINE dseupd
#ifdef HAVE_MPI
!------------------------------------------------------------------------------
!> Interface to pdnaupd from ARPACK
!!
!! Driver subroutine for non-symmetric eigenvalue problems using ARPACK's
!! Implicitly Restarted Arnoldi Method.
!------------------------------------------------------------------------------
  SUBROUTINE pdnaupd(comm,ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam, &
                     ipntr,workd,workl,lworkl,info)
  INTEGER :: comm
  INTEGER :: ido
  CHARACTER(LEN=1) :: bmat
  INTEGER :: n
  CHARACTER(LEN=2) :: which
  INTEGER :: nev
  DOUBLE PRECISION :: tol
  DOUBLE PRECISION, DIMENSION(n) :: resid
  INTEGER :: ncv
  DOUBLE PRECISION, DIMENSION(ldv,ncv) :: v
  INTEGER :: ldv
  INTEGER, DIMENSION(11) :: iparam
  INTEGER, DIMENSION(14) :: ipntr
  DOUBLE PRECISION, DIMENSION(3*n) :: workd
  DOUBLE PRECISION, DIMENSION(lworkl) :: workl
  INTEGER :: lworkl
  INTEGER :: info
  END SUBROUTINE pdnaupd
!------------------------------------------------------------------------------
!> Interface to pdsaupd from ARPACK
!------------------------------------------------------------------------------
  SUBROUTINE pdsaupd(comm,ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam, &
                     ipntr,workd,workl,lworkl,info)
  INTEGER :: comm
  INTEGER :: ido
  CHARACTER(LEN=1) :: bmat
  INTEGER :: n
  CHARACTER(LEN=2) :: which
  INTEGER :: nev
  DOUBLE PRECISION :: tol
  DOUBLE PRECISION, DIMENSION(n) :: resid
  INTEGER :: ncv
  DOUBLE PRECISION, DIMENSION(ldv,ncv) :: v
  INTEGER :: ldv
  INTEGER, DIMENSION(11) :: iparam
  INTEGER, DIMENSION(11) :: ipntr
  DOUBLE PRECISION, DIMENSION(3*n) :: workd
  DOUBLE PRECISION, DIMENSION(lworkl) :: workl
  INTEGER :: lworkl
  INTEGER :: info
  END SUBROUTINE pdsaupd
!------------------------------------------------------------------------------
!> Interface to pdneupd from ARPACK
!------------------------------------------------------------------------------
  SUBROUTINE pdneupd(comm,rvec,howmny,select,dr,di,z,ldz,sigmar,sigmai,workev, &
                    bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl, &
                    lworkl,info)
  INTEGER :: comm
  LOGICAL :: rvec
  CHARACTER(LEN=1) :: howmny
  LOGICAL, DIMENSION(ncv) :: select
  DOUBLE PRECISION, DIMENSION(nev+1) :: dr
  DOUBLE PRECISION, DIMENSION(nev+1) :: di
  DOUBLE PRECISION, DIMENSION(ldz,*) :: z
  INTEGER :: ldz
  DOUBLE PRECISION :: sigmar
  DOUBLE PRECISION :: sigmai
  DOUBLE PRECISION, DIMENSION(3*ncv) :: workev
  CHARACTER(LEN=1) :: bmat
  INTEGER :: n
  CHARACTER(LEN=2) :: which
  INTEGER :: nev
  DOUBLE PRECISION :: tol
  DOUBLE PRECISION, DIMENSION(n) :: resid
  INTEGER :: ncv
  DOUBLE PRECISION, DIMENSION(ldv,ncv) :: v
  INTEGER :: ldv
  INTEGER, DIMENSION(11) :: iparam
  INTEGER, DIMENSION(14) :: ipntr
  DOUBLE PRECISION, DIMENSION(3*n) :: workd
  DOUBLE PRECISION, DIMENSION(lworkl) :: workl
  INTEGER :: lworkl
  INTEGER :: info
  END SUBROUTINE pdneupd
!------------------------------------------------------------------------------
!> Interface to pdseupd from ARPACK
!------------------------------------------------------------------------------
  SUBROUTINE pdseupd(comm,rvec,howmny,select,d,z,ldz,sigma, &
                    bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl, &
                    lworkl,info)
  INTEGER :: comm
  LOGICAL :: rvec
  CHARACTER(LEN=1) :: howmny
  LOGICAL, DIMENSION(ncv) :: select
  DOUBLE PRECISION, DIMENSION(nev) :: d
  DOUBLE PRECISION, DIMENSION(ldz,*) :: z
  INTEGER :: ldz
  DOUBLE PRECISION :: sigma
  CHARACTER(LEN=1) :: bmat
  INTEGER :: n
  CHARACTER(LEN=2) :: which
  INTEGER :: nev
  DOUBLE PRECISION :: tol
  DOUBLE PRECISION, DIMENSION(n) :: resid
  INTEGER :: ncv
  DOUBLE PRECISION, DIMENSION(ldv,ncv) :: v
  INTEGER :: ldv
  INTEGER, DIMENSION(11) :: iparam
  INTEGER, DIMENSION(11) :: ipntr
  DOUBLE PRECISION, DIMENSION(3*n) :: workd
  DOUBLE PRECISION, DIMENSION(lworkl) :: workl
  INTEGER :: lworkl
  INTEGER :: info
  END SUBROUTINE pdseupd
#endif
END INTERFACE
CONTAINS
!------------------------------------------------------------------------------
!> Compute the eigenvalue and eigenvector of a matrix system (A*x = Lam*M*x) using
!! an Implicitly Restarted Arnoldi Iteration.
!!
!! Solver employs the ARPACK non-symmetric driver routine. The location of the
!! eigenvalue is set in the calling class.
!------------------------------------------------------------------------------
SUBROUTINE iram_eig_apply(self,u,alam)
CLASS(oft_iram_eigsolver), intent(inout) :: self
CLASS(oft_vector), intent(inout) :: u !< Guess field/Eigenvector
REAL(r8), intent(inout) :: alam !< Eigenvalue
INTEGER(i4) :: i,j,ind,nslice,neigs,comm
INTEGER(i4) :: ido,ldv,lworkl,info
INTEGER(i4), DIMENSION(:) :: iparam(11),ipntr(14)
INTEGER(i4), POINTER, DIMENSION(:) :: emap
REAL(r8) :: sigmar,sigmai
REAL(r8), POINTER, DIMENSION(:) :: vslice
REAL(r8), ALLOCATABLE, TARGET, DIMENSION(:) :: workd,workl
REAL(r8), ALLOCATABLE, DIMENSION(:) :: resid,dr,di,workev
REAL(r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: v
LOGICAL :: rvec,dist,pm_save
CHARACTER(LEN=1) :: bmat,howmny
LOGICAL, ALLOCATABLE, DIMENSION(:) :: select
CLASS(oft_vector), POINTER :: tmp1,tmp2
DEBUG_STACK_PUSH
!------------------------------------------------------------------------------
! Setup solver
!------------------------------------------------------------------------------
IF(oft_env%head_proc.AND.oft_env%pm)WRITE(*,*)'Starting IRA solver'
!---Set ARPACK parameters
IF(self%ncv==0)self%ncv=2*self%nev+1
SELECT CASE(self%mode)
CASE(1) ! Solver mode (A*x = lam*x **** w/ OP = A and B = I)
  bmat='I'
  iparam(7)=1
CASE(2) ! Solver mode (A*x = lam*M*x **** w/ OP = inv[M]*A and B = M)
  bmat='G'
  iparam(7)=2
CASE DEFAULT
  CALL oft_abort("Unknown solver mode","iram_eig_apply",__FILE__)
END SELECT
lworkl=3*(self%ncv**2)+6*self%ncv ! Length of workl array
ido=0                             ! Reverse communication variable
info=0                            ! Error reporting
iparam(1)=1                       ! Use exact shifts
IF(self%its>0)THEN                ! Maximum number of Arnoldi iterations
  iparam(3)=self%its
ELSE
  iparam(3)=self%nev*10
END IF
rvec=.TRUE.                       ! Return eigenvectors
howmny='A'                        ! Form of the basis (A = Compute Ritz vectors)
neigs=self%nev
!------------------------------------------------------------------------------
! Create local to global mapping
!------------------------------------------------------------------------------
NULLIFY(vslice)
CALL u%get_slice(vslice)
nslice=SIZE(vslice)
ldv=nslice
!---Create output arrays
IF(ASSOCIATED(self%eig_val))DEALLOCATE(self%eig_val)
IF(ASSOCIATED(self%eig_vec))DEALLOCATE(self%eig_vec)
ALLOCATE(self%eig_val(2,neigs),self%eig_vec(nslice,neigs))
self%eig_val=0.d0
self%eig_vec=0.d0
!------------------------------------------------------------------------------
! Create temporary work arrays
!------------------------------------------------------------------------------
!---APRACK workspace
ALLOCATE(resid(nslice),v(ldv,self%ncv),workd(3*nslice),workl(lworkl))
ALLOCATE(select(self%ncv),dr(self%nev+1),di(self%nev+1),workev(3*self%ncv))
!---OFT work vectors
CALL u%new(tmp1)
CALL u%new(tmp2)
!------------------------------------------------------------------------------
! Set initial guess
!------------------------------------------------------------------------------
resid=vslice
DEALLOCATE(vslice)
!---Setup slicing on work vectors
CALL tmp1%get_slice(vslice)
DEALLOCATE(vslice)
CALL tmp2%get_slice(vslice)
DEALLOCATE(vslice)
CALL u%get_slice(vslice)
DEALLOCATE(vslice)
!------------------------------------------------------------------------------
! Begin Reverse Communication
!------------------------------------------------------------------------------
#ifdef OFT_MPI_F08
comm = oft_env%comm%MPI_VAL
#else
comm = oft_env%comm
#endif
i=0
self%info=0
DO
  !---Call APRACK driver
  IF(u%ng==u%nslice)THEN
    CALL dnaupd(ido,bmat,nslice,self%which,self%nev,self%tol, &
                resid,self%ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
  ELSE
#ifdef HAVE_MPI
    CALL pdnaupd(comm,ido,bmat,nslice,self%which,self%nev,self%tol, &
                 resid,self%ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
#endif
  END IF
  i=i+1
  !IF(head_proc.AND.pm)WRITE(*,*)i,ido
  !---Handle user action request
  IF(ido==-1.OR.ido==1)THEN ! Apply y = OP*x
    !---Propogate local slice to full vector
    vslice=>workd(ipntr(1):ipntr(1)+nslice-1)
    CALL tmp1%restore_slice(vslice)
    !---Apply BC
    IF(ASSOCIATED(self%bc))CALL self%bc%apply(tmp1)
    IF(ASSOCIATED(self%orthog))CALL self%orthog%apply(tmp1)
    !---Compute y' = A*x
    pm_save=oft_env%pm; oft_env%pm=.FALSE.
    CALL self%A%apply(tmp1,tmp2)
    oft_env%pm=pm_save
    IF(ASSOCIATED(self%bc))CALL self%bc%apply(tmp2)
    !IF(ASSOCIATED(self%orthog))CALL self%orthog%apply(tmp2)
    CALL tmp1%set(0.d0)
    pm_save=oft_env%pm; oft_env%pm=.FALSE.
    !---Compute y = inv(M)*y'
    IF(bmat=='G')THEN
      CALL self%Minv%apply(tmp1,tmp2)
    ELSE
      CALL tmp1%add(0.d0,1.d0,tmp2)
    END IF
    oft_env%pm=pm_save
    !---Copy local slice into work vector
    vslice=>workd(ipntr(2):ipntr(2)+nslice-1)
    CALL tmp1%get_slice(vslice)
    CYCLE
  ELSE IF(ido==2)THEN ! Apply y = M*x
    !---Propogate local slice to full vector
    vslice=>workd(ipntr(1):ipntr(1)+nslice-1)
    CALL tmp1%restore_slice(vslice)
    !---Apply BC
    IF(ASSOCIATED(self%bc))CALL self%bc%apply(tmp1)
    IF(ASSOCIATED(self%orthog))CALL self%orthog%apply(tmp1)
    pm_save=oft_env%pm; oft_env%pm=.FALSE.
    !---Compute y = M*x
    IF(bmat=='G')THEN
      CALL self%M%apply(tmp1,tmp2)
    ELSE
      CALL tmp2%add(0.d0,1.d0,tmp1)
    END IF
    oft_env%pm=pm_save
    IF(ASSOCIATED(self%bc))CALL self%bc%apply(tmp2)
    !IF(ASSOCIATED(self%orthog))CALL self%orthog%apply(tmp2)
    !---Copy local slice into work vector
    vslice=>workd(ipntr(2):ipntr(2)+nslice-1)
    CALL tmp2%get_slice(vslice)
    CYCLE
  END IF
  !---Check for an error or completion
  IF(info<0)THEN ! Error in solve
    alam=0.d0
    self%info=info
    IF(oft_env%head_proc)WRITE(*,*)'Error in EV solve ',info
    EXIT
  ELSE ! Solver has converged
    self%info=info
    IF(oft_env%head_proc.AND.oft_env%pm)WRITE(*,*)'Fetching EVs',info
    !---Retrieve eigenvalue and eigenvector
    IF(u%ng==u%nslice)THEN
      CALL dneupd(rvec,'A',select,dr,di,v,ldv, &
                  sigmar,sigmai,workev,bmat,nslice,self%which,self%nev,self%tol, &
                  resid,self%ncv,v,ldv,iparam,ipntr,workd, &
                  workl,lworkl,info)
    ELSE
#ifdef HAVE_MPI
      CALL pdneupd(comm,rvec,'A',select,dr,di,v,ldv, &
                   sigmar,sigmai,workev,bmat,nslice,self%which,self%nev,self%tol, &
                   resid,self%ncv,v,ldv,iparam,ipntr,workd, &
                   workl,lworkl,info)
#endif
    END IF
    !---Zero non-converged eigenvalues
    IF(info==1)THEN
      IF(oft_env%head_proc)WRITE(*,*)'Only ',iparam(5),' of ',self%nev,' eigenvalues converged in IRAM solve'
      self%info=iparam(5)
      DO i=self%info+1,self%nev
        dr(i)=0.d0
        di(i)=0.d0
        v(:,i)=0.d0
      END DO
    END IF
    IF(oft_env%head_proc.AND.oft_env%pm)THEN
      WRITE(*,*)'Solve Complete:'
      WRITE(*,*)'  ',dr(1:self%nev)
      WRITE(*,*)'  ',di(1:self%nev)
    END IF
    !---
    ind=MAXLOC(dr(1:self%nev),DIM=1)
    !---Copy eigenvector to output
    vslice=>v(:,ind)
    CALL u%restore_slice(vslice)
    !---Copy eigenvalue to output
    alam=dr(ind)
    !---Save to persistent arrays
    self%eig_val(1,:)=dr(1:neigs)
    self%eig_val(2,:)=di(1:neigs)
    self%eig_vec=v(:,1:neigs)
    EXIT
  END IF
END DO
!------------------------------------------------------------------------------
! Cleanup workspace
!------------------------------------------------------------------------------
DEALLOCATE(resid,v,workd,workl)
DEALLOCATE(select,dr,di,workev)
CALL tmp1%delete
CALL tmp2%delete
DEBUG_STACK_POP
END SUBROUTINE iram_eig_apply
!------------------------------------------------------------------------------
!> Compute the largest 2 eigenvalues of a matrix system (A*x = Lam*M*x) using
!! an Implicitly Restarted Arnoldi Iteration.
!!
!! Solver employs the ARPACK non-symmetric driver routine. Currently the guess
!! field is only used to create work vectors and a random initialization is used
!! for the solver.
!------------------------------------------------------------------------------
SUBROUTINE iram_eig_max(self,u,alam)
CLASS(oft_iram_eigsolver), INTENT(inout) :: self
CLASS(oft_vector), intent(inout) :: u !< Guess field/Eigenvector
REAL(r8), intent(inout) :: alam !< Eigenvalue
INTEGER(i4) :: i,j,comm
INTEGER(i4) :: ido,nslice,nev,ncv,ldv,lworkl,info
INTEGER(i4), DIMENSION(:) :: iparam(11),ipntr(14)
REAL(r8) :: sigmar,sigmai
REAL(r8), POINTER, DIMENSION(:) :: vslice
REAL(r8), ALLOCATABLE, TARGET, DIMENSION(:) :: workd,workl
REAL(r8), ALLOCATABLE, DIMENSION(:) :: resid,dr,di,workev
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: v
LOGICAL :: rvec,pm_save
LOGICAL, ALLOCATABLE, DIMENSION(:) :: select
CHARACTER(LEN=1) :: bmat,howmny
CHARACTER(LEN=2) :: which
CLASS(oft_vector), POINTER :: tmp1,tmp2
DEBUG_STACK_PUSH
!------------------------------------------------------------------------------
! Setup solver
!------------------------------------------------------------------------------
IF(oft_env%head_proc.AND.oft_env%pm)WRITE(*,*)'Starting IRA solver'
!---Set ARPACK parameters
SELECT CASE(self%mode)
CASE(1) ! Solver mode (A*x = lam*x **** w/ OP = A and B = I)
  bmat='I'
  iparam(7)=1
CASE(2) ! Solver mode (A*x = lam*M*x **** w/ OP = inv[M]*A and B = M)
  bmat='G'
  iparam(7)=2
CASE DEFAULT
  CALL oft_abort("Unknown solver mode","",__FILE__)
END SELECT
nev=1                     ! Number of eigenvalues
ncv=nev+9                 ! Number of work vectors
ido=0                     ! Reverse communication variable
info=0                    ! Error reporting
iparam(1)=1               ! Use exact shifts
iparam(3)=50              ! Maximum number of Arnoldi iterations
which='LM'                ! Range to search for EVs
lworkl=3*(ncv**2)+6*ncv   ! Length of workl array
rvec=.FALSE.              ! Return eigenvectors?
howmny='A'                ! Form of the basis (A = Compute Ritz vectors)
!------------------------------------------------------------------------------
! Create local to global mapping
!------------------------------------------------------------------------------
NULLIFY(vslice)
CALL u%get_slice(vslice)
nslice=SIZE(vslice)
ldv=nslice
!------------------------------------------------------------------------------
! Create temporary work arrays
!------------------------------------------------------------------------------
!---APRACK workspace
ALLOCATE(resid(nslice),v(ldv,ncv),workd(3*nslice),workl(lworkl))
ALLOCATE(select(ncv),dr(nev+1),di(nev+1),workev(3*ncv))
!---OFT work vectors
CALL u%new(tmp1)
CALL u%new(tmp2)
!---Setup slicing on work vectors
DEALLOCATE(vslice)
CALL tmp1%get_slice(vslice)
DEALLOCATE(vslice)
CALL tmp2%get_slice(vslice)
DEALLOCATE(vslice)
!------------------------------------------------------------------------------
! Begin Reverse Communication
!------------------------------------------------------------------------------
#ifdef OFT_MPI_F08
comm = oft_env%comm%MPI_VAL
#else
comm = oft_env%comm
#endif
self%info=0
DO i=1,4000
  !---Call APRACK driver
  IF(u%ng==u%nslice)THEN
    CALL dnaupd(ido,bmat,nslice,which,nev,self%tol,resid,ncv,v, &
      ldv,iparam,ipntr,workd,workl,lworkl,info)
  ELSE
#ifdef HAVE_MPI
    CALL pdnaupd(comm,ido,bmat,nslice,which,nev,self%tol,resid,ncv,v, &
      ldv,iparam,ipntr,workd,workl,lworkl,info)
#endif
  END IF
  !---Handle user action request
  IF(ido==-1.OR.ido==1)THEN ! Apply y = OP*x
    !---Propogate local slice to full vector
    vslice=>workd(ipntr(1):ipntr(1)+nslice-1)
    CALL tmp1%restore_slice(vslice)
    !---Apply BC
    CALL self%bc%apply(tmp1)
    !---Compute y' = A*x
    CALL self%A%apply(tmp1,tmp2)
    CALL self%bc%apply(tmp2)
    CALL tmp1%set(0.d0)
    IF(bmat=='G')THEN
      !---Compute y = inv(M)*y'
      pm_save=oft_env%pm; oft_env%pm=.FALSE.
      CALL self%Minv%apply(tmp1,tmp2)
      oft_env%pm=pm_save
    ELSE
      CALL tmp1%add(0.d0,1.d0,tmp2)
    END IF
    !---Copy local slice into work vector
    vslice=>workd(ipntr(2):ipntr(2)+nslice-1)
    CALL tmp1%get_slice(vslice)
    CYCLE
  ELSE IF(ido==2)THEN ! Apply y = M*x
    !---Propogate local slice to full vector
    vslice=>workd(ipntr(1):ipntr(1)+nslice-1)
    CALL tmp1%restore_slice(vslice)
    !---Apply BC
    CALL self%bc%apply(tmp1)
    IF(bmat=='G')THEN
      !---Compute y = M*x
      CALL self%M%apply(tmp1,tmp2)
    ELSE
      CALL tmp2%add(0.d0,1.d0,tmp1)
    END IF
    CALL self%bc%apply(tmp2)
    !---Copy local slice into work vector
    vslice=>workd(ipntr(2):ipntr(2)+nslice-1)
    CALL tmp2%get_slice(vslice)
    CYCLE
  END IF
  !---Check for an error or completion
  IF(info<0)THEN ! Error in solve
    alam=0.d0
    self%info=info
    IF(oft_env%head_proc)WRITE(*,*)'Error in EV solve ',info
    EXIT
  ELSE ! Solver has converged
    self%info=info
    !---Retrieve eigenvalues
    IF(u%ng==u%nslice)THEN
      CALL dneupd(rvec,'A',select,dr,di,v,ldv, &
      sigmar,sigmai,workev,bmat,nslice,which,nev,self%tol, &
      resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
    ELSE
#ifdef HAVE_MPI
      CALL pdneupd(comm,rvec,'A',select,dr,di,v,ldv, &
        sigmar,sigmai,workev,bmat,nslice,which,nev,self%tol, &
        resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
#endif
    END IF
    !---Ouput solution
    alam=dr(1)
    EXIT
  END IF
END DO
!------------------------------------------------------------------------------
! Cleanup workspace
!------------------------------------------------------------------------------
DEALLOCATE(resid,v,workd,workl)
DEALLOCATE(select,dr,di,workev)
CALL tmp1%delete
CALL tmp2%delete
DEBUG_STACK_POP
END SUBROUTINE iram_eig_max
!------------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine iram_delete(self)
class(oft_iram_eigsolver), intent(inout) :: self
NULLIFY(self%A,self%M)
IF(ASSOCIATED(self%eig_val))DEALLOCATE(self%eig_val)
IF(ASSOCIATED(self%eig_vec))DEALLOCATE(self%eig_vec)
self%initialized=.FALSE.
end subroutine iram_delete
!------------------------------------------------------------------------------
!> Compute the eigenvalue and eigenvector of a matrix system (A*x = Lam*M*x) using
!! an Implicitly Restarted Lanczos Iteration.
!!
!! Solver employs the ARPACK symmetric driver routine `dsaupd`. The location of the
!! eigenvalue is set in the calling class.
!------------------------------------------------------------------------------
SUBROUTINE irlm_eig_apply(self,u,alam)
CLASS(oft_irlm_eigsolver), intent(inout) :: self
CLASS(oft_vector), intent(inout) :: u !< Guess field/Eigenvector
REAL(r8), intent(inout) :: alam !< Eigenvalue
INTEGER(i4) :: i,j,comm
INTEGER(i4) :: ido,nslice,ldv,lworkl,info
INTEGER(i4), DIMENSION(:) :: iparam(11),ipntr(11)
REAL(r8) :: sigma
REAL(r8), POINTER, DIMENSION(:) :: vslice
REAL(r8), ALLOCATABLE, TARGET, DIMENSION(:) :: workd,workl
REAL(r8), ALLOCATABLE, DIMENSION(:) :: resid,d,workev
REAL(r8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: v
LOGICAL :: rvec,pm_save
CHARACTER(LEN=1) :: bmat,howmny
LOGICAL, ALLOCATABLE, DIMENSION(:) :: select
CLASS(oft_vector), POINTER :: tmp1,tmp2
DEBUG_STACK_PUSH
!------------------------------------------------------------------------------
! Setup solver
!------------------------------------------------------------------------------
IF(oft_env%head_proc.AND.oft_env%pm)WRITE(*,*)'Starting IRL solver'
!---Set ARPACK parameters
SELECT CASE(self%mode)
CASE(1) ! Solver mode (A*x = lam*x **** w/ OP = A and B = I)
  bmat='I'
  iparam(7)=1
CASE(2) ! Solver mode (A*x = lam*M*x **** w/ OP = inv[M]*A and B = M)
  bmat='G'
  iparam(7)=2
CASE DEFAULT
  CALL oft_abort("Unknown solver mode","iram_eig_apply",__FILE__)
END SELECT
IF(self%ncv==0)self%ncv=2*self%nev+1
lworkl=3*(self%ncv**2)+6*self%ncv ! Length of workl array
ido=0                             ! Reverse communication variable
info=0                            ! Error reporting
iparam(1)=1                       ! Use exact shifts
IF(self%its>0)THEN                ! Maximum number of Lanzcos iterations
  iparam(3)=self%its
ELSE
  iparam(3)=self%nev*10
END IF
rvec=.TRUE.                       ! Return eigenvectors
howmny='A'                        ! Form of the basis (A = Compute Ritz vectors)
!------------------------------------------------------------------------------
! Create local to global mapping
!------------------------------------------------------------------------------
NULLIFY(vslice)
CALL u%get_slice(vslice)
nslice=SIZE(vslice)
ldv=nslice
!---Create output arrays
IF(ASSOCIATED(self%eig_val))DEALLOCATE(self%eig_val)
IF(ASSOCIATED(self%eig_vec))DEALLOCATE(self%eig_vec)
ALLOCATE(self%eig_val(2,self%nev),self%eig_vec(nslice,self%nev))
self%eig_val=0.d0
self%eig_vec=0.d0
!------------------------------------------------------------------------------
! Create temporary work arrays
!------------------------------------------------------------------------------
!---APRACK workspace
ALLOCATE(resid(nslice),v(ldv,self%ncv),workd(3*nslice),workl(lworkl))
ALLOCATE(select(self%ncv),d(self%nev+1),workev(3*self%ncv))
!---OFT work vectors
CALL u%new(tmp1)
CALL u%new(tmp2)
!------------------------------------------------------------------------------
! Set initial guess
!------------------------------------------------------------------------------
resid=vslice
DEALLOCATE(vslice)
!---Setup slicing on work vectors
CALL tmp1%get_slice(vslice)
DEALLOCATE(vslice)
CALL tmp2%get_slice(vslice)
DEALLOCATE(vslice)
CALL u%get_slice(vslice)
DEALLOCATE(vslice)
!------------------------------------------------------------------------------
! Begin Reverse Communication
!------------------------------------------------------------------------------
#ifdef OFT_MPI_F08
comm = oft_env%comm%MPI_VAL
#else
comm = oft_env%comm
#endif
i=0
self%info=0
DO
  !---Call APRACK driver
  IF(u%ng==u%nslice)THEN
    CALL dsaupd(ido,bmat,nslice,self%which,self%nev,self%tol,resid, &
      self%ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
  ELSE
#ifdef HAVE_MPI
    CALL pdsaupd(comm,ido,bmat,nslice,self%which,self%nev,self%tol,resid, &
      self%ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
#endif
  END IF
  i=i+1
  !IF(head_proc.AND.pm)WRITE(*,*)i,ido
  !---Handle user action request
  IF(ido==-1.OR.ido==1)THEN ! Apply y = OP*x
    !---Propogate local slice to full vector
    vslice=>workd(ipntr(1):ipntr(1)+nslice-1)
    CALL tmp1%restore_slice(vslice)
    !---Apply BC
    IF(ASSOCIATED(self%bc))CALL self%bc%apply(tmp1)
    !---Compute y' = A*x
    pm_save=oft_env%pm; oft_env%pm=.FALSE.
    CALL self%A%apply(tmp1,tmp2)
    oft_env%pm=pm_save
    IF(ASSOCIATED(self%bc))CALL self%bc%apply(tmp2)
    !---Copy local slice into work vector
    vslice=>workd(ipntr(1):ipntr(1)+nslice-1)
    CALL tmp2%get_slice(vslice)
    CALL tmp1%set(0.d0)
    pm_save=oft_env%pm; oft_env%pm=.FALSE.
    !---Compute y = inv(M)*y'
    IF(bmat=='G')THEN
      CALL self%Minv%apply(tmp1,tmp2)
    ELSE
      CALL tmp1%add(0.d0,1.d0,tmp2)
    END IF
    oft_env%pm=pm_save
    !---Copy local slice into work vector
    vslice=>workd(ipntr(2):ipntr(2)+nslice-1)
    CALL tmp1%get_slice(vslice)
    CYCLE
  ELSE IF(ido==2)THEN ! Apply y = M*x
    !---Propogate local slice to full vector
    vslice=>workd(ipntr(1):ipntr(1)+nslice-1)
    CALL tmp1%restore_slice(vslice)
    !---Apply BC
    IF(ASSOCIATED(self%bc))CALL self%bc%apply(tmp1)
    pm_save=oft_env%pm; oft_env%pm=.FALSE.
    !---Compute y = M*x
    IF(bmat=='G')THEN
      CALL self%M%apply(tmp1,tmp2)
    ELSE
      CALL tmp2%add(0.d0,1.d0,tmp1)
    END IF
    oft_env%pm=pm_save
    IF(ASSOCIATED(self%bc))CALL self%bc%apply(tmp2)
    !---Copy local slice into work vector
    vslice=>workd(ipntr(2):ipntr(2)+nslice-1)
    CALL tmp2%get_slice(vslice)
    CYCLE
  END IF
  !---Check for an error or completion
  IF(info<0)THEN ! Error in solve
    alam=0.d0
    self%info=info
    IF(oft_env%head_proc)WRITE(*,*)'Error in EV solve ',info
    EXIT
  ELSE ! Solver has converged
    self%info=info
    IF(oft_env%head_proc.AND.oft_env%pm)WRITE(*,*)'Fetching EVs'
    !---Retrieve eigenvalue and eigenvector
    IF(u%ng==u%nslice)THEN
      CALL dseupd(rvec,'A',select,d,v,ldv, &
        sigma,bmat,nslice,self%which,self%nev,self%tol, &
        resid,self%ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
    ELSE
#ifdef HAVE_MPI
      CALL pdseupd(comm,rvec,'A',select,d,v,ldv, &
        sigma,bmat,nslice,self%which,self%nev,self%tol, &
        resid,self%ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,info)
#endif
    END IF
    !---Zero non-converged eigenvalues
    IF(info==1)THEN
      IF(oft_env%head_proc)WRITE(*,*)'Only ',iparam(5),' of ',self%nev,' eigenvalues converged in IRLM solve'
      self%info=iparam(5)
      DO i=self%info+1,self%nev
        d(i)=0.d0
        v(:,i)=0.d0
      END DO
    END IF
    IF(oft_env%head_proc.AND.oft_env%pm)THEN
      WRITE(*,*)'Solve Complete:'
      WRITE(*,*)'  ',d(1:self%nev)
    END IF
    !---Copy eigenvector to output
    vslice=>v(:,1)
    CALL u%restore_slice(vslice)
    !---Copy eigenvalue to output
    alam=d(1)
    !---Save to persistent arrays
    self%eig_val(1,:)=d(1:self%nev)
    self%eig_vec=v(:,1:self%nev)
    IF(oft_env%head_proc.AND.oft_env%pm)WRITE(*,*)'Solve Complete',d(1:self%nev)
    EXIT
  END IF
END DO
!------------------------------------------------------------------------------
! Cleanup workspace
!------------------------------------------------------------------------------
DEALLOCATE(resid,v,workd,workl)
DEALLOCATE(select,d,workev)
CALL tmp1%delete
CALL tmp2%delete
DEBUG_STACK_POP
END SUBROUTINE irlm_eig_apply
!------------------------------------------------------------------------------
!> Compute the largest 2 eigenvalues of a matrix system (A*x = Lam*M*x) using
!! an Implicitly Restarted Lanczos Iteration.
!!
!! Solver employs the ARPACK symmetric driver routine. Currently the guess
!! field is only used to create work vectors and a random initialization is used
!! for the solver.
!------------------------------------------------------------------------------
SUBROUTINE irlm_eig_max(self,u,alam)
CLASS(oft_irlm_eigsolver), INTENT(inout) :: self
CLASS(oft_vector), intent(inout) :: u !< Guess field/Eigenvector
REAL(r8), intent(inout) :: alam !< Eigenvalue
INTEGER(i4) :: i,j,comm
INTEGER(i4) :: ido,nslice,nev,ncv,ldv,lworkl,info
INTEGER(i4), DIMENSION(:) :: iparam(11),ipntr(11)
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: emap
REAL(r8) :: sigma
REAL(r8), POINTER, DIMENSION(:) :: vslice
REAL(r8), ALLOCATABLE, TARGET, DIMENSION(:) :: workd,workl
REAL(r8), ALLOCATABLE, DIMENSION(:) :: resid,d,workev
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: v
LOGICAL :: rvec,dist,pm_save
LOGICAL, ALLOCATABLE, DIMENSION(:) :: select
CHARACTER(LEN=1) :: bmat,howmny
CHARACTER(LEN=2) :: which
CLASS(oft_vector), POINTER :: tmp1,tmp2
DEBUG_STACK_PUSH
!------------------------------------------------------------------------------
! Setup solver
!------------------------------------------------------------------------------
IF(oft_env%head_proc.AND.oft_env%pm)WRITE(*,*)'Starting IRL solver'
!---Set ARPACK parameters
SELECT CASE(self%mode)
CASE(1) ! Solver mode (A*x = lam*x **** w/ OP = A and B = I)
  bmat='I'
  iparam(7)=1
CASE(2) ! Solver mode (A*x = lam*M*x **** w/ OP = inv[M]*A and B = M)
  bmat='G'
  iparam(7)=2
CASE DEFAULT
  CALL oft_abort("Unknown solver mode","irlm_max_apply",__FILE__)
END SELECT
nev=1                     ! Number of eigenvalues
ncv=nev+9                 ! Number of work vectors
ido=0                     ! Reverse communication variable
info=0                    ! Error reporting
iparam(1)=1               ! Use exact shifts
iparam(3)=200             ! Maximum number of Arnoldi iterations
which='LM'                ! Range to search for EVs
lworkl=(ncv**2)+8*ncv     ! Length of workl array
rvec=.FALSE.              ! Return eigenvectors
howmny='A'                ! Form of the basis (A = Compute Ritz vectors)
!------------------------------------------------------------------------------
! Create local to global mapping
!------------------------------------------------------------------------------
NULLIFY(vslice)
CALL u%get_slice(vslice)
nslice=SIZE(vslice)
ldv=nslice
!------------------------------------------------------------------------------
! Create temporary work arrays
!------------------------------------------------------------------------------
!---APRACK workspace
ALLOCATE(resid(nslice),v(ldv,ncv),workd(3*nslice),workl(lworkl))
ALLOCATE(select(ncv),d(nev))
!---OFT work vectors
CALL u%new(tmp1)
CALL u%new(tmp2)
!---Setup slicing on work vectors
DEALLOCATE(vslice)
CALL tmp1%get_slice(vslice)
DEALLOCATE(vslice)
CALL tmp2%get_slice(vslice)
DEALLOCATE(vslice)
!------------------------------------------------------------------------------
! Begin Reverse Communication
!------------------------------------------------------------------------------
#ifdef OFT_MPI_F08
comm = oft_env%comm%MPI_VAL
#else
comm = oft_env%comm
#endif
self%info=0
DO
  !---Call APRACK driver
  IF(u%ng==u%nslice)THEN
    CALL dsaupd(ido,bmat,nslice,which,nev,self%tol,resid,ncv, &
      v,ldv,iparam,ipntr,workd,workl,lworkl,info)
  ELSE
#ifdef HAVE_MPI
    CALL pdsaupd(comm,ido,bmat,nslice,which,nev,self%tol,resid,ncv, &
      v,ldv,iparam,ipntr,workd,workl,lworkl,info)
#endif
  END IF
  !---Handle user action request
  IF(ido==-1.OR.ido==1)THEN ! Apply y = OP*x
    !---Propogate local slice to full vector
    vslice=>workd(ipntr(1):ipntr(1)+nslice-1)
    CALL tmp1%restore_slice(vslice)
    !---Apply BC
    IF(ASSOCIATED(self%bc))CALL self%bc%apply(tmp1)
    !---Compute y' = A*x
    CALL self%A%apply(tmp1,tmp2)
    IF(ASSOCIATED(self%bc))CALL self%bc%apply(tmp2)
    !IF(ASSOCIATED(self%orthog))CALL self%orthog%apply(tmp2)
    !---Copy local slice into work vector
    vslice=>workd(ipntr(1):ipntr(1)+nslice-1)
    CALL tmp2%get_slice(vslice)
    CALL tmp1%set(0.d0)
    IF(bmat=='G')THEN
      !---Compute y = inv(M)*y'
      pm_save=oft_env%pm; oft_env%pm=.FALSE.
      CALL self%Minv%apply(tmp1,tmp2)
      oft_env%pm=pm_save
    ELSE
      CALL tmp1%add(0.d0,1.d0,tmp2)
    END IF
    IF(ASSOCIATED(self%orthog))CALL self%orthog%apply(tmp1)
    !oft_env%pm=pm_save
    !---Copy local slice into work vector
    vslice=>workd(ipntr(2):ipntr(2)+nslice-1)
    CALL tmp1%get_slice(vslice)
    CYCLE
  ELSE IF(ido==2)THEN ! Apply y = M*x
    !---Propogate local slice to full vector
    vslice=>workd(ipntr(1):ipntr(1)+nslice-1)
    CALL tmp1%restore_slice(vslice)
    !---Apply BC
    IF(ASSOCIATED(self%bc))CALL self%bc%apply(tmp1)
    !---Compute y = M*x
    CALL self%M%apply(tmp1,tmp2)
    IF(bmat=='G')THEN
      !---Compute y = M*x
      CALL self%M%apply(tmp1,tmp2)
    ELSE
      CALL tmp2%add(0.d0,1.d0,tmp1)
    END IF
    IF(ASSOCIATED(self%bc))CALL self%bc%apply(tmp2)
    IF(ASSOCIATED(self%orthog))CALL self%orthog%apply(tmp2)
    !---Copy local slice into work vector
    vslice=>workd(ipntr(2):ipntr(2)+nslice-1)
    CALL tmp2%get_slice(vslice)
    CYCLE
  END IF
  !---Check for an error or completion
  IF(info<0)THEN ! Error in solve
    alam=0.d0
    self%info=info
    IF(oft_env%head_proc)WRITE(*,*)'Error in EV solve ',info
    EXIT
  ELSE ! Solver has converged
    self%info=info
    IF(info==1)THEN
      IF(oft_env%head_proc)WRITE(*,*)'Maximum number of Iterations reached'
    END IF
    !---Retrieve eigenvalues
    IF(u%ng==u%nslice)THEN
      CALL dseupd(rvec,'A',select,d,v(:,1:nev),ldv, &
        sigma,bmat,nslice,which,nev,self%tol,resid,ncv,v,ldv, &
        iparam,ipntr,workd,workl,lworkl,info)
    ELSE
#ifdef HAVE_MPI
      CALL pdseupd(comm,rvec,'A',select,d,v(:,1:nev),ldv, &
        sigma,bmat,nslice,which,nev,self%tol,resid,ncv,v,ldv, &
        iparam,ipntr,workd,workl,lworkl,info)
#endif
    END IF
    !---Ouput solution
    alam=d(1)
    EXIT
  END IF
END DO
!------------------------------------------------------------------------------
! Cleanup workspace
!------------------------------------------------------------------------------
DEALLOCATE(resid,v,workd,workl)
DEALLOCATE(select,d)
CALL tmp1%delete
CALL tmp2%delete
DEBUG_STACK_POP
END SUBROUTINE irlm_eig_max
!------------------------------------------------------------------------------
!> Destroy diagonal preconditioner and deallocate all internal storage
!------------------------------------------------------------------------------
subroutine irlm_delete(self)
class(oft_irlm_eigsolver), intent(inout) :: self
NULLIFY(self%A,self%M)
IF(ASSOCIATED(self%eig_val))DEALLOCATE(self%eig_val)
IF(ASSOCIATED(self%eig_vec))DEALLOCATE(self%eig_vec)
self%initialized=.FALSE.
end subroutine irlm_delete
#endif
END MODULE oft_arpack
