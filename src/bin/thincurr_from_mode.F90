!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file thincurr_from_mode.F90
!
!> Project DCON mode on to ThinCurr triangular grid for thin-wall studies
!!
!! **Option group:** `thincurr_from_mode_options`
!! |  Option                 |  Description  | Type [dim] |
!! |-------------------------|---------------|------------|
!! |  `filename="none"`      |  DCON mode definition file | str |
!! |  `ntheta=40`            |  Number of poloidal grid points | int |
!! |  `nphi=90`              |  Number of toroidal grid points | int |
!!
!! @authors Chris Hansen
!! @date May 2017
!! @ingroup doxy_thincurr
!---------------------------------------------------------------------------
PROGRAM thincurr_from_mode
USE oft_base
USE spline_mod
USE oft_sort, ONLY: sort_array
USE oft_mesh_type, ONLY: oft_bmesh
USE oft_trimesh_type, ONLY: oft_trimesh
USE oft_la_base, ONLY: oft_vector
USE oft_native_la, ONLY: oft_native_dense_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
USE thin_wall
IMPLICIT NONE
#include "local.h"
INTEGER(4) :: ndrivers = 0
INTEGER(4) :: nsensors = 0
TYPE(tw_type) :: tw_sim
TYPE(oft_timer) :: mytimer
INTEGER(4) :: i,j,io_unit,ierr
REAL(8) :: c1(2),c2,tmp(2)
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: vtmp,utmp,senout
REAL(8), POINTER, DIMENSION(:,:) :: bnorm
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: hole_nsets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: jumper_nsets => NULL()
TYPE(tw_sensors) :: sensors
!---
CHARACTER(LEN=OFT_PATH_SLEN) :: filename = 'none'
INTEGER(4) :: ntheta = 40
INTEGER(4) :: nphi = 90
LOGICAL :: use_spline = .TRUE.
CHARACTER(LEN=7) :: resample_type='arc_len'
NAMELIST/thincurr_from_mode_options/filename,ntheta,nphi,use_spline,resample_type
!---
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,thincurr_from_mode_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(TRIM(filename)=="none")CALL oft_abort("No mode filename specified","thincurr_from_mode",__FILE__)
!---Setup mesh
CALL setup_mode_mesh(TRIM(filename),ntheta,nphi,tw_sim%mesh,bnorm)
! Holes
ALLOCATE(hole_nsets(2))
DO i=1,2
  CALL get_torus_loop(tw_sim%mesh, i, hole_nsets(i)%v, hole_nsets(i)%n)
END DO
! Closures
tw_sim%nclosures=1
ALLOCATE(tw_sim%closures(tw_sim%nclosures))
tw_sim%closures(1)=1 !INT(tw_sim%mesh%nc/2.d0,4)
CALL tw_sim%setup(hole_nsets)
!---Setup I/0
CALL tw_sim%xdmf%setup("ThinCurr")
CALL tw_sim%mesh%setup_io(tw_sim%xdmf,1)
!---Compute face mutuals
CALL tw_compute_LmatDirect(tw_sim,tw_sim%Lmat)
!---Compute and correct perturbation
ALLOCATE(vtmp(tw_sim%nelems,2))
vtmp=0.d0
!$omp parallel do private(j)
DO i=1,tw_sim%np_active
  j=tw_sim%lpmap_inv(i)
  DO j=tw_sim%kpmap_inv(i),tw_sim%kpmap_inv(i+1)-1
    vtmp(i,:)=bnorm(tw_sim%lpmap_inv(j),:)*tw_sim%mesh%va(tw_sim%lpmap_inv(j))
  END DO
END DO
!---Save B-normal field
WRITE(*,*)'Flux chk',SUM(vtmp(:,1)),SUM(vtmp(:,2))
CALL tw_sim%mesh%save_vertex_scalar(bnorm(:,1),tw_sim%xdmf,'BSin')
CALL tw_sim%mesh%save_vertex_scalar(bnorm(:,2),tw_sim%xdmf,'BCos')
DEALLOCATE(bnorm)
WRITE(*,*)'Flux chk',SUM(vtmp(:,1)),SUM(vtmp(:,2))
!---Compute current potential
ALLOCATE(utmp(tw_sim%nelems,2))
utmp=0.d0
CALL solve_inv(tw_sim,utmp,vtmp)
CALL tw_save_pfield(tw_sim,utmp(:,1),'JSin')
CALL tw_save_pfield(tw_sim,utmp(:,2),'JCos')
!---Save to file
CALL hdf5_create_file('tCurr_mode_model.h5')
CALL hdf5_create_group('tCurr_mode_model.h5','mesh')
CALL hdf5_write(tw_sim%mesh%r,'tCurr_mode_model.h5','mesh/R')
CALL hdf5_write(tw_sim%mesh%lc,'tCurr_mode_model.h5','mesh/LC')
CALL hdf5_write(1,'tCurr_mode_model.h5','mesh/NUM_SIDESETS')
CALL hdf5_write(tw_sim%closures(1),'tCurr_mode_model.h5','mesh/SIDESET0001')
CALL hdf5_write(tw_sim%nholes,'tCurr_mode_model.h5','mesh/NUM_NODESETS')
CALL hdf5_write(tw_sim%hmesh(1)%lp,'tCurr_mode_model.h5','mesh/NODESET0001')
CALL hdf5_write(tw_sim%hmesh(2)%lp,'tCurr_mode_model.h5','mesh/NODESET0002')
CALL hdf5_create_group('tCurr_mode_model.h5','thincurr')
CALL hdf5_write(utmp,'tCurr_mode_model.h5','thincurr/driver')
! OPEN(NEWUNIT=io_unit, FILE='thincurr_mode_model.dat')
! WRITE(io_unit,*)tw_sim%mesh%np
! DO i=1,tw_sim%mesh%np
!   WRITE(io_unit,*)tw_sim%mesh%r(:,i)
! END DO
! WRITE(io_unit,*)tw_sim%closures(1)
! WRITE(io_unit,*)tw_sim%mesh%nc
! DO i=1,tw_sim%mesh%nc
!   WRITE(io_unit,*)tw_sim%mesh%lc(:,i)
! END DO
! WRITE(io_unit,*)tw_sim%nholes
! DO i=1,tw_sim%nholes
!   WRITE(io_unit,*)tw_sim%hmesh(i)%n
!   DO j=1,tw_sim%hmesh(i)%n
!     WRITE(io_unit,*)tw_sim%hmesh(i)%lp(j)
!   END DO
! END DO
! WRITE(io_unit,*)2
! DO i=1,tw_sim%nelems
!   WRITE(io_unit,*)utmp(i,:)
! END DO
! CLOSE(io_unit)
!---Save probe signals
CALL tw_load_sensors('floops.loc',tw_sim,sensors)
IF(sensors%nfloops>0)THEN
  tw_sim%n_icoils=0
  CALL tw_compute_mutuals(tw_sim,sensors%nfloops,sensors%floops)
  ALLOCATE(senout(sensors%nfloops,2))
  senout=0.d0
  CALL dgemv('N',sensors%nfloops,tw_sim%nelems,1.d0,tw_sim%Ael2sen,sensors%nfloops,utmp(:,1),1,0.d0,senout(:,1),1)
  CALL dgemv('N',sensors%nfloops,tw_sim%nelems,1.d0,tw_sim%Ael2sen,sensors%nfloops,utmp(:,2),1,0.d0,senout(:,2),1)
  OPEN(NEWUNIT=io_unit,FILE='thincurr_mode.dat')
  DO i=1,sensors%nfloops
    WRITE(io_unit,'(A,4Es24.15)')sensors%floops(i)%name,senout(i,:)
  END DO
  CLOSE(io_unit)
  DEALLOCATE(senout)
END IF
!---Cleanup
DEALLOCATE(utmp,vtmp)
DEALLOCATE(tw_sim%Lmat)
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE setup_mode_mesh
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE setup_mode_mesh(filename,nsample,nphi,bmesh,bnorm)
CHARACTER(LEN=*), INTENT(in) :: filename
INTEGER(4), INTENT(in) :: nsample,nphi
CLASS(oft_bmesh), POINTER, INTENT(out) :: bmesh
REAL(8), POINTER, INTENT(out) :: bnorm(:,:)
INTEGER(4) :: i,j,in,jn,io_unit,npts,io_stat
INTEGER(4), ALLOCATABLE :: ind_reorder(:)
REAL(8) :: theta,r0(2),rhat(2),nmode
REAL(8), ALLOCATABLE :: mode_in(:,:),mode_resample(:,:),resample_grid(:)
REAL(8), ALLOCATABLE :: phi_grid(:),mode_tmp(:,:),theta_tmp(:)
TYPE(spline_type) :: mode_spline
WRITE(*,'(2A)')oft_indent,'Loading plasma mode'
WRITE(*,'(3A)')oft_indent,'  filename = ',TRIM(filename)
!---Allocate mesh
ALLOCATE(oft_trimesh::bmesh)
CALL bmesh%setup(-1,.FALSE.)
!---Read DCON mode from modified "surf.out" file produced by match
OPEN(NEWUNIT=io_unit,FILE=TRIM(filename))
io_stat=skip_comment_lines(io_unit)
READ(io_unit,*)npts,nmode
io_stat=skip_comment_lines(io_unit)
ALLOCATE(mode_in(7,npts+2))
r0=0.d0
DO i=1,npts
  READ(io_unit,*)mode_in(1:3,i+1),mode_in(7,i+1)
  r0=r0+mode_in(1:2,i+1)
END DO
CLOSE(io_unit)
IF(MAXVAL(ABS(mode_in(:,npts+1)-mode_in(:,2)))<1.d-6)THEN
  npts=npts-1 ! First and last point are the same
  r0=r0-mode_in(1:2,npts+1) ! Remove one copy from center sum
END IF
r0=r0/REAL(npts,8)
WRITE(*,'(2A,I3)')oft_indent,'  N        = ',INT(nmode,4)
WRITE(*,'(2A,I6)')oft_indent,'  # of pts = ',npts
WRITE(*,'(2A,2ES12.4)')oft_indent,'  R0       = ',r0
!---Sort input grid to consistent ordering
ALLOCATE(theta_tmp(npts),ind_reorder(npts))
DO i=2,npts+1
  theta_tmp(i-1)=ATAN2(mode_in(2,i)-r0(2),mode_in(1,i)-r0(1))
  IF(theta_tmp(i-1)<0.d0)theta_tmp(i-1)=2.d0*pi+theta_tmp(i-1)
  in=i+1
  IF(i==npts+1)in=2
  rhat(1)=-(mode_in(2,in)-mode_in(2,i))
  rhat(2)=mode_in(1,in)-mode_in(1,i)
  rhat=rhat/SQRT(SUM(rhat**2))
  mode_in(5:6,i)=rhat
  ind_reorder(i-1)=i
END DO
IF(theta_tmp(1)>pi)theta_tmp(1)=2.d0*pi-theta_tmp(1)
CALL sort_array(theta_tmp,ind_reorder,npts)
ALLOCATE(mode_tmp(7,npts+2))
mode_tmp=mode_in
DO i=2,npts+1
  mode_in(:,i)=mode_tmp(:,ind_reorder(i-1))
  mode_in(4,i)=theta_tmp(i-1)
END DO
DEALLOCATE(theta_tmp,ind_reorder,mode_tmp)
mode_in(:,1)=mode_in(:,npts+1); mode_in(4,1)=mode_in(4,npts+1)-2.d0*pi
mode_in(:,npts+2)=mode_in(:,2); mode_in(4,npts+2)=mode_in(4,2)+2.d0*pi
!---Resample fields onto new poloidal grid
ALLOCATE(mode_resample(6,nsample),resample_grid(nsample))
mode_resample=0.d0
resample_grid=0.d0
!---Setup grid based on sampling type
SELECT CASE(resample_type)
CASE("theta")
  DO i=1,nsample
    resample_grid(i)=(i-1)*2.d0*pi/REAL(nsample,8)
  END DO
CASE("arc_len")
  mode_in(4,1)=0.d0
  DO i=1,npts+1
    mode_in(4,i+1)=mode_in(4,i)+SQRT(SUM((mode_in(1:2,i+1)-mode_in(1:2,i))**2))
  END DO
  DO i=1,nsample
    resample_grid(i)=(i-1)*mode_in(4,npts+1)/REAL(nsample,8)
  END DO
END SELECT
WRITE(*,*)
IF(use_spline)THEN
  CALL spline_alloc(mode_spline,npts+1,6)
  DO i=1,npts+2
    mode_spline%xs(i-1)=mode_in(4,i)
    mode_spline%fs(i-1,1:3)=mode_in(1:3,i)
    mode_spline%fs(i-1,4:6)=mode_in(5:7,i)
  END DO
  CALL spline_fit(mode_spline,"not-a-knot")
  DO i=1,nsample
    CALL spline_eval(mode_spline,resample_grid(i),0)
    mode_resample(:,i)=mode_spline%f
  END DO
  CALL spline_dealloc(mode_spline)
ELSE
  DO i=1,nsample
    mode_resample(1,i)=linterp(mode_in(4,:),mode_in(1,:),npts+2,resample_grid(i))
    mode_resample(2,i)=linterp(mode_in(4,:),mode_in(2,:),npts+2,resample_grid(i))
    mode_resample(3,i)=linterp(mode_in(4,:),mode_in(3,:),npts+2,resample_grid(i))
    mode_resample(4,i)=linterp(mode_in(4,:),mode_in(5,:),npts+2,resample_grid(i))
    mode_resample(5,i)=linterp(mode_in(4,:),mode_in(6,:),npts+2,resample_grid(i))
    mode_resample(6,i)=linterp(mode_in(4,:),mode_in(7,:),npts+2,resample_grid(i))
  END DO
END IF
!mode_resample(3,:)=mode_resample(3,:)-SUM(mode_resample(3,:))/REAL(nsample,8)
!mode_resample(6,:)=mode_resample(6,:)-SUM(mode_resample(6,:))/REAL(nsample,8)
WRITE(*,*)SUM(mode_resample(3,:)),SUM(mode_resample(6,:))
DEALLOCATE(mode_in)
ALLOCATE(phi_grid(nphi))
DO i=1,nphi
  phi_grid(i)=(i-1)*2.d0*pi/REAL(nphi,8)
END DO
!---Construct toroidal mesh and interpolate Bn to vertices
bmesh%np=nsample*nphi
bmesh%nc=2*nsample*nphi
ALLOCATE(bmesh%r(3,bmesh%np),bmesh%lc(3,bmesh%nc),bmesh%reg(bmesh%nc))
ALLOCATE(bnorm(bmesh%np,2))
bmesh%reg=1
DO i=1,nsample
  in=i+1
  IF(i==nsample)in=1
  DO j=1,nphi
    bmesh%r(:,(i-1)*nphi+j)=(/mode_resample(1,i)*COS(phi_grid(j)), &
      mode_resample(1,i)*SIN(phi_grid(j)),mode_resample(2,i)/)
    bnorm((i-1)*nphi+j,1)=mode_resample(3,i)*SIN(nmode*phi_grid(j)) &
      + mode_resample(6,i)*COS(nmode*phi_grid(j))
    bnorm((i-1)*nphi+j,2)=mode_resample(3,i)*COS(nmode*phi_grid(j)) &
      - mode_resample(6,i)*SIN(nmode*phi_grid(j))
    jn=j+1
    IF(j==nphi)jn=1
    bmesh%lc(:,((i-1)*nphi+j-1)*2+1)=(/(i-1)*nphi+j, (in-1)*nphi+j, (i-1)*nphi+jn/)
    bmesh%lc(:,((i-1)*nphi+j-1)*2+2)=(/(in-1)*nphi+j,(in-1)*nphi+jn,(i-1)*nphi+jn/)
  END DO
END DO
END SUBROUTINE setup_mode_mesh
!------------------------------------------------------------------------------
! SUBROUTINE solve_inv
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE solve_inv(self,u,g)
TYPE(tw_type), INTENT(in) :: self
REAL(8), TARGET, INTENT(inout) :: u(:,:),g(:,:)
!---
INTEGER(4) :: i
REAL(8) :: elapsed_time
REAL(8), POINTER :: tmp(:)
CLASS(oft_vector), POINTER :: uloc,gloc
TYPE(oft_native_dense_matrix), TARGET :: fmat
CLASS(oft_solver), POINTER :: linv
LOGICAL :: pm_save
!---Setup matrix
fmat%nr=self%nelems; fmat%nc=self%nelems
fmat%nrg=self%nelems; fmat%ncg=self%nelems
fmat%M=>self%Lmat
!---Setup local vectors
CALL self%Uloc%new(uloc)
CALL self%Uloc%new(gloc)
CALL fmat%assemble(uloc)
!---Setup L^-1 solver
CALL create_cg_solver(linv)
linv%A=>fmat
linv%its=-2
CALL create_diag_pre(linv%pre)
!---Compute surface currents for Sin/Cos phases
DO i=1,2
  CALL uloc%restore_local(u(:,i))
  CALL gloc%restore_local(g(:,i))
  WRITE(*,*)
  WRITE(*,*)'Starting linear solve'
  CALL mytimer%tick
  CALL linv%apply(uloc,gloc)
  elapsed_time=mytimer%tock()
  WRITE(*,*)'  Time = ',elapsed_time
  tmp=>u(:,i)
  CALL uloc%get_local(tmp)
  tmp=>g(:,i)
  CALL gloc%get_local(tmp)
END DO
!---Cleanup
CALL uloc%delete()
CALL gloc%delete()
DEALLOCATE(uloc,gloc)
END SUBROUTINE solve_inv
!------------------------------------------------------------------------------
! SUBROUTINE get_torus_loop
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE get_torus_loop(bmesh,dir,plist,n)
CLASS(oft_bmesh), INTENT(in) :: bmesh
INTEGER(4), INTENT(in) :: dir
INTEGER(4), POINTER, INTENT(out) :: plist(:)
INTEGER(4), INTENT(out) :: n
!---
INTEGER(4) :: i
!---
IF(dir==1)THEN
  n=nphi
  ALLOCATE(plist(n))
  DO i=1,nphi
    plist(i)=i
  END DO
ELSE IF(dir==2)THEN
  n=ntheta
  ALLOCATE(plist(n))
  DO i=1,ntheta
    plist(i)=(i-1)*nphi+1
  END DO
END IF
END SUBROUTINE get_torus_loop
! !------------------------------------------------------------------------------
! !> Needs Docs
! !------------------------------------------------------------------------------
! SUBROUTINE svd_mat(m,n,A)
! INTEGER(4), INTENT(in) :: M,N
! REAL(8), INTENT(inout) :: A(M,N)
! !---
! INTEGER(4) :: LWORK,minsize,info
! REAL(8) :: elapsed_time
! CHARACTER(LEN=1), PARAMETER :: JOBU='A',JOBVT='A'
! INTEGER(4), ALLOCATABLE, DIMENSION(:) :: IWORK
! REAL(8), ALLOCATABLE, DIMENSION(:) :: S,WORK
! REAL(8), ALLOCATABLE, DIMENSION(:,:) :: VT,U
! TYPE(oft_timer) :: stimer
! !---
! IF(oft_env%pm)THEN
!   WRITE(*,*)'Computing SVD'
!   CALL stimer%tick
! END IF
! minsize = MIN(M,N)
! ALLOCATE(S(minsize),U(M,minsize),VT(minsize,N),WORK(1),IWORK(8*minsize))
! LWORK=-1
! CALL dgesdd( JOBU, M, N, A, M, S, U, M, VT, minsize, WORK, LWORK, IWORK, info )
! LWORK=INT(WORK(1),4)
! IF(oft_debug_print(1).AND.oft_env%pm)WRITE(*,*)'  Work size = ',N*N,LWORK
! DEALLOCATE(WORK)
! ALLOCATE(WORK(LWORK))
! CALL dgesdd( JOBU, M, N, A, M, S, U, M, VT, minsize, WORK, LWORK, IWORK, info )
! IF(info/=0)WRITE(*,*)info
! DEALLOCATE(S,U,VT,WORK)
! IF(oft_env%pm)THEN
!   elapsed_time=stimer%tock()
!   WRITE(*,*)'  Time = ',elapsed_time
! END IF
! END SUBROUTINE svd_mat
END PROGRAM thincurr_from_mode