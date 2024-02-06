!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (OpenFUSIONToolkit)
!---------------------------------------------------------------------------
!> @file thincurr_eig.F90
!
!> Run thin wall eigenvalue calculations using ThinCurr
!!
!! **Option group:** `thincurr_eig_options`
!! |  Option                 |  Description  | Type [dim] |
!! |-------------------------|---------------|------------|
!! |  `mesh_file="none"`     |  Surface mesh filename (Cubit) | str |
!! |  `plot_run=F`           |  Produce plot files from stored restart files | bool |
!! |  `direct=T`             |  Use direct solver | bool |
!! |  `compute_B=F`          |  Compute B-field on cell centers | bool |
!! |  `eigs=5`               |  Number of eigenvalues/vectors to compute | ind |
!!
!! @authors Chris Hansen
!! @date Feb 2022
!! @ingroup doxy_thincurr
!---------------------------------------------------------------------------
PROGRAM thincurr_eig
USE oft_base
USE oft_sort, ONLY: sort_array
USE oft_io, ONLY: hdf5_create_timestep
USE oft_mesh_type, ONLY: smesh
USE oft_mesh_native, ONLY: native_read_nodesets, native_read_sidesets
#ifdef HAVE_NCDF
USE oft_mesh_cubit, ONLY: smesh_cubit_load, cubit_read_nodesets, cubit_read_sidesets
#endif
USE multigrid_build, ONLY: multigrid_construct_surf
!
USE oft_la_base, ONLY: oft_vector, oft_graph
USE oft_lu, ONLY: oft_lusolver, lapack_matinv, lapack_cholesky
USE oft_native_la, ONLY: oft_native_vector, oft_native_matrix, oft_native_dense_matrix, &
  partition_graph
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
#include "local.h"
INTEGER(4) :: nsensors = 0
TYPE(tw_type) :: tw_sim
TYPE(tw_sensors) :: sensors
!
INTEGER(4) :: i,n,ierr,io_unit
REAL(8), ALLOCATABLE :: eig_rval(:),eig_ival(:),eig_vec(:,:)
CHARACTER(LEN=2) :: eig_tag
TYPE(oft_timer) :: mytimer
CLASS(oft_vector), POINTER :: uio
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_nsets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_ssets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: hole_nsets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: jumper_nsets => NULL()
!
INTEGER(4) :: neigs = 5
INTEGER(4) :: jumper_start = -1
LOGICAL :: direct = .TRUE.
LOGICAL :: plot_run = .FALSE.
LOGICAL :: save_L = .FALSE.
LOGICAL :: save_Mcoil = .FALSE.
LOGICAL :: save_Msen = .FALSE.
LOGICAL :: compute_B = .FALSE.
NAMELIST/thincurr_eig_options/direct,plot_run,save_L,save_Mcoil,save_Msen,compute_B,neigs,jumper_start
!---
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit,thincurr_eig_options, IOSTAT=ierr)
CLOSE(io_unit)
if(ierr<0)call oft_abort('No "thincurr_eig_options" found in input file.','thincurr_eig',__FILE__)
if(ierr>0)call oft_abort('Error parsing "thincurr_eig_options" in input file.','thincurr_eig',__FILE__)
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
  CALL oft_abort("Unsupported mesh type","thincurr_eig",__FILE__)
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
! Eigenvalue run
!---------------------------------------------------------------------------
!---Compute inductances
WRITE(*,*)
IF(plot_run)THEN
  !---Load sensors
  CALL tw_load_sensors(tw_sim,sensors,jumper_nsets)
  tw_sim%n_icoils=0
  CALL tw_compute_mutuals(tw_sim,sensors%nfloops,sensors%floops)
  CALL plot_eig(tw_sim,sensors%nfloops,sensors%floops)
ELSE
  tw_sim%n_icoils=0
  sensors%nfloops=0
  ALLOCATE(sensors%floops(sensors%nfloops))
  IF(tw_sim%n_vcoils>0)THEN
    IF(save_Mcoil)THEN
      CALL tw_compute_Ael2dr(tw_sim,'Mcoil.save')
    ELSE
      CALL tw_compute_Ael2dr(tw_sim)
    END IF
  END IF
  IF(save_Msen)THEN
    CALL tw_compute_mutuals(tw_sim,sensors%nfloops,sensors%floops,'Msen.save')
  ELSE
    CALL tw_compute_mutuals(tw_sim,sensors%nfloops,sensors%floops)
  END IF
  IF(save_L)THEN
    CALL tw_compute_LmatDirect(tw_sim,tw_sim,tw_sim%Lmat,'Lmat.save')
  ELSE
    CALL tw_compute_LmatDirect(tw_sim,tw_sim,tw_sim%Lmat)
  END IF
  !---Setup resistivity matrix
  CALL tw_compute_Rmat(tw_sim,.TRUE.)
  !---Run eigenvalue analysis
  ALLOCATE(eig_rval(neigs),eig_ival(neigs),eig_vec(tw_sim%nelems,neigs))
  IF(direct)THEN
    CALL decay_eigenmodes(tw_sim,neigs,eig_rval,eig_ival,eig_vec)
  ELSE
    CALL run_eig(tw_sim,neigs,eig_rval,eig_ival,eig_vec)
  END IF
  !---Save plot files
  CALL tw_sim%Uloc%new(uio)
  OPEN(NEWUNIT=io_unit, FILE="thincurr_eigs.dat")
  DO i=1,neigs
    WRITE(io_unit,*)eig_rval(i),eig_ival(i)
    WRITE(eig_tag,'(I2.2)')i
    CALL uio%restore_local(eig_vec(:,i))
    CALL tw_rst_save(tw_sim,uio,'pThinCurr_eig.rst','U_'//eig_tag,append=(i/=1))
  END DO
  CLOSE(io_unit)
  CALL uio%delete()
  DEALLOCATE(uio,eig_rval,eig_ival,eig_vec)
END IF
!---
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE decay_eigenmodes
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE decay_eigenmodes(self,neigs,eig_rval,eig_ival,eig_vec)
TYPE(tw_type), INTENT(in) :: self
INTEGER(4), INTENT(in) :: neigs
REAL(8), intent(out) :: eig_rval(:),eig_ival(:),eig_vec(:,:)
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
  j = sort_tmp(N+1-i-self%nclosures)
  eig_rval(i)=WR(j)
  eig_ival(i)=WI(j)
  eig_vec(:,i)=REAL(VR(:,j),8)
END DO
DEALLOCATE(WR,WI,VR,sort_tmp)
! !---Handle closures
! DO i=1,self%nclosures
!   eig_vec(self%closures(i),1:neigs)=0.d0
! END DO
WRITE(*,*)'Eigenvalues'
WRITE(*,*)'  Real: ',eig_rval(1:MIN(5,neigs))
WRITE(*,*)'  Imag: ',eig_ival(1:MIN(5,neigs))
WRITE(*,*)
DEBUG_STACK_POP
END SUBROUTINE decay_eigenmodes
!------------------------------------------------------------------------------
! SUBROUTINE run_eig
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE run_eig(self,neigs,eig_rval,eig_ival,eig_vec)
TYPE(tw_type), INTENT(in) :: self
INTEGER(4), INTENT(in) :: neigs
REAL(8), INTENT(out) :: eig_rval(:),eig_ival(:),eig_vec(:,:)
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
TYPE(oft_iram_eigsolver) :: arsolver
TYPE(oft_timer) :: loctimer
DEBUG_STACK_PUSH
WRITE(*,*)
WRITE(*,*)'Starting eigenvalue solve'
!---Setup matrix
fmat%nr=self%nelems; fmat%nc=self%nelems
fmat%nrg=self%nelems; fmat%ncg=self%nelems
fmat%M => self%Lmat
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
arsolver%A=>fmat
arsolver%M=>self%Rmat
arsolver%tol=1.E-5_r8
arsolver%Minv=>linv
arsolver%nev=neigs+self%nclosures
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
    j = sort_tmp(arsolver%nev-self%nclosures+1-i)
    eig_rval(i)=arsolver%eig_val(1,j)
    eig_ival(i)=arsolver%eig_val(2,j)
    eig_vec(:,i)=arsolver%eig_vec(:,j)
  END DO
  DEALLOCATE(sort_tmp)
  ! !---Handle closures
  ! DO i=1,self%nclosures
  !   eig_vec(self%closures(i),1:neigs)=0.d0
  ! END DO
  WRITE(*,*)'Eigenvalues'
  WRITE(*,*)'  Real: ',eig_rval(1:MIN(5,neigs))
  WRITE(*,*)'  Imag: ',eig_ival(1:MIN(5,neigs))
  WRITE(*,*)
END IF
!---Cleanup
CALL uloc%delete()
CALL arsolver%delete()
DEALLOCATE(uloc)
#else
IF(oft_env%test_run)THEN
  WRITE(*,*)'SKIP TEST'
  CALL oft_finalize
ELSE
  CALL oft_abort("Iterative eigenvalue solve requires ARPACK", "run_eig", __FILE__)
END IF
#endif
DEBUG_STACK_POP
END SUBROUTINE run_eig
!------------------------------------------------------------------------------
! SUBROUTINE plot_eig
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE plot_eig(self,nsensors,sensors)
TYPE(tw_type), INTENT(inout) :: self
INTEGER(4), INTENT(in) :: nsensors
TYPE(floop_sensor), POINTER, INTENT(in) :: sensors(:)
INTEGER(4) :: i,j,k,jj,ntimes,ncoils,itime,io_unit
REAL(8) :: uu,t,tmp,area
REAL(8), ALLOCATABLE, DIMENSION(:) :: coil_vec
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: cc_vals,senout
REAL(8), POINTER, DIMENSION(:) :: vals
CLASS(oft_vector), POINTER :: uio
LOGICAL :: exists
CHARACTER(LEN=4) :: pltnum
CHARACTER(LEN=2) :: eig_tag
DEBUG_STACK_PUSH
WRITE(*,*)'Post-processing eigenmode run'
CALL self%Uloc%new(uio)
ALLOCATE(vals(self%nelems))
IF(nsensors>0)ALLOCATE(senout(nsensors,neigs))
IF(compute_B)THEN
  IF(.NOT.ALLOCATED(cc_vals))ALLOCATE(cc_vals(3,self%mesh%nc))
  CALL tw_compute_Bops(self)
END IF
DO i=1,neigs
  !---Load solution from file
  WRITE(eig_tag,'(I2.2)')i
  CALL tw_rst_load(uio,'pThinCurr_eig'//'.rst','U_'//eig_tag)
  CALL uio%get_local(vals)
  !---Save plot fields
  CALL tw_save_pfield(self,vals,'J_'//eig_tag)
  IF(compute_B)THEN
    !$omp parallel do private(j,jj,tmp)
    DO k=1,smesh%nc
      DO jj=1,3
        tmp=0.d0
        !$omp simd reduction(+:tmp)
        DO j=1,self%nelems
          tmp=tmp+vals(j)*self%Bel(j,k,jj)
        END DO
        cc_vals(jj,k)=tmp
      END DO
    END DO
    CALL self%mesh%save_cell_vector(cc_vals,'B_'//eig_tag)
  END IF
  !---Save sensor signals
  IF(nsensors>0)THEN
    CALL dgemv('N',nsensors,self%nelems,1.d0,self%Ael2sen,nsensors,vals,1,0.d0,senout(:,i),1)
  END IF
END DO
!---Save sensor signals
IF(nsensors>0)THEN
  OPEN(NEWUNIT=io_unit,FILE='thincurr_eig-floops.dat')
  DO i=1,nsensors
    WRITE(io_unit,*)sensors(i)%name,senout(i,:)
  END DO
  CLOSE(io_unit)
  DEALLOCATE(senout)
END IF
!---Cleanup
CALL uio%delete
DEALLOCATE(uio,vals)
IF(ALLOCATED(cc_vals))DEALLOCATE(cc_vals)
DEBUG_STACK_POP
END SUBROUTINE plot_eig
END PROGRAM thincurr_eig
