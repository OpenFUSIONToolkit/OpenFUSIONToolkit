!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
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
USE oft_io, ONLY: hdf5_create_timestep, oft_bin_file, hdf5_create_file, hdf5_write, &
  hdf5_create_group, hdf5_add_string_attribute
USE oft_mesh_type, ONLY: smesh
USE oft_mesh_native, ONLY: native_read_nodesets, native_read_sidesets
#ifdef HAVE_NCDF
USE oft_mesh_cubit, ONLY: cubit_read_nodesets, cubit_read_sidesets
#endif
USE multigrid_build, ONLY: multigrid_construct_surf
!
USE oft_la_base, ONLY: oft_vector
USE mhd_utils, ONLY: mu0
USE thin_wall
USE thin_wall_hodlr
USE thin_wall_solvers, ONLY: lr_eigenmodes_direct, lr_eigenmodes_arpack
IMPLICIT NONE
#include "local.h"
INTEGER(4) :: nsensors = 0
TYPE(tw_type), TARGET :: tw_sim
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
TYPE(oft_tw_hodlr_op), TARGET :: tw_hodlr
!
INTEGER(4) :: neigs = 5
INTEGER(4) :: jumper_start = -1
LOGICAL :: direct = .TRUE.
LOGICAL :: plot_run = .FALSE.
LOGICAL :: save_L = .FALSE.
LOGICAL :: save_Mcoil = .FALSE.
LOGICAL :: save_Msen = .FALSE.
LOGICAL :: compute_B = .FALSE.
LOGICAL :: reduce_model = .FALSE.
NAMELIST/thincurr_eig_options/direct,plot_run,save_L,save_Mcoil,save_Msen, &
  compute_B,neigs,jumper_start,reduce_model
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
  CALL tw_load_sensors('floops.loc',tw_sim,sensors,jumper_nsets)
  tw_sim%n_icoils=0
  CALL tw_compute_mutuals(tw_sim,sensors%nfloops,sensors%floops)
  CALL plot_eig(tw_sim,sensors%nfloops,sensors%floops)
ELSE
  IF(reduce_model)THEN
    CALL tw_load_sensors(tw_sim,sensors,jumper_nsets)
  ELSE
    tw_sim%n_icoils=0
    sensors%nfloops=0
    ALLOCATE(sensors%floops(sensors%nfloops))
  END IF
  IF((tw_sim%n_vcoils>0).OR.(tw_sim%n_icoils>0))THEN
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
  !---
  tw_hodlr%tw_obj=>tw_sim
  CALL tw_hodlr%setup(.FALSE.)
  IF(tw_hodlr%L_svd_tol>0.d0)THEN
    IF(direct)CALL oft_abort('HODLR compression does not support "direct=T"','thincurr_eig',__FILE__)
    IF(save_L)CALL oft_abort('HODLR compression does not support "save_L=T"','thincurr_eig',__FILE__)
    CALL tw_hodlr%compute_L()
  ELSE
    IF(save_L)THEN
      CALL tw_compute_LmatDirect(tw_sim,tw_sim%Lmat,save_file='Lmat.save')
    ELSE
      CALL tw_compute_LmatDirect(tw_sim,tw_sim%Lmat)
    END IF
  END IF
  !---Setup resistivity matrix
  CALL tw_compute_Rmat(tw_sim,.TRUE.)
  !---Run eigenvalue analysis
  ALLOCATE(eig_rval(neigs),eig_ival(neigs),eig_vec(tw_sim%nelems,neigs))
  WRITE(*,*)
  IF(direct)THEN
    CALL lr_eigenmodes_direct(tw_sim,neigs,eig_rval,eig_vec,eig_ival)
  ELSE
#if !defined(HAVE_ARPACK)
    IF(oft_env%test_run)THEN
      WRITE(*,*)'SKIP TEST'
      CALL oft_finalize
    END IF
#endif
    IF(tw_hodlr%L_svd_tol>0.d0)THEN
      CALL lr_eigenmodes_arpack(tw_sim,neigs,eig_rval,eig_vec,eig_ival=eig_ival,hodlr_op=tw_hodlr)
    ELSE
      CALL lr_eigenmodes_arpack(tw_sim,neigs,eig_rval,eig_vec,eig_ival)
    END IF
  END IF
  WRITE(*,*)
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
  DEALLOCATE(uio)
  !
  IF(reduce_model)CALL model_reduction_eig(tw_sim,neigs,eig_vec)
  !
  DEALLOCATE(eig_rval,eig_ival,eig_vec)
END IF
!---
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
! SUBROUTINE plot_eig
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE plot_eig(self,nsensors,sensors)
TYPE(tw_type), TARGET, INTENT(inout) :: self
INTEGER(4), INTENT(in) :: nsensors
TYPE(floop_sensor), POINTER, INTENT(in) :: sensors(:)
INTEGER(4) :: i,j,k,jj,ntimes,ncoils,itime,io_unit
REAL(8) :: uu,t,tmp,area
REAL(8), ALLOCATABLE, DIMENSION(:) :: coil_vec
REAL(8), ALLOCATABLE, TARGET, DIMENSION(:,:) :: cc_vals,senout
REAL(8), POINTER, DIMENSION(:) :: vals,vtmp
CLASS(oft_vector), POINTER :: uio,Bx,By,Bz
TYPE(oft_bin_file) :: floop_hist
LOGICAL :: exists
CHARACTER(LEN=4) :: pltnum
CHARACTER(LEN=2) :: eig_tag
DEBUG_STACK_PUSH
WRITE(*,*)'Post-processing eigenmode run'
CALL self%Uloc%new(uio)
ALLOCATE(vals(self%nelems))
IF(nsensors>0)ALLOCATE(senout(nsensors,neigs))
IF(compute_B)THEN
  IF(.NOT.ALLOCATED(cc_vals))ALLOCATE(cc_vals(3,self%mesh%np))
  tw_hodlr%tw_obj=>self
  CALL tw_hodlr%setup(.FALSE.)
  IF(tw_hodlr%B_svd_tol>0.d0)THEN
    CALL tw_hodlr%compute_B()
    CALL self%Uloc_pts%new(Bx)
    CALL self%Uloc_pts%new(By)
    CALL self%Uloc_pts%new(Bz)
  ELSE
    CALL tw_compute_Bops(self)
  END IF
END IF
DO i=1,neigs
  !---Load solution from file
  WRITE(eig_tag,'(I2.2)')i
  CALL tw_rst_load(uio,'pThinCurr_eig'//'.rst','U_'//eig_tag)
  CALL uio%get_local(vals)
  !---Save plot fields
  CALL tw_save_pfield(self,vals,'J_'//eig_tag)
  IF(compute_B)THEN
    IF(tw_hodlr%B_svd_tol>0.d0)THEN
      CALL tw_hodlr%apply_bop(uio,Bx,By,Bz)
      NULLIFY(vtmp)
      CALL Bx%get_local(vtmp)
      cc_vals(1,:)=vtmp
      CALL By%get_local(vtmp)
      cc_vals(2,:)=vtmp
      CALL Bz%get_local(vtmp)
      cc_vals(3,:)=vtmp
    ELSE
      !$omp parallel do private(j,jj,tmp)
      DO k=1,smesh%np
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
    CALL self%mesh%save_vertex_vector(cc_vals,'B_v_'//eig_tag)
  END IF
  !---Save sensor signals
  IF(nsensors>0)THEN
    CALL dgemv('N',nsensors,self%nelems,1.d0,self%Ael2sen,nsensors,vals,1,0.d0,senout(:,i),1)
  END IF
END DO
!---Save sensor signals
IF(nsensors>0)THEN
  !---Setup history file
  IF(oft_env%head_proc)THEN
    floop_hist%filedesc = 'ThinCurr eigenvalue flux loop coupling'
    CALL floop_hist%setup('thincurr_eig-floops.dat')
    DO i=1,nsensors
      CALL floop_hist%add_field(sensors(i)%name, 'r8')
    END DO
    CALL floop_hist%write_header
    CALL floop_hist%open
    DO i=1,neigs
      CALL floop_hist%write(data_r8=senout(:,i))
    END DO
    CALL floop_hist%close
  END IF
  DEALLOCATE(senout)
END IF
!---Cleanup
CALL uio%delete
DEALLOCATE(uio,vals)
IF(ALLOCATED(cc_vals))DEALLOCATE(cc_vals)
DEBUG_STACK_POP
END SUBROUTINE plot_eig
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE model_reduction_eig(self,neigs,eig_vec)
TYPE(tw_type), INTENT(in) :: self
INTEGER(4), INTENT(in) :: neigs
REAL(8), INTENT(out) :: eig_vec(:,:)
!---
INTEGER(4) :: i,j
REAL(8), POINTER :: vals(:)
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: Mattmp,Mat_red
CHARACTER(LEN=4) :: sub_coil_id
CLASS(oft_vector), POINTER :: atmp,btmp
character(LEN=80), dimension(1) :: description
CALL hdf5_create_file('tCurr_reduced.h5')
ALLOCATE(Mat_red(neigs,neigs),Mattmp(self%nelems,neigs))
!---Reduce L matrix
CALL dgemm('N','N',self%nelems,neigs,self%nelems,1.d0, &
  self%Lmat,self%nelems,eig_vec,self%nelems,0.d0,Mattmp,self%nelems)
CALL dgemm('T','N',neigs,neigs,self%nelems,1.d0,eig_vec, &
  self%nelems,Mattmp,self%nelems,0.d0,Mat_red,neigs)
CALL hdf5_write(Mat_red,'tCurr_reduced.h5','L')
description=["Self-inductance matrix for reduced model"]
CALL hdf5_add_string_attribute('tCurr_reduced.h5','L','description',description)
!---Reduce R matrix
NULLIFY(vals)
CALL self%Uloc%new(atmp)
CALL self%Uloc%new(btmp)
CALL atmp%get_local(vals)
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
CALL hdf5_write(Mat_red,'tCurr_reduced.h5','R')
description=["Resistance matrix for reduced model"]
CALL hdf5_add_string_attribute('tCurr_reduced.h5','R','description',description)
DEALLOCATE(Mat_red,Mattmp)
!---Reduce sensor coupling matrix
IF(sensors%nfloops>0)THEN
  ! Save sensor names and indexes
  CALL hdf5_create_group('tCurr_reduced.h5','SENSORS')
  DO i=1,sensors%nfloops
    CALL hdf5_create_group('tCurr_reduced.h5', &
      'SENSORS/'//TRIM(sensors%floops(i)%name))
    CALL hdf5_write(i,'tCurr_reduced.h5', &
      'SENSORS/'//TRIM(sensors%floops(i)%name)//'/index')
    CALL hdf5_write(sensors%floops(i)%scale_fac,'tCurr_reduced.h5', &
      'SENSORS/'//TRIM(sensors%floops(i)%name)//'/scale')
    CALL hdf5_write(sensors%floops(i)%r,'tCurr_reduced.h5', &
      'SENSORS/'//TRIM(sensors%floops(i)%name)//'/pts')
  END DO
  ! Save mutual coupling matrix
  ALLOCATE(Mat_red(sensors%nfloops,neigs))
  CALL dgemm('N','N',sensors%nfloops,neigs,self%nelems,1.d0, &
    self%Ael2sen,sensors%nfloops,eig_vec,self%nelems,0.d0,Mat_red,sensors%nfloops)
  CALL hdf5_write(Mat_red,'tCurr_reduced.h5','Ms')
  description=["Model to sensor mutual inductance matrix"]
  CALL hdf5_add_string_attribute('tCurr_reduced.h5','Ms','description',description)
  DEALLOCATE(Mat_red)
END IF
!---Reduce coil coupling matrix
IF(self%n_icoils>0)THEN
  ! Save coil names and indexes
  CALL hdf5_create_group('tCurr_reduced.h5','COILS')
  DO i=1,self%n_icoils
    CALL hdf5_create_group('tCurr_reduced.h5', &
      'COILS/'//TRIM(self%icoils(i)%name))
    CALL hdf5_write(i,'tCurr_reduced.h5', &
      'COILS/'//TRIM(self%icoils(i)%name)//'/index')
    CALL hdf5_create_group('tCurr_reduced.h5', &
      'COILS/'//TRIM(self%icoils(i)%name)//'/SUB_COILS')
    DO j=1,self%icoils(i)%ncoils
      WRITE(sub_coil_id,'(I4.4)')j
      CALL hdf5_create_group('tCurr_reduced.h5', &
        'COILS/'//TRIM(self%icoils(i)%name)//'/SUB_COILS/'//sub_coil_id)
      CALL hdf5_write(self%icoils(i)%scales(j),'tCurr_reduced.h5', &
        'COILS/'//TRIM(self%icoils(i)%name)//'/SUB_COILS/'//sub_coil_id//'/scale')
      CALL hdf5_write(self%icoils(i)%coils(j)%pts,'tCurr_reduced.h5', &
        'COILS/'//TRIM(self%icoils(i)%name)//'/SUB_COILS/'//sub_coil_id//'/pts')
    END DO
  END DO
  ! Save mutual coupling matrix
  ALLOCATE(Mat_red(neigs,self%n_icoils))
  CALL dgemm('T','N',neigs,self%n_icoils,self%nelems,1.d0, &
    eig_vec,self%nelems,self%Ael2dr,self%nelems,0.d0,Mat_red,neigs)
  CALL hdf5_write(Mat_red,'tCurr_reduced.h5','Mc')
  description=["Model to coil mutual inductance matrix"]
  CALL hdf5_add_string_attribute('tCurr_reduced.h5','Mc','description',description)
  DEALLOCATE(Mat_red)
  IF(sensors%nfloops>0)THEN
    CALL hdf5_write(self%Adr2sen,'tCurr_reduced.h5','Msc')
    description=["Coil to sensor mutual inductance matrix"]
    CALL hdf5_add_string_attribute('tCurr_reduced.h5','Msc','description',description)
  END IF
END IF
END SUBROUTINE model_reduction_eig
END PROGRAM thincurr_eig
