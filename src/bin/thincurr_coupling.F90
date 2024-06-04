!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file thincurr_coupling.F90
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
PROGRAM thincurr_coupling
USE oft_base
USE oft_io, ONLY: hdf5_field_get_sizes
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
IMPLICIT NONE
#include "local.h"
INTEGER(4) :: nsensors = 0
TYPE(tw_type) :: tw_sim,tw_sim2,mode_source
TYPE(tw_sensors) :: sensors
!
INTEGER(4) :: i,n,ierr,io_unit,ndims
integer(i4), allocatable, dimension(:) :: dim_sizes
REAL(8), POINTER, DIMENSION(:) :: vals
REAL(8), ALLOCATABLE :: eig_rval(:),eig_ival(:),eig_vec(:,:),reg_mat(:,:)
CHARACTER(LEN=2) :: eig_tag
TYPE(oft_timer) :: mytimer
CLASS(oft_vector), POINTER :: uio
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_nsets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_ssets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: hole_nsets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: jumper_nsets => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_nsets2 => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: mesh_ssets2 => NULL()
TYPE(oft_1d_int), POINTER, DIMENSION(:) :: hole_nsets2 => NULL()
LOGICAL :: success
!
INTEGER(4) :: jumper_start = 0
CHARACTER(LEN=OFT_PATH_SLEN) :: mesh_file = 'none'
LOGICAL :: plot_run = .FALSE.
NAMELIST/thincurr_coupling_options/mesh_file,plot_run
!---
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit,thincurr_coupling_options, IOSTAT=ierr)
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
CALL smesh%setup_io(1,basepath='Model1/')
IF(oft_debug_print(1))CALL tw_sim%save_debug()
!---------------------------------------------------------------------------
! Load second model
!---------------------------------------------------------------------------
ALLOCATE(oft_trimesh::tw_sim2%mesh)
CALL tw_sim2%mesh%setup(-1,.FALSE.)
CALL hdf5_field_get_sizes(TRIM(mesh_file),"mesh/R",ndims,dim_sizes)
tw_sim2%mesh%np=dim_sizes(2)
DEALLOCATE(dim_sizes)
CALL hdf5_field_get_sizes(TRIM(mesh_file),"mesh/LC",ndims,dim_sizes)
tw_sim2%mesh%nc=dim_sizes(2)
DEALLOCATE(dim_sizes)
ALLOCATE(tw_sim2%mesh%r(3,tw_sim2%mesh%np))
CALL hdf5_read(tw_sim2%mesh%r,TRIM(mesh_file),"mesh/R",success)
ALLOCATE(tw_sim2%mesh%lc(3,tw_sim2%mesh%nc))
CALL hdf5_read(tw_sim2%mesh%lc,TRIM(mesh_file),"mesh/LC",success)
ALLOCATE(tw_sim2%mesh%reg(tw_sim2%mesh%nc))
tw_sim2%mesh%reg=1.d0
! CALL native_hobase(mg_mesh2%smesh)
CALL native_read_nodesets(mesh_nsets2,native_filename=TRIM(mesh_file))
CALL native_read_sidesets(mesh_ssets2,native_filename=TRIM(mesh_file))
IF(ASSOCIATED(mesh_ssets2))THEN
  IF(mesh_ssets2(1)%n>0)THEN
    tw_sim2%nclosures=mesh_ssets2(1)%n
    ALLOCATE(tw_sim2%closures(tw_sim2%nclosures))
    tw_sim2%closures=mesh_ssets2(1)%v
  END IF
END IF
hole_nsets2=>mesh_nsets2
CALL tw_sim2%setup(hole_nsets2)
!---Setup I/0
CALL tw_sim2%mesh%setup_io(1,basepath='Model2/')
IF(oft_debug_print(1))CALL tw_sim2%save_debug()
!
IF(plot_run)THEN
  CALL tw_sim%Uloc%new(uio)
  CALL tw_rst_load(uio,'pThinCurr_coupling.rst','P')
  NULLIFY(vals)
  CALL uio%get_local(vals)
  !---Save plot fields
  CALL tw_save_pfield(tw_sim,vals,'J')
  CALL uio%delete()
  DEALLOCATE(uio)
  !
  CALL tw_sim2%Uloc%new(uio)
  CALL tw_rst_load(uio,'pThinCurr_coupling.rst','E')
  NULLIFY(vals)
  CALL uio%get_local(vals)
  !---Save plot fields
  CALL tw_sim2%mesh%save_vertex_scalar(vals(1:tw_sim2%mesh%np),'E')
  CALL uio%delete()
  DEALLOCATE(uio)
  CALL oft_finalize()
END IF
!---------------------------------------------------------------------------
! Build mutual coupling matrix
!---------------------------------------------------------------------------
CALL tw_compute_LmatDirect(tw_sim,tw_sim%Lmat,col_model=tw_sim2,save_file='Lmat_coupling.save')
!
CALL tw_get_getMat
OPEN(NEWUNIT=io_unit,FILE='Lmat_reg.save',FORM='UNFORMATTED')
WRITE(io_unit)3*tw_sim%mesh%nc,tw_sim%nelems
DO i=1,tw_sim%nelems
  WRITE(io_unit)reg_mat(:,i)
END DO
CLOSE(io_unit)
!---
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
!> Save solution vector for thin-wall model for plotting in VisIt
!------------------------------------------------------------------------------
SUBROUTINE tw_get_getMat()
INTEGER(4) :: i,j,k,jj,pt,ih,ihp,ihc
REAL(8) :: rcurr(3),ftmp(3),gop(3,3),area,norm(3)
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: ptvec,cellvec
ALLOCATE(reg_mat(3*tw_sim%mesh%nc,tw_sim%nelems))
reg_mat=0.d0
!---Avg to cells
ftmp=1.d0/3.d0
DO i=1,tw_sim%mesh%nc
  CALL tw_sim%mesh%jacobian(i,ftmp,gop,area)
  CALL tw_sim%mesh%norm(i,ftmp,norm)
  DO j=1,3
    pt=tw_sim%pmap(tw_sim%mesh%lc(j,i))
    IF(pt==0)CYCLE
    rcurr = cross_product(gop(:,j),norm)
    reg_mat((i-1)*3+1,pt) = reg_mat((i-1)*3+1,pt) + rcurr(1)*SQRT(tw_sim%mesh%ca(i))
    reg_mat((i-1)*3+2,pt) = reg_mat((i-1)*3+2,pt) + rcurr(2)*SQRT(tw_sim%mesh%ca(i))
    reg_mat((i-1)*3+3,pt) = reg_mat((i-1)*3+3,pt) + rcurr(3)*SQRT(tw_sim%mesh%ca(i))
  END DO
END DO
DO ih=1,tw_sim%nholes
DO ihp=1,tw_sim%hmesh(ih)%n
DO ihc=tw_sim%hmesh(ih)%kpc(ihp),tw_sim%hmesh(ih)%kpc(ihp+1)-1
  i=ABS(tw_sim%hmesh(ih)%lpc(ihc))
  DO j=1,3
    IF(tw_sim%mesh%lc(j,i)==tw_sim%hmesh(ih)%lp(ihp))EXIT
  END DO
  CALL tw_sim%mesh%jacobian(i,ftmp,gop,area)
  CALL tw_sim%mesh%norm(i,ftmp,norm)
  rcurr = cross_product(gop(:,j),norm)*SIGN(1,tw_sim%hmesh(ih)%lpc(ihc))
  pt = tw_sim%np_active+ih
  reg_mat((i-1)*3+1,pt) = reg_mat((i-1)*3+1,pt) + rcurr(1)*SQRT(tw_sim%mesh%ca(i))
  reg_mat((i-1)*3+2,pt) = reg_mat((i-1)*3+2,pt) + rcurr(2)*SQRT(tw_sim%mesh%ca(i))
  reg_mat((i-1)*3+3,pt) = reg_mat((i-1)*3+3,pt) + rcurr(3)*SQRT(tw_sim%mesh%ca(i))
END DO
END DO
END DO
END SUBROUTINE tw_get_getMat
END PROGRAM thincurr_coupling
