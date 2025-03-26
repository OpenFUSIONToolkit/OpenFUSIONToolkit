!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file multigrid_build.F90
!
!> Multigrid initialization using nested meshes.
!!
!! @authors George Marklin and Chris Hansen
!! @date December 2009
!! @ingroup doxy_oft_grid
!--------------------------------------------------------------------------------
module multigrid_build
use oft_base
use oft_mesh_type, only: oft_mesh, oft_bmesh, cell_is_curved
use oft_mesh_local, only: mesh_local_init, mesh_volumes, bmesh_local_init, bmesh_areas
use oft_mesh_global, only: mesh_global_init, mesh_global_link, mesh_global_decomp, &
  bmesh_global_link, smesh_global_init, bmesh_global_init
use oft_mesh_global_util, only: mesh_global_orient, mesh_global_stats, &
  mesh_global_boundary, mesh_global_save, mesh_global_igrnd, mesh_global_periodic, &
  bmesh_global_boundary, bmesh_global_periodic, bmesh_global_orient, bmesh_global_stats, &
  mesh_global_set_curved
use oft_stitching, only: oft_global_stitch
USE oft_tetmesh_type, only: oft_tetmesh
use multigrid, only: multigrid_mesh, multigrid_refine, multigrid_hybrid_base, &
  multigrid_hybrid_bmesh, multigrid_brefine, hexmesh_mg_globals, tetmesh_mg_globals, &
  multigrid_level, trimesh_mg_globals, quadmesh_mg_globals, multigrid_reffix_ho, &
  multigrid_reffix_ho_surf
use oft_mesh_native, only: native_load_vmesh, native_load_smesh, mesh_native_id, &
  native_hobase, native_set_periodic, native_bset_periodic, native_finalize_setup
use oft_mesh_t3d, only: mesh_t3d_load, mesh_t3d_cadsync, mesh_t3d_cadlink, &
  mesh_t3d_add_quad, mesh_t3d_reffix, mesh_t3d_add_quad, &
  mesh_t3d_set_periodic, smesh_t3d_load, mesh_t3d_id
#ifdef HAVE_NCDF
use oft_mesh_cubit, only: mesh_cubit_load, mesh_cubit_reffix, mesh_cubit_cadlink, &
  mesh_cubit_add_quad, mesh_cubit_hobase, mesh_cubit_set_periodic, &
  mesh_cubit_id, smesh_cubit_load, cubit_finalize_setup
#else
use oft_mesh_cubit, only: mesh_cubit_id
#endif
use oft_mesh_gmsh, only: mesh_gmsh_load, mesh_gmsh_reffix, mesh_gmsh_cadlink, &
  mesh_gmsh_add_quad, mesh_gmsh_id, gmsh_finalize_setup
use oft_mesh_sphere, only: mesh_sphere_load, mesh_sphere_reffix, mesh_sphere_cadlink, &
  mesh_sphere_add_quad, smesh_circle_load, smesh_circle_cadlink, smesh_circle_reffix, &
  smesh_circle_add_quad, mesh_sphere_id
use oft_mesh_cube, only: mesh_cube_load, mesh_cube_cadlink, &
  mesh_cube_set_periodic, mesh_cube_id, smesh_square_load, smesh_square_cadlink, &
  smesh_square_set_periodic
!---End include modules
implicit none
#define MOD_NAME "multigrid_build"
#include "local.h"
INTEGER(i4), PRIVATE :: part_meth = 1
REAL(r8), PRIVATE :: jac_ratio_tol = 10.d0
REAL(r8), PRIVATE :: cad_feature_angle = pi/4.d0
LOGICAL, PRIVATE :: fix_boundary = .TRUE.
contains
!------------------------------------------------------------------------------
!> Load in mesh and CAD information.
!! - Read in mesh points and cells
!! - Read in CAD information for refinement
!------------------------------------------------------------------------------
subroutine multigrid_load(mg_mesh,cad_type)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4), intent(in) :: cad_type !< Mesh type to load
DEBUG_STACK_PUSH
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(2A)')oft_indent,'**** Loading OFT mesh'
END IF
CALL oft_increase_indent
!---Select mesh type and load
select case(cad_type)
  case(mesh_native_id) ! Native Mesh
    CALL native_load_vmesh(mg_mesh)
    CALL mesh_global_init(mg_mesh%mesh)
    CALL native_hobase(mg_mesh%mesh)
    CALL native_set_periodic(mg_mesh%mesh)
  case(mesh_t3d_id) ! T3D Mesh
    CALL mesh_t3d_load(mg_mesh)
    CALL mesh_global_init(mg_mesh%mesh)
    CALL mesh_t3d_cadsync(mg_mesh)
    CALL mesh_t3d_cadlink(mg_mesh%mesh)
    CALL mesh_t3d_set_periodic(mg_mesh%mesh)
  case(mesh_cubit_id) ! Exodus Mesh
#ifdef HAVE_NCDF
    CALL mesh_cubit_load(mg_mesh)
    CALL mesh_global_init(mg_mesh%mesh)
    CALL mesh_cubit_cadlink(mg_mesh%mesh)
    CALL mesh_cubit_hobase(mg_mesh%mesh)
    CALL mesh_cubit_set_periodic(mg_mesh%mesh)
#else
    CALL oft_abort('CUBIT interface requires NETCDF','multigrid_load',__FILE__)
#endif
  case(mesh_gmsh_id) ! GMSH Mesh
    CALL mesh_gmsh_load(mg_mesh)
    CALL mesh_global_init(mg_mesh%mesh)
    CALL mesh_gmsh_cadlink(mg_mesh%mesh)
  case(mesh_sphere_id) ! Sphere Test Mesh
    CALL mesh_sphere_load(mg_mesh)
    CALL mesh_global_init(mg_mesh%mesh)
    CALL mesh_sphere_cadlink(mg_mesh%mesh)
  case(mesh_cube_id) ! Cube Test Mesh
    CALL mesh_cube_load(mg_mesh)
    CALL mesh_global_init(mg_mesh%mesh)
    CALL mesh_cube_cadlink(mg_mesh%mesh)
    CALL mesh_cube_set_periodic(mg_mesh%mesh)
  case default ! Invalid Mesh
    CALL oft_abort('Invalid mesh type.','multigrid_load',__FILE__)
end select
!---Mesh load complete
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine multigrid_load
!------------------------------------------------------------------------------
!> Adjust boundary points to CAD boundary following refinement.
!------------------------------------------------------------------------------
subroutine multigrid_reffix(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
INTEGER(i4) :: i,j
class(oft_mesh), pointer :: mesh
DEBUG_STACK_PUSH
IF(.NOT.fix_boundary)THEN
  IF(oft_env%head_proc)WRITE(*,'(2A)')oft_indent,'Skipping boundary corrections'
  DEBUG_STACK_POP
  RETURN
END IF
mesh=>mg_mesh%mesh
!---Always refine using mesh mapping first
CALL multigrid_reffix_ho(mg_mesh)
!---Select mesh type and adjust boundary
select case(mesh%cad_type)
  case(mesh_native_id)
    ! Do nothing
  case(mesh_t3d_id)
    call mesh_t3d_reffix(mg_mesh)
  case(mesh_cubit_id)
#ifdef HAVE_NCDF
    call mesh_cubit_reffix(mg_mesh)
#endif
  case(mesh_gmsh_id)
    call mesh_gmsh_reffix(mg_mesh)
  case(mesh_sphere_id)
    call mesh_sphere_reffix(mg_mesh)
  case(mesh_cube_id)
    ! Do nothing
  case default
    call oft_abort('Invalid mesh type.','multigrid_reffix',__FILE__)
end select
!---Propogate updates to boundary mesh
!!$omp parallel do private(j)
do i=1,mesh%bmesh%np
  j=INT(ABS(mesh%bmesh%parent%lp(i)),4)
  mesh%bmesh%r(:,i)=mesh%r(:,j)
enddo
!---Boundary adjustment complete
DEBUG_STACK_POP
end subroutine multigrid_reffix
!------------------------------------------------------------------------------
!> Add node points for quadratic elements
!------------------------------------------------------------------------------
subroutine multigrid_add_quad(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
class(oft_mesh), pointer :: mesh
INTEGER(i4) :: i,j
REAL(r8) :: pttmp(3)
DEBUG_STACK_PUSH
mesh=>mg_mesh%mesh
!---Initialize high order points with straight edges
IF(.NOT.ASSOCIATED(mesh%ho_info%r))THEN
  CALL mesh%set_order(2)
  DO i=1,mesh%ne
    mesh%ho_info%lep(1,i)=i
    mesh%ho_info%r(:,i)=(mesh%r(:,mesh%le(1,i))+mesh%r(:,mesh%le(2,i)))/2.d0
  END DO
  IF(mesh%ho_info%nfp==1)THEN
    do i=1,mesh%nf
      mesh%ho_info%lfp(1,i)=i+mesh%ne
      pttmp=0.d0
      DO j=1,mesh%face_np
        pttmp=pttmp+mesh%r(:,mesh%lf(j,i))/REAL(mesh%face_np)
      END DO
      mesh%ho_info%r(:,i+mesh%ne)=pttmp
    end do
  END IF
  IF(mesh%ho_info%ncp==1)THEN
    do i=1,mesh%nc
      mesh%ho_info%lcp(1,i)=i+mesh%ne+mesh%nf
      pttmp=0.d0
      DO j=1,mesh%cell_np
        pttmp=pttmp+mesh%r(:,mesh%lc(j,i))/REAL(mesh%cell_np)
      END DO
      mesh%ho_info%r(:,i+mesh%ne+mesh%nf)=pttmp
    end do
  END IF
ELSE
  CALL mesh%set_order(2)
END IF
CALL mesh_global_set_curved(mesh,1)
!---Select mesh type and adjust boundary
select case(mesh%cad_type)
  case(mesh_native_id)
    ! Do nothing CALL mesh_cube_add_quad
  case(mesh_t3d_id)
    call mesh_t3d_add_quad(mg_mesh)
  case(mesh_cubit_id)
#ifdef HAVE_NCDF
    call mesh_cubit_add_quad(mg_mesh)
#endif
  case(mesh_gmsh_id)
    call mesh_gmsh_add_quad(mg_mesh)
  case(mesh_sphere_id)
    call mesh_sphere_add_quad(mg_mesh%mesh)
  case(mesh_cube_id)
    ! Do nothing CALL mesh_cube_add_quad
  case default
    call oft_abort('Invalid mesh type.','multigrid_add_quad',__FILE__)
end select
CALL multigrid_check_ho(mg_mesh%mesh)
!---Propogate to boundary mesh
CALL mesh%bmesh%set_order(2)
CALL mesh_global_set_curved(mesh%bmesh,1)
!---Edge centers
!!$omp parallel do private(j)
do i=1,mesh%bmesh%ne
  mesh%bmesh%ho_info%lep(1,i)=i
  j=INT(ABS(mesh%bmesh%parent%le(i)),4)
  mesh%bmesh%ho_info%r(:,i)=mesh%ho_info%r(:,mesh%ho_info%lep(1,j))
enddo
IF(mesh%ho_info%nfp==1)THEN
  !---Face centers
  !!$omp parallel do private(j)
  do i=1,mesh%bmesh%nc
    j=INT(ABS(mesh%bmesh%parent%lf(i)),4)
    mesh%bmesh%ho_info%lcp(1,i)=i+mesh%bmesh%ne
    mesh%bmesh%ho_info%r(:,i+mesh%bmesh%ne)=mesh%ho_info%r(:,mesh%ho_info%lfp(1,j))
  enddo
END IF
!---Boundary adjustment complete
DEBUG_STACK_POP
end subroutine multigrid_add_quad
!------------------------------------------------------------------------------
!> Check high order geometric representation
!!
!! Validates the high order geometric representation and disables curvature
!! locally when negative grid jacobians are detected. Cells where ill conditioned
!! mappings are found are straightened back to their original linear spatial mapping.
!------------------------------------------------------------------------------
subroutine multigrid_check_ho(mesh)
class(oft_mesh), intent(inout) :: mesh
INTEGER(i4) :: i,j,k,cell,nbad,ed,etmp(2)
INTEGER(i4), PARAMETER :: nx = 5
REAL(r8), ALLOCATABLE, DIMENSION(:) :: eflag,fflag
REAL(r8) :: u,v,w,f(4),goptmp(3,4),vol,vol_min,vol_max
CHARACTER(LEN=70) :: error_str
DEBUG_STACK_PUSH
!---Check cells for negative jacobians
IF(oft_debug_print(2))WRITE(*,'(A)')'Checking cell curvature'
ALLOCATE(eflag(mesh%ne),fflag(mesh%nf))
eflag=0.d0
fflag=0.d0
nbad=0
cell_loop: DO cell=1,mesh%nc
  IF(.NOT.cell_is_curved(mesh,cell))CYCLE
  vol_min=1.d99
  vol_max=-1.d99
  DO i=1,nx
    u=(i-1)/REAL(nx,8)
    DO j=1,nx
      v=(j-1)*(1.d0-u)/REAL(nx,8)
      IF(v<0.d0)CYCLE
      DO k=1,nx
        w=(k-1)*(1.d0-u-v)/REAL(nx,8)
        IF((w<0.d0).OR.(1.d0-u-v-w<0.d0))CYCLE
        f=(/u,v,w,1.d0-u-v-w/)
        CALL mesh%jacobian(cell,f,goptmp,vol)
        vol_min = MIN(vol_min,vol)
        vol_max = MAX(vol_max,vol)
        IF(vol<=0.d0)THEN
          eflag(ABS(mesh%lce(:,cell)))=2.d0
          fflag(ABS(mesh%lcf(:,cell)))=2.d0
          nbad=nbad+1
          IF(oft_debug_print(2))WRITE(*,'(A,I4,I6,ES11.3)')'  Negative jacobian',oft_env%rank,cell,vol
          CYCLE cell_loop
        END IF
      END DO
    END DO
  END DO
  IF(vol_max/vol_min>jac_ratio_tol)THEN
    eflag(ABS(mesh%lce(:,cell)))=2.d0
    fflag(ABS(mesh%lcf(:,cell)))=2.d0
    nbad=nbad+1
    IF(oft_debug_print(2))WRITE(*,'(A,I4,I6,ES11.3)')'  Poorly shaped cell',oft_env%rank,cell,vol_max/vol_min
  END IF
END DO cell_loop
!---Synchornize flags
nbad=oft_mpi_sum(nbad)
IF(oft_env%head_proc.AND.nbad>0)THEN
  WRITE(error_str,'(A,I4,A)')'Removed curvature in ',nbad, &
    ' cells due to poorly conditioned Jacobians'
  CALL oft_warn(error_str)
END IF
CALL oft_global_stitch(mesh%estitch,eflag,1)
CALL oft_global_stitch(mesh%fstitch,fflag,1)
!---Disable bad curvature
IF(mesh%order==2)THEN
  !---Straighten edges
  DO i=1,mesh%ne
    IF(eflag(i)>1.d0)THEN
      mesh%ho_info%r(:,mesh%ho_info%lep(1,i))=(mesh%r(:,mesh%le(1,i))+mesh%r(:,mesh%le(2,i)))/2.d0
    END IF
  END DO
  IF(mesh%ho_info%nfp==1)THEN
    !---Straighten edges
    DO i=1,mesh%nf
      IF(fflag(i)>1.d0)THEN
        mesh%ho_info%r(:,mesh%ho_info%lfp(1,i))=0.d0
        DO j=1,mesh%face_np
          mesh%ho_info%r(:,mesh%ho_info%lfp(1,i))=mesh%ho_info%r(:,mesh%ho_info%lfp(1,i)) &
            + mesh%r(:,mesh%lf(j,i))/REAL(mesh%face_np,8)
        END DO
      END IF
    END DO
  END IF
END IF
DEALLOCATE(eflag,fflag)
!---Boundary adjustment complete
DEBUG_STACK_POP
end subroutine multigrid_check_ho
!---------------------------------------------------------------------------------
!> Set corner flags for distributed and local meshes.
!---------------------------------------------------------------------------------
subroutine multigrid_corners(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
INTEGER(i4) :: i,jc,jr,je(12)
REAL(r8) :: v,v1,v2,goptmp(3,4),f(4),ftmp(4)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: btrans
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: bvin,bvout
REAL(r8), ALLOCATABLE, DIMENSION(:,:,:) :: ct
class(oft_mesh), pointer :: mesh
DEBUG_STACK_PUSH
mesh=>mg_mesh%mesh
if(oft_debug_print(1))write(*,'(2X,A)')'Locating corner features'
ALLOCATE(mesh%cp(mesh%np),mesh%ce(mesh%ne))
mesh%cp=.FALSE.; mesh%ce=.FALSE.
!---
ALLOCATE(ct(3,2,mesh%ne))
ct=0.d0
DO i=1,mesh%nc
  !---
  IF(mesh%global%gbc(i))THEN
    je(1:mesh%cell_ne)=ABS(mesh%lce(:,i))
    DO jc=1,mesh%cell_ne
      CALL mesh%vlog(mesh%cell_ed(1,jc), f)
      CALL mesh%vlog(mesh%cell_ed(2,jc), ftmp)
      f=(f+ftmp)/2.d0
      CALL mesh%jacobian(i,f,goptmp,v)
      DO jr=1,mesh%cell_nf
        IF(ALL(ABS(mesh%cell_fe(:,jr))/=jc))CYCLE
        ! IF(f(jr)>0.d0)cycle
        IF(.NOT.mesh%global%gbf(abs(mesh%lcf(jr,i))))CYCLE
        !!$omp critical
        IF(ALL(ct(:,1,je(jc))==0.d0))THEN
          CALL mesh%snormal(i, jr, f, ct(:,1,je(jc)))
          ct(:,1,je(jc))=ct(:,1,je(jc))/SQRT(SUM(ct(:,1,je(jc))**2))
          ! ct(:,1,je(jc))=-goptmp(:,jr)/SQRT(SUM(goptmp(:,jr)**2))
        ELSE
          CALL mesh%snormal(i, jr, f, ct(:,2,je(jc)))
          ct(:,2,je(jc))=ct(:,2,je(jc))/SQRT(SUM(ct(:,2,je(jc))**2))
          ! ct(:,2,je(jc))=-goptmp(:,jr)/SQRT(SUM(goptmp(:,jr)**2))
        END IF
        !!$omp end critical
      END DO
    END DO
  END IF
END DO
!---
ALLOCATE(bvin(3,mesh%ne),bvout(3,mesh%ne))
bvin=0._r8
bvout=0._r8
!!$omp parallel do
DO i=1,mesh%ne
  IF(ALL(ct(:,2,i)==0.d0))then
    bvout(:,i)=ct(:,1,i)
  END if
END DO
CALL oft_global_stitch(mesh%estitch,bvout(1,:),2) ! Sync boundary edges across processors
CALL oft_global_stitch(mesh%estitch,bvout(2,:),2) ! Sync boundary edges across processors
CALL oft_global_stitch(mesh%estitch,bvout(3,:),2) ! Sync boundary edges across processors
!---
DO i=1,mesh%ne
  IF(ALL(ct(:,1,i)==0.d0).AND.ALL(ct(:,2,i)==0.d0))cycle
  IF(ALL(ct(:,2,i)==0.d0))ct(:,2,i)=bvout(:,i)
  v1=SQRT(DOT_PRODUCT(ct(:,1,i),ct(:,1,i)))
  v2=SQRT(DOT_PRODUCT(ct(:,2,i),ct(:,2,i)))
  IF(v1<1.E-10_r8.OR.v2<1.E-10_r8)CYCLE
  v=DOT_PRODUCT(ct(:,1,i),ct(:,2,i))/(v1*v2)
  IF(ABS(v)>1.d0)CYCLE
  v=ACOS(v)
  IF(ABS(v)>cad_feature_angle)then
    mesh%ce(i)=.TRUE.
    mesh%cp(mesh%le(:,i))=.TRUE.
  END IF
END DO
DEALLOCATE(bvin,bvout,ct)
!---
ALLOCATE(btrans(mesh%ne))
DO i=1,mesh%ne
  btrans(i)=0._r8
  if(mesh%ce(i))btrans(i)=1._r8
END DO
CALL oft_global_stitch(mesh%estitch,btrans,1) ! Sync boundary edges across processors
DO i=1,mesh%nbe
  IF(btrans(mesh%lbe(i))>0.d0)mesh%ce(mesh%lbe(i))=.TRUE. ! Edge is on boundary corner
END DO
DEALLOCATE(btrans)
!---
ALLOCATE(btrans(mesh%np))
DO i=1,mesh%np
  btrans(i)=0._r8
  if(mesh%cp(i))btrans(i)=1._r8
END DO
CALL oft_global_stitch(mesh%pstitch,btrans,1) ! Sync boundary points across processors
DO i=1,mesh%nbp
  IF(btrans(mesh%lbp(i))>0.d0)mesh%cp(mesh%lbp(i))=.TRUE. ! Point is on boundary corner
END DO
DEALLOCATE(btrans)
DEBUG_STACK_POP
end subroutine multigrid_corners
!------------------------------------------------------------------------------
!> Construct multi-level mesh.
!! - Read in mesh options
!! - Load base mesh
!! - Setup local meshes
!! - Decompose mesh
!! - Setup distributed meshes
!------------------------------------------------------------------------------
subroutine multigrid_construct(mg_mesh,grnd_pt)
type(multigrid_mesh), intent(inout) :: mg_mesh !< Multi-level mesh object
real(r8), optional, intent(in) :: grnd_pt(3) !< Reference location for global grounding point (default: [1,0,0])
integer(i4) :: i,level,io_unit
class(oft_mesh), pointer :: mesh
class(oft_bmesh), pointer :: smesh
character(LEN=20) :: meshname
integer(i4) :: ierr
integer(i4) :: cad_type = 0
integer(i4) :: nlevels = 1
integer(i4) :: nbase = 1
integer(i4) :: grid_order = 1
!---Read in mesh options
namelist/mesh_options/meshname,cad_type,nlevels,nbase,grid_order, &
fix_boundary,jac_ratio_tol,part_meth,cad_feature_angle
DEBUG_STACK_PUSH
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,mesh_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr<0)CALL oft_abort('No mesh options found in input file.','multigrid_construct',__FILE__)
!---
if(nlevels<=0)call oft_abort('Invalid number of MG levels','multigrid_construct',__FILE__)
if(nbase<=0.OR.nbase>nlevels)call oft_abort('Invalid number of MG base levels','multigrid_construct',__FILE__)
!---Setup basic mesh structure and runtime
! allocate(mg_mesh)
if(nbase==nlevels)then
  if(oft_env%nprocs/=1)call oft_abort('Requested shared run with more than one proccess.', &
    'multigrid_construct',__FILE__)
else
  if(oft_env%nprocs==1)THEN
#if !defined(HAVE_MPI)
    IF(oft_env%test_run)THEN
      WRITE(*,*)'SKIP TEST'
      CALL oft_finalize
    END IF
#endif
    call oft_abort('Requested distributed run with only one proccess.','multigrid_construct',__FILE__)
  endif
end if
!---Setup MG enviroment
mg_mesh%mgmax=nlevels
mg_mesh%nbase=nbase
oft_env%nbase=nbase
mg_mesh%mgdim=mg_mesh%mgmax
allocate(mg_mesh%inter(mg_mesh%mgdim))
allocate(mg_mesh%sinter(mg_mesh%mgdim))
mg_mesh%level=1
mg_mesh%lev=1
mg_mesh%rlevel='01'
!---Load and initialize mesh on processes
CALL multigrid_load(mg_mesh,cad_type)
mesh=>mg_mesh%mesh
smesh=>mg_mesh%smesh
!---Setup mesh
! mesh%meshname=meshname
!---Load and initialize mesh on processes
smesh%fullmesh=.TRUE.
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(2A,I2)')oft_indent,'**** Generating grid level ',mg_mesh%level
END IF
CALL oft_increase_indent
CALL mesh_global_periodic(mesh)
CALL mesh_global_link(mesh)
CALL mesh_global_boundary(mesh)
call mesh_global_save(mesh)
CALL bmesh_global_init(mesh)
CALL bmesh_global_link(smesh,mesh)
!---
smesh%global%np=smesh%np
ALLOCATE(smesh%global%lp(smesh%np))
DO i=1,smesh%np
  smesh%global%lp(i)=INT(i,8)
END DO
smesh%global%ne=smesh%ne
ALLOCATE(smesh%global%le(smesh%ne))
DO i=1,smesh%ne
  smesh%global%le(i)=INT(i,8)
END DO
smesh%global%nc=smesh%nc
ALLOCATE(smesh%global%lc(smesh%nc))
DO i=1,smesh%nc
  smesh%global%lc(i)=INT(i,8)
END DO
!---
call mesh_global_orient(mesh)
CALL bmesh_global_orient(smesh,mesh)
!---Locate corners
call multigrid_corners(mg_mesh)
!---
call mesh_global_stats(mesh)
call mesh_global_igrnd(mesh,grnd_pt)
!---
i=MAXVAL(mesh%reg)
mesh%nreg=oft_mpi_max(i)
CALL oft_decrease_indent
!---Construct OpenMP base meshes
do level=2,mg_mesh%nbase
  ! call multigrid_add
  CALL add_vol
  CALL add_surf
  CALL multigrid_level(mg_mesh,level)
  mesh=>mg_mesh%mesh
  smesh=>mg_mesh%smesh
  call multigrid_shared_level
end do
!---Construct distributed levels
if(nbase<nlevels)then
  ! call multigrid_add
  CALL add_vol
  CALL add_surf
  CALL multigrid_level(mg_mesh,mg_mesh%nbase+1)
  mesh=>mg_mesh%mesh
  smesh=>mg_mesh%smesh
  call multigrid_decomp
  do level=mg_mesh%nbase+2,mg_mesh%mgmax
    ! call multigrid_add
    CALL add_vol
    CALL add_surf
    CALL multigrid_level(mg_mesh,level)
    mesh=>mg_mesh%mesh
    smesh=>mg_mesh%smesh
    call multigrid_dist_level
  end do
end if
!---Setup high order node points
SELECT CASE(grid_order)
  CASE(1)
    !---Do nothing
  CASE(2)
    CALL multigrid_add_quad(mg_mesh)
  ! CASE(3)
  !   CALL multigrid_add_cubic
  CASE DEFAULT
    CALL oft_abort('Requested invalid grid order.','multigrid_construct',__FILE__)
END SELECT
do level=1,mg_mesh%nbase
  CALL multigrid_level(mg_mesh,level)
  mg_mesh%mesh%global%seam=>mg_mesh%seam
  mg_mesh%smesh%global%seam=>mg_mesh%seam
end do
CALL multigrid_level(mg_mesh,nlevels)
!---Finalize mesh setup for given interface
select case(mg_mesh%mesh%cad_type)
case(mesh_native_id)
  CALL native_finalize_setup
case(mesh_cubit_id)
#ifdef HAVE_NCDF
  CALL cubit_finalize_setup
#else
    CALL oft_abort('CUBIT interface requires NETCDF','multigrid_construct',__FILE__)
#endif
case(mesh_gmsh_id)
  CALL gmsh_finalize_setup
end select
IF(oft_env%head_proc)WRITE(*,*)
DEBUG_STACK_POP
CONTAINS
!
subroutine add_vol
class(oft_mesh), pointer :: fmesh,cmesh
cmesh=>mg_mesh%meshes(mg_mesh%level)   ! Get current mesh
fmesh=>mg_mesh%meshes(mg_mesh%level+1) ! Get new mesh
!---Replicate mesh to next level
fmesh%cad_type=cmesh%cad_type
fmesh%periodic%nper=cmesh%periodic%nper
fmesh%nreg=cmesh%nreg
fmesh%meshname=cmesh%meshname
fmesh%np=cmesh%np
fmesh%nc=cmesh%nc
fmesh%r=>cmesh%r
fmesh%lc=>cmesh%lc
fmesh%reg=>cmesh%reg
fmesh%global%np=cmesh%global%np
fmesh%global%ne=cmesh%global%ne
fmesh%global%nf=cmesh%global%nf
fmesh%global%nc=cmesh%global%nc
fmesh%global%lp=>cmesh%global%lp
fmesh%global%lc=>cmesh%global%lc
end subroutine add_vol
!
subroutine add_surf
class(oft_bmesh), pointer :: fmesh,cmesh
cmesh=>mg_mesh%smeshes(mg_mesh%level)   ! Get current mesh
fmesh=>mg_mesh%smeshes(mg_mesh%level+1) ! Get new mesh
!---Replicate mesh to next level
fmesh%cad_type=cmesh%cad_type
fmesh%np=cmesh%np
fmesh%nc=cmesh%nc
fmesh%r=>cmesh%r
fmesh%lc=>cmesh%lc
fmesh%reg=>cmesh%reg
fmesh%global%np=cmesh%global%np
fmesh%global%ne=cmesh%global%ne
fmesh%global%nc=cmesh%global%nc
fmesh%global%lp=>cmesh%global%lp
fmesh%global%lc=>cmesh%global%lc
end subroutine add_surf
!------------------------------------------------------------------------------
!> Build a new local mesh.
!! - Refine mesh to initialize new level
!! - Build geometry lists
!! - Adjust boundary points
!! - Orient geometry
!! - Print mesh statistics
!------------------------------------------------------------------------------
subroutine multigrid_shared_level
!---Set MG information
mg_mesh%lev=mg_mesh%level
write(mg_mesh%rlevel,'(I2.2)')mg_mesh%lev
mesh=>mg_mesh%mesh
smesh=>mg_mesh%smesh
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(2A,I2)')oft_indent,'**** Generating grid level ',mg_mesh%level
END IF
CALL oft_increase_indent
!---Refine and setup local geometry
mesh%fullmesh=.TRUE.
smesh%fullmesh=.TRUE.
call multigrid_refine(mg_mesh)
call multigrid_brefine(mg_mesh)
call mesh_local_init(mesh)
CALL bmesh_local_init(smesh,mesh)
!---Set global context information
IF(mesh%type==3)THEN
  call hexmesh_mg_globals(mg_mesh,mg_mesh%meshes(mg_mesh%level-1),mesh)
ELSE
  call tetmesh_mg_globals(mg_mesh,mg_mesh%meshes(mg_mesh%level-1),mesh)
END IF
CALL mesh_global_link(mesh)
CALL mesh_global_boundary(mesh)
CALL bmesh_global_link(smesh,mesh)
call mesh_global_save(mesh)
!---Adjust new boundary points
call multigrid_reffix(mg_mesh)
call mesh_volumes(mesh)
!---Orient geometry
call mesh_global_orient(mesh)
CALL bmesh_global_orient(smesh,mesh)
!---Locate corners
call multigrid_corners(mg_mesh)
!---Print mesh statistics
call mesh_global_stats(mesh)
call mesh_global_igrnd(mesh,grnd_pt)
CALL oft_decrease_indent
end subroutine multigrid_shared_level
!------------------------------------------------------------------------------
!> Decompose the local mesh into a distributed mesh.
!! - Decompose and scatter mesh
!! - Build geometry lists
!! - Setup hybrid transfer level
!! - Setup global mesh context
!! - Adjust boundary points
!! - Orient geometry
!! - Print mesh statistics
!------------------------------------------------------------------------------
subroutine multigrid_decomp
!---Set MG information
mg_mesh%lev=mg_mesh%level-1
write(mg_mesh%rlevel,'(I2.2)')mg_mesh%lev
mesh=>mg_mesh%mesh
smesh=>mg_mesh%smesh
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(2A,I2)')oft_indent,'**** Generating grid level ',mg_mesh%level
END IF
CALL oft_increase_indent
!---Decompose and setup local geometry
call mesh_global_decomp(mesh,part_meth)
mg_mesh%seam=>mesh%global%seam
smesh%global%seam=>mg_mesh%seam
call mesh_local_init(mesh)
mesh%fullmesh=.FALSE.
!---Setup hybrid transfer level
call multigrid_hybrid_base(mg_mesh)
!---Set global context information
call mesh_global_link(mesh)
call mesh_global_boundary(mesh)
call mesh_global_save(mesh)
!---
CALL bmesh_global_init(mesh)
smesh%fullmesh=.FALSE.
CALL bmesh_global_link(smesh,mesh)
!---
CALL multigrid_hybrid_bmesh(mg_mesh)
!---Adjust new boundary points
call multigrid_reffix(mg_mesh)
call mesh_volumes(mesh)
!---Orient geometry
call mesh_global_orient(mesh)
CALL bmesh_global_orient(smesh,mesh)
!---Locate corners
call multigrid_corners(mg_mesh)
!---Print mesh statistics
call mesh_global_stats(mesh)
call mesh_global_igrnd(mesh,grnd_pt)
CALL oft_decrease_indent
end subroutine multigrid_decomp
!------------------------------------------------------------------------------
!> Build a new distributed level.
!! - Refine mesh to initialize new level
!! - Build geometry lists
!! - Setup global mesh context
!! - Adjust boundary points
!! - Orient geometry
!! - Print mesh statistics
!------------------------------------------------------------------------------
subroutine multigrid_dist_level
!---Set MG information
mg_mesh%lev=mg_mesh%level-1
write(mg_mesh%rlevel,'(I2.2)')mg_mesh%lev
mesh=>mg_mesh%mesh
smesh=>mg_mesh%smesh
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(2A,I2)')oft_indent,'**** Generating grid level ',mg_mesh%level
END IF
CALL oft_increase_indent
!---
mesh%global%seam=>mg_mesh%seam
smesh%global%seam=>mg_mesh%seam
!---Refine and setup local geometry
mesh%fullmesh=.FALSE.
smesh%fullmesh=.FALSE.
call multigrid_refine(mg_mesh)
call multigrid_brefine(mg_mesh)
call mesh_local_init(mesh)
CALL bmesh_local_init(smesh,mesh)
!---Set global context information
IF(mesh%type==3)THEN
  call hexmesh_mg_globals(mg_mesh,mg_mesh%meshes(mg_mesh%level-1),mesh)
ELSE
  call tetmesh_mg_globals(mg_mesh,mg_mesh%meshes(mg_mesh%level-1),mesh)
END IF
call mesh_global_link(mesh)
call mesh_global_boundary(mesh)
!---
CALL bmesh_global_link(smesh,mesh)
!---
call mesh_global_save(mesh)
!---Adjust new boundary points
call multigrid_reffix(mg_mesh)
call mesh_volumes(mesh)
!---Orient geometry
call mesh_global_orient(mesh)
CALL bmesh_global_orient(smesh,mesh)
!---Locate corners
call multigrid_corners(mg_mesh)
!---Print mesh statistics
call mesh_global_stats(mesh)
call mesh_global_igrnd(mesh,grnd_pt)
CALL oft_decrease_indent
end subroutine multigrid_dist_level
end subroutine multigrid_construct
!------------------------------------------------------------------------------
!> Load in mesh and CAD information.
!! - Read in mesh points and cells
!! - Read in CAD information for refinement
!------------------------------------------------------------------------------
subroutine multigrid_load_surf(mg_mesh,cad_type)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4), intent(in) :: cad_type
integer(i4) :: i
DEBUG_STACK_PUSH
!---Select mesh type and load
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(2A)')oft_indent,'**** Loading OFT surface mesh'
END IF
CALL oft_increase_indent
select case(cad_type)
  case(mesh_native_id) ! Native Mesh
    CALL native_load_smesh(mg_mesh)
    CALL smesh_global_init(mg_mesh%smesh)
    CALL native_hobase(mg_mesh%smesh)
    CALL native_bset_periodic(mg_mesh%smesh)
  case(mesh_t3d_id) ! T3D Mesh
    CALL smesh_t3d_load(mg_mesh)
    CALL smesh_global_init(mg_mesh%smesh)
  !   CALL mesh_t3d_cadsync
  !   CALL mesh_t3d_cadlink
  !   CALL mesh_t3d_set_periodic
  case(mesh_cubit_id) ! Exodus Mesh
#ifdef HAVE_NCDF
    CALL smesh_cubit_load(mg_mesh)
    CALL smesh_global_init(mg_mesh%smesh)
    ! CALL mesh_cubit_cadlink
    CALL mesh_cubit_hobase(mg_mesh%smesh)
    ! CALL mesh_cubit_set_periodic
#else
    CALL oft_abort('CUBIT interface requires NETCDF','multigrid_load_surf',__FILE__)
#endif
  ! case(mesh_gmsh_id) ! GMSH Mesh
  !   CALL mesh_gmsh_load
  !   CALL mesh_global_init
  !   CALL mesh_gmsh_cadlink
  case(mesh_sphere_id) ! Sphere Test Mesh
    CALL smesh_circle_load(mg_mesh)
    CALL smesh_global_init(mg_mesh%smesh)
    CALL smesh_circle_cadlink(mg_mesh%smesh)
  case(mesh_cube_id) ! Cube Test Mesh
    CALL smesh_square_load(mg_mesh)
    CALL smesh_global_init(mg_mesh%smesh)
    CALL smesh_square_cadlink(mg_mesh%smesh)
    CALL smesh_square_set_periodic(mg_mesh%smesh)
  case default ! Invalid Mesh
    CALL oft_abort('Invalid mesh type.','multigrid_load_surf',__FILE__)
end select
CALL oft_decrease_indent
!---Mesh load complete
DEBUG_STACK_POP
end subroutine multigrid_load_surf
!------------------------------------------------------------------------------
!> Adjust boundary points to CAD boundary following refinement.
!------------------------------------------------------------------------------
subroutine multigrid_reffix_surf(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
DEBUG_STACK_PUSH
IF(.NOT.fix_boundary)THEN
  IF(oft_env%head_proc)WRITE(*,*)'Skipping boundary corrections'
  DEBUG_STACK_POP
  RETURN
END IF
!---Always refine using mesh mapping first
CALL multigrid_reffix_ho_surf(mg_mesh)
!---Select mesh type and adjust boundary
select case(mg_mesh%smesh%cad_type)
  case(mesh_native_id)
    ! Do nothing
  case(mesh_t3d_id)
    ! call mesh_t3d_reffix
  case(mesh_cubit_id)
#ifdef HAVE_NCDF
    ! call mesh_cubit_reffix
#endif
  case(mesh_gmsh_id)
    ! call mesh_gmsh_reffix
  case(mesh_sphere_id)
    call smesh_circle_reffix(mg_mesh)
  case(mesh_cube_id)
    ! Do nothing
  case default
    call oft_abort('Invalid mesh type.','multigrid_reffix_surf',__FILE__)
end select
!---Boundary adjustment complete
DEBUG_STACK_POP
end subroutine multigrid_reffix_surf
!------------------------------------------------------------------------------
!> Add node points for quadratic elements
!------------------------------------------------------------------------------
subroutine multigrid_add_quad_surf(smesh)
class(oft_bmesh), intent(inout) :: smesh
INTEGER(i4) :: i,j
REAL(r8) :: pttmp(3)
DEBUG_STACK_PUSH
!---Initialize high order points with straight edges
IF(.NOT.ASSOCIATED(smesh%ho_info%r))THEN
  CALL smesh%set_order(2)
  DO i=1,smesh%ne
    smesh%ho_info%lep(1,i)=i
    smesh%ho_info%r(:,i)=(smesh%r(:,smesh%le(1,i))+smesh%r(:,smesh%le(2,i)))/2.d0
  END DO
  IF(smesh%ho_info%ncp==1)THEN
    do i=1,smesh%nc
      smesh%ho_info%lcp(1,i)=i+smesh%ne
      pttmp=0.d0
      DO j=1,smesh%cell_np
        pttmp=pttmp+smesh%r(:,smesh%lc(j,i))/REAL(smesh%cell_np)
      END DO
      smesh%ho_info%r(:,i+smesh%ne)=pttmp
    end do
  END IF
ELSE
  CALL smesh%set_order(2)
END IF
CALL mesh_global_set_curved(smesh,2)
!---Select mesh type and adjust boundary
select case(smesh%cad_type)
  case(mesh_native_id)
    ! Do nothing
  case(mesh_t3d_id)
    ! call smesh_t3d_add_quad
  case(mesh_cubit_id)
#ifdef HAVE_NCDF
    ! call smesh_cubit_add_quad
#endif
  case(mesh_gmsh_id)
    ! call smesh_gmsh_add_quad
  case(mesh_sphere_id)
    call smesh_circle_add_quad(smesh)
  case(mesh_cube_id)
    ! Do nothing
  case default
    call oft_abort('Invalid mesh type.','multigrid_add_quad',__FILE__)
end select
! CALL multigrid_check_ho
!---Boundary adjustment complete
DEBUG_STACK_POP
end subroutine multigrid_add_quad_surf
!------------------------------------------------------------------------------
!> Construct multi-level surface mesh
!! - Read in mesh options
!! - Load base mesh
!! - Setup local meshes
!! - Decompose mesh
!! - Setup distributed meshes
!------------------------------------------------------------------------------
subroutine multigrid_construct_surf(mg_mesh,grnd_pt)
type(multigrid_mesh), intent(inout) :: mg_mesh !< Multi-level mesh object
real(r8), optional, intent(in) :: grnd_pt(3) !< Reference location for global grounding point (default: [1,0,0])
integer(i4) :: i,level,io_unit
class(oft_bmesh), pointer :: smesh
character(LEN=20) :: meshname
integer(i4) :: ierr
integer(i4) :: cad_type = 0
integer(i4) :: nlevels = 1
integer(i4) :: nbase = 1
integer(i4) :: grid_order = 1
!---Read in mesh options
namelist/mesh_options/meshname,cad_type,nlevels,nbase,grid_order, &
  fix_boundary,jac_ratio_tol,part_meth,cad_feature_angle
DEBUG_STACK_PUSH
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,mesh_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr<0)CALL oft_abort('No mesh options found in input file.','multigrid_construct_surf',__FILE__)
!---
if(nlevels<=0)call oft_abort('Invalid number of MG levels','multigrid_construct',__FILE__)
if(nbase<=0.OR.nbase>nlevels)call oft_abort('Invalid number of MG base levels', &
  'multigrid_construct_surf',__FILE__)
!---Setup basic mesh structure and runtime
! allocate(mg_mesh)
if(nbase==nlevels)then
  if(oft_env%nprocs/=1)call oft_abort('Requested shared run with more than one proccess.', &
    'multigrid_construct_surf',__FILE__)
else
  if(oft_env%nprocs==1)THEN
#if !defined(HAVE_MPI)
    IF(oft_env%test_run)THEN
      WRITE(*,*)'SKIP TEST'
      CALL oft_finalize
    END IF
#endif
    call oft_abort('Requested distributed run with only one proccess.', &
      'multigrid_construct_surf',__FILE__)
  endif
end if
mg_mesh%mgmax=nlevels
mg_mesh%nbase=nbase
oft_env%nbase=nbase
mg_mesh%mgdim=mg_mesh%mgmax
!---Setup MG enviroment
allocate(mg_mesh%sinter(mg_mesh%mgdim))
mg_mesh%level=1
mg_mesh%lev=1
mg_mesh%rlevel='01'
!---Load and initialize mesh on processes
CALL multigrid_load_surf(mg_mesh,cad_type)
smesh=>mg_mesh%smesh
! smesh%cad_type=cad_type
! smesh%meshname=meshname
smesh%fullmesh=.TRUE.
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(2A,I2)')oft_indent,'**** Generating surface grid level ',mg_mesh%level
END IF
CALL oft_increase_indent
! CALL bmesh_local_init(smesh)
CALL bmesh_global_periodic(smesh)
CALL bmesh_global_link(smesh)
CALL bmesh_global_boundary(smesh)
CALL mesh_global_save(smesh)
!---
call bmesh_global_orient(smesh)
!---
call bmesh_global_stats(smesh)
call mesh_global_igrnd(smesh,grnd_pt)
!---
i=MAXVAL(smesh%reg)
smesh%nreg=oft_mpi_max(i)
CALL oft_decrease_indent
!---Construct OpenMP base meshes
do level=2,mg_mesh%nbase
  CALL add_surf
  CALL multigrid_level(mg_mesh,level)
  smesh=>mg_mesh%smesh
  call multigrid_shared_level
end do
!---Construct distributed levels
if(nbase<nlevels)then
  CALL add_surf
  CALL multigrid_level(mg_mesh,mg_mesh%nbase+1)
  smesh=>mg_mesh%smesh
  call multigrid_decomp
  do level=mg_mesh%nbase+2,mg_mesh%mgmax
    CALL add_surf
    CALL multigrid_level(mg_mesh,level)
    smesh=>mg_mesh%smesh
    call multigrid_dist_level
  end do
end if
!---Setup high order node points
SELECT CASE(grid_order)
  CASE(1)
    !---Do nothing
  CASE(2)
    CALL multigrid_add_quad_surf(mg_mesh%smesh)
  CASE DEFAULT
    CALL oft_abort('Requested invalid grid order.','multigrid_construct_surf',__FILE__)
END SELECT
do level=1,mg_mesh%nbase
  CALL multigrid_level(mg_mesh,level)
  mg_mesh%smesh%global%seam=>mg_mesh%seam
end do
CALL multigrid_level(mg_mesh,nlevels)
!---Finalize mesh setup for given interface
select case(mg_mesh%smesh%cad_type)
case(mesh_native_id)
  CALL native_finalize_setup
case(mesh_cubit_id)
#ifdef HAVE_NCDF
  CALL cubit_finalize_setup
#else
    CALL oft_abort('CUBIT interface requires NETCDF','multigrid_construct_surf',__FILE__)
#endif
case(mesh_gmsh_id)
  CALL gmsh_finalize_setup
end select
IF(oft_env%head_proc)WRITE(*,*)
DEBUG_STACK_POP
CONTAINS
!
subroutine add_surf
class(oft_bmesh), pointer :: fmesh,cmesh
cmesh=>mg_mesh%smeshes(mg_mesh%level)   ! Get current mesh
fmesh=>mg_mesh%smeshes(mg_mesh%level+1) ! Get new mesh
!---Replicate mesh to next level
fmesh%cad_type=cmesh%cad_type
fmesh%periodic%nper=cmesh%periodic%nper
fmesh%nreg=cmesh%nreg
fmesh%meshname=cmesh%meshname
fmesh%dim=cmesh%dim
fmesh%np=cmesh%np
fmesh%nc=cmesh%nc
fmesh%r=>cmesh%r
fmesh%lc=>cmesh%lc
fmesh%reg=>cmesh%reg
fmesh%global%np=cmesh%global%np
fmesh%global%ne=cmesh%global%ne
fmesh%global%nc=cmesh%global%nc
fmesh%global%lp=>cmesh%global%lp
fmesh%global%lc=>cmesh%global%lc
end subroutine add_surf
!------------------------------------------------------------------------------
!> Build a new local mesh.
!! - Refine mesh to initialize new level
!! - Build geometry lists
!! - Adjust boundary points
!! - Orient geometry
!! - Print mesh statistics
!------------------------------------------------------------------------------
subroutine multigrid_shared_level
!---Set MG information
mg_mesh%lev=mg_mesh%level
write(mg_mesh%rlevel,'(I2.2)')mg_mesh%lev
smesh=>mg_mesh%smesh
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(2A,I2)')oft_indent,'**** Generating surface grid level ',mg_mesh%level
END IF
CALL oft_increase_indent
!---Refine and setup local geometry
smesh%fullmesh=.TRUE.
call multigrid_brefine(mg_mesh)
CALL bmesh_local_init(smesh)
!---Set global context information
IF(smesh%type==3)THEN
  call quadmesh_mg_globals(mg_mesh,mg_mesh%smeshes(mg_mesh%level-1),smesh)
ELSE
  CALL trimesh_mg_globals(mg_mesh,mg_mesh%smeshes(mg_mesh%level-1),smesh)
END IF
CALL bmesh_global_link(smesh)
CALL bmesh_global_boundary(smesh)
CALL mesh_global_save(smesh)
!---Adjust new boundary points
call multigrid_reffix_surf(mg_mesh)
call bmesh_areas(smesh)
!---Orient geometry
CALL bmesh_global_orient(smesh)
!---Print mesh statistics
call bmesh_global_stats(smesh)
call mesh_global_igrnd(smesh,grnd_pt)
CALL oft_decrease_indent
end subroutine multigrid_shared_level
!------------------------------------------------------------------------------
!> Decompose the local mesh into a distributed mesh.
!! - Decompose and scatter mesh
!! - Build geometry lists
!! - Setup hybrid transfer level
!! - Setup global mesh context
!! - Adjust boundary points
!! - Orient geometry
!! - Print mesh statistics
!------------------------------------------------------------------------------
subroutine multigrid_decomp
!---Set MG information
mg_mesh%lev=mg_mesh%level-1
write(mg_mesh%rlevel,'(I2.2)')mg_mesh%lev
smesh=>mg_mesh%smesh
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(2A,I2)')oft_indent,'**** Generating surface grid level ',mg_mesh%level
END IF
CALL oft_increase_indent
!---Decompose and setup local geometry
call mesh_global_decomp(smesh,part_meth)
mg_mesh%seam=>smesh%global%seam
CALL bmesh_local_init(smesh)
smesh%fullmesh=.FALSE.
!---Setup hybrid transfer level
call multigrid_hybrid_bmesh(mg_mesh)
!---Set global context information
CALL bmesh_global_link(smesh)
CALL bmesh_global_boundary(smesh)
CALL mesh_global_save(smesh)
!---
CALL multigrid_hybrid_bmesh(mg_mesh)
!---Adjust new boundary points
call multigrid_reffix_surf(mg_mesh)
call bmesh_areas(smesh)
!---Orient geometry
CALL bmesh_global_orient(smesh)
!---Print mesh statistics
call bmesh_global_stats(smesh)
call mesh_global_igrnd(smesh,grnd_pt)
CALL oft_decrease_indent
end subroutine multigrid_decomp
!------------------------------------------------------------------------------
!> Build a new distributed level.
!! - Refine mesh to initialize new level
!! - Build geometry lists
!! - Setup global mesh context
!! - Adjust boundary points
!! - Orient geometry
!! - Print mesh statistics
!------------------------------------------------------------------------------
subroutine multigrid_dist_level
!---Set MG information
mg_mesh%lev=mg_mesh%level-1
write(mg_mesh%rlevel,'(I2.2)')mg_mesh%lev
smesh=>mg_mesh%smesh
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(2A,I2)')oft_indent,'**** Generating surface grid level ',mg_mesh%level
END IF
CALL oft_increase_indent
!---
smesh%global%seam=>mg_mesh%seam
!---Refine and setup local geometry
smesh%fullmesh=.FALSE.
call multigrid_brefine(mg_mesh)
CALL bmesh_local_init(smesh)
!---Set global context information
IF(smesh%type==3)THEN
  call quadmesh_mg_globals(mg_mesh,mg_mesh%smeshes(mg_mesh%level-1),smesh)
ELSE
  CALL trimesh_mg_globals(mg_mesh,mg_mesh%smeshes(mg_mesh%level-1),smesh)
END IF
CALL bmesh_global_link(smesh)
CALL bmesh_global_boundary(smesh)
CALL mesh_global_save(smesh)
!---Adjust new boundary points
call multigrid_reffix_surf(mg_mesh)
call bmesh_areas(smesh)
!---Orient geometry
CALL bmesh_global_orient(smesh)
!---Print mesh statistics
call bmesh_global_stats(smesh)
call mesh_global_igrnd(smesh,grnd_pt)
CALL oft_decrease_indent
end subroutine multigrid_dist_level
end subroutine multigrid_construct_surf
end module multigrid_build
