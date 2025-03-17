!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file multigrid.F90
!
!> Multi-Level grid implementation using nested meshes.
!!
!! @authors George Marklin and Chris Hansen
!! @date December 2009
!! @ingroup doxy_oft_grid
!-----------------------------------------------------------------------------
module multigrid
use oft_base
USE oft_sort, ONLY: sort_array, sort_matrix
use oft_mesh_type, only: oft_mesh, oft_bmesh, mesh_seam
use oft_mesh_local_util, only: mesh_local_findedge, mesh_local_findface
use oft_mesh_global_util, only: mesh_global_set_curved
implicit none
#define MOD_NAME "multigrid"
#include "local.h"
!---------------------------------------------------------------------------
!> Level to level interface structure.
!! - Geometry daughter lists
!! - Indexing for local to global interface level
!---------------------------------------------------------------------------
type :: multigrid_inter
  integer(i4), pointer, dimension(:) :: lcdg => NULL() !< List of cell diagonals
  integer(i4), pointer, dimension(:) :: lbege => NULL() !< Edge index on base mesh (defined on interface level only)
  integer(i4), pointer, dimension(:) :: lbfgf => NULL() !< Face index on base mesh (defined on interface level only)
  integer(i4), pointer, dimension(:,:) :: lcde => NULL() !< List of cell daughter edges
  integer(i4), pointer, dimension(:,:) :: lcdf => NULL() !< List of cell daughter faces
  integer(i4), pointer, dimension(:,:) :: lfde => NULL() !< List of face daughter edges
  integer(i4), pointer, dimension(:,:) :: lfdf => NULL() !< List of face daughter faces
  integer(i4), pointer, dimension(:,:) :: lede => NULL() !< List of edge duaghter edges
end type multigrid_inter
!---------------------------------------------------------------------------
!> Multigrid meshes and ML context structure
!! - ML context information (level, nlevels, etc.)
!! - ML mesh array and current level alias
!! - ML interface structure
!---------------------------------------------------------------------------
type :: multigrid_mesh
  integer(i4) :: mgmax = 0 !< Maximum MG level
  integer(i4) :: level = 0 !< Current mesh level
  integer(i4) :: lev = 0 !< Current structure level (doesn't increment on `nbase+1`)
  integer(i4) :: mgdim = 0 !< Size of MG structure
  integer(i4) :: nbase = 0 !< Number of local base refinements
  INTEGER(i4) :: nproc_con = 0 !< Number of processor neighbors
  INTEGER(i4) :: proc_split = 0 !< Location of self in processor list
  INTEGER(i4), POINTER, DIMENSION(:) :: proc_con => NULL() !< Processor neighbor list
#ifdef OFT_MPI_F08
  TYPE(mpi_request), POINTER, DIMENSION(:) :: send_reqs => NULL() !< Asynchronous MPI Send tags
  TYPE(mpi_request), POINTER, DIMENSION(:) :: recv_reqs => NULL() !< Asynchronous MPI Recv tags
#else
  INTEGER(i4), POINTER, DIMENSION(:) :: send_reqs => NULL() !< Asynchronous MPI Send tags
  INTEGER(i4), POINTER, DIMENSION(:) :: recv_reqs => NULL() !< Asynchronous MPI Recv tags
#endif
  class(oft_mesh), pointer :: mesh => NULL() !< Structure containing current mesh
  class(oft_bmesh), pointer :: smesh => NULL() !< Structure containing current mesh
  TYPE(mesh_seam), POINTER :: seam => NULL() !< Global domain-domain connectivity information
  class(oft_mesh), pointer, dimension(:) :: meshes => NULL() !< Structure containing all meshes
  class(oft_bmesh), pointer, dimension(:) :: smeshes => NULL() !< Structure containing current mesh
  type(multigrid_inter), pointer, dimension(:) :: inter => NULL() !< Structure containing linkages
  type(multigrid_inter), pointer, dimension(:) :: sinter => NULL() !< Structure containing linkages
  character(2) :: rlevel !< Character rep of refinement level
end type multigrid_mesh
contains
!---------------------------------------------------------------------------
!> Set mesh level in ML mesh
!---------------------------------------------------------------------------
subroutine multigrid_level(mg_mesh,level)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4), intent(in) :: level !< Desired mesh level
DEBUG_STACK_PUSH
if(level<0.OR.level>mg_mesh%mgmax)call oft_abort('Requested invalid mesh level.','multigrid_level',__FILE__)
mg_mesh%level=level
if(level>mg_mesh%nbase)then
  mg_mesh%lev=level-1
else
  mg_mesh%lev=level
end if
write(mg_mesh%rlevel,'(I2.2)')mg_mesh%lev
IF(ASSOCIATED(mg_mesh%meshes))THEN
  mg_mesh%mesh=>mg_mesh%meshes(level)
  ! mesh=>mg_mesh%mesh
END IF
IF(ASSOCIATED(mg_mesh%smeshes))THEN
  mg_mesh%smesh=>mg_mesh%smeshes(level)
  ! smesh=>mg_mesh%smesh
END IF
! smesh=>mg_mesh%mesh%bmesh
DEBUG_STACK_POP
end subroutine multigrid_level
!---------------------------------------------------------------------------
!> Refine the current mesh level once
!! - Add new points at the center of each edge
!! - Update cell lists
!---------------------------------------------------------------------------
subroutine multigrid_refine(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4) :: i,j,k,reg_tmp,neskip,nfskip,ncskip
integer(i4) :: lccors(4),lecors(2),lcecors(6),lfcors(3),lfecors(3)
real(r8) :: diag(3)
class(oft_mesh), pointer :: cmesh,fmesh
DEBUG_STACK_PUSH
cmesh=>mg_mesh%meshes(mg_mesh%level-1)
fmesh=>mg_mesh%meshes(mg_mesh%level)
NULLIFY(fmesh%r,fmesh%lc)
IF(cmesh%type==3)THEN
  fmesh%global%np=cmesh%global%np+cmesh%global%ne &
    +cmesh%global%nf+cmesh%global%nc
  fmesh%global%nc=cmesh%global%nc*8
  fmesh%np=cmesh%np + cmesh%ne + cmesh%nf + cmesh%nc
  fmesh%ne=cmesh%ne*2 + cmesh%nf*4 + cmesh%nc*6
  fmesh%nf=cmesh%nf*4 + cmesh%nc*12
  fmesh%nc=cmesh%nc*8
  CALL cmesh%tessellate(fmesh%r, fmesh%lc, 2)
  fmesh%lc=fmesh%lc+1
ELSE
  fmesh%global%np=cmesh%global%np+cmesh%global%ne ! Update global point count
  fmesh%global%nc=cmesh%global%nc*8 ! Update global cell count
  fmesh%np=cmesh%np+cmesh%ne
  fmesh%ne=cmesh%ne*2+cmesh%nf*3+cmesh%nc
  fmesh%nf=cmesh%nf*4+cmesh%nc*8
  fmesh%nc=cmesh%nc*8
  CALL cmesh%tessellate(fmesh%r, fmesh%lc, 2)
  fmesh%lc=fmesh%lc+1
END IF
!---Propogate regions
ALLOCATE(fmesh%reg(fmesh%nc))
!$omp parallel do private(j,reg_tmp)
DO i=1,cmesh%nc
  reg_tmp=cmesh%reg(i)
  DO j=1,8
    fmesh%reg((i-1)*8+j)=reg_tmp
  END DO
END DO
DEBUG_STACK_POP
end subroutine multigrid_refine
!------------------------------------------------------------------------------
!> Adjust points to CAD boundary and propogate CAD linkage
!------------------------------------------------------------------------------
subroutine multigrid_reffix_ho(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4) :: i,j,k,l,cell,ho_count,ep(2),fp(4),ed,ed2,nfde,nfdf,ncde,ncdf,fc,cc
real(r8) :: f1(4),f2(4),ftmp(4),pt(3)
class(oft_mesh), pointer :: mesh
class(oft_mesh), pointer :: pmesh
mesh=>mg_mesh%mesh
pmesh=>mg_mesh%meshes(mg_mesh%level-1)
IF(.NOT.ASSOCIATED(pmesh%ho_info%r))RETURN
DEBUG_STACK_PUSH
!---If only one level do nothing
IF(mg_mesh%level==1)THEN
  DEBUG_STACK_POP
  RETURN
END IF
IF(pmesh%fullmesh.AND.(.NOT.mesh%fullmesh))THEN ! Current level is transfer level
  !---Synchronize T3D linkage to distributed mesh
  IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Copying quadratic mapping to distributed mesh'
  CALL oft_increase_indent
  !---Copy quadratic mesh if available
  IF(ASSOCIATED(pmesh%ho_info%r))THEN
    ! mesh%order=pmesh%order
    ho_count = mesh%ne*pmesh%ho_info%nep &
      + mesh%nf*pmesh%ho_info%nfp &
      + mesh%nc*pmesh%ho_info%ncp
    ALLOCATE(mesh%ho_info%r(3,ho_count))
    ! mesh%bmesh%order=pmesh%bmesh%order
    ho_count = mesh%bmesh%ne*pmesh%bmesh%ho_info%nep &
      + mesh%bmesh%nc*pmesh%bmesh%ho_info%ncp
    ALLOCATE(mesh%bmesh%ho_info%r(3,ho_count))
    IF(pmesh%ho_info%nep==1)THEN
      mesh%ho_info%nep=pmesh%ho_info%nep
      ALLOCATE(mesh%ho_info%lep(mesh%ho_info%nep,mesh%ne))
      DO i=1,mesh%ne
        mesh%ho_info%r(:,i) = pmesh%ho_info%r(:,ABS(mesh%base%le(i)))
        mesh%ho_info%lep(1,i) = i
      END DO
      mesh%bmesh%ho_info%nep=pmesh%bmesh%ho_info%nep
      ALLOCATE(mesh%bmesh%ho_info%lep(mesh%bmesh%ho_info%nep,mesh%bmesh%ne))
      DO i=1,mesh%bmesh%ne
        mesh%bmesh%ho_info%r(:,i) = mesh%ho_info%r(:,ABS(mesh%bmesh%parent%le(i)))
        mesh%bmesh%ho_info%lep(1,i) = i
      END DO
    END IF
    IF(pmesh%ho_info%nfp==1)THEN
      mesh%ho_info%nfp=pmesh%ho_info%nfp
      ALLOCATE(mesh%ho_info%lfp(mesh%ho_info%nfp,mesh%nf))
      DO i=1,mesh%nf
        mesh%ho_info%r(:,i+mesh%ne) = pmesh%ho_info%r(:,ABS(mesh%base%lf(i))+pmesh%ne)
        mesh%ho_info%lfp(1,i) = i + mesh%ne
      END DO
      mesh%bmesh%ho_info%ncp=pmesh%bmesh%ho_info%ncp
      ALLOCATE(mesh%bmesh%ho_info%lcp(mesh%bmesh%ho_info%ncp,mesh%bmesh%nc))
      DO i=1,mesh%bmesh%nc
        mesh%bmesh%ho_info%r(:,i+mesh%bmesh%ne) = pmesh%ho_info%r(:,mesh%bmesh%parent%lf(i)+pmesh%ne)
        mesh%bmesh%ho_info%lcp(1,i) = i+mesh%bmesh%ne
      END DO
    END IF
    IF(pmesh%ho_info%ncp==1)THEN
      mesh%ho_info%ncp=pmesh%ho_info%ncp
      ALLOCATE(mesh%ho_info%lcp(mesh%ho_info%ncp,mesh%nc))
      DO i=1,mesh%nc
        mesh%ho_info%r(:,i+mesh%ne+mesh%nf) = pmesh%ho_info%r(:,ABS(mesh%base%lc(i))+pmesh%ne+pmesh%nf)
        mesh%ho_info%lcp(1,i) = i + mesh%ne + mesh%nf
      END DO
    END IF
  END IF
  CALL oft_decrease_indent
  DEBUG_STACK_POP
  RETURN
END IF
!---Copy quadratic mesh if available
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Adjusting points to parent quadratic mesh'
CALL oft_increase_indent
DO i=1,pmesh%ne
  j=i+pmesh%np
  mesh%r(:,j) = pmesh%ho_info%r(:,pmesh%ho_info%lep(1,i))
END DO
IF(pmesh%ho_info%nfp==1)THEN
  DO i=1,pmesh%nf
    j=i+pmesh%ne+pmesh%np
    mesh%r(:,j) = pmesh%ho_info%r(:,pmesh%ho_info%lfp(1,i))
  END DO
END IF
IF(pmesh%ho_info%ncp==1)THEN
  DO i=1,pmesh%nc
    j=i+pmesh%nf+pmesh%ne+pmesh%np
    mesh%r(:,j) = pmesh%ho_info%r(:,pmesh%ho_info%lcp(1,i))
  END DO
END IF
!---Copy updated points to boundary mesh
do i=1,mesh%bmesh%np
  j=mesh%bmesh%parent%lp(i)
  mesh%bmesh%r(:,i)=mesh%r(:,j)
enddo
CALL pmesh%set_order(2)
CALL mesh_global_set_curved(pmesh,2)
!---Create new high order nodes
IF(pmesh%ho_info%nep==1)THEN
  mesh%ho_info%nep=pmesh%ho_info%nep
  ALLOCATE(mesh%ho_info%lep(mesh%ho_info%nep,mesh%ne))
  mesh%ho_info%lep=0
  mesh%bmesh%ho_info%nep=pmesh%bmesh%ho_info%nep
  ALLOCATE(mesh%bmesh%ho_info%lep(mesh%bmesh%ho_info%nep,mesh%bmesh%ne))
  mesh%bmesh%ho_info%lep=0
END IF
IF(pmesh%ho_info%nfp==1)THEN
  mesh%ho_info%nfp=pmesh%ho_info%nfp
  ALLOCATE(mesh%ho_info%lfp(mesh%ho_info%nfp,mesh%nf))
  mesh%ho_info%lfp=0
  mesh%bmesh%ho_info%ncp=pmesh%bmesh%ho_info%ncp
  ALLOCATE(mesh%bmesh%ho_info%lcp(mesh%bmesh%ho_info%ncp,mesh%bmesh%nc))
  mesh%bmesh%ho_info%lcp=0
END IF
IF(pmesh%ho_info%ncp==1)THEN
  mesh%ho_info%ncp=pmesh%ho_info%ncp
  ALLOCATE(mesh%ho_info%lcp(mesh%ho_info%ncp,mesh%nc))
  mesh%ho_info%lcp=0
END IF
! mesh%order=pmesh%order
ho_count = mesh%ne*mesh%ho_info%nep &
  + mesh%nf*mesh%ho_info%nfp &
  + mesh%nc*mesh%ho_info%ncp
ALLOCATE(mesh%ho_info%r(3,ho_count))
mesh%ho_info%r=0.d0
! mesh%bmesh%order=pmesh%bmesh%order
ho_count = mesh%bmesh%ne*pmesh%bmesh%ho_info%nep &
  + mesh%bmesh%nc*pmesh%bmesh%ho_info%ncp
ALLOCATE(mesh%bmesh%ho_info%r(3,ho_count))
mesh%bmesh%ho_info%r=0.d0
!
nfde=SIZE(mg_mesh%inter(mg_mesh%level-1)%lfde,DIM=1)
nfdf=SIZE(mg_mesh%inter(mg_mesh%level-1)%lfdf,DIM=1)
ncde=SIZE(mg_mesh%inter(mg_mesh%level-1)%lcde,DIM=1)
ncdf=SIZE(mg_mesh%inter(mg_mesh%level-1)%lcdf,DIM=1)
!$omp parallel private(cell,k,j,l,ep,fp,f1,f2,ftmp,ed,ed2,fc,cc)
!---Edge children
!$omp do
DO i=1,pmesh%ne
  cell = pmesh%lec(pmesh%kec(i))
  ep=0
  DO k=1,mesh%cell_np
    DO j=1,2
      IF(pmesh%lc(k,cell)==pmesh%le(j,i))ep(j)=k
    END DO
  END DO
  CALL pmesh%vlog(ep(1),f1)
  CALL pmesh%vlog(ep(2),f2)
  f2=(f1+f2)/2.d0 ! Parent edge mid point
  DO l=1,2
    ed = ABS(mg_mesh%inter(mg_mesh%level-1)%lede(l,i))
    f1=0.d0; ftmp=0.d0
    DO j=1,2
      IF(MINVAL(mesh%le(:,ed))==pmesh%le(j,i))THEN
        CALL pmesh%vlog(ep(j),f1)
      END IF
    END DO
    ftmp=(f1+f2)/2.d0
    mesh%ho_info%r(:,ed) = pmesh%log2phys(cell,ftmp)
    mesh%ho_info%lep(1,ed) = ed
  END DO
END DO
!$omp end do nowait
!---Face children
IF(nfde>0)THEN
  !$omp do
  DO i=1,pmesh%nf
    cell = pmesh%lfc(1,i)
    fp=0
    DO k=1,mesh%cell_np
      DO j=1,mesh%face_np
        IF(pmesh%lc(k,cell)==pmesh%lf(j,i))fp(j)=k
      END DO
    END DO
    !---Get location of face center
    f1=0.d0; f2=0.d0; ftmp=0.d0
    DO k=1,mesh%face_np
      CALL pmesh%vlog(fp(k),f1)
      f2=f2+f1/REAL(mesh%face_np,8)
    END DO
    !---Loop over daughter edges
    DO l=1,nfde
      ed = ABS(mg_mesh%inter(mg_mesh%level-1)%lfde(l,i))
      ftmp=0.d0
      DO j=1,2
        IF(mesh%le(j,ed)>pmesh%np+pmesh%ne)THEN
          ftmp=ftmp+f2*0.5d0
        ELSE
          !---Endpoint is mid point of edge
          ed2=mesh%le(j,ed)-pmesh%np
          DO k=1,mesh%face_np
            IF(ANY(pmesh%lf(k,i)==pmesh%le(:,ed2)))THEN
              CALL pmesh%vlog(fp(k),f1)
              ftmp=ftmp+f1/4.d0
            END IF
          END DO
        END IF
      END DO
      mesh%ho_info%r(:,ed) = pmesh%log2phys(cell,ftmp)
      mesh%ho_info%lep(1,ed) = ed
    END DO
  END DO
  !$omp end do nowait
END IF
IF(nfdf>0.AND.pmesh%ho_info%nfp==1)THEN
  !$omp do
  DO i=1,pmesh%nf
    cell = pmesh%lfc(1,i)
    fp=0
    DO k=1,mesh%cell_np
      DO j=1,mesh%face_np
        IF(pmesh%lc(k,cell)==pmesh%lf(j,i))fp(j)=k
      END DO
    END DO
    !---Get location of face center
    f1=0.d0; f2=0.d0; ftmp=0.d0
    DO k=1,mesh%face_np
      CALL pmesh%vlog(fp(k),f1)
      f2=f2+f1/REAL(mesh%face_np,8)
    END DO
    !---Loop over daughter faces
    DO l=1,nfdf
      fc = ABS(mg_mesh%inter(mg_mesh%level-1)%lfdf(l,i))
      ftmp=0.d0
      DO j=1,mesh%face_np
        IF(mesh%lf(j,fc)>pmesh%np+pmesh%ne)THEN
          ftmp=ftmp+f2
        ELSE IF(mesh%lf(j,fc)>pmesh%np)THEN
          !---Endpoint is mid point of edge
          ed2=mesh%lf(j,fc)-pmesh%np
          DO k=1,mesh%face_np
            IF(ANY(pmesh%lf(k,i)==pmesh%le(:,ed2)))THEN
              CALL pmesh%vlog(fp(k),f1)
              ftmp=ftmp+f1/2.d0
            END IF
          END DO
        ELSE
          DO k=1,mesh%face_np
            IF(pmesh%lf(k,i)==mesh%lf(j,fc))THEN
              CALL pmesh%vlog(fp(k),f1)
              ftmp=ftmp+f1
              EXIT
            END IF
          END DO
        END IF
      END DO
      ftmp=ftmp/REAL(mesh%face_np,8)
      mesh%ho_info%r(:,fc+mesh%ne) = pmesh%log2phys(cell,ftmp)
      mesh%ho_info%lfp(1,fc) = fc+mesh%ne
    END DO
  END DO
  !$omp end do nowait
END IF
!---Cell children
IF(ncde>0)THEN
  !$omp do
  DO i=1,pmesh%nc
    !---Get location of center point
    f2=0.d0
    DO k=1,mesh%cell_np
      CALL pmesh%vlog(k,f1)
      f2=f2+f1/REAL(mesh%cell_np,8)
    END DO
    !---Loop over edges
    DO l=1,ncde
      ed = ABS(mg_mesh%inter(mg_mesh%level-1)%lcde(l,i))
      ftmp=0.d0
      DO j=1,2
        IF(mesh%le(j,ed)>pmesh%np+pmesh%ne+pmesh%nf)THEN
          ftmp=ftmp+f2*0.5d0
        ELSE IF(mesh%le(j,ed)<=pmesh%np+pmesh%ne)THEN
          !---Endpoint is mid point of edge
          ed2=mesh%le(j,ed)-pmesh%np
          DO k=1,mesh%cell_np
            IF(ANY(pmesh%lc(k,i)==pmesh%le(:,ed2)))THEN
              CALL pmesh%vlog(k,f1)
              ftmp=ftmp+f1/4.d0
            END IF
          END DO
        ELSE
          !---Endpoint is mid point of face
          ed2=mesh%le(j,ed)-pmesh%np-pmesh%ne
          DO k=1,mesh%cell_np
            IF(ANY(pmesh%lc(k,i)==pmesh%lf(:,ed2)))THEN
              CALL pmesh%vlog(k,f1)
              ftmp=ftmp+f1/(2.d0*REAL(mesh%face_np,8))
            END IF
          END DO
        END IF
      END DO
      mesh%ho_info%r(:,ed) = pmesh%log2phys(i,ftmp)
      mesh%ho_info%lep(1,ed) = ed
    END DO
  END DO
  !$omp end do nowait
END IF
IF(ncdf>0.AND.pmesh%ho_info%nfp==1)THEN
  !$omp do
  DO i=1,pmesh%nc
    !---Get location of center point
    f2=0.d0
    DO k=1,mesh%cell_np
      CALL pmesh%vlog(k,f1)
      f2=f2+f1/REAL(mesh%cell_np,8)
    END DO
    !---Loop over edges
    DO l=1,ncdf
      fc = ABS(mg_mesh%inter(mg_mesh%level-1)%lcdf(l,i))
      ftmp=0.d0
      DO j=1,mesh%face_np
        IF(mesh%lf(j,fc)>pmesh%np+pmesh%ne+pmesh%nf)THEN
          ftmp=ftmp+f2
        ELSE IF(mesh%lf(j,fc)>pmesh%np+pmesh%ne)THEN
          !---Endpoint is mid point of face
          ed2=mesh%lf(j,fc)-pmesh%np-pmesh%ne
          DO k=1,mesh%cell_np
            IF(ANY(pmesh%lc(k,i)==pmesh%lf(:,ed2)))THEN
              CALL pmesh%vlog(k,f1)
              ftmp=ftmp+f1/REAL(mesh%face_np,8)
            END IF
          END DO
        ELSE IF(mesh%lf(j,fc)>pmesh%np)THEN
          !---Endpoint is mid point of edge
          ed2=mesh%lf(j,fc)-pmesh%np
          DO k=1,mesh%cell_np
            IF(ANY(pmesh%lc(k,i)==pmesh%le(:,ed2)))THEN
              CALL pmesh%vlog(k,f1)
              ftmp=ftmp+f1/2.d0
            END IF
          END DO
        ELSE
          DO k=1,mesh%cell_np
            IF(pmesh%lf(k,i)==mesh%lf(j,fc))THEN
              CALL pmesh%vlog(k,f1)
              ftmp=ftmp+f1
              EXIT
            END IF
          END DO
        END IF
      END DO
      ftmp=ftmp/REAL(mesh%face_np,8)
      mesh%ho_info%r(:,fc+mesh%ne) = pmesh%log2phys(i,ftmp)
      mesh%ho_info%lfp(1,fc) = fc+mesh%ne
    END DO
  END DO
  !$omp end do nowait
END IF
IF(pmesh%ho_info%ncp==1)THEN
  !$omp do
  DO i=1,pmesh%nc
    !---Get location of center point
    f2=0.d0
    DO k=1,mesh%cell_np
      CALL pmesh%vlog(k,f1)
      f2=f2+f1/REAL(mesh%cell_np,8)
    END DO
    !---Loop over cells
    DO l=1,8
      cc = (i-1)*8 + l
      ftmp=0.d0
      DO j=1,mesh%cell_np
        IF(mesh%lc(j,cc)>pmesh%np+pmesh%ne+pmesh%nf)THEN
          ftmp=ftmp+f2
        ELSE IF(mesh%lc(j,cc)>pmesh%np+pmesh%ne)THEN
          !---Endpoint is mid point of face
          ed2=mesh%lc(j,cc)-pmesh%np-pmesh%ne
          DO k=1,mesh%cell_np
            IF(ANY(pmesh%lc(k,i)==pmesh%lf(:,ed2)))THEN
              CALL pmesh%vlog(k,f1)
              ftmp=ftmp+f1/REAL(mesh%face_np,8)
            END IF
          END DO
        ELSE IF(mesh%lc(j,cc)>pmesh%np)THEN
          !---Endpoint is mid point of edge
          ed2=mesh%lc(j,cc)-pmesh%np
          DO k=1,mesh%cell_np
            IF(ANY(pmesh%lc(k,i)==pmesh%le(:,ed2)))THEN
              CALL pmesh%vlog(k,f1)
              ftmp=ftmp+f1/2.d0
            END IF
          END DO
        ELSE
          DO k=1,mesh%cell_np
            IF(pmesh%lc(k,i)==mesh%lc(j,cc))THEN
              CALL pmesh%vlog(k,f1)
              ftmp=ftmp+f1
              EXIT
            END IF
          END DO
        END IF
      END DO
      ftmp=ftmp/REAL(mesh%cell_np,8)
      mesh%ho_info%r(:,cc+mesh%ne+mesh%nf) = pmesh%log2phys(i,ftmp)
      mesh%ho_info%lcp(1,cc) = cc+mesh%ne+mesh%nf
    END DO
  END DO
  !$omp end do nowait
END IF
!$omp barrier
!---Propogate to boundary mesh
!$omp do
do i=1,mesh%bmesh%ne
  j=INT(ABS(mesh%bmesh%parent%le(i)),4)
  mesh%bmesh%ho_info%lep(1,i)=i
  mesh%bmesh%ho_info%r(:,i)=mesh%ho_info%r(:,mesh%ho_info%lep(1,j))
enddo
IF(pmesh%bmesh%ho_info%ncp==1)THEN
  !$omp do
  do i=1,mesh%bmesh%nc
    j=INT(ABS(mesh%bmesh%parent%lf(i)),4)
    mesh%bmesh%ho_info%lcp(1,i)=i+pmesh%bmesh%ne
    mesh%bmesh%ho_info%r(:,i+pmesh%bmesh%ne)=mesh%ho_info%r(:,mesh%ho_info%lfp(1,j))
  enddo
END IF
!$omp end parallel
CALL pmesh%set_order(1)
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine multigrid_reffix_ho
!------------------------------------------------------------------------------
!> Adjust points to CAD boundary and propogate CAD linkage
!------------------------------------------------------------------------------
subroutine multigrid_reffix_ho_surf(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4) :: i,j,k,l,cell,ho_count,ep(2),fp(4),ed,ed2,nfde,ncde,cc
real(r8) :: f1(4),f2(4),ftmp(4)
class(oft_bmesh), pointer :: smesh
class(oft_bmesh), pointer :: pmesh
smesh=>mg_mesh%smesh
pmesh=>mg_mesh%smeshes(mg_mesh%level-1)
IF(.NOT.ASSOCIATED(pmesh%ho_info%r))RETURN
DEBUG_STACK_PUSH
!---If only one level do nothing
IF(mg_mesh%level==1)THEN
  DEBUG_STACK_POP
  RETURN
END IF
IF(pmesh%fullmesh.AND.(.NOT.smesh%fullmesh))THEN ! Current level is transfer level
  !---Synchronize T3D linkage to distributed mesh
  IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Copying quadratic mapping to distributed mesh'
  CALL oft_increase_indent
  !---Copy quadratic mesh if available
  IF(ASSOCIATED(pmesh%ho_info%r))THEN
    ! smesh%order=pmesh%order
    ho_count = smesh%ne*pmesh%ho_info%nep &
      + smesh%nc*pmesh%ho_info%ncp
    ALLOCATE(smesh%ho_info%r(3,ho_count))
    IF(pmesh%ho_info%nep==1)THEN
      smesh%ho_info%nep=pmesh%ho_info%nep
      ALLOCATE(smesh%ho_info%lep(smesh%ho_info%nep,smesh%ne))
      smesh%ho_info%lep=-1
      DO i=1,smesh%ne
        smesh%ho_info%r(:,i) = pmesh%ho_info%r(:,ABS(smesh%base%le(i)))
        smesh%ho_info%lep(1,i) = i
      END DO
    END IF
    IF(pmesh%ho_info%ncp==1)THEN
      smesh%ho_info%ncp=pmesh%ho_info%ncp
      ALLOCATE(smesh%ho_info%lcp(smesh%ho_info%ncp,smesh%nc))
      smesh%ho_info%lcp=-1
      DO i=1,smesh%nc
        smesh%ho_info%r(:,i+smesh%ne) = pmesh%ho_info%r(:,ABS(smesh%base%lc(i))+pmesh%ne)
        smesh%ho_info%lcp(1,i) = i + smesh%ne
      END DO
    END IF
  END IF
  CALL oft_decrease_indent
  DEBUG_STACK_POP
  RETURN
END IF
!---Copy quadratic mesh if available
IF(oft_debug_print(1))WRITE(*,'(2A)')oft_indent,'Adjusting points to parent quadratic mesh'
CALL oft_increase_indent
DO i=1,pmesh%ne
  j=i+pmesh%np
  smesh%r(:,j) = pmesh%ho_info%r(:,pmesh%ho_info%lep(1,i))
END DO
IF(pmesh%ho_info%ncp==1)THEN
  DO i=1,pmesh%nc
    j=i+pmesh%ne+pmesh%np
    smesh%r(:,j) = pmesh%ho_info%r(:,pmesh%ho_info%lcp(1,i))
  END DO
END IF
!
! smesh%order=pmesh%order
ho_count = smesh%ne*pmesh%ho_info%nep &
  + smesh%nc*pmesh%ho_info%ncp
ALLOCATE(smesh%ho_info%r(3,ho_count))
smesh%ho_info%r=0.d0
IF(pmesh%ho_info%nep==1)THEN
  smesh%ho_info%nep=pmesh%ho_info%nep
  ALLOCATE(smesh%ho_info%lep(smesh%ho_info%nep,smesh%ne))
  smesh%ho_info%lep=-1
END IF
IF(pmesh%ho_info%ncp==1)THEN
  smesh%ho_info%ncp=pmesh%ho_info%ncp
  ALLOCATE(smesh%ho_info%lcp(smesh%ho_info%ncp,smesh%nc))
  smesh%ho_info%lcp=-1
END IF
CALL pmesh%set_order(2)
CALL mesh_global_set_curved(pmesh,2)
!---Create new high order nodes
ncde=SIZE(mg_mesh%sinter(mg_mesh%level-1)%lcde,DIM=1)
!$omp parallel private(cell,k,j,l,ep,fp,f1,f2,ftmp,ed,ed2,cc)
!---Edge children
!$omp do
DO i=1,pmesh%ne
  cell = pmesh%lec(pmesh%kec(i))
  ep=0
  DO k=1,smesh%cell_np
    DO j=1,2
      IF(pmesh%lc(k,cell)==pmesh%le(j,i))ep(j)=k
    END DO
  END DO
  f1=0.d0; f2=0.d0; ftmp=0.d0
  CALL pmesh%vlog(ep(1),f1)
  CALL pmesh%vlog(ep(2),f2)
  f2=(f1+f2)/2.d0 ! Parent edge mid point
  DO l=1,2
    ed = ABS(mg_mesh%sinter(mg_mesh%level-1)%lede(l,i))
    DO j=1,2
      IF(MINVAL(smesh%le(:,ed))==pmesh%le(j,i))THEN
        CALL pmesh%vlog(ep(j),f1)
      END IF
    END DO
    ftmp=(f1+f2)/2.d0
    smesh%ho_info%r(:,ed) = pmesh%log2phys(cell,ftmp)
    smesh%ho_info%lep(1,ed) = ed
  END DO
END DO
!$omp end do nowait
!---Cell children
IF(ncde>0)THEN
  !$omp do
  DO i=1,pmesh%nc
    !---Get location of center point
    f2=0.d0
    DO k=1,smesh%cell_np
      CALL pmesh%vlog(k,f1)
      f2=f2+f1/REAL(smesh%cell_np,8)
    END DO
    !---Loop over edges
    DO l=1,ncde
      ed = ABS(mg_mesh%sinter(mg_mesh%level-1)%lcde(l,i))
      ftmp=0.d0
      DO j=1,2
        IF(smesh%le(j,ed)>pmesh%np+pmesh%ne)THEN
          ftmp=ftmp+f2*0.5d0
        ELSE
          !---Endpoint is mid point of edge
          ed2=smesh%le(j,ed)-pmesh%np
          DO k=1,smesh%cell_np
            IF(ANY(pmesh%lc(k,i)==pmesh%le(:,ed2)))THEN
              CALL pmesh%vlog(k,f1)
              ftmp=ftmp+f1/4.d0
            END IF
          END DO
        END IF
      END DO
      smesh%ho_info%r(:,ed) = pmesh%log2phys(i,ftmp)
      smesh%ho_info%lep(1,ed) = ed
    END DO
  END DO
  !$omp end do nowait
END IF
IF(smesh%ho_info%ncp==1)THEN
  !$omp do
  DO i=1,pmesh%nc
    !---Get location of center point
    f2=0.d0
    DO k=1,smesh%cell_np
      CALL pmesh%vlog(k,f1)
      f2=f2+f1/REAL(smesh%cell_np,8)
    END DO
    !---Loop over edges
    DO l=1,4
      cc = (i-1)*4+l
      ftmp=0.d0
      DO j=1,2
        IF(smesh%lc(j,cc)>pmesh%np+pmesh%ne)THEN
          ftmp=ftmp+f2
        ELSE IF((smesh%lc(j,cc)>pmesh%np).AND.(smesh%lc(j,cc)<=pmesh%np+pmesh%ne))THEN
          !---Endpoint is mid point of edge
          ed2=smesh%lc(j,cc)-pmesh%np
          DO k=1,smesh%cell_np
            IF(ANY(pmesh%lc(k,i)==pmesh%le(:,ed2)))THEN
              CALL pmesh%vlog(k,f1)
              ftmp=ftmp+f1/2.d0
            END IF
          END DO
        ELSE
          DO k=1,smesh%cell_np
            IF(pmesh%lc(k,i)==smesh%lc(j,cc))THEN
              CALL pmesh%vlog(k,f1)
              ftmp=ftmp+f1
              EXIT
            END IF
          END DO
        END IF
      END DO
      ftmp=ftmp/REAL(smesh%cell_np,8)
      smesh%ho_info%r(:,cc+smesh%ne) = pmesh%log2phys(i,ftmp)
      smesh%ho_info%lcp(1,cc) = cc+smesh%ne
    END DO
  END DO
  !$omp end do nowait
END IF
!$omp end parallel
CALL pmesh%set_order(1)
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine multigrid_reffix_ho_surf
!---------------------------------------------------------------------------
!> Refine the current boundary mesh level once
!! - Add new points at the center of each edge
!! - Update face lists
!---------------------------------------------------------------------------
subroutine multigrid_brefine(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4) :: i,j,reg_tmp,lecors(2),lfcors(3),lfecors(3)
class(oft_bmesh), pointer :: cmesh,fmesh
DEBUG_STACK_PUSH
cmesh=>mg_mesh%smeshes(mg_mesh%level-1)!%bmesh
fmesh=>mg_mesh%smeshes(mg_mesh%level)!%bmesh
IF(cmesh%skip)THEN
  fmesh%skip=.TRUE.
  DEBUG_STACK_POP
  RETURN
END IF
NULLIFY(fmesh%r,fmesh%lc)
IF(cmesh%type==3)THEN
  ! Quadrilateral mesh
  fmesh%global%np=cmesh%global%np+cmesh%global%ne &
    +cmesh%global%nc
  fmesh%global%ne=cmesh%global%ne*2+cmesh%global%nc*4
  fmesh%global%nc=cmesh%global%nc*4
  fmesh%np=cmesh%np + cmesh%ne + cmesh%nc
  fmesh%ne=cmesh%ne*2 + cmesh%nc*4
  fmesh%nc=cmesh%nc*4
  CALL cmesh%tessellate(fmesh%r, fmesh%lc, 2)
  fmesh%lc=fmesh%lc+1
  IF(ASSOCIATED(fmesh%parent))THEN
    allocate(fmesh%parent%lp(fmesh%np))
    fmesh%parent%lp=0
    !$omp parallel do
    do i=1,cmesh%np ! re-use coarse points
      fmesh%parent%lp(i)=cmesh%parent%lp(i)
    enddo
    !$omp parallel do
    do i=1,cmesh%ne ! Loop over edges to set child edge global index
      fmesh%parent%lp(cmesh%np+i)=cmesh%parent%np &
        + ABS(cmesh%parent%le(i))
    enddo
    !$omp parallel do
    do i=1,cmesh%nc ! Loop over faces to set child edge and face global index
      fmesh%parent%lp(cmesh%np+cmesh%ne+i)=cmesh%parent%np &
        + cmesh%parent%ne + cmesh%parent%lf(i)
    enddo
  END IF
ELSE
  ! Triangular mesh
  fmesh%global%np=cmesh%global%np+cmesh%global%ne ! Update global point count
  fmesh%global%ne=cmesh%global%ne*2+cmesh%global%nc*3
  fmesh%global%nc=cmesh%global%nc*4 ! Update global cell count
  fmesh%np=cmesh%np+cmesh%ne
  fmesh%ne=cmesh%ne*2+cmesh%nc*3
  fmesh%nc=cmesh%nc*4
  CALL cmesh%tessellate(fmesh%r, fmesh%lc, 2)
  fmesh%lc=fmesh%lc+1
  !
  IF(ASSOCIATED(fmesh%parent))THEN
    allocate(fmesh%parent%lp(fmesh%np))
    !$omp parallel do
    do i=1,cmesh%np ! re-use coarse points
      fmesh%parent%lp(i)=cmesh%parent%lp(i) ! Update global point index
    enddo
    !$omp parallel do
    do i=1,cmesh%ne ! make new points & edges by dividing coarse edges
      fmesh%parent%lp(cmesh%np+i)=cmesh%parent%np + ABS(cmesh%parent%le(i)) ! Update global point index
    enddo
  END IF
END IF
!---Propogate regions
ALLOCATE(fmesh%reg(fmesh%nc))
!$omp parallel do private(j,reg_tmp)
DO i=1,cmesh%nc
  reg_tmp=cmesh%reg(i)
  DO j=1,4
    fmesh%reg((i-1)*4+j)=reg_tmp
  END DO
END DO
DEBUG_STACK_POP
end subroutine multigrid_brefine
!---------------------------------------------------------------------------
!> Generate a transfer level for local to global mapping.
!! - Populate global indexing for grid block
!---------------------------------------------------------------------------
subroutine multigrid_hybrid_base(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4) :: i,e(2),eind,find,f(4)
integer(i4), allocatable :: lptmp(:)
integer(i4), allocatable :: g_fmap(:),l_fmap(:)
class(oft_mesh), pointer :: lmesh,gmesh
DEBUG_STACK_PUSH
lmesh=>mg_mesh%meshes(mg_mesh%level-1) ! Get parent full mesh
gmesh=>mg_mesh%meshes(mg_mesh%level)   ! Get child mesh partition
!---Propogate point indices to paritioned mesh
allocate(lptmp(gmesh%np))
gmesh%base%np=lmesh%np
!$omp parallel do
do i=1,gmesh%np
  lptmp(i)=INT(gmesh%base%lp(i),4)
enddo
!---Propogate region indices to paritioned mesh
allocate(gmesh%reg(gmesh%nc))
!$omp parallel do
do i=1,gmesh%nc
  gmesh%reg(i)=lmesh%reg(gmesh%base%lc(i))
enddo
!---Index parent edges to paritioned mesh
allocate(gmesh%global%le(gmesh%ne))
allocate(gmesh%base%le(gmesh%ne))
allocate(mg_mesh%inter(mg_mesh%level-1)%lbege(gmesh%ne))
!$omp parallel do private(e,eind)
do i=1,gmesh%ne
  e=lptmp(gmesh%le(:,i))
  ! e(2)=lptmp(gmesh%le(2,i))
  eind=mesh_local_findedge(lmesh,e)
  IF(eind==0)CALL oft_abort("Bad base elink", "multigrid_hybrid_base", __FILE__)
  mg_mesh%inter(mg_mesh%level-1)%lbege(i)=eind
  gmesh%base%le(i)=eind
  gmesh%global%le(i)=abs(lmesh%global%le(abs(eind)))
enddo
!---Index parent faces to paritioned mesh
allocate(gmesh%global%lf(gmesh%nf))
allocate(gmesh%base%lf(gmesh%nf))
allocate(mg_mesh%inter(mg_mesh%level-1)%lbfgf(gmesh%nf))
!---Get face boundary mapping
allocate(g_fmap(gmesh%nf),l_fmap(lmesh%nf))
CALL get_inverse_map(gmesh%lbf,gmesh%nbf,g_fmap,gmesh%nf)
CALL get_inverse_map(lmesh%lbf,lmesh%nbf,l_fmap,lmesh%nf)
!$omp parallel do private(f,find)
do i=1,gmesh%nf
  f(1:gmesh%face_np)=lptmp(gmesh%lf(:,i))
  find=mesh_local_findface(lmesh,f(1:gmesh%face_np))
  IF(find==0)CALL oft_abort("Bad base flink", "multigrid_hybrid_base", __FILE__)
  mg_mesh%inter(mg_mesh%level-1)%lbfgf(i)=find
  gmesh%base%lf(i)=find
  gmesh%global%lf(i)=abs(lmesh%global%lf(abs(find)))
  !---
  IF(g_fmap(i)==0)CYCLE
  IF(l_fmap(ABS(find))==0)CYCLE!CALL oft_abort('Bad face link','multigrid_hybrid_base',__FILE__)
  gmesh%bfs(g_fmap(i))=lmesh%bfs(l_fmap(ABS(find)))
enddo
deallocate(lptmp,g_fmap,l_fmap)
DEBUG_STACK_POP
end subroutine multigrid_hybrid_base
!---------------------------------------------------------------------------
!> Generate a transfer level for local to global mapping for the boundary mesh.
!! - Populate global indexing for grid block
!---------------------------------------------------------------------------
subroutine multigrid_hybrid_bmesh(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4) :: i,e(2),eind
integer(i4), ALLOCATABLE :: lptmp(:),letmp(:),lftmp(:)
class(oft_mesh), pointer :: lmesh_tet,gmesh_tet
class(oft_bmesh), pointer :: lmesh,gmesh
DEBUG_STACK_PUSH
lmesh=>mg_mesh%smeshes(mg_mesh%level-1)
gmesh=>mg_mesh%smeshes(mg_mesh%level)
IF(gmesh%skip)THEN
  DEBUG_STACK_POP
  RETURN
END IF
!---
IF(ASSOCIATED(lmesh%parent))THEN
  lmesh_tet=>mg_mesh%meshes(mg_mesh%level-1) ! Get parent full mesh
  gmesh_tet=>mg_mesh%meshes(mg_mesh%level)   ! Get child mesh partition
  allocate(lptmp(lmesh_tet%np))
  ! CALL bmesh_local_lp2pp(lmesh,lptmp)
  CALL get_inverse_map(lmesh%parent%lp,lmesh%np,lptmp,lmesh%parent%np)
  ALLOCATE(gmesh%base%lp(gmesh%np))
  gmesh%base%lp=0
  !$omp parallel do
  do i=1,gmesh%np
    gmesh%base%lp(i)=lptmp(gmesh_tet%base%lp(gmesh%parent%lp(i)))
  enddo
  DEALLOCATE(lptmp)
  IF(ANY(gmesh%base%lp==0))CALL oft_abort("Error constructing base point mapping", &
    "multigrid_hybrid_bmesh",__FILE__)
  !---
  allocate(letmp(lmesh_tet%ne))
  ! CALL bmesh_local_le2pe(lmesh,letmp)
  CALL get_inverse_map(lmesh%parent%le,lmesh%ne,letmp,lmesh%parent%ne)
  ALLOCATE(gmesh%base%le(gmesh%ne))
  gmesh%base%le=0
  !$omp parallel do
  do i=1,gmesh%ne
    gmesh%base%le(i)=letmp(ABS(gmesh_tet%base%le(gmesh%parent%le(i))))
  enddo
  DEALLOCATE(letmp)
  IF(ANY(gmesh%base%le==0))CALL oft_abort("Error constructing base edge mapping", &
    "multigrid_hybrid_bmesh",__FILE__)
  !---
  allocate(lftmp(lmesh_tet%nf))
  ! CALL bmesh_local_lf2pf(lmesh,lftmp)
  CALL get_inverse_map(lmesh%parent%lf,lmesh%nc,lftmp,lmesh%parent%nf)
  ALLOCATE(gmesh%base%lc(gmesh%nc))
  gmesh%base%lc=0
  !$omp parallel do
  do i=1,gmesh%nc
    gmesh%base%lc(i)=lftmp(ABS(gmesh_tet%base%lf(gmesh%parent%lf(i))))
  enddo
  DEALLOCATE(lftmp)
  IF(ANY(gmesh%base%lp==0))CALL oft_abort("Error constructing base face mapping", &
    "multigrid_hybrid_bmesh",__FILE__)
  !---Propogate point indices to paritioned mesh
  gmesh%global%np=lmesh%global%np
  ALLOCATE(gmesh%global%lp(gmesh%np))
  !$omp parallel do
  do i=1,gmesh%np
    gmesh%global%lp(i)=lmesh%global%lp(gmesh%base%lp(i))
  enddo
  !---
  gmesh%global%ne=lmesh%global%ne
  ALLOCATE(gmesh%global%le(gmesh%ne))
  !$omp parallel do
  do i=1,gmesh%ne
    gmesh%global%le(i)=lmesh%global%le(ABS(gmesh%base%le(i)))
  enddo
  !---
  gmesh%global%nc=lmesh%global%nc
  ALLOCATE(gmesh%global%lc(gmesh%nc))
  !$omp parallel do
  do i=1,gmesh%nc
    gmesh%global%lc(i)=lmesh%global%lc(ABS(gmesh%base%lc(i)))
  enddo
  !---Propogate region indices to paritioned mesh
  allocate(gmesh%reg(gmesh%nc))
  !$omp parallel do
  do i=1,gmesh%nc
    gmesh%reg(i)=lmesh%reg(gmesh%base%lc(i))
  enddo
ELSE
  !---Propogate point indices to paritioned mesh
  allocate(lptmp(gmesh%np))
  gmesh%base%np=lmesh%np
  !$omp parallel do
  do i=1,gmesh%np
    lptmp(i)=INT(gmesh%base%lp(i),4)
  enddo
  !---Propogate region indices to paritioned mesh
  allocate(gmesh%reg(gmesh%nc))
  !$omp parallel do
  do i=1,gmesh%nc
    gmesh%reg(i)=lmesh%reg(gmesh%base%lc(i))
  enddo
  !---Index parent edges to paritioned mesh
  allocate(gmesh%global%le(gmesh%ne))
  allocate(gmesh%base%le(gmesh%ne))
  allocate(mg_mesh%sinter(mg_mesh%level-1)%lbege(gmesh%ne))
  !$omp parallel do private(e,eind)
  do i=1,gmesh%ne
    e=lptmp(gmesh%le(:,i))
    ! e(2)=lptmp(gmesh%le(2,i))
    eind=mesh_local_findedge(lmesh,e)
    IF(eind==0)CALL oft_abort("Bad base elink", "multigrid_hybrid_bmesh", __FILE__)
    mg_mesh%sinter(mg_mesh%level-1)%lbege(i)=eind
    gmesh%base%le(i)=eind
    gmesh%global%le(i)=abs(lmesh%global%le(abs(eind)))
  enddo
END IF
DEBUG_STACK_POP
end subroutine multigrid_hybrid_bmesh
!---------------------------------------------------------------------------
!> Transfer a cell based field from the distributed to base level
!! - Synchronize cell variables to the base mesh
!---------------------------------------------------------------------------
subroutine multigrid_base_pushcc(mg_mesh,bccl,bcc,n)
type(multigrid_mesh), intent(inout) :: mg_mesh
real(r8), intent(in) :: bccl(:,:) !< Cell field on local domain [n,local%nc]
real(r8), intent(out) :: bcc(:,:) !< Cell field on base mesh [n,base%nc]
integer(i4), intent(in) :: n !< Number of values per cell
real(r8), allocatable :: bcctmp(:,:)
integer(i4) :: i,j,ierr,nccors
integer(i4), pointer :: lctmp(:)
DEBUG_STACK_PUSH
if(mg_mesh%level/=mg_mesh%nbase)call oft_abort('Level is not a transfer level.','multigrid_base_pushcc',__FILE__)
lctmp=>mg_mesh%meshes(mg_mesh%nbase+1)%base%lc
nccors=mg_mesh%meshes(mg_mesh%nbase+1)%nc
allocate(bcctmp(n,mg_mesh%mesh%nc))
bcctmp=0.d0
do i=1,nccors
  DO j=1,n
    bcctmp(j,lctmp(i))=bccl(j,i)
  END DO
end do
!---Global reduction over all processors
#ifdef HAVE_MPI
call MPI_ALLREDUCE(bcctmp,bcc,n*mg_mesh%mesh%nc,OFT_MPI_R8,MPI_SUM,oft_env%COMM,ierr)
#else
bcc=bcctmp
#endif
deallocate(bcctmp)
DEBUG_STACK_POP
end subroutine multigrid_base_pushcc
!---------------------------------------------------------------------------
!> Transfer a cell based field from the base level to the distributed mesh
!! - Sample local cell variables from the base level
!---------------------------------------------------------------------------
subroutine multigrid_base_popcc(mg_mesh,bccg,bcc,n)
type(multigrid_mesh), intent(inout) :: mg_mesh
real(r8), intent(in) :: bccg(:,:) !< Cell field on base mesh [n,base%nc]
real(r8), intent(out) :: bcc(:,:) !< Cell field on local domain [n,local%nc]
integer(i4), intent(in) :: n !< Number of values per cell
integer(i4) :: i,j
DEBUG_STACK_PUSH
if(mg_mesh%level/=mg_mesh%nbase+1)call oft_abort('Level is not a transfer level.','multigrid_base_popcc',__FILE__)
do i=1,mg_mesh%mesh%nc
  DO j=1,n
    bcc(j,i)=bccg(j,mg_mesh%mesh%global%lc(i))
  END DO
end do
DEBUG_STACK_POP
end subroutine multigrid_base_popcc
!---------------------------------------------------------------------------
!> Construct multi-level mesh.
!! - Read in mesh options
!! - Load base mesh
!! - Setup local meshes
!! - Decompose mesh
!! - Setup distributed meshes
!---------------------------------------------------------------------------
subroutine multigrid_reset(mg_mesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
integer(i4) :: i,level,io_unit
NULLIFY(mg_mesh%mesh,mg_mesh%smesh)
IF(ASSOCIATED(mg_mesh%smeshes))THEN
  DO i=1,mg_mesh%mgmax
    CALL mg_mesh%smeshes(i)%delete()
  END DO
  DEALLOCATE(mg_mesh%smeshes)
END IF
IF(ASSOCIATED(mg_mesh%meshes))THEN
  DO i=1,mg_mesh%mgmax
    CALL mg_mesh%meshes(i)%delete()
  END DO
  DEALLOCATE(mg_mesh%meshes)
END IF
IF(ASSOCIATED(mg_mesh%sinter))THEN
  DO i=1,mg_mesh%mgmax
    CALL destory_inter(mg_mesh%sinter(i))
  END DO
  DEALLOCATE(mg_mesh%sinter)
END IF
IF(ASSOCIATED(mg_mesh%inter))THEN
  DO i=1,mg_mesh%mgmax
    CALL destory_inter(mg_mesh%inter(i))
  END DO
  DEALLOCATE(mg_mesh%inter)
END IF
! DEALLOCATE(mg_mesh)
!---Reset global environment info (needs to be moved to a mesh-specific object)
oft_env%nbase = -1
mg_mesh%nproc_con = 0
mg_mesh%proc_split = 0
IF(ASSOCIATED(mg_mesh%proc_con))DEALLOCATE(mg_mesh%proc_con,mg_mesh%send_reqs,mg_mesh%recv_reqs)
contains
subroutine destory_inter(obj)
type(multigrid_inter), intent(inout) :: obj
IF(ASSOCIATED(obj%lcdg))DEALLOCATE(obj%lcdg)
IF(ASSOCIATED(obj%lbege))DEALLOCATE(obj%lbege)
IF(ASSOCIATED(obj%lbfgf))DEALLOCATE(obj%lbfgf)
IF(ASSOCIATED(obj%lcde))DEALLOCATE(obj%lcde)
IF(ASSOCIATED(obj%lcdf))DEALLOCATE(obj%lcdf)
IF(ASSOCIATED(obj%lfde))DEALLOCATE(obj%lfde)
IF(ASSOCIATED(obj%lfdf))DEALLOCATE(obj%lfdf)
IF(ASSOCIATED(obj%lede))DEALLOCATE(obj%lede)
end subroutine destory_inter
end subroutine multigrid_reset
!---------------------------------------------------------------------------
!> Update global indices following refinement
!! - Populate new indices using consistent mapping
!---------------------------------------------------------------------------
subroutine tetmesh_mg_globals(mg_mesh,self,fmesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
CLASS(oft_mesh), INTENT(in) :: self
CLASS(oft_mesh), INTENT(inout) :: fmesh
integer(i4), pointer :: lfde(:,:),lede(:,:),lfdf(:,:),lcde(:,:),lcdf(:,:),lcdg(:)
integer(i4) :: i,j,k,l,lecors(2),lfecors(3),lfcors(3),lcfcors(4),lcecors(6),inds(6)
integer(8) :: matrix(2,3)
integer(i4), allocatable :: fmap(:)
DEBUG_STACK_PUSH
!---Allocate linkage arrays
ALLOCATE(mg_mesh%inter(mg_mesh%level-1)%lede(2,self%ne))
lede=>mg_mesh%inter(mg_mesh%level-1)%lede
!$omp parallel do private(lecors)
do i=1,self%ne ! loop over coarse edges & find daughter edges
  lecors=self%le(:,i) ! end points
  lede(1,i)=mesh_local_findedge(fmesh,(/lecors(1), i+self%np/))
  lede(2,i)=mesh_local_findedge(fmesh,(/lecors(2), i+self%np/))
  IF(ANY(lede(:,i)==0))THEN
    WRITE(*,*)i,lede(:,i)
    CALL oft_abort('Bad lede', "tetmesh_mg_globals", __FILE__)
  END IF
enddo
allocate(mg_mesh%inter(mg_mesh%level-1)%lfde(3,self%nf))
allocate(mg_mesh%inter(mg_mesh%level-1)%lfdf(4,self%nf))
lfde=>mg_mesh%inter(mg_mesh%level-1)%lfde
lfdf=>mg_mesh%inter(mg_mesh%level-1)%lfdf
!$omp parallel do private(j,k,l,lfecors,lfcors)
do i=1,self%nf ! loop over coarse faces & find daughter edges
  lfde(:,i)=0
  lfecors=ABS(self%lfe(:,i))+self%np ! edge points
  lfde(1,i)=mesh_local_findedge(fmesh,(/lfecors(3),lfecors(2)/))
  lfde(2,i)=mesh_local_findedge(fmesh,(/lfecors(3),lfecors(1)/))
  lfde(3,i)=mesh_local_findedge(fmesh,(/lfecors(2),lfecors(1)/))
  IF(ANY(lfde(:,i)==0))THEN
    WRITE(*,*)i,lfde(:,i)
    CALL oft_abort('Bad lfde', "tetmesh_mg_globals", __FILE__)
  END IF
  !
  lfdf(:,i)=0
  lfcors=ABS(self%lf(:,i)) ! edge points
  DO j=1,3
    outer: DO k=1,3
      DO l=1,3
        IF(k==l)CYCLE
        lfdf(j,i)=mesh_local_findface(fmesh, (/lfcors(j), lfecors(k), lfecors(l)/))
        IF(lfdf(j,i)/=0)exit outer
      END DO
    END DO outer
  END DO
  lfdf(4,i)=mesh_local_findface(fmesh, lfecors)
  IF(ANY(lfdf(:,i)==0))THEN
    WRITE(*,*)i,lfdf(:,i)
    CALL oft_abort('Bad lfdf', "tetmesh_mg_globals", __FILE__)
  END IF
enddo
allocate(mg_mesh%inter(mg_mesh%level-1)%lcdg(self%nc))
allocate(mg_mesh%inter(mg_mesh%level-1)%lcde(1,self%nc))
allocate(mg_mesh%inter(mg_mesh%level-1)%lcdf(8,self%nc))
lcdg=>mg_mesh%inter(mg_mesh%level-1)%lcdg
lcde=>mg_mesh%inter(mg_mesh%level-1)%lcde
lcdf=>mg_mesh%inter(mg_mesh%level-1)%lcdf
!$omp parallel do private(j,k,l,lcecors)
do i=1,self%nc ! loop over coarse faces & find daughter faces
  lcecors=ABS(self%lce(:,i))+self%np ! edge points
  ! Find diagonal
  lcde(:,i)=0
  DO j=1,3
    k=mesh_local_findedge(fmesh,(/lcecors(j), lcecors(j+3)/))
    IF(k/=0)THEN
      lcdg(i)=j ! diagonal placement index
      lcde(1,i)=k
      k=j
      EXIT
    END IF
  END DO
  !
  ! lcecors=ABS(self%lce(:,i))+self%np ! edge points
  ! k=mg_mesh%inter(mg_mesh%level-1)%lcdg(i) ! diagonal placement index
  ! lcde(1,i)=mesh_local_findedge(fmesh,(/lcecors(k), lcecors(k+3)/))
  IF(ANY(lcde(:,i)==0))THEN
    WRITE(*,*)i,lcde(:,i)
    CALL oft_abort('Bad lcde', "tetmesh_mg_globals", __FILE__)
  END IF
  !
  ! tet_ed(2,6)=RESHAPE((/1,4, 2,4, 3,4, 2,3, 3,1, 1,2/),(/2,6/))
  lcdf(:,i)=0
  lcdf(1,i)=mesh_local_findface(fmesh, (/lcecors(1), lcecors(5), lcecors(6)/))
  lcdf(2,i)=mesh_local_findface(fmesh, (/lcecors(2), lcecors(4), lcecors(6)/))
  lcdf(3,i)=mesh_local_findface(fmesh, (/lcecors(3), lcecors(4), lcecors(5)/))
  lcdf(4,i)=mesh_local_findface(fmesh, (/lcecors(1), lcecors(2), lcecors(3)/))
  !
  j=1
  DO l=1,6
    IF(l==k.OR.l==k+3)CYCLE
    lcdf(4+j,i)=mesh_local_findface(fmesh, (/lcecors(k), lcecors(k+3), lcecors(l)/))
    j=j+1
  END DO
  IF(ANY(lcdf(:,i)==0))THEN
    WRITE(*,*)i,lcdf(:,i)
    CALL oft_abort('Bad lcdf', "tetmesh_mg_globals", __FILE__)
  END IF
enddo
!
! Mesh globals
!
! Initialize global index arrays
allocate(fmesh%global%le(fmesh%ne),fmesh%global%lf(fmesh%nf))
allocate(fmesh%global%lp(fmesh%np),fmesh%global%lc(fmesh%nc))
fmesh%global%lp=0
fmesh%global%le=0
fmesh%global%lf=0
fmesh%global%lc=0
!$omp parallel do
do i=1,self%np ! re-use coarse points
  fmesh%global%lp(i)=self%global%lp(i)
enddo
!$omp parallel do private(j,inds)
do i=1,self%ne ! Loop over edges to set child edge global index
  fmesh%global%lp(self%np+i)=self%global%np + ABS(self%global%le(i))
  inds(1:2)=lede(:,i)
  IF(self%global%le(i)<0)inds(1:2)=inds((/2,1/))
  do j=1,2
    fmesh%global%le(ABS(inds(j)))=j+2*(ABS(self%global%le(i))-1)
  enddo
enddo
!$omp parallel do private(j,matrix,inds)
do i=1,self%nf ! Loop over faces to set child edge and face global index
  do j=1,3 ! Find edges on parent face
    inds(j)=ABS(lfde(j,i))
    matrix(:,j)=fmesh%global%lp(fmesh%le(:,inds(j)))
  enddo
  call sort_matrix(matrix,inds,3_i4) ! Sort edges by endpoint index
  do j=1,3
    fmesh%global%le(ABS(inds(j))) = j+2*self%global%ne+3*(self%global%lf(i)-1)
  enddo
  inds(1:3)=ABS(lfdf(1:3,i))
  CALL orient_listn(self%lfo(i), inds(1:3), 3_i4)
  do j=1,3
    fmesh%global%lf(inds(j)) = j+4*(self%global%lf(i)-1)
  enddo
  fmesh%global%lf(ABS(lfdf(4,i))) = 4+4*(self%global%lf(i)-1)
enddo
!$omp parallel do private(j)
do i=1,self%nc ! Loop over cells to set child edge and face global index
  fmesh%global%le(ABS(lcde(1,i))) = &
    2*self%global%ne+3*self%global%nf+self%global%lc(i)
  do j=1,8
    fmesh%global%lf(ABS(lcdf(j,i))) = j+4*self%global%nf+8*(self%global%lc(i)-1)
  enddo
  do j=1,8
    fmesh%global%lc(j+(i-1)*8) = j+8*(self%global%lc(i)-1)
  end do
enddo
!
! Boundary mesh
!
CALL trimesh_mg_globals(mg_mesh,self%bmesh,fmesh%bmesh)
!---Get edge and face boundary mapping
allocate(fmap(fmesh%nf))
CALL get_inverse_map(fmesh%lbf,fmesh%nbf,fmap,fmesh%nf) ! Get face map
!---Copy face linkage to daughter edges and faces
!$omp parallel do private(i,k)
do j=1,self%nbf
  i=self%lbf(j)
  do k=1,4
    fmesh%bfs(fmap(lfdf(k,i)))=self%bfs(j)
  end do
enddo
deallocate(fmap)
DEBUG_STACK_POP
end subroutine tetmesh_mg_globals
!---------------------------------------------------------------------------
!> Update global indices following refinement.
!! - Populate new indices using consistent mapping
!---------------------------------------------------------------------------
subroutine trimesh_mg_globals(mg_mesh,self,fmesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
CLASS(oft_bmesh), INTENT(in) :: self
CLASS(oft_bmesh), INTENT(inout) :: fmesh
integer(i4), pointer :: lede(:,:),lcde(:,:)
integer(i4) :: i,j,k,l,lecors(2),lcecors(3),inds(6)
! integer(8) :: matrix(2,3)
integer(i4), allocatable :: fmap(:)
DEBUG_STACK_PUSH
!---Allocate linkage arrays
ALLOCATE(mg_mesh%sinter(mg_mesh%level-1)%lede(2,self%ne))
lede=>mg_mesh%sinter(mg_mesh%level-1)%lede
!$omp parallel do private(lecors)
do i=1,self%ne ! loop over coarse edges & find daughter edges
  lecors=self%le(:,i) ! end points
  lede(1,i)=mesh_local_findedge(fmesh,(/lecors(1), i+self%np/))
  lede(2,i)=mesh_local_findedge(fmesh,(/lecors(2), i+self%np/))
  IF(ANY(lede(:,i)==0))THEN
    WRITE(*,*)i,lede(:,i)
    CALL oft_abort('Bad lede', "trimesh_mg_globals", __FILE__)
  END IF
enddo
allocate(mg_mesh%sinter(mg_mesh%level-1)%lcde(3,self%nc))
lcde=>mg_mesh%sinter(mg_mesh%level-1)%lcde
!$omp parallel do private(j,k,l,lcecors)
do i=1,self%nc ! loop over coarse cells & find daughter edges
  lcde(:,i)=0
  lcecors=ABS(self%lce(:,i))+self%np ! edge points
  lcde(1,i)=mesh_local_findedge(fmesh,(/lcecors(3),lcecors(2)/))
  lcde(2,i)=mesh_local_findedge(fmesh,(/lcecors(3),lcecors(1)/))
  lcde(3,i)=mesh_local_findedge(fmesh,(/lcecors(2),lcecors(1)/))
  IF(ANY(lcde(:,i)==0))THEN
    WRITE(*,*)i,lcde(:,i)
    CALL oft_abort('Bad lcde', "trimesh_mg_globals", __FILE__)
  END IF
enddo
! Initialize global index arrays
allocate(fmesh%global%le(fmesh%ne),fmesh%global%lc(fmesh%nc))
allocate(fmesh%global%lp(fmesh%np))
! allocate(fmesh%parent%lp(fmesh%np))
fmesh%global%lp=0
fmesh%global%le=0
fmesh%global%lc=0
! fmesh%parent%lp=0
!$omp parallel do
do i=1,self%np ! re-use coarse points
  fmesh%global%lp(i)=self%global%lp(i)
enddo
!$omp parallel do private(j,k,lecors)
do i=1,self%ne ! Loop over edges to set child edge global index
  fmesh%global%lp(self%np+i)=self%global%np &
    + ABS(self%global%le(i))
    lecors=self%le(:,i)
    IF(self%global%le(i)<0)lecors((/2,1/))=lecors((/1,2/))
  do j=1,2
    k=mesh_local_findedge(fmesh,(/lecors(j), self%np+i/))
    fmesh%global%le(ABS(k))=j+2*(ABS(self%global%le(i))-1)
  enddo
enddo
!$omp parallel do private(j,k,l)
do i=1,self%nc ! Loop over faces to set child edge and face global index
  do j=1,3
    l=j+1
    IF(l>3)l=1
    k=mesh_local_findedge(fmesh,(/ABS(self%lce(j,i))+self%np, &
      ABS(self%lce(l,i))+self%np/))
    IF(k==0)CALL oft_abort("Bad edge link", "tetmesh_mg_globals", __FILE__)
    fmesh%global%le(ABS(k)) = (j+2*self%ne+3*(self%global%lc(i)-1))
  enddo
  do j=1,4
    fmesh%global%lc((i-1)*4+j) = (j+4*(self%global%lc(i)-1))
  enddo
enddo
DEBUG_STACK_POP
end subroutine trimesh_mg_globals
!---------------------------------------------------------------------------
!> Update global indices following refinement.
!! - Populate new indices using consistent mapping
!---------------------------------------------------------------------------
subroutine hexmesh_mg_globals(mg_mesh,self,fmesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
CLASS(oft_mesh), INTENT(in) :: self
CLASS(oft_mesh), INTENT(inout) :: fmesh
integer(i4), pointer :: lfde(:,:),lede(:,:),lfdf(:,:),lcde(:,:),lcdf(:,:)
integer(i4) :: i,j,k,lecors(2),lfecors(4),lfcors(4),lcfcors(6),inds(4)
integer(8) :: lsort(4)
integer(i4), allocatable :: fmap(:)
DEBUG_STACK_PUSH
!---Allocate linkage arrays
ALLOCATE(mg_mesh%inter(mg_mesh%level-1)%lede(2,self%ne))
lede=>mg_mesh%inter(mg_mesh%level-1)%lede
!$omp parallel do private(lecors)
do i=1,self%ne ! loop over coarse edges & find daughter edges
  lecors=self%le(:,i) ! end points
  lede(1,i)=mesh_local_findedge(fmesh,(/lecors(1), i+self%np/))
  lede(2,i)=mesh_local_findedge(fmesh,(/lecors(2), i+self%np/))
  IF(ANY(lede(:,i)==0))THEN
    WRITE(*,*)i,lede(:,i)
    CALL oft_abort('Bad lede', "hexmesh_mg_globals", __FILE__)
  END IF
enddo
allocate(mg_mesh%inter(mg_mesh%level-1)%lfde(4,self%nf))
allocate(mg_mesh%inter(mg_mesh%level-1)%lfdf(4,self%nf))
lfde=>mg_mesh%inter(mg_mesh%level-1)%lfde
lfdf=>mg_mesh%inter(mg_mesh%level-1)%lfdf
!$omp parallel do private(lfecors,lfcors)
do i=1,self%nf ! loop over coarse faces & find daughter edges
  lfde(:,i)=0
  lfecors=ABS(self%lfe(:,i))+self%np ! edge points
  lfde(1,i)=mesh_local_findedge(fmesh,(/lfecors(1), i+self%np+self%ne/))
  lfde(2,i)=mesh_local_findedge(fmesh,(/lfecors(2), i+self%np+self%ne/))
  lfde(3,i)=mesh_local_findedge(fmesh,(/lfecors(3), i+self%np+self%ne/))
  lfde(4,i)=mesh_local_findedge(fmesh,(/lfecors(4), i+self%np+self%ne/))
  IF(ANY(lfde(:,i)==0))THEN
    WRITE(*,*)i,lfde(:,i)
    CALL oft_abort('Bad lfde', "hexmesh_mg_globals", __FILE__)
  END IF
  !
  lfdf(:,i)=0
  lfcors=ABS(self%lf(:,i)) ! edge points
  lfdf(1,i)=mesh_local_findface(fmesh, (/lfcors(1), lfecors(1), &
    i+self%np+self%ne, lfecors(4)/))
  lfdf(2,i)=mesh_local_findface(fmesh, (/lfcors(2), lfecors(2), &
    i+self%np+self%ne, lfecors(1)/))
  lfdf(3,i)=mesh_local_findface(fmesh, (/lfcors(3), lfecors(3), &
    i+self%np+self%ne, lfecors(2)/))
  lfdf(4,i)=mesh_local_findface(fmesh, (/lfcors(4), lfecors(4), &
    i+self%np+self%ne, lfecors(3)/))
  IF(ANY(lfdf(:,i)==0))THEN
    WRITE(*,*)i,lfdf(:,i)
    CALL oft_abort('Bad lfdf', "hexmesh_mg_globals", __FILE__)
  END IF
enddo
allocate(mg_mesh%inter(mg_mesh%level-1)%lcde(6,self%nc))
allocate(mg_mesh%inter(mg_mesh%level-1)%lcdf(12,self%nc))
lcde=>mg_mesh%inter(mg_mesh%level-1)%lcde
lcdf=>mg_mesh%inter(mg_mesh%level-1)%lcdf
!$omp parallel do private(j,k,lcfcors,lfecors)
do i=1,self%nc ! loop over coarse faces & find daughter faces
  lcfcors=ABS(self%lcf(:,i))
  DO j=1,6
    lcde(j,i)=mesh_local_findedge(fmesh,(/lcfcors(j)+self%ne+self%np, &
      i+self%np+self%ne+self%nf/))
  END DO
  IF(ANY(lcde(:,i)==0))THEN
    WRITE(*,*)i,lcde(:,i)
    CALL oft_abort('Bad lcde', "hexmesh_mg_globals", __FILE__)
  END IF
  DO j=1,12
    lfecors(1:2)=0
    DO k=1,6
      IF(ANY(ABS(fmesh%cell_fe(:,k))==j))THEN
        IF(lfecors(1)==0)THEN
          lfecors(1)=ABS(self%lcf(k,i)) + self%ne + self%np
        ELSE
          lfecors(2)=ABS(self%lcf(k,i)) + self%ne + self%np
        END IF
      END IF
    END DO
    IF(lfecors(2)<lfecors(1))lfecors(1:2)=lfecors((/2,1/))
    k = ABS(self%lce(j,i)) + self%np
    lcdf(j,i)=mesh_local_findface(fmesh, (/k, lfecors(1), &
      i+self%np+self%ne+self%nf, lfecors(2)/))
  END DO
  IF(ANY(lcdf(:,i)==0))THEN
    WRITE(*,*)i,lcdf(:,i)
    CALL oft_abort('Bad lcdf', "hexmesh_mg_globals", __FILE__)
  END IF
enddo
!
! Mesh globals
!
! Initialize global index arrays
allocate(fmesh%global%le(fmesh%ne),fmesh%global%lf(fmesh%nf))
allocate(fmesh%global%lp(fmesh%np),fmesh%global%lc(fmesh%nc))
fmesh%global%lp=0
fmesh%global%le=0
fmesh%global%lf=0
fmesh%global%lc=0
!$omp parallel do
do i=1,self%np ! re-use coarse points
  fmesh%global%lp(i)=self%global%lp(i)
enddo
!$omp parallel do private(j,inds)
do i=1,self%ne ! Loop over edges to set child edge global index
  fmesh%global%lp(self%np+i)=self%global%np + ABS(self%global%le(i))
  inds(1:2)=lede(:,i)
  IF(self%global%le(i)<0)inds(1:2)=inds((/2,1/))
  do j=1,2
    fmesh%global%le(ABS(inds(j)))=j+2*(ABS(self%global%le(i))-1)
  enddo
enddo
!$omp parallel do private(j,lsort,inds)
do i=1,self%nf ! Loop over faces to set child edge and face global index
  fmesh%global%lp(self%np+self%ne+i)=self%global%np + self%global%ne + ABS(self%global%lf(i))
  lsort=fmesh%global%lp(ABS(self%lfe(:,i))+self%np)
  inds(1:4)=lfde(:,i)
  CALL sort_array(lsort, inds, 4_i4)
  do j=1,4
    fmesh%global%le(ABS(inds(j))) = SIGN(1,inds(j))*(j+2*self%global%ne+4*(ABS(self%global%lf(i))-1))
  enddo
  lsort=self%global%lp(self%lf(:,i))
  inds(1:4)=lfdf(:,i)
  CALL sort_array(lsort, inds, 4_i4)
  do j=1,4
    fmesh%global%lf(ABS(inds(j))) = (j+4*(ABS(self%global%lf(i))-1))
  enddo
enddo
!$omp parallel do private(j)
do i=1,self%nc ! Loop over cells to set child edge and face global index
  fmesh%global%lp(self%np+self%ne+self%nf+i)=self%global%np + self%global%ne &
    + self%global%nf + self%global%lc(i)
  do j=1,6
    fmesh%global%le(ABS(lcde(j,i))) = SIGN(1,lcde(j,i))*(j+2*self%global%ne+4*self%global%nf+ &
      6*(self%global%lc(i)-1))
  enddo
  do j=1,12
    fmesh%global%lf(ABS(lcdf(j,i))) = (j+4*self%global%nf+12*(ABS(self%global%lc(i))-1))
  enddo
  do j=1,8
    fmesh%global%lc(j+(i-1)*8) = (j+8*(ABS(self%global%lc(i))-1))
  end do
enddo
!
! Boundary mesh
!
CALL quadmesh_mg_globals(mg_mesh,self%bmesh,fmesh%bmesh)
!---Get edge and face boundary mapping
allocate(fmap(fmesh%nf))
CALL get_inverse_map(fmesh%lbf,fmesh%nbf,fmap,fmesh%nf) ! Get face map
!---Copy face linkage to daughter edges and faces
!$omp parallel do private(i,k)
do j=1,self%nbf
  i=self%lbf(j)
  do k=1,4
    fmesh%bfs(fmap(lfdf(k,i)))=self%bfs(j)
  end do
enddo
deallocate(fmap)
DEBUG_STACK_POP
end subroutine hexmesh_mg_globals
!---------------------------------------------------------------------------
!> Update global indices following refinement.
!! - Populate new indices using consistent mapping
!---------------------------------------------------------------------------
subroutine quadmesh_mg_globals(mg_mesh,self,fmesh)
type(multigrid_mesh), intent(inout) :: mg_mesh
CLASS(oft_bmesh), INTENT(in) :: self
CLASS(oft_bmesh), INTENT(inout) :: fmesh
integer(i4), pointer :: lede(:,:),lcde(:,:)
integer(i4) :: i,j,k,lecors(2),lcecors(4),inds(4)
! integer(i8) :: lsort(4)
! integer(i4), allocatable :: fmap(:)
DEBUG_STACK_PUSH
!---Allocate linkage arrays
ALLOCATE(mg_mesh%sinter(mg_mesh%level-1)%lede(2,self%ne))
lede=>mg_mesh%sinter(mg_mesh%level-1)%lede
!$omp parallel do private(lecors)
do i=1,self%ne ! loop over coarse edges & find daughter edges
  lecors=self%le(:,i) ! end points
  lede(1,i)=mesh_local_findedge(fmesh,(/lecors(1), i+self%np/))
  lede(2,i)=mesh_local_findedge(fmesh,(/lecors(2), i+self%np/))
  IF(ANY(lede(:,i)==0))THEN
    WRITE(*,*)i,lede(:,i)
    CALL oft_abort('Bad lede', "quadmesh_mg_globals", __FILE__)
  END IF
enddo
allocate(mg_mesh%sinter(mg_mesh%level-1)%lcde(4,self%nc))
lcde=>mg_mesh%sinter(mg_mesh%level-1)%lcde
!$omp parallel do private(lcecors)
do i=1,self%nc ! loop over coarse faces & find daughter edges
  lcde(:,i)=0
  lcecors=ABS(self%lce(:,i))+self%np ! edge points
  lcde(1,i)=mesh_local_findedge(fmesh,(/lcecors(1), i+self%np+self%ne/))
  lcde(2,i)=mesh_local_findedge(fmesh,(/lcecors(2), i+self%np+self%ne/))
  lcde(3,i)=mesh_local_findedge(fmesh,(/lcecors(3), i+self%np+self%ne/))
  lcde(4,i)=mesh_local_findedge(fmesh,(/lcecors(4), i+self%np+self%ne/))
  IF(ANY(lcde(:,i)==0))THEN
    WRITE(*,*)i,lcde(:,i)
    CALL oft_abort('Bad lcde', "quadmesh_mg_globals", __FILE__)
  END IF
enddo
!
! Boundary mesh
!
! Initialize global index arrays
allocate(fmesh%global%le(fmesh%ne),fmesh%global%lc(fmesh%nc))
allocate(fmesh%global%lp(fmesh%np))
! allocate(fmesh%parent%lp(fmesh%np))
fmesh%global%lp=0
fmesh%global%le=0
fmesh%global%lc=0
! fmesh%parent%lp=0
!$omp parallel do
do i=1,self%np ! re-use coarse points
  fmesh%global%lp(i)=self%global%lp(i)
enddo
!$omp parallel do private(j,k,lecors)
do i=1,self%ne ! Loop over edges to set child edge global index
  fmesh%global%lp(self%np+i)=self%global%np &
    + ABS(self%global%le(i))
  lecors=self%le(:,i)
  IF(self%global%le(i)<0)lecors((/2,1/))=lecors((/1,2/))
  do j=1,2
    k=mesh_local_findedge(fmesh,(/lecors(j), self%np+i/))
    fmesh%global%le(ABS(k))=SIGN(1,k)*(j+2*(ABS(self%global%le(i))-1))
  enddo
enddo
!$omp parallel do private(j,k)
do i=1,self%nc ! Loop over faces to set child edge and face global index
  fmesh%global%lp(self%np+self%ne+i)=self%global%np &
    + self%global%ne + self%global%lc(i)
  do j=1,4
    k=mesh_local_findedge(fmesh,(/ABS(self%lce(j,i))+self%np, &
      self%np+self%ne+i/))
    IF(k==0)CALL oft_abort("Bad edge link", "quadmesh_mg_globals", __FILE__)
    fmesh%global%le(ABS(k)) = SIGN(1,k)*(j+2*self%global%ne+4*(ABS(self%global%lc(i))-1))
  enddo
  do j=1,4
    fmesh%global%lc((i-1)*4+j) = (j+4*(ABS(self%global%lc(i))-1))
  enddo
enddo
DEBUG_STACK_POP
end subroutine quadmesh_mg_globals
end module multigrid
