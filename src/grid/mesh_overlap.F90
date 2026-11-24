!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file mesh_overlap.F90
!
!> Mesh helpers for growing overlapping subdomains by layers of cells across processes
!!
!! The base mesh (`meshes(nbase)`) is replicated on every process and the owner of
!! each base cell is known everywhere (`meshes(nbase+1)%base%lcpart`). Each base cell
!! and all of its descendants live on exactly one process, so the base mesh bounds,
!! without communication, which remote cells a process can reach within a given
!! number of cell layers on a distributed level:
!! - @ref mesh_base_candidates finds the foreign base cells within a number of point
!!   layers of the base cells owned by this process, together with their owners.
!! - @ref mesh_cell_base_index maps the cells of a distributed level to their base cell
!!   using the deterministic global cell numbering of uniform refinement
!!   (child `j` of cell `i` has global index `j + nchild*(gl(i)-1)`).
!! - @ref mesh_base_gid_inverse inverts the global cell numbering of the base mesh.
!!
!! On shared runs (or levels `<= nbase`, `mesh%fullmesh`) every process holds the whole
!! mesh, no remote cells are ever needed and both helpers return empty results.
!!
!! @authors Chris Hansen
!! @date November 2026
!! @ingroup doxy_oft_grid
!---------------------------------------------------------------------------------
MODULE oft_mesh_overlap
USE oft_base
USE oft_mesh_type, ONLY: oft_mesh
USE multigrid, ONLY: multigrid_mesh
IMPLICIT NONE
#include "local.h"
PRIVATE
!> Number of children of a volume (tet or hex) cell under uniform refinement
INTEGER(i4), PARAMETER, PUBLIC :: mesh_overlap_nchild = 8
PUBLIC mesh_has_transfer_level, mesh_base_gid_inverse, mesh_base_candidates, mesh_cell_base_index
CONTAINS
!---------------------------------------------------------------------------------
!> Does the multi-level mesh have distributed levels (a transfer level `nbase+1`)?
!---------------------------------------------------------------------------------
FUNCTION mesh_has_transfer_level(ml_mesh) RESULT(has_transfer)
TYPE(multigrid_mesh), INTENT(in) :: ml_mesh !< Multi-level mesh
LOGICAL :: has_transfer !< True if levels above `nbase` exist and are distributed
has_transfer=(ml_mesh%nbase<ml_mesh%mgmax)
IF(has_transfer)has_transfer=(.NOT.ml_mesh%meshes(ml_mesh%nbase+1)%fullmesh)
END FUNCTION mesh_has_transfer_level
!---------------------------------------------------------------------------------
!> Invert the global cell numbering of the base mesh (`meshes(nbase)%global%lc`)
!!
!! Aborts if the numbering is not a permutation of `1..nc`.
!---------------------------------------------------------------------------------
SUBROUTINE mesh_base_gid_inverse(ml_mesh,inv)
TYPE(multigrid_mesh), INTENT(in) :: ml_mesh !< Multi-level mesh
INTEGER(i4), ALLOCATABLE, DIMENSION(:), INTENT(out) :: inv !< Base cell index of each base global cell id [nc_base]
INTEGER(i4) :: i,nc_base
INTEGER(i8) :: gl
CLASS(oft_mesh), POINTER :: bmesh
DEBUG_STACK_PUSH
bmesh=>ml_mesh%meshes(ml_mesh%nbase)
nc_base=bmesh%nc
IF(bmesh%global%nc/=INT(nc_base,i8))CALL oft_abort('Base mesh is not complete', &
  'mesh_base_gid_inverse',__FILE__)
ALLOCATE(inv(nc_base))
inv=0
DO i=1,nc_base
  gl=bmesh%global%lc(i)
  IF(gl<1_i8.OR.gl>INT(nc_base,i8))CALL oft_abort('Base cell global index out of range', &
    'mesh_base_gid_inverse',__FILE__)
  IF(inv(gl)/=0)CALL oft_abort('Duplicate base cell global index','mesh_base_gid_inverse',__FILE__)
  inv(gl)=i
END DO
DEBUG_STACK_POP
END SUBROUTINE mesh_base_gid_inverse
!---------------------------------------------------------------------------------
!> Find the foreign base cells within `layers` point layers of this process' base cells
!!
!! Runs a multi-source breadth-first search over the point adjacency (`kpc/lpc`) of the
!! replicated base mesh, seeded by the base cells owned by this process
!! (`lcpart==rank+1`). Every cell reached within `layers` layers that is owned by
!! another process is returned (sorted by owner, then by base cell index).
!!
!! Each layer of point-adjacent cells on a uniformly refined level stays within the
!! point-adjacent base cells, so `layers` fine layers of any adjacency never leave the
!! base cells found with the same number of layers here.
!!
!! No communication is needed. Returns no candidates on shared meshes.
!---------------------------------------------------------------------------------
SUBROUTINE mesh_base_candidates(ml_mesh,layers,ncand,cand,cowner)
TYPE(multigrid_mesh), INTENT(in) :: ml_mesh !< Multi-level mesh
INTEGER(i4), INTENT(in) :: layers !< Number of point layers to grow by
INTEGER(i4), INTENT(out) :: ncand !< Number of candidate base cells
INTEGER(i4), ALLOCATABLE, DIMENSION(:), INTENT(out) :: cand !< Candidate base cell indices (in `meshes(nbase)`) [ncand]
INTEGER(i4), ALLOCATABLE, DIMENSION(:), INTENT(out) :: cowner !< Owner rank of each candidate [ncand]
INTEGER(i4) :: i,j,k,l,c,nc_base,me,nfront,nnew,q
INTEGER(i4), POINTER, DIMENSION(:) :: lcpart
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: lvl,front,newf,ncnt
CLASS(oft_mesh), POINTER :: bmesh
DEBUG_STACK_PUSH
ncand=0
IF(.NOT.mesh_has_transfer_level(ml_mesh))THEN
  ALLOCATE(cand(0),cowner(0))
  DEBUG_STACK_POP
  RETURN
END IF
IF(layers<0)CALL oft_abort('Negative number of layers','mesh_base_candidates',__FILE__)
me=oft_env%rank
bmesh=>ml_mesh%meshes(ml_mesh%nbase)
nc_base=bmesh%nc
lcpart=>ml_mesh%meshes(ml_mesh%nbase+1)%base%lcpart
IF(.NOT.ASSOCIATED(lcpart))CALL oft_abort('Base cell partition not available','mesh_base_candidates',__FILE__)
IF(SIZE(lcpart)/=nc_base)CALL oft_abort('Base cell partition size mismatch','mesh_base_candidates',__FILE__)
!---Multi-source BFS seeded by owned base cells
ALLOCATE(lvl(nc_base),front(nc_base),newf(nc_base))
lvl=-1
nfront=0
DO i=1,nc_base
  IF(lcpart(i)==me+1)THEN
    lvl(i)=0
    nfront=nfront+1
    front(nfront)=i
  END IF
END DO
DO l=1,layers
  nnew=0
  DO i=1,nfront
    c=front(i)
    DO j=1,bmesh%cell_np
      k=bmesh%lc(j,c)
      DO q=bmesh%kpc(k),bmesh%kpc(k+1)-1
        IF(lvl(bmesh%lpc(q))<0)THEN
          lvl(bmesh%lpc(q))=l
          nnew=nnew+1
          newf(nnew)=bmesh%lpc(q)
        END IF
      END DO
    END DO
  END DO
  IF(nnew==0)EXIT
  nfront=nnew
  front(1:nnew)=newf(1:nnew)
END DO
!---Foreign cells reached, grouped by owner
ALLOCATE(ncnt(0:oft_env%nprocs))
ncnt=0
DO i=1,nc_base
  IF(lvl(i)>0)THEN
    IF(lcpart(i)<1.OR.lcpart(i)>oft_env%nprocs)CALL oft_abort('Invalid base cell owner', &
      'mesh_base_candidates',__FILE__)
    ncnt(lcpart(i))=ncnt(lcpart(i))+1
  END IF
END DO
ncand=SUM(ncnt)
DO q=1,oft_env%nprocs
  ncnt(q)=ncnt(q)+ncnt(q-1)
END DO
ALLOCATE(cand(ncand),cowner(ncand))
DO i=1,nc_base
  IF(lvl(i)>0)THEN
    q=lcpart(i)-1
    ncnt(q)=ncnt(q)+1
    cand(ncnt(q))=i
    cowner(ncnt(q))=q
  END IF
END DO
DEALLOCATE(lvl,front,newf,ncnt)
DEBUG_STACK_POP
END SUBROUTINE mesh_base_candidates
!---------------------------------------------------------------------------------
!> Map the cells of a level mesh to their base cell (index in `meshes(nbase)`)
!!
!! For a distributed level the number of refinements `nref` above the transfer level
!! is found from `mesh%global%nc == nchild**nref * nc_base` (aborts if no such `nref`),
!! and the base cell global index of a cell with global index `gl` is
!! `(gl-1)/nchild**nref + 1`. At the transfer level (`nref==0`) the result is
!! cross-checked against `mesh%base%lc`. Every cell must lie in a base cell owned by
!! this process.
!!
!! Shared meshes (`mesh%fullmesh`) have no base mapping: `nref=-1` and `cb=0`.
!---------------------------------------------------------------------------------
SUBROUTINE mesh_cell_base_index(ml_mesh,mesh,cb,nref)
TYPE(multigrid_mesh), INTENT(in) :: ml_mesh !< Multi-level mesh
CLASS(oft_mesh), INTENT(in) :: mesh !< Level mesh (one of `ml_mesh%meshes`)
INTEGER(i4), ALLOCATABLE, DIMENSION(:), INTENT(out) :: cb !< Base cell index of each cell [mesh%nc]
INTEGER(i4), INTENT(out) :: nref !< Number of refinements above the transfer level (-1 on shared meshes)
INTEGER(i4) :: i,nc_base
INTEGER(i8) :: nc_ref,nfac,gb
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: inv
INTEGER(i4), POINTER, DIMENSION(:) :: lcpart
CLASS(oft_mesh), POINTER :: bmesh
DEBUG_STACK_PUSH
ALLOCATE(cb(mesh%nc))
cb=0
nref=-1
IF(mesh%fullmesh)THEN
  DEBUG_STACK_POP
  RETURN
END IF
IF(.NOT.mesh_has_transfer_level(ml_mesh))CALL oft_abort('Distributed mesh without transfer level', &
  'mesh_cell_base_index',__FILE__)
bmesh=>ml_mesh%meshes(ml_mesh%nbase)
nc_base=bmesh%nc
lcpart=>ml_mesh%meshes(ml_mesh%nbase+1)%base%lcpart
CALL mesh_base_gid_inverse(ml_mesh,inv)
!---Refinement depth from the global cell counts (uniform refinement)
nref=0
nfac=1
nc_ref=INT(nc_base,i8)
DO WHILE(nc_ref<mesh%global%nc)
  nref=nref+1
  nfac=nfac*mesh_overlap_nchild
  nc_ref=nc_ref*mesh_overlap_nchild
END DO
IF(nc_ref/=mesh%global%nc)CALL oft_abort('Level is not a uniform refinement of the base mesh', &
  'mesh_cell_base_index',__FILE__)
IF(nref>ml_mesh%mgmax-ml_mesh%nbase-1)CALL oft_abort('Refinement depth exceeds distributed levels', &
  'mesh_cell_base_index',__FILE__)
DO i=1,mesh%nc
  gb=(mesh%global%lc(i)-1_i8)/nfac+1_i8
  IF(gb<1_i8.OR.gb>INT(nc_base,i8))CALL oft_abort('Base cell global index out of range', &
    'mesh_cell_base_index',__FILE__)
  cb(i)=inv(gb)
  IF(lcpart(cb(i))/=oft_env%rank+1)CALL oft_abort('Cell lies in a base cell owned by another process', &
    'mesh_cell_base_index',__FILE__)
END DO
!---Cross-check at the transfer level
IF(nref==0)THEN
  IF(.NOT.ASSOCIATED(mesh%base%lc))CALL oft_abort('Transfer level without base cell map', &
    'mesh_cell_base_index',__FILE__)
  DO i=1,mesh%nc
    IF(ABS(mesh%base%lc(i))/=cb(i))CALL oft_abort('Base cell map mismatch at transfer level', &
      'mesh_cell_base_index',__FILE__)
  END DO
END IF
DEALLOCATE(inv)
DEBUG_STACK_POP
END SUBROUTINE mesh_cell_base_index
END MODULE oft_mesh_overlap
