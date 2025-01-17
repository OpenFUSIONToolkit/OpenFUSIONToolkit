!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!------------------------------------------------------------------------------
!> @file oft_mesh_type.F90
!
!> @defgroup doxy_oft_grid Meshing
!! Modules containing the Open FUSION Toolkit meshing framework.
!
!> Tetrahedral mesh structure definitions
!! - MPI seam
!! - MPI I/O index
!! - MPI global context
!! - Base mesh linkage
!! - MPI preallocated Send/Recv
!! - High order tet representation
!! - Mesh container
!!
!! Global Tet variables
!! - Tetrahedra edge and face lists
!! - Grounding point position
!!
!! @author Chris Hansen
!! @date June 2010
!! @ingroup doxy_oft_grid
!------------------------------------------------------------------------------
MODULE oft_mesh_type
USE oft_base
USE oft_stitching, ONLY: oft_seam, destory_seam
USE oft_io, ONLY: hdf5_create_files, oft_hdf5_write_dump, oft_hdf5_add_dump, &
  hdf5_proc_str, hdf5_ts_str, hdf5_write
USE oft_quadrature
IMPLICIT NONE
#include "local.h"
PRIVATE
PUBLIC cell_is_curved, mesh_findcell, mesh_findcell2, bmesh_findcell
!------------------------------------------------------------------------------
! TYPE mesh_per
!------------------------------------------------------------------------------
!> Global mesh information and indicies
!!
!! Contains global mexh context information.
!! - Global counts
!! - Global element indices
!! - Global boundary flags
!------------------------------------------------------------------------------
TYPE, PUBLIC :: mesh_per
  INTEGER(i4) :: nper = 0 !< Number of periodic directions
  INTEGER(i4), POINTER, DIMENSION(:) :: lp => NULL() !< Global index of points (np)
  INTEGER(i4), POINTER, DIMENSION(:) :: le => NULL() !< Global index of edges (ne)
  INTEGER(i4), POINTER, DIMENSION(:) :: lf => NULL() !< Global index of faces (nf)
END TYPE mesh_per
!------------------------------------------------------------------------------
! TYPE mesh_save_index
!------------------------------------------------------------------------------
!> MPI global index information (For I/O Only)
!!
!! Contains global indices for HDF5 I/O.
!! - Local starting and closing index
!! - Maximum counts
!------------------------------------------------------------------------------
TYPE, PUBLIC :: mesh_save_index
  INTEGER(i8) :: npb = 0 !< Index of first point on proc
  INTEGER(i8) :: neb = 0 !< Index of first edge on proc
  INTEGER(i8) :: nfb = 0 !< Index of first face on proc
  INTEGER(i8) :: ncb = 0 !< Index of first cell on proc
  INTEGER(i8) :: npl = 0 !< Index of last point on proc
  INTEGER(i8) :: nel = 0 !< Index of last edge on proc
  INTEGER(i8) :: nfl = 0 !< Index of last face on proc
  INTEGER(i8) :: nbf = 0 !< Number of global boundary faces on proc
  INTEGER(i8) :: npmax = 0 !< Max # of points on one proc
  INTEGER(i8) :: nemax = 0 !< Max # of edges on one proc
  INTEGER(i8) :: nfmax = 0 !< Max # of faces on one proc
  INTEGER(i8) :: ncmax = 0 !< Max # of cells on one proc
END TYPE mesh_save_index
!------------------------------------------------------------------------------
! TYPE mpi_global
!------------------------------------------------------------------------------
!> Global mesh information and indicies
!!
!! Contains global mexh context information.
!! - Global counts
!! - Global element indices
!! - Global boundary flags
!------------------------------------------------------------------------------
TYPE, PUBLIC :: mesh_global
  INTEGER(i8) :: np = 0 !< Global point count
  INTEGER(i8) :: ne = 0 !< Global edge count
  INTEGER(i8) :: nf = 0 !< Global face count
  INTEGER(i8) :: nc = 0 !< Global cell count
  LOGICAL, POINTER, DIMENSION(:) :: gbp => NULL() !< Global boundary point flag (np)
  LOGICAL, POINTER, DIMENSION(:) :: gbe => NULL() !< Global boundary edge flag (ne)
  LOGICAL, POINTER, DIMENSION(:) :: gbf => NULL() !< Global boundary face flag (nf)
  LOGICAL, POINTER, DIMENSION(:) :: gbc => NULL() !< Global boundary cell flag (nc)
  INTEGER(i8), POINTER, DIMENSION(:) :: lp => NULL() !< Global index of points (np)
  INTEGER(i8), POINTER, DIMENSION(:) :: le => NULL() !< Global index of edges (ne) [oriented]
  INTEGER(i8), POINTER, DIMENSION(:) :: lf => NULL() !< Global index of faces (nf)
  INTEGER(i8), POINTER, DIMENSION(:) :: lc => NULL() !< Global index of cells (nc)
END TYPE mesh_global
!------------------------------------------------------------------------------
! TYPE mpi_base
!------------------------------------------------------------------------------
!> Base mesh information and indicies
!!
!! Contains global mexh context information.
!! - Global counts
!! - Global element indices
!! - Global boundary flags
!------------------------------------------------------------------------------
TYPE, PUBLIC :: mesh_base
  INTEGER(i4) :: np = 0 !< Global point count
  INTEGER(i4) :: nc = 0 !< Global cell count
  INTEGER(i4), POINTER, DIMENSION(:) :: lp => NULL() !< Global index of points (np)
  INTEGER(i4), POINTER, DIMENSION(:) :: le => NULL() !< Global index of points (np)
  INTEGER(i4), POINTER, DIMENSION(:) :: lf => NULL() !< Global index of points (np)
  INTEGER(i4), POINTER, DIMENSION(:) :: lc => NULL() !< Global index of cells (nc)
  INTEGER(i4), POINTER, DIMENSION(:) :: lcpart => NULL() !< Global index of cells (nc)
END TYPE mesh_base
!------------------------------------------------------------------------------
! TYPE ho_mesh
!------------------------------------------------------------------------------
!> High order tet geometry information
!!
!! Contains additional data for high order tetrahedra.
!! - Number of additional nodes per geometry primative
!! - List of high-order node locations
!! - List of node points
!------------------------------------------------------------------------------
TYPE, PUBLIC :: ho_mesh
  INTEGER(i4) :: nep = 0 !< Number of nodes per edge
  INTEGER(i4) :: nfp = 0 !< Number of nodes per face
  INTEGER(i4) :: ncp = 0 !< Number of nodes per cell
  LOGICAL, POINTER, DIMENSION(:) :: is_curved => NULL() !< Flag indicating cell is curved
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lep => NULL() !< List of edge points
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lfp => NULL() !< List of face points
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lcp => NULL() !< List of cell points
  REAL(r8), POINTER, DIMENSION(:,:) :: r => NULL() !< List of high-order points
END TYPE ho_mesh
!------------------------------------------------------------------------------
!> Abstrac mesh type (surface or volume)
!!
!! Contains geometry information for the computational grid.
!! - Entity counts
!! - Mesh type and order
!! - Global mesh information
!! - Linkage of geometric primatives
!------------------------------------------------------------------------------
TYPE, PUBLIC, ABSTRACT :: oft_amesh
  LOGICAL :: fullmesh = .TRUE. !< Local mesh flag (False if distributed)
  INTEGER(i4) :: cad_type = 1 !< Type of CAD geometry
  INTEGER(i4) :: type = 0 !< Mesh type
  INTEGER(i4) :: order = 1 !< order of boundary tets (default=linear)
  INTEGER(i4) :: tess_order = 1 !< order of boundary tets (default=linear)
  INTEGER(i4) :: cell_np = 0 !< Number of points per cell
  INTEGER(i4) :: cell_ne = 0 !< Number of edged per cell
  INTEGER(i4) :: np = 0 !< Number of points
  INTEGER(i4) :: ne = 0 !< Number of edges
  INTEGER(i4) :: nc = 0 !< Number of cells
  INTEGER(i4) :: npc = 0 !< Number of point to cell interactions
  INTEGER(i4) :: nec = 0 !< Number of edge to cell interactions
  INTEGER(i4) :: npp = 0 !< Number of point to point interactions
  INTEGER(i4) :: nee = 0 !< Number of edge to edge interactions
  INTEGER(i4) :: npe = 0 !< Number of point to edge interactions
  INTEGER(i4) :: nbp = 0 !< Number of boundary points
  INTEGER(i4) :: nbe = 0 !< Number of boundary edges
  INTEGER(i4) :: nbc = 0 !< Number of boundary cells
  INTEGER(i4) :: igrnd(2) = 0  !< Index of ground point
  INTEGER(i4) :: nreg = 1  !< Number of mesh regions
  INTEGER(i4) :: nparts = 0 !< Number of local partitions
  REAL(r8) :: hmin = 0.d0 !< Minimum edge length
  REAL(r8) :: hrms = 0.d0 !< Mean edge length
  REAL(r8) :: hmax = 0.d0 !< Maximum edge length
  LOGICAL, POINTER, DIMENSION(:) :: bp => NULL() !< Boundary point flag
  LOGICAL, POINTER, DIMENSION(:) :: be => NULL() !< Boundary edge flag
  LOGICAL, POINTER, DIMENSION(:) :: bc => NULL() !< Boundary cell flag
  LOGICAL, POINTER, DIMENSION(:) :: cp => NULL() !< Corner point flag
  LOGICAL, POINTER, DIMENSION(:) :: ce => NULL() !< Corner edge flag
  INTEGER(i4), POINTER, DIMENSION(:) :: reg => NULL() !< Cell list of region ids
  INTEGER(i4), POINTER, DIMENSION(:) :: kpc => NULL() !< Pointer to point shared cell list
  INTEGER(i4), POINTER, DIMENSION(:) :: lpc => NULL() !< Linkage of points to shared cells
  INTEGER(i4), POINTER, DIMENSION(:) :: kec => NULL() !< Pointer to edge shared cell list
  INTEGER(i4), POINTER, DIMENSION(:) :: lec => NULL() !< Linkage of edges to shared cells
  INTEGER(i4), POINTER, DIMENSION(:) :: kpp => NULL() !< Pointer to point shared point list
  INTEGER(i4), POINTER, DIMENSION(:) :: lpp => NULL() !< Linkage of points to shared points
  INTEGER(i4), POINTER, DIMENSION(:) :: kee => NULL() !< Pointer to edge shared edge list
  INTEGER(i4), POINTER, DIMENSION(:) :: lee => NULL() !< Linkage of edges to shared edges
  INTEGER(i4), POINTER, DIMENSION(:) :: lbp => NULL() !< List of boundary points
  INTEGER(i4), POINTER, DIMENSION(:) :: lbe => NULL() !< List of boundary edges
  INTEGER(i4), POINTER, DIMENSION(:) :: lbc => NULL() !< List of boundary cells
  INTEGER(i4), POINTER, DIMENSION(:) :: kpe => NULL() !< Pointer to point shared edge list
  INTEGER(i4), POINTER, DIMENSION(:) :: lpe => NULL() !< Linkage of points to shared edges
  INTEGER(i4), POINTER, DIMENSION(:) :: klpe => NULL() !< Pointer to low point shared edge list
  INTEGER(i4), POINTER, DIMENSION(:) :: llpe => NULL() !< Linkage of low points to shared edge
  INTEGER(i4), POINTER, DIMENSION(:,:) :: cell_ed => NULL() !< Cell edge end points (local numbering)
  INTEGER(i4), POINTER, DIMENSION(:,:) :: le => NULL() !< List of edge points
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lc => NULL() !< List of cell points
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lce => NULL() !< List of cell edges (oriented)
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lcc => NULL() !< List of cell neighbor cells
  REAL(r8), POINTER :: r(:,:) => NULL() !< List of node locations
  CHARACTER(LEN=OFT_PATH_SLEN) :: meshname = 'none' !< Meshname for mesh
  CHARACTER(LEN=OFT_PATH_SLEN) :: filename = 'none' !< Filename for mesh
  CHARACTER(LEN=OFT_PATH_SLEN) :: io_path = '' !< Base path for I/O (relative to working directory)
  TYPE(mesh_global) :: global !< Global mesh information
  TYPE(mesh_base) :: base !< Base mesh information
  TYPE(mesh_save_index) :: save !< Processor to processor linkage information
  TYPE(oft_seam) :: pstitch
  TYPE(oft_seam) :: estitch
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: tloc_p => NULL() !< Point thread ownership
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: tloc_e => NULL() !< Edge thread ownership
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: tloc_c => NULL() !< Cell thread ownership
  TYPE(mesh_per) :: periodic !< Periodic information
  TYPE(ho_mesh) :: ho_info !< High order geometry information
END TYPE oft_amesh
!------------------------------------------------------------------------------
!> Global mesh information and indicies
!!
!! Contains global mexh context information.
!! - Global counts
!! - Global element indices
!! - Global boundary flags
!------------------------------------------------------------------------------
TYPE, PUBLIC :: bmesh_parent
  INTEGER(i4) :: np = 0 !< Global point count
  INTEGER(i4) :: ne = 0 !< Global edge count
  INTEGER(i4) :: nf = 0 !< Global face count
  INTEGER(i4), POINTER, DIMENSION(:) :: lp => NULL() !< Global index of points (np)
  INTEGER(i4), POINTER, DIMENSION(:) :: le => NULL() !< Global index of edges (ne) [oriented]
  INTEGER(i4), POINTER, DIMENSION(:) :: lf => NULL() !< Global index of faces (nf)
END TYPE bmesh_parent
!------------------------------------------------------------------------------
!> Surface mesh type
!!
!! Contains geometry information for the computational grid.
!! - Entity counts
!! - Mesh type and order
!! - Global mesh information
!! - Linkage of geometric primatives
!------------------------------------------------------------------------------
TYPE, PUBLIC, ABSTRACT, EXTENDS(oft_amesh) :: oft_bmesh
  LOGICAL :: skip = .FALSE. !< Skip mesh on this processor
  INTEGER(i4) :: dim = 3 !< Spatial dimension of grid (2 or 3)
  INTEGER(i4), POINTER, DIMENSION(:) :: lco => NULL() !< List of face orientations
  INTEGER(i4), POINTER, DIMENSION(:) :: bes => NULL() !< Surface ID of boundary edges
  REAL(r8), POINTER, DIMENSION(:) :: ca => NULL() !< Cell volumes
  REAL(r8), POINTER, DIMENSION(:) :: va => NULL() !< Node vertex volumes
  TYPE(bmesh_parent), POINTER :: parent => NULL() !< Parent mesh information
CONTAINS
  !> Setup mesh
  PROCEDURE(bmesh_setup), DEFERRED :: setup
  !> Save mesh to simple text format
  PROCEDURE(bmesh_save), DEFERRED :: save_to_file
  !> Save mesh to simple text format
  PROCEDURE(bmesh_load), DEFERRED :: load_from_file
  !> Set geometric mapping order
  PROCEDURE(bmesh_set_order), DEFERRED :: set_order
  !> Invert the sense of a given cell
  PROCEDURE(bmesh_invert_face), DEFERRED :: invert_face
  !> Convert logical to physical coordinates
  PROCEDURE(bmesh_log2phys), DEFERRED :: log2phys
  !> Convert physical to logical coordinates
  PROCEDURE(bmesh_phys2log), DEFERRED :: phys2log
  !> Get logical to physical jacobian matrix
  PROCEDURE(bmesh_jacobian), DEFERRED :: jacobian
  !> Get logical to physical hessian matrices
  PROCEDURE(bmesh_hessian), DEFERRED :: hessian
  !> Get surface unit normal
  PROCEDURE(bmesh_norm), DEFERRED :: norm
  !> Get surface tangent basis 
  PROCEDURE(bmesh_tang), DEFERRED :: tang
  !> Needs docs
  PROCEDURE(bmesh_vlog), DEFERRED :: vlog
  !> Needs docs
  PROCEDURE(bmesh_in_cell), DEFERRED :: in_cell
  !> Get quadrature rule for surface mesh
  PROCEDURE(bmesh_quad_rule), DEFERRED :: quad_rule
  !> Tessellate mesh
  PROCEDURE(bmesh_tessellate), DEFERRED :: tessellate
  !> Get vertex and cell counts for tessallated I/O mesh
  PROCEDURE(bmesh_get_io_sizes), DEFERRED :: get_io_sizes
  !> Setup I/O files for surface mesh
  PROCEDURE :: setup_io => bmesh_setup_io
  !> Save cell-centered scalar field
  PROCEDURE :: save_cell_scalar => bmesh_save_cell_scalar
  !> Save cell-centered vector field
  PROCEDURE :: save_cell_vector => bmesh_save_cell_vector
  !> Save vertex-centered scalar field
  PROCEDURE :: save_vertex_scalar => bmesh_save_vertex_scalar
  !> Save vertex-centered vector field
  PROCEDURE :: save_vertex_vector => bmesh_save_vertex_vector
  !> Compute surface area of mesh
  PROCEDURE :: area => bmesh_area
  !> Delete mesh object
  PROCEDURE :: delete => bmesh_destroy
END TYPE oft_bmesh
! Class procedure interfaces
ABSTRACT INTERFACE
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
  SUBROUTINE bmesh_setup(self,cad_type,has_parent)
  IMPORT oft_bmesh, i4
  CLASS(oft_bmesh), INTENT(inout) :: self
  INTEGER(i4), INTENT(in) :: cad_type
  LOGICAL, INTENT(in) :: has_parent
  END SUBROUTINE bmesh_setup
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
  SUBROUTINE bmesh_save(self,filename)
  IMPORT oft_bmesh
  CLASS(oft_bmesh), INTENT(in) :: self
  CHARACTER(LEN=*), INTENT(in) :: filename
  END SUBROUTINE bmesh_save
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
  SUBROUTINE bmesh_load(self,filename)
  IMPORT oft_bmesh
  CLASS(oft_bmesh), INTENT(inout) :: self
  CHARACTER(LEN=*), INTENT(in) :: filename
  END SUBROUTINE bmesh_load
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
  SUBROUTINE bmesh_set_order(self,order)
  IMPORT oft_bmesh, i4
  CLASS(oft_bmesh), INTENT(inout) :: self
  INTEGER(i4), INTENT(in) :: order
  END SUBROUTINE bmesh_set_order
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
  SUBROUTINE bmesh_invert_face(self,i)
  IMPORT oft_bmesh, i4
  CLASS(oft_bmesh), INTENT(inout) :: self
  INTEGER(i4), INTENT(in) :: i
  END SUBROUTINE bmesh_invert_face
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
  FUNCTION bmesh_log2phys(self,cell,f) RESULT(pt)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: cell
  REAL(r8), INTENT(in) :: f(:)
  REAL(r8) :: pt(3)
  END FUNCTION bmesh_log2phys
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
  SUBROUTINE bmesh_phys2log(self,cell,pt,f)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: cell
  REAL(r8), INTENT(in) :: pt(3)
  REAL(r8), INTENT(out) :: f(:)
  END SUBROUTINE bmesh_phys2log
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
  SUBROUTINE bmesh_jacobian(self,cell,f,gop,j)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: cell
  REAL(r8), INTENT(in) :: f(:)
  REAL(r8), INTENT(out) :: gop(:,:)
  REAL(r8), INTENT(out) :: j
  END SUBROUTINE bmesh_jacobian
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
  SUBROUTINE bmesh_hessian(self,cell,f,g2op,K)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: cell
  REAL(r8), INTENT(in) :: f(:)
  REAL(r8), INTENT(out) :: g2op(:,:)
  REAL(r8), INTENT(out) :: K(:,:)
  END SUBROUTINE bmesh_hessian
!---------------------------------------------------------------------------
!> Get unit normal for surface at a given point
!!
!! @param[in] self Mesh object
!! @param[in] i Cell containing point
!! @param[in] f Logical coordinates in cell
!! @param[out] n Unit normal [3]
!---------------------------------------------------------------------------
  SUBROUTINE bmesh_norm(self,i,f,n)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), TARGET, INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: i
  REAL(r8), INTENT(in) :: f(:)
  REAL(r8), INTENT(out) :: n(3)
  END SUBROUTINE bmesh_norm
!---------------------------------------------------------------------------
!> Get tangent basis set for surface at a given point
!!
!! @param[in] self Mesh object
!! @param[in] i Cell containing point
!! @param[in] f Logical coordinates in cell
!! @param[out] t Unit tangent basis set [3,2]
!---------------------------------------------------------------------------
  SUBROUTINE bmesh_tang(self,i,f,t)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), TARGET, INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: i
  REAL(r8), INTENT(in) :: f(:)
  REAL(r8), INTENT(out) :: t(3,2)
  END SUBROUTINE bmesh_tang
!---------------------------------------------------------------------------
!> Get quadrature rule for a given order in logical coordinates
!!
!! @param[in] self Mesh object
!! @param[in] order Desired order of quadrature rule
!! @param[out] quad_rule Quadrature rule
!---------------------------------------------------------------------------
  SUBROUTINE bmesh_quad_rule(self,order,quad_rule)
  IMPORT oft_bmesh, oft_quad_type, i4
  CLASS(oft_bmesh), INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: order
  TYPE(oft_quad_type), INTENT(out) :: quad_rule
  END SUBROUTINE bmesh_quad_rule
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
  SUBROUTINE bmesh_vlog(self,i,f)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: i
  REAL(r8), INTENT(out) :: f(:)
  END SUBROUTINE bmesh_vlog
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
  FUNCTION bmesh_in_cell(self,f,tol) RESULT(eedge)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), INTENT(in) :: self
  REAL(r8), INTENT(in) :: f(:)
  REAL(r8), INTENT(in) :: tol
  INTEGER(i4) :: eedge
  END FUNCTION bmesh_in_cell
!---------------------------------------------------------------------------
!> Tessellate grid onto Lagrange node points for a given order
!!
!! @param[in] self Mesh to tessellate
!! @param[out] rtmp Vertices for tessellation [3,:]
!! @param[out] lftmp Cell list for tessellation [face_np,:]
!! @param[in] order Order of tessellation
!---------------------------------------------------------------------------
  SUBROUTINE bmesh_tessellate(self,rtmp,lctmp,order)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), INTENT(in) :: self
  REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp
  INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp
  INTEGER(i4), INTENT(in) :: order
  END SUBROUTINE bmesh_tessellate
!---------------------------------------------------------------------------
!> Get variable sizes following tessellation
!---------------------------------------------------------------------------
  function bmesh_get_io_sizes(self) result(sizes)
  IMPORT oft_bmesh, i4
  CLASS(oft_bmesh), INTENT(in) :: self
  integer(i4) :: sizes(2)
  end function bmesh_get_io_sizes
END INTERFACE
!------------------------------------------------------------------------------
! TYPE oft_mesh
!------------------------------------------------------------------------------
!> Tetrahedral Mesh type
!!
!! Contains geometry information for the computational grid.
!! - Entity counts
!! - Mesh type and order
!! - Global mesh information
!! - Linkage of geometric primatives
!------------------------------------------------------------------------------
TYPE, PUBLIC, ABSTRACT, EXTENDS(oft_amesh) :: oft_mesh
  INTEGER(i4) :: cell_nf = 0 !< Number of faces per cell
  INTEGER(i4) :: face_np = 0 !< Number of points per face
  INTEGER(i4) :: nf = 0 !< Number of faces
  INTEGER(i4) :: nbf = 0 !< Number of boundary faces
  LOGICAL, POINTER, DIMENSION(:) :: bf => NULL() !< Boundary face flag
  INTEGER(i4), POINTER, DIMENSION(:) :: lbf => NULL() !< List of boundary faces
  INTEGER(i4), POINTER, DIMENSION(:) :: klef => NULL() !< Pointer to low edge shared face list
  INTEGER(i4), POINTER, DIMENSION(:) :: llef => NULL() !< Linkage of low edges to shared
  INTEGER(i4), POINTER, DIMENSION(:) :: lfo => NULL() !< List of cell faces orientations
  INTEGER(i4), POINTER, DIMENSION(:) :: bfs => NULL() !< Surface ID of boundary faces
  INTEGER(i4), POINTER, DIMENSION(:,:) :: cell_fc => NULL() !< Cell face corner points (local numbering)
  INTEGER(i4), POINTER, DIMENSION(:,:) :: cell_fe => NULL() !< Cell face edges (local numbering)
  INTEGER(i4), POINTER, DIMENSION(:,:) :: face_ed => NULL() !< Face edge indices
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lf => NULL() !< List of face points
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lfe => NULL() !< List of face edges
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lfc => NULL() !< List of faces neighbor cells
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lcf => NULL() !< List of cell faces
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lcfo => NULL() !< List of cell faces orientations
  REAL(r8), POINTER :: cv(:) => NULL() !< Cell volumes
  REAL(r8), POINTER :: vv(:) => NULL() !< Node vertex volumes
  TYPE(oft_seam) :: fstitch
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: tloc_f => NULL() !< Face thread ownership
  CLASS(oft_bmesh), POINTER :: bmesh => NULL() !< Boundary mesh
CONTAINS
  PROCEDURE(mesh_setup), DEFERRED :: setup
  PROCEDURE(mesh_set_order), DEFERRED :: set_order
  PROCEDURE(mesh_invert_cell), DEFERRED :: invert_cell
  PROCEDURE(mesh_log2phys), DEFERRED :: log2phys
  PROCEDURE(mesh_phys2log), DEFERRED :: phys2log
  PROCEDURE(mesh_jacobian), DEFERRED :: jacobian
  PROCEDURE(mesh_hessian), DEFERRED :: hessian
  PROCEDURE(mesh_snormal), DEFERRED :: snormal
  PROCEDURE(mesh_ctang), DEFERRED :: ctang
  PROCEDURE(mesh_get_surf_map), DEFERRED :: get_surf_map
  PROCEDURE(mesh_surf_to_vol), DEFERRED :: surf_to_vol
  PROCEDURE(mesh_vlog), DEFERRED :: vlog
  PROCEDURE(mesh_in_cell), DEFERRED :: in_cell
  PROCEDURE(mesh_quad_rule), DEFERRED :: quad_rule
  PROCEDURE(mesh_tessellate), DEFERRED :: tessellate
  PROCEDURE(mesh_get_io_sizes), DEFERRED :: get_io_sizes
  PROCEDURE :: setup_io => mesh_setup_io
  PROCEDURE :: save_cell_scalar => mesh_save_cell_scalar
  PROCEDURE :: save_cell_vector => mesh_save_cell_vector
  PROCEDURE :: save_vertex_scalar => mesh_save_vertex_scalar
  PROCEDURE :: save_vertex_vector => mesh_save_vertex_vector
  PROCEDURE :: volume => mesh_volume
  !> Delete mesh object
  PROCEDURE :: delete => mesh_destroy
END TYPE oft_mesh
!
ABSTRACT INTERFACE
  !
  SUBROUTINE mesh_setup(self,cad_type)
  IMPORT oft_mesh, i4
  CLASS(oft_mesh), INTENT(inout) :: self
  INTEGER(i4), INTENT(in) :: cad_type
  END SUBROUTINE mesh_setup
  !
  SUBROUTINE mesh_set_order(self,order)
  IMPORT oft_mesh, i4
  CLASS(oft_mesh), INTENT(inout) :: self
  INTEGER(i4), INTENT(in) :: order
  END SUBROUTINE mesh_set_order
  !
  SUBROUTINE mesh_invert_cell(self,i)
  IMPORT oft_mesh, i4
  CLASS(oft_mesh), INTENT(inout) :: self
  INTEGER(i4), INTENT(in) :: i
  END SUBROUTINE mesh_invert_cell
  !
  FUNCTION mesh_log2phys(self,cell,f) RESULT(pt)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: cell
  REAL(r8), INTENT(in) :: f(:)
  REAL(r8) :: pt(3)
  END FUNCTION mesh_log2phys
  !
  SUBROUTINE mesh_phys2log(self,i,pt,f)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: i
  REAL(r8), INTENT(in) :: pt(3)
  REAL(r8), INTENT(out) :: f(:)
  END SUBROUTINE mesh_phys2log
  !
  SUBROUTINE mesh_jacobian(self,cell,f,gop,j)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: cell
  REAL(r8), INTENT(in) :: f(:)
  REAL(r8), INTENT(out) :: gop(:,:)
  REAL(r8), INTENT(out) :: j
  END SUBROUTINE mesh_jacobian
  !
  SUBROUTINE mesh_hessian(self,i,f,g2op,K)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: i
  REAL(r8), INTENT(in) :: f(:)
  REAL(r8), INTENT(out) :: g2op(:,:)
  REAL(r8), INTENT(out) :: K(:,:)
  END SUBROUTINE mesh_hessian
  !
  SUBROUTINE mesh_snormal(self,i,ind,f,norm)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: i
  INTEGER(i4), INTENT(in) :: ind
  REAL(r8), INTENT(in) :: f(:)
  REAL(r8), INTENT(out) :: norm(3)
  END SUBROUTINE mesh_snormal
  !
  SUBROUTINE mesh_ctang(self,i,ind,f,tang)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: i
  INTEGER(i4), INTENT(in) :: ind
  REAL(r8), INTENT(in) :: f(:)
  REAL(r8), INTENT(out) :: tang(3)
  END SUBROUTINE mesh_ctang
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_get_surf_map
!
!> Get mapping between boundary and volume logical coordinates
!------------------------------------------------------------------------------
  SUBROUTINE mesh_get_surf_map(self,face,cell,lmap)
  IMPORT oft_mesh, i4
  CLASS(oft_mesh), INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: face !< Index of face on boundary mesh
  INTEGER(i4), INTENT(out) :: cell !< Cell containing face
  INTEGER(i4), INTENT(out) :: lmap(3) !< Coordinate mapping
  END SUBROUTINE mesh_get_surf_map
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_surf_to_vol
!
!> Map between surface and volume logical coordinates
!------------------------------------------------------------------------------
  SUBROUTINE mesh_surf_to_vol(self,fsurf,lmap,fvol)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self
  REAL(r8), INTENT(in) :: fsurf(:) !< Surface coordinates [3]
  INTEGER(i4), INTENT(in) :: lmap(3) !< Coordinate mapping
  REAL(r8), INTENT(out) :: fvol(:) !< Volume coordinates [4]
  END SUBROUTINE mesh_surf_to_vol
  !
  SUBROUTINE mesh_vlog(self,i,f)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: i
  REAL(r8), INTENT(out) :: f(:)
  END SUBROUTINE mesh_vlog
  !
  FUNCTION mesh_in_cell(self,f,tol) RESULT(eface)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self
  REAL(r8), INTENT(in) :: f(:)
  REAL(r8), INTENT(in) :: tol
  INTEGER(i4) :: eface
  END FUNCTION mesh_in_cell
  !
  SUBROUTINE mesh_quad_rule(self,order,quad_rule)
  IMPORT oft_mesh, oft_quad_type, i4
  CLASS(oft_mesh), INTENT(in) :: self
  INTEGER(i4), INTENT(in) :: order
  TYPE(oft_quad_type), INTENT(out) :: quad_rule
  END SUBROUTINE mesh_quad_rule
  !
  SUBROUTINE mesh_tessellate(self,rtmp,lctmp,order)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self
  REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp
  INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp
  INTEGER(i4), INTENT(in) :: order
  END SUBROUTINE mesh_tessellate
  !---------------------------------------------------------------------------
  ! FUNCTION: mesh_get_io_sizes
  !---------------------------------------------------------------------------
  !> Get variable sizes following tessellation
  !---------------------------------------------------------------------------
  function mesh_get_io_sizes(self) result(sizes)
  IMPORT oft_mesh, i4
  CLASS(oft_mesh), INTENT(in) :: self
  integer(i4) :: sizes(2)
  end function mesh_get_io_sizes
END INTERFACE
!---
REAL(r8), PUBLIC :: rgrnd(3)=(/1.d0,0.d0,0.d0/) !< Grounding point position
CLASS(oft_mesh), PUBLIC, POINTER :: mesh => NULL()
CLASS(oft_bmesh), PUBLIC, POINTER :: smesh => NULL()
INTEGER(i4), PRIVATE, PARAMETER :: ho_find_retry=20 !< Number of retry attempts during high order find_cell
#ifdef OFT_PLOT_DOUBLE
LOGICAL, PARAMETER :: PLOT_R4_FLAG=.FALSE.
#else
LOGICAL, PARAMETER :: PLOT_R4_FLAG=.TRUE.
#endif
CONTAINS
!------------------------------------------------------------------------------
!> Checks if a global mesh cell is curved or not
!!
!! @param[in] self Mesh containing cell
!! @param[in] cell Index of cell to check
!! @result (T/F) cell is curved?
!------------------------------------------------------------------------------
function cell_is_curved(self,cell) result(curved)
class(oft_amesh), intent(in) :: self
integer(i4), intent(in) :: cell
integer(i4) :: k,i
logical :: curved
DEBUG_STACK_PUSH
IF(self%type==3)THEN
  curved=.TRUE.
ELSE
  if(self%order==1)then
    curved=.FALSE.
  else
    curved=self%ho_info%is_curved(cell)
  end if
END IF
DEBUG_STACK_POP
end function cell_is_curved
!---------------------------------------------------------------------------
! FUNCTION: mesh_volume
!---------------------------------------------------------------------------
!> Estimate mesh volume
!---------------------------------------------------------------------------
function mesh_volume(self) result(vol)
class(oft_mesh), INTENT(IN) :: self
REAL(r8) :: vol
INTEGER(i4) :: i,j,ierr
REAL(r8) :: gop(3,4),vtmp
TYPE(oft_quad_type) :: quad
DEBUG_STACK_PUSH
!---Integerate over the volume
vol = 0.d0
CALL self%quad_rule(self%order*2, quad)
!$omp parallel do private(gop,vtmp) reduction(+:vol)
do i=1,self%nc ! Loop over cells
  DO j=1,quad%np ! Loop over quadrature points
    CALL self%jacobian(i,quad%pts(:,j),gop,vtmp)
    vol = vol + vtmp*quad%wts(j)
  END DO
end do
IF(.NOT.self%fullmesh)vol=oft_mpi_sum(vol)
CALL quad%delete()
DEBUG_STACK_POP
end function mesh_volume
!---------------------------------------------------------------------------
! SUBROUTINE: mesh_setup_io
!---------------------------------------------------------------------------
!> Estimate mesh volume
!---------------------------------------------------------------------------
SUBROUTINE mesh_setup_io(self,tess_order,basepath)
class(oft_mesh), INTENT(inout) :: self
integer(i4), intent(in) :: tess_order
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: basepath
integer(i4) :: i,j,k,id,error,two=2,io_sizes(2),dims(2)
integer(i4), POINTER :: lctmp(:,:),fmap(:)
real(r8), POINTER :: rtmp(:,:),reg_tmp(:)
DEBUG_STACK_PUSH
CALL hdf5_create_files(basepath=self%io_path)
!---Initialize HDF5 and open mesh file
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Writing mesh to plot files'
CALL oft_increase_indent
self%io_path=''
IF(PRESENT(basepath))THEN
  self%io_path=basepath
  call execute_command_line('mkdir -p '//TRIM(self%io_path), exitstat=error)
  IF(error/=0)CALL oft_abort('Failed to create output directory: '//TRIM(self%io_path), &
    "mesh_setup_io", __FILE__)
END IF
self%tess_order=tess_order
!---Get grid tessellation
CALL self%tessellate(rtmp, lctmp, self%tess_order)
!---Write out point list
CALL hdf5_write(rtmp,TRIM(self%io_path)//'mesh.'//hdf5_proc_str()//'.h5','R_vol',single_prec=PLOT_R4_FLAG)
deallocate(rtmp)
!---Write out cell point list
CALL hdf5_write(lctmp,TRIM(self%io_path)//'mesh.'//hdf5_proc_str()//'.h5','LC_vol')
dims=SHAPE(lctmp)
k=dims(2)/self%nc
deallocate(lctmp)
!---Create "dump.dat" file
CALL oft_hdf5_write_dump(self%type,self%get_io_sizes(),self%bmesh%get_io_sizes(),basepath=self%io_path)
!---Setup I/O for boundary mesh
self%bmesh%io_path=self%io_path
CALL self%bmesh%setup_io(self%tess_order,append_files=.TRUE.)
!---Write regions
ALLOCATE(reg_tmp(dims(2)))
!$omp parallel do private(j,id)
DO i=1,self%nc
  id=self%reg(i)
  DO j=1,k
  reg_tmp((i-1)*k+j)=REAL(id,8)
  END DO
END DO
CALL self%save_cell_scalar(reg_tmp,'REG_vol')
DEALLOCATE(reg_tmp)
!---Write surface IDs
IF(self%bmesh%nc>0)THEN
  ALLOCATE(fmap(self%nf))
  fmap=0
  !$omp parallel do
  do i=1,self%nbf
    fmap(self%lbf(i))=i
  end do
  io_sizes=self%bmesh%get_io_sizes()
  k=io_sizes(2)/self%bmesh%nc
  ALLOCATE(reg_tmp(io_sizes(2)))
  !$omp parallel do private(j,id)
  DO i=1,self%bmesh%nc
    id=self%bfs(fmap(ABS(self%bmesh%parent%lf(i))))
    DO j=1,k
      reg_tmp((i-1)*k+j)=REAL(id,8)
    END DO
  END DO
  CALL self%bmesh%save_cell_scalar(reg_tmp,'SID')
  DEALLOCATE(reg_tmp,fmap)
END IF
CALL oft_decrease_indent
DEBUG_STACK_POP
END SUBROUTINE mesh_setup_io
!---------------------------------------------------------------------------
! SUBROUTINE: mesh_save_vertex_scalar
!---------------------------------------------------------------------------
!> Write scalar vertex data out to file
!!
!! @param[in] p Vertex data
!! @param[in] tag Name of the output field
!---------------------------------------------------------------------------
subroutine mesh_save_vertex_scalar(self,p,path)
class(oft_mesh), INTENT(IN) :: self
real(r8), intent(in) :: p(:)
character(LEN=*), intent(in) :: path
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving scalar plot field: ',TRIM(path)
sizes=self%get_io_sizes()
IF(SIZE(p,DIM=1)/=sizes(1))CALL oft_abort("Incorrect array size","mesh_save_vertex_scalar",__FILE__)
CALL hdf5_write(p,TRIM(self%io_path)//"scalar_dump."//hdf5_proc_str()//".h5",path//hdf5_ts_str(),PLOT_R4_FLAG)
CALL oft_hdf5_add_dump(path,11,basepath=self%io_path)
DEBUG_STACK_POP
end subroutine mesh_save_vertex_scalar
!---------------------------------------------------------------------------
! SUBROUTINE: mesh_save_cell_scalar
!---------------------------------------------------------------------------
!> Write scalar cell data out to file
!!
!! @param[in] p Cell data
!! @param[in] tag Name of the output field
!---------------------------------------------------------------------------
subroutine mesh_save_cell_scalar(self,p,path)
class(oft_mesh), INTENT(IN) :: self
real(r8), intent(in) :: p(:)
character(LEN=*), intent(in) :: path
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving scalar plot field: ',TRIM(path)
sizes=self%get_io_sizes()
IF(SIZE(p,DIM=1)/=sizes(2))CALL oft_abort("Incorrect array size","mesh_save_cell_scalar",__FILE__)
CALL hdf5_write(p,TRIM(self%io_path)//"scalar_dump."//hdf5_proc_str()//".h5",path//hdf5_ts_str(),PLOT_R4_FLAG)
CALL oft_hdf5_add_dump(path,12,basepath=self%io_path)
DEBUG_STACK_POP
end subroutine mesh_save_cell_scalar
!---------------------------------------------------------------------------
! SUBROUTINE: mesh_save_vertex_vector
!---------------------------------------------------------------------------
!> Write vector vertex data out to file
!!
!! @param[in] bv Vertex data
!! @param[in] tag Name of the output field
!---------------------------------------------------------------------------
subroutine mesh_save_vertex_vector(self,bv,path)
class(oft_mesh), INTENT(IN) :: self
real(r8), intent(in) :: bv(:,:)
character(LEN=*), intent(in) :: path
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving vector plot field: ',TRIM(path)
IF(SIZE(bv,DIM=1)/=3)CALL oft_abort("Output array is not 3 vector","mesh_save_vertex_vector",__FILE__)
sizes=self%get_io_sizes()
IF(SIZE(bv,DIM=2)/=sizes(1))CALL oft_abort("Incorrect array size","mesh_save_vertex_vector",__FILE__)
CALL hdf5_write(bv,TRIM(self%io_path)//"vector_dump."//hdf5_proc_str()//".h5",path//hdf5_ts_str(),PLOT_R4_FLAG)
CALL oft_hdf5_add_dump(path,21,basepath=self%io_path)
DEBUG_STACK_POP
end subroutine mesh_save_vertex_vector
!---------------------------------------------------------------------------
! SUBROUTINE: mesh_save_cell_vector
!---------------------------------------------------------------------------
!> Write vector cell data out to file
!!
!! @param[in] bcc Cell data
!! @param[in] tag Name of the output field
!---------------------------------------------------------------------------
subroutine mesh_save_cell_vector(self,bcc,path)
class(oft_mesh), INTENT(IN) :: self
real(r8), intent(in) :: bcc(:,:)
character(LEN=*), intent(in) :: path
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving vector plot field: ',TRIM(path)
IF(SIZE(bcc,DIM=1)/=3)CALL oft_abort("Output array is not 3 vector","mesh_save_cell_vector",__FILE__)
sizes=self%get_io_sizes()
IF(SIZE(bcc,DIM=2)/=sizes(2))CALL oft_abort("Incorrect array size","mesh_save_cell_vector",__FILE__)
CALL hdf5_write(bcc,TRIM(self%io_path)//"vector_dump."//hdf5_proc_str()//".h5",path//hdf5_ts_str(),PLOT_R4_FLAG)
CALL oft_hdf5_add_dump(path,22,basepath=self%io_path)
DEBUG_STACK_POP
end subroutine mesh_save_cell_vector
!---------------------------------------------------------------------------
!> Destroy mesh object
!!
!! @note Should only be used via class \ref tri_mesh or children
!---------------------------------------------------------------------------
SUBROUTINE mesh_destroy(self)
CLASS(oft_mesh), INTENT(inout) :: self
INTEGER(i4) :: i
DEBUG_STACK_PUSH
CALL amesh_destroy(self)
self%nf = 0
!---Deallocate integer arrays
IF(ASSOCIATED(self%bf))DEALLOCATE(self%bf)
IF(ASSOCIATED(self%lbf))DEALLOCATE(self%lbf)
IF(ASSOCIATED(self%klef))DEALLOCATE(self%klef)
IF(ASSOCIATED(self%llef))DEALLOCATE(self%llef)
IF(ASSOCIATED(self%lfo))DEALLOCATE(self%lfo)
IF(ASSOCIATED(self%bfs))DEALLOCATE(self%bfs)
IF(ASSOCIATED(self%cell_fc))DEALLOCATE(self%cell_fc)
IF(ASSOCIATED(self%cell_fe))DEALLOCATE(self%cell_fe)
IF(ASSOCIATED(self%face_ed))DEALLOCATE(self%face_ed)
IF(ASSOCIATED(self%lfe))DEALLOCATE(self%lfe)
IF(ASSOCIATED(self%lfc))DEALLOCATE(self%lfc)
IF(ASSOCIATED(self%lcf))DEALLOCATE(self%lcf)
IF(ASSOCIATED(self%lcfo))DEALLOCATE(self%lcfo)
IF(ASSOCIATED(self%cv))DEALLOCATE(self%cv)
IF(ASSOCIATED(self%vv))DEALLOCATE(self%vv)
!---
CALL destory_seam(self%fstitch)
IF(ASSOCIATED(self%tloc_f))THEN
  DO i=1,self%nparts
    DEALLOCATE(self%tloc_f(i)%v)
  END DO
  DEALLOCATE(self%tloc_f)
END IF
!---
NULLIFY(self%bmesh)
DEBUG_STACK_POP
END SUBROUTINE mesh_destroy
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_findcell
!------------------------------------------------------------------------------
!> Find physical point in mesh.
!!
!! For high order grids an approximate guess is first computed using only the linear
!! portion of the mesh representation. This guess is then refined using the full mesh
!! representation and the @ref tetmesh_mapping::tetmesh_phys2logho "tetmesh_phys2logho"
!! subroutine. The maximum number of cell searches during this non-linear phase is controlled
!! by the module variable @ref tetmesh_mapping::ho_find_retry "ho_find_retry". For
!! more information see the documentation for @ref tetmesh_mapping::tetmesh_phys2logho
!! "tetmesh_phys2logho"
!!
!! @param[in] self Mesh to search
!! @param[in,out] cell Cell containing point on output, guess on input
!! @param[in] pt Coordinates to locate [3]
!! @param[out] fout Logical coordinates of point in cell (optional)
!------------------------------------------------------------------------------
subroutine mesh_findcell(self,cell,pt,fout)
class(oft_mesh), target, intent(inout) :: self
integer(i4), intent(inout) :: cell
real(r8), intent(in) :: pt(3)
real(r8), optional, intent(out) :: fout(4)
real(r8) :: d2,d2min,rcc(3),f(4),ftmp(4),fmin,fmax,tol=1.d-10
integer(i4) :: next,i,ii,minf
IF(SUM(pt**2)>1.d90)THEN
  cell=0
  IF(PRESENT(fout))fout=1.d2
  RETURN
END IF
DEBUG_STACK_PUSH
if((cell<=0).OR.(cell>self%nc))then
  !---Find closest point in mesh
  d2min=1.d99
  do i=1,self%np
    d2=SUM((pt-self%r(:,i))**2)
    if(d2<d2min)then
      d2min=d2
      next=i
    end if
  end do
  !---Search surounding cells for closest
  d2min=1.d99
  f=0.d0
  DO i=1,self%cell_np
    CALL self%vlog(i,ftmp)
    f=f+ftmp
  END DO
  f=f/REAL(self%cell_np,8)
  do ii=self%kpc(next),self%kpc(next+1)-1
    i=self%lpc(ii)
    rcc=self%log2phys(i,f)
    d2=SUM((pt-rcc)**2)
    if(d2<d2min)then
      d2min=d2
      cell=i
    end if
  end do
end if
next=cell
do i=1,self%nc
  cell=next
  CALL self%phys2log(cell,pt,f)
  minf=self%in_cell(f, tol)
  IF(minf==0)EXIT
  next=self%lcc(minf,cell)
  IF(next==0)EXIT ! pt off mesh but closest to cell
end do
if(i>=self%nc)then
  cell=0
  f=-1.d0
end if
if(PRESENT(fout))fout=f
DEBUG_STACK_POP
end subroutine mesh_findcell
!------------------------------------------------------------------------------
! SUBROUTINE: mesh_findcell2
!------------------------------------------------------------------------------
!> Find physical point in mesh using robust method
!!
!! First, a collection of the `nclosest` mesh vertices is found. Then the surrounding
!! cells are searched for the specified point. If the point is not found the
!! standard @ref tetmesh_mapping::tetmesh_findcell "tetmesh_findcell" subroutine is
!! used as a fallback.
!!
!! @param[in] self Mesh to search
!! @param[in,out] cell Cell containing point on output, guess on input
!! @param[in] pt Coordinates to locate [3]
!! @param[in] nclosest Number of candidate vertices to use for search
!! @param[out] fout Logical coordinates of point in cell (optional)
!------------------------------------------------------------------------------
subroutine mesh_findcell2(self,cell,pt,nclosest,fout)
class(oft_mesh), target, intent(inout) :: self
integer(i4), intent(inout) :: cell
real(r8), intent(in) :: pt(3)
integer(i4), intent(in) :: nclosest
real(r8), optional, intent(out) :: fout(4)
real(r8) :: d2,f(4),fmin,fmax,tol=1.d-10,d2mins(nclosest)
integer(i4) :: i,ii,k,imins(nclosest),minf
IF(SUM(pt**2)>1.d90)THEN
  cell=0
  IF(PRESENT(fout))fout=1.d2
  RETURN
END IF
DEBUG_STACK_PUSH
if((cell<=0).OR.(cell>self%nc))then
   d2mins=1.d99
   imins=0
   do i=1,self%np
     d2=SUM((pt-self%r(:,i))**2)
     do ii=1,nclosest
       if(d2<d2mins(ii))then
         if(ii<nclosest)then
           d2mins(ii+1:nclosest)=d2mins(ii:nclosest-1)
           imins(ii+1:nclosest)=imins(ii:nclosest-1)
         endif
         d2mins(ii)=d2
         imins(ii)=i
         exit
       end if
     end do
   end do
   cell=0
   do k=1,nclosest
     if(imins(k)==0)exit
     do ii=self%kpc(imins(k)),self%kpc(imins(k)+1)-1
       i=self%lpc(ii)
       CALL self%phys2log(i,pt,f)
       minf=self%in_cell(f, tol)
       IF(minf==0)THEN
         cell=i
         EXIT
       END IF
     end do
     if(cell>0)exit
   end do
   ! Fallback to old method
   if(cell==0)call mesh_findcell(self,cell,pt,f)
else
  ! Use old method
  call mesh_findcell(self,cell,pt,f)
end if
if(PRESENT(fout))fout=f
DEBUG_STACK_POP
end subroutine mesh_findcell2
!---------------------------------------------------------------------------
!> Estimate mesh area
!---------------------------------------------------------------------------
function bmesh_area(self) result(area)
class(oft_bmesh), intent(in) :: self
REAL(r8) :: area
INTEGER(i4) :: i,j,ierr
REAL(r8) :: gop(3,3),vtmp
TYPE(oft_quad_type) :: quad
DEBUG_STACK_PUSH
!---Integerate over the surface
area = 0.d0
CALL self%quad_rule(self%order*2, quad)
!$omp parallel do private(gop,vtmp) reduction(+:area)
do i=1,self%nc
  DO j=1,quad%np
    call self%jacobian(i,quad%pts(:,j),gop,vtmp)
    area = area + vtmp*quad%wts(j)
  END DO
end do
IF(.NOT.self%fullmesh)area=oft_mpi_sum(area)
CALL quad%delete()
DEBUG_STACK_POP
end function bmesh_area
!------------------------------------------------------------------------------
!> Find physical point in mesh.
!!
!! @warning Only works for 2D meshes and uses linear logical mapping
!!
!! @param[in] self Mesh to search
!! @param[in,out] cell Cell containing point on output, guess on input
!! @param[in] pt Coordinates to locate [2]
!! @param[out] fout Logical coordinates of point in cell (optional)
!------------------------------------------------------------------------------
subroutine bmesh_findcell(self,cell,pt,fout)
class(oft_bmesh), intent(in) :: self
integer(i4), intent(inout) :: cell
real(r8), intent(in) :: pt(2)
real(r8), optional, intent(out) :: fout(3)
real(r8) :: d2,d2min,rcc(2),f(3),fmin,fmax,tol=1.d-10,pttmp(3)
integer(i4) :: next,i,ii,mine
IF(self%dim/=2)CALL oft_abort("Only supported for dim=2","bmesh_findcell",__FILE__)
IF(self%type/=1)CALL oft_abort("Only supported for type=1","bmesh_findcell",__FILE__)
pttmp(1:2)=pt
pttmp(3)=0.d0
d2min=1.d99
if(( cell<=0 ).OR.( cell>self%nc ))then
  do ii=1,self%nbc
    i=self%lbc(ii)
    rcc=(self%r(1:2,self%lc(1,i)) &
      + self%r(1:2,self%lc(2,i)) &
      + self%r(1:2,self%lc(3,i)))
    d2=sum((3.d0*pt-rcc)**2)
    if( d2<d2min )then
      d2min=d2
      cell=i
    end if
  end do
end if
next=cell
do i=1,self%nc
  cell=next
  CALL self%phys2log(cell,pttmp,f)
  fmin=minval(f)
  fmax=maxval(f)
  if (( fmax<=1.d0+tol ).AND.( fmin>=-tol ))exit
  mine=minloc(f,DIM=1)
  next=self%lcc(mine,cell)
  if(next==0)exit ! pt off mesh but closest to cell
end do
if(present(fout))fout=f
if(i>=self%nc)then
  cell=0
  fout=-1.d0
endif
end subroutine bmesh_findcell
!---------------------------------------------------------------------------
!> Destroy trimesh object
!!
!! @note Should only be used via class \ref tri_mesh or children
!---------------------------------------------------------------------------
SUBROUTINE bmesh_setup_io(self,tess_order,append_files,basepath)
CLASS(oft_bmesh), INTENT(inout) :: self
integer(i4), intent(in) :: tess_order
logical, optional, intent(in) :: append_files
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: basepath
logical :: create_files
integer(i4) :: i,k,j,id,error,dims(2)
integer(i4), POINTER :: lftmp(:,:),fmap(:)
real(r8), POINTER :: rtmp(:,:),reg_tmp(:)
DEBUG_STACK_PUSH
!---Setup I/O
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Writing boundary mesh to plot files'
CALL oft_increase_indent
self%io_path=''
IF(PRESENT(basepath))THEN
  self%io_path=basepath
  call execute_command_line('mkdir -p '//TRIM(self%io_path), exitstat=error)
  IF(error/=0)CALL oft_abort('Failed to create output directory: '//TRIM(self%io_path), &
    "bmesh_setup_io", __FILE__)
END IF
self%tess_order=tess_order
create_files=.TRUE.
IF(PRESENT(append_files))create_files=(.NOT.append_files)
IF(create_files)THEN
  CALL hdf5_create_files(basepath=self%io_path)
  CALL oft_hdf5_write_dump(self%type,[0,0],self%get_io_sizes(),basepath=self%io_path)
END IF
IF(self%nc==0)THEN
  CALL oft_decrease_indent
  DEBUG_STACK_POP
  RETURN
END IF
!---Get grid tessellation
CALL self%tessellate(rtmp, lftmp, self%tess_order)
!---Write out point list
CALL hdf5_write(rtmp,TRIM(self%io_path)//'mesh.'//hdf5_proc_str()//'.h5', "R_surf",single_prec=PLOT_R4_FLAG)
deallocate(rtmp)
!---Write out cell point list
CALL hdf5_write(lftmp,TRIM(self%io_path)//'mesh.'//hdf5_proc_str()//'.h5', "LC_surf")
dims=SHAPE(lftmp)
k=INT(dims(2),4)/self%nc
deallocate(lftmp)
!---Write regions
ALLOCATE(reg_tmp(dims(2)))
IF(ASSOCIATED(self%reg))THEN
  !$omp parallel do private(j,id)
  DO i=1,self%nc
    id=self%reg(i)
    DO j=1,k
      reg_tmp((i-1)*k+j)=REAL(id,8)
    END DO
  END DO
ELSE
  reg_tmp=0.d0
END IF
CALL self%save_cell_scalar(reg_tmp,'REG_surf')
DEALLOCATE(reg_tmp)
CALL oft_decrease_indent
DEBUG_STACK_POP
END SUBROUTINE bmesh_setup_io
!---------------------------------------------------------------------------
!> Write scalar vertex data out to file
!!
!! @param[in] p Vertex data
!! @param[in] tag Name of the output field
!---------------------------------------------------------------------------
subroutine bmesh_save_vertex_scalar(self,p,path)
class(oft_bmesh), INTENT(IN) :: self
real(r8), intent(in) :: p(:)
character(LEN=*), intent(in) :: path
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving scalar plot field: ',TRIM(path)
sizes=self%get_io_sizes()
IF(SIZE(p,DIM=1)/=sizes(1))CALL oft_abort("Incorrect array size","bmesh_save_vertex_scalar",__FILE__)
CALL hdf5_write(p,TRIM(self%io_path)//"scalar_dump."//hdf5_proc_str()//".h5",path//hdf5_ts_str(),PLOT_R4_FLAG)
CALL oft_hdf5_add_dump(path,31,basepath=self%io_path)
DEBUG_STACK_POP
end subroutine bmesh_save_vertex_scalar
!---------------------------------------------------------------------------
!> Write scalar cell data out to file
!!
!! @param[in] p Cell data
!! @param[in] tag Name of the output field
!---------------------------------------------------------------------------
subroutine bmesh_save_cell_scalar(self,p,path)
class(oft_bmesh), INTENT(IN) :: self
real(r8), intent(in) :: p(:)
character(LEN=*), intent(in) :: path
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving scalar plot field: ',TRIM(path)
sizes=self%get_io_sizes()
IF(SIZE(p,DIM=1)/=sizes(2))CALL oft_abort("Incorrect array size","bmesh_save_cell_scalar",__FILE__)
CALL hdf5_write(p,TRIM(self%io_path)//"scalar_dump."//hdf5_proc_str()//".h5",path//hdf5_ts_str(),PLOT_R4_FLAG)
CALL oft_hdf5_add_dump(path,32,basepath=self%io_path)
DEBUG_STACK_POP
end subroutine bmesh_save_cell_scalar
!---------------------------------------------------------------------------
!> Write vector vertex data out to file
!!
!! @param[in] bv Vertex data
!! @param[in] tag Name of the output field
!---------------------------------------------------------------------------
subroutine bmesh_save_vertex_vector(self,bv,path)
class(oft_bmesh), INTENT(IN) :: self
real(r8), intent(in) :: bv(:,:)
character(LEN=*), intent(in) :: path
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving vector plot field: ',TRIM(path)
IF(SIZE(bv,DIM=1)/=3)CALL oft_abort("Output array is not 3 vector","bmesh_save_vertex_vector",__FILE__)
sizes=self%get_io_sizes()
IF(SIZE(bv,DIM=2)/=sizes(1))CALL oft_abort("Incorrect array size","bmesh_save_vertex_vector",__FILE__)
CALL hdf5_write(bv,TRIM(self%io_path)//"vector_dump."//hdf5_proc_str()//".h5",path//hdf5_ts_str(),PLOT_R4_FLAG)
CALL oft_hdf5_add_dump(path,41,basepath=self%io_path)
DEBUG_STACK_POP
end subroutine bmesh_save_vertex_vector
!---------------------------------------------------------------------------
!> Write vector cell data out to file
!!
!! @param[in] bcc Cell data
!! @param[in] tag Name of the output field
!---------------------------------------------------------------------------
subroutine bmesh_save_cell_vector(self,bcc,path)
class(oft_bmesh), INTENT(IN) :: self
real(r8), intent(in) :: bcc(:,:)
character(LEN=*), intent(in) :: path
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving vector plot field: ',TRIM(path)
IF(SIZE(bcc,DIM=1)/=3)CALL oft_abort("Output array is not 3 vector","bmesh_save_cell_vector",__FILE__)
sizes=self%get_io_sizes()
IF(SIZE(bcc,DIM=2)/=sizes(2))CALL oft_abort("Incorrect array size","bmesh_save_cell_vector",__FILE__)
CALL hdf5_write(bcc,TRIM(self%io_path)//"vector_dump."//hdf5_proc_str()//".h5",path//hdf5_ts_str(),PLOT_R4_FLAG)
CALL oft_hdf5_add_dump(path,42,basepath=self%io_path)
DEBUG_STACK_POP
end subroutine bmesh_save_cell_vector
!---------------------------------------------------------------------------
!> Destroy surface mesh object
!!
!! @note Should only be used via class \ref tri_mesh or children
!---------------------------------------------------------------------------
SUBROUTINE bmesh_destroy(self)
CLASS(oft_bmesh), INTENT(inout) :: self
DEBUG_STACK_PUSH
CALL amesh_destroy(self)
!---Deallocate class-specific arrays
IF(ASSOCIATED(self%lco))DEALLOCATE(self%lco)
IF(ASSOCIATED(self%bes))DEALLOCATE(self%bes)
IF(ASSOCIATED(self%ca))DEALLOCATE(self%ca)
IF(ASSOCIATED(self%va))DEALLOCATE(self%va)
!---Deallocated parent mesh info
IF(ASSOCIATED(self%parent))THEN
  IF(ASSOCIATED(self%parent%lp))DEALLOCATE(self%parent%lp)
  IF(ASSOCIATED(self%parent%le))DEALLOCATE(self%parent%le)
  IF(ASSOCIATED(self%parent%lf))DEALLOCATE(self%parent%lf)
END IF
DEBUG_STACK_POP
END SUBROUTINE bmesh_destroy
!---------------------------------------------------------------------------
!> Destroy mesh object
!!
!! @note Should only be used via class \ref tri_mesh or children
!---------------------------------------------------------------------------
SUBROUTINE amesh_destroy(self)
CLASS(oft_amesh), INTENT(inout) :: self
INTEGER(i4) :: i
DEBUG_STACK_PUSH
!---Reset a few counters
self%np = 0
self%ne = 0
self%nc = 0
!---Deallocate logical arrays
IF(ASSOCIATED(self%bp))DEALLOCATE(self%bp)
IF(ASSOCIATED(self%be))DEALLOCATE(self%be)
IF(ASSOCIATED(self%bc))DEALLOCATE(self%bc)
IF(ASSOCIATED(self%cp))DEALLOCATE(self%cp)
IF(ASSOCIATED(self%ce))DEALLOCATE(self%ce)
!---Deallocate integer arrays
IF(ASSOCIATED(self%reg))DEALLOCATE(self%reg)
IF(ASSOCIATED(self%kpc))DEALLOCATE(self%kpc)
IF(ASSOCIATED(self%lpc))DEALLOCATE(self%lpc)
IF(ASSOCIATED(self%kec))DEALLOCATE(self%kec)
IF(ASSOCIATED(self%lec))DEALLOCATE(self%lec)
IF(ASSOCIATED(self%kpp))DEALLOCATE(self%kpp)
IF(ASSOCIATED(self%lpp))DEALLOCATE(self%lpp)
IF(ASSOCIATED(self%kee))DEALLOCATE(self%kee)
IF(ASSOCIATED(self%lee))DEALLOCATE(self%lee)
IF(ASSOCIATED(self%lbp))DEALLOCATE(self%lbp)
IF(ASSOCIATED(self%lbe))DEALLOCATE(self%lbe)
IF(ASSOCIATED(self%lbc))DEALLOCATE(self%lbc)
IF(ASSOCIATED(self%kpe))DEALLOCATE(self%kpe)
IF(ASSOCIATED(self%lpe))DEALLOCATE(self%lpe)
IF(ASSOCIATED(self%klpe))DEALLOCATE(self%klpe)
IF(ASSOCIATED(self%llpe))DEALLOCATE(self%llpe)
IF(ASSOCIATED(self%cell_ed))DEALLOCATE(self%cell_ed)
IF(ASSOCIATED(self%le))DEALLOCATE(self%le)
IF(ASSOCIATED(self%lc))DEALLOCATE(self%lc)
IF(ASSOCIATED(self%lce))DEALLOCATE(self%lce)
IF(ASSOCIATED(self%lcc))DEALLOCATE(self%lcc)
!---Deallocate real arrays
IF(ASSOCIATED(self%r))DEALLOCATE(self%r)
!---Destory thread partitioning pointers
IF(ASSOCIATED(self%tloc_p))THEN
  DO i=1,self%nparts
    DEALLOCATE(self%tloc_p(i)%v)
    DEALLOCATE(self%tloc_e(i)%v)
    DEALLOCATE(self%tloc_c(i)%v)
  END DO
  DEALLOCATE(self%tloc_p,self%tloc_e,self%tloc_c)
END IF
!---
NULLIFY(self%pstitch%be,self%pstitch%lbe)
CALL destory_seam(self%pstitch)
NULLIFY(self%estitch%be,self%estitch%lbe)
CALL destory_seam(self%estitch)
!---
IF(ASSOCIATED(self%global%lp))DEALLOCATE(self%global%lp)
IF(ASSOCIATED(self%global%le))DEALLOCATE(self%global%le)
IF(ASSOCIATED(self%global%lf))DEALLOCATE(self%global%lf)
IF(ASSOCIATED(self%global%lc))DEALLOCATE(self%global%lc)
IF(ASSOCIATED(self%global%gbp))DEALLOCATE(self%global%gbp)
IF(ASSOCIATED(self%global%gbe))DEALLOCATE(self%global%gbe)
IF(ASSOCIATED(self%global%gbf))DEALLOCATE(self%global%gbf)
IF(ASSOCIATED(self%global%gbc))DEALLOCATE(self%global%gbc)
!---
IF(ASSOCIATED(self%base%lp))DEALLOCATE(self%base%lp)
IF(ASSOCIATED(self%base%le))DEALLOCATE(self%base%le)
IF(ASSOCIATED(self%base%lf))DEALLOCATE(self%base%lf)
IF(ASSOCIATED(self%base%lc))DEALLOCATE(self%base%lc)
IF(ASSOCIATED(self%base%lcpart))DEALLOCATE(self%base%lcpart)
!---
IF(ASSOCIATED(self%periodic%lp))DEALLOCATE(self%periodic%lp)
IF(ASSOCIATED(self%periodic%le))DEALLOCATE(self%periodic%le)
IF(ASSOCIATED(self%periodic%lf))DEALLOCATE(self%periodic%lf)
!---
IF(ASSOCIATED(self%ho_info%is_curved))DEALLOCATE(self%ho_info%is_curved)
IF(ASSOCIATED(self%ho_info%lep))DEALLOCATE(self%ho_info%lep)
IF(ASSOCIATED(self%ho_info%lfp))DEALLOCATE(self%ho_info%lfp)
IF(ASSOCIATED(self%ho_info%lcp))DEALLOCATE(self%ho_info%lcp)
IF(ASSOCIATED(self%ho_info%r))DEALLOCATE(self%ho_info%r)
DEBUG_STACK_POP
END SUBROUTINE amesh_destroy
END MODULE oft_mesh_type
