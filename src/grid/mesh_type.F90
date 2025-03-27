!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------
MODULE oft_mesh_type
USE oft_base
USE oft_stitching, ONLY: oft_seam, destory_seam
USE oft_io, ONLY: hdf5_write, hdf5_field_exist, &
  hdf5_create_group, xdmf_plot_file
USE oft_quadrature
IMPLICIT NONE
#include "local.h"
PRIVATE
PUBLIC cell_is_curved, oft_init_seam, mesh_findcell, mesh_findcell2, bmesh_findcell
!---------------------------------------------------------------------------------
!> Global mesh information and indicies
!!
!! Contains global mexh context information.
!! - Global counts
!! - Global element indices
!! - Global boundary flags
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: mesh_per
  INTEGER(i4) :: nper = 0 !< Number of periodic directions
  INTEGER(i4), POINTER, DIMENSION(:) :: lp => NULL() !< Global index of points (np)
  INTEGER(i4), POINTER, DIMENSION(:) :: le => NULL() !< Global index of edges (ne)
  INTEGER(i4), POINTER, DIMENSION(:) :: lf => NULL() !< Global index of faces (nf)
END TYPE mesh_per
!---------------------------------------------------------------------------------
!> MPI global index information (For I/O Only)
!!
!! Contains global indices for HDF5 I/O.
!! - Local starting and closing index
!! - Maximum counts
!---------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------
!> Processor-processor connectivity information for mesh
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: mesh_seam
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
END TYPE mesh_seam
!---------------------------------------------------------------------------------
!> Global mesh information and indicies
!!
!! Contains global mexh context information.
!! - Global counts
!! - Global element indices
!! - Global boundary flags
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: mesh_global
  INTEGER(i8) :: np = 0 !< Global point count
  INTEGER(i8) :: ne = 0 !< Global edge count
  INTEGER(i8) :: nf = 0 !< Global face count
  INTEGER(i8) :: nc = 0 !< Global cell count
  INTEGER(i8) :: nbp = 0 !< Global boundary point count
  INTEGER(i8) :: nbe = 0 !< Global boundary edge count
  INTEGER(i8) :: nbf = 0 !< Global boundary face count
  INTEGER(i8) :: nbc = 0 !< Global boundary cell count
  LOGICAL, POINTER, DIMENSION(:) :: gbp => NULL() !< Global boundary point flag (np)
  LOGICAL, POINTER, DIMENSION(:) :: gbe => NULL() !< Global boundary edge flag (ne)
  LOGICAL, POINTER, DIMENSION(:) :: gbf => NULL() !< Global boundary face flag (nf)
  LOGICAL, POINTER, DIMENSION(:) :: gbc => NULL() !< Global boundary cell flag (nc)
  INTEGER(i8), POINTER, DIMENSION(:) :: lp => NULL() !< Global index of points (np)
  INTEGER(i8), POINTER, DIMENSION(:) :: le => NULL() !< Global index of edges (ne) [oriented]
  INTEGER(i8), POINTER, DIMENSION(:) :: lf => NULL() !< Global index of faces (nf)
  INTEGER(i8), POINTER, DIMENSION(:) :: lc => NULL() !< Global index of cells (nc)
  TYPE(mesh_seam), POINTER :: seam => NULL() !< Global domain-domain connectivity information
END TYPE mesh_global
!---------------------------------------------------------------------------------
!> Base mesh information and indicies
!!
!! Contains global mexh context information.
!! - Global counts
!! - Global element indices
!! - Global boundary flags
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: mesh_base
  INTEGER(i4) :: np = 0 !< Global point count
  INTEGER(i4) :: nc = 0 !< Global cell count
  INTEGER(i4), POINTER, DIMENSION(:) :: lp => NULL() !< Global index of points (np)
  INTEGER(i4), POINTER, DIMENSION(:) :: le => NULL() !< Global index of points (np)
  INTEGER(i4), POINTER, DIMENSION(:) :: lf => NULL() !< Global index of points (np)
  INTEGER(i4), POINTER, DIMENSION(:) :: lc => NULL() !< Global index of cells (nc)
  INTEGER(i4), POINTER, DIMENSION(:) :: lcpart => NULL() !< Global index of cells (nc)
END TYPE mesh_base
!---------------------------------------------------------------------------------
!> High order tet geometry information
!!
!! Contains additional data for high order tetrahedra.
!! - Number of additional nodes per geometry primative
!! - List of high-order node locations
!! - List of node points
!---------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------
!> Abstrac mesh type (surface or volume)
!!
!! Contains geometry information for the computational grid.
!! - Entity counts
!! - Mesh type and order
!! - Global mesh information
!! - Linkage of geometric primatives
!---------------------------------------------------------------------------------
TYPE, PUBLIC, ABSTRACT :: oft_amesh
  LOGICAL :: fullmesh = .TRUE. !< Local mesh flag (False if distributed)
  INTEGER(i4) :: cad_type = 1 !< Type of CAD geometry
  INTEGER(i4) :: type = 0 !< Mesh type
  INTEGER(i4) :: order = 1 !< order of boundary tets (default=linear)
  INTEGER(i4) :: tess_order = 0 !< order of boundary tets (default=linear)
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
  TYPE(oft_seam) :: pstitch !< Needs docs
  TYPE(oft_seam) :: estitch !< Needs docs
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: tloc_p => NULL() !< Point thread ownership
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: tloc_e => NULL() !< Edge thread ownership
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: tloc_c => NULL() !< Cell thread ownership
  TYPE(mesh_per) :: periodic !< Periodic information
  TYPE(ho_mesh) :: ho_info !< High order geometry information
END TYPE oft_amesh
!---------------------------------------------------------------------------------
!> Global mesh information and indicies
!!
!! Contains global mexh context information.
!! - Global counts
!! - Global element indices
!! - Global boundary flags
!---------------------------------------------------------------------------------
TYPE, PUBLIC :: bmesh_parent
  INTEGER(i4) :: np = 0 !< Global point count
  INTEGER(i4) :: ne = 0 !< Global edge count
  INTEGER(i4) :: nf = 0 !< Global face count
  INTEGER(i4), POINTER, DIMENSION(:) :: lp => NULL() !< Global index of points (np)
  INTEGER(i4), POINTER, DIMENSION(:) :: le => NULL() !< Global index of edges (ne) [oriented]
  INTEGER(i4), POINTER, DIMENSION(:) :: lf => NULL() !< Global index of faces (nf)
END TYPE bmesh_parent
!---------------------------------------------------------------------------------
!> Surface mesh type
!!
!! Contains geometry information for the computational grid.
!! - Entity counts
!! - Mesh type and order
!! - Global mesh information
!! - Linkage of geometric primatives
!---------------------------------------------------------------------------------
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
  PROCEDURE(bmesh_invert_cell), DEFERRED :: invert_face
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
  PROCEDURE(bmesh_tessellated_sizes), DEFERRED :: tessellated_sizes
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
  !------------------------------------------------------------------------------
  !> Setup mesh with implementation specifics (`cell_np`, `cell_ne`, etc.)
  !------------------------------------------------------------------------------
  SUBROUTINE bmesh_setup(self,cad_type,has_parent)
  IMPORT oft_bmesh, i4
  CLASS(oft_bmesh), INTENT(inout) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cad_type !< CAD/mesh interface ID number
  LOGICAL, INTENT(in) :: has_parent !< Is this mesh the/a surface of a volume mesh?
  END SUBROUTINE bmesh_setup
  !------------------------------------------------------------------------------
  !> Save mesh to transfer file
  !------------------------------------------------------------------------------
  SUBROUTINE bmesh_save(self,filename)
  IMPORT oft_bmesh
  CLASS(oft_bmesh), INTENT(in) :: self !< Mesh object
  CHARACTER(LEN=*), INTENT(in) :: filename !< File to save mesh to
  END SUBROUTINE bmesh_save
  !------------------------------------------------------------------------------
  !> Load mesh from transfer file
  !------------------------------------------------------------------------------
  SUBROUTINE bmesh_load(self,filename)
  IMPORT oft_bmesh
  CLASS(oft_bmesh), INTENT(inout) :: self !< Mesh object
  CHARACTER(LEN=*), INTENT(in) :: filename !< File to load mesh from
  END SUBROUTINE bmesh_load
  !------------------------------------------------------------------------------
  !> Set maximum order of spatial mapping
  !------------------------------------------------------------------------------
  SUBROUTINE bmesh_set_order(self,order)
  IMPORT oft_bmesh, i4
  CLASS(oft_bmesh), INTENT(inout) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: order !< Maximum order of spatial mapping
  END SUBROUTINE bmesh_set_order
  !------------------------------------------------------------------------------
  !> Turn cell "inside out", used to ensure consistent orientations
  !------------------------------------------------------------------------------
  SUBROUTINE bmesh_invert_cell(self,cell)
  IMPORT oft_bmesh, i4
  CLASS(oft_bmesh), INTENT(inout) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cell !< Index of cell to invert
  END SUBROUTINE bmesh_invert_cell
  !------------------------------------------------------------------------------
  !> Map from logical to physical coordinates in a given cell
  !------------------------------------------------------------------------------
  FUNCTION bmesh_log2phys(self,cell,f) RESULT(pt)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cell !< Index of cell for evaulation
  REAL(r8), INTENT(in) :: f(:) !< Logical coordinate in cell [4]
  REAL(r8) :: pt(3) !< Physical position [3]
  END FUNCTION bmesh_log2phys
  !------------------------------------------------------------------------------
  !> Map from physical to logical coordinates in a given cell
  !------------------------------------------------------------------------------
  SUBROUTINE bmesh_phys2log(self,cell,pt,f)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cell !< Index of cell for evaulation
  REAL(r8), INTENT(in) :: pt(3) !< Physical position [3]
  REAL(r8), INTENT(out) :: f(:) !< Logical coordinates within the cell [4]
  END SUBROUTINE bmesh_phys2log
  !------------------------------------------------------------------------------
  !> Compute the spatial jacobian matrix and its determinant for a given cell at a given logical position
  !------------------------------------------------------------------------------
  SUBROUTINE bmesh_jacobian(self,cell,f,gop,j)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cell !< Index of cell for evaulation
  REAL(r8), INTENT(in) :: f(:) !< Logical coordinate in cell [3]
  REAL(r8), INTENT(out) :: gop(:,:) !< Jacobian matrix \f$ (\frac{\partial x_i}{\partial \lambda_j})^{-1} \f$ [3,4]
  REAL(r8), INTENT(out) :: j !< Jacobian of transformation from logical to physical coordinates
  END SUBROUTINE bmesh_jacobian
  !------------------------------------------------------------------------------
  !> Compute the spatial hessian matrices for a given cell at a given logical position
  !------------------------------------------------------------------------------
  SUBROUTINE bmesh_hessian(self,cell,f,g2op,K)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cell !< Index of cell for evaulation
  REAL(r8), INTENT(in) :: f(:) !< Logical coordinate in cell [4]
  REAL(r8), INTENT(out) :: g2op(:,:) !< Second order Jacobian matrix
  !! \f$ (\frac{\partial x_i}{\partial \lambda_l} \frac{\partial x_j}{\partial \lambda_k})^{-1} \f$
  REAL(r8), INTENT(out) :: K(:,:) !< Gradient correction matrix
  !! \f$ \frac{\partial^2 x_i}{\partial \lambda_k \partial \lambda_l}\f$ [10,3]
  END SUBROUTINE bmesh_hessian
  !------------------------------------------------------------------------------
  !> Get unit normal for surface at a given point in a given cell
  !------------------------------------------------------------------------------
  SUBROUTINE bmesh_norm(self,cell,f,n)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), TARGET, INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cell !< Cell containing point
  REAL(r8), INTENT(in) :: f(:) !< Logical coordinates in cell
  REAL(r8), INTENT(out) :: n(3) !< Unit normal [3]
  END SUBROUTINE bmesh_norm
  !------------------------------------------------------------------------------
  !> Get tangent basis set for surface at a given point in a given cell
  !------------------------------------------------------------------------------
  SUBROUTINE bmesh_tang(self,cell,f,t)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), TARGET, INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cell !< Cell containing point
  REAL(r8), INTENT(in) :: f(:) !< Logical coordinates in cell
  REAL(r8), INTENT(out) :: t(3,2) !< Unit tangent basis set [3,2]
  END SUBROUTINE bmesh_tang
  !------------------------------------------------------------------------------
  !> Retrieve suitable quadrature rule for mesh with given order
  !------------------------------------------------------------------------------
  SUBROUTINE bmesh_quad_rule(self,order,quad_rule)
  IMPORT oft_bmesh, oft_quad_type, i4
  CLASS(oft_bmesh), INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: order !< Desired order of quadrature rule
  TYPE(oft_quad_type), INTENT(out) :: quad_rule !< Resulting quadrature rule
  END SUBROUTINE bmesh_quad_rule
  !------------------------------------------------------------------------------
  !> Get position in logical space of vertex `i`
  !------------------------------------------------------------------------------
  SUBROUTINE bmesh_vlog(self,i,f)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: i !< Vertex to locate
  REAL(r8), INTENT(out) :: f(:) !< Logical coordinates of vertex `i`
  END SUBROUTINE bmesh_vlog
  !------------------------------------------------------------------------------
  !> Test if logical position lies within the base cell
  !!
  !! @returns Position `f` is inside the base cell?
  !------------------------------------------------------------------------------
  FUNCTION bmesh_in_cell(self,f,tol) RESULT(eedge)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), INTENT(in) :: self !< Mesh object
  REAL(r8), INTENT(in) :: f(:) !< Logical coordinate to evaluate
  REAL(r8), INTENT(in) :: tol !< Tolerance for test
  INTEGER(i4) :: eedge
  END FUNCTION bmesh_in_cell
  !------------------------------------------------------------------------------
  !> Tessellate mesh onto lagrange FE nodes of specified order (usually for plotting)
  !!
  !! @note The maximum tessellation order currently supported is 4
  !! (may be lower for certain mesh types).
  !!
  !! @warning Cell lists are returned with zero based indexing
  !------------------------------------------------------------------------------
  SUBROUTINE bmesh_tessellate(self,rtmp,lctmp,order)
  IMPORT oft_bmesh, i4, r8
  CLASS(oft_bmesh), INTENT(in) :: self !< Mesh object
  REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp !< Vertices for tessellation [3,:]
  INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp !< Cell list for tessellation [self%cell_np,:]
  INTEGER(i4), INTENT(in) :: order !< Tessellation order
  END SUBROUTINE bmesh_tessellate
  !------------------------------------------------------------------------------
  !> Get variable sizes following tessellation
  !------------------------------------------------------------------------------
  function bmesh_tessellated_sizes(self) result(sizes)
  IMPORT oft_bmesh, i4
  CLASS(oft_bmesh), INTENT(in) :: self !< Mesh object
  integer(i4) :: sizes(2) !< Array sizes following tessellation [np_tess,nc_tess]
  end function bmesh_tessellated_sizes
END INTERFACE
!---------------------------------------------------------------------------------
!> Tetrahedral Mesh type
!!
!! Contains geometry information for the computational grid.
!! - Entity counts
!! - Mesh type and order
!! - Global mesh information
!! - Linkage of geometric primatives
!---------------------------------------------------------------------------------
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
  TYPE(oft_seam) :: fstitch !< Needs docs
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: tloc_f => NULL() !< Face thread ownership
  CLASS(oft_bmesh), POINTER :: bmesh => NULL() !< Boundary mesh
CONTAINS
  !> Needs Docs
  PROCEDURE(mesh_setup), DEFERRED :: setup
  !> Needs Docs
  PROCEDURE(mesh_set_order), DEFERRED :: set_order
  !> Needs Docs
  PROCEDURE(mesh_invert_cell), DEFERRED :: invert_cell
  !> Needs Docs
  PROCEDURE(mesh_log2phys), DEFERRED :: log2phys
  !> Needs Docs
  PROCEDURE(mesh_phys2log), DEFERRED :: phys2log
  !> Needs Docs
  PROCEDURE(mesh_jacobian), DEFERRED :: jacobian
  !> Needs Docs
  PROCEDURE(mesh_hessian), DEFERRED :: hessian
  !> Needs Docs
  PROCEDURE(mesh_snormal), DEFERRED :: snormal
  !> Needs Docs
  PROCEDURE(mesh_ctang), DEFERRED :: ctang
  !> Needs Docs
  PROCEDURE(mesh_get_surf_map), DEFERRED :: get_surf_map
  !> Needs Docs
  PROCEDURE(mesh_surf_to_vol), DEFERRED :: surf_to_vol
  !> Needs Docs
  PROCEDURE(mesh_vlog), DEFERRED :: vlog
  !> Needs Docs
  PROCEDURE(mesh_in_cell), DEFERRED :: in_cell
  !> Needs Docs
  PROCEDURE(mesh_quad_rule), DEFERRED :: quad_rule
  !> Needs Docs
  PROCEDURE(mesh_tessellate), DEFERRED :: tessellate
  !> Needs Docs
  PROCEDURE(mesh_tessellated_sizes), DEFERRED :: tessellated_sizes
  !> Needs Docs
  PROCEDURE :: setup_io => mesh_setup_io
  !> Needs Docs
  PROCEDURE :: save_cell_scalar => mesh_save_cell_scalar
  !> Needs Docs
  PROCEDURE :: save_cell_vector => mesh_save_cell_vector
  !> Needs Docs
  PROCEDURE :: save_vertex_scalar => mesh_save_vertex_scalar
  !> Needs Docs
  PROCEDURE :: save_vertex_vector => mesh_save_vertex_vector
  !> Needs Docs
  PROCEDURE :: volume => mesh_volume
  !> Delete mesh object
  PROCEDURE :: delete => mesh_destroy
END TYPE oft_mesh
!
ABSTRACT INTERFACE
  !------------------------------------------------------------------------------
  !> Setup mesh with implementation specifics (`cell_np`, `cell_ne`, etc.)
  !------------------------------------------------------------------------------
  SUBROUTINE mesh_setup(self,cad_type)
  IMPORT oft_mesh, i4
  CLASS(oft_mesh), INTENT(inout) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cad_type !< CAD/mesh interface ID number
  END SUBROUTINE mesh_setup
  !------------------------------------------------------------------------------
  !> Set maximum order of spatial mapping
  !------------------------------------------------------------------------------
  SUBROUTINE mesh_set_order(self,order)
  IMPORT oft_mesh, i4
  CLASS(oft_mesh), INTENT(inout) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: order !< Maximum order of spatial mapping
  END SUBROUTINE mesh_set_order
  !------------------------------------------------------------------------------
  !> Turn cell "inside out", used to ensure consistent orientations
  !------------------------------------------------------------------------------
  SUBROUTINE mesh_invert_cell(self,cell)
  IMPORT oft_mesh, i4
  CLASS(oft_mesh), INTENT(inout) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cell !< Index of cell to invert
  END SUBROUTINE mesh_invert_cell
  !------------------------------------------------------------------------------
  !> Map from logical to physical coordinates in a given cell
  !------------------------------------------------------------------------------
  FUNCTION mesh_log2phys(self,cell,f) RESULT(pt)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cell !< Index of cell for evaulation
  REAL(r8), INTENT(in) :: f(:) !< Logical coordinate in cell [4]
  REAL(r8) :: pt(3) !< Physical position [3]
  END FUNCTION mesh_log2phys
  !------------------------------------------------------------------------------
  !> Map from physical to logical coordinates in a given cell
  !------------------------------------------------------------------------------
  SUBROUTINE mesh_phys2log(self,cell,pt,f)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cell !< Index of cell for evaulation
  REAL(r8), INTENT(in) :: pt(3) !< Physical position [3]
  REAL(r8), INTENT(out) :: f(:) !< Logical coordinates within the cell [4]
  END SUBROUTINE mesh_phys2log
  !------------------------------------------------------------------------------
  !> Compute the spatial jacobian matrix and its determinant for a given cell at a given logical position
  !------------------------------------------------------------------------------
  SUBROUTINE mesh_jacobian(self,cell,f,gop,j)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cell !< Index of cell for evaulation
  REAL(r8), INTENT(in) :: f(:) !< Logical coordinate in cell [4]
  REAL(r8), INTENT(out) :: gop(:,:) !< Jacobian matrix \f$ (\frac{\partial x_i}{\partial \lambda_j})^{-1} \f$ [3,4]
  REAL(r8), INTENT(out) :: j !< Jacobian of transformation from logical to physical coordinates
  END SUBROUTINE mesh_jacobian
  !------------------------------------------------------------------------------
  !> Compute the spatial hessian matrices for a given cell at a given logical position
  !------------------------------------------------------------------------------
  SUBROUTINE mesh_hessian(self,cell,f,g2op,K)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cell !< Index of cell for evaulation
  REAL(r8), INTENT(in) :: f(:) !< Logical coordinate in cell [4]
  REAL(r8), INTENT(out) :: g2op(:,:) !< Second order Jacobian matrix
  !! \f$ (\frac{\partial x_i}{\partial \lambda_l} \frac{\partial x_j}{\partial \lambda_k})^{-1} \f$
  REAL(r8), INTENT(out) :: K(:,:) !< Gradient correction matrix
  !! \f$ \frac{\partial^2 x_i}{\partial \lambda_k \partial \lambda_l}\f$ [10,3]
  END SUBROUTINE mesh_hessian
  !------------------------------------------------------------------------------
  !> Compute the surface normal vector for a given face on a cell
  !!
  !! If face is not a global boundary face the function returns with `norm = 0`
  !!
  !! @note The logical position in the cell must be on the chosen face for this
  !! subroutine, else an error will be thrown
  !------------------------------------------------------------------------------
  SUBROUTINE mesh_snormal(self,cell,ind,f,norm)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cell !< Index of cell
  INTEGER(i4), INTENT(in) :: ind !< Index of edge within cell
  REAL(r8), INTENT(in) :: f(:) !< Logical coordinate in cell [4]
  REAL(r8), INTENT(out) :: norm(3) !< Unit vector normal to the face [3]
  END SUBROUTINE mesh_snormal
  !------------------------------------------------------------------------------
  !> Compute the curve tangent vector for a given edge on a cell
  !!
  !! If edge is not a global boundary edge the function returns with `tang = 0`
  !!
  !! @note The logical position in the cell must be on the chosen edge for this
  !! subroutine to return a meaningful result
  !------------------------------------------------------------------------------
  SUBROUTINE mesh_ctang(self,cell,ind,f,tang)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: cell !< Index of cell
  INTEGER(i4), INTENT(in) :: ind !< Index of edge within cell
  REAL(r8), INTENT(in) :: f(:) !< Logical coordinate in cell [4]
  REAL(r8), INTENT(out) :: tang(3) !< Unit vector tangent to the edge [3]
  END SUBROUTINE mesh_ctang
  !---------------------------------------------------------------------------------
  !> Get mapping between boundary and volume logical coordinates
  !---------------------------------------------------------------------------------
  SUBROUTINE mesh_get_surf_map(self,face,cell,lmap)
  IMPORT oft_mesh, i4
  CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: face !< Index of face on boundary mesh
  INTEGER(i4), INTENT(out) :: cell !< Cell containing face
  INTEGER(i4), INTENT(out) :: lmap(3) !< Coordinate mapping
  END SUBROUTINE mesh_get_surf_map
  !---------------------------------------------------------------------------------
  !> Map between surface and volume logical coordinates
  !---------------------------------------------------------------------------------
  SUBROUTINE mesh_surf_to_vol(self,fsurf,lmap,fvol)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
  REAL(r8), INTENT(in) :: fsurf(:) !< Surface coordinates [3]
  INTEGER(i4), INTENT(in) :: lmap(3) !< Coordinate mapping
  REAL(r8), INTENT(out) :: fvol(:) !< Volume coordinates [4]
  END SUBROUTINE mesh_surf_to_vol
  !------------------------------------------------------------------------------
  !> Get position in logical space of vertex `i`
  !------------------------------------------------------------------------------
  SUBROUTINE mesh_vlog(self,i,f)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: i !< Vertex to locate
  REAL(r8), INTENT(out) :: f(:) !< Logical coordinates of vertex `i`
  END SUBROUTINE mesh_vlog
  !------------------------------------------------------------------------------
  !> Test if logical position lies within the base cell
  !!
  !! @returns Position `f` is inside the base cell?
  !------------------------------------------------------------------------------
  FUNCTION mesh_in_cell(self,f,tol) RESULT(eface)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
  REAL(r8), INTENT(in) :: f(:) !< Logical coordinate to evaluate
  REAL(r8), INTENT(in) :: tol !< Tolerance for test
  INTEGER(i4) :: eface
  END FUNCTION mesh_in_cell
  !------------------------------------------------------------------------------
  !> Retrieve suitable quadrature rule for mesh with given order
  !------------------------------------------------------------------------------
  SUBROUTINE mesh_quad_rule(self,order,quad_rule)
  IMPORT oft_mesh, oft_quad_type, i4
  CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
  INTEGER(i4), INTENT(in) :: order !< Desired order of quadrature rule
  TYPE(oft_quad_type), INTENT(out) :: quad_rule !< Resulting quadrature rule
  END SUBROUTINE mesh_quad_rule
  !------------------------------------------------------------------------------
  !> Tessellate mesh onto lagrange FE nodes of specified order (usually for plotting)
  !!
  !! @note The maximum tessellation order currently supported is 4
  !! (may be lower for certain mesh types).
  !!
  !! @warning Cell lists are returned with zero based indexing
  !------------------------------------------------------------------------------
  SUBROUTINE mesh_tessellate(self,rtmp,lctmp,order)
  IMPORT oft_mesh, i4, r8
  CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
  REAL(r8), POINTER, DIMENSION(:,:), INTENT(out) :: rtmp !< Tessellated point list [3,:]
  INTEGER(i4), POINTER, DIMENSION(:,:), INTENT(out) :: lctmp !< Tessellated cell list [self%ncp,:]
  INTEGER(i4), INTENT(in) :: order !< Tessellation order
  END SUBROUTINE mesh_tessellate
  !------------------------------------------------------------------------------
  !> Get sizes of arrays returned by @ref mesh_tessellate
  !------------------------------------------------------------------------------
  function mesh_tessellated_sizes(self) result(sizes)
  IMPORT oft_mesh, i4
  CLASS(oft_mesh), INTENT(in) :: self !< Mesh object
  integer(i4) :: sizes(2) !< Array sizes following tessellation [np_tess,nc_tess]
  end function mesh_tessellated_sizes
END INTERFACE
!---
INTEGER(i4), PRIVATE, PARAMETER :: ho_find_retry=20 !< Number of retry attempts during high order find_cell
#ifdef OFT_PLOT_DOUBLE
LOGICAL, PARAMETER :: PLOT_R4_FLAG=.FALSE.
#else
LOGICAL, PARAMETER :: PLOT_R4_FLAG=.TRUE.
#endif
CONTAINS
!---------------------------------------------------------------------------------
!> Checks if a global mesh cell is curved or not
!!
!! @result (T/F) cell is curved?
!---------------------------------------------------------------------------------
function cell_is_curved(self,cell) result(curved)
class(oft_amesh), intent(in) :: self !< Mesh containing cell
integer(i4), intent(in) :: cell !< Index of cell to check
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
!---------------------------------------------------------------------------------
!> Create @ref oft_seam object from mesh connectivity
!---------------------------------------------------------------------------------
subroutine oft_init_seam(self,seam_obj)
class(oft_amesh), intent(in) :: self !< Mesh containing cell
type(oft_seam), intent(out) :: seam_obj !< Resulting seam object
DEBUG_STACK_PUSH
IF(ASSOCIATED(self%global%seam))THEN
  seam_obj%nproc_con=self%global%seam%nproc_con
  seam_obj%proc_split=self%global%seam%proc_split
  seam_obj%proc_con=>self%global%seam%proc_con
  seam_obj%send_reqs=>self%global%seam%send_reqs
  seam_obj%recv_reqs=>self%global%seam%recv_reqs
ELSE
  seam_obj%nproc_con=0
END IF
DEBUG_STACK_POP
end subroutine oft_init_seam
!------------------------------------------------------------------------------
!> Estimate mesh volume
!------------------------------------------------------------------------------
function mesh_volume(self) result(vol)
class(oft_mesh), INTENT(IN) :: self !< Needs docs
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
!------------------------------------------------------------------------------
!> Estimate mesh volume
!------------------------------------------------------------------------------
SUBROUTINE mesh_setup_io(self,xdmf_obj,tess_order)
class(oft_mesh), INTENT(inout) :: self !< Needs docs
class(xdmf_plot_file), intent(inout) :: xdmf_obj !< Needs docs
integer(i4), intent(in) :: tess_order !< Needs docs
integer(i4) :: i,j,k,id,error,two=2,io_sizes(2),dims(2)
integer(i4), POINTER :: lctmp(:,:),fmap(:)
real(r8), POINTER :: rtmp(:,:),reg_tmp(:)
DEBUG_STACK_PUSH
!---Initialize HDF5 and open mesh file
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Writing mesh to plot files'
CALL oft_increase_indent
self%tess_order=tess_order
self%bmesh%tess_order=tess_order
!---Get grid tessellation
CALL self%tessellate(rtmp, lctmp, self%tess_order)
IF(TRIM(self%meshname)=='none')self%meshname='vmesh'
CALL xdmf_obj%add_mesh(30+self%type,rtmp,lctmp,self%meshname)
deallocate(rtmp)
dims=SHAPE(lctmp)
k=dims(2)/self%nc
deallocate(lctmp)
!---Setup I/O for boundary mesh
CALL self%bmesh%setup_io(xdmf_obj,self%tess_order)
!---Write regions
ALLOCATE(reg_tmp(dims(2)))
!$omp parallel do private(j,id)
DO i=1,self%nc
  id=self%reg(i)
  DO j=1,k
  reg_tmp((i-1)*k+j)=REAL(id,8)
  END DO
END DO
CALL self%save_cell_scalar(reg_tmp,xdmf_obj,'REG_vol')
DEALLOCATE(reg_tmp)
!---Write surface IDs
IF(self%bmesh%nc>0)THEN
  ALLOCATE(fmap(self%nf))
  fmap=0
  !$omp parallel do
  do i=1,self%nbf
    fmap(self%lbf(i))=i
  end do
  io_sizes=self%bmesh%tessellated_sizes()
  k=io_sizes(2)/self%bmesh%nc
  ALLOCATE(reg_tmp(io_sizes(2)))
  !$omp parallel do private(j,id)
  DO i=1,self%bmesh%nc
    id=self%bfs(fmap(ABS(self%bmesh%parent%lf(i))))
    DO j=1,k
      reg_tmp((i-1)*k+j)=REAL(id,8)
    END DO
  END DO
  CALL self%bmesh%save_cell_scalar(reg_tmp,xdmf_obj,'SID')
  DEALLOCATE(reg_tmp,fmap)
END IF
CALL oft_decrease_indent
DEBUG_STACK_POP
END SUBROUTINE mesh_setup_io
!------------------------------------------------------------------------------
!> Write scalar vertex data out to file
!------------------------------------------------------------------------------
subroutine mesh_save_vertex_scalar(self,p,xdmf_obj,path)
class(oft_mesh), INTENT(IN) :: self
real(r8), intent(in) :: p(:) !< Vertex data [np]
class(xdmf_plot_file), intent(in) :: xdmf_obj !< XDMF save object
character(LEN=*), intent(in) :: path !< Name of the output field
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving scalar plot field: ',TRIM(path)
sizes=self%tessellated_sizes()
IF(SIZE(p,DIM=1)/=sizes(1))CALL oft_abort("Incorrect array size","mesh_save_vertex_scalar",__FILE__)
CALL xdmf_obj%write(p,self%meshname,path,1,PLOT_R4_FLAG)
DEBUG_STACK_POP
end subroutine mesh_save_vertex_scalar
!------------------------------------------------------------------------------
!> Write scalar cell data out to file
!------------------------------------------------------------------------------
subroutine mesh_save_cell_scalar(self,p,xdmf_obj,path)
class(oft_mesh), INTENT(IN) :: self
real(r8), intent(in) :: p(:) !< Cell data [nc]
class(xdmf_plot_file), intent(in) :: xdmf_obj !< XDMF save object
character(LEN=*), intent(in) :: path !< Name of the output field
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving scalar plot field: ',TRIM(path)
sizes=self%tessellated_sizes()
IF(SIZE(p,DIM=1)/=sizes(2))CALL oft_abort("Incorrect array size","mesh_save_cell_scalar",__FILE__)
CALL xdmf_obj%write(p,self%meshname,path,2,PLOT_R4_FLAG)
DEBUG_STACK_POP
end subroutine mesh_save_cell_scalar
!------------------------------------------------------------------------------
!> Write vector vertex data out to file
!------------------------------------------------------------------------------
subroutine mesh_save_vertex_vector(self,bv,xdmf_obj,path)
class(oft_mesh), INTENT(IN) :: self
real(r8), intent(in) :: bv(:,:) !< Vertex data [3,np]
class(xdmf_plot_file), intent(in) :: xdmf_obj !< XDMF save object
character(LEN=*), intent(in) :: path !< Name of the output field
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving vector plot field: ',TRIM(path)
IF(SIZE(bv,DIM=1)/=3)CALL oft_abort("Output array is not 3 vector","mesh_save_vertex_vector",__FILE__)
sizes=self%tessellated_sizes()
IF(SIZE(bv,DIM=2)/=sizes(1))CALL oft_abort("Incorrect array size","mesh_save_vertex_vector",__FILE__)
CALL xdmf_obj%write(bv,self%meshname,path,1,PLOT_R4_FLAG)
DEBUG_STACK_POP
end subroutine mesh_save_vertex_vector
!------------------------------------------------------------------------------
!> Write vector cell data out to file
!------------------------------------------------------------------------------
subroutine mesh_save_cell_vector(self,bcc,xdmf_obj,path)
class(oft_mesh), INTENT(IN) :: self
real(r8), intent(in) :: bcc(:,:) !< Cell data [3,nc]
class(xdmf_plot_file), intent(in) :: xdmf_obj !< XDMF save object
character(LEN=*), intent(in) :: path !< Name of the output field
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving vector plot field: ',TRIM(path)
IF(SIZE(bcc,DIM=1)/=3)CALL oft_abort("Output array is not 3 vector","mesh_save_cell_vector",__FILE__)
sizes=self%tessellated_sizes()
IF(SIZE(bcc,DIM=2)/=sizes(2))CALL oft_abort("Incorrect array size","mesh_save_cell_vector",__FILE__)
CALL xdmf_obj%write(bcc,self%meshname,path,2,PLOT_R4_FLAG)
DEBUG_STACK_POP
end subroutine mesh_save_cell_vector
!------------------------------------------------------------------------------
!> Destroy mesh object
!------------------------------------------------------------------------------
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
NULLIFY(self%fstitch%be,self%fstitch%lbe)
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
!---------------------------------------------------------------------------------
!> Find physical point in mesh.
!!
!! For high order grids an approximate guess is first computed using only the linear
!! portion of the mesh representation. This guess is then refined using the full mesh
!! representation and the @ref tetmesh_mapping::tetmesh_phys2logho "tetmesh_phys2logho"
!! subroutine. The maximum number of cell searches during this non-linear phase is controlled
!! by the module variable @ref tetmesh_mapping::ho_find_retry "ho_find_retry". For
!! more information see the documentation for @ref tetmesh_mapping::tetmesh_phys2logho
!! "tetmesh_phys2logho"
!---------------------------------------------------------------------------------
subroutine mesh_findcell(self,cell,pt,fout)
class(oft_mesh), target, intent(inout) :: self !< Mesh to search
integer(i4), intent(inout) :: cell !< Cell containing point on output, guess on input
real(r8), intent(in) :: pt(3) !< Coordinates to locate [3]
real(r8), optional, intent(out) :: fout(4) !< Logical coordinates of point in cell (optional)
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
!---------------------------------------------------------------------------------
!> Find physical point in mesh using robust method
!!
!! First, a collection of the `nclosest` mesh vertices is found. Then the surrounding
!! cells are searched for the specified point. If the point is not found the
!! standard @ref tetmesh_mapping::tetmesh_findcell "tetmesh_findcell" subroutine is
!! used as a fallback.
!---------------------------------------------------------------------------------
subroutine mesh_findcell2(self,cell,pt,nclosest,fout)
class(oft_mesh), target, intent(inout) :: self !< Mesh to search
integer(i4), intent(inout) :: cell !< Cell containing point on output, guess on input
real(r8), intent(in) :: pt(3) !< Coordinates to locate [3]
integer(i4), intent(in) :: nclosest !< Number of candidate vertices to use for search
real(r8), optional, intent(out) :: fout(4) !< Logical coordinates of point in cell (optional)
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
!------------------------------------------------------------------------------
!> Estimate mesh area
!------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------
!> Find physical point in mesh.
!!
!! @warning Only works for 2D meshes and uses linear logical mapping
!---------------------------------------------------------------------------------
subroutine bmesh_findcell(self,cell,pt,fout)
class(oft_bmesh), intent(in) :: self !< Mesh to search
integer(i4), intent(inout) :: cell !< Cell containing point on output, guess on input
real(r8), intent(in) :: pt(2) !< Coordinates to locate [2]
real(r8), optional, intent(out) :: fout(3) !< Logical coordinates of point in cell (optional)
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
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE bmesh_setup_io(self,xdmf_obj,tess_order)
CLASS(oft_bmesh), INTENT(inout) :: self !< Needs docs
class(xdmf_plot_file), intent(inout) :: xdmf_obj !< Needs docs
integer(i4), intent(in) :: tess_order !< Needs docs
logical :: create_files
integer(i4) :: i,k,j,id,error,dims(2)
integer(i4), POINTER :: lftmp(:,:),fmap(:)
real(r8), POINTER :: rtmp(:,:),reg_tmp(:)
DEBUG_STACK_PUSH
!---Setup I/O
if(oft_debug_print(1))write(*,'(2A)')oft_indent,'Writing boundary mesh to plot files'
CALL oft_increase_indent
self%tess_order=tess_order
IF(self%nc==0)THEN
  CALL oft_decrease_indent
  DEBUG_STACK_POP
  RETURN
END IF
!---Get grid tessellation
CALL self%tessellate(rtmp, lftmp, self%tess_order)
IF(TRIM(self%meshname)=='none')self%meshname='smesh'
CALL xdmf_obj%add_mesh(20+self%type,rtmp,lftmp,self%meshname)
deallocate(rtmp)
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
CALL self%save_cell_scalar(reg_tmp,xdmf_obj,'REG_surf')
DEALLOCATE(reg_tmp)
CALL oft_decrease_indent
DEBUG_STACK_POP
END SUBROUTINE bmesh_setup_io
!------------------------------------------------------------------------------
!> Write scalar vertex data out to file
!------------------------------------------------------------------------------
subroutine bmesh_save_vertex_scalar(self,p,xdmf_obj,path)
class(oft_bmesh), INTENT(IN) :: self
real(r8), intent(in) :: p(:) !< Vertex data [np]
class(xdmf_plot_file), intent(in) :: xdmf_obj !< XDMF save object
character(LEN=*), intent(in) :: path !< Name of the output field
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving scalar plot field: ',TRIM(path)
sizes=self%tessellated_sizes()
IF(SIZE(p,DIM=1)/=sizes(1))CALL oft_abort("Incorrect array size","bmesh_save_vertex_scalar",__FILE__)
CALL xdmf_obj%write(p,self%meshname,path,1,PLOT_R4_FLAG)
DEBUG_STACK_POP
end subroutine bmesh_save_vertex_scalar
!------------------------------------------------------------------------------
!> Write scalar cell data out to file
!------------------------------------------------------------------------------
subroutine bmesh_save_cell_scalar(self,p,xdmf_obj,path)
class(oft_bmesh), INTENT(IN) :: self
real(r8), intent(in) :: p(:) !< Cell data [nc]
class(xdmf_plot_file), intent(in) :: xdmf_obj !< XDMF save object
character(LEN=*), intent(in) :: path !< Name of the output field
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving scalar plot field: ',TRIM(path)
sizes=self%tessellated_sizes()
IF(SIZE(p,DIM=1)/=sizes(2))CALL oft_abort("Incorrect array size","bmesh_save_cell_scalar",__FILE__)
CALL xdmf_obj%write(p,self%meshname,path,2,PLOT_R4_FLAG)
DEBUG_STACK_POP
end subroutine bmesh_save_cell_scalar
!------------------------------------------------------------------------------
!> Write vector vertex data out to file
!------------------------------------------------------------------------------
subroutine bmesh_save_vertex_vector(self,bv,xdmf_obj,path)
class(oft_bmesh), INTENT(IN) :: self
real(r8), intent(in) :: bv(:,:) !< Vertex data [3,np]
class(xdmf_plot_file), intent(in) :: xdmf_obj !< XDMF save object
character(LEN=*), intent(in) :: path !< Name of the output field
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving vector plot field: ',TRIM(path)
IF(SIZE(bv,DIM=1)/=3)CALL oft_abort("Output array is not 3 vector","bmesh_save_vertex_vector",__FILE__)
sizes=self%tessellated_sizes()
IF(SIZE(bv,DIM=2)/=sizes(1))CALL oft_abort("Incorrect array size","bmesh_save_vertex_vector",__FILE__)
CALL xdmf_obj%write(bv,self%meshname,path,1,PLOT_R4_FLAG)
DEBUG_STACK_POP
end subroutine bmesh_save_vertex_vector
!------------------------------------------------------------------------------
!> Write vector cell data out to file
!------------------------------------------------------------------------------
subroutine bmesh_save_cell_vector(self,bcc,xdmf_obj,path)
class(oft_bmesh), INTENT(IN) :: self
real(r8), intent(in) :: bcc(:,:) !< Cell data [3,nc]
class(xdmf_plot_file), intent(in) :: xdmf_obj !< XDMF save object
character(LEN=*), intent(in) :: path !< Name of the output field
integer(i4) :: sizes(2)
DEBUG_STACK_PUSH
IF(oft_debug_print(1))WRITE(*,'(3A)')oft_indent,'Saving vector plot field: ',TRIM(path)
IF(SIZE(bcc,DIM=1)/=3)CALL oft_abort("Output array is not 3 vector","bmesh_save_cell_vector",__FILE__)
sizes=self%tessellated_sizes()
IF(SIZE(bcc,DIM=2)/=sizes(2))CALL oft_abort("Incorrect array size","bmesh_save_cell_vector",__FILE__)
CALL xdmf_obj%write(bcc,self%meshname,path,2,PLOT_R4_FLAG)
DEBUG_STACK_POP
end subroutine bmesh_save_cell_vector
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
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
