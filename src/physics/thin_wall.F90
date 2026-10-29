!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file thin_wall.F90
!
!> Module for thin-wall modeling on 3D triangular meshes
!!
!!
!! @authors Chris Hansen
!! @date May 2017
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
MODULE thin_wall
USE, INTRINSIC :: iso_c_binding, only: c_loc
USE oft_base
USE oft_sort, ONLY: sort_array, search_array
USE oft_io, ONLY: xdmf_plot_file, hdf5_field_exist, hdf5_rst, hdf5_write, hdf5_read, &
  hdf5_rst_destroy, hdf5_create_file, hdf5_field_get_sizes, hdf5_create_group
USE oft_quadrature
USE oft_mesh_type, ONLY: oft_bmesh
USE oft_mesh_local, ONLY: bmesh_local_init
USE oft_mesh_local_util, ONLY: mesh_local_findedge
USE oft_mesh_global_util, ONLY: mesh_global_resolution
USE oft_tetmesh_type, ONLY: oft_tetmesh
USE oft_trimesh_type, ONLY: oft_trimesh
USE oft_tet_quadrature, ONLY: set_quad_2d
USE oft_la_base, ONLY: oft_vector, oft_cvector
USE oft_la_utils, ONLY: csr_remove_redundant
USE oft_deriv_matrices, ONLY: oft_noop_matrix
USE oft_native_la, ONLY: oft_native_vector, oft_native_matrix, native_vector_cast, &
  native_vector_slice_push, native_vector_slice_pop, oft_native_dense_matrix
USE oft_solver_base, ONLY: oft_solver, oft_csolver
USE oft_lu, ONLY: lapack_matinv, lapack_cholesky
USE axi_green, ONLY: green, grad_green, axi_coil_set
USE mhd_utils, ONLY: mu0
#define MAX_EDGE_CONN 6
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------------
!> Structure containing definition of "hole" elements for multiply connected
!! surfaces
!---------------------------------------------------------------------------------
TYPE :: hole_mesh
  INTEGER(i4) :: i0 = 0 !< Starting point for point chain
  INTEGER(i4) :: n = 0 !< Number of points in chain
  INTEGER(i4), POINTER, DIMENSION(:) :: lp => NULL() !< List of points in chain
  INTEGER(i4), POINTER, DIMENSION(:) :: lp_sorted => NULL() !< Sorted list of points in chain
  INTEGER(i4), POINTER, DIMENSION(:) :: sort_ind => NULL() !< Index from sorted list
  INTEGER(i4), POINTER, DIMENSION(:) :: kpc => NULL() !< Pointer to point-cell linkage
  INTEGER(i4), POINTER, DIMENSION(:) :: lpc => NULL() !< List of cells tied to each point
  ! REAL(r8), POINTER, DIMENSION(:) :: fsign => NULL() !< Sign of
  REAL(r8) :: ptcc(3) = 0.d0
END TYPE hole_mesh
!---------------------------------------------------------------------------------
!> Structure containing definition of a flux loop sensor
!---------------------------------------------------------------------------------
TYPE :: floop_sensor
  INTEGER(i4) :: np = 0 !< Number of points in loop
  REAL(r8) :: scale_fac = 1.d0 !< Scale factor to apply to signal
  REAL(r8), POINTER, DIMENSION(:,:) :: r => NULL() !< List of points [3,np]
  CHARACTER(LEN=40) :: name = '' !< Name of sensor
END TYPE floop_sensor
!---------------------------------------------------------------------------------
!> Structure containing definition of a current jumper sensor
!---------------------------------------------------------------------------------
TYPE :: jumper_sensor
  INTEGER(i4) :: np = 0 !< Number of points on jumper
  CHARACTER(LEN=40) :: name = '' !< Name of sensor
  INTEGER(i4), POINTER, DIMENSION(:) :: points => NULL() !< List of points on jumper
  REAL(r8), POINTER, DIMENSION(:) :: hole_facs => NULL() !< Coupling weight to "holes"
CONTAINS
  !> Compute current [A] through jumper from thin-wall solution vector
  PROCEDURE :: compute_current => jumper_compute_current
END TYPE jumper_sensor
!---------------------------------------------------------------------------------
!> Structure containing sensor sets
!---------------------------------------------------------------------------------
TYPE :: tw_sensors
  INTEGER(i4) :: nfloops = 0 !< Number of flux loops
  INTEGER(i4) :: njumpers = 0 !< Number of current jumpers
  TYPE(floop_sensor), POINTER, DIMENSION(:) :: floops => NULL() !< List of flux loops
  TYPE(jumper_sensor), POINTER, DIMENSION(:) :: jumpers => NULL() !< List of current jumpers
END TYPE tw_sensors
!---------------------------------------------------------------------------------
!> Junction data of a model whose mesh is split at junctions (edges shared by more than two sheets)
!!
!! The model mesh is the split, manifold mesh: each junction vertex has one copy per fan of cells, and the copies
!! are joined by the junction map P (current conservation across junction edges), \f$ x_{exp} = P x \f$. Quantities
!! with the suffix `_exp` ("expanded") are on the split mesh, before this reduction, and only exist for models with
!! junctions.
!---------------------------------------------------------------------------------
TYPE :: tw_junctions
  INTEGER(i4) :: nelems_exp = 0 !< Number of elements on the split mesh (before the junction reduction)
  INTEGER(i4) :: nred = 0 !< Number of reduced surface unknowns
  INTEGER(i4) :: nholes = 0 !< Number of explicit junction holes, after the mesh file holes
  INTEGER(i4), POINTER, DIMENSION(:) :: p_kr => NULL() !< Junction map P, row pointer [nelems_exp+1] (V-coil rows appended at setup)
  INTEGER(i4), POINTER, DIMENSION(:) :: p_lc => NULL() !< Junction map P, column indices
  REAL(r8), POINTER, DIMENSION(:) :: p_val => NULL() !< Junction map P, values
  INTEGER(i4), POINTER, DIMENSION(:,:) :: table => NULL() !< Junction table (vertex, vertex, cell, sign) [4,n]
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: holes => NULL() !< Explicit junction holes, (cell, corner) pairs
  TYPE(oft_native_matrix), POINTER :: Rmat_exp => NULL() !< Resistance matrix before the junction reduction
END TYPE tw_junctions
!---------------------------------------------------------------------------------
!> Structure containing filament coil definition
!---------------------------------------------------------------------------------
TYPE :: tw_gen_coil
  LOGICAL :: is_ring = .FALSE. !< Circular coil about the Z axis, defined by its radius and height?
  INTEGER(i4) :: npts = 0 !< Number of points in coil
  REAL(r8) :: ring_rz(2) = 0.d0 !< Radius and height of a circular coil
  REAL(r8), POINTER, DIMENSION(:,:) :: pts => NULL() !< Points [3,npts]
END TYPE tw_gen_coil
!---------------------------------------------------------------------------------
!> Structure containing a coil sets composed of one or more individual filament coils
!---------------------------------------------------------------------------------
TYPE :: tw_coil_set
  LOGICAL :: sens_mask = .FALSE. !< Mask from sensor output?
  LOGICAL :: is_vcoil = .FALSE. !< Is this a voltage-specified coil set?
  INTEGER(i4) :: ncoils = 0 !< Number of coils in set
  REAL(r8) :: Lself = -1.d0 !< Self inductance of coil set (if required)
  REAL(r8) :: Rself = -1.d0 !< Resistance of coil set (if required)
  REAL(r8), POINTER, DIMENSION(:) :: scales => NULL() !< Scale factor for each coil
  REAL(r8), POINTER, DIMENSION(:) :: res_per_len => NULL() !< Resistance/length of each coil (if required)
  REAL(r8), POINTER, DIMENSION(:) :: radius => NULL() !< Effective radius of each coil (for calculation of Lself)
  ! REAL(r8), POINTER, DIMENSION(:,:) :: axi_pt => NULL() !< Coil definitions for axisymmetric coils
  CHARACTER(LEN=40) :: name = '' !< Name of coil set
  TYPE(tw_gen_coil), POINTER, DIMENSION(:) :: coils => NULL() !< List of coils
END TYPE tw_coil_set
!---------------------------------------------------------------------------------
!> Class for thin-wall simulation
!---------------------------------------------------------------------------------
TYPE :: tw_type
  LOGICAL :: model_frozen = .FALSE. !< Flag indicating if model is frozen (no more shape changes allowed)
  INTEGER(i4) :: nelems = 0 !< Number of elements in model (np_active+nholes+n_vcoils)
  INTEGER(i4) :: np_active = 0 !< Number of active vertices in model
  INTEGER(i4) :: n_coils = 0 !< Number of general coils in model
  INTEGER(i4) :: n_vcoils = 0 !< Number of voltage-specified coils in model
  INTEGER(i4) :: n_icoils = 0 !< Number of current-specified coils in model
  INTEGER(i4) :: nholes = 0 !< Number of "holes" in model
  INTEGER(i4) :: holes_version = 0 !< Version of the hole construction (0: unknown, >0: from "OFT_ThinCurr_holes")
  TYPE(tw_junctions), POINTER :: junc => NULL() !< Junction data (only for meshes split at junctions)
  INTEGER(i4) :: nclosures = 0 !< Number of "closures" in model
  INTEGER(i4) :: nfh = 0 !< Number of face-hole interactions
  REAL(r8) :: B_dx = 1.d-3 !< Spatial step size for finite difference evaluation of B-field
  LOGICAL, POINTER, DIMENSION(:) :: sens_mask => NULL() !< Mask array for sensors [nreg]
  INTEGER(i4), POINTER, DIMENSION(:) :: pmap => NULL() !< Map from mesh vertices to active vertices [mesh%np]
  INTEGER(i4), POINTER, DIMENSION(:) :: kpmap_inv => NULL() !< Map from active vertices to mesh vertices [np_active]
  INTEGER(i4), POINTER, DIMENSION(:) :: lpmap_inv => NULL() !< Map from active vertices to mesh vertices [mesh%np]
  INTEGER(i4), POINTER, DIMENSION(:) :: closures => NULL() !< List of closure vertices [nclosures]
  INTEGER(i4), POINTER, DIMENSION(:) :: kfh => NULL() !< Pointer to face-hole interaction list [mesh%nc+1]
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lfh => NULL() !< List of face-hole interactions [nfh]
  REAL(r8), POINTER, DIMENSION(:) :: Eta_vol => NULL() !< Resistivity values for each region, divided by mu0 [nreg].
  REAL(r8), POINTER, DIMENSION(:) :: Thickness => NULL() !< Thickness values for each region [nreg]. Units: meters
  REAL(r8), POINTER, DIMENSION(:) :: Eta_surf => NULL() !< Surface resistivity (eta/thickness) values for each region, divided by mu0 [nreg].
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: Ael2dr => NULL() !< Element to driver (icoils) coupling matrix
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: Ael2coil => NULL() !< Element to coil (vcoils+icoils) coupling matrix
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: Ael2sen => NULL() !< Element to sensor coupling matrix
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: Acoil2coil => NULL() !< Coil to coil coupling matrix
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: Lmat => NULL() !< Full inductance matrix
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: Adr2sen => NULL() !< Driver (icoils) to sensor coupling matrix
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: Bel => NULL() !< Element to vertex B-field reconstruction operator [np,nelems,3]
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: Bdr => NULL() !< Driver (icoils) to vertex B-field reconstruction operator [np,n_icoils,3]
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: qbasis => NULL() !< Basis function pre-evaluated at cell centers
  CHARACTER(LEN=OFT_PATH_SLEN) :: rst_prefix = '' !< Path prefix for .rst and .hist output files
  TYPE(xdmf_plot_file) :: xdmf
  CLASS(oft_vector), POINTER :: Uloc => NULL() !< FE vector for thin-wall model
  CLASS(oft_vector), POINTER :: Uloc_pts => NULL() !< Needs docs
  TYPE(oft_native_matrix), POINTER :: Rmat => NULL() !< Resistivity matrix for thin-wall model
  CLASS(oft_bmesh), POINTER :: mesh => NULL() !< Underlying surface mesh
  TYPE(hole_mesh), POINTER, DIMENSION(:) :: hmesh => NULL() !< Hole definitions
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: jumper_nsets => NULL() !< Jumper definitions
  TYPE(tw_coil_set), POINTER, DIMENSION(:) :: coils => NULL() !< List of general coils
  TYPE(tw_coil_set), POINTER, DIMENSION(:) :: vcoils => NULL() !< List of voltage-specified coils
  TYPE(tw_coil_set), POINTER, DIMENSION(:) :: icoils => NULL() !< List of current-specified coils
  TYPE(xml_node) :: xml
CONTAINS
  !> Setup thin-wall model
  PROCEDURE :: setup => tw_setup
  !> Save debug information for model
  PROCEDURE :: save_debug => tw_save_debug
END TYPE tw_type
!---
! REAL(r8) :: mag_dx = 1.d-5
REAL(r8) :: quad_tols(3) = [0.75d0, 0.95d0, 0.995d0] !< Distance tolerances for quadrature order selection
INTEGER(i4) :: quad_orders(3) = [18, 10, 6] !< Quadrature order for each tolerance
REAL(r8), PARAMETER :: target_err = 1.d-8
!> Quadrature target of the inductance cell-pair kernel (@ref tw_cell_pair_kernel). Looser than `target_err`: the
!! order formula assumes the error falls as (2h/d)^order, which is very pessimistic for 1/r between separated
!! triangles. With 1E-6 the change to L is ~1E-9 relative to the matrix (eigenvalues ~1E-9), well below the HODLR/ACA
!! compression tolerances, while assembly is ~1.75x faster than with 1E-8.
REAL(r8), PARAMETER :: tw_L_quad_err = 1.d-6
INTEGER(i4), PARAMETER :: tw_quad_far_max = 10 !< Highest quadrature order of far-field cell-pair integrals (closer pairs integrate one cell analytically)
INTEGER(i4), PARAMETER :: tw_quad_far_np = 25 !< Number of points of the order `tw_quad_far_max` triangle rule (sizes cell-pair work arrays)
REAL(r8), PARAMETER :: coil_min_rad = 1.d-6
TYPE(oft_quad_type), TARGET, SAVE :: tw_quads(18) !< Triangle quadrature rules of orders 1-18 (shared)
LOGICAL, SAVE :: tw_quads_ready = .FALSE. !< Have tw_quads been set up?
integer(i4), public, parameter :: tw_idx_ver=1 !< File version for array indexing
integer(i4), public, parameter :: tw_cache_ver=2 !< Version of cached matrix files (increment when a change invalidates existing caches)
character(LEN=16), public, parameter :: tw_idx_path="ThinCurr_Version" !< HDF5 field name
CONTAINS
!---------------------------------------------------------------------------------
!> Set up the shared triangle quadrature rules (done once, by the first @ref tw_setup)
!---------------------------------------------------------------------------------
SUBROUTINE tw_quads_setup()
INTEGER(i4) :: i
IF(tw_quads_ready)RETURN
DO i=1,18
  CALL set_quad_2d(tw_quads(i),i)
END DO
IF(tw_quads(tw_quad_far_max)%np>tw_quad_far_np)CALL oft_abort('Far-field quadrature rule exceeds "tw_quad_far_np"', &
  'tw_quads_setup',__FILE__)
tw_quads_ready=.TRUE.
END SUBROUTINE tw_quads_setup
!---------------------------------------------------------------------------------
!> Hash of the coil definitions of a model, seeded with the cache file version, for keys of cached matrices
!!
!! With `vcoils_only`, only what the self-inductance matrix depends on: the V-coils, without their resistances.
!---------------------------------------------------------------------------------
FUNCTION tw_coil_hash(self,vcoils_only) RESULT(hash)
TYPE(tw_type), INTENT(in) :: self !< ThinCurr object
LOGICAL, OPTIONAL, INTENT(in) :: vcoils_only !< Hash only the V-coils, without resistances (default: all coils)
INTEGER(4) :: hash
INTEGER(4) :: i
LOGICAL :: lonly
lonly=.FALSE.
IF(PRESENT(vcoils_only))lonly=vcoils_only
hash=tw_cache_ver
IF(.NOT.lonly)THEN
  DO i=1,self%n_icoils
    hash=coil_set_hash(self%icoils(i),hash)
  END DO
END IF
DO i=1,self%n_vcoils
  hash=coil_set_hash(self%vcoils(i),hash)
END DO
CONTAINS
!> Combine one coil set's definition (count, scales, radii, resistances unless vcoils_only, each coil's points or R, Z)
FUNCTION coil_set_hash(coil_set,hash_in) RESULT(hash_out)
TYPE(tw_coil_set), INTENT(in) :: coil_set
INTEGER(4), INTENT(in) :: hash_in
INTEGER(4) :: hash_out,j
INTEGER(4), TARGET :: itmp(2)
REAL(8), TARGET :: rtmp(2)
itmp=[coil_set%ncoils,0]
hash_out=oft_simple_hash(C_LOC(itmp),INT(4*2,8),starting_value=hash_in)
IF(ASSOCIATED(coil_set%scales))hash_out=oft_simple_hash(C_LOC(coil_set%scales),INT(8*coil_set%ncoils,8),starting_value=hash_out)
IF(ASSOCIATED(coil_set%radius))hash_out=oft_simple_hash(C_LOC(coil_set%radius),INT(8*coil_set%ncoils,8),starting_value=hash_out)
IF(ASSOCIATED(coil_set%res_per_len).AND.(.NOT.lonly))hash_out=oft_simple_hash(C_LOC(coil_set%res_per_len),INT(8*coil_set%ncoils,8),starting_value=hash_out)
DO j=1,coil_set%ncoils
  itmp=[coil_set%coils(j)%npts,0]
  IF(coil_set%coils(j)%is_ring)itmp(2)=1
  hash_out=oft_simple_hash(C_LOC(itmp),INT(4*2,8),starting_value=hash_out)
  IF(coil_set%coils(j)%is_ring)THEN
    rtmp=coil_set%coils(j)%ring_rz
    hash_out=oft_simple_hash(C_LOC(rtmp),INT(8*2,8),starting_value=hash_out)
  ELSE
    hash_out=oft_simple_hash(C_LOC(coil_set%coils(j)%pts),INT(8*3*coil_set%coils(j)%npts,8),starting_value=hash_out)
  END IF
END DO
END FUNCTION coil_set_hash
END FUNCTION tw_coil_hash
!---------------------------------------------------------------------------------
!> Hash of the mesh and coil definitions of a model, for keys of cached coil and sensor matrices
!---------------------------------------------------------------------------------
FUNCTION tw_model_hash(self) RESULT(hash)
TYPE(tw_type), INTENT(in) :: self !< ThinCurr object
INTEGER(4) :: hash
hash=tw_coil_hash(self)
hash=oft_simple_hash(C_LOC(self%mesh%lc),INT(4*3*self%mesh%nc,8),starting_value=hash)
hash=oft_simple_hash(C_LOC(self%mesh%r),INT(8*3*self%mesh%np,8),starting_value=hash)
END FUNCTION tw_model_hash
!---------------------------------------------------------------------------------
!> Number of elements on the model mesh before the junction reduction
!!
!! @ref tw_junctions::nelems_exp for meshes split at junctions, otherwise the number of elements (`nelems`)
!---------------------------------------------------------------------------------
PURE FUNCTION tw_nelems_exp(self) RESULT(n)
CLASS(tw_type), INTENT(in) :: self !< Thin-wall model object
INTEGER(i4) :: n
IF(ASSOCIATED(self%junc))THEN
  n=self%junc%nelems_exp
ELSE
  n=self%nelems
END IF
END FUNCTION tw_nelems_exp
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE tw_setup(self,hole_ns,error_string,filepath)
CLASS(tw_type), INTENT(INOUT) :: self !< Thin-wall model object
TYPE(oft_1d_int), POINTER, INTENT(inout) :: hole_ns(:) !< Hole nodesets
CHARACTER(LEN=:), ALLOCATABLE, INTENT(out) :: error_string
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: filepath !< Path to mesh file
LOGICAL :: coils_xml
INTEGER(4) :: i,j,k,l,face,ioffset,ed,error_flag,njholes
INTEGER(4), ALLOCATABLE :: kfh_tmp(:),np_inverse(:)
REAL(8) :: f(3),rgop(3,3),area_i,norm_i(3)
TYPE(xml_node) :: coil_element
!---
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Creating thin-wall model'
CALL oft_increase_indent
CALL tw_quads_setup
CALL mesh_global_resolution(self%mesh)
CALL bmesh_local_init(self%mesh,sync_normals=.TRUE.)
!---Load ThinCurr-specific mesh information
coils_xml=.TRUE.
IF(PRESENT(filepath))THEN
  CALL tw_load_hdf5(self,filepath,hole_ns,error_string)
  IF(self%n_coils>0)coils_xml=.FALSE.
END IF
IF(ALLOCATED(error_string))RETURN
!---Load coils from XML for backwards compatibility
IF(.NOT.self%xml%associated().AND.ASSOCIATED(oft_env%xml))CALL xml_get_element(oft_env%xml,"thincurr",self%xml,error_flag)
IF(coils_xml)THEN
  !---Read V(t) driver coils from XML file
  CALL xml_get_element(self%xml,"vcoils",coil_element,error_flag)
  IF(error_flag==0)THEN
    WRITE(*,'(2A)')oft_indent,'Loading V(t) driver coils from XML file'
    CALL tw_load_coils_xml(coil_element,0,self%n_vcoils,self%vcoils)
  ELSE
    WRITE(*,'(2A)')oft_indent,'No V(t) driver coils found'
  END IF
  DO i=1,self%n_vcoils
    IF(ANY(self%vcoils(i)%res_per_len<0.d0))CALL oft_abort("Invalid resistivity for passive coil", &
      "tw_setup", __FILE__)
    IF(ANY(self%vcoils(i)%radius<coil_min_rad))CALL oft_abort("Invalid radius for passive coil", &
      "tw_setup", __FILE__)
  END DO
  !---Read I(t) driver coils from XML file
  CALL xml_get_element(self%xml,"icoils",coil_element,error_flag)
  IF(error_flag==0)THEN
    WRITE(*,'(2A)')oft_indent,'Loading I(t) driver coils from XML file'
    CALL tw_load_coils_xml(coil_element,self%n_vcoils,self%n_icoils,self%icoils)
  ELSE
    WRITE(*,'(2A)')oft_indent,'No I(t) driver coils found'
  END IF
  DO i=1,self%n_icoils
    DO j=1,self%icoils(i)%ncoils
      self%icoils(i)%radius(j)=MAX(coil_min_rad,self%icoils(i)%radius(j)) ! Remove dummy radius on Icoils
    END DO
  END DO
  !---TEMPORARY to match new behavior of storing all coils together
  self%n_coils=self%n_vcoils+self%n_icoils
  IF(self%n_coils>0)ALLOCATE(self%coils(self%n_coils))
  DO i=1,self%n_vcoils
    CALL tw_copy_coil(self%vcoils(i),self%coils(i))
    self%coils(i)%is_vcoil=.TRUE.
  END DO
  DO i=1,self%n_icoils
    CALL tw_copy_coil(self%icoils(i),self%coils(i+self%n_vcoils))
    self%coils(i+self%n_vcoils)%is_vcoil=.FALSE.
  END DO
ELSE
  CALL xml_get_element(self%xml,"vcoils",coil_element,error_flag)
  IF(error_flag==0)THEN
    IF(self%n_vcoils>0)THEN
      CALL oft_warn("V-coils specified in mesh and XML files, ignoring XML definitions")
    ELSE
      CALL oft_warn("V-coils specified in XML file, but not in mesh file, ignoring XML definitions")
    END IF
  END IF
  CALL xml_get_element(self%xml,"icoils",coil_element,error_flag)
  IF(error_flag==0)THEN
    IF(self%n_icoils>0)THEN
      CALL oft_warn("I-coils specified in mesh and XML files, ignoring XML definitions")
    ELSE
      CALL oft_warn("I-coils specified in XML file, but not in mesh file, ignoring XML definitions")
    END IF
  END IF
END IF
!---Analyze mesh to construct holes
njholes=0
IF(ASSOCIATED(self%junc))njholes=self%junc%nholes
self%nholes=njholes
IF(ASSOCIATED(hole_ns))self%nholes=SIZE(hole_ns)+njholes
WRITE(*,'(2A)')oft_indent,'Creating hole elements'
ALLOCATE(self%hmesh(self%nholes))
ALLOCATE(self%kfh(self%mesh%nc+1))
self%kfh=0
IF(self%nholes>0)THEN
  DO i=1,self%nholes-njholes
    IF(oft_debug_print(1))THEN
      WRITE(*,'(2A,I12)')oft_indent,'  Hole: ',i
      WRITE(*,'(2A,I12)')oft_indent,'    Nset size = ',hole_ns(i)%n
    END IF
    IF(hole_ns(i)%n==1)THEN
      self%hmesh(i)%i0 = hole_ns(i)%v(1)
      CALL get_hole_pseq(self%hmesh(i)%i0, self%hmesh(i)%lp, self%hmesh(i)%n)
    ELSE
      self%hmesh(i)%n=hole_ns(i)%n
      ALLOCATE(self%hmesh(i)%lp(self%hmesh(i)%n))
      CALL order_hole_list(hole_ns(i)%v,self%hmesh(i)%lp,self%hmesh(i)%n)
    END IF
    CALL tw_setup_hole(self%mesh, self%hmesh(i))
    IF(oft_debug_print(1))THEN
      WRITE(*,'(2A,I12)')oft_indent,'    np        = ',self%hmesh(i)%n
      WRITE(*,'(2A,3ES12.4)')oft_indent,'    center    = ',self%hmesh(i)%ptcc
    END IF
    IF(oft_debug_print(3))THEN
      WRITE(*,'(2A)')oft_indent,'    points:'
      CALL oft_increase_indent
      DO j=1,self%hmesh(i)%n
        WRITE(*,'(A,3ES12.4)')oft_indent,self%mesh%r(:,self%hmesh(i)%lp(j))
      END DO
      CALL oft_decrease_indent
    END IF
  END DO
  DO i=1,njholes
    CALL tw_setup_jhole(self%mesh,self%junc%holes(i),self%hmesh(self%nholes-njholes+i))
  END DO
  !---Build cell list
  DO i=1,self%nholes
    DO j=1,self%hmesh(i)%n
      DO k=self%hmesh(i)%kpc(j),self%hmesh(i)%kpc(j+1)-1
        face=ABS(self%hmesh(i)%lpc(k))
        self%kfh(face)=self%kfh(face)+1
      END DO
    END DO
  END DO
  self%nfh=SUM(self%kfh)
  self%kfh(self%mesh%nc+1)=self%nfh+1
  do i=self%mesh%nc,1,-1 ! cumulative count
    self%kfh(i) = self%kfh(i+1) - self%kfh(i)
  end do
  if(self%kfh(1)/=1)call oft_abort('Bad element to element count','tw_setup',__FILE__)
  ALLOCATE(self%lfh(2,self%nfh),kfh_tmp(self%mesh%nc+1))
  self%lfh=0
  kfh_tmp=0
  DO i=1,self%nholes
    DO j=1,self%hmesh(i)%n
      DO k=self%hmesh(i)%kpc(j),self%hmesh(i)%kpc(j+1)-1
        face=ABS(self%hmesh(i)%lpc(k))
        DO l=1,3
          IF(self%mesh%lc(l,face)==self%hmesh(i)%lp(j))EXIT
        END DO
        self%lfh(:,self%kfh(face)+kfh_tmp(face))=[SIGN(i,self%hmesh(i)%lpc(k)),l]
        kfh_tmp(face)=kfh_tmp(face)+1
      END DO
    END DO
  END DO
  DEALLOCATE(kfh_tmp)
  !---Check for duplicates (holes from "OFT_ThinCurr_holes" are independent by construction; junction holes differ
  !---by cell and sign)
  DO i=1,MERGE(0,self%nholes-njholes,self%holes_version>0)
    DO j=1,self%nholes-njholes
      IF(i==j.OR.self%hmesh(i)%n/=self%hmesh(j)%n)CYCLE
      IF(ALL(self%hmesh(i)%lp_sorted==self%hmesh(j)%lp_sorted))THEN
        WRITE(*,*)"Duplicate hole detected: ",i,j
        CALL oft_abort("Duplicate hole detected","tw_setup",__FILE__)
      END IF
    END DO
  END DO
END IF
!---Build mapping from mesh vertices to active vertices
IF(.NOT.ASSOCIATED(self%pmap))THEN
  ALLOCATE(self%pmap(self%mesh%np))
  self%pmap=0
  self%np_active=0
  DO i=1,self%mesh%np
    IF(self%mesh%bp(i))CYCLE
    self%np_active=self%np_active+1
    self%pmap(i)=self%np_active
  END DO
  ! Convert closure cells to vertices and mark
  DO i=1,self%nclosures
    l=-1
    j=1
    DO k=1,3
      IF(self%pmap(self%mesh%lc(k,self%closures(i)))<=0)CYCLE
      IF(self%mesh%kpc(self%mesh%lc(k,self%closures(i))+1)-self%mesh%kpc(self%mesh%lc(k,self%closures(i)))>l)THEN
        l=self%mesh%kpc(self%mesh%lc(k,self%closures(i))+1)-self%mesh%kpc(self%mesh%lc(k,self%closures(i)))
        j=k
      END IF
    END DO
    j=self%mesh%lc(j,self%closures(i))
    IF(self%pmap(j)==0)CALL oft_abort("Error getting closure vertex","tw_setup",__FILE__)
    self%pmap(j)=-i
    self%closures(i)=j
  END DO
  self%nclosures=0 !-self%nclosures
  ! Reindex, removing closure vertices
  self%np_active=0
  DO i=1,self%mesh%np
    IF(self%pmap(i)>0)THEN
      self%np_active=self%np_active+1
      self%pmap(i)=self%np_active
    ELSE
      self%pmap(i)=0
    END IF
  END DO
ELSE
  self%np_active=MAXVAL(self%pmap)
END IF
!---Operator assembly assumes a non-empty model
IF(self%mesh%nc<1.OR.self%np_active<1)THEN
  error_string="Thin-wall model must have at least one cell and one active vertex"
  RETURN
END IF
!
ALLOCATE(self%kpmap_inv(self%np_active+1))
self%kpmap_inv=0
DO i=1,self%mesh%np
  IF(self%pmap(i)/=0)self%kpmap_inv(self%pmap(i))=self%kpmap_inv(self%pmap(i))+1
END DO
self%kpmap_inv(self%np_active+1)=SUM(self%kpmap_inv)+1
DO i=self%np_active,1,-1
  self%kpmap_inv(i)=self%kpmap_inv(i+1)-self%kpmap_inv(i)
END DO
IF(self%kpmap_inv(1)/=1)CALL oft_abort("Invalid inverse point linkage","tw_setup",__FILE__)
ALLOCATE(self%lpmap_inv(self%mesh%np),np_inverse(self%np_active))
np_inverse=0
self%lpmap_inv=0
DO i=1,self%mesh%np
  IF(self%pmap(i)/=0)THEN
    self%lpmap_inv(self%kpmap_inv(self%pmap(i))+np_inverse(self%pmap(i)))=i
    np_inverse(self%pmap(i))=np_inverse(self%pmap(i))+1
  END IF
END DO
DEALLOCATE(np_inverse)
self%nelems = self%np_active + self%nholes + self%n_vcoils
IF(ASSOCIATED(self%junc))THEN
  self%junc%nelems_exp = self%nelems
  CALL tw_setup_junctions(self)
END IF
!---Build basis array
ALLOCATE(self%qbasis(3,3,self%mesh%nc))
f=1.d0/3.d0
DO i=1,self%mesh%nc
  CALL self%mesh%jacobian(i,f,rgop,area_i)
  CALL self%mesh%norm(i,f,norm_i)
  DO j=1,3
    self%qbasis(:,j,i)=cross_product(rgop(:,j),norm_i)
  END DO
END DO
!---Load resistivity
IF(.NOT.ASSOCIATED(self%Eta_surf))CALL tw_load_eta(self)
!---Create local vector
ALLOCATE(oft_native_vector::self%Uloc)
SELECT TYPE(this=>self%Uloc)
  CLASS IS(oft_native_vector)
    this%n=self%nelems; this%ng=self%nelems
    this%nslice=this%n
    ALLOCATE(this%v(this%n))
    ALLOCATE(this%stitch_info)
    ALLOCATE(this%stitch_info%be(this%n))
    this%stitch_info%full=.TRUE.
    this%stitch_info%nbe=0
    this%stitch_info%be=.FALSE.
  CLASS DEFAULT
    CALL oft_abort('Failed to allocate "Uloc" vector','tw_setup',__FILE__)
END SELECT
!---Create point vector
ALLOCATE(oft_native_vector::self%Uloc_pts)
SELECT TYPE(this=>self%Uloc_pts)
  CLASS IS(oft_native_vector)
    this%n=self%mesh%np; this%ng=self%mesh%np
    this%nslice=this%n
    ALLOCATE(this%v(this%n))
    ALLOCATE(this%stitch_info)
    ALLOCATE(this%stitch_info%be(this%n))
    this%stitch_info%full=.TRUE.
    this%stitch_info%nbe=0
    this%stitch_info%be=.FALSE.
  CLASS DEFAULT
    CALL oft_abort('Failed to allocate "Uloc_pts" vector','tw_setup',__FILE__)
END SELECT
!---Print summary
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Setup complete:'
CALL oft_increase_indent
WRITE(*,'(2A,I12)')oft_indent,'# of points    = ',self%mesh%np
WRITE(*,'(2A,I12)')oft_indent,'# of edges     = ',self%mesh%ne
WRITE(*,'(2A,I12)')oft_indent,'# of cells     = ',self%mesh%nc
WRITE(*,'(2A,I12)')oft_indent,'# of holes     = ',self%nholes
WRITE(*,'(2A,I12)')oft_indent,'# of closures  = ',self%nclosures
IF(oft_debug_print(1))THEN
  CALL oft_increase_indent
  DO i=1,self%nclosures
    WRITE(*,'(A,I8,3ES14.5)')oft_indent,self%closures(i),self%mesh%r(:,self%closures(i))
  END DO
  CALL oft_decrease_indent
END IF
WRITE(*,'(2A,I12)')oft_indent,'# of Vcoils    = ',self%n_vcoils
WRITE(*,'(2A,I12)')oft_indent,'# of Icoils    = ',self%n_icoils
CALL oft_decrease_indent
CALL oft_decrease_indent
CONTAINS
!---------------------------------------------------------------------------------
!> Find connected chain for a boundary hole from a starting vertex
!---------------------------------------------------------------------------------
SUBROUTINE get_hole_pseq(i0,plist,n)
INTEGER(4), INTENT(in) :: i0 !< Starting vertex (must be on boundary)
INTEGER(4), POINTER, INTENT(out) :: plist(:) !< List of vertices forming chain
INTEGER(4), INTENT(out) :: n !< Number of vertices in chain
INTEGER(4) :: ii,jj,k,ipt,eprev,ed
INTEGER(4), ALLOCATABLE :: lloop_tmp(:)
!---Find loop points and edges
ALLOCATE(lloop_tmp(self%mesh%nbp))
lloop_tmp=0
!
IF(.NOT.self%mesh%bp(i0))CALL oft_abort('Hole starting vertex is not on boundary', &
  'tw_setup::get_hole_pseq', __FILE__)
ipt=i0
k=1
lloop_tmp(k)=ipt
eprev=0
DO jj=1,self%mesh%nbe
  DO ii=self%mesh%kpe(ipt),self%mesh%kpe(ipt+1)-1
    ed = self%mesh%lpe(ii)
    IF(ed==eprev.OR.(.NOT.self%mesh%be(ed)))CYCLE
    ipt=SUM(self%mesh%le(:,ed))-ipt
    k=k+1
    lloop_tmp(k)=ipt
    eprev=ed
    EXIT
  END DO
  IF(ipt==i0)EXIT
END DO
IF(ipt/=i0)CALL oft_abort('Error building hole mesh, could not find periodic path', &
  'tw_setup::get_hole_pseq', __FILE__)
!---
n=k-1
ALLOCATE(plist(n))
plist=lloop_tmp(1:n)
DEALLOCATE(lloop_tmp)
END SUBROUTINE get_hole_pseq
!---------------------------------------------------------------------------------
!> Reorder hole vertices into a sequential chain
!---------------------------------------------------------------------------------
SUBROUTINE order_hole_list(list_in,list_out,n)
INTEGER(4), INTENT(in) :: n !< Number of points in list
INTEGER(4), INTENT(in) :: list_in(n) !< Input vertex list
INTEGER(4), INTENT(out) :: list_out(n) !< Reordered list
INTEGER(4) :: ii,jj,k,l,nlinks,ipt,eprev,ed,ed2,ptp2,ptp,candidate,candidate2,last_item(3),i0
INTEGER(4), ALLOCATABLE :: lloop_tmp(:),flag_list(:)
!---Check for a valid point list
IF(n<3)CALL oft_abort('Error building hole mesh, hole has fewer than 3 points', &
  'tw_setup::order_hole_list', __FILE__)
ALLOCATE(lloop_tmp(n),flag_list(n))
flag_list=0
lloop_tmp=list_in
CALL sort_array(lloop_tmp, n)
IF(ANY(lloop_tmp(2:n)==lloop_tmp(1:n-1)))CALL oft_abort('Error building hole mesh, hole has repeated points', &
  'tw_setup::order_hole_list', __FILE__)
!---Use the given order if consecutive points (including last -> first) are all joined by mesh edges,
!---starting at the lowest point and going towards its lower neighbor (same as the walk below)
DO jj=1,n
  IF(mesh_local_findedge(self%mesh,[list_in(jj),list_in(MOD(jj,n)+1)])==0)EXIT
END DO
IF(jj>n)THEN
  i0=MINLOC(list_in,DIM=1)
  IF(list_in(MOD(i0,n)+1)<list_in(MOD(i0+n-2,n)+1))THEN
    DO k=1,n
      list_out(k)=list_in(MOD(i0+k-2,n)+1)
    END DO
  ELSE
    DO k=1,n
      list_out(k)=list_in(MOD(i0-k+n,n)+1)
    END DO
  END IF
  DEALLOCATE(lloop_tmp,flag_list)
  RETURN
END IF
IF(self%holes_version>0)CALL oft_abort('Hole points from "OFT_ThinCurr_holes" are not an ordered loop', &
  'tw_setup::order_hole_list', __FILE__)
CALL oft_warn('Hole points are not stored as an ordered loop, using legacy ordering (may fail or flip the hole ' &
  //'direction for loops with chords); recreate holes with "OFT_ThinCurr_holes"')
!---Find loop points and edges
!
! DO i0=1,n
!   nlinks=0
!   DO l=self%mesh%kpe(lloop_tmp(i0)),self%mesh%kpe(lloop_tmp(i0)+1)-1
!     ed2 = self%mesh%lpe(l)
!     ptp2 = SUM(self%mesh%le(:,ed2))-lloop_tmp(i0)
!     candidate2=search_array(ptp2, lloop_tmp, n)
!     IF(candidate2==0)CYCLE
!     nlinks=nlinks+1
!   END DO
!   IF(nlinks==2)EXIT
! END DO
i0=1
flag_list(i0)=1
ipt=lloop_tmp(i0)
k=1
list_out(k)=ipt
eprev=0
DO jj=1,n
  IF(jj==n-2)flag_list(1)=0 ! Unflag first point for last link (needs to be 2 offset for link count check below)
  last_item=-1
  DO ii=self%mesh%kpe(ipt),self%mesh%kpe(ipt+1)-1
    ed = self%mesh%lpe(ii)
    IF(ed==eprev)CYCLE
    ptp = SUM(self%mesh%le(:,ed))-ipt
    candidate=search_array(ptp, lloop_tmp, n)
    IF(candidate==0)THEN
      CYCLE
    ELSE
      IF(flag_list(candidate)==1)CYCLE
      nlinks=0
      DO l=self%mesh%kpe(ptp),self%mesh%kpe(ptp+1)-1
        ed2 = self%mesh%lpe(l)
        ptp2 = SUM(self%mesh%le(:,ed2))-ptp
        candidate2=search_array(ptp2, lloop_tmp, n)
        IF(candidate2==0)CYCLE
        IF(flag_list(candidate2)==1)CYCLE
        nlinks=nlinks+1
      END DO
      last_item=[ptp,candidate,ed]
      IF(nlinks>1)CYCLE
      last_item=-1
    END IF
    flag_list(candidate)=1
    ipt=ptp
    IF(jj<n)THEN
      k=k+1
      list_out(k)=ipt
      eprev=ed
    END IF
    EXIT
  END DO
  IF(last_item(1)>0)THEN
    flag_list(last_item(2))=1
    ipt=last_item(1)
    IF(jj<n)THEN
      k=k+1
      list_out(k)=ipt
      eprev=last_item(3)
    END IF
  END IF
END DO
IF(jj<n+1)CALL oft_abort('Error building hole mesh, unmatched points exist', &
  'tw_setup::order_hole_list', __FILE__)
IF(ipt/=lloop_tmp(i0))CALL oft_abort('Error building hole mesh, path is not periodic', &
  'tw_setup::order_hole_list', __FILE__)
DEALLOCATE(lloop_tmp,flag_list)
END SUBROUTINE order_hole_list
END SUBROUTINE tw_setup
!---------------------------------------------------------------------------------
!> Save debug plotting information for thin-wall model
!---------------------------------------------------------------------------------
SUBROUTINE tw_save_debug(self)
CLASS(tw_type), INTENT(INOUT) :: self !< Thin-wall model object
INTEGER(4) :: i,j,k
REAL(8) :: ftmp(3)
REAl(8), ALLOCATABLE, DIMENSION(:,:) :: normals
CHARACTER(LEN=4) :: plt_tag
!---Save normals
ftmp=1.d0/3.d0
ALLOCATE(normals(3,self%mesh%nc))
DO i=1,self%mesh%nc
  CALL self%mesh%norm(i,ftmp,normals(:,i))
END DO
CALL self%mesh%save_cell_vector(normals,self%xdmf,'Nhat')
!---Save hole info
DO i=1,self%nholes
  normals=0.d0
  DO j=1,self%hmesh(i)%n
    DO k=self%hmesh(i)%kpc(j),self%hmesh(i)%kpc(j+1)-1
      normals(1,ABS(self%hmesh(i)%lpc(k)))=REAL(SIGN(j,self%hmesh(i)%lpc(k)),8)
    END DO
  END DO
  WRITE(plt_tag,'(I4.4)')i
  CALL self%mesh%save_cell_scalar(normals(1,:),self%xdmf,'Ho_'//plt_tag)
END DO
CALL tw_save_hole_debug(self)
DEALLOCATE(normals)
!--Save jumper info
ALLOCATE(normals(self%mesh%np,1))
DO i=1,SIZE(self%jumper_nsets)
  normals=0.d0
  DO j=1,self%jumper_nsets(i)%n
    IF(self%jumper_nsets(i)%v(j)>0)normals(self%jumper_nsets(i)%v(j),1)=1.d0
  END DO
  WRITE(plt_tag,'(I4.4)')i
  CALL self%mesh%save_vertex_scalar(normals(:,1),self%xdmf,'Jumper_'//plt_tag)
END DO
END SUBROUTINE tw_save_debug
!---------------------------------------------------------------------------------
!> Define a circular coil about the Z axis from its radius and height
!!
!! The coil is flagged as circular, so couplings use the exact ring formulas. Its points (traversed
!! counter-clockwise about +Z) are still generated: the polygon is integrated for elements too close to the ring for
!! the exact formula and for mutuals with coils given as points, and it is used for plotting.
!---------------------------------------------------------------------------------
SUBROUTINE tw_coil_set_ring(coil,rz,npts)
TYPE(tw_gen_coil), INTENT(inout) :: coil !< Coil
REAL(r8), INTENT(in) :: rz(2) !< Radius and height
INTEGER(i4), INTENT(in) :: npts !< Number of points for the coil polygon (181 if <= 0)
INTEGER(i4) :: k
REAL(r8) :: theta
coil%is_ring=.TRUE.
coil%ring_rz=rz
coil%npts=npts
IF(coil%npts<=0)coil%npts=181
ALLOCATE(coil%pts(3,coil%npts))
DO k=1,coil%npts
  theta=(k-1)*2.d0*pi/REAL(coil%npts-1,8)
  coil%pts(:,k)=[rz(1)*COS(theta),rz(1)*SIN(theta),rz(2)]
END DO
END SUBROUTINE tw_coil_set_ring
!---------------------------------------------------------------------------------
!> Exact vector potential kernel \f$ \oint d\mathbf{l}'/|\mathbf{r}-\mathbf{r}'| \f$ of a circular coil about the Z axis
!!
!! For a unit current traversed counter-clockwise about +Z: \f$ A_\phi = -4 \pi G(R,Z;R_c,Z_c)/R \f$ (@ref axi_green::green)
!---------------------------------------------------------------------------------
FUNCTION tw_ring_avec(rz,pt) RESULT(avec)
REAL(r8), INTENT(in) :: rz(2) !< Radius and height of the coil
REAL(r8), INTENT(in) :: pt(3) !< Evaluation point
REAL(r8) :: avec(3) !< Vector potential kernel
REAL(r8) :: rp,aphi
avec=0.d0
rp=SQRT(pt(1)**2+pt(2)**2)
IF(rp<1.d-12)RETURN
aphi=-4.d0*pi*green(rp,pt(3),rz(1),rz(2))/rp
avec(1:2)=[-aphi*pt(2)/rp,aphi*pt(1)/rp]
END FUNCTION tw_ring_avec
!---------------------------------------------------------------------------------
!> Exact Biot-Savart kernel \f$ \oint d\mathbf{l}' \times (\mathbf{r}-\mathbf{r}')/|\mathbf{r}-\mathbf{r}'|^3 \f$ of a circular coil
!! about the Z axis
!!
!! For a unit current traversed counter-clockwise about +Z, the curl of @ref tw_ring_avec:
!! \f$ B_R = 4 \pi (\partial G/\partial Z)/R \f$, \f$ B_Z = -4 \pi (\partial G/\partial R)/R \f$
!---------------------------------------------------------------------------------
FUNCTION tw_ring_bvec(rz,pt) RESULT(bvec)
REAL(r8), INTENT(in) :: rz(2) !< Radius and height of the coil
REAL(r8), INTENT(in) :: pt(3) !< Evaluation point
REAL(r8) :: bvec(3) !< Biot-Savart kernel
REAL(r8) :: rp,fg,gg(2),br
rp=SQRT(pt(1)**2+pt(2)**2)
IF(rp<1.d-12)THEN ! On axis
  bvec=[0.d0,0.d0,2.d0*pi*rz(1)**2/(rz(1)**2+(pt(3)-rz(2))**2)**1.5d0]
  RETURN
END IF
CALL grad_green(rp,pt(3),rz(1),rz(2),fg,gg)
br=4.d0*pi*gg(2)/rp
bvec=[br*pt(1)/rp,br*pt(2)/rp,-4.d0*pi*gg(1)/rp]
END FUNCTION tw_ring_bvec
!---------------------------------------------------------------------------------
!> Compute element to driver (Icoils) coupling matrix
!!
!! Circular coils about the Z axis are integrated with the exact vector potential of a ring (elliptic integrals)
!! where they are well separated from an element; otherwise the coil polygon is integrated point by point.
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_Ael2dr_exp(tw_obj,save_file)
TYPE(tw_type), INTENT(inout) :: tw_obj !< Thin-wall model object
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: save_file
LOGICAL :: exists
INTEGER(4) :: i,ii,j,jj,k,kk,ik,ncoils_tot,ih,ihp,ihc,file_counts(4),ierr,io_unit,iquad
INTEGER(i4), ALLOCATABLE :: nrad_cross(:)
REAL(8) :: tmp(3),cvec(3),cpt(3),pt_i(3),evec_i(3,3),pts_i(3,3),elapsed_time
REAL(8) :: pot_tmp,rgop(3,3),area_i,norm_i(3),f(3),dl_min,dl_max,pot_last
REAL(8) :: avec(3),dring
REAL(8), allocatable :: atmp(:,:),Ael2coil_tmp(:,:),ring_rz(:,:,:)
LOGICAL, ALLOCATABLE :: ring_ok(:,:)
CLASS(oft_bmesh), POINTER :: bmesh
TYPE(tw_coil_set), POINTER, DIMENSION(:) :: coils_tot
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
!
IF(PRESENT(save_file))THEN
  IF(TRIM(save_file)/='none')THEN
    INQUIRE(FILE=TRIM(save_file),EXIST=exists)
    IF(exists)THEN
      WRITE(*,'(2A)')oft_indent,'Reading coil mutual matrices'
      CALL oft_increase_indent
      OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
      READ(io_unit, IOSTAT=ierr)file_counts
      IF((ierr/=0).OR.ANY(file_counts/=[tw_nelems_exp(tw_obj),tw_obj%n_vcoils,tw_obj%n_icoils,tw_model_hash(tw_obj)]))THEN
        exists=.FALSE.
        WRITE(*,'(2A)')oft_indent,'Ignoring stored matrix: Model, coils or cache version do not match'
      END IF
      IF(exists)THEN
        ALLOCATE(tw_obj%Ael2coil(tw_nelems_exp(tw_obj),tw_obj%n_vcoils))
        READ(io_unit, IOSTAT=ierr)tw_obj%Ael2coil
        IF(ierr/=0)THEN
          WRITE(*,'(2A)')oft_indent,'Error reading matrix from file'
          DEALLOCATE(tw_obj%Ael2coil)
          exists=.FALSE.
        END IF
      END IF
      IF(exists)THEN
        ALLOCATE(tw_obj%Ael2dr(tw_nelems_exp(tw_obj),tw_obj%n_icoils))
        READ(io_unit, IOSTAT=ierr)tw_obj%Ael2dr
        IF(ierr/=0)THEN
          WRITE(*,'(2A)')oft_indent,'Error reading matrix from file'
          DEALLOCATE(tw_obj%Ael2coil,tw_obj%Ael2dr)
          exists=.FALSE.
        END IF
      END IF
      !---Coil<->coil mutuals.  Without these the cache is unusable for a model with Vcoils:
      IF(exists.AND.(tw_obj%n_vcoils>0))THEN
        IF(ASSOCIATED(tw_obj%Acoil2coil))DEALLOCATE(tw_obj%Acoil2coil)
        ALLOCATE(tw_obj%Acoil2coil(tw_obj%n_vcoils,tw_obj%n_vcoils))
        READ(io_unit, IOSTAT=ierr)tw_obj%Acoil2coil
        IF(ierr/=0)THEN
          WRITE(*,'(2A)')oft_indent,'Stored matrix predates coil<->coil caching, rebuilding'
          DEALLOCATE(tw_obj%Ael2coil,tw_obj%Ael2dr,tw_obj%Acoil2coil)
          exists=.FALSE.
        ELSE
          DO i=1,tw_obj%n_vcoils
            tw_obj%vcoils(i)%Lself=tw_obj%Acoil2coil(i,i)
            WRITE(*,"(2A,1X,I4,A,ES12.4)")oft_indent,"Vcoil",i,": L [H] = ", &
              tw_obj%vcoils(i)%Lself*1.d-7
          END DO
        END IF
      END IF
      CLOSE(io_unit)
      CALL oft_decrease_indent
    END IF
    IF(exists)THEN
      DEBUG_STACK_POP
      RETURN
    END IF
  END IF
END IF
!---Setup quadrature
!---Build temporary coil set from all filaments in model
ncoils_tot=tw_obj%n_vcoils+tw_obj%n_icoils
ALLOCATE(coils_tot(ncoils_tot))
DO i=1,tw_obj%n_vcoils
  CALL tw_copy_coil(tw_obj%vcoils(i),coils_tot(i))
END DO
DO i=1,tw_obj%n_icoils
  CALL tw_copy_coil(tw_obj%icoils(i),coils_tot(i+tw_obj%n_vcoils))
END DO
!
bmesh=>tw_obj%mesh
WRITE(*,'(2A)')oft_indent,'Building coil<->element inductance matrices'
CALL oft_increase_indent
CALL mytimer%tick
ALLOCATE(Ael2coil_tmp(tw_nelems_exp(tw_obj),ncoils_tot),nrad_cross(ncoils_tot))
nrad_cross=0
!---Circular coils (rings about the Z axis): radius, height and orientation
k=0
DO j=1,ncoils_tot
  k=MAX(k,coils_tot(j)%ncoils)
END DO
ALLOCATE(ring_ok(k,ncoils_tot),ring_rz(3,k,ncoils_tot))
ring_ok=.FALSE.
DO j=1,ncoils_tot
  DO k=1,coils_tot(j)%ncoils
    ring_ok(k,j)=coils_tot(j)%coils(k)%is_ring
    ring_rz(:,k,j)=[coils_tot(j)%coils(k)%ring_rz(1),coils_tot(j)%coils(k)%ring_rz(2),1.d0]
  END DO
END DO
Ael2coil_tmp=0.d0
f=1.d0/3.d0
!$omp parallel private(ii,ik,j,jj,k,kk,pts_i,pt_i,cvec, &
!$omp cpt,tmp,atmp,evec_i,pot_tmp,pot_last,i,ih,ihp, &
!$omp ihc,rgop,area_i,norm_i,dl_min,dl_max,iquad,avec,dring)
ALLOCATE(atmp(ncoils_tot,3))
!$omp do
DO i=1,bmesh%nc
  ! CALL bmesh%jacobian(i,f,rgop,area_i)
  ! CALL bmesh%norm(i,f,norm_i)
  area_i=bmesh%ca(i)
  DO ii=1,3
    pts_i(:,ii)=bmesh%r(:,bmesh%lc(ii,i))
    ! evec_i(:,ii)=cross_product(rgop(:,ii),norm_i)
    evec_i(:,ii)=tw_obj%qbasis(:,ii,i)
  END DO
  !---Compute driver contributions
  atmp=0.d0
  DO j=1,ncoils_tot
    DO k=1,coils_tot(j)%ncoils
      !---Ring well separated from the element: integrate its exact vector potential over the element
      IF(ring_ok(k,j))THEN
        dl_min=1.d99
        dl_max=SQRT(area_i*2.d0)
        DO ii=1,3
          dring=SQRT((SQRT(pts_i(1,ii)**2+pts_i(2,ii)**2)-ring_rz(1,k,j))**2+(pts_i(3,ii)-ring_rz(2,k,j))**2)
          dl_min=MIN(dl_min,dring)
          dl_max=MAX(dl_max,dring)
        END DO
        iquad=18
        IF(dl_min>1.d-8)iquad=MAX(4,MIN(18,ABS(INT(LOG(target_err)/LOG(1.d0-dl_min/dl_max)))))
        IF(iquad<=10)THEN
          !---Elements close to the ring use the polygon path below, which counts intersections itself
          IF(dl_min<coils_tot(j)%radius(k))THEN
            !$omp atomic
            nrad_cross(j)=nrad_cross(j)+1
          END IF
          avec=0.d0
          DO ii=1,tw_quads(iquad)%np
            pt_i = tw_quads(iquad)%pts(1,ii)*pts_i(:,1) &
              + tw_quads(iquad)%pts(2,ii)*pts_i(:,2) &
              + tw_quads(iquad)%pts(3,ii)*pts_i(:,3)
            avec=avec+tw_quads(iquad)%wts(ii)*tw_ring_avec(ring_rz(1:2,k,j),pt_i)
          END DO
          avec=avec*area_i*ring_rz(3,k,j)
          DO jj=1,3
            atmp(j,jj)=atmp(j,jj)+coils_tot(j)%scales(k)*DOT_PRODUCT(evec_i(:,jj),avec)
          END DO
          CYCLE
        END IF
      END IF
      tmp=0.d0
      pot_last=0.d0
      DO kk=1,coils_tot(j)%coils(k)%npts
        cpt = coils_tot(j)%coils(k)%pts(:,kk)
        !---Compute minimum separation
        dl_min=1.d99
        dl_max=SQRT(area_i*2.d0)
        DO ii=1,3
          dl_min=MIN(dl_min,SQRT(SUM((pts_i(:,ii)-cpt)**2)))
          dl_max=MAX(dl_max,SQRT(SUM((pts_i(:,ii)-cpt)**2)))
        END DO
        IF(dl_min<coils_tot(j)%radius(k))THEN
          !$omp atomic
          nrad_cross(j)=nrad_cross(j)+1
        END IF
        IF(dl_min<1.d-8)THEN
          iquad = 18
        ELSE
          iquad = MAX(4,MIN(18,ABS(INT(LOG(target_err)/LOG(1.d0-dl_min/dl_max)))))
        END IF
        IF(iquad>10)THEN
          pot_tmp=tw_compute_phipot(pts_i,cpt)
        ELSE
          pot_tmp=0.d0
          !$omp simd private(pt_i) reduction(+:pot_tmp)
          DO ii=1,tw_quads(iquad)%np
            pt_i = tw_quads(iquad)%pts(1,ii)*pts_i(:,1) &
              + tw_quads(iquad)%pts(2,ii)*pts_i(:,2) &
              + tw_quads(iquad)%pts(3,ii)*pts_i(:,3)
            pot_tmp = pot_tmp + tw_quads(iquad)%wts(ii)/SQRT(SUM((pt_i-cpt)**2))
          END DO
          pot_tmp=pot_tmp*area_i
        END IF
        IF(kk>1)THEN
          cvec = coils_tot(j)%coils(k)%pts(:,kk)-coils_tot(j)%coils(k)%pts(:,kk-1)
          DO jj=1,3
            tmp(jj) = tmp(jj) + DOT_PRODUCT(evec_i(:,jj),cvec)*(pot_tmp+pot_last)/2.d0
          END DO
        END IF
        pot_last=pot_tmp
      END DO
      DO jj=1,3
        atmp(j,jj)=atmp(j,jj)+coils_tot(j)%scales(k)*tmp(jj)
      END DO
    END DO
  END DO
  DO ii=1,3
    ik=tw_obj%pmap(bmesh%lc(ii,i))
    IF(ik==0)CYCLE
    DO j=1,ncoils_tot
      !$omp atomic
      Ael2coil_tmp(ik,j) = Ael2coil_tmp(ik,j) + atmp(j,ii)
    END DO
  END DO
  DO ii=tw_obj%kfh(i),tw_obj%kfh(i+1)-1
    ik=ABS(tw_obj%lfh(1,ii))+tw_obj%np_active
    DO j=1,ncoils_tot
      tmp(1)=SIGN(1,tw_obj%lfh(1,ii))*atmp(j,tw_obj%lfh(2,ii))
      !$omp atomic
      Ael2coil_tmp(ik,j) = Ael2coil_tmp(ik,j) + tmp(1)
    END DO
  END DO
END DO
DEALLOCATE(atmp)
!$omp end parallel
IF(ANY(nrad_cross(1:tw_obj%n_vcoils)>0))THEN
  CALL oft_warn("One or more elements intersect a Vcoil within its radius, which may lead to invalid inductance values.")
  DO j=1,tw_obj%n_vcoils
    IF(nrad_cross(j)>0)WRITE(*,'(A,I8,A,I6)')'  ',nrad_cross(j),' intersecting elements for coil ',j
  END DO
END IF
DEALLOCATE(nrad_cross)
!---Unpack passive and driver coils
ALLOCATE(tw_obj%Ael2coil(tw_nelems_exp(tw_obj),tw_obj%n_vcoils))
DO i=1,tw_obj%n_vcoils
  tw_obj%Ael2coil(:,i)=Ael2coil_tmp(:,i)
END DO
ALLOCATE(tw_obj%Ael2dr(tw_nelems_exp(tw_obj),tw_obj%n_icoils))
DO i=1,tw_obj%n_icoils
  tw_obj%Ael2dr(:,i)=Ael2coil_tmp(:,i+tw_obj%n_vcoils)
END DO
!---Cleanup temporaries
DEALLOCATE(Ael2coil_tmp,coils_tot,ring_ok,ring_rz)
elapsed_time=mytimer%tock()
WRITE(*,'(3A)')oft_indent,'Time = ',time_to_string(elapsed_time)
CALL oft_decrease_indent
!
CALL tw_compute_Lmat_coils(tw_obj)
!
IF(PRESENT(save_file))THEN
  IF(TRIM(save_file)/='none')THEN
    OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
    WRITE(io_unit)tw_nelems_exp(tw_obj),tw_obj%n_vcoils,tw_obj%n_icoils,tw_model_hash(tw_obj)
    WRITE(io_unit)tw_obj%Ael2coil
    WRITE(io_unit)tw_obj%Ael2dr
    IF(tw_obj%n_vcoils>0)WRITE(io_unit)tw_obj%Acoil2coil
    CLOSE(io_unit)
  END IF
END IF
DEBUG_STACK_POP
END SUBROUTINE tw_compute_Ael2dr_exp
!---------------------------------------------------------------------------------
!> Compute coupling from thin-wall model elements to flux loop sensors
!!
!! @note The asymptotic form of self-inductance for thin circular coils derived
!! by Hurwitz and Landreman [arXiv:2310.09313 (2023)] is used for all inductance
!! calculations. Note that this is not strictly valid for mutual inductances,
!! but is instead used to avoid integration challenges with very-closely-spaced coils.
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_Lmat_coils(tw_obj)
TYPE(tw_type), INTENT(inout) :: tw_obj !< Thin-wall model object
LOGICAL :: exists
INTEGER(4) :: i,ii,j,jj,k,kk,l,ik,jk,ih,ihp,ihc,file_counts(4),ierr,io_unit,iquad
REAL(8) :: tmp,dl,pt_i(3),pt_j(3),evec_i(3,3),evec_j(3,3),pts_i(3,3),pt_i_last(3)
REAL(8) :: rvec_i(3),r1,rmag,rvec_j(3),z1,cvec(3),cpt(3),coil_thickness,sqrt_e
REAL(8) :: rgop(3,3),norm_i(3),area_i,dl_min,dl_max,pot_last,pot_tmp,elapsed_time,tmp_seg
REAL(8), allocatable :: atmp(:,:),Acoil2sen_tmp(:,:),Acoil2coil_tmp(:,:),ring_rz(:,:,:)
LOGICAL, ALLOCATABLE :: ring_ok(:,:)
CLASS(oft_bmesh), POINTER :: bmesh
INTEGER(4) :: ncoils_tot
TYPE(tw_coil_set), POINTER, DIMENSION(:) :: coils_tot
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
!
!---Build temporary coil set from all filaments in model
ncoils_tot=tw_obj%n_vcoils
IF(ASSOCIATED(tw_obj%Ael2dr))ncoils_tot=ncoils_tot+tw_obj%n_icoils
ALLOCATE(coils_tot(ncoils_tot))
DO i=1,tw_obj%n_vcoils
  CALL tw_copy_coil(tw_obj%vcoils(i),coils_tot(i))
END DO
IF(ASSOCIATED(tw_obj%Ael2dr))THEN
  DO i=1,tw_obj%n_icoils
    CALL tw_copy_coil(tw_obj%icoils(i),coils_tot(i+tw_obj%n_vcoils))
  END DO
END IF
bmesh=>tw_obj%mesh
!---Circular coils (rings about the Z axis): radius, height and orientation
k=0
DO j=1,ncoils_tot
  k=MAX(k,coils_tot(j)%ncoils)
END DO
ALLOCATE(ring_ok(k,ncoils_tot),ring_rz(3,k,ncoils_tot))
ring_ok=.FALSE.
DO j=1,ncoils_tot
  DO k=1,coils_tot(j)%ncoils
    ring_ok(k,j)=coils_tot(j)%coils(k)%is_ring
    ring_rz(:,k,j)=[coils_tot(j)%coils(k)%ring_rz(1),coils_tot(j)%coils(k)%ring_rz(2),1.d0]
  END DO
END DO
!---Compute coupling between coils
ALLOCATE(Acoil2coil_tmp(tw_obj%n_vcoils,ncoils_tot))
Acoil2coil_tmp=0.d0
IF((tw_obj%n_vcoils>0).AND.(ncoils_tot>0))THEN
  WRITE(*,'(2A)')oft_indent,'Building coil<->coil inductance matrix'
  CALL oft_increase_indent
  CALL mytimer%tick()
  sqrt_e=SQRT(EXP(1.d0))
  !$omp parallel private(i,ii,j,k,kk,tmp,tmp_seg,pt_i,rvec_i,atmp,cvec,cpt,coil_thickness,pot_tmp,pot_last,pt_i_last)
  ALLOCATE(atmp(ncoils_tot,1))
  !$omp do
  DO l=1,tw_obj%n_vcoils
    atmp=0.d0
    DO i=1,coils_tot(l)%ncoils
      coil_thickness=(coils_tot(l)%radius(i)**2)/sqrt_e
      DO j=1,ncoils_tot
        DO k=1,coils_tot(j)%ncoils
          IF(ring_ok(i,l).AND.ring_ok(k,j))THEN
            !---Two coaxial rings: the softened kernel 1/sqrt(|r-r'|^2+t) is the exact ring-ring mutual with the
            !---height difference dz replaced by sqrt(dz^2+t), so sum dl.dl'/... = -8*pi^2*green
            tmp=-8.d0*pi*pi*green(ring_rz(1,i,l),SQRT((ring_rz(2,i,l)-ring_rz(2,k,j))**2+coil_thickness), &
              ring_rz(1,k,j),0.d0)*ring_rz(3,i,l)*ring_rz(3,k,j)
          ELSE
            tmp=0.d0
            pt_i_last = coils_tot(l)%coils(i)%pts(:,1)
            DO ii=2,coils_tot(l)%coils(i)%npts
              rvec_i = coils_tot(l)%coils(i)%pts(:,ii)-coils_tot(l)%coils(i)%pts(:,ii-1)
              pt_i = coils_tot(l)%coils(i)%pts(:,ii)
              pot_last=0.d0
              tmp_seg=0.d0
              DO kk=1,coils_tot(j)%coils(k)%npts
                cpt = coils_tot(j)%coils(k)%pts(:,kk)
                pot_tmp = (1.d0/SQRT(SUM((pt_i-cpt)**2) + coil_thickness) + 1.d0/SQRT(SUM((pt_i_last-cpt)**2) + coil_thickness))/2.d0
                IF(kk>1)THEN
                  cvec = coils_tot(j)%coils(k)%pts(:,kk)-coils_tot(j)%coils(k)%pts(:,kk-1)
                  tmp_seg = tmp_seg + DOT_PRODUCT(rvec_i,cvec)*(pot_tmp+pot_last)/2.d0
                END IF
                pot_last = pot_tmp
              END DO
              tmp = tmp + tmp_seg
              pt_i_last = pt_i
            END DO
          END IF
          atmp(j,1)=atmp(j,1)+coils_tot(j)%scales(k)*coils_tot(l)%scales(i)*tmp
        END DO
      END DO
    END DO
    DO j=1,ncoils_tot
      Acoil2coil_tmp(l,j) = atmp(j,1)
    END DO
  END DO
  DEALLOCATE(atmp)
  !$omp end parallel
  elapsed_time=mytimer%tock()
  WRITE(*,'(3A)')oft_indent,'Time = ',time_to_string(elapsed_time)
END IF
DEALLOCATE(coils_tot,ring_ok,ring_rz)
!---Unpack passive and driver coils
IF(ASSOCIATED(tw_obj%Acoil2coil))DEALLOCATE(tw_obj%Acoil2coil)
ALLOCATE(tw_obj%Acoil2coil(tw_obj%n_vcoils,tw_obj%n_vcoils))
DO i=1,tw_obj%n_vcoils
  tw_obj%Acoil2coil(:,i)=Acoil2coil_tmp(:,i)
  tw_obj%vcoils(i)%Lself=Acoil2coil_tmp(i,i)
  WRITE(*,"(2A,1X,I4,A,ES12.4)")oft_indent,"Vcoil",i,": L [H] = ",tw_obj%vcoils(i)%Lself*1.d-7
END DO
!---Compute coupling between elements and drivers
IF(ASSOCIATED(tw_obj%Ael2dr))THEN
  ! WRITE(*,*)'Building driver->element inductance matrices'
  IF(tw_obj%n_vcoils>0)THEN
    !$omp parallel do private(jj)
    DO i=1,tw_obj%n_vcoils
      DO jj=1,tw_obj%n_icoils
        tw_obj%Ael2dr(tw_obj%np_active+tw_obj%nholes+i,jj) = Acoil2coil_tmp(i,jj+tw_obj%n_vcoils)
      END DO
    END DO
  END IF
  tw_obj%Ael2dr = tw_obj%Ael2dr*mu0/(4.d0*pi)
END IF
DEALLOCATE(Acoil2coil_tmp)
!
CALL oft_decrease_indent
DEBUG_STACK_POP
END SUBROUTINE tw_compute_Lmat_coils
!---------------------------------------------------------------------------------
!> Compute mutual inductance matrix between two thin-wall models
!!
!! For models with junctions each element pair contribution is scattered through the junction maps as it is
!! computed, giving the reduced matrix \f$ P_c^T L P_r \f$ without storing the matrix over expanded elements.
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_LmatDirect(row_model,Lmat,col_model,save_file)
TYPE(tw_type), TARGET, INTENT(in) :: row_model !< Thin-wall model object for rows
REAL(8), CONTIGUOUS, POINTER, INTENT(inout) :: Lmat(:,:) !< Mutual inductance matrix [col%nelems,row%nelems]
TYPE(tw_type), TARGET, OPTIONAL, INTENT(in) :: col_model !< Thin-wall model object for columns
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: save_file !< Cache file
INTEGER(4) :: i,ii,j,jj,k,ik,jk,imin,jmax,io_unit,ierr,nsurf,nhash,nbuf,nbuf_max
INTEGER(4) :: file_counts(8),hash_tmp(8)
INTEGER(4), ALLOCATABLE :: rkr(:),rlc(:),ckr(:),clc(:),ebuf(:,:)
INTEGER(8) :: counts(4)
REAL(8) :: tmp,dl_min,f(3),elapsed_time
REAL(8) :: area_i,evec_i(3,3),pts_i(3,3)
REAL(8) :: area_j,evec_j(3,3),pts_j(3,3)
REAL(8), ALLOCATABLE :: rval(:),cval(:),vbuf(:)
REAL(8), CONTIGUOUS, POINTER :: Ael2coil(:,:)
LOGICAL :: exists,Lself,jscatter
CLASS(oft_bmesh), POINTER :: rmesh,cmesh
TYPE(tw_type), POINTER :: row_obj,col_obj
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
row_obj=>row_model
IF(PRESENT(col_model))THEN
  Lself=.FALSE.
  col_obj=>col_model
ELSE
  Lself=.TRUE.
  col_obj=>row_obj
END IF
jscatter=ASSOCIATED(row_obj%junc).OR.ASSOCIATED(col_obj%junc)
IF(ASSOCIATED(Lmat))THEN
  IF(SIZE(Lmat,DIM=1)/=col_obj%nelems.OR.SIZE(Lmat,DIM=2)/=row_obj%nelems)THEN
    CALL oft_abort("Input size incorrect","tw_compute_LmatDirect",__FILE__)
  END IF
ELSE
  ALLOCATE(Lmat(col_obj%nelems,row_obj%nelems))
END IF
!
IF(PRESENT(save_file))THEN
  !---Cache key: element counts and mesh hashes, plus the junction maps (marked by a negative count) if present
  !---Self inductance depends on the V-coils only (no I-coils or resistances); mutual inductance on no coils
  nhash=6
  hash_tmp(3) = oft_simple_hash(C_LOC(row_obj%mesh%lc),INT(4*3*row_obj%mesh%nc,8), &
    starting_value=MERGE(tw_coil_hash(row_obj,vcoils_only=.TRUE.),tw_cache_ver,Lself))
  hash_tmp(4) = oft_simple_hash(C_LOC(col_obj%mesh%lc),INT(4*3*col_obj%mesh%nc,8))
  hash_tmp(5) = oft_simple_hash(C_LOC(row_obj%mesh%r),INT(8*3*row_obj%mesh%np,8))
  hash_tmp(6) = oft_simple_hash(C_LOC(col_obj%mesh%r),INT(8*3*col_obj%mesh%np,8))
  IF(Lself)THEN
    hash_tmp(1) = row_obj%nelems
    hash_tmp(2) = row_obj%mesh%nc
    IF(jscatter)THEN
      nhash=7
      hash_tmp(1) = -hash_tmp(1)
      hash_tmp(7) = tw_junction_hash(row_obj)
    END IF
  ELSE
    hash_tmp(1) = col_obj%nelems
    hash_tmp(2) = row_obj%nelems
    IF(jscatter)THEN
      nhash=8
      hash_tmp(1) = -hash_tmp(1)
      hash_tmp(7) = tw_junction_hash(row_obj)
      hash_tmp(8) = tw_junction_hash(col_obj)
    END IF
  END IF
  IF(TRIM(save_file)/='none')THEN
    INQUIRE(FILE=TRIM(save_file),EXIST=exists)
    IF(exists)THEN
      IF(Lself)THEN
        WRITE(*,'(2A)')oft_indent,'Reading element<->element self inductance matrix'
        CALL oft_increase_indent
        OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
        READ(io_unit, IOSTAT=ierr)file_counts(1:nhash)
        IF((ierr==0).AND.ALL(file_counts(1:nhash)==hash_tmp(1:nhash)))THEN
          DO i=1,row_obj%nelems
            READ(io_unit, IOSTAT=ierr)Lmat(i,i:row_obj%nelems)
            IF(ierr/=0)THEN
              WRITE(*,'(2A)')oft_indent,'Error reading matrix from file'
              exists=.FALSE.
              EXIT
            END IF
          END DO
          CLOSE(io_unit)
          IF(exists)THEN
            DO i=1,row_obj%nelems
              DO j=i+1,row_obj%nelems
                Lmat(j,i)=Lmat(i,j)
              END DO
            END DO
          END IF
        ELSE
          WRITE(*,'(2A)')oft_indent,'Ignoring stored matrix: Model hashes do not match'
          CLOSE(io_unit)
          exists=.FALSE.
        END IF
        CALL oft_decrease_indent
      ELSE
        WRITE(*,'(2A)')oft_indent,'Reading element<->element mutual inductance matrix'
        CALL oft_increase_indent
        OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
        READ(io_unit, IOSTAT=ierr)file_counts(1:nhash)
        IF((ierr==0).AND.ALL(file_counts(1:nhash)==hash_tmp(1:nhash)))THEN
          DO i=1,row_obj%nelems
            READ(io_unit, IOSTAT=ierr)Lmat(:,i)
            IF(ierr/=0)THEN
              WRITE(*,'(2A)')oft_indent,'Error reading matrix from file'
              exists=.FALSE.
              EXIT
            END IF
          END DO
          CLOSE(io_unit)
        ELSE
          WRITE(*,'(2A)')oft_indent,'Ignoring stored matrix: Model hashes do not match'
          CLOSE(io_unit)
          exists=.FALSE.
        END IF
        CALL oft_decrease_indent
      END IF
    END IF
    IF(exists)THEN
      DEBUG_STACK_POP
      RETURN
    END IF
  END IF
END IF
!---Contributions of each cell pair are buffered, then added together (size of buffer bounded by the vertex
!---and hole counts of the two cells)
nbuf_max=(3+MAXVAL(row_obj%kfh(2:)-row_obj%kfh(1:row_obj%mesh%nc)))* &
  (3+MAXVAL(col_obj%kfh(2:)-col_obj%kfh(1:col_obj%mesh%nc)))
!---Junction maps (identity for a model without junctions)
IF(jscatter)THEN
  CALL tw_junction_map(row_obj,rkr,rlc,rval)
  CALL tw_junction_map(col_obj,ckr,clc,cval)
END IF
!---Element to V-coil couplings, with rows reduced by the junction map (V-coil rows of P are the identity)
!---(only computed, and only needed, when the model has V-coils)
NULLIFY(Ael2coil)
IF(Lself.AND.row_obj%n_vcoils>0)THEN
  IF(ASSOCIATED(row_obj%junc))THEN
    ALLOCATE(Ael2coil(tw_nelems_exp(row_obj),row_obj%n_vcoils))
    Ael2coil=row_obj%Ael2coil
    CALL tw_junction_reduce_rows(row_obj,Ael2coil)
  ELSE
    Ael2coil=>row_obj%Ael2coil
  END IF
END IF
!
IF(Lself)THEN
  WRITE(*,'(2A)')oft_indent,'Building element<->element self inductance matrix'
ELSE
  WRITE(*,'(2A)')oft_indent,'Building element<->element mutual inductance matrix'
END IF
CALL oft_increase_indent
Lmat=0.d0
CALL mytimer%tick
!---Setup quadrature
!
rmesh=>row_obj%mesh
cmesh=>col_obj%mesh
!---
counts=0
f=1.d0/3.d0
!$omp parallel private(ii,j,jj,ik,jk,pts_i,pts_j,tmp, &
!$omp evec_i,evec_j,dl_min,i,jmax,imin, &
!$omp area_i,area_j,k,nbuf,ebuf,vbuf) reduction(+:counts)
ALLOCATE(ebuf(2,nbuf_max),vbuf(nbuf_max))
!$omp do schedule(dynamic,100)
DO i=1,rmesh%nc
  ! CALL rmesh%jacobian(i,f,rgop,area_i)
  ! CALL rmesh%norm(i,f,norm_i)
  area_i=rmesh%ca(i)
  DO ii=1,3
    pts_i(:,ii)=rmesh%r(:,rmesh%lc(ii,i))
    ! evec_i(:,ii)=cross_product(rgop(:,ii),norm_i)
    evec_i(:,ii)=row_obj%qbasis(:,ii,i)
  END DO
  IF(Lself)THEN
    imin=MINVAL(row_obj%pmap(rmesh%lc(:,i)))
    DO ii=row_obj%kfh(i),row_obj%kfh(i+1)-1
      imin=MIN(imin,ABS(row_obj%lfh(1,ii))+row_obj%np_active)
    END DO
  END IF
  !---Compute inter-edge inductances
  DO j=1,cmesh%nc
    IF(Lself)THEN
      jmax=MAXVAL(col_obj%pmap(cmesh%lc(:,j)))
      DO jj=col_obj%kfh(j),col_obj%kfh(j+1)-1
        jmax=MAX(jmax,ABS(col_obj%lfh(1,jj))+col_obj%np_active)
      END DO
      IF(jmax<imin)CYCLE
    END IF
    ! CALL cmesh%jacobian(j,f,cgop,area_j)
    ! CALL cmesh%norm(j,f,norm_j)
    area_j=cmesh%ca(j)
    DO jj=1,3
      pts_j(:,jj)=cmesh%r(:,cmesh%lc(jj,j))
      ! evec_j(:,jj)=cross_product(cgop(:,jj),norm_j)
      evec_j(:,jj)=col_obj%qbasis(:,jj,j)
    END DO
    !---Cell-pair integral of 1/|r-r'|
    tmp=tw_cell_pair_kernel(pts_i,area_i,pts_j,area_j)
    !
    nbuf=0
    DO ii=1,3
      ik=row_obj%pmap(rmesh%lc(ii,i))
      IF(ik==0)CYCLE
      DO jj=1,3
        jk=col_obj%pmap(cmesh%lc(jj,j))
        IF((Lself.AND.jk<ik).OR.jk<=0)CYCLE
        dl_min = DOT_PRODUCT(evec_i(:,ii),evec_j(:,jj))*tmp
        nbuf=nbuf+1
        ebuf(:,nbuf)=[jk,ik]
        vbuf(nbuf)=dl_min
      END DO
      DO jj=col_obj%kfh(j),col_obj%kfh(j+1)-1
        jk=ABS(col_obj%lfh(1,jj))+col_obj%np_active
        dl_min = SIGN(1,col_obj%lfh(1,jj))* &
          DOT_PRODUCT(evec_i(:,ii),evec_j(:,col_obj%lfh(2,jj)))*tmp
        nbuf=nbuf+1
        ebuf(:,nbuf)=[jk,ik]
        vbuf(nbuf)=dl_min
      END DO
    END DO
    DO ii=row_obj%kfh(i),row_obj%kfh(i+1)-1
      ik=ABS(row_obj%lfh(1,ii))+row_obj%np_active
      DO jj=col_obj%kfh(j),col_obj%kfh(j+1)-1
        jk=ABS(col_obj%lfh(1,jj))+col_obj%np_active
        IF(Lself.AND.jk<ik)CYCLE
        dl_min = SIGN(1,row_obj%lfh(1,ii))*SIGN(1,col_obj%lfh(1,jj))* &
          DOT_PRODUCT(evec_i(:,row_obj%lfh(2,ii)),evec_j(:,col_obj%lfh(2,jj)))*tmp
        nbuf=nbuf+1
        ebuf(:,nbuf)=[jk,ik]
        vbuf(nbuf)=dl_min
      END DO
      !
      IF(.NOT.Lself)THEN
        DO jj=1,3
          jk=col_obj%pmap(cmesh%lc(jj,j))
          IF(jk<=0)CYCLE
          dl_min = SIGN(1,row_obj%lfh(1,ii))* &
            DOT_PRODUCT(evec_i(:,row_obj%lfh(2,ii)),evec_j(:,jj))*tmp
          nbuf=nbuf+1
          ebuf(:,nbuf)=[jk,ik]
          vbuf(nbuf)=dl_min
        END DO
      END IF
    END DO
    !---Add this cell pair's contributions, through the junction maps if present
    IF(jscatter)THEN
      CALL tw_junction_scatter_add(Lmat,ckr,clc,cval,rkr,rlc,rval,nbuf,ebuf,vbuf,Lself)
    ELSE
      DO k=1,nbuf
        !$omp atomic
        Lmat(ebuf(1,k),ebuf(2,k)) = Lmat(ebuf(1,k),ebuf(2,k)) + vbuf(k)
      END DO
    END IF
  END DO
END DO
!$omp end do nowait
DEALLOCATE(ebuf,vbuf)
IF(Lself)THEN
  !---Add passive coils to model
  nsurf=row_obj%nelems-row_obj%n_vcoils
  !$omp do
  DO i=1,nsurf
    DO j=1,row_obj%n_vcoils
      jj=nsurf+j
      Lmat(jj,i) = Ael2coil(i,j)
    END DO
  END DO
  !$omp end do nowait
  !$omp do
  DO i=1,row_obj%n_vcoils
    ii=nsurf+i
    DO j=i,row_obj%n_vcoils
      jj=nsurf+j
      Lmat(jj,ii) = row_obj%Acoil2coil(i,j)
    END DO
  END DO
  !$omp do private(j)
  DO i=1,row_obj%nelems
    DO j=1,i-1
      Lmat(j,i)=Lmat(i,j)
    END DO
  END DO
END IF
!$omp end parallel
Lmat = Lmat/(4.d0*pi)
elapsed_time=mytimer%tock()
WRITE(*,'(3A)')oft_indent,'Time = ',time_to_string(elapsed_time)
IF(jscatter)DEALLOCATE(rkr,rlc,rval,ckr,clc,cval)
IF(Lself.AND.row_obj%n_vcoils>0.AND.ASSOCIATED(row_obj%junc))DEALLOCATE(Ael2coil)
IF(PRESENT(save_file))THEN
  IF(TRIM(save_file)/='none')THEN
    OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
    WRITE(io_unit)hash_tmp(1:nhash)
    IF(Lself)THEN
      DO i=1,row_obj%nelems
        WRITE(io_unit)Lmat(i,i:row_obj%nelems)
      END DO
    ELSE
      DO i=1,row_obj%nelems
        WRITE(io_unit)Lmat(:,i)
      END DO
    END IF
    CLOSE(io_unit)
  END IF
END IF
CALL oft_decrease_indent
DEBUG_STACK_POP
END SUBROUTINE tw_compute_LmatDirect
!---------------------------------------------------------------------------------
!> Compute mutual inductance matrix between two thin-wall models
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_Lmat_MF_exp(row_obj,col_obj,nrhs,a,b)
TYPE(tw_type), INTENT(in) :: row_obj !< Thin-wall model object for rows
TYPE(tw_type), INTENT(in) :: col_obj !< Thin-wall model object for columns
INTEGER(4), INTENT(in) :: nrhs
REAL(8), CONTIGUOUS, POINTER, INTENT(in) :: a(:,:) !< Mutual inductance matrix
REAL(8), CONTIGUOUS, POINTER, INTENT(inout) :: b(:,:) !< Mutual inductance matrix
INTEGER(4) :: i,ii,j,jj,ik,jk,imin,jmax,irhs
INTEGER(8) :: counts(4)
REAL(8) :: tmp,dl_min,dl_max,f(3),elapsed_time
REAL(8) :: pt_i(3),rgop(3,3),area_i,norm_i(3),evec_i(3,3),pts_i(3,3)
REAL(8) :: pt_j(3),cgop(3,3),area_j,norm_j(3),evec_j(3,3),pts_j(3,3)
LOGICAL :: vvclose_flag, close_flag
CLASS(oft_bmesh), POINTER :: rmesh,cmesh
TYPE(oft_quad_type) :: quad,quad2
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
CALL mytimer%tick
!---Setup quadrature
CALL set_quad_2d(quad,quad_orders(2))
CALL set_quad_2d(quad2,quad_orders(3))
!
rmesh=>row_obj%mesh
cmesh=>col_obj%mesh
!---
WRITE(*,'(2A)')oft_indent,'Applying MF element<->element inductance matrix'
CALL oft_increase_indent
b=0.d0
counts=0
f=1.d0/3.d0
!$omp parallel private(ii,j,jj,ik,jk,pts_i,pts_j,tmp,irhs, &
!$omp pt_i,pt_j,evec_i,evec_j,dl_min,dl_max,i,jmax,imin, &
!$omp rgop,area_i,norm_i,cgop,area_j,norm_j,close_flag,vvclose_flag) reduction(+:counts)
!$omp do schedule(dynamic,100)
DO i=1,rmesh%nc
  ! CALL rmesh%jacobian(i,f,rgop,area_i)
  ! CALL rmesh%norm(i,f,norm_i)
  area_i=rmesh%ca(i)
  DO ii=1,3
    pts_i(:,ii)=rmesh%r(:,rmesh%lc(ii,i))
    ! evec_i(:,ii)=cross_product(rgop(:,ii),norm_i)
    evec_i(:,ii)=row_obj%qbasis(:,ii,i)
  END DO
  !---Compute inter-edge inductances
  DO j=1,cmesh%nc
    ! CALL cmesh%jacobian(j,f,cgop,area_j)
    ! CALL cmesh%norm(j,f,norm_j)
    area_j=cmesh%ca(j)
    DO jj=1,3
      pts_j(:,jj)=cmesh%r(:,cmesh%lc(jj,j))
      ! evec_j(:,jj)=cross_product(cgop(:,jj),norm_j)
      evec_j(:,jj)=col_obj%qbasis(:,jj,j)
    END DO
    !---Compute minimum separation
    close_flag=.FALSE.
    vvclose_flag=.FALSE.
    dl_max=-1.d99
    !!$omp simd collapse(1) reduction(max:dl_max) reduction(min:dl_min)
    DO ii=1,3
      pt_j=pts_j(:,1)-pts_i(:,ii)
      tmp=magnitude(pt_j)
      IF(tmp<1.d-10)THEN
        close_flag=.TRUE.
        EXIT
      END IF
      pt_j=pt_j/tmp
      pt_i=pts_j(:,2)-pts_i(:,ii)
      tmp=magnitude(pt_i)
      IF(tmp<1.d-10)THEN
        close_flag=.TRUE.
        EXIT
      END IF
      pt_i=pt_i/tmp
      dl_max=MAX(dl_max,ABS(ACOS(DOT_PRODUCT(pt_j,pt_i))))
      !
      pt_i=pts_j(:,3)-pts_i(:,ii)
      tmp=magnitude(pt_i)
      IF(tmp<1.d-10)THEN
        close_flag=.TRUE.
        EXIT
      END IF
      pt_i=pt_i/tmp
      dl_max=MAX(dl_max,ABS(ACOS(DOT_PRODUCT(pt_j,pt_i))))
    END DO
    IF(dl_max>pi/8.d0)THEN
      close_flag=.TRUE.
      IF(dl_max>pi/4.d0)vvclose_flag=.TRUE.
    ELSE
      dl_min=1.d99
      dl_max=-1.d99
      !!$omp simd collapse(1) reduction(max:dl_max) reduction(min:dl_min)
      DO ii=1,3
        DO jj=1,3
          dl_min=MIN(dl_min,SQRT(SUM((pts_i(:,ii)-pts_j(:,jj))**2)))
          dl_max=MAX(dl_max,SQRT(SUM((pts_i(:,ii)-pts_j(:,jj))**2)))
        END DO
      END DO
      IF(dl_min/dl_max<quad_tols(2))close_flag=.TRUE.
      IF(dl_min/dl_max<quad_tols(1))vvclose_flag=.TRUE.
    END IF
    !---Chose quadrature order based on distance
    tmp=0.d0
    IF(close_flag)THEN
      ! counts(1)=counts(1)+1
      !IF(i==j.OR.ANY(cmesh%lcc(:,j)==i).OR.vvclose_flag)THEN
      IF(vvclose_flag)THEN
        DO jj=1,quad%np
          pt_j = quad%pts(1,jj)*pts_j(:,1) &
            + quad%pts(2,jj)*pts_j(:,2) &
            + quad%pts(3,jj)*pts_j(:,3)
          tmp = tmp + quad%wts(jj)*tw_compute_phipot(pts_i,pt_j)
        END DO
        tmp=tmp*area_j
      ELSE
        !!$omp simd collapse(1) private(pt_i,pt_j) reduction(+:tmp)
        DO ii=1,quad%np
          pt_i = quad%pts(1,ii)*pts_i(:,1) &
            + quad%pts(2,ii)*pts_i(:,2) &
            + quad%pts(3,ii)*pts_i(:,3)
          DO jj=1,quad%np
            pt_j = quad%pts(1,jj)*pts_j(:,1) &
              + quad%pts(2,jj)*pts_j(:,2) &
              + quad%pts(3,jj)*pts_j(:,3)
            tmp = tmp + quad%wts(jj)*quad%wts(ii)/SQRT(SUM((pt_i-pt_j)**2))
          END DO
        END DO
        tmp=tmp*area_i*area_j
      END IF
    ELSE !IF(dl_min/dl_max<quad_tols(2))THEN
      ! counts(2)=counts(2)+1
      !!$omp simd collapse(1) private(pt_i,pt_j) reduction(+:tmp)
      DO ii=1,quad2%np
        pt_i = quad2%pts(1,ii)*pts_i(:,1) &
          + quad2%pts(2,ii)*pts_i(:,2) &
          + quad2%pts(3,ii)*pts_i(:,3)
        DO jj=1,quad2%np
          pt_j = quad2%pts(1,jj)*pts_j(:,1) &
            + quad2%pts(2,jj)*pts_j(:,2) &
            + quad2%pts(3,jj)*pts_j(:,3)
          tmp = tmp + quad2%wts(jj)*quad2%wts(ii)/SQRT(SUM((pt_i-pt_j)**2))
        END DO
      END DO
      tmp=tmp*area_i*area_j
    ! ELSE
    !   ! counts(3)=counts(3)+1
    !   pt_i = SUM(pts_i,DIM=2)/3.d0
    !   pt_j = SUM(pts_j,DIM=2)/3.d0
    !   tmp=area_i*area_j/SQRT(SUM((pt_i-pt_j)**2))
    END IF
    !
    DO ii=1,3
      ik=row_obj%pmap(rmesh%lc(ii,i))
      IF(ik==0)CYCLE
      DO jj=1,3
        jk=col_obj%pmap(cmesh%lc(jj,j))
        ! IF((row_obj%nelems==col_obj%nelems).AND.jk<ik.OR.jk<=0)CYCLE
        IF(jk==0)CYCLE
        dl_min = DOT_PRODUCT(evec_i(:,ii),evec_j(:,jj))*tmp
        DO irhs=1,nrhs
          !$omp atomic
          b(jk,irhs) = b(jk,irhs) + dl_min*a(ik,irhs)
        END DO
      END DO
      DO jj=col_obj%kfh(j),col_obj%kfh(j+1)-1
        jk=ABS(col_obj%lfh(1,jj))+col_obj%np_active
        dl_min = SIGN(1,col_obj%lfh(1,jj))* &
          DOT_PRODUCT(evec_i(:,ii),evec_j(:,col_obj%lfh(2,jj)))*tmp
        DO irhs=1,nrhs
          !$omp atomic
          b(jk,irhs) = b(jk,irhs) + dl_min*a(ik,irhs)
        END DO
      END DO
    END DO
    DO ii=row_obj%kfh(i),row_obj%kfh(i+1)-1
      ik=ABS(row_obj%lfh(1,ii))+row_obj%np_active
      DO jj=1,3
        jk=col_obj%pmap(cmesh%lc(jj,j))
        ! IF((row_obj%nelems==col_obj%nelems).AND.jk<ik.OR.jk<=0)CYCLE
        IF(jk==0)CYCLE
        dl_min = SIGN(1,row_obj%lfh(1,ii))*DOT_PRODUCT(evec_i(:,row_obj%lfh(2,ii)),evec_j(:,jj))*tmp
        DO irhs=1,nrhs
          !$omp atomic
          b(jk,irhs) = b(jk,irhs) + dl_min*a(ik,irhs)
        END DO
      END DO
      DO jj=col_obj%kfh(j),col_obj%kfh(j+1)-1
        jk=ABS(col_obj%lfh(1,jj))+col_obj%np_active
        ! IF((row_obj%nelems==col_obj%nelems).AND.jk<ik.OR.jk<=0)CYCLE
        dl_min = SIGN(1,row_obj%lfh(1,ii))*SIGN(1,col_obj%lfh(1,jj))* &
          DOT_PRODUCT(evec_i(:,row_obj%lfh(2,ii)),evec_j(:,col_obj%lfh(2,jj)))*tmp
        DO irhs=1,nrhs
          !$omp atomic
          b(jk,irhs) = b(jk,irhs) + dl_min*a(ik,irhs)
        END DO
      END DO
    END DO
  END DO
END DO
!$omp end do nowait
IF((col_obj%n_vcoils>0).OR.(row_obj%n_vcoils>0))CALL oft_warn("V-coil contributions were not computed.")
! !---Add passive coils to model
! !$omp do
! DO i=1,col_obj%n_vcoils
!   ii=col_obj%np_active+col_obj%nholes+i
!   DO j=1,row_obj%np_active+row_obj%nholes
!     b(ii,:) = b(ii,:) + row_obj%Ael2coil(j,i)*a(j,:)
!   END DO
! END DO
! !$omp end do nowait
! !$omp do
! DO i=1,col_obj%n_vcoils
!   ii=col_obj%np_active+col_obj%nholes+i
!   DO j=1,row_obj%n_vcoils
!     jj=row_obj%np_active+row_obj%nholes+j
!     b(ii,:) = b(ii,:) + row_obj%Acoil2coil(j,i)*a(jj,:)
!   END DO
! END DO
!$omp end parallel
b = b/(4.d0*pi)
CALL quad%delete()
CALL quad2%delete()
elapsed_time=mytimer%tock()
WRITE(*,'(3A)')oft_indent,'Time = ',time_to_string(elapsed_time)
CALL oft_decrease_indent
DEBUG_STACK_POP
END SUBROUTINE tw_compute_Lmat_MF_exp
!---------------------------------------------------------------------------------
!> Compute coupling from thin-wall model elements to flux loop sensors
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_mutuals_exp(tw_obj,nsensors,sensors,save_file)
TYPE(tw_type), INTENT(inout) :: tw_obj !< Thin-wall model object
INTEGER(4), INTENT(in) :: nsensors !< Number of flux loops
TYPE(floop_sensor), POINTER, INTENT(in) :: sensors(:) !< List of flux loops
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: save_file
LOGICAL :: exists
INTEGER(4) :: i,ii,j,jj,k,kk,l,ik,jk,ih,ihp,ihc,file_counts(5),ierr,io_unit,iquad,sens_hash
REAL(8) :: tmp,dl,pt_i(3),pt_j(3),evec_i(3,3),evec_j(3,3),pts_i(3,3),pt_i_last(3)
REAL(8) :: rvec_i(3),r1,rmag,rvec_j(3),z1,cvec(3),cpt(3),coil_thickness,elapsed_time
REAL(8) :: rgop(3,3),norm_i(3),area_i,f(3),dl_min,dl_max,pot_last,pot_tmp
REAL(8), allocatable :: atmp(:,:),Acoil2sen_tmp(:,:),Acoil2coil_tmp(:,:)
CLASS(oft_bmesh), POINTER :: bmesh
INTEGER(4) :: ncoils_tot
TYPE(tw_coil_set), POINTER, DIMENSION(:) :: coils_tot
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
!
!---Key for cached matrices: mesh, coils, sensors and cache version
sens_hash=tw_model_hash(tw_obj)
DO i=1,nsensors
  sens_hash=oft_simple_hash(C_LOC(sensors(i)%r),INT(8*3*sensors(i)%np,8),starting_value=sens_hash)
END DO
IF(PRESENT(save_file))THEN
  IF(TRIM(save_file)/='none')THEN
    INQUIRE(FILE=TRIM(save_file),EXIST=exists)
    IF(exists)THEN
      WRITE(*,'(2A)')oft_indent,'Reading sensor mutual matrices'
      CALL oft_increase_indent
      OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
      READ(io_unit, IOSTAT=ierr)file_counts
      IF((ierr/=0).OR.ANY(file_counts/=[tw_nelems_exp(tw_obj),tw_obj%n_vcoils,tw_obj%n_icoils,nsensors,sens_hash]))THEN
        exists=.FALSE.
        WRITE(*,'(2A)')oft_indent,'Ignoring stored matrix: Sizes do not match'
      END IF
      IF(exists)THEN
        ALLOCATE(tw_obj%Ael2sen(nsensors,tw_nelems_exp(tw_obj)))
        READ(io_unit, IOSTAT=ierr)tw_obj%Ael2sen
        IF(ierr/=0)THEN
          WRITE(*,'(2A)')oft_indent,'Error reading matrix from file'
          DEALLOCATE(tw_obj%Ael2sen)
          exists=.FALSE.
        END IF
      END IF
      ! IF(exists)THEN
      !   ALLOCATE(tw_obj%Acoil2sen(tw_obj%n_vcoils,nsensors))
      !   READ(io_unit, IOSTAT=ierr)tw_obj%Acoil2sen
      !   IF(ierr/=0)THEN
      !     WRITE(*,'(2A)')oft_indent,'Error reading matrix from file'
      !     DEALLOCATE(tw_obj%Ael2sen,tw_obj%Acoil2sen)
      !     exists=.FALSE.
      !   END IF
      ! END IF
      IF(exists)THEN
        ALLOCATE(tw_obj%Adr2sen(nsensors,tw_obj%n_icoils))
        READ(io_unit, IOSTAT=ierr)tw_obj%Adr2sen
        IF(ierr/=0)THEN
          WRITE(*,'(2A)')oft_indent,'Error reading matrix from file'
          DEALLOCATE(tw_obj%Ael2sen,tw_obj%Adr2sen)
          exists=.FALSE.
        END IF
      END IF
      CLOSE(io_unit)
      CALL oft_decrease_indent
    END IF
    IF(exists)THEN
      DEBUG_STACK_POP
      RETURN
    END IF
  END IF
END IF
!
bmesh=>tw_obj%mesh
!---Compute coupling between vertices and sensors
f=1.d0/3.d0
WRITE(*,'(2A)')oft_indent,'Building element->sensor inductance matrix'
CALL oft_increase_indent
IF(ASSOCIATED(tw_obj%Ael2sen))DEALLOCATE(tw_obj%Ael2sen)
ALLOCATE(tw_obj%Ael2sen(nsensors,tw_nelems_exp(tw_obj)))
tw_obj%Ael2sen=0.d0
IF(nsensors>0)THEN
  CALL mytimer%tick()
  !$omp parallel private(ii,j,jj,ik,pts_i,tmp,pt_i,pt_j,evec_i,rvec_j, &
  !$omp atmp,i,ih,ihp,ihc,rgop,area_i,norm_i,dl_min,dl_max,pot_last,iquad)
  ALLOCATE(atmp(3,nsensors))
  !$omp do
  DO i=1,bmesh%nc
    IF(tw_obj%sens_mask(bmesh%reg(i)))CYCLE
    ! CALL bmesh%jacobian(i,f,rgop,area_i)
    ! CALL bmesh%norm(i,f,norm_i)
    area_i=bmesh%ca(i)
    DO ii=1,3
      pts_i(:,ii)=bmesh%r(:,bmesh%lc(ii,i))
      ! evec_i(:,ii)=cross_product(rgop(:,ii),norm_i)
      evec_i(:,ii)=tw_obj%qbasis(:,ii,i)
    END DO
    !---Compute sensor couplings
    atmp=0.d0
    DO j=1,nsensors
      pot_last = 0.d0
      DO jj=1,sensors(j)%np
        pt_j=sensors(j)%r(:,jj)
        !---Compute minimum separation
        dl_min=1.d99
        dl_max=SQRT(area_i*2.d0)
        DO ii=1,3
          dl_min=MIN(dl_min,SQRT(SUM((pts_i(:,ii)-pt_j)**2)))
          dl_max=MAX(dl_max,SQRT(SUM((pts_i(:,ii)-pt_j)**2)))
        END DO
        IF(dl_min<1.d-8)THEN
          iquad = 18
        ELSE
          iquad = MAX(4,MIN(18,ABS(INT(LOG(target_err)/LOG(1.d0-dl_min/dl_max)))))
        END IF
        IF(iquad>10)THEN
          tmp=tw_compute_phipot(pts_i,pt_j)
        ELSE
          tmp=0.d0
          !$omp simd private(pt_i) reduction(+:tmp)
          DO ii=1,tw_quads(iquad)%np
            pt_i = tw_quads(iquad)%pts(1,ii)*pts_i(:,1) &
              + tw_quads(iquad)%pts(2,ii)*pts_i(:,2) &
              + tw_quads(iquad)%pts(3,ii)*pts_i(:,3)
            tmp = tmp + tw_quads(iquad)%wts(ii)/SQRT(SUM((pt_i-pt_j)**2))
          END DO
          tmp=tmp*area_i
        END IF
        IF(jj>1)THEN
          rvec_j=sensors(j)%r(:,jj)-sensors(j)%r(:,jj-1)
          DO ik=1,3
            atmp(ik,j)=atmp(ik,j)+DOT_PRODUCT(rvec_j,evec_i(:,ik))*(tmp+pot_last)/2.d0
          END DO
        END IF
        pot_last=tmp
      END DO
    END DO
    DO ii=1,3
      ik=tw_obj%pmap(bmesh%lc(ii,i))
      IF(ik==0)CYCLE
      DO j=1,nsensors
        !$omp atomic
        tw_obj%Ael2sen(j,ik)=tw_obj%Ael2sen(j,ik)+atmp(ii,j)
      END DO
    END DO
    DO ii=tw_obj%kfh(i),tw_obj%kfh(i+1)-1
      ik=ABS(tw_obj%lfh(1,ii))+tw_obj%np_active
      DO j=1,nsensors
        tmp=SIGN(1,tw_obj%lfh(1,ii))*atmp(tw_obj%lfh(2,ii),j)
        !$omp atomic
        tw_obj%Ael2sen(j,ik)=tw_obj%Ael2sen(j,ik)+tmp
      END DO
    END DO
  END DO
  DEALLOCATE(atmp)
  !$omp end parallel
  DO j=1,nsensors
    tw_obj%Ael2sen(j,:)=tw_obj%Ael2sen(j,:)*sensors(j)%scale_fac
  END DO
  elapsed_time=mytimer%tock()
  WRITE(*,'(3A)')oft_indent,'Time = ',time_to_string(elapsed_time)
ELSE
  WRITE(*,'(2A)')oft_indent,'No magnetic sensors, skipping...'
END IF
CALL oft_decrease_indent
!---Build temporary coil set from all filaments in model
ncoils_tot=tw_obj%n_vcoils+tw_obj%n_icoils
ALLOCATE(coils_tot(ncoils_tot))
DO i=1,tw_obj%n_vcoils
  CALL tw_copy_coil(tw_obj%vcoils(i),coils_tot(i))
END DO
DO i=1,tw_obj%n_icoils
  CALL tw_copy_coil(tw_obj%icoils(i),coils_tot(i+tw_obj%n_vcoils))
END DO
!---Compute coupling between coils and sensors
ALLOCATE(Acoil2sen_tmp(nsensors,ncoils_tot))
Acoil2sen_tmp=0.d0
WRITE(*,'(2A)')oft_indent,'Building coil->sensor inductance matrix'
CALL oft_increase_indent
IF(nsensors>0.AND.ncoils_tot>0)THEN
  CALL mytimer%tick()
  !$omp parallel private(ii,j,k,kk,tmp,pt_i,rvec_i,rmag,atmp,cvec,cpt,pot_tmp,pt_i_last,pot_last)
  ALLOCATE(atmp(ncoils_tot,1))
  !$omp do
  DO i=1,nsensors
    atmp=0.d0
    DO ii=1,sensors(i)%np
      pt_i=sensors(i)%r(:,ii)
      !
      IF(ii>1)THEN
        rvec_i=sensors(i)%r(:,ii)-sensors(i)%r(:,ii-1)
        DO j=1,ncoils_tot
          IF(coils_tot(j)%sens_mask)CYCLE
          DO k=1,coils_tot(j)%ncoils
            IF(coils_tot(j)%coils(k)%is_ring)THEN
              tmp=DOT_PRODUCT(rvec_i,tw_ring_avec(coils_tot(j)%coils(k)%ring_rz,pt_i) &
                + tw_ring_avec(coils_tot(j)%coils(k)%ring_rz,pt_i_last))/2.d0
              atmp(j,1)=atmp(j,1)+coils_tot(j)%scales(k)*tmp
              CYCLE
            END IF
            tmp=0.d0
            pot_last = 0.d0
            DO kk=1,coils_tot(j)%coils(k)%npts
              cpt = coils_tot(j)%coils(k)%pts(:,kk)
              pot_tmp = (1.d0/SQRT(SUM((pt_i-cpt)**2)) + 1.d0/SQRT(SUM((pt_i_last-cpt)**2)))/2.d0
              IF(kk>1)THEN
                cvec = coils_tot(j)%coils(k)%pts(:,kk)-coils_tot(j)%coils(k)%pts(:,kk-1)
                tmp = tmp + DOT_PRODUCT(rvec_i,cvec)*(pot_tmp+pot_last)/2.d0
              END IF
              pot_last = pot_tmp
            END DO
            atmp(j,1)=atmp(j,1)+coils_tot(j)%scales(k)*tmp
          END DO
        END DO
      END IF
      pt_i_last=pt_i
    END DO
    DO j=1,ncoils_tot
      Acoil2sen_tmp(i,j) = atmp(j,1)
    END DO
  END DO
  DEALLOCATE(atmp)
  !$omp end parallel
  DO j=1,nsensors
    Acoil2sen_tmp(j,:)=Acoil2sen_tmp(j,:)*sensors(j)%scale_fac
  END DO
  elapsed_time=mytimer%tock()
  WRITE(*,'(3A)')oft_indent,'Time = ',time_to_string(elapsed_time)
ELSE
  WRITE(*,'(2A)')oft_indent,'No magnetic sensors or coils, skipping...'
END IF
DEALLOCATE(coils_tot)
CALL oft_decrease_indent
!---Unpack passive and driver coils
! IF(ASSOCIATED(tw_obj%Acoil2sen))DEALLOCATE(tw_obj%Acoil2sen)
! ALLOCATE(tw_obj%Acoil2sen(nsensors,tw_obj%n_vcoils))
! DO i=1,tw_obj%n_vcoils
!   tw_obj%Acoil2sen(:,i)=Acoil2sen_tmp(:,i)
! END DO
IF(ASSOCIATED(tw_obj%Adr2sen))DEALLOCATE(tw_obj%Adr2sen)
ALLOCATE(tw_obj%Adr2sen(nsensors,tw_obj%n_icoils))
DO i=1,tw_obj%n_icoils
  tw_obj%Adr2sen(:,i)=Acoil2sen_tmp(:,i+tw_obj%n_vcoils)
END DO
tw_obj%Adr2sen=tw_obj%Adr2sen*mu0/(4.d0*pi)
!---Copy coupling between passive coils and sensors
IF(nsensors>0)THEN
  !$omp parallel private(j,ii,jj,ik,jk)
  !$omp do
  DO i=1,tw_obj%n_vcoils
    ! Sensors
    DO jj=1,nsensors
      tw_obj%Ael2sen(jj,tw_obj%np_active+tw_obj%nholes+i)=Acoil2sen_tmp(jj,i)! tw_obj%Acoil2sen(jj,i)
    END DO
  END DO
  !$omp end parallel
  tw_obj%Ael2sen = tw_obj%Ael2sen/(4.d0*pi)
END IF
DEALLOCATE(Acoil2sen_tmp)
!
!---
IF(PRESENT(save_file))THEN
  IF(TRIM(save_file)/='none')THEN
    OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
    WRITE(io_unit)tw_nelems_exp(tw_obj),tw_obj%n_vcoils,tw_obj%n_icoils,nsensors,sens_hash
    WRITE(io_unit)tw_obj%Ael2sen
    ! WRITE(io_unit)tw_obj%Acoil2sen
    WRITE(io_unit)tw_obj%Adr2sen
    CLOSE(io_unit)
  END IF
END IF
DEBUG_STACK_POP
END SUBROUTINE tw_compute_mutuals_exp
!---------------------------------------------------------------------------------
!> Compute resistivity matrix for thin-wall model
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_Rmat_exp(tw_obj,keep_closures)
TYPE(tw_type), INTENT(inout) :: tw_obj !< Thin-wall model object
LOGICAL, INTENT(in) :: keep_closures !< Keep diagonal entries (1) for closure elements
INTEGER(4) :: i,ii,j,jj,jj2,k,kk,kr,ed,ninteract,pt,ih,ihp,ihc,max_felems
INTEGER(4), ALLOCATABLE :: j_add(:),face_interact(:,:,:),krtmp(:)
REAL(8) :: eta1,eta2,eta_eff,rsign,rmag,rcurr(3),dl,f(3),gop(3,3),area,norm(3)
REAL(8), ALLOCATABLE :: curr_signs(:),eta_add(:,:),evec(:,:)
LOGICAL, ALLOCATABLE :: hole_mat(:,:),hole_mark(:)
CLASS(oft_bmesh), POINTER :: bmesh
TYPE(oft_native_matrix), POINTER :: Rmat
DEBUG_STACK_PUSH
bmesh=>tw_obj%mesh
WRITE(*,'(2A)')oft_indent,'Building resistivity matrix'
CALL oft_increase_indent
IF(ALL(tw_obj%Eta_surf<0.d0))THEN
  CALL oft_warn('Resistivity not set, using "eta=mu0" for all regions')
  tw_obj%Eta_surf=1.d0
END IF
ALLOCATE(tw_obj%Rmat)
Rmat=>tw_obj%Rmat
Rmat%nr=tw_nelems_exp(tw_obj); Rmat%nrg=tw_nelems_exp(tw_obj)
Rmat%nc=tw_nelems_exp(tw_obj); Rmat%ncg=tw_nelems_exp(tw_obj)
!---Build face element lists
max_felems=bmesh%cell_np+2*tw_obj%nholes
ALLOCATE(face_interact(2,max_felems,bmesh%nc))
face_interact=0
DO i=1,bmesh%nc
  k=0
  DO j=1,bmesh%cell_np
    IF(tw_obj%pmap(bmesh%lc(j,i))>0)THEN
      k=k+1
      face_interact(:,k,i)=[tw_obj%pmap(bmesh%lc(j,i)),j]
    END IF
  END DO
END DO
! Add holes
DO ih=1,tw_obj%nholes
DO ihp=1,tw_obj%hmesh(ih)%n
DO ihc=tw_obj%hmesh(ih)%kpc(ihp),tw_obj%hmesh(ih)%kpc(ihp+1)-1
i=ABS(tw_obj%hmesh(ih)%lpc(ihc))
  DO k=1,bmesh%cell_np
    IF(bmesh%lc(k,i)==tw_obj%hmesh(ih)%lp(ihp))EXIT
  END DO
  DO j=1,max_felems
    IF(face_interact(1,j,i)==0)THEN
      face_interact(:,j,i)=[tw_obj%np_active+ih,SIGN(k,tw_obj%hmesh(ih)%lpc(ihc))]
      EXIT
    END IF
  END DO
  IF(j>max_felems)CALL oft_abort("Exceeded max interaction size","tw_compute_Rmat",__FILE__)
END DO
END DO
END DO
! Find hole to hole interactions
ALLOCATE(hole_mat(tw_obj%nholes,tw_obj%nholes))
hole_mat=.FALSE.
DO i=1,bmesh%nc
  DO j=1,max_felems
    IF(face_interact(1,j,i)==0)EXIT
    IF(face_interact(1,j,i)>tw_obj%np_active)THEN
      DO k=1,max_felems
        IF(face_interact(1,k,i)==0)EXIT
        IF(face_interact(1,k,i)>tw_obj%np_active)THEN
          hole_mat(face_interact(1,j,i)-tw_obj%np_active, &
            face_interact(1,k,i)-tw_obj%np_active)=.TRUE.
        END IF
      END DO
    END IF
  END DO
END DO
!
ALLOCATE(Rmat%kr(Rmat%nr+1),krtmp(Rmat%nr+1))
Rmat%kr=0
! !$omp parallel private(ii,j,k,ed,hole_mark)
ALLOCATE(hole_mark(tw_obj%nholes))
! !$omp do
DO i=1,bmesh%np
  ii=tw_obj%pmap(i)
  IF(ii==0)CYCLE
  DO j=bmesh%kpp(i),bmesh%kpp(i+1)-1
    IF(tw_obj%pmap(bmesh%lpp(j))>0)Rmat%kr(ii) = Rmat%kr(ii) + 1
  END DO
  !---Add holes
  hole_mark=.FALSE.
  DO j=bmesh%kpc(i),bmesh%kpc(i+1)-1
    DO k=1,max_felems
      ih=face_interact(1,k,bmesh%lpc(j))
      IF(ih==0)EXIT
      IF(ih>tw_obj%np_active)hole_mark(ih-tw_obj%np_active)=.TRUE.
    END DO
  END DO
  DO j=1,tw_obj%nholes
    IF(hole_mark(j))THEN
      Rmat%kr(ii) = Rmat%kr(ii) + 1
      Rmat%kr(j+tw_obj%np_active) = Rmat%kr(j+tw_obj%np_active) + 1 ! Reverse interaction
    END IF
  END DO
END DO
DEALLOCATE(hole_mark)
!---Add hole to hole interactions
! !$omp do
DO i=1,tw_obj%nholes
  DO j=1,tw_obj%nholes
    IF(hole_mat(i,j))Rmat%kr(tw_obj%np_active+i)=Rmat%kr(tw_obj%np_active+i)+1
  END DO
END DO
!---Add passive coils
! !$omp do
DO i=1,tw_obj%n_vcoils
  Rmat%kr(tw_obj%np_active+tw_obj%nholes+i)=1
END DO
! !$omp end parallel
Rmat%nnz=SUM(Rmat%kr)
Rmat%kr(Rmat%nr+1)=Rmat%nnz+1
do i=Rmat%nr,1,-1 ! cumulative point to point count
  Rmat%kr(i) = Rmat%kr(i+1) - Rmat%kr(i)
end do
if(Rmat%kr(1)/=1)call oft_abort('Bad element to element count','tw_compute_Rmat',__FILE__)
krtmp=Rmat%kr-1
ALLOCATE(Rmat%lc(Rmat%nnz))
! !$omp parallel private(ii,j,k,kr,ed,hole_mark)
ALLOCATE(hole_mark(tw_obj%nholes))
! !$omp do
DO i=1,bmesh%np
  ii=tw_obj%pmap(i)
  IF(ii==0)CYCLE
  DO j=bmesh%kpp(i),bmesh%kpp(i+1)-1
    IF(tw_obj%pmap(bmesh%lpp(j))>0)THEN
      krtmp(ii)=krtmp(ii)+1
      Rmat%lc(krtmp(ii)) = tw_obj%pmap(bmesh%lpp(j))
    END IF
  END DO
  !---Add holes
  hole_mark=.FALSE.
  DO j=bmesh%kpc(i),bmesh%kpc(i+1)-1
    DO k=1,max_felems
      ih=face_interact(1,k,bmesh%lpc(j))
      IF(ih==0)EXIT
      IF(ih>tw_obj%np_active)hole_mark(ih-tw_obj%np_active)=.TRUE.
    END DO
  END DO
  DO j=1,tw_obj%nholes
    IF(hole_mark(j))THEN
      krtmp(ii)=krtmp(ii)+1
      Rmat%lc(krtmp(ii)) = j+tw_obj%np_active
      krtmp(j+tw_obj%np_active)=krtmp(j+tw_obj%np_active)+1
      Rmat%lc(krtmp(j+tw_obj%np_active)) = ii ! Reverse interaction
    END IF
  END DO
  CALL sort_array(Rmat%lc(Rmat%kr(ii):krtmp(ii)),krtmp(ii)+1-Rmat%kr(ii))
END DO
DEALLOCATE(hole_mark)
!---Add hole to hole interactions
! !$omp do
DO i=1,tw_obj%nholes
  DO j=1,tw_obj%nholes
    IF(hole_mat(i,j))THEN
      krtmp(tw_obj%np_active+i)=krtmp(tw_obj%np_active+i)+1
      Rmat%lc(krtmp(tw_obj%np_active+i))=tw_obj%np_active+j
    END IF
  END DO
  CALL sort_array(Rmat%lc(Rmat%kr(tw_obj%np_active+i):krtmp(tw_obj%np_active+i)), &
    krtmp(tw_obj%np_active+i)+1-Rmat%kr(tw_obj%np_active+i))
END DO
DEALLOCATE(hole_mat)
!---Add passive coils
! !$omp do
DO i=1,tw_obj%n_vcoils
  Rmat%lc(Rmat%kr(tw_obj%np_active+tw_obj%nholes+i))=tw_obj%np_active+tw_obj%nholes+i
END DO
! !$omp end parallel
CALL csr_remove_redundant(Rmat%nr,Rmat%kr,Rmat%nnz,Rmat%lc)
!---Compute resistivity matrix
ALLOCATE(Rmat%M(Rmat%nnz))
Rmat%M=0.d0
ALLOCATE(j_add(max_felems),evec(3,max_felems))
ALLOCATE(eta_add(max_felems,max_felems))
f=1.d0/3.d0
DO i=1,bmesh%nc
  j_add=0
  ninteract=0
  DO j=1,max_felems
    IF(face_interact(1,j,i)==0)EXIT
    ninteract=ninteract+1
    j_add(ninteract)=face_interact(1,j,i)
  END DO
  eta_eff=tw_obj%Eta_surf(bmesh%reg(i))
  !---Compute mass coupling
  CALL bmesh%jacobian(i,f,gop,area)
  CALL bmesh%norm(i,f,norm)
  DO j=1,ninteract
    evec(:,j)=cross_product(gop(:,ABS(face_interact(2,j,i))),norm)*SIGN(1,face_interact(2,j,i))
  END DO
  DO j=1,ninteract
    DO k=1,ninteract
      eta_add(j,k)=eta_eff*DOT_PRODUCT(evec(:,j),evec(:,k))*area
    END DO
  END DO
  !---Add values to matrix
  CALL Rmat%add_values(j_add(1:ninteract), j_add(1:ninteract), &
                       eta_add(1:ninteract,1:ninteract), ninteract, ninteract)
END DO
DEALLOCATE(j_add,eta_add,evec)
!---Add passive coils
ALLOCATE(j_add(1),eta_add(1,1))
DO i=1,tw_obj%n_vcoils
  ! Compute resistivity from length
  tw_obj%vcoils(i)%Rself=0.d0
  DO j=1,tw_obj%vcoils(i)%ncoils
    dl=0.d0
    IF(tw_obj%vcoils(i)%coils(j)%is_ring)THEN
      dl = 2.d0*pi*tw_obj%vcoils(i)%coils(j)%ring_rz(1) ! Exact length of a circular coil
    ELSE
      DO k=2,tw_obj%vcoils(i)%coils(j)%npts
        dl = dl + SQRT(SUM((tw_obj%vcoils(i)%coils(j)%pts(:,k)-tw_obj%vcoils(i)%coils(j)%pts(:,k-1))**2))
      END DO
    END IF
    tw_obj%vcoils(i)%Rself = tw_obj%vcoils(i)%Rself + tw_obj%vcoils(i)%res_per_len(j)*dl
  END DO
  WRITE(*,"(2A,1X,I4,A,ES12.4)")oft_indent,"Vcoil",i,": R [Ohm] = ",tw_obj%vcoils(i)%Rself !*pi*4.d-7
  tw_obj%vcoils(i)%Rself = tw_obj%vcoils(i)%Rself/mu0 ! Convert to magnetic units
  !
  eta_add=tw_obj%vcoils(i)%Rself
  j_add(1)=tw_obj%np_active+tw_obj%nholes+i
  CALL Rmat%add_values(j_add(1:1), j_add(1:1), eta_add(1:1,1:1), 1, 1)
END DO
DEALLOCATE(j_add,eta_add)
! CALL Rmat%assemble()
!---Handle closures
! IF(tw_obj%nclosures>0)THEN
!   CALL Rmat%zero_rows(tw_obj%nclosures,tw_obj%closures)
!   CALL Rmat%zero_cols(tw_obj%nclosures,tw_obj%closures)
!   IF(keep_closures)THEN
!     ALLOCATE(j_add(1),eta_add(1,1))
!     eta_add=1.d0
!     DO i=1,tw_obj%nclosures
!       j_add(1)=tw_obj%closures(i)
!       CALL Rmat%add_values(j_add(1:1), j_add(1:1), eta_add(1:1,1:1), 1, 1)
!     END DO
!     DEALLOCATE(j_add,eta_add)
!   END IF
! END IF
CALL oft_decrease_indent
DEBUG_STACK_POP
END SUBROUTINE tw_compute_Rmat_exp
!---------------------------------------------------------------------------------
!> Compute \f$ \int \int 1/|r-r'| dA dA' \f$ for a pair of triangles
!!
!! The quadrature order is chosen from the pair's separation; close pairs integrate the inner triangle analytically
!! (@ref tw_compute_phipot). Uses the shared quadrature rules (@ref tw_quads_setup, done by @ref tw_setup).
!---------------------------------------------------------------------------------
FUNCTION tw_cell_pair_kernel(pts_i,area_i,pts_j,area_j) RESULT(tmp)
REAL(8), INTENT(in) :: pts_i(3,3) !< Vertices of first triangle
REAL(8), INTENT(in) :: area_i !< Area of first triangle
REAL(8), INTENT(in) :: pts_j(3,3) !< Vertices of second triangle
REAL(8), INTENT(in) :: area_j !< Area of second triangle
REAL(8) :: tmp
INTEGER(4) :: ii,jj,iquad,np
REAL(8) :: dl_min,dl_max,d2,s,pt_i(3),pt_j(3),xj(tw_quad_far_np),yj(tw_quad_far_np),zj(tw_quad_far_np)
!---Compute minimum separation (squared distances)
dl_min=1.d99
dl_max=MAX(area_i,area_j)*2.d0
DO ii=1,3
  DO jj=1,3
    d2=SUM((pts_i(:,ii)-pts_j(:,jj))**2)
    dl_min=MIN(dl_min,d2)
    dl_max=MAX(dl_max,d2)
  END DO
END DO
dl_min=SQRT(dl_min)
dl_max=SQRT(dl_max)
!---Chose quadrature order based on distance
IF(dl_min<1.d-8)THEN
  iquad = 18
ELSE
  iquad = MAX(4,MIN(18,ABS(INT(LOG(tw_L_quad_err)/LOG(1.d0-dl_min/dl_max)))))
END IF
tmp=0.d0
np=tw_quads(iquad)%np
IF(iquad>tw_quad_far_max)THEN
  DO jj=1,np
    pt_j = tw_quads(iquad)%pts(1,jj)*pts_j(:,1) &
      + tw_quads(iquad)%pts(2,jj)*pts_j(:,2) &
      + tw_quads(iquad)%pts(3,jj)*pts_j(:,3)
    tmp = tmp + tw_quads(iquad)%wts(jj)*tw_compute_phipot(pts_i,pt_j)
  END DO
  tmp=tmp*area_j
ELSE
  !---Map the second triangle's points once, then a vectorizable 1/r loop
  DO jj=1,np
    pt_j = tw_quads(iquad)%pts(1,jj)*pts_j(:,1) &
      + tw_quads(iquad)%pts(2,jj)*pts_j(:,2) &
      + tw_quads(iquad)%pts(3,jj)*pts_j(:,3)
    xj(jj)=pt_j(1); yj(jj)=pt_j(2); zj(jj)=pt_j(3)
  END DO
  DO ii=1,np
    pt_i = tw_quads(iquad)%pts(1,ii)*pts_i(:,1) &
      + tw_quads(iquad)%pts(2,ii)*pts_i(:,2) &
      + tw_quads(iquad)%pts(3,ii)*pts_i(:,3)
    s=0.d0
    !$omp simd reduction(+:s)
    DO jj=1,np
      s = s + tw_quads(iquad)%wts(jj)/SQRT((pt_i(1)-xj(jj))**2+(pt_i(2)-yj(jj))**2+(pt_i(3)-zj(jj))**2)
    END DO
    tmp = tmp + tw_quads(iquad)%wts(ii)*s
  END DO
  tmp=tmp*area_i*area_j
END IF
END FUNCTION tw_cell_pair_kernel
!---------------------------------------------------------------------------------
!> Compute the B-field at a point from the three basis functions of a triangle (without the 1/(4 pi) factor)
!!
!! Far points use quadrature of the Biot-Savart kernel; near points use the analytic gradient of the cell potential
!! (@ref tw_compute_phigrad). Points on the cell are evaluated `10*B_dx` inside and outside the face and averaged. Uses the shared quadrature rules (@ref tw_quads_setup, done by @ref tw_setup).
!---------------------------------------------------------------------------------
FUNCTION tw_cell_point_Bkernel(pts_i,evec_i,area_i,norm_i,pt,va,B_dx) RESULT(atmp)
REAL(8), INTENT(in) :: pts_i(3,3) !< Vertices of triangle
REAL(8), INTENT(in) :: evec_i(3,3) !< Basis function current directions
REAL(8), INTENT(in) :: area_i !< Area of triangle
REAL(8), INTENT(in) :: norm_i(3) !< Unit normal of triangle
REAL(8), INTENT(in) :: pt(3) !< Evaluation point
REAL(8), INTENT(in) :: va !< Vertex area of evaluation point
REAL(8), INTENT(in) :: B_dx !< Length scale for points on the cell (evaluated 10*B_dx off the face)
REAL(8) :: atmp(3,3) !< B-field from each basis function [3,3]
INTEGER(4) :: ii,ik,iquad
REAL(8) :: dl_min,dl_max,d2,pt_j(3),pt_i(3),diffvec(3),gvec(3)
LOGICAL :: is_neighbor
!---Compute minimum separation (squared distances)
dl_min=1.d99
dl_max=MAX(area_i,va/(pi**2))
DO ii=1,3
  d2=SUM((pts_i(:,ii)-pt)**2)
  dl_min=MIN(dl_min,d2)
  dl_max=MAX(dl_max,d2)
END DO
dl_min=SQRT(dl_min)
dl_max=SQRT(dl_max)
!---Chose quadrature order based on distance
IF(dl_min<1.d-8)THEN
  iquad = 18
  is_neighbor=.TRUE.
ELSE
  iquad = MAX(4,MIN(18,ABS(INT(LOG(target_err)/LOG(1.d0-dl_min/dl_max)))))
  is_neighbor=.FALSE.
END IF
!
pt_j=pt
IF(iquad>10)THEN
  IF(is_neighbor)THEN
    !---Point on the cell: average just inside and just outside the face (the tangential field jumps across the sheet)
    diffvec=(tw_compute_phigrad(pts_i,pt-norm_i*10.d0*B_dx)+tw_compute_phigrad(pts_i,pt+norm_i*10.d0*B_dx))/2.d0
  ELSE
    diffvec=tw_compute_phigrad(pts_i,pt)
  END IF
  DO ik=1,3
    atmp(1,ik) = diffvec(2)*evec_i(3,ik) - diffvec(3)*evec_i(2,ik)
    atmp(2,ik) = diffvec(3)*evec_i(1,ik) - diffvec(1)*evec_i(3,ik)
    atmp(3,ik) = diffvec(1)*evec_i(2,ik) - diffvec(2)*evec_i(1,ik)
  END DO
ELSE
  !---sum_q w (e x d)/|d|^3 = e x sum_q w d/|d|^3: one pass over the points for all three basis functions
  gvec=0.d0
  DO ii=1,tw_quads(iquad)%np
    pt_i = pt_j - (tw_quads(iquad)%pts(1,ii)*pts_i(:,1) &
      + tw_quads(iquad)%pts(2,ii)*pts_i(:,2) &
      + tw_quads(iquad)%pts(3,ii)*pts_i(:,3))
    d2=SUM(pt_i**2)
    gvec = gvec + pt_i*(tw_quads(iquad)%wts(ii)/(d2*SQRT(d2)))
  END DO
  gvec=gvec*area_i
  DO ik=1,3
    atmp(1,ik) = evec_i(2,ik)*gvec(3) - evec_i(3,ik)*gvec(2)
    atmp(2,ik) = evec_i(3,ik)*gvec(1) - evec_i(1,ik)*gvec(3)
    atmp(3,ik) = evec_i(1,ik)*gvec(2) - evec_i(2,ik)*gvec(1)
  END DO
END IF
END FUNCTION tw_cell_point_Bkernel
!---------------------------------------------------------------------------------
!> Compute \f$ \int 1/(r-r') dA' \f$ for a triangle
!---------------------------------------------------------------------------------
FUNCTION tw_compute_phipot(pt_cell,pt) RESULT(phi)
REAL(8), INTENT(in) :: pt_cell(3,3) !< Vertices defining triangle
REAL(8), INTENT(in) :: pt(3) !< Observation point (r)
REAL(8) :: phi
REAL(8) :: nhat(3),a(3),b(3),r(3,3),rmag(3),e(3),elen,num,den,omega
INTEGER(4) :: i,ip
DEBUG_STACK_PUSH
!---Get unit normal
a=pt_cell(:,2)-pt_cell(:,1)
b=pt_cell(:,3)-pt_cell(:,2)
nhat=[a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
nhat=nhat/SQRT(SUM(nhat**2))
!---Compute local vectors
DO i=1,3
  r(:,i)=pt_cell(:,i)-pt
  rmag(i)=SQRT(SUM(r(:,i)**2))
END DO
!---Compute solid angle
num=r(1,1)*(r(2,2)*r(3,3)-r(3,2)*r(2,3)) + r(2,1)*(r(3,2)*r(1,3)-r(1,2)*r(3,3)) &
  + r(3,1)*(r(1,2)*r(2,3)-r(2,2)*r(1,3))
den=rmag(1)*rmag(2)*rmag(3) + SUM(r(:,1)*r(:,2))*rmag(3) &
  + SUM(r(:,1)*r(:,3))*rmag(2) + SUM(r(:,2)*r(:,3))*rmag(1)
omega=2.d0*ATAN2(num,den)
!---Edge terms; nhat.(c_i x c_i+1) = nhat.(r_i x r_i+1) since the in-plane parts c differ from r by the same normal offset
phi=-DOT_PRODUCT(nhat,r(:,1))*omega
DO i=1,3
  ip=MOD(i,3)+1
  e=r(:,ip)-r(:,i)
  elen=SQRT(SUM(e**2))
  num=rmag(ip)*elen + SUM(r(:,ip)*e)
  den=rmag(i)*elen + SUM(r(:,i)*e)
  IF(ABS(den)<1.d-14.OR.elen<1.d-14)CYCLE
  phi=phi+(nhat(1)*(r(2,i)*r(3,ip)-r(3,i)*r(2,ip)) + nhat(2)*(r(3,i)*r(1,ip)-r(1,i)*r(3,ip)) &
    + nhat(3)*(r(1,i)*r(2,ip)-r(2,i)*r(1,ip)))*LOG(num/den)/elen
END DO
DEBUG_STACK_POP
END FUNCTION tw_compute_phipot
!---------------------------------------------------------------------------------
!> Compute \f$ \nabla \int 1/|r-r'| dA' \f$ for a triangle (gradient with respect to the observation point)
!!
!! \f$ \nabla \phi = -\sum_i \hat{m}_i \int_{e_i} dl'/|r-r'| + \hat{n} \Omega \f$, with \f$ \hat{m}_i \f$ the outward
!! in-plane normal of edge \f$ e_i \f$ and \f$ \Omega \f$ the signed solid angle (as in @ref tw_compute_phipot)
!---------------------------------------------------------------------------------
FUNCTION tw_compute_phigrad(pt_cell,pt) RESULT(grad)
REAL(8), INTENT(in) :: pt_cell(3,3) !< Vertices defining triangle
REAL(8), INTENT(in) :: pt(3) !< Observation point (r)
REAL(8) :: grad(3)
REAL(8) :: nhat(3),a(3),b(3),r(3,3),rmag(3),e(3),elen,num,den
INTEGER(4) :: i,ip
DEBUG_STACK_PUSH
!---Get unit normal
a=pt_cell(:,2)-pt_cell(:,1)
b=pt_cell(:,3)-pt_cell(:,2)
nhat=[a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
nhat=nhat/SQRT(SUM(nhat**2))
!---Compute local vectors
DO i=1,3
  r(:,i)=pt_cell(:,i)-pt
  rmag(i)=SQRT(SUM(r(:,i)**2))
END DO
!---Normal part: signed solid angle
num=r(1,1)*(r(2,2)*r(3,3)-r(3,2)*r(2,3)) + r(2,1)*(r(3,2)*r(1,3)-r(1,2)*r(3,3)) &
  + r(3,1)*(r(1,2)*r(2,3)-r(2,2)*r(1,3))
den=rmag(1)*rmag(2)*rmag(3) + SUM(r(:,1)*r(:,2))*rmag(3) &
  + SUM(r(:,1)*r(:,3))*rmag(2) + SUM(r(:,2)*r(:,3))*rmag(1)
grad=nhat*2.d0*ATAN2(num,den)
!---In-plane part: outward edge normals times the edge integrals of 1/|r-r'|
DO i=1,3
  ip=MOD(i,3)+1
  e=r(:,ip)-r(:,i)
  elen=SQRT(SUM(e**2))
  num=rmag(ip)*elen + SUM(r(:,ip)*e)
  den=rmag(i)*elen + SUM(r(:,i)*e)
  IF(ABS(den)<1.d-14.OR.elen<1.d-14)CYCLE
  grad=grad-[e(2)*nhat(3)-e(3)*nhat(2), e(3)*nhat(1)-e(1)*nhat(3), e(1)*nhat(2)-e(2)*nhat(1)]*LOG(num/den)/elen
END DO
DEBUG_STACK_POP
END FUNCTION tw_compute_phigrad
!---------------------------------------------------------------------------------
!> Compute the magnetic field reconstruction operators at the mesh vertices
!!
!! For models with junctions the element operators are scattered through the junction map as they are
!! computed, giving the reduced operator (`Bel(:,e,:)` P) without storing the operator over expanded elements.
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_Bops(self,B_dx,save_file)
TYPE(tw_type), INTENT(inout) :: self !< Thin-wall model object
REAL(r8), INTENT(in) :: B_dx !< Distance for finite-difference field evaluation
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: save_file !< Cache file
REAL(r8) :: evec_i(3,3),pts_i(3,3),pt_j(3),diffvec(3),ecc(3)
REAL(r8) :: cvec(3),cpt(3),tmp,area_i,norm_j(3),f(3)
REAL(r8), ALLOCATABLE :: atmp(:,:,:)
INTEGER(4) :: i,ii,j,jj,ik,k,kk,nsurf
LOGICAL :: load_success
CLASS(oft_bmesh), POINTER :: bmesh
IF(PRESENT(save_file))THEN
  load_success=load_from_file()
  IF(load_success)RETURN
END IF
!
bmesh=>self%mesh
ALLOCATE(self%Bel(bmesh%np,self%nelems,3))
self%Bel=0.d0
f=1.d0/3.d0
WRITE(*,'(2A)')oft_indent,'Building element->element magnetic reconstruction operator'
!$omp parallel private(ii,j,jj,ik,pts_i,tmp,pt_j,evec_i, &
!$omp atmp,i,area_i,norm_j,diffvec)
ALLOCATE(atmp(3,3,bmesh%np))
!$omp do schedule(dynamic,100)
DO i=1,bmesh%nc
  ! CALL bmesh%jacobian(i,f,rgop,area_i)
  ! CALL bmesh%norm(i,f,norm_i)
  area_i=bmesh%ca(i)
  DO ii=1,3
    pts_i(:,ii)=bmesh%r(:,bmesh%lc(ii,i))
    ! evec_i(:,ii)=cross_product(rgop(:,ii),norm_i)
    evec_i(:,ii)=self%qbasis(:,ii,i)
  END DO
  CALL bmesh%norm(i,f,norm_j)
  !---Compute sensor couplings
  atmp=0.d0
  DO j=1,bmesh%np
    !---Field from the cell's three basis functions
    atmp(:,:,j)=tw_cell_point_Bkernel(pts_i,evec_i,area_i,norm_j,bmesh%r(:,j),bmesh%va(j),B_dx)
  END DO
  DO ii=1,3
    ik=self%pmap(bmesh%lc(ii,i))
    IF(ik==0)CYCLE
    CALL add_elem(atmp,ii,ik,1.d0)
  END DO
  DO ii=self%kfh(i),self%kfh(i+1)-1
    ik=ABS(self%lfh(1,ii))+self%np_active
    CALL add_elem(atmp,self%lfh(2,ii),ik,REAL(SIGN(1,self%lfh(1,ii)),r8))
  END DO
END DO
DEALLOCATE(atmp)
!$omp end parallel
!
WRITE(*,'(2A)')oft_indent,'Building vcoil->element magnetic reconstruction operator'
nsurf=self%nelems-self%n_vcoils ! V-coil rows of the junction map are the identity
!$omp parallel do private(ii,j,k,kk,pt_j,ecc,diffvec,cvec,cpt)
DO i=1,bmesh%np
  pt_j=bmesh%r(:,i)
  !---Compute driver contributions
  DO j=1,self%n_vcoils
    ecc = 0.d0
    IF(self%vcoils(j)%sens_mask)CYCLE
    DO k=1,self%vcoils(j)%ncoils
      IF(self%vcoils(j)%coils(k)%is_ring)THEN
        diffvec=tw_ring_bvec(self%vcoils(j)%coils(k)%ring_rz,pt_j)
      ELSE
        diffvec=0.d0
        DO kk=2,self%vcoils(j)%coils(k)%npts
          cvec = self%vcoils(j)%coils(k)%pts(:,kk)-self%vcoils(j)%coils(k)%pts(:,kk-1)
          cpt = (self%vcoils(j)%coils(k)%pts(:,kk)+self%vcoils(j)%coils(k)%pts(:,kk-1))/2.d0
          diffvec = diffvec + cross_product(cvec,pt_j-cpt)/SUM((pt_j-cpt)**2)**1.5d0
        END DO
      END IF
      ecc=ecc+self%vcoils(j)%scales(k)*diffvec
    END DO
    DO jj=1,3
      !$omp atomic
      self%Bel(i,nsurf+j,jj) = self%Bel(i,nsurf+j,jj) + ecc(jj)
    END DO
  END DO
END DO
self%Bel=self%Bel/(4.d0*pi)
!
WRITE(*,'(2A)')oft_indent,'Building icoil->element magnetic reconstruction operator'
ALLOCATE(self%Bdr(bmesh%np,self%n_icoils,3))
self%Bdr=0.d0
!$omp parallel do private(ii,j,k,kk,pt_j,ecc,diffvec,cvec,cpt)
DO i=1,bmesh%np
  pt_j=bmesh%r(:,i)
  !---Compute driver contributions
  DO j=1,self%n_icoils
    ecc = 0.d0
    IF(self%icoils(j)%sens_mask)CYCLE
    DO k=1,self%icoils(j)%ncoils
      IF(self%icoils(j)%coils(k)%is_ring)THEN
        diffvec=tw_ring_bvec(self%icoils(j)%coils(k)%ring_rz,pt_j)
      ELSE
        diffvec=0.d0
        DO kk=2,self%icoils(j)%coils(k)%npts
          cvec = self%icoils(j)%coils(k)%pts(:,kk)-self%icoils(j)%coils(k)%pts(:,kk-1)
          cpt = (self%icoils(j)%coils(k)%pts(:,kk)+self%icoils(j)%coils(k)%pts(:,kk-1))/2.d0
          diffvec = diffvec + cross_product(cvec,pt_j-cpt)/SUM((pt_j-cpt)**2)**1.5d0
        END DO
      END IF
      ecc=ecc+self%icoils(j)%scales(k)*diffvec
    END DO
    DO jj=1,3
      !$omp atomic
      self%Bdr(i,j,jj) = self%Bdr(i,j,jj) + ecc(jj)
    END DO
  END DO
END DO
self%Bdr=self%Bdr*mu0/(4.d0*pi)
!
CALL oft_increase_indent
IF(PRESENT(save_file))CALL save_to_file()
CALL oft_decrease_indent
CONTAINS
!> Add the field of one basis function of a cell to the operator of an element, scattered through the junction
!! map if present
SUBROUTINE add_elem(atmp,ib,ik,s)
REAL(r8), INTENT(in) :: atmp(:,:,:) !< Field of the cell's basis functions at each vertex [3,3,np]
INTEGER(4), INTENT(in) :: ib !< Basis function of the cell
INTEGER(4), INTENT(in) :: ik !< Element (expanded if the model has junctions)
REAL(r8), INTENT(in) :: s !< Sign of the element's basis function
INTEGER(4) :: k
IF(ASSOCIATED(self%junc))THEN
  DO k=self%junc%p_kr(ik),self%junc%p_kr(ik+1)-1
    CALL add_col(atmp,ib,self%junc%p_lc(k),s*self%junc%p_val(k))
  END DO
ELSE
  CALL add_col(atmp,ib,ik,s)
END IF
END SUBROUTINE add_elem
!> Add a scaled basis function field to one column of the operator
SUBROUTINE add_col(atmp,ib,l,w)
REAL(r8), INTENT(in) :: atmp(:,:,:) !< Field of the cell's basis functions at each vertex [3,3,np]
INTEGER(4), INTENT(in) :: ib !< Basis function of the cell
INTEGER(4), INTENT(in) :: l !< Column (element) of the operator
REAL(r8), INTENT(in) :: w !< Weight
INTEGER(4) :: j,jj
REAL(r8) :: tmp
DO jj=1,3
  DO j=1,SIZE(atmp,DIM=3)
    tmp=w*atmp(jj,ib,j)
    !$omp atomic
    self%Bel(j,l,jj) = self%Bel(j,l,jj) + tmp
  END DO
END DO
END SUBROUTINE add_col
!> Load the operators from the cache file if its key matches the model
FUNCTION load_from_file() RESULT(exists)
INTEGER(4) :: hash_tmp(4),file_counts(4),jhash(1)
LOGICAL :: exists
IF(TRIM(save_file)/='none')THEN
  INQUIRE(FILE=TRIM(save_file),EXIST=exists)
  IF(exists)THEN
    CALL cache_key(hash_tmp)
    WRITE(*,'(3A)')oft_indent,'Loading B-field operator from file: ',TRIM(save_file)
    CALL hdf5_read(file_counts,TRIM(save_file),'MODEL_hash',success=exists)
    IF(exists)exists=ALL(file_counts==hash_tmp)
    IF(exists.AND.ASSOCIATED(self%junc))THEN
      CALL hdf5_read(jhash,TRIM(save_file),'JUNCTION_hash',success=exists)
      IF(exists)exists=(jhash(1)==tw_junction_hash(self))
    END IF
    IF(exists)THEN
      ALLOCATE(self%Bel(self%mesh%np,self%nelems,3))
      CALL hdf5_read(self%Bel(:,:,1),TRIM(save_file),'Bel_X',success=exists)
      IF(exists)CALL hdf5_read(self%Bel(:,:,2),TRIM(save_file),'Bel_Y',success=exists)
      IF(exists)CALL hdf5_read(self%Bel(:,:,3),TRIM(save_file),'Bel_Z',success=exists)
      IF(exists)THEN
        ALLOCATE(self%Bdr(self%mesh%np,self%n_icoils,3))
        CALL hdf5_read(self%Bdr(:,:,1),TRIM(save_file),'Bdr_X',success=exists)
        IF(exists)CALL hdf5_read(self%Bdr(:,:,2),TRIM(save_file),'Bdr_Y',success=exists)
        IF(exists)CALL hdf5_read(self%Bdr(:,:,3),TRIM(save_file),'Bdr_Z',success=exists)
      END IF
    END IF
    IF(exists)RETURN
    IF(ASSOCIATED(self%Bel))DEALLOCATE(self%Bel)
    IF(ASSOCIATED(self%Bdr))DEALLOCATE(self%Bdr)
  END IF
ELSE
  exists=.FALSE.
END IF
END FUNCTION load_from_file
!
!> Save the operators to the cache file
SUBROUTINE save_to_file()
INTEGER(4) :: hash_tmp(4)
IF(TRIM(save_file)/='none')THEN
  CALL cache_key(hash_tmp)
  WRITE(*,'(2A)')oft_indent,'Saving B-field operator to file: ',TRIM(save_file)
  CALL hdf5_create_file(TRIM(save_file))
  CALL hdf5_write(hash_tmp,TRIM(save_file),'MODEL_hash')
  IF(ASSOCIATED(self%junc))CALL hdf5_write([tw_junction_hash(self)],TRIM(save_file),'JUNCTION_hash')
  CALL hdf5_write(self%Bel(:,:,1),TRIM(save_file),'Bel_X')
  CALL hdf5_write(self%Bel(:,:,2),TRIM(save_file),'Bel_Y')
  CALL hdf5_write(self%Bel(:,:,3),TRIM(save_file),'Bel_Z')
  CALL hdf5_write(self%Bdr(:,:,1),TRIM(save_file),'Bdr_X')
  CALL hdf5_write(self%Bdr(:,:,2),TRIM(save_file),'Bdr_Y')
  CALL hdf5_write(self%Bdr(:,:,3),TRIM(save_file),'Bdr_Z')
END IF
END SUBROUTINE save_to_file
!> Cache key: element count and mesh hashes, with the count negated for a model with junctions (whose
!! map is keyed separately by the `JUNCTION_hash` dataset)
SUBROUTINE cache_key(hash_tmp)
INTEGER(4), INTENT(out) :: hash_tmp(4) !< Cache key
hash_tmp(1) = self%nelems
IF(ASSOCIATED(self%junc))hash_tmp(1) = -self%nelems
hash_tmp(2) = self%mesh%nc
hash_tmp(3) = oft_simple_hash(C_LOC(self%mesh%lc),INT(4*3*self%mesh%nc,8),starting_value=tw_coil_hash(self))
hash_tmp(4) = oft_simple_hash(C_LOC(self%mesh%r),INT(8*3*self%mesh%np,8))
END SUBROUTINE cache_key
END SUBROUTINE tw_compute_Bops
!---------------------------------------------------------------------------------
!> Setup hole definition for ordered chain of vertices
!---------------------------------------------------------------------------------
SUBROUTINE tw_setup_hole(bmesh,hmesh)
CLASS(oft_bmesh), INTENT(inout) :: bmesh !< Surface mesh containing hole
TYPE(hole_mesh), INTENT(inout) :: hmesh !< Hole definition
INTEGER(4) :: i,j,k,l,face,nneg,nposs,prev_orient
LOGICAL :: all_oriented
REAL(8) :: ptcc(3),evec(3),ecc(3),norm(3),f(3)
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: face_orientation,point_orient
DEBUG_STACK_PUSH
hmesh%ptcc = 0.d0
DO i=1,hmesh%n
  hmesh%ptcc = hmesh%ptcc + bmesh%r(:,hmesh%lp(i))
END DO
hmesh%ptcc = hmesh%ptcc/hmesh%n
ALLOCATE(hmesh%lp_sorted(hmesh%n),hmesh%sort_ind(hmesh%n))
hmesh%lp_sorted=hmesh%lp
hmesh%sort_ind=[(i,i=1,hmesh%n)]
CALL sort_array(hmesh%lp_sorted,hmesh%sort_ind,hmesh%n)
!---Determine face orientation
f=1.d0/3.d0
ALLOCATE(face_orientation(bmesh%nc),point_orient(hmesh%n))
point_orient=0
face_orientation=0
DO i=1,hmesh%n
  IF(i==hmesh%n)THEN
    k = ABS(mesh_local_findedge(bmesh, [hmesh%lp(i), hmesh%lp(1)]))
    evec=bmesh%r(:,hmesh%lp(1))-bmesh%r(:,hmesh%lp(i))
    ecc=(bmesh%r(:,hmesh%lp(1))+bmesh%r(:,hmesh%lp(i)))/2.d0
  ELSE
    k = ABS(mesh_local_findedge(bmesh, [hmesh%lp(i), hmesh%lp(i+1)]))
    evec=bmesh%r(:,hmesh%lp(i+1))-bmesh%r(:,hmesh%lp(i))
    ecc=(bmesh%r(:,hmesh%lp(i+1))+bmesh%r(:,hmesh%lp(i)))/2.d0
  END IF
  IF(k==0)CALL oft_abort('Could not find edge','tw_setup_hole',__FILE__)
  DO j=bmesh%kec(k),bmesh%kec(k+1)-1
    IF(face_orientation(bmesh%lec(j))/=0)CYCLE
    ptcc=(bmesh%r(:,bmesh%lc(1,bmesh%lec(j))) &
      + bmesh%r(:,bmesh%lc(2,bmesh%lec(j))) &
      + bmesh%r(:,bmesh%lc(3,bmesh%lec(j))))/3.d0
    CALL bmesh%norm(bmesh%lec(j),f,norm)
    face_orientation(bmesh%lec(j))=INT(SIGN(1.d0,DOT_PRODUCT(cross_product(ptcc-ecc,evec),norm)),4)
  END DO
  IF(bmesh%be(k))THEN
    IF(i==hmesh%n)THEN
      point_orient([1,i])=face_orientation(bmesh%lec(bmesh%kec(k)))
    ELSE
      point_orient(i:i+1)=face_orientation(bmesh%lec(bmesh%kec(k)))
    END IF
  END IF
END DO
IF(ALL(point_orient>=0))THEN
  point_orient=1.d0
ELSE IF(ALL(point_orient<0))THEN
  point_orient=-1.d0
ELSE
  prev_orient=0
  DO i=1,hmesh%n
    IF(point_orient(i)==0)THEN
      IF(prev_orient/=0)point_orient(i)=prev_orient
    ELSE
      prev_orient=point_orient(i)
    END IF
  END DO
  ! Wrap if necessary
  DO i=1,hmesh%n
    IF(point_orient(i)/=0)EXIT
    point_orient(i)=prev_orient
  END DO
END IF
!---Propogate
DO k=1,10
  all_oriented=.TRUE.
  DO i=1,hmesh%n
    DO j=bmesh%kpc(hmesh%lp(i)),bmesh%kpc(hmesh%lp(i)+1)-1
      IF(face_orientation(bmesh%lpc(j))==0)THEN
        DO l=1,3
          face=bmesh%lcc(l,bmesh%lpc(j))
          IF(face==0)CYCLE
          ! WRITE(*,*)i,j,k,l,face_orientation(face)
          IF(face_orientation(face)/=0)THEN
            face_orientation(bmesh%lpc(j))=face_orientation(face)
            EXIT
          END IF
        END DO
        IF(face_orientation(bmesh%lpc(j))==0)all_oriented=.FALSE.
      END IF
    END DO
  END DO
  IF(all_oriented)EXIT
END DO
IF(.NOT.all_oriented)CALL oft_abort("Error orienting cells on hole","tw_setup_hole",__FILE__)
!---Save orientations
ALLOCATE(hmesh%kpc(hmesh%n+1))
hmesh%kpc=0
DO i=1,hmesh%n
  DO j=bmesh%kpc(hmesh%lp(i)),bmesh%kpc(hmesh%lp(i)+1)-1
    !IF(ABS(face_orientation(bmesh%lpc(j)))>0)hmesh%kpc(i)=hmesh%kpc(i)+1
    IF(face_orientation(bmesh%lpc(j))==point_orient(i))hmesh%kpc(i)=hmesh%kpc(i)+1
  END DO
END DO
nneg=SUM(hmesh%kpc)
hmesh%kpc(hmesh%n+1)=nneg+1
do i=hmesh%n,1,-1 ! cumulative count
  hmesh%kpc(i) = hmesh%kpc(i+1) - hmesh%kpc(i)
end do
if(hmesh%kpc(1)/=1)call oft_abort('Bad element to element count','tw_setup_hole',__FILE__)
ALLOCATE(hmesh%lpc(nneg))
nneg=0
DO i=1,hmesh%n
  DO j=bmesh%kpc(hmesh%lp(i)),bmesh%kpc(hmesh%lp(i)+1)-1
    !IF(ABS(face_orientation(bmesh%lpc(j)))>0)THEN
    IF(face_orientation(bmesh%lpc(j))==point_orient(i))THEN
      nneg=nneg+1
      ! hmesh%lpc(nneg)=SIGN(bmesh%lpc(j),face_orientation(bmesh%lpc(j)))
      hmesh%lpc(nneg)=bmesh%lpc(j)*face_orientation(bmesh%lpc(j))
    END IF
  END DO
END DO
DEALLOCATE(face_orientation,point_orient)
DEBUG_STACK_POP
END SUBROUTINE tw_setup_hole
!------------------------------------------------------------------------------
!> Read coil sets for HDF5 mesh file
!------------------------------------------------------------------------------
subroutine tw_load_hdf5(self,filepath,hole_ns,error_string)
TYPE(tw_type), INTENT(inout) :: self !< ThinCurr object
CHARACTER(LEN=*), INTENT(in) :: filepath !< Path to mesh file
TYPE(oft_1d_int), POINTER, INTENT(inout) :: hole_ns(:)
CHARACTER(LEN=:), ALLOCATABLE, INTENT(out) :: error_string
integer(4) :: ncoils,ncoil_sets,nread,coil_type,nmasked
LOGICAL :: success
INTEGER(4) :: i,j,k,io_unit,ierr,id,cell,ipath,ndims,n_jumpers,npts_tmp
INTEGER(4), ALLOCATABLE :: dim_sizes(:),itmp2(:,:)
INTEGER(4) :: psize(2)
REAL(8) :: res_per_len,radius,dl,theta,rz_tmp(2)
REAL(8), POINTER :: pts_tmp(:)
CHARACTER(LEN=4) :: blknum,blknum2
CHARACTER(LEN=16) :: coil_ind
CHARACTER(LEN=OFT_PATH_SLEN) :: coil_path
CHARACTER(LEN=:), ALLOCATABLE :: str_tmp
TYPE(tw_coil_set), POINTER :: coil_tmp
nmasked=0
CALL oft_increase_indent
!---Count coil sets
IF(hdf5_field_exist(TRIM(filepath),'thincurr/coils'))THEN
  WRITE(*,'(2A)')oft_indent,'Loading driver coils from mesh file'
  CALL hdf5_read(ncoil_sets,TRIM(filepath),'thincurr/coils/NCOILSETS',success=success)
  self%n_coils=ncoil_sets
  ALLOCATE(self%coils(self%n_coils))
  self%n_icoils=0
  self%n_vcoils=0
  DO i=1,ncoil_sets
    WRITE(blknum,'(I4.4)')i
    !
    CALL hdf5_read(ncoils,TRIM(filepath),'thincurr/coils/coilset'//blknum//'/NCOILS',success=success)
    IF(.NOT.success)CALL oft_abort('Failed to read coils for coilset '//blknum,'tw_load_coils_hdf5',__FILE__)
    coil_tmp=>self%coils(i)
    coil_tmp%ncoils=ncoils
    ALLOCATE(coil_tmp%coils(coil_tmp%ncoils))
    ALLOCATE(coil_tmp%scales(coil_tmp%ncoils))
    ALLOCATE(coil_tmp%res_per_len(coil_tmp%ncoils))
    ALLOCATE(coil_tmp%radius(coil_tmp%ncoils))
    coil_tmp%scales=1.d0
    coil_tmp%res_per_len=-1.d0
    coil_tmp%radius=-1.d0
    coil_tmp%Rself=0.d0
    !
    IF(hdf5_field_exist(TRIM(filepath),'thincurr/coils/coilset'//blknum//'/NAME'))THEN
      CALL hdf5_read(str_tmp,TRIM(filepath),'thincurr/coils/coilset'//blknum//'/NAME',success)
      IF(.NOT.success)CALL oft_abort('Failed to read name for coilset '//blknum,'tw_load_coils_hdf5',__FILE__)
      IF(LEN(str_tmp)>LEN(coil_tmp%name))THEN
        WRITE(*,'(3A)')oft_indent,'Coil set name "',str_tmp,'" too long for coil set truncating'
        coil_tmp%name=str_tmp(1:LEN(coil_tmp%name))
      ELSE
        coil_tmp%name=str_tmp
      END IF
    ELSE
      coil_tmp%name='COIL_'//blknum
    END IF
    IF(hdf5_field_exist(TRIM(filepath),'thincurr/coils/coilset'//blknum//'/RADIUS'))THEN
      CALL hdf5_read(coil_tmp%radius,TRIM(filepath),'thincurr/coils/coilset'//blknum//'/RADIUS',success)
      IF(.NOT.success)CALL oft_abort('Failed to read radius for coilset '//blknum,'tw_load_coils_hdf5',__FILE__)
    END IF
    IF(hdf5_field_exist(TRIM(filepath),'thincurr/coils/coilset'//blknum//'/RES_PER_LEN'))THEN
      CALL hdf5_read(coil_tmp%res_per_len,TRIM(filepath),'thincurr/coils/coilset'//blknum//'/RES_PER_LEN',success)
      IF(.NOT.success)CALL oft_abort('Failed to read res_per_len for coilset '//blknum,'tw_load_coils_hdf5',__FILE__)
    END IF
    IF(hdf5_field_exist(TRIM(filepath),'thincurr/coils/coilset'//blknum//'/SCALES'))THEN
      CALL hdf5_read(coil_tmp%scales,TRIM(filepath),'thincurr/coils/coilset'//blknum//'/SCALES',success)
      IF(.NOT.success)CALL oft_abort('Failed to read scales for coilset '//blknum,'tw_load_coils_hdf5',__FILE__)
    END IF
    IF(hdf5_field_exist(TRIM(filepath),'thincurr/coils/coilset'//blknum//'/SENS_MASK'))THEN
      CALL hdf5_read(coil_tmp%sens_mask,TRIM(filepath),'thincurr/coils/coilset'//blknum//'/SENS_MASK',success)
      IF(.NOT.success)CALL oft_abort('Failed to read sens_mask for coilset '//blknum,'tw_load_coils_hdf5',__FILE__)
      IF(coil_tmp%sens_mask)THEN
        nmasked=nmasked+1
        IF(oft_debug_print(2))WRITE(*,'(2A,I6,A)')oft_indent,'Masked coil set ',i,' from sensors'
      END IF
    END IF
    !
    DO j=1,ncoils
      WRITE(blknum2,'(I4.4)')j
      coil_path='thincurr/coils/coilset'//blknum//'/coil'//blknum2
      !---Circular coil: radius and height (and optionally the number of polygon points)
      IF(hdf5_field_exist(TRIM(filepath),TRIM(coil_path)//'/RZ'))THEN
        CALL hdf5_read(rz_tmp,TRIM(filepath),TRIM(coil_path)//'/RZ',success=success)
        npts_tmp=0
        IF(success.AND.hdf5_field_exist(TRIM(filepath),TRIM(coil_path)//'/NPTS'))THEN
          CALL hdf5_read(npts_tmp,TRIM(filepath),TRIM(coil_path)//'/NPTS',success=success)
        END IF
        IF(.NOT.success)THEN
          WRITE(coil_ind,'(I6,2X,I6)')i,j
          CALL oft_abort('Failed to read circular coil for coilset/coil '//coil_ind,'tw_load_coils_hdf5',__FILE__)
        END IF
        CALL tw_coil_set_ring(coil_tmp%coils(j),rz_tmp,npts_tmp)
        CYCLE
      END IF
      CALL hdf5_field_get_sizes(TRIM(filepath),'thincurr/coils/coilset'//blknum//'/coil'//blknum2//'/PTS',ndims,dim_sizes)
      IF(ndims==2)THEN
        coil_tmp%coils(j)%npts=dim_sizes(2)
        ALLOCATE(coil_tmp%coils(j)%pts(3,coil_tmp%coils(j)%npts))
        CALL hdf5_read(coil_tmp%coils(j)%pts,TRIM(filepath),'thincurr/coils/coilset'//blknum//'/coil'//blknum2//'/PTS',success)
        IF(.NOT.success)ndims=-1
      END IF
      IF(ndims/=2)THEN
        WRITE(coil_ind,'(I6,2X,I6)')i,j
        CALL oft_abort('Failed to read coil points for coilset/coil '//coil_ind,'tw_load_coils_hdf5',__FILE__)
      END IF
    END DO
    !---TEMPORARY to mirror old behavior of separating vcoils and icoils
    IF(ALL(coil_tmp%res_per_len>0.d0).AND.ALL(coil_tmp%radius>0.d0))THEN
      self%n_vcoils=self%n_vcoils+1
    ELSE
      self%n_icoils=self%n_icoils+1
    END IF
  END DO
  !---TEMPORARY to mirror old behavior of separating vcoils and icoils
  IF(self%n_vcoils>0)ALLOCATE(self%vcoils(self%n_vcoils))
  IF(self%n_icoils>0)ALLOCATE(self%icoils(self%n_icoils))
  self%n_icoils=0
  self%n_vcoils=0
  DO i=1,ncoil_sets
    IF(ALL(self%coils(i)%res_per_len>0.d0).AND.ALL(self%coils(i)%radius>0.d0))THEN
      self%n_vcoils=self%n_vcoils+1
      self%coils(i)%is_vcoil=.TRUE.
      CALL tw_copy_coil(self%coils(i),self%vcoils(self%n_vcoils))
    ELSE
      self%n_icoils=self%n_icoils+1
      self%coils(i)%is_vcoil=.FALSE.
      CALL tw_copy_coil(self%coils(i),self%icoils(self%n_icoils))
    END IF
  END DO
END IF
!---
IF(nmasked>0)WRITE(*,'(2A,I6,A)')oft_indent,'Masked ',nmasked,' coils from sensors'
IF(oft_debug_print(1))THEN
  IF(self%n_vcoils>0)THEN
    WRITE(*,*)
    WRITE(*,'(2A)')oft_indent,'V-Coils set definitions'
    WRITE(*,'(2A)')oft_indent,'========================='
    WRITE(*,'(2A,I6,A)')oft_indent,'Found ',self%n_vcoils,' coil sets'
    DO i=1,self%n_vcoils
      WRITE(*,'(2A,I6,A)')oft_indent,'  Set  : ',i
      WRITE(*,'(2A,I6,A)')oft_indent,'    nCoils  : ',self%vcoils(i)%ncoils
    END DO
    WRITE(*,'(2A)')oft_indent,'========================='
    WRITE(*,*)
  END IF
  IF(self%n_icoils>0)THEN
    WRITE(*,*)
    WRITE(*,'(2A)')oft_indent,'I-Coils set definitions'
    WRITE(*,'(2A)')oft_indent,'========================='
    WRITE(*,'(2A,I6,A)')oft_indent,'Found ',self%n_icoils,' coil sets'
    DO i=1,self%n_icoils
      WRITE(*,'(2A,I6,A)')oft_indent,'  Set  : ',i
      WRITE(*,'(2A,I6,A)')oft_indent,'    nCoils  : ',self%icoils(i)%ncoils
    END DO
    WRITE(*,'(2A)')oft_indent,'========================='
    WRITE(*,*)
  END IF
END IF
CALL oft_decrease_indent
ALLOCATE(self%sens_mask(self%mesh%nreg))
self%sens_mask=.FALSE.
!---Read periodicity information
IF(hdf5_field_exist(TRIM(filepath),'thincurr/periodicity/PMAP'))THEN
  IF(ASSOCIATED(self%pmap))THEN
    CALL oft_warn("Point map specified directly and in mesh file, ignoring mesh file definition")
  ELSE
    ALLOCATE(self%pmap(self%mesh%np))
    CALL hdf5_read(self%pmap,TRIM(filepath),'thincurr/periodicity/PMAP',success)
    IF(.NOT.success)THEN
      error_string='Error reading periodicity information from mesh file'
      RETURN
    END IF
  END IF
END IF
!---Read surface resistivity
IF(hdf5_field_exist(TRIM(filepath),'thincurr/ETA_SURF'))THEN
  ALLOCATE(self%Eta_surf(self%mesh%nreg))
  CALL hdf5_read(self%Eta_surf,TRIM(filepath),'thincurr/ETA_SURF',success)
  IF(.NOT.success)THEN
    error_string='Error reading surface resistivity from mesh file'
    RETURN
  END IF
  self%Eta_surf=self%Eta_surf/mu0
END IF
!---Read volume resistivity
IF(hdf5_field_exist(TRIM(filepath),'thincurr/ETA_VOL'))THEN
  IF(ASSOCIATED(self%Eta_surf))THEN
    error_string='Volume resistivity and surface resistivity cannot both be present in mesh file'
    RETURN
  END IF
  ALLOCATE(self%Eta_vol(self%mesh%nreg))
  CALL hdf5_read(self%Eta_vol,TRIM(filepath),'thincurr/ETA_VOL',success)
  IF(.NOT.success)THEN
    error_string='Error reading volume resistivity from mesh file'
    RETURN
  END IF
  self%Eta_vol=self%Eta_vol/mu0
END IF
!---Read thickness
IF(hdf5_field_exist(TRIM(filepath),'thincurr/THICKNESS'))THEN
  ALLOCATE(self%thickness(self%mesh%nreg))
  CALL hdf5_read(self%thickness,TRIM(filepath),'thincurr/THICKNESS',success)
  IF(.NOT.success)THEN
    error_string='Error reading thickness from mesh file'
    RETURN
  END IF
  IF(ASSOCIATED(self%Eta_vol))THEN
    ALLOCATE(self%Eta_surf(self%mesh%nreg))
    self%Eta_surf=self%Eta_vol/self%thickness
  END IF
ELSE
  IF(ASSOCIATED(self%Eta_vol))THEN
    error_string='Volume resistivity requires thickness to be present in mesh file'
    RETURN
  END IF
END IF
!---
IF(hdf5_field_exist(TRIM(filepath),'thincurr/holes'))THEN
  IF(ASSOCIATED(hole_ns))THEN
    CALL oft_warn("Holes specified in mesh and XML files, ignoring XML definitions")
    DO i=1,SIZE(hole_ns)
      IF(ASSOCIATED(hole_ns(i)%v))DEALLOCATE(hole_ns(i)%v)
    END DO
    DEALLOCATE(hole_ns)
  END IF
  CALL hdf5_read(self%nholes,TRIM(filepath),'thincurr/holes/NHOLES',success)
  IF(.NOT.success)THEN
    error_string='Error reading number of holes from mesh file'
    RETURN
  END IF
  ALLOCATE(hole_ns(self%nholes))
  DO i=1,self%nholes
    WRITE(blknum,'(I4.4)')i
    CALL hdf5_field_get_sizes(TRIM(filepath),'thincurr/holes/HOLE'//blknum,ndims,dim_sizes)
    IF(ndims==1)THEN
      hole_ns(i)%n=dim_sizes(1)
      ALLOCATE(hole_ns(i)%v(hole_ns(i)%n))
      CALL hdf5_read(hole_ns(i)%v,TRIM(filepath),'thincurr/holes/HOLE'//blknum,success)
      IF(.NOT.success)ndims=-1
    ELSE
      ndims=-1
    END IF
    IF(ndims<0)THEN
      error_string='Error reading hole nodeset from mesh file'
      RETURN
    END IF
  END DO
END IF
!---Junctions
IF(hdf5_field_exist(TRIM(filepath),'thincurr/junctions'))THEN
  ALLOCATE(self%junc)
  IF(ASSOCIATED(self%pmap))THEN
    error_string='Periodic point maps ("pmap") are not supported for meshes split at junctions'
    RETURN
  END IF
  ALLOCATE(self%pmap(self%mesh%np))
  CALL hdf5_read(self%pmap,TRIM(filepath),'thincurr/junctions/PMAP',success)
  IF(.NOT.success)THEN
    error_string='Error reading junction point map from mesh file'
    RETURN
  END IF
  CALL hdf5_read(self%junc%nholes,TRIM(filepath),'thincurr/junctions/NHOLES',success)
  ALLOCATE(self%junc%holes(self%junc%nholes))
  DO i=1,self%junc%nholes
    WRITE(blknum,'(I4.4)')i
    CALL hdf5_field_get_sizes(TRIM(filepath),'thincurr/junctions/HOLE'//blknum,ndims,dim_sizes)
    ALLOCATE(itmp2(dim_sizes(1),dim_sizes(2)))
    CALL hdf5_read(itmp2,TRIM(filepath),'thincurr/junctions/HOLE'//blknum,success)
    IF(.NOT.success)THEN
      error_string='Error reading junction hole from mesh file'
      RETURN
    END IF
    IF(SIZE(itmp2,DIM=1)/=3)THEN
      error_string='Junction holes must list (cell, vertex, sign) triples'
      RETURN
    END IF
    self%junc%holes(i)%n=SIZE(itmp2,DIM=2)
    ALLOCATE(self%junc%holes(i)%v(3*self%junc%holes(i)%n))
    self%junc%holes(i)%v=RESHAPE(itmp2,[3*self%junc%holes(i)%n])
    DEALLOCATE(itmp2)
  END DO
  CALL hdf5_read(psize,TRIM(filepath),'thincurr/junctions/P_SIZE',success)
  self%junc%nred=psize(2)
  ALLOCATE(self%junc%p_kr(psize(1)+1))
  CALL hdf5_read(self%junc%p_kr,TRIM(filepath),'thincurr/junctions/P_KR',success)
  ALLOCATE(self%junc%p_lc(self%junc%p_kr(psize(1)+1)-1),self%junc%p_val(self%junc%p_kr(psize(1)+1)-1))
  CALL hdf5_read(self%junc%p_lc,TRIM(filepath),'thincurr/junctions/P_LC',success)
  CALL hdf5_read(self%junc%p_val,TRIM(filepath),'thincurr/junctions/P_VAL',success)
  CALL hdf5_field_get_sizes(TRIM(filepath),'thincurr/junctions/LJE',ndims,dim_sizes)
  ALLOCATE(self%junc%table(dim_sizes(1),dim_sizes(2)))
  CALL hdf5_read(self%junc%table,TRIM(filepath),'thincurr/junctions/LJE',success)
  IF(.NOT.success)THEN
    error_string='Error reading junction map from mesh file'
    RETURN
  END IF
END IF
!---
IF(hdf5_field_exist(TRIM(filepath),'thincurr/jumpers'))THEN
  IF(ASSOCIATED(self%jumper_nsets))THEN
    CALL oft_warn("Jumpers specified in mesh and XML files, ignoring XML definitions")
    DO i=1,SIZE(self%jumper_nsets)
      IF(ASSOCIATED(self%jumper_nsets(i)%v))DEALLOCATE(self%jumper_nsets(i)%v)
    END DO
    DEALLOCATE(self%jumper_nsets)
  END IF
  CALL hdf5_read(n_jumpers,TRIM(filepath),'thincurr/jumpers/NJUMPERS',success)
  IF(.NOT.success)THEN
    error_string='Error reading number of jumpers from mesh file'
    RETURN
  END IF
  ALLOCATE(self%jumper_nsets(n_jumpers))
  DO i=1,n_jumpers
    WRITE(blknum,'(I4.4)')i
    CALL hdf5_field_get_sizes(TRIM(filepath),'thincurr/jumpers/JUMPER'//blknum,ndims,dim_sizes)
    IF(ndims==1)THEN
      self%jumper_nsets(i)%n=dim_sizes(1)
      ALLOCATE(self%jumper_nsets(i)%v(self%jumper_nsets(i)%n))
      CALL hdf5_read(self%jumper_nsets(i)%v,TRIM(filepath),'thincurr/jumpers/JUMPER'//blknum,success)
      IF(.NOT.success)ndims=-1
    ELSE
      ndims=-1
    END IF
    IF(ndims<0)THEN
      error_string='Error reading jumper nodeset from mesh file'
      RETURN
    END IF
  END DO
END IF
!---Read hole construction version
IF(hdf5_field_exist(TRIM(filepath),'thincurr/HOLES_VERSION'))THEN
  CALL hdf5_read(self%holes_version,TRIM(filepath),'thincurr/HOLES_VERSION',success)
  IF(.NOT.success)THEN
    error_string='Error reading hole version from mesh file'
    RETURN
  END IF
END IF
!---Read closures
IF(hdf5_field_exist(TRIM(filepath),'thincurr/CLOSURES'))THEN
  CALL hdf5_field_get_sizes(TRIM(filepath),'thincurr/CLOSURES',ndims,dim_sizes)
  IF(ndims==1)THEN
    IF(ASSOCIATED(self%closures))THEN
      CALL oft_warn("Closures specified in mesh and XML files, ignoring XML definitions")
      DEALLOCATE(self%closures)
    END IF
    self%nclosures=dim_sizes(1)
    ALLOCATE(self%closures(self%nclosures))
    CALL hdf5_read(self%closures,TRIM(filepath),'thincurr/CLOSURES',success)
    IF(.NOT.success)ndims=-1
  ELSE
    ndims=-1
  END IF
  IF(ndims<0)THEN
    error_string='Error reading closures from mesh file'
    RETURN
  END IF
END IF
end subroutine tw_load_hdf5
!------------------------------------------------------------------------------
!> Read coil sets for "oft_in.xml" file
!------------------------------------------------------------------------------
subroutine tw_load_coils_xml(group_node,name_offset,ncoils,coils)
TYPE(xml_node), INTENT(IN) :: group_node !< XML node relative to base `<thincurr>` node
INTEGER(i4), INTENT(IN) :: name_offset !< Offset for coil set names
INTEGER(4), INTENT(out) :: ncoils !< Number of coil sets found
TYPE(tw_coil_set), POINTER, INTENT(out) :: coils(:) !< List of coil sets
!---XML solver fields
integer(4) :: ncoil_sets,nread,nmasked,npts_tmp
TYPE(xml_node) :: doc,thincurr_group
TYPE(xml_node), POINTER :: coil_set,coil
TYPE(xml_nodelist) :: coil_sets,coil_list
!---
LOGICAL :: success
INTEGER(4) :: i,j,io_unit,ierr,id,cell,ipath,ndims
INTEGER(4), ALLOCATABLE :: dim_sizes(:)
REAL(8) :: res_per_len,radius,dl
REAL(8), POINTER :: pts2_tmp(:,:)
CHARACTER(LEN=16) :: coil_ind
CHARACTER(LEN=OFT_PATH_SLEN) :: coil_path
CHARACTER(LEN=:), ALLOCATABLE :: str_tmp
TYPE(tw_coil_set), POINTER :: coil_tmp
!---Count coil sets
ncoils=0
CALL xml_get_element(group_node,"coil_set",coil_sets,ierr)
IF(ierr==0)ncoils=coil_sets%n
ALLOCATE(coils(ncoils))
IF(ncoils==0)RETURN
nmasked=0
CALL oft_increase_indent
!---Setup coil sets
DO i=1,ncoils
  coil_tmp=>coils(i)
  coil_set=>coil_sets%nodes(i)!%this
  !
  CALL xml_get_element(coil_set,"coil",coil_list,ierr)
  IF(ierr/=0)CYCLE
  coil_tmp%ncoils=coil_list%n
  ALLOCATE(coil_tmp%scales(coil_tmp%ncoils))
  ALLOCATE(coil_tmp%res_per_len(coil_tmp%ncoils))
  ALLOCATE(coil_tmp%radius(coil_tmp%ncoils))
  coil_tmp%scales=1.d0
  coil_tmp%res_per_len=-1.d0
  coil_tmp%radius=-1.d0
  coil_tmp%Rself=0.d0
  !---Get coil set name
  IF(xml_hasAttribute(coil_set,"name"))THEN
    CALL xml_read_attribute(coil_set,"name",str_tmp,iostat=ierr)
    IF(ierr/=0)THEN
      WRITE(coil_ind,'(I6)')i
      CALL oft_xml_abort('Error reading "name" for coil set '//coil_ind,'tw_load_coils_xml',__FILE__)
    END IF
    IF(LEN(str_tmp)>LEN(coil_tmp%name))THEN
      WRITE(*,'(2A,I6,A)')oft_indent,'Coil set name too long for coil set ',i,', truncating'
      coil_tmp%name=str_tmp(1:LEN(coil_tmp%name))
    ELSE
      coil_tmp%name=str_tmp
    END IF
    DEALLOCATE(str_tmp)
  ELSE
    WRITE(coil_tmp%name,'(A5,I4.4)')'COIL_',i+name_offset
  END IF
  !---Get coil set resistivity per unit length (can be overriden)
  IF(xml_hasAttribute(coil_set,"res_per_len"))THEN
    CALL xml_read_attribute(coil_set,"res_per_len",res_per_len,iostat=ierr)
    IF(ierr/=0)THEN
      WRITE(coil_ind,'(I6)')i
      CALL oft_xml_abort('Error reading "res_per_len" for coil set '//coil_ind,'tw_load_coils_xml',__FILE__)
    END IF
    coil_tmp%res_per_len=res_per_len
  END IF
  !---Get coil set radius (can be overriden)
  IF(xml_hasAttribute(coil_set,"radius"))THEN
    CALL xml_read_attribute(coil_set,"radius",radius,iostat=ierr)
    IF(ierr/=0)THEN
      WRITE(coil_ind,'(I6)')i
      CALL oft_xml_abort('Error reading "radius" for coil set '//coil_ind,'tw_load_coils_xml',__FILE__)
    END IF
    coil_tmp%radius=radius
  END IF
  !---Get sensor flag
  IF(xml_hasAttribute(coil_set,"sens_mask"))THEN
    CALL xml_read_attribute(coil_set,"sens_mask",coil_tmp%sens_mask,iostat=ierr)
    IF(ierr/=0)THEN
      WRITE(coil_ind,'(I6)')i
      CALL oft_xml_abort('Error reading "sens_mask" for coil set '//coil_ind,'tw_load_coils_xml',__FILE__)
    END IF
    IF(coil_tmp%sens_mask)THEN
      IF(oft_debug_print(2))WRITE(*,'(2A,I6,A)')oft_indent,'Masking coil ',i,' from sensors'
      nmasked=nmasked+1
    END IF
  END IF
  ALLOCATE(coil_tmp%coils(coil_tmp%ncoils))
  DO j=1,coil_tmp%ncoils
    coil=>coil_list%nodes(j)!%this
    !---Look for HDF5 path
    IF(xml_hasAttribute(coil,"path"))THEN
      CALL xml_read_attribute(coil,"path",str_tmp,iostat=ierr)
      IF(ierr/=0)THEN
        WRITE(coil_ind,'(I6,2X,I6)')i,j
        CALL oft_xml_abort('Error reading "path" in coil '//coil_ind,'tw_load_coils_xml',__FILE__)
      END IF
      ipath=INDEX(str_tmp,":")
      IF(ipath==0)THEN
        WRITE(coil_ind,'(I6,2X,I6)')i,j
        CALL oft_abort('Misformatted "path" attribute in coil '//coil_ind,'tw_load_coils_xml',__FILE__)
      END IF
      CALL hdf5_field_get_sizes(str_tmp(1:ipath-1),str_tmp(ipath+1:),ndims,dim_sizes)
      IF(ndims<0)THEN
        WRITE(coil_ind,'(I6,2X,I6)')i,j
        CALL oft_abort('Failed to read HDF5 data sizes for coil '//coil_ind,'tw_load_coils_xml',__FILE__)
      END IF
      IF(dim_sizes(1)/=3)THEN
        WRITE(coil_ind,'(I6,2X,I6)')i,j
        CALL oft_abort('Incorrect first dimension of HDF5 dataset for coil '//coil_ind,'tw_load_coils_xml',__FILE__)
      END IF
      coil_tmp%coils(j)%npts=dim_sizes(2)
      DEALLOCATE(dim_sizes)
      ALLOCATE(coil_tmp%coils(j)%pts(3,coil_tmp%coils(j)%npts))
      CALL hdf5_read(coil_tmp%coils(j)%pts,str_tmp(1:ipath-1),str_tmp(ipath+1:),success=success)
      IF(.NOT.success)THEN
        WRITE(coil_ind,'(I6,2X,I6)')i,j
        CALL oft_abort('Failed to read HDF5 data for coil '//coil_ind,'tw_load_coils_xml',__FILE__)
      END IF
      DEALLOCATE(str_tmp)
    ELSE
      !---Read number of points (optional; polygon resolution for a circular coil)
      npts_tmp=0
      IF(xml_hasAttribute(coil,"npts"))THEN
        CALL xml_read_attribute(coil,"npts",npts_tmp,iostat=ierr)
        IF(ierr/=0)THEN
          WRITE(coil_ind,'(I6,2X,I6)')i,j
          CALL oft_xml_abort('Error reading "npts" for coil '//coil_ind,'tw_load_coils_xml',__FILE__)
        END IF
      END IF
      !---Points: a single (R, Z) pair for a circular coil, otherwise rows of (x, y, z)
      CALL xml_read_content(coil,pts2_tmp,iostat=ierr)
      IF(ierr/=0)THEN
        WRITE(coil_ind,'(I6,2X,I6)')i,j
        CALL oft_xml_abort('Error reading coil '//coil_ind,'tw_load_coils_xml',__FILE__)
      END IF
      IF(SIZE(pts2_tmp)==2)THEN
        CALL tw_coil_set_ring(coil_tmp%coils(j),RESHAPE(pts2_tmp,[2]),npts_tmp)
        DEALLOCATE(pts2_tmp)
      ELSE
        IF(SIZE(pts2_tmp,1)/=3)THEN
          WRITE(coil_ind,'(I6,2X,I6)')i,j
          CALL oft_abort('Incorrect first dimension of coil points in coil '//coil_ind,'tw_load_coils_xml',__FILE__)
        END IF
        IF((npts_tmp>0).AND.(SIZE(pts2_tmp,2)/=npts_tmp))THEN
          WRITE(coil_ind,'(I6,2X,I6)')i,j
          CALL oft_abort('Incorrect second dimension of coil points in coil '//coil_ind,'tw_load_coils_xml',__FILE__)
        END IF
        coil_tmp%coils(j)%npts=SIZE(pts2_tmp,2)
        coil_tmp%coils(j)%pts=>pts2_tmp
      END IF
    END IF
    !---Get scale factor
    IF(xml_hasAttribute(coil,"scale"))THEN
      CALL xml_read_attribute(coil,"scale",coil_tmp%scales(j),iostat=ierr)
      IF(ierr/=0)THEN
        WRITE(coil_ind,'(I6,2X,I6)')i,j
        CALL oft_xml_abort('Error reading "scale" for coil '//coil_ind,'tw_load_coils_xml',__FILE__)
      END IF
    END IF
    !---Get coil resistivity per unit length
    IF(xml_hasAttribute(coil,"res_per_len"))THEN
      CALL xml_read_attribute(coil,"res_per_len",coil_tmp%res_per_len(j),iostat=ierr)
      IF(ierr/=0)THEN
        WRITE(coil_ind,'(I6,2X,I6)')i,j
        CALL oft_xml_abort('Error reading "res_per_len" for coil '//coil_ind,'tw_load_coils_xml',__FILE__)
      END IF
    END IF
    !---Get coil radius
    IF(xml_hasAttribute(coil,"radius"))THEN
      CALL xml_read_attribute(coil,"radius",coil_tmp%radius(j),iostat=ierr)
      IF(ierr/=0)THEN
        WRITE(coil_ind,'(I6,2X,I6)')i,j
        CALL oft_xml_abort('Error reading "radius" for coil '//coil_ind,'tw_load_coils_xml',__FILE__)
      END IF
    END IF
  END DO
  IF(ASSOCIATED(coil_list%nodes))DEALLOCATE(coil_list%nodes)
END DO
IF(ASSOCIATED(coil_sets%nodes))DEALLOCATE(coil_sets%nodes)
!---
IF(nmasked>0)WRITE(*,'(2A,I6,A)')oft_indent,'Masked ',nmasked,' coils from sensors'
IF(oft_debug_print(1))THEN
  WRITE(*,*)
  WRITE(*,'(2A)')oft_indent,'Coils set definitions'
  WRITE(*,'(2A)')oft_indent,'========================='
  WRITE(*,'(2A,I6,A)')oft_indent,'Found ',ncoils,' coil sets'
  DO i=1,ncoils
    WRITE(*,'(2A,I6,A)')oft_indent,'  Set  : ',i
    WRITE(*,'(2A,I6,A)')oft_indent,'    nCoils  : ',coils(i)%ncoils
  END DO
  WRITE(*,'(2A)')oft_indent,'========================='
  WRITE(*,*)
END IF
CALL oft_decrease_indent
end subroutine tw_load_coils_xml
!------------------------------------------------------------------------------
!> Create a copy (by reference) of a coil set
!------------------------------------------------------------------------------
SUBROUTINE tw_copy_coil(coil_in,coil_out)
TYPE(tw_coil_set), INTENT(in) :: coil_in !< Source coil set
TYPE(tw_coil_set), INTENT(inout) :: coil_out !< Copy of source
! coil_out%type=coil_in%type
coil_out%name=coil_in%name
coil_out%ncoils=coil_in%ncoils
coil_out%Lself=coil_in%Lself
coil_out%Rself=coil_in%Rself
coil_out%scales=>coil_in%scales
coil_out%res_per_len=>coil_in%res_per_len
coil_out%radius=>coil_in%radius
! coil_out%axi_pt=>coil_in%axi_pt
coil_out%coils=>coil_in%coils
coil_out%sens_mask=coil_in%sens_mask
END SUBROUTINE tw_copy_coil
!------------------------------------------------------------------------------
!> Load sensors from "floops.loc" and build jumpers from nodesets
!------------------------------------------------------------------------------
subroutine tw_load_sensors(filename,self,sensors)
CHARACTER(LEN=*), INTENT(in) :: filename !< Thin-wall model object
CLASS(tw_type), INTENT(inout) :: self !< Thin-wall model object
TYPE(tw_sensors), INTENT(inout) :: sensors !< Sensor container
INTEGER(4) :: i,j,k,io_unit,ierr,id,vert,cell,p1,p2
INTEGER(4), ALLOCATABLE :: ed_mark(:),list_out(:),hole_tmp(:,:)
REAL(8) :: location(2),norm(3),ed_norm(3),cell_norm(3),f(3)
REAL(8), ALLOCATABLE :: hole_facs(:)
LOGICAL :: exists
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Loading sensor information'
CALL oft_increase_indent
!---Load flux loops
sensors%nfloops=0
INQUIRE(FILE=TRIM(filename), EXIST=exists)
IF(exists)THEN
  WRITE(*,'(2A)')oft_indent,'Loading flux loops from file: '//TRIM(filename)
  OPEN(NEWUNIT=io_unit, FILE=TRIM(filename))
  ierr=skip_comment_lines(io_unit)
  READ(io_unit,*)sensors%nfloops
  WRITE(*,'(2A,I6)')oft_indent,'  # of floops =',sensors%nfloops
  ALLOCATE(sensors%floops(sensors%nfloops))
  DO i=1,sensors%nfloops
    READ(io_unit,*)
    READ(io_unit,*)sensors%floops(i)%np,sensors%floops(i)%scale_fac,sensors%floops(i)%name
    IF(oft_debug_print(1))WRITE(*,'(2A,I6)')oft_indent,'    # of floops pts =',sensors%floops(i)%np
    ALLOCATE(sensors%floops(i)%r(3,sensors%floops(i)%np))
    DO j=1,sensors%floops(i)%np
      READ(io_unit,*)sensors%floops(i)%r(:,j)
    END DO
  END DO
  CLOSE(io_unit)
END IF
!---Load jumpers
sensors%njumpers=0
IF(ASSOCIATED(self%jumper_nsets))THEN
  WRITE(*,'(2A)')oft_indent,'Setting up jumpers'
  OPEN(NEWUNIT=io_unit, FILE='jumpers_orient.dat')
  sensors%njumpers=SIZE(self%jumper_nsets)
  ALLOCATE(sensors%jumpers(sensors%njumpers))
  ALLOCATE(hole_facs(self%nholes))
  DO i=1,sensors%njumpers
    WRITE(sensors%jumpers(i)%name,'(A,I4.4)')'JUMPER_',i
    sensors%jumpers(i)%np=self%jumper_nsets(i)%n
    IF(ASSOCIATED(self%junc))THEN
      ! Already ordered on the split mesh, with segments separated by 0 where the jumper crosses a junction
      ALLOCATE(sensors%jumpers(i)%points(sensors%jumpers(i)%np))
      sensors%jumpers(i)%points=self%jumper_nsets(i)%v
    ELSE
      ALLOCATE(sensors%jumpers(i)%points(sensors%jumpers(i)%np+1))
      CALL order_jumper_list(self%jumper_nsets(i)%v,sensors%jumpers(i)%points,self%jumper_nsets(i)%n)
      IF(sensors%jumpers(i)%points(sensors%jumpers(i)%np+1)>0)sensors%jumpers(i)%np=sensors%jumpers(i)%np+1
    END IF
    ! Setup hole couplings
    hole_facs=0.d0
    DO j=1,sensors%jumpers(i)%np-1
      IF(sensors%jumpers(i)%points(j)==0.OR.sensors%jumpers(i)%points(j+1)==0)CYCLE
      id=mesh_local_findedge(self%mesh,[sensors%jumpers(i)%points(j),sensors%jumpers(i)%points(j+1)])
      IF(id==0)THEN
        WRITE(*,*)"No Edge matching points: ",sensors%jumpers(i)%points(j),sensors%jumpers(i)%points(j+1),' on jumper ',i
        CALL oft_abort("No matching edge for jumper points","tw_load_sensors",__FILE__)
      END IF
      cell=self%mesh%lec(self%mesh%kec(ABS(id)))
      DO k=self%kfh(cell),self%kfh(cell+1)-1
        vert=self%mesh%lc(self%lfh(2,k),cell)
        IF(vert==sensors%jumpers(i)%points(j))THEN
          hole_facs(ABS(self%lfh(1,k)))=hole_facs(ABS(self%lfh(1,k)))-SIGN(1,self%lfh(1,k))
        ELSE IF(vert==sensors%jumpers(i)%points(j+1))THEN
          hole_facs(ABS(self%lfh(1,k)))=hole_facs(ABS(self%lfh(1,k)))+SIGN(1,self%lfh(1,k))
        END IF
      END DO
    END DO
    !
    ALLOCATE(sensors%jumpers(i)%hole_facs(self%nholes))
    sensors%jumpers(i)%hole_facs=hole_facs
    ! Save orientation
    p1=sensors%jumpers(i)%points(1)
    p2=sensors%jumpers(i)%points(2)
    id=mesh_local_findedge(self%mesh,[p1,p2])
    cell=self%mesh%lec(self%mesh%kec(ABS(id)))
    DO vert=1,self%mesh%cell_np
      IF(self%mesh%lc(vert,cell)==sensors%jumpers(i)%points(1))EXIT
    END DO
    norm=-self%qbasis(:,vert,cell)
    DO vert=1,self%mesh%cell_np
      IF(self%mesh%lc(vert,cell)==sensors%jumpers(i)%points(2))EXIT
    END DO
    norm=norm+self%qbasis(:,vert,cell)
    ! Project into edge perpendicular direction and normalize
    f=1.d0/3.d0
    CALL self%mesh%norm(cell,f,cell_norm)
    ed_norm=cross_product(self%mesh%r(:,p2)-self%mesh%r(:,p1),cell_norm)
    norm=DOT_PRODUCT(norm,ed_norm)*ed_norm
    norm=norm/magnitude(norm)
    WRITE(io_unit,'(3ES24.15)')norm
    IF(oft_debug_print(1))WRITE(*,'(A,I8,3ES24.15)')oft_indent,i,norm
    ! Close loop, remove points and only use holes
    IF((sensors%jumpers(i)%points(1)==sensors%jumpers(i)%points(sensors%jumpers(i)%np)) &
      .AND.ALL(sensors%jumpers(i)%points(1:sensors%jumpers(i)%np)>0))THEN
      DEALLOCATE(sensors%jumpers(i)%points)
      sensors%jumpers(i)%np=0
    END IF
  END DO
  DEALLOCATE(hole_facs)
  CLOSE(io_unit)
END IF
CALL oft_decrease_indent
CONTAINS
!---------------------------------------------------------------------------------
!> Order jumper list into sequential chain
!---------------------------------------------------------------------------------
SUBROUTINE order_jumper_list(list_in,list_out,n)
INTEGER(4), INTENT(in) :: n !< Number of vertices in jumper definition
INTEGER(4), INTENT(in) :: list_in(n) !< Input vertex list
INTEGER(4), INTENT(out) :: list_out(n+1) !< Reordered vertex list
INTEGER(4) :: ii,jj,k,l,nlinks,ipt,eprev,ed,ed2,ptp2,ptp,candidate,candidate2,last_item(3),i0
INTEGER(4), ALLOCATABLE :: lloop_tmp(:),flag_list(:)
!---Find loop points and edges
ALLOCATE(lloop_tmp(n),flag_list(n))
flag_list=0
lloop_tmp=list_in
CALL sort_array(lloop_tmp, n)
!---Find suitable starting point
i0=1
DO jj=1,n
  k=0
  DO ii=self%mesh%kpe(lloop_tmp(jj)),self%mesh%kpe(lloop_tmp(jj)+1)-1
    ed = self%mesh%lpe(ii)
    ptp = SUM(self%mesh%le(:,ed))-lloop_tmp(jj)
    IF(search_array(ptp, lloop_tmp, n)==0)CYCLE
    k=k+1
    IF(k>1)EXIT
  END DO
  IF(k==1)THEN
    i0=jj
    EXIT
  END IF
END DO
!---Create oriented list of points
flag_list(i0)=1
ipt=lloop_tmp(i0)
k=1
list_out(k)=ipt
eprev=0
DO jj=1,n
  IF(jj==n-2)flag_list(i0)=0 ! Unflag first point for last link (needs to be 2 offset for link count check below)
  last_item=-1
  DO ii=self%mesh%kpe(ipt),self%mesh%kpe(ipt+1)-1
    ed = self%mesh%lpe(ii)
    IF(ed==eprev)CYCLE
    ptp = SUM(self%mesh%le(:,ed))-ipt
    candidate=search_array(ptp, lloop_tmp, n)
    IF(candidate==0)THEN
      CYCLE
    ELSE
      IF(flag_list(candidate)==1)CYCLE
      nlinks=0
      DO l=self%mesh%kpe(ptp),self%mesh%kpe(ptp+1)-1
        ed2 = self%mesh%lpe(l)
        ptp2 = SUM(self%mesh%le(:,ed2))-ptp
        candidate2=search_array(ptp2, lloop_tmp, n)
        IF(candidate2==0)CYCLE
        IF(flag_list(candidate2)==1)CYCLE
        nlinks=nlinks+1
      END DO
      last_item=[ptp,candidate,ed]
      IF(nlinks>1)CYCLE
      last_item=-1
    END IF
    flag_list(candidate)=1
    ipt=ptp
    IF(jj<n)THEN
      k=k+1
      list_out(k)=ipt
      eprev=ed
    END IF
    EXIT
  END DO
  IF(last_item(1)>0)THEN
    flag_list(last_item(2))=1
    ipt=last_item(1)
    IF(jj<n)THEN
      k=k+1
      list_out(k)=ipt
      eprev=last_item(3)
    END IF
  END IF
END DO
IF((flag_list(i0)==1).AND.(n>2))THEN
  list_out(n+1)=lloop_tmp(i0)
ELSE
  list_out(n+1)=0
  flag_list(i0)=1 ! Ok to not be a full loop
END IF
IF(ANY(flag_list==0))THEN
  WRITE(*,'(2A)')oft_indent,"No matching edges for following points:"
  CALL oft_increase_indent
  DO jj=1,n
    IF(flag_list(jj)==0)WRITE(*,'(A,I8,3ES14.5)')oft_indent,lloop_tmp(jj),self%mesh%r(:,lloop_tmp(jj))
  END DO
  CALL oft_decrease_indent
  CALL oft_abort('Error building jumper mesh, unmatched points exist', &
    'tw_load_sensors::order_jumper_list', __FILE__)
END IF
DEALLOCATE(lloop_tmp,flag_list)
END SUBROUTINE order_jumper_list
end subroutine tw_load_sensors
!------------------------------------------------------------------------------
!> Compute current [A] through a jumper given a ThinCurr solution vector
!------------------------------------------------------------------------------
FUNCTION jumper_compute_current_exp(self,tw_obj,vals) RESULT(jumper_current)
CLASS(jumper_sensor), INTENT(IN) :: self !< Jumper sensor object
CLASS(tw_type), INTENT(IN) :: tw_obj !< Thin-wall model object
REAL(r8), INTENT(IN) :: vals(:) !< Solution vector
REAL(r8) :: jumper_current !< Computed current through jumper  [A]
INTEGER(i4) :: k,k0
!---Change along each segment (segments are separated by 0 where a jumper crosses a junction)
jumper_current=0.d0
k0=1
DO k=1,self%np
  IF(self%points(k)==0)THEN
    k0=k+1
  ELSE IF(k==self%np.OR.self%points(MIN(k+1,self%np))==0)THEN
    jumper_current=jumper_current+point_value(self%points(k))-point_value(self%points(k0))
  END IF
END DO
DO k=1,tw_obj%nholes
  jumper_current=jumper_current+vals(tw_obj%np_active+k)*self%hole_facs(k)
END DO
jumper_current=jumper_current/mu0
CONTAINS
!> Solution value at a vertex (zero if inactive)
FUNCTION point_value(ipt) RESULT(val)
INTEGER(i4), INTENT(in) :: ipt !< Vertex index
REAL(r8) :: val
val=0.d0
IF(tw_obj%pmap(ipt)>0)val=vals(tw_obj%pmap(ipt))
END FUNCTION point_value
END FUNCTION jumper_compute_current_exp
!------------------------------------------------------------------------------
!> Load resistivity and sensor mask from "oft_in.xml" file
!------------------------------------------------------------------------------
subroutine tw_load_eta(self)
TYPE(tw_type), INTENT(inout) :: self !< Thin-wall model object
!---XML solver fields
integer(4) :: nshells,nreg_mesh,nread
TYPE(xml_node) :: sens_node,eta_group,eta_surf_group,eta_vol_group,thincurr_group,thickness_group
!---
INTEGER(4) :: i,j,io_unit,ierr,id,cell
REAL(8) :: location(2)
LOGICAL :: has_eta_surf,has_eta_vol,has_thickness
REAL(r8), POINTER :: eta_tmp(:),eta_vol_tmp(:),thickness_tmp(:)
LOGICAL, POINTER :: sens_mask_tmp(:)
nreg_mesh=MAXVAL(self%mesh%reg)
!--- Deallocate if already set
IF(ASSOCIATED(self%Eta_vol))DEALLOCATE(self%Eta_vol)
IF(ASSOCIATED(self%Thickness))DEALLOCATE(self%Thickness)
IF(ASSOCIATED(self%Eta_surf))DEALLOCATE(self%Eta_surf)
ALLOCATE(self%Eta_vol(nreg_mesh))
ALLOCATE(self%Thickness(nreg_mesh))
ALLOCATE(self%Eta_surf(nreg_mesh))
!---Default to negative values that are obviously wrong
self%Thickness=-1.d0
self%Eta_vol=-1.d0
self%Eta_surf=-1.d0
ALLOCATE(self%sens_mask(nreg_mesh))
self%sens_mask=.FALSE.
has_eta_surf=.FALSE.
has_eta_vol=.FALSE.
has_thickness=.FALSE.
IF(.NOT.self%xml%associated())THEN
  CALL oft_warn('No "thincurr" XML node specified. Ignore this warning if an XML node does not need to be specified.')
  RETURN
END IF
! Read surface resistivity values
CALL xml_get_element(self%xml,"eta",eta_group,ierr)
IF(ierr==0)THEN
  WRITE(*,*)
  WRITE(*,'(2A)')oft_indent,'Loading region surface resistivity:'
  CALL xml_read_content(eta_group,eta_tmp,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort('Error reading eta values','tw_load_eta',__FILE__)
  IF(SIZE(eta_tmp)/=SIZE(self%Eta_surf))CALL oft_abort('Eta size mismatch','tw_load_eta',__FILE__)
  IF(ANY(eta_tmp<=0.d0))CALL oft_abort('All "eta" values must be > 0','tw_load_eta',__FILE__)
  DO i=1,nreg_mesh
    WRITE(*,'(A,I4,ES12.4)')oft_indent,i,eta_tmp(i)
    self%Eta_surf(i)=eta_tmp(i)/mu0 ! Convert to magnetic units
  END DO
  DEALLOCATE(eta_tmp)
  has_eta_surf=.TRUE.
ELSE
  CALL xml_get_element(self%xml,"eta_surf",eta_surf_group,ierr)
  IF(ierr==0)THEN
    WRITE(*,*)
    WRITE(*,'(2A)')oft_indent,'Loading region surface resistivity:'
    CALL xml_read_content(eta_surf_group,eta_tmp,iostat=ierr)
    IF(ierr/=0)CALL oft_abort('Error reading eta values','tw_load_eta',__FILE__)
    IF(SIZE(eta_tmp)/=SIZE(self%Eta_surf))CALL oft_abort('Eta size mismatch','tw_load_eta',__FILE__)
    IF(ANY(eta_tmp<=0.d0))CALL oft_abort('All "eta" values must be > 0','tw_load_eta',__FILE__)
    DO i=1,nreg_mesh
      WRITE(*,'(A,I4,ES12.4)')oft_indent,i,eta_tmp(i)
      self%Eta_surf(i)=eta_tmp(i)/mu0 ! Convert to magnetic units
    END DO
    DEALLOCATE(eta_tmp)
    has_eta_surf=.TRUE.
  END IF
END IF
! Read volumetric resistivity values
CALL xml_get_element(self%xml,"eta_vol",eta_vol_group,ierr)
IF(ierr==0)THEN
  WRITE(*,*)
  WRITE(*,'(2A)')oft_indent,'Loading region volumetric resistivity:'
  CALL xml_read_content(eta_vol_group,eta_vol_tmp,iostat=ierr)
  IF(ierr/=0)CALL oft_abort('Error reading eta_vol values','tw_load_eta',__FILE__)
  IF(SIZE(eta_vol_tmp)/=SIZE(self%Eta_vol))CALL oft_abort('Eta_vol size mismatch','tw_load_eta',__FILE__)
  IF(ANY(eta_vol_tmp<=0.d0))CALL oft_abort('All "eta_vol" values must be > 0','tw_load_eta',__FILE__)
  DO i=1,nreg_mesh
    WRITE(*,'(A,I4,ES12.4)')oft_indent,i,eta_vol_tmp(i)
    self%Eta_vol(i)=eta_vol_tmp(i)/mu0 ! Convert to magnetic units
  END DO
  DEALLOCATE(eta_vol_tmp)
  has_eta_vol=.TRUE.
END IF
! Read thickness values
CALL xml_get_element(self%xml,"thickness",thickness_group,ierr)
IF(ierr==0)THEN
  WRITE(*,'(2A)')oft_indent,'Loading region thickness:'
  CALL xml_read_content(thickness_group,thickness_tmp,iostat=ierr)
  IF(ierr/=0)CALL oft_abort('Error reading thickness values','tw_load_eta',__FILE__)
  IF(SIZE(thickness_tmp)/=SIZE(self%Thickness))CALL oft_abort('Thickness size mismatch','tw_load_eta',__FILE__)
  IF(ANY(thickness_tmp<=0.d0))CALL oft_abort('All "thickness" values must be > 0','tw_load_eta',__FILE__)
  DO i=1,nreg_mesh
    WRITE(*,'(A,I4,ES12.4)')oft_indent,i,thickness_tmp(i)
    self%Thickness(i)=thickness_tmp(i)
  END DO
  DEALLOCATE(thickness_tmp)
  has_thickness=.TRUE.
END IF
IF(has_eta_surf)THEN
  IF(has_eta_vol.AND.has_thickness)THEN
    CALL oft_warn('"eta_surf" (or "eta"), "eta_vol," and "thickness" all provided from XML. Recalculating eta_surf from eta_vol')
    self%Eta_surf=self%Eta_vol/self%Thickness
  ELSE IF((.NOT.has_eta_vol).AND.has_thickness)THEN
    self%Eta_vol=self%Eta_surf*self%Thickness
  END IF
ELSE IF(has_eta_vol.AND.has_thickness)THEN
  self%Eta_surf=self%Eta_vol/self%Thickness
ELSE IF(has_eta_vol.AND.(.NOT.has_thickness))THEN
  CALL oft_warn('"eta_vol" specified without "thickness" nor "eta_surf" in XML. Please specify "eta_surf" or both "eta_vol" and "thickness".')
ELSE
  CALL oft_warn('Cannot gather or infer surface resisitivity from XML. Ignore this warning if resistivity is specified later or is not needed.')
END IF
! Read sensor mask
CALL xml_get_element(self%xml,"sens_mask",sens_node,ierr)
IF(ierr==0)THEN
  WRITE(*,'(2A)')oft_indent,'Loading sensor mask:'
  CALL xml_read_content(sens_node,sens_mask_tmp,iostat=ierr)
  IF(ierr/=0)CALL oft_xml_abort('Error reading sensor mask values','tw_load_eta',__FILE__)
  IF(SIZE(sens_mask_tmp)/=SIZE(self%sens_mask))CALL oft_abort('Sensor mask size mismatch','tw_load_eta',__FILE__)
  DO i=1,nreg_mesh
    WRITE(*,'(A,I4,L1)')oft_indent,i,sens_mask_tmp(i)
    self%sens_mask(i)=sens_mask_tmp(i)
  END DO
  DEALLOCATE(sens_mask_tmp)
  ! WRITE(*,*)'  Sens mask = ',self%sens_mask
END IF
end subroutine tw_load_eta
!---------------------------------------------------------------------------------
!> Save solution vector for thin-wall model for plotting in VisIt
!---------------------------------------------------------------------------------
SUBROUTINE tw_recon_curr_exp(self,pot,curr)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
REAL(8), INTENT(in) :: pot(:) !< Solution values [self%nelems]
REAL(8), INTENT(out) :: curr(:,:) !< Solution values [3,self%mesh%nelems]
INTEGER(4) :: i,j,k,jj,pt,ih,ihp,ihc
REAL(8) :: rcurr(3),ftmp(3),gop(3,3),area,norm(3)
DEBUG_STACK_PUSH
!---Avg to cells
ftmp=1.d0/3.d0
DO i=1,self%mesh%nc
  curr(:,i)=0.d0
  CALL self%mesh%jacobian(i,ftmp,gop,area)
  CALL self%mesh%norm(i,ftmp,norm)
  DO j=1,3
    pt=self%pmap(self%mesh%lc(j,i))
    IF(pt==0)CYCLE
    curr(:,i) = curr(:,i) + pot(pt)*cross_product(gop(:,j),norm)
  END DO
END DO
DO ih=1,self%nholes
DO ihp=1,self%hmesh(ih)%n
DO ihc=self%hmesh(ih)%kpc(ihp),self%hmesh(ih)%kpc(ihp+1)-1
i=ABS(self%hmesh(ih)%lpc(ihc))
  DO j=1,3
    IF(self%mesh%lc(j,i)==self%hmesh(ih)%lp(ihp))EXIT
  END DO
  CALL self%mesh%jacobian(i,ftmp,gop,area)
  CALL self%mesh%norm(i,ftmp,norm)
  curr(:,i) = curr(:,i) &
    + pot(self%np_active+ih)*cross_product(gop(:,j),norm)*SIGN(1,self%hmesh(ih)%lpc(ihc))
END DO
END DO
END DO
DEBUG_STACK_POP
END SUBROUTINE tw_recon_curr_exp
!---------------------------------------------------------------------------------
!> Save solution vector for thin-wall model for plotting in VisIt
!---------------------------------------------------------------------------------
SUBROUTINE tw_save_pfield_exp(self,a,tag)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
REAL(8), INTENT(in) :: a(:) !< Solution values [self%nelems]
CHARACTER(LEN=*), INTENT(in) :: tag !< Path to save vector in HDF5 plot files
INTEGER(4) :: i,j,k,jj,pt,ih,ihp,ihc,ic
REAL(8) :: rcurr(3),ftmp(3),gop(3,3),area,norm(3)
REAL(8), ALLOCATABLE, DIMENSION(:) :: thickness_cell
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: ptvec,cellvec,jvol_cellvec
DEBUG_STACK_PUSH
!---Avg to cells
ALLOCATE(cellvec(3,self%mesh%nc))
CALL tw_recon_curr_exp(self,a,cellvec)
CALL self%mesh%save_cell_vector(cellvec/mu0,self%xdmf,TRIM(tag)) ! Convert back to Amps
!---Automatically save volumetric current density when thickness is valid
IF(ASSOCIATED(self%Thickness))THEN
  IF(ALL(self%Thickness>0.d0))THEN
    ALLOCATE(thickness_cell(self%mesh%nc))
    ALLOCATE(jvol_cellvec(3,self%mesh%nc))
    jvol_cellvec=cellvec
    DO ic=1,self%mesh%nc
      thickness_cell(ic)=self%Thickness(self%mesh%reg(ic))
      jvol_cellvec(:,ic)=jvol_cellvec(:,ic)/(mu0*thickness_cell(ic))
    END DO
    CALL self%mesh%save_cell_scalar(thickness_cell,self%xdmf,'thickness')
    CALL self%mesh%save_cell_vector(jvol_cellvec,self%xdmf,'J_vol')
    DEALLOCATE(thickness_cell,jvol_cellvec)
  END IF
END IF
!---Avg to points
ALLOCATE(ptvec(3,self%mesh%np))
DO i=1,self%mesh%np
  ptvec(:,i)=0.d0
  DO j=self%mesh%kpc(i),self%mesh%kpc(i+1)-1
    ptvec(:,i) = ptvec(:,i) + cellvec(:,self%mesh%lpc(j))*self%mesh%ca(self%mesh%lpc(j))/3.d0
  END DO
  ptvec(:,i) = ptvec(:,i)/self%mesh%va(i)
END DO
CALL self%mesh%save_vertex_vector(ptvec/mu0,self%xdmf,TRIM(tag)//'_v') ! Convert back to Amps
!
ptvec=0.d0
DO i=1,self%mesh%np
  IF(self%pmap(i)>0)ptvec(1,i)=a(self%pmap(i))
END DO
CALL self%mesh%save_vertex_scalar(ptvec(1,:),self%xdmf,TRIM(tag)//'_p')
DEALLOCATE(ptvec,cellvec)
DEBUG_STACK_POP
END SUBROUTINE tw_save_pfield_exp
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE tw_save_hole_debug(self)
TYPE(tw_type), INTENT(in) :: self
INTEGER(4) :: i,j,k,jj,pt,ih,ihp,ihc
REAL(8) :: rcurr(3),ftmp(3),gop(3,3),area,norm(3)
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: cellvec
CHARACTER(LEN=3) :: pltnum
!---Avg to cells
ALLOCATE(cellvec(3,self%mesh%nc))
cellvec=0.d0
DO ih=1,self%nholes
!cellvec=0.d0
DO ihp=1,self%hmesh(ih)%n
DO ihc=self%hmesh(ih)%kpc(ihp),self%hmesh(ih)%kpc(ihp+1)-1
i=ABS(self%hmesh(ih)%lpc(ihc))
  DO j=1,3
    IF(self%mesh%lc(j,i)==self%hmesh(ih)%lp(ihp))EXIT
  END DO
  CALL self%mesh%jacobian(i,ftmp,gop,area)
  CALL self%mesh%norm(i,ftmp,norm)
  cellvec(:,i) = cellvec(:,i) &
    + cross_product(gop(:,j),norm)*SIGN(1,self%hmesh(ih)%lpc(ihc))
END DO
END DO
!WRITE(*,*)'Saving hole',ih
!WRITE(pltnum,'(I3.3)')ih
!CALL self%mesh%save_cell_vector(cellvec,"hole_"//pltnum)
END DO
CALL self%mesh%save_cell_vector(cellvec,self%xdmf,"hole_vec")
DEALLOCATE(cellvec)
END SUBROUTINE tw_save_hole_debug
!------------------------------------------------------------------------------
!> Save Thin-wall solution vector to a restart file
!------------------------------------------------------------------------------
subroutine tw_rst_save(self,u,filename,path,append)
class(tw_type), intent(in) :: self !< Thin-wall model object
class(oft_vector), pointer, intent(inout) :: u !< Solution to save
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to store solution vector in file
logical, optional, intent(in) :: append !< Append to file? (optional)
integer(8) :: i
integer(8), allocatable, dimension(:) :: global_le
real(8), allocatable :: pcoil_vals(:)
class(oft_native_vector), pointer :: outvec
type(hdf5_rst) :: rst_info
logical :: do_append
DEBUG_STACK_PUSH
do_append=.FALSE.
IF(PRESENT(append))do_append=append
IF(.NOT.do_append)THEN
  IF(oft_env%head_proc)THEN
    CALL hdf5_create_file(filename)
    CALL hdf5_write(tw_idx_ver,filename,tw_idx_path)
  END IF
END IF
!---
IF(oft_env%head_proc.AND.oft_env%pm)WRITE(*,'(6A)')oft_indent,'Writing "',TRIM(path), &
  '" to restart file "',TRIM(filename),'"'
IF(.NOT.native_vector_cast(outvec,u))CALL oft_abort('Failed to cast "source".', &
  'tw_rst_save',__FILE__)
!---
ALLOCATE(global_le(u%n))
DO i=1,u%n
  global_le(i) = i
END DO
CALL native_vector_slice_push(outvec,global_le,rst_info)
CALL hdf5_write(rst_info,filename,path)
IF(self%n_vcoils>0)THEN
  allocate(pcoil_vals(self%n_vcoils))
  pcoil_vals=outvec%v(self%nelems-self%n_vcoils+1:self%nelems)
  CALL hdf5_write(pcoil_vals,filename,path//"_Vcoils")
  deallocate(pcoil_vals)
END IF
CALL hdf5_rst_destroy(rst_info)
DEALLOCATE(global_le)
DEBUG_STACK_POP
end subroutine tw_rst_save
!------------------------------------------------------------------------------
!> Load Thin-wall solution vector from a restart file
!------------------------------------------------------------------------------
subroutine tw_rst_load(u,filename,path)
class(oft_vector), target, intent(inout) :: u !< Solution to load
character(LEN=*), intent(in) :: filename !< Name of restart file
character(LEN=*), intent(in) :: path !< Path to load solution vector in file
integer(8) :: i
integer(8), allocatable, dimension(:) :: global_le
class(oft_native_vector), pointer :: invec
type(hdf5_rst) :: rst_info
DEBUG_STACK_PUSH
IF(oft_env%head_proc.AND.oft_env%pm)WRITE(*,*)'Reading "',TRIM(path), &
  '" from restart file "',TRIM(filename),'"'
IF(.NOT.native_vector_cast(invec,u))CALL oft_abort('Failed to cast "source".', &
  'tw_rst_load',__FILE__)
!
ALLOCATE(global_le(u%n))
DO i=1,u%n
  global_le(i) = i
END DO
!---
CALL native_vector_slice_push(invec,global_le,rst_info,alloc_only=.TRUE.)
CALL hdf5_read(rst_info,filename,path)
CALL native_vector_slice_pop(invec,global_le,rst_info)
!---
CALL hdf5_rst_destroy(rst_info)
DEALLOCATE(global_le)
DEBUG_STACK_POP
end subroutine tw_rst_load
!---------------------------------------------------------------------------------
!> Set up an explicit hole from its (cell, vertex, sign) triples (value sign at that corner of the cell)
!---------------------------------------------------------------------------------
SUBROUTINE tw_setup_jhole(bmesh,pairs,hmesh)
CLASS(oft_bmesh), INTENT(in) :: bmesh !< Surface mesh
TYPE(oft_1d_int), INTENT(in) :: pairs !< (cell, vertex, sign) triples [3*n]
TYPE(hole_mesh), INTENT(inout) :: hmesh !< Hole definition
INTEGER(4) :: i,j,n
INTEGER(4), ALLOCATABLE :: pts(:),cnt(:)
n=pairs%n
ALLOCATE(pts(n))
DO i=1,n
  pts(i)=pairs%v(3*i-1)
END DO
CALL sort_array(pts,n)
hmesh%n=1
DO i=2,n
  IF(pts(i)/=pts(hmesh%n))THEN
    hmesh%n=hmesh%n+1
    pts(hmesh%n)=pts(i)
  END IF
END DO
ALLOCATE(hmesh%lp(hmesh%n),hmesh%lp_sorted(hmesh%n),hmesh%sort_ind(hmesh%n))
hmesh%lp=pts(1:hmesh%n)
hmesh%lp_sorted=hmesh%lp
hmesh%sort_ind=[(i,i=1,hmesh%n)]
ALLOCATE(cnt(hmesh%n+1),hmesh%kpc(hmesh%n+1),hmesh%lpc(n))
cnt=0
DO i=1,n
  j=search_array(pairs%v(3*i-1),hmesh%lp,hmesh%n)
  cnt(j)=cnt(j)+1
END DO
hmesh%kpc(1)=1
DO i=1,hmesh%n
  hmesh%kpc(i+1)=hmesh%kpc(i)+cnt(i)
END DO
cnt=0
DO i=1,n
  j=search_array(pairs%v(3*i-1),hmesh%lp,hmesh%n)
  hmesh%lpc(hmesh%kpc(j)+cnt(j))=SIGN(pairs%v(3*i-2),pairs%v(3*i))
  cnt(j)=cnt(j)+1
END DO
hmesh%ptcc=0.d0
DO i=1,hmesh%n
  hmesh%ptcc=hmesh%ptcc+bmesh%r(:,hmesh%lp(i))
END DO
hmesh%ptcc=hmesh%ptcc/REAL(hmesh%n,8)
DEALLOCATE(pts,cnt)
END SUBROUTINE tw_setup_jhole
!---------------------------------------------------------------------------------
!> Check the junction map against the model and set the reduced element count
!---------------------------------------------------------------------------------
SUBROUTINE tw_setup_junctions(self)
TYPE(tw_type), INTENT(inout) :: self !< Thin-wall model object
INTEGER(4) :: i,j,s
IF(SIZE(self%junc%p_kr)-1/=self%np_active+self%nholes)THEN
  WRITE(*,*)SIZE(self%junc%p_kr)-1,self%np_active,self%nholes
  CALL oft_abort("Junction map does not match the active vertices and holes","tw_setup_junctions",__FILE__)
END IF
DO i=1,SIZE(self%junc%table,DIM=2)
  DO j=1,3
    IF(self%mesh%lc(j,self%junc%table(3,i))==self%junc%table(1,i))EXIT
  END DO
  s=-1
  IF(self%mesh%lc(MOD(j,3)+1,self%junc%table(3,i))==self%junc%table(2,i))s=1
  IF(s/=self%junc%table(4,i))CALL oft_abort("Junction signs do not match the mesh orientation", &
    "tw_setup_junctions",__FILE__)
END DO
CALL tw_junction_vcoil_rows(self)
WRITE(*,'(2A,I8,A,I8)')oft_indent,'Junctions: unknowns reduced from ',self%junc%nelems_exp,' to ',self%nelems
END SUBROUTINE tw_setup_junctions
!---------------------------------------------------------------------------------
!> Set the V-coil rows of the junction map (identity) for the current V-coils, so that P covers all expanded
!! elements, and the reduced element count. Called at setup and whenever the coil types change.
!---------------------------------------------------------------------------------
SUBROUTINE tw_junction_vcoil_rows(self)
TYPE(tw_type), INTENT(inout) :: self !< Thin-wall model object
INTEGER(4) :: i,ns,nnz
INTEGER(i4), POINTER, DIMENSION(:) :: kr_tmp,lc_tmp
REAL(r8), POINTER, DIMENSION(:) :: val_tmp
!---Keep the vertex and hole rows, drop any previous V-coil rows
ns=self%np_active+self%nholes
nnz=self%junc%p_kr(ns+1)-1
ALLOCATE(kr_tmp(ns+self%n_vcoils+1),lc_tmp(nnz+self%n_vcoils),val_tmp(nnz+self%n_vcoils))
kr_tmp(1:ns+1)=self%junc%p_kr(1:ns+1)
lc_tmp(1:nnz)=self%junc%p_lc(1:nnz)
val_tmp(1:nnz)=self%junc%p_val(1:nnz)
DO i=1,self%n_vcoils
  kr_tmp(ns+i+1)=nnz+i+1
  lc_tmp(nnz+i)=self%junc%nred+i
  val_tmp(nnz+i)=1.d0
END DO
DEALLOCATE(self%junc%p_kr,self%junc%p_lc,self%junc%p_val)
self%junc%p_kr=>kr_tmp
self%junc%p_lc=>lc_tmp
self%junc%p_val=>val_tmp
self%nelems=self%junc%nred+self%n_vcoils
END SUBROUTINE tw_junction_vcoil_rows
!---------------------------------------------------------------------------------
!> Expand a reduced vector, x_exp = P x
!---------------------------------------------------------------------------------
SUBROUTINE tw_junction_expand(self,x,x_exp)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
REAL(8), INTENT(in) :: x(:) !< Reduced vector [nelems]
REAL(8), INTENT(out) :: x_exp(:) !< Expanded vector [nelems_exp]
INTEGER(4) :: e,k
DO e=1,self%junc%nelems_exp
  x_exp(e)=0.d0
  DO k=self%junc%p_kr(e),self%junc%p_kr(e+1)-1
    x_exp(e)=x_exp(e)+self%junc%p_val(k)*x(self%junc%p_lc(k))
  END DO
END DO
END SUBROUTINE tw_junction_expand
!---------------------------------------------------------------------------------
!> Junction map P in CSR form, or the identity for a model without junctions
!---------------------------------------------------------------------------------
SUBROUTINE tw_junction_map(self,kr,lc,val)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
INTEGER(4), ALLOCATABLE, INTENT(out) :: kr(:) !< Row pointer [nelems_exp+1]
INTEGER(4), ALLOCATABLE, INTENT(out) :: lc(:) !< Reduced element of each entry
REAL(8), ALLOCATABLE, INTENT(out) :: val(:) !< Coefficients
INTEGER(4) :: e,n
n=tw_nelems_exp(self)
IF(ASSOCIATED(self%junc))THEN
  ALLOCATE(kr(n+1),lc(self%junc%p_kr(n+1)-1),val(self%junc%p_kr(n+1)-1))
  kr=self%junc%p_kr(1:n+1)
  lc=self%junc%p_lc(1:SIZE(lc))
  val=self%junc%p_val(1:SIZE(val))
ELSE
  ALLOCATE(kr(n+1),lc(n),val(n))
  kr=[(e,e=1,n+1)]
  lc=[(e,e=1,n)]
  val=1.d0
END IF
END SUBROUTINE tw_junction_map
!---------------------------------------------------------------------------------
!> Add the contributions of a list of expanded element pairs to a reduced inductance matrix
!!
!! For each pair (`e_c`,`e_r`) with value `v`, adds \f$ P_c(e_c,a) P_r(e_r,b) v \f$ to `A(a,b)` for each entry
!! of row `e_c` of the column model map \f$ P_c \f$ and row `e_r` of the row model map \f$ P_r \f$. If `lsym`
!! is set the maps are the same, only expanded pairs with `e_c >= e_r` are given (the mirror pair is implied)
!! and the lower triangle of `A` is accumulated. The diagonal of `A` then receives the mirror pair's
!! contribution too, while reduced pairs from a diagonal expanded entry are only added once.
!!
!! @note Updates are atomic, so this may be called from within an OpenMP parallel region
!---------------------------------------------------------------------------------
SUBROUTINE tw_junction_scatter_add(A,ckr,clc,cval,rkr,rlc,rval,n,epairs,v,lsym)
REAL(8), INTENT(inout) :: A(:,:) !< Reduced matrix [nelems_c,nelems_r]
INTEGER(4), INTENT(in) :: ckr(:) !< Column model map row pointer
INTEGER(4), INTENT(in) :: clc(:) !< Column model map column indices
REAL(8), INTENT(in) :: cval(:) !< Column model map values
INTEGER(4), INTENT(in) :: rkr(:) !< Row model map row pointer
INTEGER(4), INTENT(in) :: rlc(:) !< Row model map column indices
REAL(8), INTENT(in) :: rval(:) !< Row model map values
INTEGER(4), INTENT(in) :: n !< Number of expanded element pairs
INTEGER(4), INTENT(in) :: epairs(:,:) !< Expanded elements (`e_c`,`e_r`) of each pair [2,n]
REAL(8), INTENT(in) :: v(:) !< Contribution of each pair over expanded elements [n]
LOGICAL, INTENT(in) :: lsym !< Accumulate the lower triangle of a symmetric matrix?
INTEGER(4) :: i,k,l,ec,er,ia,ib
REAL(8) :: w
DO i=1,n
  ec=epairs(1,i)
  er=epairs(2,i)
  DO k=ckr(ec),ckr(ec+1)-1
    DO l=rkr(er),rkr(er+1)-1
      w=cval(k)*rval(l)*v(i)
      ia=clc(k)
      ib=rlc(l)
      IF(lsym)THEN
        IF(ec==er)THEN
          IF(ia<ib)CYCLE
        ELSE IF(ia==ib)THEN
          w=2.d0*w
        END IF
        IF(ia<ib)THEN
          ia=rlc(l)
          ib=clc(k)
        END IF
      END IF
      !$omp atomic
      A(ia,ib) = A(ia,ib) + w
    END DO
  END DO
END DO
END SUBROUTINE tw_junction_scatter_add
!---------------------------------------------------------------------------------
!> Hash of the junction map P, used to key cache files (0 for a model without junctions)
!!
!! Only the vertex and hole rows are hashed. The V-coil rows are the identity, and the coils are keyed
!! separately where a matrix depends on them, so caches stay valid when only the V-coils change.
!---------------------------------------------------------------------------------
FUNCTION tw_junction_hash(self) RESULT(hash)
TYPE(tw_type), TARGET, INTENT(in) :: self !< Thin-wall model object
INTEGER(4) :: hash !< Hash of the map
INTEGER(4) :: ns,nnz
IF(.NOT.ASSOCIATED(self%junc))THEN
  hash=0
  RETURN
END IF
ns=self%np_active+self%nholes
nnz=self%junc%p_kr(ns+1)-1
hash=oft_simple_hash(C_LOC(self%junc%p_kr),INT(4*(ns+1),8))
hash=oft_simple_hash(C_LOC(self%junc%p_lc),INT(4*nnz,8),starting_value=hash)
hash=oft_simple_hash(C_LOC(self%junc%p_val),INT(8*nnz,8),starting_value=hash)
END FUNCTION tw_junction_hash
!---------------------------------------------------------------------------------
!> Transpose of the junction map, P^T in CSR form
!---------------------------------------------------------------------------------
SUBROUTINE tw_junction_transpose(self,pt_kr,pt_lc,pt_val)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
INTEGER(4), ALLOCATABLE, INTENT(out) :: pt_kr(:) !< Row pointer [nelems+1]
INTEGER(4), ALLOCATABLE, INTENT(out) :: pt_lc(:) !< Expanded element of each entry
REAL(8), ALLOCATABLE, INTENT(out) :: pt_val(:) !< Coefficients
INTEGER(4) :: a,e,k
INTEGER(4), ALLOCATABLE :: cnt(:)
ALLOCATE(pt_kr(self%nelems+1),cnt(self%nelems))
cnt=0
DO k=1,self%junc%p_kr(self%junc%nelems_exp+1)-1
  cnt(self%junc%p_lc(k))=cnt(self%junc%p_lc(k))+1
END DO
pt_kr(1)=1
DO a=1,self%nelems
  pt_kr(a+1)=pt_kr(a)+cnt(a)
END DO
ALLOCATE(pt_lc(pt_kr(self%nelems+1)-1),pt_val(pt_kr(self%nelems+1)-1))
cnt=0
DO e=1,self%junc%nelems_exp
  DO k=self%junc%p_kr(e),self%junc%p_kr(e+1)-1
    a=self%junc%p_lc(k)
    pt_lc(pt_kr(a)+cnt(a))=e
    pt_val(pt_kr(a)+cnt(a))=self%junc%p_val(k)
    cnt(a)=cnt(a)+1
  END DO
END DO
DEALLOCATE(cnt)
END SUBROUTINE tw_junction_transpose
!---------------------------------------------------------------------------------
!> Reduce the resistance matrix in place, R = P^T R_exp P
!---------------------------------------------------------------------------------
SUBROUTINE tw_junction_reduce_sparse(self)
TYPE(tw_type), INTENT(inout) :: self !< Thin-wall model object
INTEGER(4) :: a,e,i,j,k,l,m,nr,nnz
INTEGER(4), ALLOCATABLE :: pt_kr(:),pt_lc(:),kr(:),lc(:),mark(:),ind(:),tmp(:)
REAL(8), ALLOCATABLE :: pt_val(:),mvals(:),acc(:)
TYPE(oft_native_matrix), POINTER :: Rexp,Rred
Rexp=>self%Rmat
nr=self%nelems
CALL tw_junction_transpose(self,pt_kr,pt_lc,pt_val)
ALLOCATE(tmp(nr))
!---Sparsity of P^T R P
ALLOCATE(mark(nr),kr(nr+1))
mark=0
kr(1)=1
DO a=1,nr
  m=0
  DO i=pt_kr(a),pt_kr(a+1)-1
    e=pt_lc(i)
    DO j=Rexp%kr(e),Rexp%kr(e+1)-1
      DO k=self%junc%p_kr(Rexp%lc(j)),self%junc%p_kr(Rexp%lc(j)+1)-1
        IF(mark(self%junc%p_lc(k))/=a)THEN
          mark(self%junc%p_lc(k))=a
          m=m+1
        END IF
      END DO
    END DO
  END DO
  kr(a+1)=kr(a)+m
END DO
nnz=kr(nr+1)-1
ALLOCATE(lc(nnz),mvals(nnz),acc(nr),ind(nr))
mark=0
acc=0.d0
DO a=1,nr
  m=0
  DO i=pt_kr(a),pt_kr(a+1)-1
    e=pt_lc(i)
    DO j=Rexp%kr(e),Rexp%kr(e+1)-1
      DO k=self%junc%p_kr(Rexp%lc(j)),self%junc%p_kr(Rexp%lc(j)+1)-1
        IF(mark(self%junc%p_lc(k))/=a)THEN
          mark(self%junc%p_lc(k))=a
          m=m+1
          tmp(m)=self%junc%p_lc(k)
        END IF
        acc(self%junc%p_lc(k))=acc(self%junc%p_lc(k))+pt_val(i)*Rexp%M(j)*self%junc%p_val(k)
      END DO
    END DO
  END DO
  ind(1:m)=[(l,l=1,m)]
  CALL sort_array(tmp(1:m),ind(1:m),m)
  DO l=1,m
    lc(kr(a)+l-1)=tmp(l)
    mvals(kr(a)+l-1)=acc(tmp(l))
    acc(tmp(l))=0.d0
  END DO
END DO
ALLOCATE(Rred)
Rred%nr=nr; Rred%nrg=nr
Rred%nc=nr; Rred%ncg=nr
Rred%nnz=nnz
Rred%kr=>NULL(); Rred%lc=>NULL(); Rred%M=>NULL()
ALLOCATE(Rred%kr(nr+1),Rred%lc(nnz),Rred%M(nnz))
Rred%kr=kr; Rred%lc=lc; Rred%M=mvals
self%junc%Rmat_exp=>Rexp
self%Rmat=>Rred
DEALLOCATE(pt_kr,pt_lc,pt_val,tmp,mark,kr,lc,mvals,acc,ind)
END SUBROUTINE tw_junction_reduce_sparse
!---------------------------------------------------------------------------------
!> Compute the resistance matrix (reduced to the junction unknowns if present)
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_Rmat(tw_obj,keep_closures)
TYPE(tw_type), INTENT(inout) :: tw_obj !< Thin-wall model object
LOGICAL, INTENT(in) :: keep_closures !< Keep diagonal entries (1) for closure elements
CALL tw_compute_Rmat_exp(tw_obj,keep_closures)
IF(ASSOCIATED(tw_obj%junc))CALL tw_junction_reduce_sparse(tw_obj)
END SUBROUTINE tw_compute_Rmat
!---------------------------------------------------------------------------------
!> Reconstruct the surface current on each cell (expanding junction unknowns if present)
!---------------------------------------------------------------------------------
SUBROUTINE tw_recon_curr(self,pot,curr)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
REAL(8), INTENT(in) :: pot(:) !< Solution values [self%nelems]
REAL(8), INTENT(out) :: curr(:,:) !< Solution values [3,self%mesh%nelems]
REAL(8), ALLOCATABLE :: pot_exp(:)
IF(.NOT.ASSOCIATED(self%junc))THEN
  CALL tw_recon_curr_exp(self,pot,curr)
  RETURN
END IF
ALLOCATE(pot_exp(tw_nelems_exp(self)))
CALL tw_junction_expand(self,pot,pot_exp)
CALL tw_recon_curr_exp(self,pot_exp,curr)
DEALLOCATE(pot_exp)
END SUBROUTINE tw_recon_curr
!---------------------------------------------------------------------------------
!> Reduce the rows of a matrix over expanded elements, A = P^T A_exp (A_exp is replaced)
!---------------------------------------------------------------------------------
SUBROUTINE tw_junction_reduce_rows(self,A)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
REAL(8), CONTIGUOUS, POINTER, INTENT(inout) :: A(:,:) !< Matrix [nelems_exp,m] -> [nelems,m]
INTEGER(4) :: e,k
REAL(8), CONTIGUOUS, POINTER :: Ared(:,:)
ALLOCATE(Ared(self%nelems,SIZE(A,DIM=2)))
Ared=0.d0
DO e=1,self%junc%nelems_exp
  DO k=self%junc%p_kr(e),self%junc%p_kr(e+1)-1
    Ared(self%junc%p_lc(k),:)=Ared(self%junc%p_lc(k),:)+self%junc%p_val(k)*A(e,:)
  END DO
END DO
DEALLOCATE(A)
A=>Ared
END SUBROUTINE tw_junction_reduce_rows
!---------------------------------------------------------------------------------
!> Reduce the columns of a matrix over expanded elements, A = A_exp P (A_exp is replaced)
!---------------------------------------------------------------------------------
SUBROUTINE tw_junction_reduce_cols(self,A)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
REAL(8), CONTIGUOUS, POINTER, INTENT(inout) :: A(:,:) !< Matrix [m,nelems_exp] -> [m,nelems]
INTEGER(4) :: e,k
REAL(8), CONTIGUOUS, POINTER :: Ared(:,:)
ALLOCATE(Ared(SIZE(A,DIM=1),self%nelems))
Ared=0.d0
DO e=1,self%junc%nelems_exp
  DO k=self%junc%p_kr(e),self%junc%p_kr(e+1)-1
    Ared(:,self%junc%p_lc(k))=Ared(:,self%junc%p_lc(k))+self%junc%p_val(k)*A(:,e)
  END DO
END DO
DEALLOCATE(A)
A=>Ared
END SUBROUTINE tw_junction_reduce_cols
!---------------------------------------------------------------------------------
!> Compute the element to driver coupling (reduced to the junction unknowns if present)
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_Ael2dr(tw_obj,save_file)
TYPE(tw_type), INTENT(inout) :: tw_obj !< Thin-wall model object
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: save_file !< Cache file
CALL tw_compute_Ael2dr_exp(tw_obj,save_file)
IF(ASSOCIATED(tw_obj%junc).AND.ASSOCIATED(tw_obj%Ael2dr))CALL tw_junction_reduce_rows(tw_obj,tw_obj%Ael2dr)
END SUBROUTINE tw_compute_Ael2dr
!---------------------------------------------------------------------------------
!> Compute the element to sensor coupling (reduced to the junction unknowns if present)
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_mutuals(tw_obj,nsensors,sensors,save_file)
TYPE(tw_type), INTENT(inout) :: tw_obj !< Thin-wall model object
INTEGER(4), INTENT(in) :: nsensors !< Number of sensors
TYPE(floop_sensor), POINTER, INTENT(in) :: sensors(:) !< Sensors
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: save_file !< Cache file
CALL tw_compute_mutuals_exp(tw_obj,nsensors,sensors,save_file)
IF(.NOT.ASSOCIATED(tw_obj%junc))RETURN
CALL tw_junction_reduce_cols(tw_obj,tw_obj%Ael2sen)
END SUBROUTINE tw_compute_mutuals
!---------------------------------------------------------------------------------
!> Current through a jumper (expanding junction unknowns if present)
!---------------------------------------------------------------------------------
FUNCTION jumper_compute_current(self,tw_obj,vals) RESULT(jumper_current)
CLASS(jumper_sensor), INTENT(IN) :: self !< Jumper sensor object
CLASS(tw_type), INTENT(IN) :: tw_obj !< Thin-wall model object
REAL(r8), INTENT(IN) :: vals(:) !< Solution vector
REAL(r8) :: jumper_current !< Computed current through jumper  [A]
REAL(r8), ALLOCATABLE :: vals_exp(:)
SELECT TYPE(tw_obj)
TYPE IS(tw_type)
  IF(ASSOCIATED(tw_obj%junc))THEN
    ALLOCATE(vals_exp(tw_nelems_exp(tw_obj)))
    CALL tw_junction_expand(tw_obj,vals,vals_exp)
    jumper_current=jumper_compute_current_exp(self,tw_obj,vals_exp)
    DEALLOCATE(vals_exp)
    RETURN
  END IF
END SELECT
jumper_current=jumper_compute_current_exp(self,tw_obj,vals)
END FUNCTION jumper_compute_current
!---------------------------------------------------------------------------------
!> Value of hole element i in a solution vector (expanding junction unknowns if present)
!---------------------------------------------------------------------------------
FUNCTION tw_hole_value(self,vals,i) RESULT(val)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
REAL(r8), INTENT(in) :: vals(:) !< Solution vector [nelems]
INTEGER(4), INTENT(in) :: i !< Hole index
REAL(r8) :: val
INTEGER(4) :: k
IF(.NOT.ASSOCIATED(self%junc))THEN
  val=vals(self%np_active+i)
  RETURN
END IF
val=0.d0
DO k=self%junc%p_kr(self%np_active+i),self%junc%p_kr(self%np_active+i+1)-1
  val=val+self%junc%p_val(k)*vals(self%junc%p_lc(k))
END DO
END FUNCTION tw_hole_value
!---------------------------------------------------------------------------------
!> Save the current and potential to plot files (expanding junction unknowns if present)
!---------------------------------------------------------------------------------
SUBROUTINE tw_save_pfield(self,a,tag)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
REAL(8), INTENT(in) :: a(:) !< Solution values [self%nelems]
CHARACTER(LEN=*), INTENT(in) :: tag !< Path to save vector in HDF5 plot files
REAL(8), ALLOCATABLE :: a_exp(:)
IF(.NOT.ASSOCIATED(self%junc))THEN
  CALL tw_save_pfield_exp(self,a,tag)
  RETURN
END IF
ALLOCATE(a_exp(tw_nelems_exp(self)))
CALL tw_junction_expand(self,a,a_exp)
CALL tw_save_pfield_exp(self,a_exp,tag)
DEALLOCATE(a_exp)
END SUBROUTINE tw_save_pfield
!---------------------------------------------------------------------------------
!> Reduce an expanded vector, x = P^T x_exp
!---------------------------------------------------------------------------------
SUBROUTINE tw_junction_reduce_vec(self,x_exp,x)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
REAL(8), INTENT(in) :: x_exp(:) !< Expanded vector [nelems_exp]
REAL(8), INTENT(out) :: x(:) !< Reduced vector [nelems]
INTEGER(4) :: e,k
x=0.d0
DO e=1,self%junc%nelems_exp
  DO k=self%junc%p_kr(e),self%junc%p_kr(e+1)-1
    x(self%junc%p_lc(k))=x(self%junc%p_lc(k))+self%junc%p_val(k)*x_exp(e)
  END DO
END DO
END SUBROUTINE tw_junction_reduce_vec
!---------------------------------------------------------------------------------
!> Expand a reduced complex vector, x_exp = P x
!---------------------------------------------------------------------------------
SUBROUTINE tw_junction_expand_c(self,x,x_exp)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
COMPLEX(8), INTENT(in) :: x(:) !< Reduced vector [nelems]
COMPLEX(8), INTENT(out) :: x_exp(:) !< Expanded vector [nelems_exp]
INTEGER(4) :: e,k
DO e=1,self%junc%nelems_exp
  x_exp(e)=(0.d0,0.d0)
  DO k=self%junc%p_kr(e),self%junc%p_kr(e+1)-1
    x_exp(e)=x_exp(e)+self%junc%p_val(k)*x(self%junc%p_lc(k))
  END DO
END DO
END SUBROUTINE tw_junction_expand_c
!---------------------------------------------------------------------------------
!> Reduce an expanded complex vector, x = P^T x_exp
!---------------------------------------------------------------------------------
SUBROUTINE tw_junction_reduce_vec_c(self,x_exp,x)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
COMPLEX(8), INTENT(in) :: x_exp(:) !< Expanded vector [nelems_exp]
COMPLEX(8), INTENT(out) :: x(:) !< Reduced vector [nelems]
INTEGER(4) :: e,k
x=(0.d0,0.d0)
DO e=1,self%junc%nelems_exp
  DO k=self%junc%p_kr(e),self%junc%p_kr(e+1)-1
    x(self%junc%p_lc(k))=x(self%junc%p_lc(k))+self%junc%p_val(k)*x_exp(e)
  END DO
END DO
END SUBROUTINE tw_junction_reduce_vec_c
!---------------------------------------------------------------------------------
!> Diagonal of P^T D P for a diagonal D over expanded elements (cross terms dropped)
!---------------------------------------------------------------------------------
SUBROUTINE tw_junction_reduce_diag(self,d_exp,d)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
REAL(8), INTENT(in) :: d_exp(:) !< Expanded diagonal [nelems_exp]
REAL(8), INTENT(out) :: d(:) !< Reduced diagonal [nelems]
INTEGER(4) :: e,k
d=0.d0
DO e=1,self%junc%nelems_exp
  DO k=self%junc%p_kr(e),self%junc%p_kr(e+1)-1
    d(self%junc%p_lc(k))=d(self%junc%p_lc(k))+self%junc%p_val(k)**2*d_exp(e)
  END DO
END DO
END SUBROUTINE tw_junction_reduce_diag
!---------------------------------------------------------------------------------
!> Apply the mutual inductance between two models, b = M a (reduced to the junction unknowns if present)
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_Lmat_MF(row_obj,col_obj,nrhs,a,b)
TYPE(tw_type), INTENT(in) :: row_obj !< Thin-wall model object for rows
TYPE(tw_type), INTENT(in) :: col_obj !< Thin-wall model object for columns
INTEGER(4), INTENT(in) :: nrhs !< Number of vectors
REAL(8), CONTIGUOUS, POINTER, INTENT(in) :: a(:,:) !< Currents on row_obj [nelems,nrhs]
REAL(8), CONTIGUOUS, POINTER, INTENT(inout) :: b(:,:) !< Fluxes on col_obj [nelems,nrhs]
INTEGER(4) :: i
REAL(8), CONTIGUOUS, POINTER :: aexp(:,:),bexp(:,:)
IF(.NOT.(ASSOCIATED(row_obj%junc).OR.ASSOCIATED(col_obj%junc)))THEN
  CALL tw_compute_Lmat_MF_exp(row_obj,col_obj,nrhs,a,b)
  RETURN
END IF
ALLOCATE(aexp(tw_nelems_exp(row_obj),nrhs),bexp(tw_nelems_exp(col_obj),nrhs))
DO i=1,nrhs
  IF(ASSOCIATED(row_obj%junc))THEN
    CALL tw_junction_expand(row_obj,a(:,i),aexp(:,i))
  ELSE
    aexp(:,i)=a(:,i)
  END IF
END DO
bexp=0.d0
CALL tw_compute_Lmat_MF_exp(row_obj,col_obj,nrhs,aexp,bexp)
DO i=1,nrhs
  IF(ASSOCIATED(col_obj%junc))THEN
    CALL tw_junction_reduce_vec(col_obj,bexp(:,i),b(:,i))
  ELSE
    b(:,i)=bexp(:,i)
  END IF
END DO
DEALLOCATE(aexp,bexp)
END SUBROUTINE tw_compute_Lmat_MF
END MODULE thin_wall
