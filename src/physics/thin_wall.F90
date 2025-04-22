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
USE oft_io, ONLY: xdmf_plot_file
USE oft_quadrature
USE oft_mesh_type, ONLY: oft_bmesh
USE oft_mesh_local, ONLY: bmesh_local_init
USE oft_mesh_local_util, ONLY: mesh_local_findedge
USE oft_io, ONLY: hdf5_rst, hdf5_write, hdf5_read, hdf5_rst_destroy, hdf5_create_file, &
  hdf5_field_get_sizes, hdf5_create_group
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
USE axi_green, ONLY: green, axi_coil_set
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
!> Structure containing filament coil definition
!---------------------------------------------------------------------------------
TYPE :: tw_gen_coil
  INTEGER(i4) :: npts = 0 !< Number of points in coil
  REAL(r8), POINTER, DIMENSION(:,:) :: pts => NULL() !< Points [3,npts]
END TYPE tw_gen_coil
!---------------------------------------------------------------------------------
!> Structure containing a coil sets composed of one or more individual filament coils
!---------------------------------------------------------------------------------
TYPE :: tw_coil_set
  LOGICAL :: sens_mask = .FALSE. !< Mask from sensor output?
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
  INTEGER(i4) :: nelems = 0 !< Number of elements in model (np_active+nholes+n_vcoils)
  INTEGER(i4) :: np_active = 0 !< Number of active vertices in model
  INTEGER(i4) :: n_vcoils = 0 !< Number of voltage-specified coils in model
  INTEGER(i4) :: n_icoils = 0 !< Number of current-specified coils in model
  INTEGER(i4) :: nholes = 0 !< Number of "holes" in model
  INTEGER(i4) :: nclosures = 0 !< Number of "closures" in model
  INTEGER(i4) :: nfh = 0 !< Number of face-hole interactions
  LOGICAL, POINTER, DIMENSION(:) :: sens_mask => NULL() !< Mask array for sensors [nreg]
  INTEGER(i4), POINTER, DIMENSION(:) :: pmap => NULL() !< Map from mesh vertices to active vertices [mesh%np]
  INTEGER(i4), POINTER, DIMENSION(:) :: kpmap_inv => NULL() !< Map from active vertices to mesh vertices [np_active]
  INTEGER(i4), POINTER, DIMENSION(:) :: lpmap_inv => NULL() !< Map from active vertices to mesh vertices [mesh%np]
  INTEGER(i4), POINTER, DIMENSION(:) :: closures => NULL() !< List of closure vertices [nclosures]
  INTEGER(i4), POINTER, DIMENSION(:) :: kfh => NULL() !< Pointer to face-hole interaction list [mesh%nc+1]
  INTEGER(i4), POINTER, DIMENSION(:,:) :: lfh => NULL() !< List of face-hole interactions [nfh]
  REAL(r8), POINTER, DIMENSION(:) :: Eta_reg => NULL() !< Resistivity*thickness values for each region [nreg]
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: Ael2dr => NULL() !< Element to driver (icoils) coupling matrix
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: Ael2coil => NULL() !< Element to coil (vcoils+icoils) coupling matrix
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: Ael2sen => NULL() !< Element to sensor coupling matrix
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: Acoil2coil => NULL() !< Coil to coil coupling matrix
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: Lmat => NULL() !< Full inductance matrix
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:) :: Adr2sen => NULL() !< Driver (icoils) to sensor coupling matrix
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: Bel => NULL()
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: Bdr => NULL()
  REAL(r8), POINTER, CONTIGUOUS, DIMENSION(:,:,:) :: qbasis => NULL() !< Basis function pre-evaluated at cell centers
  TYPE(xdmf_plot_file) :: xdmf
  CLASS(oft_vector), POINTER :: Uloc => NULL() !< FE vector for thin-wall model
  CLASS(oft_vector), POINTER :: Uloc_pts => NULL() !< Needs docs
  TYPE(oft_native_matrix), POINTER :: Rmat => NULL() !< Resistivity matrix for thin-wall model
  CLASS(oft_bmesh), POINTER :: mesh => NULL() !< Underlying surface mesh
  TYPE(hole_mesh), POINTER, DIMENSION(:) :: hmesh => NULL() !< Hole definitions
  TYPE(oft_1d_int), POINTER, DIMENSION(:) :: jumper_nsets => NULL() !< Jumper definitions
  TYPE(tw_coil_set), POINTER, DIMENSION(:) :: vcoils => NULL() !< List of Vcoils
  TYPE(tw_coil_set), POINTER, DIMENSION(:) :: icoils => NULL() !< List of Icoils
  TYPE(xml_node), POINTER :: xml => NULL()
CONTAINS
  !> Setup thin-wall model
  PROCEDURE :: setup => tw_setup
  !> Save debug information for model
  PROCEDURE :: save_debug => tw_save_debug
END TYPE tw_type

!---------------------------------------------------------------------------------
!> Class for thin-wall simulation
!---------------------------------------------------------------------------------
TYPE :: tw_plasma_boozer
  CLASS(tw_type), POINTER :: wall => NULL() !< Thin-wall model for structures
  CLASS(tw_type), POINTER :: plasma => NULL() !< Thin-wall model for plasma
  REAL(8) :: s = 0.d0
  REAL(8) :: alpha = 0.d0
  COMPLEX(8), POINTER, DIMENSION(:,:) :: rho => NULL()
  COMPLEX(8), POINTER, DIMENSION(:,:) :: L_w => NULL()
  COMPLEX(8), POINTER, DIMENSION(:,:) :: L_wc => NULL()
  COMPLEX(8), POINTER, DIMENSION(:,:) :: L_wd => NULL()
  COMPLEX(8), POINTER, DIMENSION(:,:) :: L_cw => NULL()
  COMPLEX(8), POINTER, DIMENSION(:,:) :: L_c => NULL()
  COMPLEX(8), POINTER, DIMENSION(:,:) :: L_cd => NULL()
  COMPLEX(8), POINTER, DIMENSION(:,:) :: L_dw => NULL()
  COMPLEX(8), POINTER, DIMENSION(:,:) :: L_dc => NULL()
  COMPLEX(8), POINTER, DIMENSION(:,:) :: L_d => NULL()
CONTAINS
  !> Setup thin-wall model
  PROCEDURE :: setup => tw_build_boozer
END TYPE tw_plasma_boozer
! REAL(r8) :: mag_dx = 1.d-5
REAL(r8) :: quad_tols(3) = [0.75d0, 0.95d0, 0.995d0] !< Distance tolerances for quadrature order selection
INTEGER(i4) :: quad_orders(3) = [18, 10, 6] !< Quadrature order for each tolerance
REAL(r8), PARAMETER :: target_err = 1.d-8
REAL(r8), PARAMETER :: coil_min_rad = 1.d-6
integer(i4), public, parameter :: tw_idx_ver=1 !< File version for array indexing
character(LEN=16), public, parameter :: tw_idx_path="ThinCurr_Version" !< HDF5 field name
CONTAINS
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE tw_setup(self,hole_ns)
CLASS(tw_type), INTENT(INOUT) :: self !< Thin-wall model object
TYPE(oft_1d_int), POINTER, INTENT(IN) :: hole_ns(:) !< Hole nodesets
INTEGER(4) :: i,j,k,l,face,ioffset,ed,error_flag
INTEGER(4), ALLOCATABLE :: kfh_tmp(:),np_inverse(:)
REAL(8) :: f(3),rgop(3,3),area_i,norm_i(3)
#ifdef HAVE_XML
TYPE(xml_node), POINTER :: coil_element
#endif
!
IF(ASSOCIATED(hole_ns))self%nholes=SIZE(hole_ns)
!---
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Creating thin-wall model'
CALL oft_increase_indent
CALL bmesh_local_init(self%mesh,sync_normals=.TRUE.)
#ifdef HAVE_XML
!---Load coils
IF(.NOT.ASSOCIATED(self%xml))THEN
  CALL xml_get_element(oft_env%xml,"thincurr",self%xml,error_flag)
  IF(error_flag/=0)CALL oft_warn('Unable to find "thincurr" XML node')
END IF
WRITE(*,'(2A)')oft_indent,'Loading V(t) driver coils'
CALL xml_get_element(self%xml,"vcoils",coil_element,error_flag)
IF(error_flag==0)CALL tw_load_coils(coil_element,self%n_vcoils,self%vcoils)
DO i=1,self%n_vcoils
  IF(ANY(self%vcoils(i)%res_per_len<0.d0))CALL oft_abort("Invalid resistivity for passive coil", &
    "tw_setup", __FILE__)
  IF(ANY(self%vcoils(i)%radius<coil_min_rad))CALL oft_abort("Invalid radius for passive coil", &
    "tw_setup", __FILE__)
END DO
WRITE(*,'(2A)')oft_indent,'Loading I(t) driver coils'
CALL xml_get_element(self%xml,"icoils",coil_element,error_flag)
IF(error_flag==0)CALL tw_load_coils(coil_element,self%n_icoils,self%icoils)
DO i=1,self%n_icoils
  DO j=1,self%icoils(i)%ncoils
    self%icoils(i)%radius(j)=MAX(coil_min_rad,self%icoils(i)%radius(j)) ! Remove dummy radius on Icoils
  END DO
END DO
#endif
WRITE(*,*)
! WRITE(*,'(2A)')oft_indent,'Thin-wall model loaded:'
! WRITE(*,*)'  filename    = ',TRIM(meshfile)
WRITE(*,'(2A,2I12)')oft_indent,'# of points    = ',self%mesh%np
WRITE(*,'(2A,I12)')oft_indent,'# of edges     = ',self%mesh%ne
WRITE(*,'(2A,I12)')oft_indent,'# of cells     = ',self%mesh%nc
WRITE(*,'(2A,I12)')oft_indent,'# of holes     = ',self%nholes
WRITE(*,'(2A,I12)')oft_indent,'# of Vcoils    = ',self%n_vcoils
WRITE(*,'(2A,I12)')oft_indent,'# of closures  = ',self%nclosures
IF(oft_debug_print(1))WRITE(*,*)oft_indent,'  Closures: ',self%closures
WRITE(*,'(2A,I12)')oft_indent,'# of Icoils    = ',self%n_icoils
!---Analyze mesh to construct holes
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Building holes'
ALLOCATE(self%hmesh(self%nholes))
ALLOCATE(self%kfh(self%mesh%nc+1))
self%kfh=0
IF(self%nholes>0)THEN
  DO i=1,self%nholes
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
  !---Check for duplicates
  DO i=1,self%nholes
    DO j=1,self%nholes
      IF(i==j.OR.self%hmesh(i)%n/=self%hmesh(j)%n)CYCLE
      IF(ALL(self%hmesh(i)%lp_sorted==self%hmesh(j)%lp_sorted))THEN
        WRITE(*,*)"Duplicate hole detected: ",i,j
        CALL oft_abort("Duplicate hole detected","tw_setup",__FILE__)
      END IF
    END DO
  END DO
END IF
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
CALL tw_load_eta(self)
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
END SELECT
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
INTEGER(4), INTENT(in) :: list_in(n) !< Input vertex list
INTEGER(4), INTENT(out) :: list_out(n) !< Reordered list
INTEGER(4), INTENT(in) :: n !< Number of points in list
INTEGER(4) :: ii,jj,k,l,nlinks,ipt,eprev,ed,ed2,ptp2,ptp,candidate,candidate2,last_item(3),i0
INTEGER(4), ALLOCATABLE :: lloop_tmp(:),flag_list(:)
!---Find loop points and edges
ALLOCATE(lloop_tmp(n),flag_list(n))
flag_list=0
lloop_tmp=list_in
CALL sort_array(lloop_tmp, n)
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
END SUBROUTINE tw_save_debug
!---------------------------------------------------------------------------------
!> Compute element to driver (Icoils) coupling matrix
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_Ael2dr(tw_obj,save_file)
TYPE(tw_type), INTENT(inout) :: tw_obj !< Thin-wall model object
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: save_file
LOGICAL :: exists
INTEGER(4) :: i,ii,j,jj,k,kk,ik,ncoils_tot,ih,ihp,ihc,file_counts(3),ierr,io_unit,iquad
REAL(8) :: tmp(3),cvec(3),cpt(3),pt_i(3),evec_i(3,3),pts_i(3,3),elapsed_time
REAL(8) :: pot_tmp,rgop(3,3),area_i,norm_i(3),f(3),dl_min,dl_max,pot_last
REAL(8), allocatable :: atmp(:,:),Ael2coil_tmp(:,:)
CLASS(oft_bmesh), POINTER :: bmesh
TYPE(oft_quad_type), ALLOCATABLE :: quads(:)
TYPE(tw_coil_set), POINTER, DIMENSION(:) :: coils_tot
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
!
IF(PRESENT(save_file))THEN
  IF(TRIM(save_file)/='none')THEN
    INQUIRE(FILE=TRIM(save_file),EXIST=exists)
    IF(exists)THEN
      WRITE(*,*)'Reading coil mutual matrices'
      OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
      READ(io_unit, IOSTAT=ierr)file_counts
      IF((ierr/=0).OR.ANY(file_counts/=[tw_obj%nelems,tw_obj%n_vcoils,tw_obj%n_icoils]))THEN
        exists=.FALSE.
        WRITE(*,*)'  Ignoring stored matrix: Sizes do not match'
      END IF
      IF(exists)THEN
        ALLOCATE(tw_obj%Ael2coil(tw_obj%nelems,tw_obj%n_vcoils))
        READ(io_unit, IOSTAT=ierr)tw_obj%Ael2coil
        IF(ierr/=0)THEN
          WRITE(*,*)'  Error reading matrix from file'
          DEALLOCATE(tw_obj%Ael2coil)
          exists=.FALSE.
        END IF
      END IF
      IF(exists)THEN
        ALLOCATE(tw_obj%Ael2dr(tw_obj%nelems,tw_obj%n_icoils))
        READ(io_unit, IOSTAT=ierr)tw_obj%Ael2dr
        IF(ierr/=0)THEN
          WRITE(*,*)'  Error reading matrix from file'
          DEALLOCATE(tw_obj%Ael2coil,tw_obj%Ael2dr)
          exists=.FALSE.
        END IF
      END IF
      CLOSE(io_unit)
    END IF
    IF(exists)THEN
      DEBUG_STACK_POP
      RETURN
    END IF
  END IF
END IF
!---Setup quadrature
ALLOCATE(quads(18))
DO i=1,18
  CALL set_quad_2d(quads(i),i)
END DO
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
WRITE(*,*)'Building coil<->element inductance matrices'
CALL mytimer%tick
ALLOCATE(Ael2coil_tmp(tw_obj%nelems,ncoils_tot))
Ael2coil_tmp=0.d0
f=1.d0/3.d0
!$omp parallel private(ii,ik,j,jj,k,kk,pts_i,pt_i,cvec, &
!$omp cpt,tmp,atmp,evec_i,pot_tmp,pot_last,i,ih,ihp, &
!$omp ihc,rgop,area_i,norm_i,dl_min,dl_max,iquad)
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
          DO ii=1,quads(iquad)%np
            pt_i = quads(iquad)%pts(1,ii)*pts_i(:,1) &
              + quads(iquad)%pts(2,ii)*pts_i(:,2) &
              + quads(iquad)%pts(3,ii)*pts_i(:,3)
            pot_tmp = pot_tmp + quads(iquad)%wts(ii)/SQRT(SUM((pt_i-cpt)**2))
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
!---Unpack passive and driver coils
ALLOCATE(tw_obj%Ael2coil(tw_obj%nelems,tw_obj%n_vcoils))
DO i=1,tw_obj%n_vcoils
  tw_obj%Ael2coil(:,i)=Ael2coil_tmp(:,i)
END DO
ALLOCATE(tw_obj%Ael2dr(tw_obj%nelems,tw_obj%n_icoils))
DO i=1,tw_obj%n_icoils
  tw_obj%Ael2dr(:,i)=Ael2coil_tmp(:,i+tw_obj%n_vcoils)
END DO
!---Cleanup temporaries
DEALLOCATE(Ael2coil_tmp,coils_tot)
DO i=1,18
  CALL quads(i)%delete()
END DO
DEALLOCATE(quads)
elapsed_time=mytimer%tock()
WRITE(*,'(5X,2A)')'Time = ',time_to_string(elapsed_time)
!
CALL tw_compute_Lmat_coils(tw_obj)
!
IF(PRESENT(save_file))THEN
  IF(TRIM(save_file)/='none')THEN
    OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
    WRITE(io_unit)tw_obj%nelems,tw_obj%n_vcoils,tw_obj%n_icoils
    WRITE(io_unit)tw_obj%Ael2coil
    WRITE(io_unit)tw_obj%Ael2dr
    CLOSE(io_unit)
  END IF
END IF
DEBUG_STACK_POP
END SUBROUTINE tw_compute_Ael2dr
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
REAL(8) :: rgop(3,3),norm_i(3),area_i,f(3),dl_min,dl_max,pot_last,pot_tmp
REAL(8), allocatable :: atmp(:,:),Acoil2sen_tmp(:,:),Acoil2coil_tmp(:,:)
CLASS(oft_bmesh), POINTER :: bmesh
INTEGER(4) :: ncoils_tot
TYPE(tw_coil_set), POINTER, DIMENSION(:) :: coils_tot
TYPE(oft_quad_type), ALLOCATABLE :: quads(:)
DEBUG_STACK_PUSH
!
ALLOCATE(quads(18))
DO i=1,18
  CALL set_quad_2d(quads(i),i)
END DO
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
!---Compute coupling between vertices and sensors
f=1.d0/3.d0
!---Compute coupling between coils
ALLOCATE(Acoil2coil_tmp(tw_obj%n_vcoils,ncoils_tot))
Acoil2coil_tmp=0.d0
WRITE(*,*)'Building coil<->coil inductance matrix'
IF(tw_obj%n_vcoils>0.AND.ncoils_tot>0)THEN
  sqrt_e=SQRT(EXP(1.d0))
  !$omp parallel private(i,ii,j,k,kk,tmp,pt_i,rvec_i,atmp,cvec,cpt,coil_thickness)
  ALLOCATE(atmp(ncoils_tot,1))
  !$omp do
  DO l=1,tw_obj%n_vcoils
    atmp=0.d0
    DO i=1,coils_tot(l)%ncoils
      DO ii=2,coils_tot(l)%coils(i)%npts
        rvec_i = coils_tot(l)%coils(i)%pts(:,ii)-coils_tot(l)%coils(i)%pts(:,ii-1)
        pt_i = (coils_tot(l)%coils(i)%pts(:,ii)+coils_tot(l)%coils(i)%pts(:,ii-1))/2.d0
        DO j=1,ncoils_tot
          DO k=1,coils_tot(j)%ncoils
            coil_thickness=(MAX(coils_tot(l)%radius(i),coils_tot(j)%radius(k))**2)/sqrt_e
            tmp=0.d0
            DO kk=2,coils_tot(j)%coils(k)%npts
              cvec = coils_tot(j)%coils(k)%pts(:,kk)-coils_tot(j)%coils(k)%pts(:,kk-1)
              cpt = (coils_tot(j)%coils(k)%pts(:,kk)+coils_tot(j)%coils(k)%pts(:,kk-1))/2.d0
              tmp = tmp + DOT_PRODUCT(rvec_i,cvec)/SQRT(SUM((pt_i-cpt)**2) + coil_thickness)
            END DO
            atmp(j,1)=atmp(j,1)+coils_tot(j)%scales(k)*coils_tot(l)%scales(i)*tmp
          END DO
        END DO
      END DO
    END DO
    DO j=1,ncoils_tot
      Acoil2coil_tmp(l,j) = atmp(j,1)
    END DO
  END DO
  DEALLOCATE(atmp)
  !$omp end parallel
END IF
DEALLOCATE(coils_tot)
!---Unpack passive and driver coils
IF(ASSOCIATED(tw_obj%Acoil2coil))DEALLOCATE(tw_obj%Acoil2coil)
ALLOCATE(tw_obj%Acoil2coil(tw_obj%n_vcoils,tw_obj%n_vcoils))
DO i=1,tw_obj%n_vcoils
  tw_obj%Acoil2coil(:,i)=Acoil2coil_tmp(:,i)
  tw_obj%vcoils(i)%Lself=Acoil2coil_tmp(i,i)
  WRITE(*,"(A,1X,I4,A,ES12.4)")"  Vcoil",i,": L [H] = ",tw_obj%vcoils(i)%Lself*1.d-7
END DO
!---Compute coupling between elements and drivers
IF(ASSOCIATED(tw_obj%Ael2dr))THEN
  ! WRITE(*,*)'Building driver->element inductance matrices'
  IF(tw_obj%n_vcoils>0)THEN
    !$omp parallel private(j,ii,jj,ik,jk)
    !$omp do
    DO i=1,tw_obj%n_vcoils
      DO jj=1,tw_obj%n_icoils
        tw_obj%Ael2dr(tw_obj%np_active+tw_obj%nholes+i,jj) = Acoil2coil_tmp(i,jj+tw_obj%n_vcoils)
      END DO
    END DO
    !$omp end parallel
  END IF
  tw_obj%Ael2dr = tw_obj%Ael2dr*mu0/(4.d0*pi)
END IF
DEALLOCATE(Acoil2coil_tmp)
!
DO i=1,18
  CALL quads(i)%delete()
END DO
DEALLOCATE(quads)
DEBUG_STACK_POP
END SUBROUTINE tw_compute_Lmat_coils
!---------------------------------------------------------------------------------
!> Compute mutual inductance matrix between two thin-wall models
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_LmatDirect(row_model,Lmat,col_model,save_file)
TYPE(tw_type), TARGET, INTENT(in) :: row_model !< Thin-wall model object for rows
REAL(8), CONTIGUOUS, POINTER, INTENT(inout) :: Lmat(:,:) !< Mutual inductance matrix
TYPE(tw_type), TARGET, OPTIONAL, INTENT(in) :: col_model !< Thin-wall model object for columns
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: save_file
INTEGER(4) :: i,ii,j,jj,ik,jk,imin,jmax,io_unit,ierr,file_counts(6),hash_tmp(6),iquad
INTEGER(8) :: counts(4)
REAL(8) :: tmp,dl_min,dl_max,f(3),elapsed_time
REAL(8) :: pt_i(3),rgop(3,3),area_i,norm_i(3),evec_i(3,3),pts_i(3,3)
REAL(8) :: pt_j(3),cgop(3,3),area_j,norm_j(3),evec_j(3,3),pts_j(3,3)
LOGICAL :: vvclose_flag,close_flag,exists,Lself
CLASS(oft_bmesh), POINTER :: rmesh,cmesh
TYPE(tw_type), POINTER :: row_obj,col_obj
TYPE(oft_quad_type), ALLOCATABLE :: quads(:)
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
IF(ASSOCIATED(Lmat))THEN
  IF(SIZE(Lmat,DIM=1)/=col_obj%nelems.OR.SIZE(Lmat,DIM=2)/=row_obj%nelems)THEN
    CALL oft_abort("Input size incorrect","tw_compute_LmatDirect",__FILE__)
  END IF
ELSE
  ALLOCATE(Lmat(col_obj%nelems,row_obj%nelems))
END IF
!
IF(PRESENT(save_file))THEN
  hash_tmp(3) = oft_simple_hash(C_LOC(row_obj%mesh%lc),INT(4*3*row_obj%mesh%nc,8))
  hash_tmp(4) = oft_simple_hash(C_LOC(col_obj%mesh%lc),INT(4*3*col_obj%mesh%nc,8))
  hash_tmp(5) = oft_simple_hash(C_LOC(row_obj%mesh%r),INT(8*3*row_obj%mesh%np,8))
  hash_tmp(6) = oft_simple_hash(C_LOC(col_obj%mesh%r),INT(8*3*col_obj%mesh%np,8))
  IF(TRIM(save_file)/='none')THEN
    INQUIRE(FILE=TRIM(save_file),EXIST=exists)
    IF(exists)THEN
      IF(Lself)THEN
        WRITE(*,*)'Reading element<->element self inductance matrix'
        OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
        READ(io_unit, IOSTAT=ierr)file_counts
        hash_tmp(1) = row_obj%nelems
        hash_tmp(2) = row_obj%mesh%nc
        IF((ierr==0).AND.ALL(file_counts==hash_tmp))THEN
          DO i=1,row_obj%nelems
            READ(io_unit, IOSTAT=ierr)Lmat(i,i:row_obj%nelems)
            IF(ierr/=0)THEN
              WRITE(*,*)'  Error reading matrix from file'
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
          WRITE(*,*)'  Ignoring stored matrix: Model hashes do not match'
          CLOSE(io_unit)
          exists=.FALSE.
        END IF
      ELSE
        WRITE(*,*)'Reading element<->element mutual inductance matrix'
        OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
        READ(io_unit, IOSTAT=ierr)file_counts
        hash_tmp(1) = col_obj%nelems
        hash_tmp(2) = row_obj%nelems
        IF((ierr==0).AND.ALL(file_counts==hash_tmp))THEN
          DO i=1,row_obj%nelems
            READ(io_unit, IOSTAT=ierr)Lmat(:,i)
            IF(ierr/=0)THEN
              WRITE(*,*)'  Error reading matrix from file'
              exists=.FALSE.
              EXIT
            END IF
          END DO
          CLOSE(io_unit)
        ELSE
          WRITE(*,*)'  Ignoring stored matrix: Model hashes do not match'
          CLOSE(io_unit)
          exists=.FALSE.
        END IF
      END IF
    END IF
    IF(exists)THEN
      DEBUG_STACK_POP
      RETURN
    END IF
  END IF
END IF
!
IF(Lself)THEN
  WRITE(*,*)'Building element<->element self inductance matrix'
ELSE
  WRITE(*,*)'Building element<->element mutual inductance matrix'
END IF
Lmat=0.d0
CALL mytimer%tick
!---Setup quadrature
ALLOCATE(quads(18))
DO i=1,18
  CALL set_quad_2d(quads(i),i)
END DO
!
rmesh=>row_obj%mesh
cmesh=>col_obj%mesh
!---
counts=0
f=1.d0/3.d0
!$omp parallel private(ii,j,jj,ik,jk,pts_i,pts_j,tmp,iquad, &
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
    !---Compute minimum separation
    dl_min=1.d99
    dl_max=SQRT(MAX(area_i,area_j)*2.d0) !-1.d99
    !$omp simd collapse(1) reduction(max:dl_max) reduction(min:dl_min)
    DO ii=1,3
      DO jj=1,3
        dl_min=MIN(dl_min,SQRT(SUM((pts_i(:,ii)-pts_j(:,jj))**2)))
        dl_max=MAX(dl_max,SQRT(SUM((pts_i(:,ii)-pts_j(:,jj))**2)))
      END DO
    END DO
    !---Chose quadrature order based on distance
    IF(dl_min<1.d-8)THEN
      iquad = 18
    ELSE
      iquad = MAX(4,MIN(18,ABS(INT(LOG(target_err)/LOG(1.d0-dl_min/dl_max)))))
    END IF
    tmp=0.d0
    IF(iquad>10)THEN
      DO jj=1,quads(iquad)%np
        pt_j = quads(iquad)%pts(1,jj)*pts_j(:,1) &
          + quads(iquad)%pts(2,jj)*pts_j(:,2) &
          + quads(iquad)%pts(3,jj)*pts_j(:,3)
        tmp = tmp + quads(iquad)%wts(jj)*tw_compute_phipot(pts_i,pt_j)
      END DO
      tmp=tmp*area_j
    ELSE
      !$omp simd collapse(1) private(pt_i,pt_j) reduction(+:tmp)
      DO ii=1,quads(iquad)%np
        pt_i = quads(iquad)%pts(1,ii)*pts_i(:,1) &
          + quads(iquad)%pts(2,ii)*pts_i(:,2) &
          + quads(iquad)%pts(3,ii)*pts_i(:,3)
        DO jj=1,quads(iquad)%np
          pt_j = quads(iquad)%pts(1,jj)*pts_j(:,1) &
            + quads(iquad)%pts(2,jj)*pts_j(:,2) &
            + quads(iquad)%pts(3,jj)*pts_j(:,3)
          tmp = tmp + quads(iquad)%wts(jj)*quads(iquad)%wts(ii)/SQRT(SUM((pt_i-pt_j)**2))
        END DO
      END DO
      tmp=tmp*area_i*area_j
    END IF
    !
    DO ii=1,3
      ik=row_obj%pmap(rmesh%lc(ii,i))
      IF(ik==0)CYCLE
      DO jj=1,3
        jk=col_obj%pmap(cmesh%lc(jj,j))
        IF((Lself.AND.jk<ik).OR.jk<=0)CYCLE
        dl_min = DOT_PRODUCT(evec_i(:,ii),evec_j(:,jj))*tmp
        !$omp atomic
        Lmat(jk,ik) = Lmat(jk,ik) + dl_min
      END DO
      DO jj=col_obj%kfh(j),col_obj%kfh(j+1)-1
        jk=ABS(col_obj%lfh(1,jj))+col_obj%np_active
        dl_min = SIGN(1,col_obj%lfh(1,jj))* &
          DOT_PRODUCT(evec_i(:,ii),evec_j(:,col_obj%lfh(2,jj)))*tmp
        !$omp atomic
        Lmat(jk,ik) = Lmat(jk,ik) + dl_min
      END DO
    END DO
    DO ii=row_obj%kfh(i),row_obj%kfh(i+1)-1
      ik=ABS(row_obj%lfh(1,ii))+row_obj%np_active
      DO jj=col_obj%kfh(j),col_obj%kfh(j+1)-1
        jk=ABS(col_obj%lfh(1,jj))+col_obj%np_active
        IF(Lself.AND.jk<ik)CYCLE
        dl_min = SIGN(1,row_obj%lfh(1,ii))*SIGN(1,col_obj%lfh(1,jj))* &
          DOT_PRODUCT(evec_i(:,row_obj%lfh(2,ii)),evec_j(:,col_obj%lfh(2,jj)))*tmp
        !$omp atomic
        Lmat(jk,ik) = Lmat(jk,ik) + dl_min
      END DO
      !
      IF(.NOT.Lself)THEN
        DO jj=1,3
          jk=col_obj%pmap(cmesh%lc(jj,j))
          IF(jk<=0)CYCLE
          dl_min = SIGN(1,row_obj%lfh(1,ii))* &
            DOT_PRODUCT(evec_i(:,row_obj%lfh(2,ii)),evec_j(:,jj))*tmp
          !$omp atomic
          Lmat(jk,ik) = Lmat(jk,ik) + dl_min
        END DO
      END IF
    END DO
  END DO
END DO
!$omp end do nowait
IF(Lself)THEN
  !---Add passive coils to model
  !$omp do
  DO i=1,row_obj%np_active+row_obj%nholes
    DO j=1,row_obj%n_vcoils
      jj=row_obj%np_active+row_obj%nholes+j
      Lmat(jj,i) = row_obj%Ael2coil(i,j)
    END DO
  END DO
  !$omp end do nowait
  !$omp do
  DO i=1,row_obj%n_vcoils
    ii=row_obj%np_active+row_obj%nholes+i
    DO j=i,row_obj%n_vcoils
      jj=row_obj%np_active+row_obj%nholes+j
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
DO i=1,18
  CALL quads(i)%delete()
END DO
DEALLOCATE(quads)
elapsed_time=mytimer%tock()
WRITE(*,'(5X,2A)')'Time = ',time_to_string(elapsed_time)
IF(PRESENT(save_file))THEN
  IF(TRIM(save_file)/='none')THEN
    IF(Lself)THEN
      OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
      hash_tmp(1) = row_obj%nelems
      hash_tmp(2) = row_obj%mesh%nc
      WRITE(io_unit)hash_tmp
      DO i=1,row_obj%nelems
        WRITE(io_unit)row_obj%Lmat(i,i:row_obj%nelems)
      END DO
      CLOSE(io_unit)
    ELSE
      OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
      hash_tmp(1) = col_obj%nelems
      hash_tmp(2) = row_obj%nelems
      WRITE(io_unit)hash_tmp
      DO i=1,row_obj%nelems
        WRITE(io_unit)Lmat(:,i)
      END DO
      CLOSE(io_unit)
    END IF
  END IF
END IF
DEBUG_STACK_POP
END SUBROUTINE tw_compute_LmatDirect
!---------------------------------------------------------------------------------
!> Compute mutual inductance matrix between two thin-wall models
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_Lmat_MF(row_obj,col_obj,nrhs,a,b)
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
WRITE(*,*)'Applying MF element<->element inductance matrix'
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
WRITE(*,'(5X,2A)')'Time = ',time_to_string(elapsed_time)
DEBUG_STACK_POP
END SUBROUTINE tw_compute_Lmat_MF
!---------------------------------------------------------------------------------
!> Compute coupling from thin-wall model elements to flux loop sensors
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_mutuals(tw_obj,nsensors,sensors,save_file)
TYPE(tw_type), INTENT(inout) :: tw_obj !< Thin-wall model object
INTEGER(4), INTENT(in) :: nsensors !< Number of flux loops
TYPE(floop_sensor), POINTER, INTENT(in) :: sensors(:) !< List of flux loops
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: save_file
LOGICAL :: exists
INTEGER(4) :: i,ii,j,jj,k,kk,l,ik,jk,ih,ihp,ihc,file_counts(4),ierr,io_unit,iquad
REAL(8) :: tmp,dl,pt_i(3),pt_j(3),evec_i(3,3),evec_j(3,3),pts_i(3,3),pt_i_last(3)
REAL(8) :: rvec_i(3),r1,rmag,rvec_j(3),z1,cvec(3),cpt(3),coil_thickness
REAL(8) :: rgop(3,3),norm_i(3),area_i,f(3),dl_min,dl_max,pot_last,pot_tmp
REAL(8), allocatable :: atmp(:,:),Acoil2sen_tmp(:,:),Acoil2coil_tmp(:,:)
CLASS(oft_bmesh), POINTER :: bmesh
INTEGER(4) :: ncoils_tot
TYPE(tw_coil_set), POINTER, DIMENSION(:) :: coils_tot
TYPE(oft_quad_type), ALLOCATABLE :: quads(:)
DEBUG_STACK_PUSH
!
IF(PRESENT(save_file))THEN
  IF(TRIM(save_file)/='none')THEN
    INQUIRE(FILE=TRIM(save_file),EXIST=exists)
    IF(exists)THEN
      WRITE(*,*)'Reading sensor mutual matrices'
      OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
      READ(io_unit, IOSTAT=ierr)file_counts
      IF((ierr/=0).OR.ANY(file_counts/=[tw_obj%nelems,tw_obj%n_vcoils,tw_obj%n_icoils,nsensors]))THEN
        exists=.FALSE.
        WRITE(*,*)'  Ignoring stored matrix: Sizes do not match'
      END IF
    IF(exists)THEN
      ALLOCATE(tw_obj%Ael2sen(tw_obj%nelems,nsensors))
      READ(io_unit, IOSTAT=ierr)tw_obj%Ael2sen
      IF(ierr/=0)THEN
        WRITE(*,*)'  Error reading matrix from file'
        DEALLOCATE(tw_obj%Ael2sen)
        exists=.FALSE.
      END IF
    END IF
    ! IF(exists)THEN
    !   ALLOCATE(tw_obj%Acoil2sen(tw_obj%n_vcoils,nsensors))
    !   READ(io_unit, IOSTAT=ierr)tw_obj%Acoil2sen
    !   IF(ierr/=0)THEN
    !     WRITE(*,*)'  Error reading matrix from file'
    !     DEALLOCATE(tw_obj%Ael2sen,tw_obj%Acoil2sen)
    !     exists=.FALSE.
    !   END IF
    ! END IF
    IF(exists)THEN
      ALLOCATE(tw_obj%Adr2sen(tw_obj%n_icoils,nsensors))
      READ(io_unit, IOSTAT=ierr)tw_obj%Adr2sen
      IF(ierr/=0)THEN
        WRITE(*,*)'  Error reading matrix from file'
        DEALLOCATE(tw_obj%Ael2sen,tw_obj%Adr2sen)
        exists=.FALSE.
      END IF
    END IF
      CLOSE(io_unit)
    END IF
    IF(exists)THEN
      DEBUG_STACK_POP
      RETURN
    END IF
  END IF
END IF
!
ALLOCATE(quads(18))
DO i=1,18
  CALL set_quad_2d(quads(i),i)
END DO
!---Build temporary coil set from all filaments in model
ncoils_tot=tw_obj%n_vcoils+tw_obj%n_icoils
ALLOCATE(coils_tot(ncoils_tot))
DO i=1,tw_obj%n_vcoils
  CALL tw_copy_coil(tw_obj%vcoils(i),coils_tot(i))
END DO
DO i=1,tw_obj%n_icoils
  CALL tw_copy_coil(tw_obj%icoils(i),coils_tot(i+tw_obj%n_vcoils))
END DO
bmesh=>tw_obj%mesh
!---Compute coupling between vertices and sensors
f=1.d0/3.d0
WRITE(*,*)'Building element->sensor inductance matrix'
IF(ASSOCIATED(tw_obj%Ael2sen))DEALLOCATE(tw_obj%Ael2sen)
ALLOCATE(tw_obj%Ael2sen(nsensors,tw_obj%nelems))
tw_obj%Ael2sen=0.d0
IF(nsensors>0)THEN
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
          DO ii=1,quads(iquad)%np
            pt_i = quads(iquad)%pts(1,ii)*pts_i(:,1) &
              + quads(iquad)%pts(2,ii)*pts_i(:,2) &
              + quads(iquad)%pts(3,ii)*pts_i(:,3)
            tmp = tmp + quads(iquad)%wts(ii)/SQRT(SUM((pt_i-pt_j)**2))
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
END IF
!---Compute coupling between coils and sensors
ALLOCATE(Acoil2sen_tmp(nsensors,ncoils_tot))
Acoil2sen_tmp=0.d0
WRITE(*,*)'Building coil->sensor inductance matrix'
IF(nsensors>0.AND.ncoils_tot>0)THEN
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
END IF
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
DO i=1,18
  CALL quads(i)%delete()
END DO
DEALLOCATE(quads)
!---
IF(PRESENT(save_file))THEN
  IF(TRIM(save_file)/='none')THEN
    OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
    WRITE(io_unit)tw_obj%nelems,tw_obj%n_vcoils,tw_obj%n_icoils,nsensors
    WRITE(io_unit)tw_obj%Ael2sen
    ! WRITE(io_unit)tw_obj%Acoil2sen
    WRITE(io_unit)tw_obj%Adr2sen
    CLOSE(io_unit)
  END IF
END IF
DEBUG_STACK_POP
END SUBROUTINE tw_compute_mutuals
!---------------------------------------------------------------------------------
!> Compute resistivity matrix for thin-wall model
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_Rmat(tw_obj,keep_closures)
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
WRITE(*,*)'Building resistivity matrix'
ALLOCATE(tw_obj%Rmat)
Rmat=>tw_obj%Rmat
Rmat%nr=tw_obj%nelems; Rmat%nrg=tw_obj%nelems
Rmat%nc=tw_obj%nelems; Rmat%ncg=tw_obj%nelems
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
  eta_eff=tw_obj%Eta_reg(bmesh%reg(i))
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
    DO k=2,tw_obj%vcoils(i)%coils(j)%npts
      dl = dl + SQRT(SUM((tw_obj%vcoils(i)%coils(j)%pts(:,k)-tw_obj%vcoils(i)%coils(j)%pts(:,k-1))**2))
    END DO
    tw_obj%vcoils(i)%Rself = tw_obj%vcoils(i)%Rself + tw_obj%vcoils(i)%res_per_len(j)*dl
  END DO
  WRITE(*,"(A,1X,I4,A,ES12.4)")"  Vcoil",i,": R [Ohm] = ",tw_obj%vcoils(i)%Rself !*pi*4.d-7
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
DEBUG_STACK_POP
END SUBROUTINE tw_compute_Rmat
!---------------------------------------------------------------------------------
!> Compute \f$ \int 1/(r-r') dA' \f$ for a triangle
!---------------------------------------------------------------------------------
FUNCTION tw_compute_phipot(pt_cell,pt) RESULT(phi)
REAL(8), INTENT(in) :: pt_cell(3,3) !< Vertices defining triangle
REAL(8), INTENT(in) :: pt(3) !< Observation point (r)
REAL(8) :: phi
REAL(8) :: gam(3),num,den,omega,nhat(3),c(3,3),r(3,3),rmag(3),tmp,yp(3),yp1(3)
INTEGER(4) :: i
DEBUG_STACK_PUSH
!---Get unit normal
nhat = cross_product(pt_cell(:,2)-pt_cell(:,1),pt_cell(:,3)-pt_cell(:,2))
nhat = nhat/SQRT(DOT_PRODUCT(nhat,nhat))
!---Compute local vectors
DO i=1,3
  r(:,i)=pt_cell(:,i)-pt
  rmag(i)=SQRT(DOT_PRODUCT(r(:,i),r(:,i)))
  c(:,i)=r(:,i)-DOT_PRODUCT(nhat,r(:,i))*nhat
END DO
!---Compute solid angle
num=DOT_PRODUCT(r(:,1),cross_product(r(:,2),r(:,3)))
den=PRODUCT(rmag) + DOT_PRODUCT(r(:,1),r(:,2))*rmag(3) &
  + DOT_PRODUCT(r(:,1),r(:,3))*rmag(2) + DOT_PRODUCT(r(:,2),r(:,3))*rmag(1)
omega=2.d0*ATAN2(num,den)
!---Compute gammas
DO i=1,3
  yp=r(:,i)
  IF(i==3)THEN
    yp1=r(:,1)
  ELSE
    yp1=r(:,i+1)
  END IF
  tmp=magnitude(yp1-yp)
  num=magnitude(yp1)*tmp + DOT_PRODUCT(yp1,yp1-yp)
  den=magnitude(yp)*tmp + DOT_PRODUCT(yp,yp1-yp)
  IF(ABS(den)<1.d-14.OR.tmp<1.d-14)THEN
    gam(i)=0.d0
  ELSE
    gam(i)=LOG(num/den)/tmp
  END IF
END DO
!---Evalute final function
phi=0.d0
DO i=1,3
  yp=c(:,i)
  IF(i==3)THEN
    yp1=c(:,1)
  ELSE
    yp1=c(:,i+1)
  END IF
  phi=phi+DOT_PRODUCT(nhat,cross_product(yp,yp1))*gam(i)
END DO
phi=phi-DOT_PRODUCT(nhat,r(:,1))*omega
DEBUG_STACK_POP
END FUNCTION tw_compute_phipot
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE tw_compute_Bops(self,save_file)
TYPE(tw_type), INTENT(inout) :: self
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: save_file
REAL(r8) :: evec_i(3,3),evec_j(3),pts_i(3,3),pt_i(3),pt_j(3),diffvec(3),ecc(3)
REAL(r8) :: r1,z1,rmag,cvec(3),cpt(3),tmp,area_i,dl_min,dl_max,norm_j(3),f(3),pot_tmp,pot_last
REAL(r8), ALLOCATABLE :: atmp(:,:,:)
REAL(8), PARAMETER :: B_dx = 1.d-6
INTEGER(4) :: i,ii,j,jj,ik,jk,k,kk,iquad,hash_tmp(4),file_counts(4)
LOGICAL :: is_neighbor,exists
CLASS(oft_bmesh), POINTER :: bmesh
TYPE(oft_quad_type), ALLOCATABLE :: quads(:)
IF(TRIM(save_file)/='none')THEN
  INQUIRE(FILE=TRIM(save_file),EXIST=exists)
  IF(exists)THEN
    hash_tmp(1) = self%nelems
    hash_tmp(2) = self%mesh%nc
    hash_tmp(3) = oft_simple_hash(C_LOC(self%mesh%lc),INT(4*3*self%mesh%nc,8))
    hash_tmp(4) = oft_simple_hash(C_LOC(self%mesh%r),INT(8*3*self%mesh%np,8))
    WRITE(*,*)'  Loading B-field operator from file: ',TRIM(save_file)
    CALL hdf5_read(file_counts,TRIM(save_file),'MODEL_hash')
    IF(exists.AND.ALL(file_counts==hash_tmp))THEN
      ALLOCATE(self%Bel(self%nelems,self%mesh%np,3))
      CALL hdf5_read(self%Bel(:,:,1),TRIM(save_file),'Bel_X',success=exists)
      IF(exists)CALL hdf5_read(self%Bel(:,:,2),TRIM(save_file),'Bel_Y',success=exists)
      IF(exists)CALL hdf5_read(self%Bel(:,:,3),TRIM(save_file),'Bel_Z',success=exists)
      IF(exists)THEN
        ALLOCATE(self%Bdr(self%mesh%np,self%n_icoils,3))
        CALL hdf5_read(self%Bdr(:,:,1),TRIM(save_file),'Bdr_X',success=exists)
        IF(exists)CALL hdf5_read(self%Bdr(:,:,2),TRIM(save_file),'Bdr_X',success=exists)
        IF(exists)CALL hdf5_read(self%Bdr(:,:,3),TRIM(save_file),'Bdr_X',success=exists)
      END IF
    END IF
  END IF
  IF(exists)RETURN
  DEALLOCATE(self%Bel,self%Bdr)
END IF
!
bmesh=>self%mesh
ALLOCATE(quads(18))
DO i=1,18
  CALL set_quad_2d(quads(i),i)
END DO
ALLOCATE(self%Bel(self%nelems,bmesh%np,3))
self%Bel=0.d0
f=1.d0/3.d0
WRITE(*,*)'Building element->element magnetic reconstruction operator'
!$omp parallel private(ii,j,jj,ik,pts_i,tmp,pt_i,pt_j,evec_i, &
!$omp atmp,i,area_i,dl_min,dl_max,norm_j,diffvec,is_neighbor,iquad)
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
    pt_j=bmesh%r(:,j)
    !---Compute minimum separation
    dl_min=1.d99
    dl_max=SQRT(MAX(area_i,bmesh%va(j)/(pi**2))) !-1.d99
    !!$omp simd collapse(1) reduction(max:dl_max) reduction(min:dl_min)
    DO ii=1,3
      dl_min=MIN(dl_min,SQRT(SUM((pts_i(:,ii)-pt_j)**2)))
      dl_max=MAX(dl_max,SQRT(SUM((pts_i(:,ii)-pt_j)**2)))
    END DO
    !---Chose quadrature order based on distance
    IF(dl_min<1.d-8)THEN
      iquad = 18
      is_neighbor=.TRUE.
    ELSE
      iquad = MAX(4,MIN(18,ABS(INT(LOG(target_err)/LOG(1.d0-dl_min/dl_max)))))
      is_neighbor=.FALSE.
    END IF
    !
    IF(iquad>10)THEN
      IF(is_neighbor)THEN
        ! pt_j = (SUM(pts_i,DIM=2)/3.d0 + 19.d0*pt_j)/20.d0 ! Move incrementally toward cell center
        pt_j = pt_j - norm_j*10.d0*B_dx ! Sample just inside face
      END IF
      diffvec=0.d0
      DO ik=1,2
        IF(ik==2)pt_j=pt_j+norm_j*20.d0*B_dx ! Sample just outside face
        DO jj=1,3
          ! Forward step
          pt_j(jj)=pt_j(jj)+B_dx
          tmp=tw_compute_phipot(pts_i,pt_j)
          diffvec(jj)=diffvec(jj)+tmp/(2.d0*B_dx)
          ! Backward step
          pt_j(jj)=pt_j(jj)-2.d0*B_dx
          tmp=tw_compute_phipot(pts_i,pt_j)
          diffvec(jj)=diffvec(jj)-tmp/(2.d0*B_dx)
          pt_j(jj)=pt_j(jj)+B_dx ! Reset point
        END DO
        IF(.NOT.is_neighbor)EXIT
      END DO
      IF(is_neighbor)diffvec=diffvec/2.d0
      DO ik=1,3
        atmp(1,ik,j) = diffvec(2)*evec_i(3,ik) - diffvec(3)*evec_i(2,ik)
        atmp(2,ik,j) = diffvec(3)*evec_i(1,ik) - diffvec(1)*evec_i(3,ik)
        atmp(3,ik,j) = diffvec(1)*evec_i(2,ik) - diffvec(2)*evec_i(1,ik)
      END DO
    ELSE
      DO ik=1,3
        diffvec=0.d0
        ! !$omp simd private(pt_i) reduction(+:diffvec)
        DO ii=1,quads(iquad)%np
          pt_i = pt_j - (quads(iquad)%pts(1,ii)*pts_i(:,1) &
            + quads(iquad)%pts(2,ii)*pts_i(:,2) &
            + quads(iquad)%pts(3,ii)*pts_i(:,3))
          diffvec = diffvec + cross_product(evec_i(:,ik),pt_i)*quads(iquad)%wts(ii)/(SUM(pt_i**2))**1.5d0
        END DO
        atmp(:,ik,j) = diffvec*area_i
      END DO
    END IF
  END DO
  DO ii=1,3
    ik=self%pmap(bmesh%lc(ii,i))
    IF(ik==0)CYCLE
    DO j=1,bmesh%np
      DO jj=1,3
        !$omp atomic
        self%Bel(ik,j,jj) = self%Bel(ik,j,jj) + atmp(jj,ii,j)
      END DO
    END DO
  END DO
  DO ii=self%kfh(i),self%kfh(i+1)-1
    ik=ABS(self%lfh(1,ii))+self%np_active
    DO j=1,bmesh%np
      DO jj=1,3
        tmp=SIGN(1,self%lfh(1,ii))*atmp(jj,self%lfh(2,ii),j)
        !$omp atomic
        self%Bel(ik,j,jj) = self%Bel(ik,j,jj) + tmp
      END DO
    END DO
  END DO
END DO
DEALLOCATE(atmp)
!$omp end parallel
DO i=1,18
  CALL quads(i)%delete()
END DO
DEALLOCATE(quads)
!
WRITE(*,*)'Building vcoil->element magnetic reconstruction operator'
!$omp parallel do private(ii,j,k,kk,pt_j,ecc,diffvec,cvec,cpt,pot_tmp,pot_last)
DO i=1,bmesh%np
  pt_j=bmesh%r(:,i)
  !---Compute driver contributions
  DO j=1,self%n_vcoils
    ecc = 0.d0
    IF(self%vcoils(j)%sens_mask)CYCLE
    DO k=1,self%vcoils(j)%ncoils
      diffvec=0.d0
      DO kk=2,self%vcoils(j)%coils(k)%npts
        cvec = self%vcoils(j)%coils(k)%pts(:,kk)-self%vcoils(j)%coils(k)%pts(:,kk-1)
        cpt = (self%vcoils(j)%coils(k)%pts(:,kk)+self%vcoils(j)%coils(k)%pts(:,kk-1))/2.d0
        diffvec = diffvec + cross_product(cvec,pt_j-cpt)/SUM((pt_j-cpt)**2)**1.5d0
      END DO
      ecc=ecc+self%vcoils(j)%scales(k)*diffvec
    END DO
    DO jj=1,3
      !$omp atomic
      self%Bel(self%np_active+self%nholes+j,i,jj) = &
        self%Bel(self%np_active+self%nholes+j,i,jj) + ecc(jj)
    END DO
  END DO
END DO
self%Bel=self%Bel/(4.d0*pi)
!
WRITE(*,*)'Building icoil->element magnetic reconstruction operator'
ALLOCATE(self%Bdr(bmesh%np,self%n_icoils,3))
self%Bdr=0.d0
!$omp parallel do private(ii,j,k,kk,pt_j,ecc,diffvec,cvec,cpt,pot_tmp,pot_last)
DO i=1,bmesh%np
  pt_j=bmesh%r(:,i)
  !---Compute driver contributions
  DO j=1,self%n_icoils
    ecc = 0.d0
    IF(self%icoils(j)%sens_mask)CYCLE
    DO k=1,self%icoils(j)%ncoils
      diffvec=0.d0
      DO kk=2,self%icoils(j)%coils(k)%npts
        cvec = self%icoils(j)%coils(k)%pts(:,kk)-self%icoils(j)%coils(k)%pts(:,kk-1)
        cpt = (self%icoils(j)%coils(k)%pts(:,kk)+self%icoils(j)%coils(k)%pts(:,kk-1))/2.d0
        diffvec = diffvec + cross_product(cvec,pt_j-cpt)/SUM((pt_j-cpt)**2)**1.5d0
      END DO
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
IF(TRIM(save_file)/='none')THEN
  hash_tmp(1) = self%nelems
  hash_tmp(2) = self%mesh%nc
  hash_tmp(3) = oft_simple_hash(C_LOC(self%mesh%lc),INT(4*3*self%mesh%nc,8))
  hash_tmp(4) = oft_simple_hash(C_LOC(self%mesh%r),INT(8*3*self%mesh%np,8))
  WRITE(*,*)'  Saving B-field operator to file: ',TRIM(save_file)
  CALL hdf5_create_file(TRIM(save_file))
  CALL hdf5_write(hash_tmp,TRIM(save_file),'MODEL_hash')
  IF(exists.AND.ALL(file_counts==hash_tmp))THEN
    CALL hdf5_write(self%Bel(:,:,1),TRIM(save_file),'Bel_X')
    CALL hdf5_write(self%Bel(:,:,2),TRIM(save_file),'Bel_Y')
    CALL hdf5_write(self%Bel(:,:,3),TRIM(save_file),'Bel_Z')
    IF(exists)THEN
      CALL hdf5_write(self%Bdr(:,:,1),TRIM(save_file),'Bdr_X')
      CALL hdf5_write(self%Bdr(:,:,2),TRIM(save_file),'Bdr_Y')
      CALL hdf5_write(self%Bdr(:,:,3),TRIM(save_file),'Bdr_Z')
    END IF
  END IF
END IF
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
IF(.NOT.all_oriented)CALL oft_abort("Error orienting cells","",__FILE__)
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
!> Read coil sets for "oft_in.xml" file
!------------------------------------------------------------------------------
subroutine tw_load_coils(group_node,ncoils,coils)
TYPE(xml_node), POINTER, INTENT(IN) :: group_node !< XML node relative to base `<thincurr>` node
INTEGER(4), INTENT(out) :: ncoils !< Number of coil sets found
TYPE(tw_coil_set), POINTER, INTENT(out) :: coils(:) !< List of coil sets
!---XML solver fields
integer(4) :: ncoil_sets,nread,coil_type
TYPE(xml_node), POINTER :: doc,coil_set,coil,thincurr_group
TYPE(xml_nodelist) :: coil_sets,coil_list
!---
LOGICAL :: success
INTEGER(4) :: i,j,k,io_unit,ierr,id,cell,ipath,ndims
INTEGER(4), ALLOCATABLE :: dim_sizes(:)
REAL(8) :: pts_tmp(2),res_per_len,radius,dl,theta
CHARACTER(LEN=10) :: coil_ind
CHARACTER(LEN=OFT_PATH_SLEN) :: coil_path
TYPE(tw_coil_set), POINTER :: coil_tmp
!---Count coil sets
ncoils=0
CALL xml_get_element(group_node,"coil_set",coil_sets,ierr)
IF(ierr==0)ncoils=coil_sets%n
ALLOCATE(coils(ncoils))
IF(ncoils==0)RETURN
!---Setup coil sets
DO i=1,ncoils
  coil_tmp=>coils(i)
  coil_set=>coil_sets%nodes(i)%this
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
    CALL xml_extractDataAttribute(coil_set,"name",coil_tmp%name,num=nread,iostat=ierr)
  ELSE
    WRITE(coil_tmp%name,'(A8,I5.5)')'UNKNOWN_',i
  END IF
  !---Get coil set resistivity per unit length (can be overriden)
  IF(xml_hasAttribute(coil_set,"res_per_len"))THEN
    CALL xml_extractDataAttribute(coil_set,"res_per_len",res_per_len,num=nread,iostat=ierr)
    coil_tmp%res_per_len=res_per_len
  END IF
  !---Get coil set radius (can be overriden)
  IF(xml_hasAttribute(coil_set,"radius"))THEN
    CALL xml_extractDataAttribute(coil_set,"radius",radius,num=nread,iostat=ierr)
    coil_tmp%radius=radius
  END IF
  !---Get sensor flag
  IF(xml_hasAttribute(coil_set,"sens_mask"))THEN
    CALL xml_extractDataAttribute(coil_set,"sens_mask",coil_tmp%sens_mask,num=nread,iostat=ierr)
    IF(coil_tmp%sens_mask)WRITE(*,'(2A,I4,A)')oft_indent,'Masking coil ',i,' from sensors'
  END IF
  ALLOCATE(coil_tmp%coils(coil_tmp%ncoils))
  DO j=1,coil_tmp%ncoils
    coil=>coil_list%nodes(j)%this
    !---Look for HDF5 path
    IF(xml_hasAttribute(coil,"path"))THEN
      CALL xml_extractDataAttribute(coil,"path",coil_path,num=nread,iostat=ierr)
      IF(ierr/=0)THEN
        WRITE(coil_ind,'(I4,2X,I4)')i,j
        CALL oft_abort('Error reading "path" in coil '//coil_ind,'tw_load_coils',__FILE__)
      END IF
      ipath=INDEX(coil_path,":")
      IF(ipath==0)THEN
        WRITE(coil_ind,'(I4,2X,I4)')i,j
        CALL oft_abort('Misformatted "path" attribute in coil '//coil_ind,'tw_load_coils',__FILE__)
      END IF
      CALL hdf5_field_get_sizes(coil_path(1:ipath-1),coil_path(ipath+1:OFT_PATH_SLEN),ndims,dim_sizes)
      IF(ndims<0)THEN
        WRITE(coil_ind,'(I4,2X,I4)')i,j
        CALL oft_abort('Failed to read HDF5 data sizes for coil '//coil_ind,'tw_load_coils',__FILE__)
      END IF
      IF(dim_sizes(1)/=3)THEN
        WRITE(coil_ind,'(I4,2X,I4)')i,j
        CALL oft_abort('Incorrect first dimension of HDF5 dataset for coil '//coil_ind,'tw_load_coils',__FILE__)
      END IF
      coil_tmp%coils(j)%npts=dim_sizes(2)
      DEALLOCATE(dim_sizes)
      ALLOCATE(coil_tmp%coils(j)%pts(3,coil_tmp%coils(j)%npts))
      CALL hdf5_read(coil_tmp%coils(j)%pts,coil_path(1:ipath-1),coil_path(ipath+1:OFT_PATH_SLEN),success=success)
      IF(.NOT.success)THEN
        WRITE(coil_ind,'(I4,2X,I4)')i,j
        CALL oft_abort('Failed to read HDF5 data for coil '//coil_ind,'tw_load_coils',__FILE__)
      END IF
    ELSE
      !---Read number of points
      IF(xml_hasAttribute(coil,"npts"))THEN
        CALL xml_extractDataAttribute(coil,"npts",coil_tmp%coils(j)%npts,num=nread,iostat=ierr)
        coil_type=2
      ELSE
        coil_type=1
      END IF
      SELECT CASE(coil_type)
        CASE(1)
          CALL xml_extractDataContent(coil,pts_tmp,num=nread,iostat=ierr)
          IF(ierr/=0)THEN
            WRITE(coil_ind,'(I4,2X,I4)')i,j
            CALL oft_abort('Error reading circular coil '//coil_ind,'tw_load_coils',__FILE__)
          END IF
          IF(coil_tmp%coils(j)%npts==0)coil_tmp%coils(j)%npts=181
          ALLOCATE(coil_tmp%coils(j)%pts(3,coil_tmp%coils(j)%npts))
          DO k=1,coil_tmp%coils(j)%npts
            theta=(k-1)*2.d0*pi/REAL(coil_tmp%coils(j)%npts-1,8)
            coil_tmp%coils(j)%pts(:,k)=[pts_tmp(1)*COS(theta),pts_tmp(1)*SIN(theta),pts_tmp(2)]
          END DO
        CASE(2)
          ALLOCATE(coil_tmp%coils(j)%pts(3,coil_tmp%coils(j)%npts))
          CALL xml_extractDataContent(coil,coil_tmp%coils(j)%pts,num=nread,iostat=ierr)
          IF(ierr/=0)THEN
            WRITE(coil_ind,'(I4,2X,I4)')i,j
            CALL oft_abort('Error reading coil '//coil_ind,'tw_load_coils',__FILE__)
          END IF
      END SELECT
    END IF
    !---Get scale factor
    IF(xml_hasAttribute(coil,"scale"))CALL xml_extractDataAttribute(coil,"scale",coil_tmp%scales(j),num=nread,iostat=ierr)
    !---Get coil resistivity per unit length
    IF(xml_hasAttribute(coil,"res_per_len"))CALL xml_extractDataAttribute(coil,"res_per_len",coil_tmp%res_per_len(j),num=nread,iostat=ierr)
    !---Get coil radius
    IF(xml_hasAttribute(coil,"radius"))CALL xml_extractDataAttribute(coil,"radius",coil_tmp%radius(j),num=nread,iostat=ierr)
  END DO
  IF(ASSOCIATED(coil_list%nodes))DEALLOCATE(coil_list%nodes)
END DO
IF(ASSOCIATED(coil_sets%nodes))DEALLOCATE(coil_sets%nodes)
!---
IF(oft_debug_print(1))THEN
  WRITE(*,*)
  WRITE(*,*)'Coils set definitions'
  WRITE(*,*)'========================='
  WRITE(*,*)'Found ',INT(ncoils,2),' coil sets'
  DO i=1,ncoils
    WRITE(*,*)'  Set  : ',i
    WRITE(*,*)'    nCoils  : ',coils(i)%ncoils
    ! DO j=1,coils(i)%ncoils
    !   WRITE(*,*)'    Position, Scale  : ',REAL(coils(i)%pt(:,j),4),REAL(coils(i)%scales(j),4)
    ! END DO
  END DO
  WRITE(*,*)'========================='
  WRITE(*,*)
END IF
end subroutine tw_load_coils
!------------------------------------------------------------------------------
!> Create a copy (by reference) of a coil set
!------------------------------------------------------------------------------
SUBROUTINE tw_copy_coil(coil_in,coil_out)
TYPE(tw_coil_set), INTENT(in) :: coil_in !< Source coil set
TYPE(tw_coil_set), INTENT(inout) :: coil_out !< Copy of source
! coil_out%type=coil_in%type
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
!---Load flux loops
sensors%nfloops=0
INQUIRE(FILE=TRIM(filename), EXIST=exists)
IF(exists)THEN
  WRITE(*,*)
  WRITE(*,*)'Loading floop information:'
  OPEN(NEWUNIT=io_unit, FILE=TRIM(filename))
  ierr=skip_comment_lines(io_unit)
  READ(io_unit,*)sensors%nfloops
  WRITE(*,*)'  # of floops =',sensors%nfloops
  ALLOCATE(sensors%floops(sensors%nfloops))
  DO i=1,sensors%nfloops
    READ(io_unit,*)
    READ(io_unit,*)sensors%floops(i)%np,sensors%floops(i)%scale_fac,sensors%floops(i)%name
    IF(oft_debug_print(1))WRITE(*,*)'    # of floops pts =',sensors%floops(i)%np
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
  WRITE(*,*)
  WRITE(*,*)'Setting jumper information:'
  OPEN(NEWUNIT=io_unit, FILE='jumpers_orient.dat')
  sensors%njumpers=SIZE(self%jumper_nsets)
  ALLOCATE(sensors%jumpers(sensors%njumpers))
  ALLOCATE(hole_facs(self%nholes))
  DO i=1,sensors%njumpers
    WRITE(sensors%jumpers(i)%name,'(A,I4.4)')'JUMPER_',i
    sensors%jumpers(i)%np=self%jumper_nsets(i)%n
    ALLOCATE(sensors%jumpers(i)%points(sensors%jumpers(i)%np))
    CALL order_jumper_list(self%jumper_nsets(i)%v,sensors%jumpers(i)%points,self%jumper_nsets(i)%n)
    ! Setup hole couplings
    hole_facs=0.d0
    DO j=1,sensors%jumpers(i)%np-1
      id=mesh_local_findedge(self%mesh,[sensors%jumpers(i)%points(j),sensors%jumpers(i)%points(j+1)])
      IF(id==0)CALL oft_abort("No matching edge for jumper points","tw_load_sensors",__FILE__)
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
    WRITE(io_unit,'(3Es24.15)')norm
    WRITE(*,'(I8,3Es24.15)')i,norm
  END DO
  DEALLOCATE(hole_facs)
  CLOSE(io_unit)
END IF
CONTAINS
!---------------------------------------------------------------------------------
!> Order jumper list into sequential chain
!---------------------------------------------------------------------------------
SUBROUTINE order_jumper_list(list_in,list_out,n)
INTEGER(4), INTENT(in) :: list_in(n) !< Input vertex list
INTEGER(4), INTENT(out) :: list_out(n) !< Reordered vertex list
INTEGER(4), INTENT(in) :: n !< Number of vertices in jumper definition
!---
INTEGER(4) :: ii,jj,k,ipt,eprev,ed,ptp
INTEGER(4), ALLOCATABLE :: lloop_tmp(:)
!---Find loop points and edges
ALLOCATE(lloop_tmp(n))
lloop_tmp=list_in
CALL sort_array(lloop_tmp, n)
!
DO jj=1,n
  IF(self%mesh%bp(lloop_tmp(jj)))THEN
    ipt=lloop_tmp(jj)
    EXIT
  END IF
  ! ipt=lloop_tmp(1)
END DO
k=1
list_out(k)=ipt
eprev=0
DO jj=1,n
  DO ii=self%mesh%kpe(ipt),self%mesh%kpe(ipt+1)-1
    ed = self%mesh%lpe(ii)
    IF(ed==eprev)CYCLE
    ptp = SUM(self%mesh%le(:,ed))-ipt
    IF(search_array(ptp, lloop_tmp, n)==0)CYCLE
    ipt=ptp
    IF(jj==n)EXIT
    k=k+1
    list_out(k)=ipt
    eprev=ed
    EXIT
  END DO
END DO
DEALLOCATE(lloop_tmp)
END SUBROUTINE order_jumper_list
end subroutine tw_load_sensors
!------------------------------------------------------------------------------
!> Load resistivity and sensor mask from "oft_in.xml" file
!------------------------------------------------------------------------------
subroutine tw_load_eta(self)
TYPE(tw_type), INTENT(inout) :: self !< Thin-wall model object
!---XML solver fields
integer(4) :: nshells,nreg_mesh,nread
TYPE(xml_node), POINTER :: sens_node,eta_group,thincurr_group
!---
INTEGER(4) :: i,j,io_unit,ierr,id,cell
REAL(8) :: location(2)
nreg_mesh=MAXVAL(self%mesh%reg)
ALLOCATE(self%Eta_reg(nreg_mesh))
self%Eta_reg=1.d0
ALLOCATE(self%sens_mask(nreg_mesh))
self%sens_mask=.FALSE.
IF(.NOT.ASSOCIATED(self%xml))THEN
  CALL oft_warn('No "thincurr" XML node, using "eta=mu0" for all regions')
  RETURN
END IF
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Loading region resistivity:'
!
CALL xml_get_element(self%xml,"eta",eta_group,ierr)
CALL xml_extractDataContent(eta_group,self%Eta_reg,num=nread,iostat=ierr)
IF(nread/=nreg_mesh)CALL oft_abort('Eta size mismatch','tw_load_eta',__FILE__)
! WRITE(*,'(2A)')oft_indent,'  Eta = ',REAL(self%Eta_reg,4)
DO i=1,nreg_mesh
  WRITE(*,'(A,I4,ES12.4)')oft_indent,i,self%Eta_reg(i)
  self%Eta_reg(i)=self%Eta_reg(i)/mu0 ! Convert to magnetic units
END DO
! Load sensor mask
CALL xml_get_element(self%xml,"sens_mask",sens_node,ierr)
IF(ierr==0)THEN
  WRITE(*,'(2A)')oft_indent,'Loading sensor mask:'
  CALL xml_extractDataContent(sens_node,self%sens_mask,num=nread,iostat=ierr)
  IF(nread/=nreg_mesh)CALL oft_abort('Sensor mask size mismatch','tw_load_eta',__FILE__)
  DO i=1,nreg_mesh
    WRITE(*,'(A,I4,L)')oft_indent,i,self%sens_mask(i)
  END DO
  ! WRITE(*,*)'  Sens mask = ',self%sens_mask
END IF
end subroutine tw_load_eta
!------------------------------------------------------------------------------
!> Load forcing mesh and fields for an MHD-style mode 
!------------------------------------------------------------------------------
subroutine tw_load_mode(filename,self,driver)
CHARACTER(LEN=*) :: filename !< Filename containing mode definition
TYPE(tw_type), INTENT(inout) :: self !< Thin-wall model of mode
REAL(8), POINTER, INTENT(out) :: driver(:,:) !< Sin/Cos pair corresponding to mode
!---
INTEGER(4) :: i,j,k,l,face,io_unit,ierr,ndim
INTEGER(4), ALLOCATABLE :: kfh_tmp(:),dim_sizes(:)
REAL(8) :: rgop(3,3),norm_i(3),f(3),area_i
WRITE(*,*)
WRITE(*,'(2A)')oft_indent,'Loading mode forcing mesh:'
WRITE(*,'(3A)')oft_indent,'  filename    = ',TRIM(filename)
!---Save to file
ALLOCATE(oft_trimesh::self%mesh)
CALL self%mesh%setup(-1,.FALSE.)
!
CALL hdf5_field_get_sizes('tCurr_mode_model.h5','mesh/R',ndim,dim_sizes)
self%mesh%np=dim_sizes(2)
DEALLOCATE(dim_sizes)
ALLOCATE(self%mesh%r(3,self%mesh%np))
CALL hdf5_read(self%mesh%r,'tCurr_mode_model.h5','mesh/R')
!
CALL hdf5_field_get_sizes('tCurr_mode_model.h5','mesh/LC',ndim,dim_sizes)
self%mesh%nc=dim_sizes(2)
DEALLOCATE(dim_sizes)
ALLOCATE(self%mesh%lc(3,self%mesh%nc))
CALL hdf5_read(self%mesh%lc,'tCurr_mode_model.h5','mesh/LC')
!
self%nclosures=1
ALLOCATE(self%closures(self%nclosures))
CALL hdf5_read(self%closures(1),'tCurr_mode_model.h5','mesh/SIDESET0001')
! OPEN(NEWUNIT=io_unit, FILE=TRIM(filename))
! READ(io_unit,*)self%mesh%np
! ALLOCATE(self%mesh%r(3,self%mesh%np))
! DO i=1,self%mesh%np
!   READ(io_unit,*)self%mesh%r(:,i)
! END DO
! self%nclosures=1
! ALLOCATE(self%closures(self%nclosures))
! READ(io_unit,*)self%closures(1)
! !
! READ(io_unit,*)self%mesh%nc
! ALLOCATE(self%mesh%lc(3,self%mesh%nc),self%mesh%reg(self%mesh%nc))
! self%mesh%reg=1
! DO i=1,self%mesh%nc
!   READ(io_unit,*)self%mesh%lc(:,i)
! END DO
CALL bmesh_local_init(self%mesh)
!---
self%nholes=2
! READ(io_unit,*)self%nholes
ALLOCATE(self%hmesh(self%nholes))
CALL hdf5_field_get_sizes('tCurr_mode_model.h5','mesh/NODESET0001',ndim,dim_sizes)
self%hmesh(1)%n=dim_sizes(1)
DEALLOCATE(dim_sizes)
ALLOCATE(self%hmesh(1)%lp(self%hmesh(1)%n))
CALL hdf5_read(self%hmesh(1)%lp,'tCurr_mode_model.h5','mesh/SIDESET0001')
CALL hdf5_field_get_sizes('tCurr_mode_model.h5','mesh/NODESET0002',ndim,dim_sizes)
self%hmesh(2)%n=dim_sizes(1)
DEALLOCATE(dim_sizes)
ALLOCATE(self%hmesh(2)%lp(self%hmesh(2)%n))
CALL hdf5_read(self%hmesh(2)%lp,'tCurr_mode_model.h5','mesh/SIDESET0002')
ALLOCATE(self%kfh(self%mesh%nc+1))
self%kfh=0
DO i=1,self%nholes
  ! READ(io_unit,*)self%hmesh(i)%n
  ! ALLOCATE(self%hmesh(i)%lp(self%hmesh(i)%n))
  ! DO j=1,self%hmesh(i)%n
  !   READ(io_unit,*)self%hmesh(i)%lp(j)
  ! END DO
  CALL tw_setup_hole(self%mesh, self%hmesh(i))
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
self%nelems = self%mesh%np + self%nholes - 1
!
! READ(io_unit,*)
ALLOCATE(driver(self%nelems,2))
CALL hdf5_read(driver,'tCurr_mode_model.h5','thincurr/driver')
! DO i=1,self%nelems
!   READ(io_unit,*)driver(i,:)
! END DO
! CLOSE(io_unit)
WRITE(*,'(2A,I8)')oft_indent,'  # of points = ',self%mesh%np
WRITE(*,'(2A,I8)')oft_indent,'  # of edges  = ',self%mesh%ne
WRITE(*,'(2A,I8)')oft_indent,'  # of faces  = ',self%mesh%nc
WRITE(*,'(2A,I8)')oft_indent,'  # of holes  = ',self%nholes
!
ALLOCATE(self%pmap(self%mesh%np))
self%pmap=0
self%np_active=0
DO i=1,self%mesh%np
  IF(self%mesh%bp(i))CYCLE
  self%np_active=self%np_active+1
  self%pmap(i)=self%np_active
END DO
! Mark closure vertex for removal
DO i=1,self%nclosures
  self%pmap(self%closures(i))=-i
END DO
self%nclosures=-self%nclosures
! Reindex, removing closure vertex
self%np_active=0
DO i=1,self%mesh%np
  IF(self%pmap(i)>0)THEN
    self%np_active=self%np_active+1
    self%pmap(i)=self%np_active
  ELSE
    self%pmap(i)=0
  END IF
END DO
!
ALLOCATE(self%kpmap_inv(self%np_active+1),self%lpmap_inv(self%np_active))
self%kpmap_inv=[(i,i=1,self%np_active+1)]
DO i=1,self%mesh%np
  IF(self%pmap(i)/=0)self%lpmap_inv(self%pmap(i))=i
END DO
self%nelems = self%np_active + self%nholes
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
ALLOCATE(self%sens_mask(1))
self%sens_mask=.FALSE.
end subroutine tw_load_mode
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE tw_build_boozer(self,s,alpha)
CLASS(tw_plasma_boozer), INTENT(INOUT) :: self !< Thin-wall model for structures
REAL(8), INTENT(in) :: s,alpha
!
INTEGER(4) :: i,j,k,l,info,nwall,nplasma,ncoils
COMPLEX(8) :: tmp
REAL(8), CONTIGUOUS, POINTER, DIMENSION(:,:) :: M_wp,M_pw
COMPLEX(8), CONTIGUOUS, POINTER, DIMENSION(:,:) :: tmp_pc,tmp_pw
!---Build helper matrices
nwall=self%wall%nelems-self%wall%n_vcoils
nplasma=self%plasma%nelems-self%plasma%n_vcoils
ncoils=self%plasma%n_vcoils
NULLIFY(M_wp,M_pw)
CALL tw_compute_LmatDirect(self%plasma,M_wp,col_model=self%wall)
ALLOCATE(M_pw(nplasma,nwall))
M_pw=TRANSPOSE(M_wp)
! !---Build rho matrix
! ALLOCATE(self%rho(nplasma,nplasma))
! self%rho=self%plasma%Lmat(1:nplasma,1:nplasma)*(1.d0,0.d0)
! CALL lapack_matinv(nplasma,self%rho,info)
! DO i=1,nplasma
!   DO j=1,nplasma
!     ! self%rho(i,j)=self%rho(i,j)*(-1.d0/(s,alpha)-1.d0)
!   END DO
! END DO
! !---Build l_w = L_w + M_wp * rho * M_pw
! ALLOCATE(self%L_w(nwall,nwall))
! self%L_w=self%wall%Lmat(1:nwall,1:nwall) &
!   + MATMUL(M_wp,MATMUL(rho,M_pw))

! !---Build l_wc = M_wc + M_wp * rho * M_pc
! ALLOCATE(self%L_wc(nwall,self%wall%n_vcoils))
! self%L_wc=self%wall%Lmat(1:nwall,nwall+1:nwall+ncoils) &
!   + MATMUL(M_wp,MATMUL(rho,M_pc))

! !---Build l_wd = M_wp + M_wp * rho * L_p
! ALLOCATE(self%L_wd(nwall,nplasma))
! self%L_wd=M_wp

! !---Build l_cw = M_cw + M_cp * rho * M_pw
! ALLOCATE(self%L_cw(ncoils,nwall))
! self%L_cw=self%wall%Lmat(nwall+1:nwall+ncoils,1:nwall)

! !---Build l_c = L_c + M_cp * rho * M_pc
! ALLOCATE(self%L_c(ncoils,ncoils))
! self%L_c=self%wall%Lmat(nwall+1:nwall+ncoils,nwall+1:nwall+ncoils)

! !---Build l_cd = M_cp + M_cp * rho * L_p
! ALLOCATE(self%L_cd(ncoils,nplasma))
! self%L_cd=self%plasma%Lmat(nplasma+1:nplasma+ncoils,1:nplasma)

! !---Build l_dw = M_pw + L_p * rho * M_pw
! ALLOCATE(self%L_dw(nplasma,nwall))
! self%L_dw=M_pw

! !---Build l_dc = M_pc + L_p * rho * M_pc
! ALLOCATE(self%L_dc(nplasma,ncoils))
! self%L_dc=self%plasma%Lmat(1:nplasma,nplasma+1:nplasma+ncoils)

! !---Build l_d = L_p + L_p * rho * L_p
! ALLOCATE(self%L_d(nplasma,nplasma))
! self%L_d=self%plasma%Lmat(1:nplasma,1:nplasma)

END SUBROUTINE tw_build_boozer
!---------------------------------------------------------------------------------
!> Save solution vector for thin-wall model for plotting in VisIt
!---------------------------------------------------------------------------------
SUBROUTINE tw_recon_curr(self,pot,curr)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
real(8), intent(in) :: pot(:) !< Solution values [self%nelems]
real(8), intent(out) :: curr(:,:) !< Solution values [3,self%mesh%nelems]
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
END SUBROUTINE tw_recon_curr
!---------------------------------------------------------------------------------
!> Save solution vector for thin-wall model for plotting in VisIt
!---------------------------------------------------------------------------------
SUBROUTINE tw_save_pfield(self,a,tag)
TYPE(tw_type), INTENT(in) :: self !< Thin-wall model object
real(8), intent(in) :: a(:) !< Solution values [self%nelems]
character(LEN=*), intent(in) :: tag !< Path to save vector in HDF5 plot files
INTEGER(4) :: i,j,k,jj,pt,ih,ihp,ihc
REAL(8) :: rcurr(3),ftmp(3),gop(3,3),area,norm(3)
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: ptvec,cellvec
DEBUG_STACK_PUSH
!---Avg to cells
ALLOCATE(cellvec(3,self%mesh%nc))
CALL tw_recon_curr(self,a,cellvec)
CALL self%mesh%save_cell_vector(cellvec/mu0,self%xdmf,TRIM(tag)) ! Convert back to Amps
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
END SUBROUTINE tw_save_pfield
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
IF(oft_env%head_proc.AND.oft_env%pm)WRITE(*,'(6A)')oft_indent,'Writting "',TRIM(path), &
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
  pcoil_vals=outvec%v(self%np_active+self%nholes:self%nelems)
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
END MODULE thin_wall
