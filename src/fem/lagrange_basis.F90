!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_lag_basis.F90
!
!> @defgroup doxy_oft_lag Lagrange
!! Lagrange finite element implementation for the Open FUSION Toolkit
!! @ingroup doxy_oft_fem
!
!> Base Lagrange FE class and basis evaluation
!! - FE Construction
!! - Basis evaluation
!!   - Interpolation
!!   - Gradient
!!
!! @authors Chris Hansen
!! @date August 2011
!! @ingroup doxy_oft_lag
!---------------------------------------------------------------------------
MODULE oft_lag_basis
USE, INTRINSIC :: iso_c_binding, only: c_int, c_double
USE oft_base
USE oft_lag_poly
USE oft_mesh_type, ONLY: oft_mesh, oft_bmesh
USE oft_mesh_local_util, ONLY: mesh_local_orient, oriented_cell, oriented_edges, &
  oriented_faces
USE oft_stitching, ONLY: oft_global_stitch
USE oft_trimesh_type, ONLY: tri_2d_grid
USE oft_quadmesh_type, ONLY: quad_2d_grid, quad_grid_orient
USE oft_tetmesh_type, ONLY: tet_3d_grid
USE oft_hexmesh_type, ONLY: hex_3d_grid, hex_grid_forient
USE multigrid, ONLY: multigrid_mesh, multigrid_level
USE oft_la_base, ONLY: oft_matrix, oft_graph
USE fem_base, ONLY: oft_afem_type, oft_fem_type, oft_ml_fem_type, oft_bfem_type, &
  fem_delete, bfem_delete
USE fem_composite, ONLY: oft_fem_comp_type, oft_ml_fem_comp_type
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------
!> Lagrange operator container
!---------------------------------------------------------------------------
type, extends(oft_fem_type) :: oft_scalar_fem
  INTEGER(i4) :: inodesp(3,8) !< Needs docs
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:,:,:) :: inodese => NULL() !< Needs docs
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:,:,:) :: inodesf => NULL() !< Needs docs
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:,:) :: inodesc => NULL() !< Needs docs
  REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:) :: xnodes => NULL() !< Needs docs
  REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:,:) :: sn => NULL() !< Surface normal vector
  ! REAL(r8), POINTER, DIMENSION(:,:) :: ct => NULL() !< Curve tangent vector
  TYPE(oft_graph), POINTER :: interp_graph => NULL() !< Interpolation graph
  CLASS(oft_matrix), POINTER :: interp => NULL() !< Interpolation matrix
  CLASS(oft_matrix), POINTER :: vinterp => NULL() !< Vector interpolation matrix
contains
  !> Destory FE type
  PROCEDURE :: delete => scalar_fem_delete
end type oft_scalar_fem
!---------------------------------------------------------------------------
!> Lagrange operator container
!---------------------------------------------------------------------------
type, extends(oft_bfem_type) :: oft_scalar_bfem
  INTEGER(i4) :: inodesp(2,4) !< Needs docs
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:,:,:) :: inodese => NULL() !< Needs docs
  INTEGER(i4), CONTIGUOUS, POINTER, DIMENSION(:,:) :: inodesf => NULL() !< Needs docs
  REAL(r8), CONTIGUOUS, POINTER, DIMENSION(:) :: xnodes => NULL() !< Needs docs
contains
  !> Destory FE type
  PROCEDURE :: delete => scalar_bfem_delete
end type oft_scalar_bfem
!---Global Variables
integer(i4), parameter :: oft_lagrange_id = 1 !< FE type ID
contains
!------------------------------------------------------------------------------
!> Cast abstract FE type to 3D lagrange finite element type
!!
!! The source matrix must be @ref oft_scalar_fem or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION oft_3D_lagrange_cast(self,source) RESULT(success)
CLASS(oft_scalar_fem), POINTER, INTENT(out) :: self !< Reference to source object with desired class
CLASS(oft_afem_type), TARGET, INTENT(in) :: source !< Source object to reference
LOGICAL :: success !< Cast success flag
DEBUG_STACK_PUSH
SELECT TYPE(source)
  CLASS IS(oft_scalar_fem)
    self=>source
    success=.TRUE.
  CLASS DEFAULT
    NULLIFY(self)
    success=.FALSE.
END SELECT
DEBUG_STACK_POP
END FUNCTION oft_3D_lagrange_cast
!------------------------------------------------------------------------------
!> Cast abstract FE type to 2D lagrange finite element type
!!
!! The source matrix must be @ref oft_scalar_bfem or a child class, otherwise
!! pointer will be returned as `null` and `success == .FALSE.`
!------------------------------------------------------------------------------
FUNCTION oft_2D_lagrange_cast(self,source) RESULT(success)
CLASS(oft_scalar_bfem), POINTER, INTENT(out) :: self !< Reference to source object with desired class
CLASS(oft_afem_type), TARGET, INTENT(in) :: source !< Source object to reference
LOGICAL :: success !< Cast success flag
DEBUG_STACK_PUSH
SELECT TYPE(source)
  CLASS IS(oft_scalar_bfem)
    self=>source
    success=.TRUE.
  CLASS DEFAULT
    NULLIFY(self)
    success=.FALSE.
END SELECT
DEBUG_STACK_POP
END FUNCTION oft_2D_lagrange_cast
!---------------------------------------------------------------------------
!> Compute surface normals for use in boundary conditions
!---------------------------------------------------------------------------
SUBROUTINE oft_lag_boundary(lag_rep)
class(oft_scalar_fem), intent(inout) :: lag_rep
INTEGER(i4), PARAMETER :: i2(3)=(/2,3,1/),i3(3)=(/3,1,2/)
INTEGER(i4), PARAMETER :: fc1(4)=(/2,3,1,1/),fc2(4)=(/3,1,2,2/),fc3(4)=(/4,4,4,3/)
INTEGER(i4), PARAMETER :: ed1(6)=(/1,2,3,2,3,1/),ed2(6)=(/4,4,4,3,1,2/)
INTEGER(i4) :: i,k,l,m,el,jc
INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: j,emap
REAL(r8) :: v,goptmp(3,4),norm(3)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: pc
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: f,ct
REAL(r8), POINTER, DIMENSION(:,:) :: sn
class(oft_mesh), pointer :: mesh
DEBUG_STACK_PUSH
mesh=>lag_rep%mesh
!---
ALLOCATE(lag_rep%sn(3,lag_rep%ne),lag_rep%bc(lag_rep%ne))
sn=>lag_rep%sn
sn=0._r8
lag_rep%bc=0_i4
!---
ALLOCATE(f(4,lag_rep%nce))
!---
ALLOCATE(j(lag_rep%nce))
DO i=1,mesh%nbc
  !---Check if boundary cell
  l=mesh%lbc(i)
  IF(.NOT.mesh%global%gbc(l))CYCLE
  !---Get local to global DOF mapping
  call lag_rep%ncdofs(l,j)
  !---
  DO m=1,lag_rep%nce
    CALL oft_lag_npos(lag_rep,l,m,f(:,m))
    k=lag_rep%cmap(m)%el
    SELECT CASE(lag_rep%cmap(m)%type)
      CASE(2)
        el=ABS(mesh%lce(k,l))
        IF(.NOT.mesh%global%gbe(el))CYCLE
        !---
        IF(mesh%ce(el))THEN
          CALL mesh%ctang(l,k,f(:,m),norm)
          sn(:,j(m))=norm
          lag_rep%bc(j(m))=2._i4
        ELSE
          DO k=1,mesh%cell_nf
            ! IF(f(k,m)/=0._r8)CYCLE
            IF(.NOT.mesh%global%gbf(ABS(mesh%lcf(k,l))))CYCLE
            IF(ALL(lag_rep%cmap(m)%el/=ABS(mesh%cell_fe(:,k))))CYCLE
            CALL mesh%snormal(l,k,f(:,m),norm)
            sn(:,j(m))=sn(:,j(m))+norm
            lag_rep%bc(j(m))=3._i4
          END DO
        END IF
      CASE(3)
        el=mesh%lcf(k,l)
        IF(.NOT.mesh%global%gbf(el))CYCLE
        !---
        CALL mesh%snormal(l,k,f(:,m),norm)
        sn(:,j(m))=norm
        lag_rep%bc(j(m))=3._i4
      CASE DEFAULT
        CYCLE
      END SELECT
  END DO
END DO
!---
CALL oft_global_stitch(lag_rep%linkage,sn(1,:),1)
CALL oft_global_stitch(lag_rep%linkage,sn(2,:),1)
CALL oft_global_stitch(lag_rep%linkage,sn(3,:),1)
!---
DO i=1,lag_rep%ne
  v=SQRT(SUM(sn(:,i)**2))
  IF(v>1.E-10_r8)sn(:,i)=sn(:,i)/v
END DO
!---
DO i=1,mesh%nbe
  jc=mesh%lbe(i)
  IF(.NOT.mesh%global%gbe(jc))CYCLE
  l=mesh%lec(mesh%kec(jc))
  !---Get local to global DOF mapping
  call lag_rep%ncdofs(l,j)
  !---
  DO m=1,lag_rep%nce
    k=lag_rep%cmap(m)%el
    IF(lag_rep%cmap(m)%type/=2)CYCLE
    el=ABS(mesh%lce(k,l))
    IF(el/=jc)CYCLE
    !---
    IF(mesh%ce(el))THEN
      lag_rep%bc(j(m))=2._i4
    ELSE
      lag_rep%bc(j(m))=3._i4
    END IF
  END DO
END DO
!---
ALLOCATE(pc(mesh%np),emap(mesh%ne))
CALL get_inverse_map(mesh%lbe,mesh%nbe,emap,mesh%ne)
pc=0._r8
DO i=1,mesh%nbp
  l=mesh%lbp(i)
  IF(.NOT.mesh%global%gbp(l))CYCLE
  !---
  DO jc=mesh%kpe(l),mesh%kpe(l+1)-1
    k=mesh%lpe(jc)
    IF(.NOT.mesh%global%gbe(k))CYCLE
    IF(mesh%ce(k).AND.emap(k)>0)THEN
      ! IF(mesh%linkage%leo(emap(k)))pc(l)=pc(l)+1_i4
      IF(mesh%estitch%leo(emap(k)))pc(l)=pc(l)+1_i4
    END IF
  END DO
END DO
DEALLOCATE(emap)
CALL oft_global_stitch(mesh%pstitch,pc,1)
ALLOCATE(ct(3,mesh%np))
ct=0._r8
DO i=1,mesh%nbp
  l=mesh%lbp(i)
  IF(.NOT.mesh%global%gbp(l))CYCLE
  IF(pc(l)>2)then
    sn(:,l)=0._r8
    lag_rep%bc(l)=1._i4
    CYCLE
  END IF
  !---
  DO jc=mesh%kpc(l),mesh%kpc(l+1)-1
    k=mesh%lpc(jc)
    !IF(.NOT.mesh%global%gbc(k))CYCLE
    DO m=1,mesh%cell_np
      IF(mesh%lc(m,k)==l)EXIT
    END DO
    el=m
    ! f(:,1)=0._r8; f(el,1)=1._r8
    CALL mesh%vlog(el, f(:,1))
    !---
    IF(pc(l)==0)THEN
      lag_rep%bc(l)=3._i4
      IF(.NOT.mesh%global%gbc(k))CYCLE
      DO m=1,mesh%cell_nf
        IF(.NOT.mesh%global%gbf(mesh%lcf(m,k)))CYCLE
        ! IF(fc1(m)==el.OR.fc2(m)==el.OR.fc3(m)==el)THEN
        IF(ANY(mesh%cell_fc(:,m)==el))THEN
          CALL mesh%snormal(k,m,f(:,1),norm)
          sn(:,l)=sn(:,l)+norm
        END IF
      END DO
    ELSE IF(pc(l)>0)THEN
      lag_rep%bc(l)=2._i4
      IF(.NOT.mesh%global%gbc(k))CYCLE
      DO m=1,mesh%cell_ne
        IF(.NOT.mesh%global%gbe(ABS(mesh%lce(m,k))))CYCLE
        ! IF(ed1(m)==el.OR.ed2(m)==el)THEN
        IF(ANY(mesh%cell_ed(:,m)==el))THEN
          IF(.NOT.mesh%ce(ABS(mesh%lce(m,k))))CYCLE
          CALL mesh%ctang(k,m,f(:,1),norm)
          IF(ALL(norm==0._r8))CYCLE
          ct(:,l)=ct(:,l)+norm*SIGN(1._r8,DOT_PRODUCT(ct(:,l),norm)+1.E-10_r8)
        END IF
      END DO
    END IF
  END DO
END DO
!---
CALL oft_global_stitch(mesh%pstitch,sn(1,1:mesh%np),1)
CALL oft_global_stitch(mesh%pstitch,sn(2,1:mesh%np),1)
CALL oft_global_stitch(mesh%pstitch,sn(3,1:mesh%np),1)
!---
DO i=1,mesh%nbp
  l=mesh%lbp(i)
  IF(lag_rep%bc(l)/=2._i4)CYCLE
  sn(:,l)=ct(:,l)
  v=SQRT(SUM(sn(:,l)**2))
  IF(v>1.E-10_r8)sn(:,l)=sn(:,l)/v
END DO
!---
CALL oft_global_stitch(mesh%pstitch,ct(1,1:mesh%np),2)
CALL oft_global_stitch(mesh%pstitch,ct(2,1:mesh%np),2)
CALL oft_global_stitch(mesh%pstitch,ct(3,1:mesh%np),2)
!---
DO i=1,mesh%nbp
  l=mesh%lbp(i)
  IF(lag_rep%bc(l)/=2._i4)CYCLE
  sn(:,l)=sn(:,l)+ct(:,l)*SIGN(1._r8,DOT_PRODUCT(ct(:,l),sn(:,l))+1.E-10_r8)
  v=SQRT(SUM(sn(:,l)**2))
  IF(v>1.E-10_r8)sn(:,l)=sn(:,l)/v
END DO
!---
DO i=1,mesh%np
  v=SQRT(SUM(sn(:,i)**2))
  IF(v>1.E-10_r8)sn(:,i)=sn(:,i)/v
END DO
DEBUG_STACK_POP
END SUBROUTINE oft_lag_boundary
!---------------------------------------------------------------------------
!> Construct lagrange scalar FE on each mesh level
!!
!! @note Highest supported representation is quartic
!---------------------------------------------------------------------------
subroutine  oft_lag_setup(mg_mesh,order,ML_lag_obj,ML_blag_obj,ML_vlag_obj,minlev)
type(multigrid_mesh), target, intent(inout) :: mg_mesh
integer(i4), intent(in) :: order !< Order of representation desired
TYPE(oft_ml_fem_type), target, optional, INTENT(inout) :: ML_lag_obj
TYPE(oft_ml_fem_type), optional, INTENT(inout) :: ML_blag_obj
TYPE(oft_ml_fem_comp_type), optional, INTENT(inout) :: ML_vlag_obj
integer(i4), optional, intent(in) :: minlev !< Lowest level to construct
integer(i4) :: i,nlevels,minlev_out
DEBUG_STACK_PUSH
minlev_out=1
IF(PRESENT(minlev))minlev_out=minlev
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(2A)')oft_indent,'**** Creating Lagrange FE space'
  WRITE(*,'(A,2X,A,I4)')oft_indent,'Order  = ',order
  WRITE(*,'(A,2X,A,I4)')oft_indent,'Minlev = ',minlev_out
END IF
CALL oft_increase_indent
!---Allocate multigrid operators
nlevels=mg_mesh%mgdim+(order-1)
IF(minlev_out<0)minlev_out=nlevels
IF(PRESENT(ML_lag_obj))THEN
  IF(.NOT.ASSOCIATED(mg_mesh%meshes))CALL oft_abort("No volume mesh available for 3D elements","oft_lag_setup",__FILE__)
  ML_lag_obj%nlevels=nlevels
  ML_lag_obj%minlev=minlev_out
  ML_lag_obj%ml_mesh=>mg_mesh
ELSE
  IF(PRESENT(ML_vlag_obj))CALL oft_abort("Vector FE space requires scalar space","oft_lag_setup",__FILE__)
END IF
IF(PRESENT(ML_blag_obj))THEN
  IF(.NOT.ASSOCIATED(mg_mesh%smeshes))CALL oft_abort("No surface mesh available for 2D elements","oft_lag_setup",__FILE__)
  ML_blag_obj%ml_mesh=>mg_mesh
  ML_blag_obj%nlevels=nlevels
  ML_blag_obj%minlev=minlev_out
END IF
!---Set linear elements
do i=1,mg_mesh%mgdim-1
  IF(i<minlev_out)CYCLE
  CALL multigrid_level(mg_mesh,i)
  IF(PRESENT(ML_lag_obj))THEN
    CALL oft_lag_setup_vol(ML_lag_obj%levels(i)%fe,mg_mesh%mesh,1)
    IF(mg_mesh%level==mg_mesh%nbase)ML_lag_obj%blevel=i
    CALL ML_lag_obj%set_level(i)
    SELECT TYPE(this=>ML_lag_obj%current_level)
      CLASS IS(oft_scalar_fem)
        CALL oft_lag_boundary(this)
      CLASS DEFAULT
        CALL oft_abort("Error casting FE obj","oft_lag_setup",__FILE__)
    END SELECT
    IF(mg_mesh%level==mg_mesh%nbase)ML_lag_obj%blevel=i
  END IF
  IF(PRESENT(ML_blag_obj))THEN
    CALL oft_lag_setup_bmesh(ML_blag_obj%levels(i)%fe,mg_mesh%smesh,1)
    IF(mg_mesh%level==mg_mesh%nbase)ML_blag_obj%blevel=i
  END IF
end do
call multigrid_level(mg_mesh,mg_mesh%mgdim)
!---Set high order elements
do i=1,order
  IF(i>1.AND.mg_mesh%mgdim+i-1<minlev_out)CYCLE
  IF(PRESENT(ML_lag_obj))THEN
    CALL oft_lag_setup_vol(ML_lag_obj%levels(mg_mesh%mgdim+i-1)%fe,mg_mesh%mesh,i)
    CALL ML_lag_obj%set_level(mg_mesh%mgdim+i-1)
    SELECT TYPE(this=>ML_lag_obj%current_level)
      CLASS IS(oft_scalar_fem)
        CALL oft_lag_boundary(this)
      CLASS DEFAULT
        CALL oft_abort("Error casting FE obj","oft_lag_setup",__FILE__)
    END SELECT
  END IF
  IF(PRESENT(ML_blag_obj))THEN
    CALL oft_lag_setup_bmesh(ML_blag_obj%levels(mg_mesh%mgdim+i-1)%fe,mg_mesh%smesh,i)
  END IF
end do
!---Setup vector rep
IF(PRESENT(ML_vlag_obj))THEN
  ML_vlag_obj%nlevels=ML_lag_obj%nlevels
  ML_vlag_obj%nfields=3
  ALLOCATE(ML_vlag_obj%ml_fields(ML_vlag_obj%nfields))
  ALLOCATE(ML_vlag_obj%field_tags(ML_vlag_obj%nfields))
  ML_vlag_obj%ml_fields(1)%ml=>ML_lag_obj
  ML_vlag_obj%field_tags(1)='x'
  ML_vlag_obj%ml_fields(2)%ml=>ML_lag_obj
  ML_vlag_obj%field_tags(2)='y'
  ML_vlag_obj%ml_fields(3)%ml=>ML_lag_obj
  ML_vlag_obj%field_tags(3)='z'
  CALL ML_vlag_obj%setup
END IF
!---Set to highest level when done
IF(PRESENT(ML_lag_obj))CALL ML_lag_obj%set_level(ML_lag_obj%nlevels)
IF(PRESENT(ML_vlag_obj))CALL ML_vlag_obj%set_level(ML_lag_obj%nlevels)
IF(PRESENT(ML_blag_obj))CALL ML_blag_obj%set_level(ML_blag_obj%nlevels)
IF(oft_env%head_proc)WRITE(*,*)
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine oft_lag_setup
!---------------------------------------------------------------------------
!> Construct lagrange scalar FE for a given order
!!
!! @note Highest supported representation is quartic
!---------------------------------------------------------------------------
subroutine oft_lag_setup_vol(self,tmesh,order)
class(oft_afem_type), pointer, intent(out) :: self !< Needs docs
class(oft_mesh), target, intent(in) :: tmesh !< Needs docs
integer(i4), intent(in) :: order !< Order of representation desired
DEBUG_STACK_PUSH
IF(oft_debug_print(1))THEN
  WRITE(*,'(2A)')oft_indent,'Creating 3D Lagrange FE space'
  WRITE(*,'(A,2X,A,I4)')oft_indent,'Order  = ',order
END IF
CALL oft_increase_indent
!---
ALLOCATE(oft_scalar_fem::self)
SELECT TYPE(self)
CLASS IS(oft_scalar_fem)
  self%mesh=>tmesh
  self%order=order
  self%dim=1
  self%type=oft_lagrange_id
  IF(self%mesh%type==3)THEN
    CALL hex_3d_grid(self%order,self%xnodes,self%inodesp, &
      self%inodese,self%inodesf,self%inodesc)
    select case(self%order)
      case(1)
        self%gstruct=(/1,0,0,0/)
      case(2)
        self%gstruct=(/1,1,1,1/)
      case(3)
        self%gstruct=(/1,2,4,8/)
      case(4)
        self%gstruct=(/1,3,9,27/)
      case default
        call oft_abort('Invalid polynomial degree (npmax=4)','oft_lag_setup_vol',__FILE__)
    end select
  ELSE
    CALL tet_3d_grid(self%order,self%xnodes,self%inodesf,self%inodesc)
    select case(self%order)
      case(1)
        self%gstruct=(/1,0,0,0/)
      case(2)
        self%gstruct=(/1,1,0,0/)
      case(3)
        self%gstruct=(/1,2,1,0/)
      case(4)
        self%gstruct=(/1,3,3,1/)
      case default
        call oft_abort('Invalid polynomial degree (npmax=4)','oft_lag_setup_vol',__FILE__)
    end select
  END IF
CLASS DEFAULT
  CALL oft_abort("Error allocate Lagrange FE object","oft_lag_setup_vol",__FILE__)
END SELECT
call self%setup(self%order*2+1)
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine oft_lag_setup_vol
!---------------------------------------------------------------------------
!> Construct lagrange scalar FE for a given order
!!
!! @note Highest supported representation is quartic
!---------------------------------------------------------------------------
subroutine oft_lag_setup_bmesh(self,tmesh,order)
class(oft_afem_type), pointer, intent(out) :: self !< Needs docs
class(oft_bmesh), target, intent(in) :: tmesh !< Needs docs
integer(i4), intent(in) :: order !< Order of representation desired
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(tmesh%parent))THEN
  IF(oft_debug_print(1))THEN
    WRITE(*,'(2A)')oft_indent,'Creating 2D Lagrange FE space'
    WRITE(*,'(A,2X,A,I4)')oft_indent,'Order  = ',order
  END IF
END IF
CALL oft_increase_indent
!---
ALLOCATE(oft_scalar_bfem::self)
SELECT TYPE(self)
CLASS IS(oft_scalar_bfem)
  self%mesh=>tmesh
  self%order=order
  self%dim=1
  self%type=oft_lagrange_id
  IF(self%mesh%type==3)THEN
    CALL quad_2d_grid(order,self%xnodes,self%inodesp,self%inodese,self%inodesf)
    select case(self%order)
      case(1)
        self%gstruct=(/1,0,0/)
      case(2)
        self%gstruct=(/1,1,1/)
      case(3)
        self%gstruct=(/1,2,4/)
      case(4)
        self%gstruct=(/1,3,9/)
      case default
        call oft_abort('Invalid polynomial degree (npmax=4)','oft_lag_setup_bmesh',__FILE__)
    end select
  ELSE
    CALL tri_2d_grid(order,self%xnodes,self%inodesf)
    select case(self%order)
      case(1)
        self%gstruct=(/1,0,0/)
      case(2)
        self%gstruct=(/1,1,0/)
      case(3)
        self%gstruct=(/1,2,1/)
      case(4)
        self%gstruct=(/1,3,3/)
      case default
        call oft_abort('Invalid polynomial degree (npmax=4)','oft_lag_setup_bmesh',__FILE__)
    end select
  END IF
CLASS DEFAULT
  CALL oft_abort("Error allocate Lagrange FE object","oft_lag_setup_bmesh",__FILE__)
END SELECT
call self%setup(order*2+1)
CALL oft_decrease_indent
DEBUG_STACK_POP
end subroutine oft_lag_setup_bmesh
!---------------------------------------------------------------------------
!> Evaluate lagrange interpolation function
!!
!! @note Evaluation is performed in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_lag_eval(self,cell,dof,f,val)
class(oft_scalar_fem), intent(in) :: self !< Lagrange type for evaluation
integer(i4), intent(in) :: cell !< Cell for evaluation
integer(i4), intent(in) :: dof !< Element to evaluate
real(r8), intent(in) :: f(:) !< Position in cell in logical space
real(r8), intent(out) :: val !< Value of interpolation function (dof) at point (f)
integer(i4) :: ed,etmp(2),fc,ftmp(3),ind,finds(16)
DEBUG_STACK_PUSH
IF(self%mesh%type==3)THEN
  select case(self%cmap(dof)%type)
    case(1)
      val=lag_3d(self%inodesp(:,self%cmap(dof)%el),f, &
        self%xnodes,self%order+1)
    case(2)
      ind=self%cmap(dof)%ind
      IF(self%mesh%lce(self%cmap(dof)%el,cell)<0)ind=self%order-ind
      val=lag_3d(self%inodese(:,ind,self%cmap(dof)%el),f, &
        self%xnodes,self%order+1)
    case(3)
      CALL hex_grid_forient(self%mesh%lcfo(self%cmap(dof)%el,cell),self%order,finds)
      val=lag_3d(self%inodesf(:,finds(self%cmap(dof)%ind),self%cmap(dof)%el),f, &
        self%xnodes,self%order+1)
    case(4)
      val=lag_3d(self%inodesc(:,self%cmap(dof)%ind),f, &
        self%xnodes,self%order+1)
  end select
ELSE
  IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
  select case(self%cmap(dof)%type)
    case(1)
      val=lag_1d(self%order+1,f(self%cmap(dof)%el),self%xnodes, &
         self%order+1)
    case(2)
      etmp=oriented_edges(:,self%cmap(dof)%el)
      val=lag_1d_bary(self%cmap(dof)%ind,f(etmp),self%xnodes, &
        self%order+1)
    case(3)
      ftmp=oriented_faces(:,self%cmap(dof)%el)
      val=lag_2d_bary(self%inodesf(:,self%cmap(dof)%ind,1), &
        f(ftmp),self%xnodes,self%order+1)
    case(4)
      val=lag_3d_bary(self%inodesc(:,self%cmap(dof)%ind), &
        f,self%xnodes,self%order+1)
  end select
END IF
DEBUG_STACK_POP
end subroutine oft_lag_eval
!---------------------------------------------------------------------------
!> Evaluate lagrange interpolation function
!!
!! @note Evaluation is performed in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_blag_eval(self,face,dof,f,val)
class(oft_scalar_bfem), intent(in) :: self !< Lagrange type for evaluation
integer(i4), intent(in) :: face !< Cell for evaluation
integer(i4), intent(in) :: dof !< Element to evaluate
real(r8), intent(in) :: f(:) !< Position in cell in logical space
real(r8), intent(out) :: val !< Value of interpolation function (dof) at point (f)
integer(i4) :: ed,etmp(2),fc,ftmp(3),finds(16),ind
DEBUG_STACK_PUSH
IF(self%mesh%type==3)THEN
  select case(self%cmap(dof)%type)
    case(1)
      val=lag_2d(self%inodesp(:,self%cmap(dof)%el),f, &
        self%xnodes,self%order+1)
    case(2)
      ind=self%cmap(dof)%ind
      IF(self%mesh%lce(self%cmap(dof)%el,face)<0)ind=self%order-ind
      val=lag_2d(self%inodese(:,ind,self%cmap(dof)%el),f, &
        self%xnodes,self%order+1)
    case(3)
      CALL quad_grid_orient(self%mesh%lco(face),self%order,finds)
      val=lag_2d(self%inodesf(:,finds(self%cmap(dof)%ind)),f, &
        self%xnodes,self%order+1)
  end select
ELSE
  select case(self%cmap(dof)%type)
    case(1)
      val=lag_1d(self%order+1,f(self%cmap(dof)%el),self%xnodes, &
         self%order+1)
    case(2)
      ed=self%mesh%lce(self%cmap(dof)%el,face)
      etmp=(/self%mesh%cell_ed(1,self%cmap(dof)%el),self%mesh%cell_ed(2,self%cmap(dof)%el)/)
      call orient_list2(ed,etmp)
      val=lag_1d_bary(self%cmap(dof)%ind,f(etmp),self%xnodes, &
        self%order+1)
    case(3)
      fc=self%mesh%lco(face)
      ftmp=(/1,2,3/)
      call orient_listn(fc,ftmp,3_i4)
      val=lag_2d_bary(self%inodesf(:,self%cmap(dof)%ind), &
        f(ftmp),self%xnodes,self%order+1)
  end select
END IF
DEBUG_STACK_POP
end subroutine oft_blag_eval
!---------------------------------------------------------------------------
!> Evaluate all lagrange interpolation functions
!!
!! @note Evaluation is performed in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_lag_eval_all(self,cell,f,rop)
class(oft_scalar_fem), intent(in) :: self !< Lagrange type for evaluation
integer(i4), intent(in) :: cell !< Cell for evaluation
real(r8), intent(in) :: f(:) !< Position in cell in logical space
real(r8), contiguous, intent(out) :: rop(:) !< Value of interpolation functions at point (f) [ncdofs]
integer(i4) :: i,j,offset,ind,inc,etmp(2),ftmp(3)
integer(i4), ALLOCATABLE, DIMENSION(:) :: finds
real(r8) :: pnorm
real(r8), ALLOCATABLE, DIMENSION(:,:) :: grid_1d
DEBUG_STACK_PUSH
IF(self%mesh%type==3)THEN
  ALLOCATE(grid_1d(self%order+1,3))
  DO i=1,self%order+1
    grid_1d(i,1) = lag_1d(i,f(1),self%xnodes,self%order+1)
    grid_1d(i,2) = lag_1d(i,f(2),self%xnodes,self%order+1)
    grid_1d(i,3) = lag_1d(i,f(3),self%xnodes,self%order+1)
  END DO
  DO i=1,8
    rop(i)=grid_1d(self%inodesp(1,i),1)*grid_1d(self%inodesp(2,i),2)*grid_1d(self%inodesp(3,i),3)
  END DO
  IF(self%order>1)THEN
    !---Edges
    DO i=1,12
      offset=(i-1)*self%gstruct(2)+8
      ind=1; inc=1
      IF(self%mesh%lce(i,cell)<0)THEN
        ind=self%gstruct(2); inc=-1
      END IF
      DO j=1,self%gstruct(2)
        rop(offset+j)=grid_1d(self%inodese(1,ind,i),1) &
          *grid_1d(self%inodese(2,ind,i),2)*grid_1d(self%inodese(3,ind,i),3)
        ind=ind+inc
      END DO
    END DO
    !---Faces
    ALLOCATE(finds(self%gstruct(3)))
    DO i=1,6
      offset=(i-1)*self%gstruct(3)+12*self%gstruct(2)+8
      CALL hex_grid_forient(self%mesh%lcfo(i,cell),self%order,finds)
      DO j=1,self%gstruct(3)
        rop(offset+j)=grid_1d(self%inodesf(1,finds(j),i),1) &
          *grid_1d(self%inodesf(2,finds(j),i),2)*grid_1d(self%inodesf(3,finds(j),i),3)
      END DO
    END DO
    DEALLOCATE(finds)
    !---Cell
    offset=6*self%gstruct(3)+12*self%gstruct(2)+8
    DO j=1,self%gstruct(4)
      rop(offset+j)=grid_1d(self%inodesc(1,j),1) &
        *grid_1d(self%inodesc(2,j),2)*grid_1d(self%inodesc(3,j),3)
    END DO
  END IF
  DEALLOCATE(grid_1d)
ELSE
  SELECT CASE(self%order)
  CASE(1)
    DO i=1,4
      rop(i)=f(i)
    END DO
  CASE(2)
    CALL tet_eval_all2()!self,cell,f,rop)
  CASE(3)
    CALL tet_eval_all3()!self,cell,f,rop)
  CASE(4)
    CALL tet_eval_all4()!self,cell,f,rop)
  CASE DEFAULT
    IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
    DO i=1,4
      rop(i)=lag_1d(self%order+1,f(i),self%xnodes,self%order+1)
    END DO
    IF(self%order>1)THEN
      !---Edges
      DO i=1,6
        offset=(i-1)*self%gstruct(2)+4
        etmp=oriented_edges(:,i)
        DO j=1,self%gstruct(2)
          rop(offset+j)=lag_1d_bary(j,f(etmp),self%xnodes,self%order+1)
        END DO
      END DO
      !---Faces
      DO i=1,4
        offset=(i-1)*self%gstruct(3)+6*self%gstruct(2)+4
        ftmp=oriented_faces(:,i)
        DO j=1,self%gstruct(3)
          rop(offset+j)=lag_2d_bary(self%inodesf(:,j,1),f(ftmp),self%xnodes,self%order+1)
        END DO
      END DO
      !---Cell
      offset=4*self%gstruct(3)+6*self%gstruct(2)+4
      DO j=1,self%gstruct(4)
        rop(offset+j)=lag_3d_bary(self%inodesc(:,j),f,self%xnodes,self%order+1)
      END DO
    END IF
  END SELECT
END IF
DEBUG_STACK_POP
contains
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine tet_eval_all2()!self,cell,f,rop)
! class(oft_scalar_fem), intent(in) :: self
! integer(i4), intent(in) :: cell
! real(r8), intent(in) :: f(:)
! real(r8), intent(out) :: rop(10)
integer(i4) :: i,offset,etmp(2)
real(r8) :: pnorm,nodes1(4)
DEBUG_STACK_PUSH
!---Vertices
pnorm = 1.d0/PRODUCT(self%xnodes(3)-self%xnodes(1:2))
DO i=1,4
  nodes1(i) = f(i) - self%xnodes(1)
  rop(i) = nodes1(i)*(f(i) - self%xnodes(2))*pnorm
END DO
!---Edges
pnorm = 1.d0/((self%xnodes(2)-self%xnodes(1))**2)
DO i=1,6
  offset=(i-1)*1 + 4
  etmp=self%mesh%cell_ed(:,i)
  rop(offset+1) = nodes1(etmp(1))*nodes1(etmp(2))*pnorm
END DO
DEBUG_STACK_POP
end subroutine tet_eval_all2
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine tet_eval_all3()!self,cell,f,rop)
! class(oft_scalar_fem), intent(in) :: self
! integer(i4), intent(in) :: cell
! real(r8), intent(in) :: f(:)
! real(r8), intent(out) :: rop(20)
integer(i4) :: i,offset,etmp(2),ftmp(3)
real(r8) :: pnorm,enorm,fnorm,nodes1(4),nodes2(4)
DEBUG_STACK_PUSH
IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
!---Vertices
pnorm = 1.d0/PRODUCT(self%xnodes(4)-self%xnodes(1:3))
DO i=1,4
  nodes1(i) = f(i)-self%xnodes(1)
  nodes2(i) = nodes1(i)*(f(i)-self%xnodes(2))
  rop(i) =    nodes2(i)*(f(i)-self%xnodes(3))*pnorm
END DO
!---Edges
enorm = 1.d0/((self%xnodes(2)-self%xnodes(1))*PRODUCT(self%xnodes(3)-self%xnodes(1:2)))
DO i=1,6
  offset=(i-1)*2 + 4
  etmp=oriented_edges(:,i)
  rop(offset+1) = nodes1(etmp(1))*nodes2(etmp(2))*enorm
  rop(offset+2) = nodes2(etmp(1))*nodes1(etmp(2))*enorm
END DO
!---Faces
fnorm = 1.d0/((self%xnodes(2)-self%xnodes(1))**3)
DO i=1,4
  offset=(i-1)*1 + 6*2 + 4
  ftmp=oriented_faces(:,i)
  rop(offset+1) = nodes1(ftmp(1))*nodes1(ftmp(2))*nodes1(ftmp(3))*fnorm
END DO
DEBUG_STACK_POP
end subroutine tet_eval_all3
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine tet_eval_all4()!self,cell,f,rop)
! class(oft_scalar_fem), intent(in) :: self
! integer(i4), intent(in) :: cell
! real(r8), intent(in) :: f(:)
! real(r8), intent(out) :: rop(35)
integer(i4) :: i,offset,etmp(2),ftmp(3)
real(r8) :: pnorm,enorm(2),fnorm,nodes1(4),nodes2(4),nodes3(4)
DEBUG_STACK_PUSH
IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
!---Vertices
pnorm = 1.d0/PRODUCT(self%xnodes(5)-self%xnodes(1:4))
DO i=1,4
  nodes1(i) = f(i)-self%xnodes(1)
  nodes2(i) = nodes1(i)*(f(i)-self%xnodes(2))
  nodes3(i) = nodes2(i)*(f(i)-self%xnodes(3))
  rop(i) =    nodes3(i)*(f(i)-self%xnodes(4))*pnorm
END DO
!---Edges
enorm(1) = 1.d0/((self%xnodes(2)-self%xnodes(1))*PRODUCT(self%xnodes(4)-self%xnodes(1:3)))
enorm(2) = 1.d0/(PRODUCT(self%xnodes(3)-self%xnodes(1:2))**2)
DO i=1,6
  offset=(i-1)*3 + 4
  etmp=oriented_edges(:,i)
  rop(offset+1) = nodes1(etmp(1))*nodes3(etmp(2))*enorm(1)
  rop(offset+2) = nodes2(etmp(1))*nodes2(etmp(2))*enorm(2)
  rop(offset+3) = nodes3(etmp(1))*nodes1(etmp(2))*enorm(1)
END DO
!---Faces
fnorm = 1.d0/(((self%xnodes(2)-self%xnodes(1))**2)*PRODUCT(self%xnodes(3)-self%xnodes(1:2)))
DO i=1,4
  offset=(i-1)*3 + 6*3+4
  ftmp=oriented_faces(:,i)
  rop(offset+1) = nodes1(ftmp(1))*nodes1(ftmp(2))*nodes2(ftmp(3))*fnorm
  rop(offset+2) = nodes1(ftmp(1))*nodes2(ftmp(2))*nodes1(ftmp(3))*fnorm
  rop(offset+3) = nodes2(ftmp(1))*nodes1(ftmp(2))*nodes1(ftmp(3))*fnorm
END DO
!---Cell
offset=4*3 + 6*3 + 4
rop(offset+1) = PRODUCT(f-self%xnodes(1))/((self%xnodes(2)-self%xnodes(1))**4)
DEBUG_STACK_POP
end subroutine tet_eval_all4
end subroutine oft_lag_eval_all
!---------------------------------------------------------------------------
!> Evaluate lagrange gradient function
!!
!! @note Evaluation is performed in logical coordinates with the resulting
!! gradient with respect to physical coordinates
!---------------------------------------------------------------------------
subroutine oft_lag_geval(self,cell,dof,f,val,gop)
class(oft_scalar_fem), intent(in) :: self !< Lagrange type for evaluation
integer(i4), intent(in) :: cell !< Cell for evaluation
integer(i4), intent(in) :: dof !< Element to evaluate
real(r8), intent(in) :: f(:) !< Position in cell in logical space
real(r8), intent(in) :: gop(:,:) !< Gradient of lagrange element (dof) at point (f) [3]
real(r8), intent(out) :: val(3) !< Cell Jacobian matrix at point (f) [3,4]
real(r8) :: cofs(4),vtmp(3)
integer(i4) :: ed,etmp(2),fc,ftmp(3),i,ind,finds(16)
DEBUG_STACK_PUSH
IF(self%mesh%type==3)THEN
  val=0.d0
  select case(self%cmap(dof)%type)
    case(1)
      vtmp=dlag_3d(self%inodesp(:,self%cmap(dof)%el),f, &
          self%xnodes,self%order+1)
      DO i=1,3
        val=val+gop(:,i)*vtmp(i)
      END DO
    case(2)
      ind=self%cmap(dof)%ind
      IF(self%mesh%lce(self%cmap(dof)%el,cell)<0)ind=self%order-ind
      vtmp=dlag_3d(self%inodese(:,ind,self%cmap(dof)%el),f, &
          self%xnodes,self%order+1)
      DO i=1,3
        val=val+gop(:,i)*vtmp(i)
      END DO
    case(3)
      CALL hex_grid_forient(self%mesh%lcfo(self%cmap(dof)%el,cell),self%order,finds)
      vtmp=dlag_3d(self%inodesf(:,finds(self%cmap(dof)%ind), &
          self%cmap(dof)%el),f,self%xnodes,self%order+1)
      DO i=1,3
        val=val+gop(:,i)*vtmp(i)
      END DO
    case(4)
      vtmp=dlag_3d(self%inodesc(:,self%cmap(dof)%ind),f, &
          self%xnodes,self%order+1)
      DO i=1,3
        val=val+gop(:,i)*vtmp(i)
      END DO
  end select
ELSE
  IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
  cofs=0.d0
  select case(self%cmap(dof)%type)
    case(1)
      cofs(self%cmap(dof)%el)=dlag_1d(self%order+1,f(self%cmap(dof)%el), &
         self%xnodes,self%order+1)
    case(2)
      etmp=oriented_edges(:,self%cmap(dof)%el)
      cofs(etmp)=dlag_1d_bary(self%cmap(dof)%ind,f(etmp),self%xnodes, &
        self%order+1)
    case(3)
      ftmp=oriented_faces(:,self%cmap(dof)%el)
      cofs(ftmp)=dlag_2d_bary(self%inodesf(:,self%cmap(dof)%ind,1), &
        f(ftmp),self%xnodes,self%order+1)
    case(4)
      cofs=dlag_3d_bary(self%inodesc(:,self%cmap(dof)%ind), &
        f,self%xnodes,self%order+1)
  end select
  !---Sum contributions
  val=0.d0
  do i=1,4
    val=val+gop(:,i)*cofs(i)
  end do
END IF
DEBUG_STACK_POP
end subroutine oft_lag_geval
!---------------------------------------------------------------------------
!> Evaluate lagrange gradient function
!!
!! @note Evaluation is performed in logical coordinates with the resulting
!! gradient with respect to physical coordinates
!---------------------------------------------------------------------------
subroutine oft_blag_geval(self,face,dof,f,val,gop)
class(oft_scalar_bfem), intent(in) :: self !< Lagrange type for evaluation
integer(i4), intent(in) :: face !< Cell for evaluation
integer(i4), intent(in) :: dof !< Element to evaluate
real(r8), intent(in) :: f(:) !< Position in cell in logical space
real(r8), intent(out) :: val(3) !< Gradient of lagrange element (dof) at point (f) [3]
real(r8), optional, intent(in) :: gop(3,3) !< Cell Jacobian matrix at point (f) [3,4]
real(r8) :: grads(3,3),cofs(3)
integer(i4) :: ed,etmp(2),fc,ftmp(3),i,finds(16),ind
DEBUG_STACK_PUSH
grads=1.d0
if(present(gop))grads=gop
IF(self%mesh%type==3)THEN
  val=0.d0
  select case(self%cmap(dof)%type)
    case(1)
      DO i=1,2
        val=val+grads(:,i)*dlag_2d(self%inodesp(:,self%cmap(dof)%el),i,f, &
            self%xnodes,self%order+1)
      END DO
    case(2)
      ind=self%cmap(dof)%ind
      IF(self%mesh%lce(self%cmap(dof)%el,face)<0)ind=self%order-ind
      DO i=1,2
        val=val+grads(:,i)*dlag_2d(self%inodese(:,ind,self%cmap(dof)%el),i,f, &
            self%xnodes,self%order+1)
      END DO
    case(3)
      CALL quad_grid_orient(self%mesh%lco(face),self%order,finds)
      DO i=1,2
        val=val+grads(:,i)*dlag_2d(self%inodesf(:,finds(self%cmap(dof)%ind)), &
            i,f,self%xnodes,self%order+1)
      END DO
  end select
ELSE
  cofs=0.d0
  select case(self%cmap(dof)%type)
    case(1)
      cofs(self%cmap(dof)%el)=dlag_1d(self%order+1,f(self%cmap(dof)%el), &
         self%xnodes,self%order+1)
    case(2)
      ed=self%mesh%lce(self%cmap(dof)%el,face)
      etmp=(/self%mesh%cell_ed(1,self%cmap(dof)%el),self%mesh%cell_ed(2,self%cmap(dof)%el)/)
      call orient_list2(ed,etmp)
      cofs(etmp)=dlag_1d_bary(self%cmap(dof)%ind,f(etmp),self%xnodes, &
        self%order+1)
    case(3)
      fc=self%mesh%lco(face)
      ftmp=(/1,2,3/)
      call orient_listn(fc,ftmp,3_i4)
      cofs(ftmp)=dlag_2d_bary(self%inodesf(:,self%cmap(dof)%ind), &
        f(ftmp),self%xnodes,self%order+1)
  end select
  !---Sum contributions
  val=0.d0
  do i=1,3
    val=val+grads(:,i)*cofs(i)
  end do
END IF
DEBUG_STACK_POP
end subroutine oft_blag_geval
!---------------------------------------------------------------------------
!> Evaluate all lagrange interpolation functions
!!
!! @note Evaluation is performed in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_lag_geval_all(self,cell,f,rop,gop)
class(oft_scalar_fem), intent(in) :: self !< Lagrange type for evaluation
integer(i4), intent(in) :: cell !< Cell for evaluation
real(r8), intent(in) :: f(:) !< Position in cell in logical space
real(r8), contiguous, intent(out) :: rop(:,:) !< Value of interpolation functions at point (f) [3,ncdofs]
real(r8), intent(in) :: gop(:,:) !< Cell Jacobian matrix at point (f) [3,4]
integer(i4) :: i,j,k,offset,ind,inc,etmp(2),ftmp(3)
integer(i4), ALLOCATABLE, DIMENSION(:) :: finds
real(r8) :: cofs(4),val(3)
REAL(r8), ALLOCATABLE, DIMENSION(:,:) :: grid_1d,dgrid_1d
DEBUG_STACK_PUSH
IF(self%mesh%type==3)THEN
  ALLOCATE(grid_1d(self%order+1,3),dgrid_1d(self%order+1,3))
  DO i=1,self%order+1
    grid_1d(i,1) = lag_1d(i,f(1),self%xnodes,self%order+1)
    grid_1d(i,2) = lag_1d(i,f(2),self%xnodes,self%order+1)
    grid_1d(i,3) = lag_1d(i,f(3),self%xnodes,self%order+1)
    dgrid_1d(i,1) = dlag_1d(i,f(1),self%xnodes,self%order+1)
    dgrid_1d(i,2) = dlag_1d(i,f(2),self%xnodes,self%order+1)
    dgrid_1d(i,3) = dlag_1d(i,f(3),self%xnodes,self%order+1)
  END DO
  DO i=1,8
    val=0.d0
    val = val + gop(:,1)*dgrid_1d(self%inodesp(1,i),1) &
      *grid_1d(self%inodesp(2,i),2)*grid_1d(self%inodesp(3,i),3)
    val = val + gop(:,2)*grid_1d(self%inodesp(1,i),1) &
      *dgrid_1d(self%inodesp(2,i),2)*grid_1d(self%inodesp(3,i),3)
    val = val + gop(:,3)*grid_1d(self%inodesp(1,i),1) &
      *grid_1d(self%inodesp(2,i),2)*dgrid_1d(self%inodesp(3,i),3)
    rop(:,i)=val
  END DO
  IF(self%order>1)THEN
    !---Edges
    DO i=1,12
      offset=(i-1)*self%gstruct(2)+8
      ind=1; inc=1
      IF(self%mesh%lce(i,cell)<0)THEN
        ind=self%gstruct(2); inc=-1
      END IF
      DO j=1,self%gstruct(2)
        val=0.d0
        val = val + gop(:,1)*dgrid_1d(self%inodese(1,ind,i),1) &
          *grid_1d(self%inodese(2,ind,i),2)*grid_1d(self%inodese(3,ind,i),3)
        val = val + gop(:,2)*grid_1d(self%inodese(1,ind,i),1) &
          *dgrid_1d(self%inodese(2,ind,i),2)*grid_1d(self%inodese(3,ind,i),3)
        val = val + gop(:,3)*grid_1d(self%inodese(1,ind,i),1) &
          *grid_1d(self%inodese(2,ind,i),2)*dgrid_1d(self%inodese(3,ind,i),3)
        rop(:,offset+j)=val
        ind=ind+inc
      END DO
    END DO
    !---Faces
    ALLOCATE(finds(self%gstruct(3)))
    DO i=1,6
      offset=(i-1)*self%gstruct(3)+12*self%gstruct(2)+8
      CALL hex_grid_forient(self%mesh%lcfo(i,cell),self%order,finds)
      DO j=1,self%gstruct(3)
        val=0.d0
        val = val + gop(:,1)*dgrid_1d(self%inodesf(1,finds(j),i),1) &
          *grid_1d(self%inodesf(2,finds(j),i),2)*grid_1d(self%inodesf(3,finds(j),i),3)
        val = val + gop(:,2)*grid_1d(self%inodesf(1,finds(j),i),1) &
          *dgrid_1d(self%inodesf(2,finds(j),i),2)*grid_1d(self%inodesf(3,finds(j),i),3)
        val = val + gop(:,3)*grid_1d(self%inodesf(1,finds(j),i),1) &
          *grid_1d(self%inodesf(2,finds(j),i),2)*dgrid_1d(self%inodesf(3,finds(j),i),3)
        rop(:,offset+j)=val
      END DO
    END DO
    DEALLOCATE(finds)
    !---Cell
    offset=6*self%gstruct(3)+12*self%gstruct(2)+8
    DO j=1,self%gstruct(4)
      val=0.d0
      val = val + gop(:,1)*dgrid_1d(self%inodesc(1,j),1) &
        *grid_1d(self%inodesc(2,j),2)*grid_1d(self%inodesc(3,j),3)
      val = val + gop(:,2)*grid_1d(self%inodesc(1,j),1) &
        *dgrid_1d(self%inodesc(2,j),2)*grid_1d(self%inodesc(3,j),3)
      val = val + gop(:,3)*grid_1d(self%inodesc(1,j),1) &
        *grid_1d(self%inodesc(2,j),2)*dgrid_1d(self%inodesc(3,j),3)
      rop(:,offset+j)=val
    END DO
  END IF
  DEALLOCATE(grid_1d,dgrid_1d)
ELSE
  SELECT CASE(self%order)
  CASE(1)
    DO i=1,4
      rop(:,i)=gop(:,i)
    END DO
  CASE(2)
    CALL tet_geval_all2()!self,cell,f,gop,rop)
  CASE(3)
    CALL tet_geval_all3()!self,cell,f,gop,rop)
  CASE(4)
    CALL tet_geval_all4()!self,cell,f,gop,rop)
  CASE DEFAULT
    IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
    DO i=1,4
      rop(:,i)=gop(:,i)*dlag_1d(self%order+1,f(i),self%xnodes,self%order+1)
    END DO
    IF(self%order>1)THEN
      !---Edges
      DO i=1,6
        offset=(i-1)*self%gstruct(2)+4
        etmp=oriented_edges(:,i)
        DO j=1,self%gstruct(2)
          cofs(1:2)=dlag_1d_bary(j,f(etmp),self%xnodes,self%order+1)
          val=0.d0
          DO k=1,2
            val=val+gop(:,etmp(k))*cofs(k)
          END DO
          rop(:,offset+j)=val
        END DO
      END DO
      !---Faces
      DO i=1,4
        offset=(i-1)*self%gstruct(3)+6*self%gstruct(2)+4
        ftmp=oriented_faces(:,i)
        DO j=1,self%gstruct(3)
          cofs(1:3)=dlag_2d_bary(self%inodesf(:,j,1),f(ftmp),self%xnodes,self%order+1)
          val=0.d0
          DO k=1,3
            val=val+gop(:,ftmp(k))*cofs(k)
          END DO
          rop(:,offset+j)=val
        END DO
      END DO
      !---Cell
      offset=4*self%gstruct(3)+6*self%gstruct(2)+4
      DO j=1,self%gstruct(4)
        cofs=dlag_3d_bary(self%inodesc(:,j),f,self%xnodes,self%order+1)
        val=0.d0
        DO k=1,4
          val=val+gop(:,k)*cofs(k)
        END DO
        rop(:,offset+j)=val
      END DO
    END IF
  END SELECT
END IF
DEBUG_STACK_POP
contains
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine tet_geval_all2()!self,cell,f,gop,rop)
! class(oft_scalar_fem), intent(in) :: self
! integer(i4), intent(in) :: cell
! real(r8), intent(in) :: f(:)
! real(r8), intent(in) :: gop(3,4)
! real(r8), intent(out) :: rop(3,10)
! integer(i4) :: i,j,k,offset,ind,inc,etmp(2),ftmp(3)
! integer(i4), ALLOCATABLE, DIMENSION(:) :: finds
! real(r8) :: cof1,cofs2(2),pnorm,nodes1(4)
integer(i4) :: i,offset,etmp(2)
real(r8) :: cof1,cofs2(2),pnorm,nodes1(4)
DEBUG_STACK_PUSH
!---Vertices
pnorm = 1.d0/PRODUCT(self%xnodes(3)-self%xnodes(1:2))
DO i=1,4
  nodes1(i) = f(i) - self%xnodes(1)
  rop(:,i) = gop(:,i)*(nodes1(i) + (f(i) - self%xnodes(2)))*pnorm
END DO
!---Edges
pnorm = 1.d0/((self%xnodes(2)-self%xnodes(1))**2)
DO i=1,6
  offset=(i-1)*1 + 4
  etmp=self%mesh%cell_ed(:,i)
  rop(:,offset+1) = (gop(:,etmp(1))*nodes1(etmp(2)) &
    + gop(:,etmp(2))*nodes1(etmp(1)))*pnorm
END DO
DEBUG_STACK_POP
end subroutine tet_geval_all2
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine tet_geval_all3()!self,cell,f,gop,rop)
! class(oft_scalar_fem), intent(in) :: self
! integer(i4), intent(in) :: cell
! real(r8), intent(in) :: f(:)
! real(r8), intent(in) :: gop(3,4)
! real(r8), intent(out) :: rop(3,20)
integer(i4) :: i,offset,etmp(2),ftmp(3)
real(r8) :: pnorm,enorm,fnorm,nodes1(4),nodes2(4),dnodes2(4)
DEBUG_STACK_PUSH
IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
!---Vertices
pnorm = 1.d0/PRODUCT(self%xnodes(4)-self%xnodes(1:3))
DO i=1,4
  nodes1(i) = f(i) - self%xnodes(1)
  nodes2(i) = nodes1(i)*(f(i) - self%xnodes(2))
  dnodes2(i) = (f(i) - self%xnodes(1)) + (f(i) - self%xnodes(2))
  rop(:,i) = gop(:,i)*(nodes2(i) + dnodes2(i)*(f(i) - self%xnodes(3)))*pnorm
END DO
!---Edges
enorm = 1.d0/((self%xnodes(2)-self%xnodes(1))*PRODUCT(self%xnodes(3)-self%xnodes(1:2)))
DO i=1,6
  offset=(i-1)*2 + 4
  etmp=oriented_edges(:,i)
  rop(:,offset+1) = (gop(:,etmp(1))*nodes2(etmp(2)) &
    + gop(:,etmp(2))*nodes1(etmp(1))*dnodes2(etmp(2)))*enorm
  rop(:,offset+2) = (gop(:,etmp(1))*dnodes2(etmp(1))*nodes1(etmp(2)) &
    + gop(:,etmp(2))*nodes2(etmp(1)))*enorm
END DO
!---Faces
fnorm = 1.d0/((self%xnodes(2)-self%xnodes(1))**3)
DO i=1,4
  offset=(i-1)*1 + 6*2 + 4
  ftmp=oriented_faces(:,i)
  rop(:,offset+1) = (gop(:,ftmp(1))*nodes1(ftmp(2))*nodes1(ftmp(3)) &
    + gop(:,ftmp(2))*nodes1(ftmp(1))*nodes1(ftmp(3)) &
    + gop(:,ftmp(3))*nodes1(ftmp(1))*nodes1(ftmp(2)))*fnorm
END DO
DEBUG_STACK_POP
end subroutine tet_geval_all3
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
subroutine tet_geval_all4()!self,cell,f,gop,rop)
! class(oft_scalar_fem), intent(in) :: self
! integer(i4), intent(in) :: cell
! real(r8), intent(in) :: f(:)
! real(r8), intent(in) :: gop(3,4)
! real(r8), intent(out) :: rop(3,35)
! integer(i4) :: i,offset,etmp(2),ftmp(3)
integer(i4) :: i,offset,etmp(2),ftmp(3)
real(r8) :: pnorm,enorm(2),fnorm,cofs4(4)
real(r8) :: nodes1(4),nodes2(4),nodes3(4),dnodes2(4),dnodes3(4)
DEBUG_STACK_PUSH
IF(oriented_cell/=cell)CALL mesh_local_orient(self%mesh,cell)
!---Vertices
pnorm = 1.d0/PRODUCT(self%xnodes(5)-self%xnodes(1:4))
DO i=1,4
  nodes1(i) = f(i) - self%xnodes(1)
  nodes2(i) = nodes1(i)*(f(i) - self%xnodes(2))
  nodes3(i) = nodes2(i)*(f(i) - self%xnodes(3))
  dnodes2(i) = (f(i) - self%xnodes(1)) + (f(i) - self%xnodes(2))
  dnodes3(i) = nodes2(i) + dnodes2(i)*(f(i) - self%xnodes(3))
  rop(:,i) = gop(:,i)*(nodes3(i) + dnodes3(i)*(f(i) - self%xnodes(4)))*pnorm
END DO
!---Edges
enorm(1) = 1.d0/((self%xnodes(2)-self%xnodes(1))*PRODUCT(self%xnodes(4)-self%xnodes(1:3)))
enorm(2) = 1.d0/(PRODUCT(self%xnodes(3)-self%xnodes(1:2))**2)
DO i=1,6
  offset=(i-1)*3 + 4
  etmp=oriented_edges(:,i)
  rop(:,offset+1) = (gop(:,etmp(1))*nodes3(etmp(2)) &
    + gop(:,etmp(2))*nodes1(etmp(1))*dnodes3(etmp(2)))*enorm(1)
  rop(:,offset+2) = (gop(:,etmp(1))*dnodes2(etmp(1))*nodes2(etmp(2)) &
    + gop(:,etmp(2))*nodes2(etmp(1))*dnodes2(etmp(2)))*enorm(2)
  rop(:,offset+3) = (gop(:,etmp(1))*dnodes3(etmp(1))*nodes1(etmp(2)) &
    + gop(:,etmp(2))*nodes3(etmp(1)))*enorm(1)
END DO
!---Faces
fnorm = 1.d0/(((self%xnodes(2)-self%xnodes(1))**2)*PRODUCT(self%xnodes(3)-self%xnodes(1:2)))
DO i=1,4
  offset=(i-1)*3 + 6*3 + 4
  ftmp=oriented_faces(:,i)
  rop(:,offset+1) = (gop(:,ftmp(1))*nodes1(ftmp(2))*nodes2(ftmp(3)) &
    + gop(:,ftmp(2))*nodes1(ftmp(1))*nodes2(ftmp(3)) &
    + gop(:,ftmp(3))*nodes1(ftmp(1))*nodes1(ftmp(2))*dnodes2(ftmp(3)))*fnorm
  rop(:,offset+2) = (gop(:,ftmp(1))*nodes2(ftmp(2))*nodes1(ftmp(3)) &
    + gop(:,ftmp(2))*nodes1(ftmp(1))*dnodes2(ftmp(2))*nodes1(ftmp(3)) &
    + gop(:,ftmp(3))*nodes1(ftmp(1))*nodes2(ftmp(2)))*fnorm
  rop(:,offset+3) = (gop(:,ftmp(1))*dnodes2(ftmp(1))*nodes1(ftmp(2))*nodes1(ftmp(3)) &
    + gop(:,ftmp(2))*nodes2(ftmp(1))*nodes1(ftmp(3)) &
    + gop(:,ftmp(3))*nodes2(ftmp(1))*nodes1(ftmp(2)))*fnorm
END DO
!---Cell
offset=4*3 + 6*3 + 4
cofs4 = (/PRODUCT(nodes1(2:4)),PRODUCT(nodes1((/1,3,4/))), &
  PRODUCT(nodes1((/1,2,4/))),PRODUCT(nodes1(1:3))/)
rop(:,offset+1) = MATMUL(gop,cofs4)/((self%xnodes(2)-self%xnodes(1))**4)
DEBUG_STACK_POP
end subroutine tet_geval_all4
end subroutine oft_lag_geval_all
!---------------------------------------------------------------------------
!> Evaluate lagrange gradient function
!!
!! @note Evaluation is performed in logical coordinates with the resulting
!! gradient with respect to physical coordinates
!---------------------------------------------------------------------------
subroutine oft_lag_d2eval(self,cell,dof,f,val,g2op)
class(oft_scalar_fem), intent(in) :: self !< Lagrange type for evaluation
integer(i4), intent(in) :: cell !< Cell for evaluation
integer(i4), intent(in) :: dof !< Element to evaluate
real(r8), intent(in) :: f(:) !< Position in cell in logical space
real(r8), intent(out) :: val(6) !< Second derivatives of function [6]
real(r8), optional, intent(in) :: g2op(6,10) !< Grid Hessian [6,6]
real(r8) :: grad(2),cofs(10),d2etmp(3),d2ftmp(6)
integer(i4) :: ed,etmp(2),fc,ftmp(3),i,j,k,ind,finds(16)
integer(i4), parameter :: pmap(4)=(/1,5,8,10/)
integer(i4), parameter :: emap(6)=(/4,7,9,6,3,2/)
integer(i4), parameter :: fmap(4,4)=RESHAPE((/1,2,3,4,2,5,6,7,3,6,8,9,4,7,9,10/),(/4,4/))
cofs=0.d0
IF(self%mesh%type==3)THEN
  val=0.d0
  select case(self%cmap(dof)%type)
    case(1)
      k=0
      DO i=1,3
        DO j=i,3
          k=k+1
          val=val+g2op(:,k)*d2lag_3d(self%inodesp(:,self%cmap(dof)%el),i,j,f, &
            self%xnodes,self%order+1)
        END DO
      END DO
    case(2)
      ind=self%cmap(dof)%ind
      IF(self%mesh%lce(self%cmap(dof)%el,cell)<0)ind=self%order-ind
      k=0
      DO i=1,3
        DO j=i,3
          k=k+1
          val=val+g2op(:,k)*d2lag_3d(self%inodese(:,ind,self%cmap(dof)%el),i,j,f, &
            self%xnodes,self%order+1)
        END DO
      END DO
    case(3)
      CALL hex_grid_forient(self%mesh%lcfo(self%cmap(dof)%el,cell),self%order,finds)
      k=0
      DO i=1,3
        DO j=i,3
          k=k+1
          val=val+g2op(:,k)*d2lag_3d(self%inodesf(:,finds(self%cmap(dof)%ind), &
              self%cmap(dof)%el),i,j,f,self%xnodes,self%order+1)
        END DO
      END DO
    case(4)
      k=0
      DO i=1,3
        DO j=i,3
          k=k+1
          val=val+g2op(:,k)*d2lag_3d(self%inodesc(:,self%cmap(dof)%ind), &
              i,j,f,self%xnodes,self%order+1)
        END DO
      END DO
  end select
ELSE
  select case(self%cmap(dof)%type)
    case(1)
      cofs(pmap(self%cmap(dof)%el))=d2lag_1d(self%order+1,f(self%cmap(dof)%el), &
          self%xnodes,self%order+1)
    case(2)
      ed=self%mesh%lce(self%cmap(dof)%el,cell)
      etmp=self%mesh%cell_ed(:,self%cmap(dof)%el)
      call orient_list2(ed,etmp)
      d2etmp=d2lag_1d_bary(self%cmap(dof)%ind,f(etmp),self%xnodes, &
        self%order+1)
      cofs(pmap(etmp(1)))=d2etmp(1)
      cofs(emap(self%cmap(dof)%el))=d2etmp(2)
      cofs(pmap(etmp(2)))=d2etmp(3)
    case(3)
      fc=self%mesh%lcfo(self%cmap(dof)%el,cell)
      ftmp=self%mesh%cell_fc(:,self%cmap(dof)%el)
      call orient_listn(fc,ftmp,3_i4)
      d2ftmp=d2lag_2d_bary(self%inodesf(:,self%cmap(dof)%ind,1), &
        f(ftmp),self%xnodes,self%order+1)
      cofs(pmap(ftmp(1)))=d2ftmp(1)
      cofs(fmap(ftmp(1),ftmp(2)))=d2ftmp(2)
      cofs(fmap(ftmp(1),ftmp(3)))=d2ftmp(3)
      cofs(pmap(ftmp(2)))=d2ftmp(4)
      cofs(fmap(ftmp(2),ftmp(3)))=d2ftmp(5)
      cofs(pmap(ftmp(3)))=d2ftmp(6)
    case(4)
      cofs=d2lag_3d_bary(self%inodesc(:,self%cmap(dof)%ind), &
        f,self%xnodes,self%order+1)
  end select
  !---Sum contributions
  val=0.d0
  do i=1,10
    val=val+g2op(:,i)*cofs(i)
  end do
END IF
end subroutine oft_lag_d2eval
!---------------------------------------------------------------------------
!> Evaluate lagrange gradient function
!!
!! @note Evaluation is performed in logical coordinates with the resulting
!! gradient with respect to physical coordinates
!---------------------------------------------------------------------------
subroutine oft_blag_d2eval(self,cell,dof,f,val,g2op)
class(oft_scalar_bfem), intent(in) :: self !<  Lagrange type for evaluation
integer(i4), intent(in) :: cell !< Cell for evaluation
integer(i4), intent(in) :: dof !< Element to evaluate
real(r8), intent(in) :: f(:) !< Position in cell in logical space
real(r8), intent(out) :: val(6) !< Second derivatives of function [6]
real(r8), optional, intent(in) :: g2op(6,6) !< Grid Hessian [6,6]
real(r8) :: grad(2),cofs(10),d2etmp(3),d2ftmp(6)
integer(i4) :: ed,etmp(2),fc,ftmp(3),i,j,k,ind,finds(16)
integer(i4), parameter :: pmap(3)=(/1,4,6/)
integer(i4), parameter :: emap(3)=(/5,3,2/)
integer(i4), parameter :: fmap(3,3)=RESHAPE((/1,2,3,2,4,5,3,4,6/),(/3,3/))
IF(self%mesh%type==3)THEN
  val=0.d0
  select case(self%cmap(dof)%type)
    case(1)
      k=0
      DO i=1,3
        DO j=i,3
          k=k+1
          val=val+g2op(:,k)*d2lag_2d(self%inodesp(:,self%cmap(dof)%el),i,j,f, &
            self%xnodes,self%order+1)
        END DO
      END DO
    case(2)
      ind=self%cmap(dof)%ind
      IF(self%mesh%lce(self%cmap(dof)%el,cell)<0)ind=self%order-ind
      k=0
      DO i=1,3
        DO j=i,3
          k=k+1
          val=val+g2op(:,k)*d2lag_2d(self%inodese(:,ind,self%cmap(dof)%el),i,j,f, &
            self%xnodes,self%order+1)
        END DO
      END DO
    case(3)
      CALL quad_grid_orient(self%mesh%lco(cell),self%order,finds)
      k=0
      DO i=1,3
        DO j=i,3
          k=k+1
          val=val+g2op(:,k)*d2lag_2d(self%inodesf(:,finds(self%cmap(dof)%ind)), &
            i,j,f,self%xnodes,self%order+1)
        END DO
      END DO
  end select
ELSE
  cofs=0.d0
  select case(self%cmap(dof)%type)
    case(1)
      cofs(pmap(self%cmap(dof)%el))=d2lag_1d(self%order+1,f(self%cmap(dof)%el), &
          self%xnodes,self%order+1)
    case(2)
      ed=self%mesh%lce(self%cmap(dof)%el,cell)
      etmp=self%mesh%cell_ed(:,self%cmap(dof)%el)
      call orient_list2(ed,etmp)
      d2etmp=d2lag_1d_bary(self%cmap(dof)%ind,f(etmp),self%xnodes, &
        self%order+1)
      cofs(pmap(etmp(1)))=d2etmp(1)
      cofs(emap(self%cmap(dof)%el))=d2etmp(2)
      cofs(pmap(etmp(2)))=d2etmp(3)
    case(3)
      fc=self%mesh%lco(cell)
      ftmp=[1,2,3]
      call orient_listn(fc,ftmp,3_i4)
      d2ftmp=d2lag_2d_bary(self%inodesf(:,self%cmap(dof)%ind), &
        f(ftmp),self%xnodes,self%order+1)
      cofs(pmap(ftmp(1)))=d2ftmp(1)
      cofs(fmap(ftmp(1),ftmp(2)))=d2ftmp(2)
      cofs(fmap(ftmp(1),ftmp(3)))=d2ftmp(3)
      cofs(pmap(ftmp(2)))=d2ftmp(4)
      cofs(fmap(ftmp(2),ftmp(3)))=d2ftmp(5)
      cofs(pmap(ftmp(3)))=d2ftmp(6)
  end select
  !---Sum contributions
  val=0.d0
  do i=1,6
    val=val+g2op(:,i)*cofs(i)
  end do
END IF
end subroutine oft_blag_d2eval
!---------------------------------------------------------------------------
!> Retrieve lagrange node locations in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_lag_npos(self,cell,dof,f)
class(oft_scalar_fem), intent(in) :: self !< Lagrange type for evaluation
integer(i4), intent(in) :: cell !< Cell for evaluation
integer(i4), intent(in) :: dof !< Element to locate
real(r8), intent(out) :: f(:) !< Position of node in logical space
integer(i4) :: i,ed,etmp(2),fc,ftmp(3),ind,finds(16)
DEBUG_STACK_PUSH
!---
IF(self%mesh%type==3)THEN
  select case(self%cmap(dof)%type)
    case(1)
      f(1:3)=self%xnodes(self%inodesp(:,self%cmap(dof)%el))
    case(2)
      ind=self%cmap(dof)%ind
      IF(self%mesh%lce(self%cmap(dof)%el,cell)<0)ind=self%order-ind
      f(1:3)=self%xnodes(self%inodese(:,ind,self%cmap(dof)%el))
    case(3)
      CALL hex_grid_forient(self%mesh%lcfo(self%cmap(dof)%el,cell),self%order,finds)
      ind=finds(self%cmap(dof)%ind)
      f(1:3)=self%xnodes(self%inodesf(:,ind,self%cmap(dof)%el))
    case(4)
      f(1:3)=self%xnodes(self%inodesc(:,self%cmap(dof)%ind))
  end select
ELSE
  f=0._r8
  select case(self%cmap(dof)%type)
    case(1)
      f(self%cmap(dof)%el)=1.d0
    case(2)
      ed=self%mesh%lce(self%cmap(dof)%el,cell)
      etmp=self%mesh%cell_ed(:,self%cmap(dof)%el)
      call orient_list2(ed,etmp)
      f(etmp(1))=self%xnodes(self%cmap(dof)%ind+1)
      f(etmp(2))=1.d0-f(etmp(1))
    case(3)
      fc=self%mesh%lcfo(self%cmap(dof)%el,cell)
      ftmp=self%mesh%cell_fc(:,self%cmap(dof)%el)
      call orient_listn(fc,ftmp,3_i4)
      f(ftmp(1))=self%xnodes(self%inodesf(1,self%cmap(dof)%ind,1)+1)
      f(ftmp(2))=self%xnodes(self%inodesf(2,self%cmap(dof)%ind,1)+1)
      f(ftmp(3))=1.d0-f(ftmp(1))-f(ftmp(2))
    case(4)
      f(1)=self%xnodes(self%inodesc(1,self%cmap(dof)%ind)+1)
      f(2)=self%xnodes(self%inodesc(2,self%cmap(dof)%ind)+1)
      f(3)=self%xnodes(self%inodesc(2,self%cmap(dof)%ind)+1)
      f(4)=1.d0-SUM(f(1:3))
  end select
END IF
DEBUG_STACK_POP
end subroutine oft_lag_npos
!---------------------------------------------------------------------------
!> Retrieve all lagrange node locations in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_lag_nodes(order,ed_nodes,fc_nodes,c_nodes)
integer(i4), intent(in) :: order !< Needs docs
real(r8), pointer, intent(out) :: ed_nodes(:,:) !< Needs docs
real(r8), pointer, intent(out) :: fc_nodes(:,:) !< Needs docs
real(r8), pointer, intent(out) :: c_nodes(:,:) !< Needs docs
INTEGER(i4) :: i,n
INTEGER(i4), POINTER, DIMENSION(:,:,:) :: inodesf
INTEGER(i4), POINTER, DIMENSION(:,:) :: inodesc
REAL(r8), POINTER, DIMENSION(:) :: xnodes
DEBUG_STACK_PUSH
!
IF(order>1)THEN
  NULLIFY(inodesf,inodesc,xnodes)
  CALL tet_3d_grid(order,xnodes,inodesf,inodesc)
  !
  ALLOCATE(ed_nodes(2,order-1))
  ed_nodes(1,:)=xnodes(2:order)
  ed_nodes(2,:)=1.d0-ed_nodes(1,:)
  !
  IF(ASSOCIATED(inodesf))THEN
    n=SIZE(inodesf, DIM=2)
    ALLOCATE(fc_nodes(3,n))
    DO i=1,n
      fc_nodes(1:2,i)=xnodes(inodesf(:,i,1)+1)
      fc_nodes(3,i)=1.d0-SUM(fc_nodes(1:2,i))
    END DO
    DEALLOCATE(inodesf)
  ELSE
    ALLOCATE(fc_nodes(1,1))
  END IF
  !
  IF(ASSOCIATED(inodesc))THEN
    n=SIZE(inodesc, DIM=2)
    ALLOCATE(c_nodes(4,n))
    DO i=1,n
      c_nodes(1:3,i)=xnodes(inodesc(:,i)+1)
      c_nodes(4,i)=1.d0-SUM(c_nodes(1:3,i))
    END DO
    DEALLOCATE(inodesc)
  ELSE
    ALLOCATE(c_nodes(1,1))
  END IF
ELSE
  ALLOCATE(ed_nodes(1,1))
END IF
! select case(order)
!   case(1)
!     ALLOCATE(ed_nodes(1,1),fc_nodes(1,1),c_nodes(1,1))
!   case(2)
!     ALLOCATE(fc_nodes(1,1),c_nodes(1,1))
!     ALLOCATE(ed_nodes(2,1))
!     ed_nodes(:,1)=(/.5d0,.5d0/)
!   case(3)
!     ALLOCATE(c_nodes(1,1))
!     ALLOCATE(ed_nodes(2,2))
!     ed_nodes(:,1)=(/1.d0,2.d0/)/3
!     ed_nodes(:,2)=(/2.d0,1.d0/)/3
!     ALLOCATE(fc_nodes(3,1))
!     fc_nodes(:,1)=(/1.d0,1.d0,1.d0/)/3
!   case(4)
!     ALLOCATE(ed_nodes(2,3))
!     ed_nodes(:,1)=(/1.d0,3.d0/)/4
!     ed_nodes(:,2)=(/2.d0,2.d0/)/4
!     ed_nodes(:,3)=(/3.d0,1.d0/)/4
!     ALLOCATE(fc_nodes(3,3))
!     fc_nodes(:,1)=(/1.d0,1.d0,2.d0/)/4
!     fc_nodes(:,2)=(/1.d0,2.d0,1.d0/)/4
!     fc_nodes(:,3)=(/2.d0,1.d0,1.d0/)/4
!     ALLOCATE(c_nodes(4,1))
!     c_nodes(:,1)=(/1.d0,1.d0,1.d0,1.d0/)/4
! end select
DEBUG_STACK_POP
end subroutine oft_lag_nodes
!---------------------------------------------------------------------------
!> Retrieve lagrange node locations in logical coordinates
!---------------------------------------------------------------------------
subroutine oft_blag_npos(self,cell,dof,f)
class(oft_scalar_bfem), intent(in) :: self !< Lagrange type for evaluation
integer(i4), intent(in) :: cell !< Cell for evaluation
integer(i4), intent(in) :: dof !< Element to locate
real(r8), intent(out) :: f(:) !< Position of node in logical space
integer(i4) :: i,ed,etmp(2),fc,ftmp(3),ind,finds(16)
DEBUG_STACK_PUSH
!---
IF(self%mesh%type==3)THEN
  select case(self%cmap(dof)%type)
    case(1)
      f(1:2)=self%xnodes(self%inodesp(1:2,self%cmap(dof)%el))
    case(2)
      ind=self%cmap(dof)%ind
      IF(self%mesh%lce(self%cmap(dof)%el,cell)<0)ind=self%order-ind
      f(1:2)=self%xnodes(self%inodese(1:2,ind,self%cmap(dof)%el))
    case(3)
      CALL quad_grid_orient(self%mesh%lco(cell),self%order,finds)
      ind=finds(self%cmap(dof)%ind)
      f(1:2)=self%xnodes(self%inodesf(:,ind))
  end select
ELSE
  f=0._r8
  select case(self%cmap(dof)%type)
    case(1)
      f(self%cmap(dof)%el)=1.d0
    case(2)
      ed=self%mesh%lce(self%cmap(dof)%el,cell)
      etmp=self%mesh%cell_ed(:,self%cmap(dof)%el)
      call orient_list2(ed,etmp)
      f(etmp(1))=self%xnodes(self%cmap(dof)%ind+1)
      f(etmp(2))=1.d0-f(etmp(1))
    case(3)
      fc=self%mesh%lco(cell)
      ftmp=[1,2,3]
      call orient_listn(fc,ftmp,3_i4)
      f(ftmp(1))=self%xnodes(self%inodesf(1,self%cmap(dof)%ind)+1)
      f(ftmp(2))=self%xnodes(self%inodesf(2,self%cmap(dof)%ind)+1)
      f(ftmp(3))=1.d0-f(ftmp(1))-f(ftmp(2))
  end select
END IF
DEBUG_STACK_POP
end subroutine oft_blag_npos
!---------------------------------------------------------------------------
!> Destroy boundary FE object
!---------------------------------------------------------------------------
SUBROUTINE scalar_bfem_delete(self)
CLASS(oft_scalar_bfem), INTENT(inout) :: self
DEBUG_STACK_PUSH
CALL bfem_delete(self)
IF(ASSOCIATED(self%inodese))DEALLOCATE(self%inodese)
IF(ASSOCIATED(self%inodesf))DEALLOCATE(self%inodesf)
IF(ASSOCIATED(self%xnodes))DEALLOCATE(self%xnodes)
DEBUG_STACK_POP
END SUBROUTINE scalar_bfem_delete
!---------------------------------------------------------------------------
!> Destroy FE object
!---------------------------------------------------------------------------
SUBROUTINE scalar_fem_delete(self)
CLASS(oft_scalar_fem), INTENT(inout) :: self
DEBUG_STACK_PUSH
CALL fem_delete(self)
IF(ASSOCIATED(self%inodese))DEALLOCATE(self%inodese)
IF(ASSOCIATED(self%inodesf))DEALLOCATE(self%inodesf)
IF(ASSOCIATED(self%inodesc))DEALLOCATE(self%inodesc)
IF(ASSOCIATED(self%xnodes))DEALLOCATE(self%xnodes)
IF(ASSOCIATED(self%sn))DEALLOCATE(self%sn)
IF(ASSOCIATED(self%interp_graph))THEN
  IF(ASSOCIATED(self%interp_graph%kr))DEALLOCATE(self%interp_graph%kr)
  IF(ASSOCIATED(self%interp_graph%lc))DEALLOCATE(self%interp_graph%lc)
  DEALLOCATE(self%interp_graph)
END IF
IF(ASSOCIATED(self%interp))THEN
  CALL self%interp%delete()
  DEALLOCATE(self%interp)
END IF
IF(ASSOCIATED(self%vinterp))THEN
  CALL self%vinterp%delete()
  DEALLOCATE(self%vinterp)
END IF
DEBUG_STACK_POP
END SUBROUTINE scalar_fem_delete
end module oft_lag_basis
