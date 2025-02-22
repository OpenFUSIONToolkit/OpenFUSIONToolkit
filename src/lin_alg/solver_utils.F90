!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file oft_solver_utils.F90
!
!> Matrix and vector management routines
!!
!! @authors Chris Hansen
!! @date December 2012
!! @ingroup doxy_oft_lin_alg
!---------------------------------------------------------------------------
MODULE oft_solver_utils
USE oft_base
USE oft_stitching, ONLY: oft_seam, seam_list, oft_global_stitch, &
  oft_stitch_check
USE oft_la_base, ONLY: oft_vector, oft_map, map_list, oft_matrix, oft_matrix_ptr, &
  oft_graph, oft_graph_ptr, oft_matrix_map
USE oft_native_la, ONLY: oft_native_vector, native_vector_cast, oft_native_matrix
USE oft_solver_base, ONLY: oft_solver, oft_bc_proto
USE oft_native_solvers, ONLY: oft_native_cg_solver, native_cg_solver_cast, &
  oft_native_gmres_solver, native_gmres_solver_cast, oft_jblock_precond, &
  oft_ml_precond, oft_veccreate_proto, oft_interp_proto, jblock_precond_cast, &
  ml_precond_cast, oft_diag_scale, diag_scale_cast, oft_bjprecond
USE oft_lu, ONLY: oft_lusolver, oft_ilusolver
#ifdef HAVE_PETSC
USE oft_petsc_la, ONLY: oft_petsc_vector, oft_petsc_vector_cast, oft_petsc_matrix, oft_petsc_matrix_cast
USE oft_petsc_solvers, ONLY: oft_petsc_sjacobi_solver, oft_petsc_cg_solver, oft_petsc_gmres_solver, &
oft_petsc_diagprecond, petsc_cg_solver_cast, petsc_gmres_solver_cast, petsc_sjacobi_solver_cast, &
petsc_diagprecond_cast, oft_petsc_asprecond, oft_petsc_direct_solver, petsc_direct_solver_cast, &
oft_petsc_pre_solver, oft_petsc_luprecond
#endif
IMPLICIT NONE
#include "local.h"
CONTAINS
!---------------------------------------------------------------------------
! SUBROUTINE: create_mlpre
!---------------------------------------------------------------------------
!> Construct Multi-Grid preconditioner
!!
!! This subroutine is a wrapper around specific subroutines for construction
!! from an XML specification or standard native/PETSc preconditioners.
!!
!! @param[out] pre Preconditioner object
!! @param[in] Mats Operator matrices [nlevels]
!! @param[in] levels List of level indices [nlevels]
!! @param[in] nlevels Number of levels
!! @param[in] create_vec Vector creation subroutine
!! @param[in] interp Interpolation subroutine
!! @param[in] inject Restriction subroutine
!! param[in] bc Bondary condition subroutine (optional)
!! @param[in] stype Smoother type (optional)
!! @param[in] df Smoother damping factors [nlevels] (optional)
!! @param[in] nu Number of smoother iterations [nlevels] (optional)
!! @param[in] xml_root Preconditioner definition node (optional)
!---------------------------------------------------------------------------
subroutine create_mlpre(pre,Mats,levels,nlevels,create_vec,interp, &
  inject,bc,stype,df,nu,xml_root)
class(oft_solver), pointer, intent(out) :: pre
TYPE(oft_matrix_ptr), INTENT(in) :: Mats(:)
integer(i4), intent(in) :: levels(:)
integer(i4), intent(in) :: nlevels
procedure(oft_veccreate_proto) :: create_vec
procedure(oft_interp_proto) :: interp
procedure(oft_interp_proto) :: inject
procedure(oft_bc_proto), optional :: bc
integer(i4), optional, intent(in) :: stype
real(r8), optional, intent(in) :: df(:)
integer(i4), optional, intent(in) :: nu(:)
TYPE(xml_node), optional, pointer, intent(in) :: xml_root
!---
NULLIFY(pre)
#ifdef HAVE_XML
IF(ASSOCIATED(oft_env%xml).AND.PRESENT(xml_root))THEN
  IF(ASSOCIATED(xml_root))THEN
    !---Create preconditioner
    CALL create_ml_xml(pre,Mats,levels,nlevels=nlevels, &
      create_vec=create_vec,interp=interp,inject=inject,&
      pre_node=xml_root,bc=bc)
  END IF
END IF
#endif
IF(.NOT.ASSOCIATED(pre))THEN
  IF(use_petsc)THEN
    CALL create_petsc_mlpre(pre,Mats,levels,nlevels=nlevels, &
      create_vec=create_vec,interp=interp,inject=inject, &
      bc=bc,stype=stype,nu=nu,df=df)
  ELSE
    CALL create_native_mlpre(pre,Mats,levels,nlevels=nlevels, &
      create_vec=create_vec,interp=interp,inject=inject, &
      bc=bc,stype=stype,nu=nu,df=df)
  END IF
END IF
end subroutine create_mlpre
!---------------------------------------------------------------------------
! SUBROUTINE: create_native_mlpre
!---------------------------------------------------------------------------
!> Construct native Multi-Grid preconditioner
!!
!! @param[out] pre Preconditioner object
!! @param[in] Mats Operator matrices [nlevels]
!! @param[in] levels List of level indices [nlevels]
!! @param[in] nlevels Number of levels
!! @param[in] create_vec Vector creation subroutine
!! @param[in] interp Interpolation subroutine
!! @param[in] inject Restriction subroutine
!! @param[in] bc Bondary condition subroutine (optional)
!! @param[in] stype Smoother type (optional)
!! @param[in] df Smoother damping factors [nlevels] (optional)
!! @param[in] nu Number of smoother iterations [nlevels] (optional)
!---------------------------------------------------------------------------
subroutine create_native_mlpre(pre,Mats,levels,nlevels,create_vec,interp,inject,bc,stype,df,nu)
class(oft_solver), pointer, intent(out) :: pre
TYPE(oft_matrix_ptr), INTENT(in) :: Mats(:)
integer(i4), intent(in) :: levels(:)
integer(i4), intent(in) :: nlevels
procedure(oft_veccreate_proto) :: create_vec
procedure(oft_interp_proto) :: interp
procedure(oft_interp_proto) :: inject
procedure(oft_bc_proto), optional :: bc
integer(i4), optional, intent(in) :: stype
real(r8), optional, intent(in) :: df(:)
integer(i4), optional, intent(in) :: nu(:)
integer(i4) :: i,smoother
class(oft_ml_precond), pointer :: this_ml
type(oft_jblock_precond), pointer :: this_jac
type(oft_native_cg_solver), pointer :: this_cg
type(oft_native_gmres_solver), pointer :: this_gmres
type(oft_diag_scale), pointer :: this_diag
DEBUG_STACK_PUSH
!---
smoother=1 ! Default to jacobi
IF(PRESENT(stype))smoother=stype
SELECT CASE(smoother)
  CASE(1)
    IF(.NOT.(PRESENT(df).AND.PRESENT(nu)))THEN
      CALL oft_abort('Jacobi smoother requires (df, nu)!','create_native_mlpre',__FILE__)
    END IF
    IF(oft_debug_print(1))WRITE(*,*)'Creating MG-Jacobi',df,nu
  CASE(2)
    IF(.NOT.PRESENT(nu))THEN
      CALL oft_abort('GMRES smoother requires (nu)!','create_native_mlpre',__FILE__)
    END IF
    IF(oft_debug_print(1))WRITE(*,*)'Creating MG-GMRES',nu
  CASE(3)
    IF(.NOT.PRESENT(nu))THEN
      CALL oft_abort('CG smoother requires (nu)!','create_native_mlpre',__FILE__)
    END IF
    IF(oft_debug_print(1))WRITE(*,*)'Creating MG-CG',nu
  CASE DEFAULT
    CALL oft_abort('Unknown smoother type!','create_native_mlpre',__FILE__)
END SELECT
!---Set smoother
ALLOCATE(oft_ml_precond::pre)
SELECT TYPE(pre)
  TYPE IS(oft_ml_precond)
    this_ml=>pre
  CLASS DEFAULT
    CALL oft_abort('Error in precon allocation!','create_native_mlpre',__FILE__)
END SELECT
!---Setup MG smoother
this_ml%interp=>interp
this_ml%inject=>inject
this_ml%vec_create=>create_vec
this_ml%level=ABS(levels(nlevels))
!---Allocate top-level smoother
SELECT CASE(smoother)
  CASE(1) ! Jacobi
    ALLOCATE(oft_jblock_precond::this_ml%smooth_up)
    !---Setup top-level smoother
    IF(jblock_precond_cast(this_jac,this_ml%smooth_up)<0) &
      CALL oft_abort('Failed to allocate "this_jac".','create_native_mlpre',__FILE__)
    this_jac%A=>Mats(nlevels)%M
    IF(PRESENT(bc))this_jac%bc=>bc
    this_jac%df=df(nlevels)
    this_jac%its=nu(nlevels)
  CASE(2) ! FGMRES
    ALLOCATE(oft_native_gmres_solver::this_ml%smooth_up)
    !---Setup top-level smoother
    IF(native_gmres_solver_cast(this_gmres,this_ml%smooth_up)<0) &
      CALL oft_abort('Failed to allocate "this_gmres".','create_native_mlpre',__FILE__)
    this_gmres%A=>Mats(nlevels)%M
    IF(PRESENT(bc))this_gmres%bc=>bc
    this_gmres%nrits=nu(nlevels)
    this_gmres%its=nu(nlevels)
    !---
    ALLOCATE(oft_diag_scale::this_gmres%pre)
    IF(diag_scale_cast(this_diag,this_gmres%pre)<0) &
      CALL oft_abort('Failed to allocate "this_diag".','create_native_mlpre',__FILE__)
    this_diag%A=>Mats(nlevels)%M
  CASE(3) ! CG
    ALLOCATE(oft_native_cg_solver::this_ml%smooth_up)
    !---Setup top-level smoother
    IF(native_cg_solver_cast(this_cg,this_ml%smooth_up)<0) &
      CALL oft_abort('Failed to allocate "this_cg".','create_native_mlpre',__FILE__)
    this_cg%A=>Mats(nlevels)%M
    IF(PRESENT(bc))this_cg%bc=>bc
    this_cg%its=nu(nlevels)
END SELECT
this_ml%smooth_down=>this_ml%smooth_up
this_ml%symmetric=.TRUE.
!---Build ML preconditioner
DO i=nlevels-1,1,-1
  !---Set smoother
  ALLOCATE(oft_ml_precond::this_ml%base_solve)
  IF(ml_precond_cast(this_ml,this_ml%base_solve)<0) &
    CALL oft_abort('Failed to allocate "this_ml".','create_native_mlpre',__FILE__)
  !---Setup MG smoother
  this_ml%interp=>interp
  this_ml%inject=>inject
  this_ml%vec_create=>create_vec
  this_ml%level=ABS(levels(i))
  IF(levels(i)<0)CYCLE
  !---Allocate level smoother
  SELECT CASE(smoother)
    CASE(1)
      ALLOCATE(oft_jblock_precond::this_ml%smooth_up)
      !---Setup top-level smoother
      IF(jblock_precond_cast(this_jac,this_ml%smooth_up)<0) &
        CALL oft_abort('Failed to allocate "this_jac".','create_native_mlpre',__FILE__)
      this_jac%A=>Mats(i)%M
      IF(PRESENT(bc))this_jac%bc=>bc
      this_jac%df=df(i)
      this_jac%its=nu(i)
      this_ml%smooth_down=>this_ml%smooth_up
      this_ml%symmetric=.TRUE.
    CASE(2) ! FGMRES
      ALLOCATE(oft_native_gmres_solver::this_ml%smooth_up)
      !---Setup top-level smoother
      IF(native_gmres_solver_cast(this_gmres,this_ml%smooth_up)<0) &
        CALL oft_abort('Failed to allocate "this_gmres".','create_native_mlpre',__FILE__)
      this_gmres%A=>Mats(i)%M
      IF(PRESENT(bc))this_gmres%bc=>bc
      this_gmres%nrits=MIN(20,nu(i))
      this_gmres%its=nu(i)
      !---
      ALLOCATE(oft_diag_scale::this_gmres%pre)
      IF(diag_scale_cast(this_diag,this_gmres%pre)<0) &
        CALL oft_abort('Failed to allocate "this_diag".','create_native_mlpre',__FILE__)
      this_diag%A=>Mats(i)%M
      IF(i>1)THEN
        this_ml%smooth_down=>this_ml%smooth_up
        this_ml%symmetric=.TRUE.
      END IF
    CASE(3) ! CG
      ALLOCATE(oft_native_cg_solver::this_ml%smooth_up)
      !---Setup top-level smoother
      IF(native_cg_solver_cast(this_cg,this_ml%smooth_up)<0) &
        CALL oft_abort('Failed to allocate "this_cg".','create_native_mlpre',__FILE__)
      this_cg%A=>Mats(i)%M
      IF(PRESENT(bc))this_cg%bc=>bc
      this_cg%its=nu(i)
      IF(i>1)THEN
        this_ml%smooth_down=>this_ml%smooth_up
        this_ml%symmetric=.TRUE.
      END IF
  END SELECT
END DO
DEBUG_STACK_POP
end subroutine create_native_mlpre
!---------------------------------------------------------------------------
! SUBROUTINE: create_petsc_mlpre
!---------------------------------------------------------------------------
!> Construct PETSc Multi-Grid preconditioner using native mechanics
!!
!! @param[out] pre Preconditioner object
!! @param[in] Mats Operator matrices [nlevels]
!! @param[in] levels List of level indices [nlevels]
!! @param[in] nlevels Number of levels
!! @param[in] create_vec Vector creation subroutine
!! @param[in] interp Interpolation subroutine
!! @param[in] inject Restriction subroutine
!! @param[in] bc Bondary condition subroutine (optional)
!! @param[in] stype Smoother type (optional)
!! @param[in] df Smoother damping factors [nlevels] (optional)
!! @param[in] nu Number of smoother iterations [nlevels] (optional)
!---------------------------------------------------------------------------
subroutine create_petsc_mlpre(pre,Mats,levels,nlevels,create_vec,interp,inject,bc,stype,df,nu)
class(oft_solver), pointer, intent(out) :: pre
TYPE(oft_matrix_ptr), INTENT(in) :: Mats(:)
integer(i4), intent(in) :: levels(:)
integer(i4), intent(in) :: nlevels
procedure(oft_veccreate_proto) :: create_vec
procedure(oft_interp_proto) :: interp
procedure(oft_interp_proto) :: inject
procedure(oft_bc_proto), optional :: bc
integer(i4), optional, intent(in) :: stype
real(r8), optional, intent(in) :: df(:)
integer(i4), optional, intent(in) :: nu(:)
#ifdef HAVE_PETSC
!---
integer(i4) :: i,smoother
!---
class(oft_ml_precond), pointer :: this_ml
type(oft_petsc_sjacobi_solver), pointer :: this_jac
type(oft_petsc_cg_solver), pointer :: this_cg
type(oft_petsc_gmres_solver), pointer :: this_gmres
type(oft_petsc_diagprecond), pointer :: this_diag
type(oft_petsc_direct_solver), pointer :: this_direct
DEBUG_STACK_PUSH
!---
smoother=1 ! Default to jacobi
IF(PRESENT(stype))smoother=ABS(stype)
SELECT CASE(smoother)
  CASE(1)
    IF(.NOT.(PRESENT(df).AND.PRESENT(nu)))THEN
      CALL oft_abort('Jacobi smoother requires (df, nu)!','create_petsc_mlpre',__FILE__)
    END IF
    IF(oft_debug_print(1))WRITE(*,*)'Creating MG-Jacobi',df,nu
  CASE(2)
    IF(.NOT.PRESENT(nu))THEN
      CALL oft_abort('GMRES smoother requires (nu)!','create_petsc_mlpre',__FILE__)
    END IF
    IF(oft_debug_print(1))WRITE(*,*)'Creating MG-GMRES',nu
  CASE(3)
    IF(.NOT.PRESENT(nu))THEN
      CALL oft_abort('CG smoother requires (nu)!','create_petsc_mlpre',__FILE__)
    END IF
    IF(oft_debug_print(1))WRITE(*,*)'Creating MG-CG',nu
  CASE DEFAULT
    CALL oft_abort('Unknown smoother type!','create_petsc_mlpre',__FILE__)
END SELECT
!---Set smoother
allocate(oft_ml_precond::pre)
SELECT TYPE(pre)
  TYPE IS(oft_ml_precond)
    this_ml=>pre
  CLASS DEFAULT
    CALL oft_abort('Error in precon allocation!','create_petsc_mlpre',__FILE__)
END SELECT
!---Setup MG smoother
this_ml%interp=>interp
this_ml%inject=>inject
this_ml%vec_create=>create_vec
this_ml%level=ABS(levels(nlevels))
IF(PRESENT(bc))this_ml%bc=>bc
!---Allocate top-level smoother
SELECT CASE(smoother)
  CASE(1) ! Jacobi
    ALLOCATE(oft_petsc_sjacobi_solver::this_ml%smooth_up)
    !---Setup top-level smoother
    IF(petsc_sjacobi_solver_cast(this_jac,this_ml%smooth_up)<0) &
      CALL oft_abort('Failed to allocate "this_jac".','create_petsc_mlpre',__FILE__)
    this_jac%A=>Mats(nlevels)%M
    this_jac%df=df(nlevels)
    this_jac%its=nu(nlevels)
    IF(PRESENT(bc))this_jac%bc=>bc
  CASE(2) ! FGMRES
    ALLOCATE(oft_petsc_gmres_solver::this_ml%smooth_up)
    !---Setup top-level smoother
    IF(petsc_gmres_solver_cast(this_gmres,this_ml%smooth_up)<0) &
      CALL oft_abort('Failed to allocate "this_gmres".','create_petsc_mlpre',__FILE__)
    this_gmres%A=>Mats(nlevels)%M
    this_gmres%nrits=nu(nlevels)
    this_gmres%its=nu(nlevels)
    !---
    ALLOCATE(oft_petsc_diagprecond::this_gmres%pre)
  CASE(3) ! CG
    ALLOCATE(oft_petsc_cg_solver::this_ml%smooth_up)
    !---Setup top-level smoother
    IF(petsc_cg_solver_cast(this_cg,this_ml%smooth_up)<0) &
      CALL oft_abort('Failed to allocate "this_cg".','create_petsc_mlpre',__FILE__)
    this_cg%A=>Mats(nlevels)%M
    this_cg%its=nu(nlevels)
    !---
    ALLOCATE(oft_petsc_diagprecond::this_cg%pre)
    IF(petsc_diagprecond_cast(this_diag,this_cg%pre)<0) &
      CALL oft_abort('Failed to allocate "this_diag".','create_petsc_mlpre',__FILE__)
END SELECT
this_ml%smooth_down=>this_ml%smooth_up
this_ml%symmetric=.TRUE.
!---Build ML preconditioner
DO i=nlevels-1,1,-1
  !---Set smoother
  ALLOCATE(oft_ml_precond::this_ml%base_solve)
  IF(ml_precond_cast(this_ml,this_ml%base_solve)<0) &
    CALL oft_abort('Failed to allocate "this_ml".','create_petsc_mlpre',__FILE__)
  !---Setup MG smoother
  this_ml%interp=>interp
  this_ml%inject=>inject
  this_ml%vec_create=>create_vec
  this_ml%level=ABS(levels(i))
  IF(levels(i)<0)CYCLE
  IF(PRESENT(bc))this_ml%bc=>bc
  IF(i==1.AND.stype<0)THEN
    ALLOCATE(oft_petsc_direct_solver::this_ml%smooth_up)
    !---Setup top-level smoother
    IF(petsc_direct_solver_cast(this_direct,this_ml%smooth_up)<0) &
      CALL oft_abort('Failed to allocate "this_direct".','create_petsc_mlpre',__FILE__)
    this_direct%A=>Mats(i)%M
    this_direct%its=nu(i)
  ELSE
  !---Allocate level smoother
  SELECT CASE(smoother)
    CASE(1)
      ALLOCATE(oft_petsc_sjacobi_solver::this_ml%smooth_up)
      !---Setup top-level smoother
      IF(petsc_sjacobi_solver_cast(this_jac,this_ml%smooth_up)<0) &
        CALL oft_abort('Failed to allocate "this_jac".','create_petsc_mlpre',__FILE__)
      this_jac%A=>Mats(i)%M
      this_jac%df=df(i)
      this_jac%its=nu(i)
      this_ml%smooth_down=>this_ml%smooth_up
      this_ml%symmetric=.TRUE.
      IF(PRESENT(bc))this_jac%bc=>bc
    CASE(2) ! FGMRES
      ALLOCATE(oft_petsc_gmres_solver::this_ml%smooth_up)
      !---Setup top-level smoother
      IF(petsc_gmres_solver_cast(this_gmres,this_ml%smooth_up)<0) &
        CALL oft_abort('Failed to allocate "this_gmres".','create_petsc_mlpre',__FILE__)
      this_gmres%A=>Mats(i)%M
      this_gmres%nrits=MIN(20,nu(i))
      this_gmres%its=nu(i)
      !---
      ALLOCATE(oft_petsc_diagprecond::this_gmres%pre)
      IF(i>1)THEN
        this_ml%smooth_down=>this_ml%smooth_up
        this_ml%symmetric=.TRUE.
      END IF
    CASE(3) ! CG
      ALLOCATE(oft_petsc_cg_solver::this_ml%smooth_up)
      !---Setup top-level smoother
      IF(petsc_cg_solver_cast(this_cg,this_ml%smooth_up)<0) &
        CALL oft_abort('Failed to allocate "this_cg".','create_petsc_mlpre',__FILE__)
      this_cg%A=>Mats(i)%M
      this_cg%its=nu(i)
      !---
      ALLOCATE(oft_petsc_diagprecond::this_cg%pre)
      IF(petsc_diagprecond_cast(this_diag,this_cg%pre)<0) &
        CALL oft_abort('Failed to allocate "this_diag".','create_petsc_mlpre',__FILE__)
      IF(i>1)THEN
        this_ml%smooth_down=>this_ml%smooth_up
        this_ml%symmetric=.TRUE.
      END IF
  END SELECT
  END IF
END DO
DEBUG_STACK_POP
#else
CALL oft_abort("Not compiled with PETSc", "create_petsc_mlpre", __FILE__)
#endif
end subroutine create_petsc_mlpre
!---------------------------------------------------------------------------
! SUBROUTINE: create_native_solver
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
SUBROUTINE create_native_solver(solver,solver_type)
CLASS(oft_solver), POINTER, INTENT(out) :: solver
CHARACTER(LEN=*), INTENT(in) :: solver_type
DEBUG_STACK_PUSH
SELECT CASE(TRIM(solver_type))
  CASE("cg")
    ALLOCATE(oft_native_cg_solver::solver)
  CASE("gmres")
    ALLOCATE(oft_native_gmres_solver::solver)
  CASE("sjacobi")
    ALLOCATE(oft_jblock_precond::solver)
  CASE("lu")
    ALLOCATE(oft_lusolver::solver)
  CASE("ilu")
    ALLOCATE(oft_ilusolver::solver)
  CASE DEFAULT
    CALL oft_abort("Invalid solver type.","create_native_solver",__FILE__)
END SELECT
DEBUG_STACK_POP
end subroutine create_native_solver
!---------------------------------------------------------------------------
! SUBROUTINE: create_petsc_solver
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
SUBROUTINE create_petsc_solver(solver,solver_type)
CLASS(oft_solver), POINTER, INTENT(out) :: solver
CHARACTER(LEN=*), INTENT(in) :: solver_type
#ifdef HAVE_PETSC
DEBUG_STACK_PUSH
SELECT CASE(TRIM(solver_type))
  CASE("cg")
    ALLOCATE(oft_petsc_cg_solver::solver)
  CASE("gmres")
    ALLOCATE(oft_petsc_gmres_solver::solver)
  CASE("sjacobi")
    ALLOCATE(oft_petsc_sjacobi_solver::solver)
  CASE("lu")
    ALLOCATE(oft_petsc_direct_solver::solver)
  CASE("ilu")
    CALL oft_abort("ILU0 not currently supported via PETSc.","create_petsc_solver",__FILE__)
  CASE DEFAULT
    CALL oft_abort("Invalid solver type.","create_petsc_solver",__FILE__)
END SELECT
DEBUG_STACK_POP
#else
  CALL oft_abort("Not compiled with PETSc", "create_petsc_solver", __FILE__)
#endif
end subroutine create_petsc_solver
!---------------------------------------------------------------------------
! SUBROUTINE: create_cg_solver
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
SUBROUTINE create_cg_solver(solver,force_native)
CLASS(oft_solver), POINTER, INTENT(out) :: solver
LOGICAL, OPTIONAL, INTENT(in) :: force_native
LOGICAL :: native_solver
DEBUG_STACK_PUSH
native_solver=.FALSE.
IF(PRESENT(force_native))native_solver=force_native
IF(use_petsc.AND.(.NOT.native_solver))THEN
  CALL create_petsc_solver(solver,"cg")
ELSE
  CALL create_native_solver(solver,"cg")
END IF
DEBUG_STACK_POP
end subroutine create_cg_solver
!---------------------------------------------------------------------------
! SUBROUTINE: create_cg_solver
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
SUBROUTINE create_gmres_solver(solver,nrits,force_native)
CLASS(oft_solver), POINTER, INTENT(out) :: solver
INTEGER(i4), OPTIONAL, INTENT(in) :: nrits
LOGICAL, OPTIONAL, INTENT(in) :: force_native
LOGICAL :: native_solver
DEBUG_STACK_PUSH
native_solver=.FALSE.
IF(PRESENT(force_native))native_solver=force_native
IF(use_petsc.AND.(.NOT.native_solver))THEN
#ifdef HAVE_PETSC
  ALLOCATE(oft_petsc_gmres_solver::solver)
  IF(PRESENT(nrits))THEN
    SELECT TYPE(this=>solver)
      CLASS IS(oft_petsc_gmres_solver)
        this%nrits=nrits
    END SELECT
  END IF
#endif
ELSE
  ALLOCATE(oft_native_gmres_solver::solver)
  IF(PRESENT(nrits))THEN
    SELECT TYPE(this=>solver)
      CLASS IS(oft_native_gmres_solver)
        this%nrits=nrits
    END SELECT
  END IF
END IF
DEBUG_STACK_POP
end subroutine create_gmres_solver
!---------------------------------------------------------------------------
! SUBROUTINE: create_solver_xml
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
RECURSIVE SUBROUTINE create_solver_xml(solver,solver_node,level)
CLASS(oft_solver), POINTER, INTENT(out) :: solver
TYPE(xml_node), POINTER, INTENT(in) :: solver_node
INTEGER(i4), OPTIONAL, INTENT(in) :: level
#ifdef HAVE_XML
!---
INTEGER(i4) :: nread,nnodes,ierr
TYPE(xml_node), POINTER :: pre_node
!---
integer(i4) :: i,val_level
CHARACTER(LEN=20) :: solver_type,temp_string
LOGICAL :: force_native,native_solver,petsc_solver
DEBUG_STACK_PUSH
val_level=1
IF(PRESENT(level))val_level=level
native_solver=.FALSE.
!---
CALL xml_extractDataAttribute(solver_node,"type",solver_type,iostat=ierr)
IF(oft_debug_print(2))WRITE(*,*)'Found solver: ',solver_type
force_native=.FALSE.
IF(xml_hasAttribute(solver_node,"native"))THEN
  CALL xml_extractDataAttribute(solver_node,"native",temp_string,iostat=ierr)
  force_native=((temp_string(1:1)=='t').OR.(temp_string(1:1)=='T'))
END IF
!---
petsc_solver=use_petsc.AND.(.NOT.force_native)
SELECT CASE(TRIM(solver_type))
  CASE("cg")
    IF(petsc_solver)THEN
      CALL create_petsc_solver(solver,"cg")
    ELSE
      CALL create_native_solver(solver,"cg")
      native_solver=.TRUE.
    END IF
  CASE("gmres")
    IF(petsc_solver)THEN
      CALL create_petsc_solver(solver,"gmres")
    ELSE
      CALL create_native_solver(solver,"gmres")
      native_solver=.TRUE.
    END IF
  CASE("sjacobi")
    IF(petsc_solver)THEN
      CALL create_petsc_solver(solver,"sjacobi")
    ELSE
      CALL create_native_solver(solver,"sjacobi")
      native_solver=.TRUE.
    END IF
  CASE("lu")
    IF(use_petsc)THEN
      CALL create_petsc_solver(solver,"lu")
    ELSE
      CALL create_native_solver(solver,"lu")
      native_solver=.TRUE.
    END IF
  CASE("ilu")
    IF(use_petsc)THEN
      CALL create_petsc_solver(solver,"ilu")
    ELSE
      CALL create_native_solver(solver,"ilu")
      native_solver=.TRUE.
    END IF
  CASE DEFAULT
    CALL oft_abort("Invalid solver type.","create_solver_xml",__FILE__)
END SELECT
!---
CALL solver%setup_from_xml(solver_node,val_level)
!---
CALL xml_get_element(solver_node,"pre",pre_node,ierr)
IF(ierr==0)THEN
  CALL create_pre_xml(solver%pre,pre_node,native_solver,val_level)
END IF
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','create_solver_xml',__FILE__)
#endif
end subroutine create_solver_xml
!---------------------------------------------------------------------------
! SUBROUTINE: create_native_pre
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
SUBROUTINE create_native_pre(pre,pre_type)
CLASS(oft_solver), POINTER, INTENT(out) :: pre
CHARACTER(LEN=*), INTENT(in) :: pre_type
DEBUG_STACK_PUSH
SELECT CASE(TRIM(pre_type))
  CASE("jacobi")
    ALLOCATE(oft_diag_scale::pre)
  CASE("block_jacobi")
    ALLOCATE(oft_bjprecond::pre)
  CASE("lu")
    ALLOCATE(oft_lusolver::pre)
  CASE("ilu")
    ALLOCATE(oft_ilusolver::pre)
  CASE DEFAULT
    CALL oft_abort("Invalid preconditioner type.","create_native_pre",__FILE__)
END SELECT
DEBUG_STACK_POP
end subroutine create_native_pre
!---------------------------------------------------------------------------
! SUBROUTINE: create_petsc_pre
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
SUBROUTINE create_petsc_pre(pre,pre_type)
CLASS(oft_solver), POINTER, INTENT(out) :: pre
CHARACTER(LEN=*), INTENT(in) :: pre_type
#ifdef HAVE_PETSC
DEBUG_STACK_PUSH
SELECT CASE(TRIM(pre_type))
  CASE("jacobi")
    ALLOCATE(oft_petsc_diagprecond::pre)
  CASE("block_jacobi")
    ALLOCATE(oft_petsc_asprecond::pre)
    SELECT TYPE(this=>pre)
      CLASS IS(oft_petsc_asprecond)
        this%overlap=0
    END SELECT
  CASE("add_schwarz")
    ALLOCATE(oft_petsc_asprecond::pre)
  CASE("lu")
    ALLOCATE(oft_petsc_luprecond::pre)
  CASE("ilu")
    CALL oft_abort("ILU0 not currently supported via PETSc.","create_petsc_pre",__FILE__)
  !   ALLOCATE(oft_petsc_luprecond::pre)
  CASE DEFAULT
    CALL oft_abort("Invalid preconditioner type.","create_petsc_pre",__FILE__)
END SELECT
DEBUG_STACK_POP
#else
  CALL oft_abort("Not compiled with PETSc", "create_petsc_pre", __FILE__)
#endif
end subroutine create_petsc_pre
!---------------------------------------------------------------------------
! SUBROUTINE: create_diag_pre
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
SUBROUTINE create_diag_pre(pre)
CLASS(oft_solver), POINTER, INTENT(out) :: pre
LOGICAL :: native_solver
DEBUG_STACK_PUSH
IF(use_petsc)THEN
  CALL create_petsc_pre(pre,"jacobi")
ELSE
  CALL create_native_pre(pre,"jacobi")
END IF
DEBUG_STACK_POP
end subroutine create_diag_pre
!---------------------------------------------------------------------------
!> Create ILU(0) preconditioner (native or MKL)
!---------------------------------------------------------------------------
SUBROUTINE create_ilu_pre(pre)
CLASS(oft_solver), POINTER, INTENT(out) :: pre
LOGICAL :: native_solver
DEBUG_STACK_PUSH
IF(use_petsc)THEN
  CALL create_petsc_pre(pre,"ilu")
ELSE
  CALL create_native_pre(pre,"ilu")
END IF
DEBUG_STACK_POP
end subroutine create_ilu_pre
!---------------------------------------------------------------------------
! SUBROUTINE: create_bjacobi_pre
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
SUBROUTINE create_bjacobi_pre(pre,nlocal)
CLASS(oft_solver), POINTER, INTENT(out) :: pre
INTEGER(i4), INTENT(in) :: nlocal
DEBUG_STACK_PUSH
IF(use_petsc)THEN
#ifdef HAVE_PETSC
  ALLOCATE(oft_petsc_asprecond::pre)
  SELECT TYPE(this=>pre)
    CLASS IS(oft_petsc_asprecond)
      this%overlap=0
      this%n_local=nlocal
  END SELECT
#endif
ELSE
  ALLOCATE(oft_bjprecond::pre)
  SELECT TYPE(this=>pre)
    CLASS IS(oft_bjprecond)
      this%nlocal=nlocal
  END SELECT
  ALLOCATE(oft_lusolver::pre%pre)
END IF
DEBUG_STACK_POP
end subroutine create_bjacobi_pre
!---------------------------------------------------------------------------
! SUBROUTINE: create_pre_xml
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
RECURSIVE SUBROUTINE create_pre_xml(pre,pre_node,native_solver,level)
CLASS(oft_solver), POINTER, INTENT(out) :: pre
TYPE(xml_node), POINTER, INTENT(in) :: pre_node
LOGICAL, INTENT(in) :: native_solver
INTEGER(i4), OPTIONAL, INTENT(in) :: level
#ifdef HAVE_XML
!---
INTEGER(i4) :: nread,nnodes,ierr
TYPE(xml_node), POINTER :: solver_node
!---
integer(i4) :: i,val_level,smoother
logical :: switch
CHARACTER(LEN=20) :: pre_type,temp_string
DEBUG_STACK_PUSH
val_level=1
IF(PRESENT(level))val_level=level
!---
CALL xml_extractDataAttribute(pre_node,"type",pre_type,iostat=ierr)
IF(oft_debug_print(2))WRITE(*,*)'Found preconditioner: ',pre_type
!---
SELECT CASE(TRIM(pre_type))
  CASE("jacobi")
    IF(use_petsc)THEN
      CALL create_petsc_pre(pre,"jacobi")
    ELSE
      CALL create_native_pre(pre,"jacobi")
    END IF
  CASE("block_jacobi")
    IF(use_petsc)THEN
      CALL create_petsc_pre(pre,"block_jacobi")
    ELSE
      CALL create_native_pre(pre,"block_jacobi")
    END IF
  CASE("add_schwarz")
    IF(use_petsc)THEN
      CALL create_petsc_pre(pre,"add_schwarz")
    ELSE
      CALL oft_abort("Additive-Schwarz precon requires PETSc interface","create_pre_xml",__FILE__)
    END IF
  CASE("lu")
    IF(use_petsc)THEN
      IF(native_solver)THEN
        CALL oft_abort("LU precon requires PETSc parent solver.","create_pre_xml",__FILE__)
      END IF
      CALL create_petsc_pre(pre,"lu")
    ELSE
      CALL create_native_pre(pre,"lu")
    END IF
  CASE("ilu")
    IF(use_petsc)THEN
      IF(native_solver)THEN
        CALL oft_abort("ILU precon requires PETSc parent solver.","create_pre_xml",__FILE__)
      END IF
      CALL create_petsc_pre(pre,"ilu")
    ELSE
      CALL create_native_pre(pre,"ilu")
    END IF
  CASE DEFAULT
    CALL oft_abort("Invalid precon type.","create_pre_xml",__FILE__)
END SELECT
!---
CALL pre%setup_from_xml(pre_node,val_level)
!---
CALL xml_get_element(pre_node,"solver",solver_node,ierr)
IF(ierr==0)THEN
  CALL create_solver_xml(pre%pre,solver_node,val_level)
END IF
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','create_pre_xml',__FILE__)
#endif
end subroutine create_pre_xml
!---------------------------------------------------------------------------
! SUBROUTINE: create_ml_xml
!---------------------------------------------------------------------------
!> Construct PETSc Multi-Grid preconditioner using native mechanics
!!
!! @param[out] pre Preconditioner object
!! @param[in] Mats Operator matrices [nlevels]
!! @param[in] levels List of level indices [nlevels]
!! @param[in] nlevels Number of levels
!! @param[in] create_vec Vector creation subroutine
!! @param[in] interp Interpolation subroutine
!! @param[in] inject Restriction subroutine
!! @param[in] pre_node
!---------------------------------------------------------------------------
subroutine create_ml_xml(pre,Mats,levels,nlevels,create_vec,interp,inject,pre_node,bc)
class(oft_solver), pointer, intent(out) :: pre
TYPE(oft_matrix_ptr), INTENT(in) :: Mats(:)
integer(i4), intent(in) :: levels(:)
integer(i4), intent(in) :: nlevels
procedure(oft_veccreate_proto) :: create_vec
procedure(oft_interp_proto) :: interp
procedure(oft_interp_proto) :: inject
TYPE(xml_node), POINTER, INTENT(in) :: pre_node
procedure(oft_bc_proto), optional :: bc
#ifdef HAVE_XML
!---
integer(i4) :: i,ierr,nnodes
class(oft_ml_precond), pointer :: this_ml
LOGICAL :: symmetric,up_present,down_present,coarse_present
CHARACTER(LEN=20) :: dir_type
TYPE(xml_node), POINTER :: up_node,down_node,coarse_node,current_node,solver_node
TYPE(xml_nodelist) :: current_nodes
DEBUG_STACK_PUSH
!---
symmetric=.FALSE.
up_present=.FALSE.
down_present=.FALSE.
coarse_present=.FALSE.
IF(oft_debug_print(1))WRITE(*,*)'Creating MG smoother'
!---
CALL xml_get_element(pre_node,"smoother",current_nodes,ierr)
IF(current_nodes%n==0)CALL oft_abort("Object contains no smoother definitions.","create_ml_xml",__FILE__)
DO i=1,current_nodes%n
  current_node=>current_nodes%nodes(i)%this
  CALL xml_extractDataAttribute(current_node,"direction",dir_type,iostat=ierr)
  IF(oft_debug_print(2))WRITE(*,*)'Found smoother: ',dir_type
  SELECT CASE(TRIM(dir_type))
    CASE("up")
      up_present=.TRUE.
      CALL xml_get_element(current_node,"solver",up_node,ierr)
    CASE("down")
      down_present=.TRUE.
      CALL xml_get_element(current_node,"solver",down_node,ierr)
    CASE("both")
      symmetric=.TRUE.
      up_present=.TRUE.
      CALL xml_get_element(current_node,"solver",up_node,ierr)
      EXIT
    CASE DEFAULT
      CALL oft_abort("Invalid smoother direction.","create_ml_xml",__FILE__)
  END SELECT
END DO
IF(ASSOCIATED(current_nodes%nodes))DEALLOCATE(current_nodes%nodes)
!---
CALL xml_get_element(pre_node,"coarse",current_node,ierr)
IF(ierr==0)THEN
  IF(oft_debug_print(2))WRITE(*,*)'Found coarse solver'
  coarse_present=.TRUE.
  CALL xml_get_element(current_node,"solver",coarse_node,ierr)
END IF
!---Set smoother
ALLOCATE(oft_ml_precond::pre)
SELECT TYPE(pre)
  TYPE IS(oft_ml_precond)
    this_ml=>pre
  CLASS DEFAULT
    CALL oft_abort('Error in precon allocation!','create_petsc_mlpre',__FILE__)
END SELECT
!---Setup MG smoother
this_ml%interp=>interp
this_ml%inject=>inject
this_ml%vec_create=>create_vec
this_ml%level=ABS(levels(nlevels))
IF(PRESENT(bc))this_ml%bc=>bc
!---Allocate top-level smoother
IF(up_present)THEN
  CALL create_solver_xml(this_ml%smooth_up,up_node,this_ml%level)
  this_ml%smooth_up%A=>Mats(nlevels)%M
  IF(PRESENT(bc))this_ml%smooth_up%bc=>bc
  IF(ASSOCIATED(this_ml%smooth_up%pre))this_ml%smooth_up%pre%A=>Mats(nlevels)%M
END IF
IF(down_present)THEN
  CALL create_solver_xml(this_ml%smooth_down,down_node,this_ml%level)
  this_ml%smooth_down%A=>Mats(nlevels)%M
  IF(PRESENT(bc))this_ml%smooth_down%bc=>bc
  IF(ASSOCIATED(this_ml%smooth_down%pre))this_ml%smooth_down%pre%A=>Mats(nlevels)%M
END IF
IF(symmetric)THEN
  this_ml%symmetric=.TRUE.
  this_ml%smooth_down=>this_ml%smooth_up
END IF
!---Build ML preconditioner
DO i=nlevels-1,1,-1
  !---Set smoother
  ALLOCATE(oft_ml_precond::this_ml%base_solve)
  IF(ml_precond_cast(this_ml,this_ml%base_solve)<0) &
    CALL oft_abort('Failed to allocate "this_ml".','create_petsc_mlpre',__FILE__)
  !---Setup MG smoother
  this_ml%interp=>interp
  this_ml%inject=>inject
  this_ml%vec_create=>create_vec
  this_ml%level=ABS(levels(i))
  IF(levels(i)<0)CYCLE
  IF(i==1.AND.coarse_present)THEN
    CALL create_solver_xml(this_ml%smooth_up,coarse_node,this_ml%level)
    this_ml%smooth_up%A=>Mats(i)%M
    IF(PRESENT(bc))this_ml%smooth_up%bc=>bc
    IF(ASSOCIATED(this_ml%smooth_up%pre))this_ml%smooth_up%pre%A=>Mats(i)%M
  ELSE
    IF(up_present)THEN
      CALL create_solver_xml(this_ml%smooth_up,up_node,this_ml%level)
      this_ml%smooth_up%A=>Mats(i)%M
      IF(PRESENT(bc))this_ml%smooth_up%bc=>bc
      IF(ASSOCIATED(this_ml%smooth_up%pre))this_ml%smooth_up%pre%A=>Mats(i)%M
    END IF
    IF(down_present)THEN
      CALL create_solver_xml(this_ml%smooth_down,down_node,this_ml%level)
      this_ml%smooth_down%A=>Mats(i)%M
      IF(PRESENT(bc))this_ml%smooth_down%bc=>bc
      IF(ASSOCIATED(this_ml%smooth_down%pre))this_ml%smooth_down%pre%A=>Mats(i)%M
    END IF
    IF(symmetric)THEN
      this_ml%symmetric=.TRUE.
      this_ml%smooth_down=>this_ml%smooth_up
    END IF
  END IF
END DO
DEBUG_STACK_POP
#else
CALL oft_abort('OFT not compiled with xml support.','create_ml_xml',__FILE__)
#endif
end subroutine create_ml_xml
END MODULE oft_solver_utils
