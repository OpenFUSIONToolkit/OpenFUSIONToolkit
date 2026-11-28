!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_hcurl_grad.F90
!
!> Regression tests for vector H(Curl) + Grad(H^1) finite elements. Tests are performed
!! on a unit cube at different polynomial orders.
!!
!! The current test cases are:
!! - Solve the equation \f$ M B = M \hat{I} \f$
!! - L2 projection of an analytic field with vanishing tangential boundary
!!   component using the `'zerob'` mass matrix
!!
!! @authors Chris Hansen
!! @date April 2013
!! @ingroup testing
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Analytic vector field with vanishing tangential component on the boundary of
! the unit cube, used to test projection with the 'zerob' mass matrix
!------------------------------------------------------------------------------
MODULE test_hcurl_grad_helpers
USE oft_base
USE fem_utils, ONLY: fem_interp
IMPLICIT NONE
TYPE, EXTENDS(fem_interp) :: tang_zero_field
CONTAINS
  PROCEDURE :: interp => tang_zero_field_interp
END TYPE tang_zero_field
CONTAINS
SUBROUTINE tang_zero_field_interp(self,cell,f,gop,val)
CLASS(tang_zero_field), INTENT(inout) :: self
INTEGER(i4), INTENT(in) :: cell
REAL(r8), INTENT(in) :: f(:)
REAL(r8), INTENT(in) :: gop(3,4)
REAL(r8), INTENT(out) :: val(:)
REAL(r8) :: pt(3),s(3)
pt=self%mesh%log2phys(cell,f)
s=SIN(pi*pt)
val(1)=s(2)*s(3)*(1.d0+pt(1))
val(2)=s(1)*s(3)*(1.d0+pt(2)**2)
val(3)=s(1)*s(2)*COS(pt(3))
END SUBROUTINE tang_zero_field_interp
END MODULE test_hcurl_grad_helpers
program test_hcurl_grad
USE oft_base
USE oft_mesh_cube, ONLY: mesh_cube_id
USE oft_mesh_type, ONLY: cell_is_curved
! USE oft_io
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---
USE oft_la_base, ONLY: oft_vector, oft_matrix, oft_matrix_ptr
USE oft_solver_base, ONLY: oft_solver
USE oft_native_solvers, ONLY: oft_native_gmres_solver
USE oft_solver_utils, ONLY: create_native_mlpre, create_cg_solver, create_diag_pre
!---
USE fem_base, ONLY: oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
!---Grad(H^1) FE space
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp
!---H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup, oft_hcurl_fem
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_getmop, hcurl_grad_setup_interp, hcurl_grad_getmop_pre, hcurl_grad_mloptions, &
  oft_hcurl_grad_rinterp, oft_hcurl_grad_gzerop, oft_hcurl_grad_zerob, oft_hcurl_grad_project, oft_ml_hcurl_grad_vecspace
USE test_hcurl_grad_helpers, ONLY: tang_zero_field
IMPLICIT NONE
INTEGER(i4) :: order,ierr,io_unit,minlev
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_h1,ML_h1grad
TYPE(oft_ml_fem_type), TARGET :: ML_oft_hcurl
TYPE(oft_ml_fem_comp_type), TARGET :: ML_hcurl_grad
TYPE(oft_hcurl_grad_gzerop), TARGET :: hcurl_grad_gzerop
TYPE(oft_hcurl_grad_zerob), TARGET :: hcurl_grad_zerob
LOGICAL :: mg_test,zerob_test,adjoint_test
NAMELIST/test_hcurl_grad_options/order,mg_test,zerob_test,adjoint_test
!---Initialize enviroment
CALL oft_init
!---Read in options
mg_test=.FALSE.
zerob_test=.FALSE.
adjoint_test=.FALSE.
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_hcurl_grad_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Setup grid
CALL multigrid_construct(mg_mesh)
IF(mg_mesh%mesh%cad_type/=mesh_cube_id)CALL oft_abort('Wrong mesh type, test for CUBE only.','main',__FILE__)
!------------------------------------------------------------------------------
! Build FE structures
!------------------------------------------------------------------------------
minlev=2
IF(oft_env%nprocs>1)minlev=mg_mesh%nbase+1
IF(mg_mesh%mesh%type==3)minlev=mg_mesh%mgdim
!--- Grad(H^1) subspace
CALL oft_h1_setup(mg_mesh,order+1,ML_oft_h1,minlev=minlev)
IF(mg_test.OR.adjoint_test)CALL h1_setup_interp(ML_oft_h1)
!--- H(Curl) subspace
CALL oft_hcurl_setup(mg_mesh,order,ML_oft_hcurl,minlev=minlev)
IF(mg_test.OR.adjoint_test)CALL hcurl_setup_interp(ML_oft_hcurl)
!--- Full H(Curl) space
CALL oft_hcurl_grad_setup(ML_oft_hcurl,ML_oft_h1,ML_hcurl_grad,ML_h1grad,minlev=minlev)
hcurl_grad_gzerop%ML_hcurl_grad_rep=>ML_hcurl_grad
hcurl_grad_zerob%ML_hcurl_grad_rep=>ML_hcurl_grad
IF(mg_test.OR.adjoint_test)THEN
  CALL hcurl_grad_setup_interp(ML_hcurl_grad,ML_oft_h1,create_full=.TRUE.)
END IF
IF(mg_test)CALL hcurl_grad_mloptions
!---Run tests
oft_env%pm=.FALSE.
IF(adjoint_test)THEN
  CALL test_mg_adjoint
ELSE IF(zerob_test)THEN
  CALL test_mop_zerob
ELSE IF(mg_test)THEN
  CALL test_mopmg
ELSE
  CALL test_mop
END IF
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!---------------------------------------------------------------------------------
!> Solve the equation \f$ \nabla \times \nabla \times B = K \hat{I} \f$, where
!! \f$ K \f$ is the helicity matrix and \f$ \hat{I} \f$ is the identity vector
!! using H(Curl) elements.
!---------------------------------------------------------------------------------
SUBROUTINE test_mop
!---Create solver objects
CLASS(oft_solver), POINTER :: winv => NULL()
!---Local variables
REAL(r8) :: uu
CLASS(oft_vector), POINTER :: u,v
CLASS(oft_matrix), POINTER :: mop => NULL()
!---Set FE level
CALL ML_hcurl_grad%set_level(ML_hcurl_grad%nlevels,propogate=.TRUE.)
!---Create solver fields
CALL ML_hcurl_grad%vec_create(u)
CALL ML_hcurl_grad%vec_create(v)
!---Get FE operators
CALL hcurl_grad_getmop(ML_hcurl_grad%current_level,mop,'none')
!---Setup matrix solver
CALL create_cg_solver(winv)
winv%A=>mop
winv%its=-3
winv%rtol=1.d-10
CALL create_diag_pre(winv%pre)
!---Solve
CALL u%set(1.d0)
CALL hcurl_grad_gzerop%apply(u)
CALL mop%apply(u,v)
CALL u%set(0.d0)
CALL winv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='hcurl_grad.results')
  WRITE(io_unit,*)winv%cits
  WRITE(io_unit,*)uu/(u%ng-mg_mesh%mesh%global%np)
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
END SUBROUTINE test_mop
!---------------------------------------------------------------------------------
!> Solve the equation \f$ \nabla \times \nabla \times B = K \hat{I} \f$, where
!! \f$ K \f$ is the helicity matrix and \f$ \hat{I} \f$ is the identity vector
!! using H(Curl) elements.
!---------------------------------------------------------------------------------
SUBROUTINE test_mopmg
!---Create solver objects
! TYPE(oft_native_gmres_solver) :: winv
CLASS(oft_solver), POINTER :: winv => NULL()
CLASS(oft_solver), POINTER :: pre
TYPE(oft_matrix_ptr), POINTER :: mats(:) => NULL()
!---Local variables
INTEGER(i4) :: i,nlevels
REAL(r8) :: uu,r,rsq
REAL(r8), POINTER :: bvout(:,:)
REAL(r8), POINTER :: vals(:) => NULL()
CLASS(oft_vector), POINTER :: u,v,x,y
CLASS(oft_matrix), POINTER :: mop => NULL()
TYPE(oft_hcurl_grad_rinterp) :: field
!---Set FE level
CALL ML_hcurl_grad%set_level(ML_hcurl_grad%nlevels,propogate=.TRUE.)
!---Create solver fields
CALL ML_hcurl_grad%vec_create(u)
CALL ML_hcurl_grad%vec_create(v)
!---Setup matrix solver
nlevels=ML_hcurl_grad%nlevels-ML_hcurl_grad%minlev
CALL create_cg_solver(winv,force_native=.TRUE.)
winv%its=-3
winv%rtol=1.d-10
CALL hcurl_grad_getmop_pre(ML_hcurl_grad,winv%pre,mats,nlevels=nlevels)
mop=>mats(nlevels)%M
winv%A=>mop
!---Solve
CALL u%set(1.d0)
CALL hcurl_grad_gzerop%apply(u)
CALL mop%apply(u,v)
CALL u%set(0.d0)
CALL winv%apply(u,v)
uu=u%dot(u)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='hcurl_grad.results')
  WRITE(io_unit,*)winv%cits
  WRITE(io_unit,*)uu/(u%ng-mg_mesh%mesh%global%np)
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL u%delete
CALL v%delete
END SUBROUTINE test_mopmg
!---------------------------------------------------------------------------------
! L2 project an analytic field with vanishing tangential boundary component
! using the 'zerob' mass matrix and compute the relative projection error
!---------------------------------------------------------------------------------
SUBROUTINE test_mop_zerob
CLASS(oft_solver), POINTER :: winv => NULL()
INTEGER(i4) :: i,m
REAL(r8) :: res_norm,rhs_norm,err_sq,ref_sq,vol,det,goptmp(3,4),bh(3),bex(3)
LOGICAL :: curved
CLASS(oft_vector), POINTER :: u,v,w,rhs
CLASS(oft_matrix), POINTER :: mop => NULL()
TYPE(tang_zero_field) :: field_ex
TYPE(oft_hcurl_grad_rinterp) :: field_h
CLASS(oft_hcurl_fem), POINTER :: hcurl_rep
!---Set FE level
CALL ML_hcurl_grad%set_level(ML_hcurl_grad%nlevels,propogate=.TRUE.)
!---Create solver fields
CALL ML_hcurl_grad%vec_create(u)
CALL ML_hcurl_grad%vec_create(v)
CALL ML_hcurl_grad%vec_create(w)
CALL ML_hcurl_grad%vec_create(rhs)
!---Get FE operators
CALL hcurl_grad_getmop(ML_hcurl_grad%current_level,mop,'zerob')
!---Setup matrix solver
CALL create_cg_solver(winv)
winv%A=>mop
winv%its=-3
winv%rtol=1.d-12
CALL create_diag_pre(winv%pre)
!---Build RHS with boundary and vertex Grad(H^1) entries zeroed
field_ex%mesh=>mg_mesh%mesh
CALL oft_hcurl_grad_project(ML_hcurl_grad%current_level,field_ex,v)
CALL hcurl_grad_zerob%apply(v)
CALL hcurl_grad_gzerop%apply(v)
rhs_norm=v%norm(2)
CALL rhs%add(0.d0,1.d0,v)
!---Solve
CALL u%set(0.d0)
CALL winv%apply(u,v)
!---Check true residual
CALL mop%apply(u,w)
CALL w%add(1.d0,-1.d0,rhs)
res_norm=w%norm(2)/rhs_norm
!---Compute relative L2 projection error
field_h%u=>u
CALL field_h%setup(ML_hcurl_grad%current_level)
SELECT TYPE(this=>ML_oft_hcurl%current_level)
  CLASS IS(oft_hcurl_fem)
    hcurl_rep=>this
  CLASS DEFAULT
    CALL oft_abort('Invalid H(Curl) space','test_mop_zerob',__FILE__)
END SELECT
err_sq=0.d0; ref_sq=0.d0
DO i=1,mg_mesh%mesh%nc
  curved=cell_is_curved(mg_mesh%mesh,i)
  DO m=1,hcurl_rep%quad%np
    IF(curved.OR.m==1)CALL mg_mesh%mesh%jacobian(i,hcurl_rep%quad%pts(:,m),goptmp,vol)
    det=vol*hcurl_rep%quad%wts(m)
    CALL field_h%interp(i,hcurl_rep%quad%pts(:,m),goptmp,bh)
    CALL field_ex%interp(i,hcurl_rep%quad%pts(:,m),goptmp,bex)
    err_sq=err_sq+SUM((bh-bex)**2)*det
    ref_sq=ref_sq+SUM(bex**2)*det
  END DO
END DO
err_sq=oft_mpi_sum(err_sq)
ref_sq=oft_mpi_sum(ref_sq)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='hcurl_grad.results')
  WRITE(io_unit,*)winv%cits
  WRITE(io_unit,*)res_norm
  WRITE(io_unit,*)SQRT(err_sq/ref_sq)
  CLOSE(io_unit)
END IF
!---Destroy vectors
CALL field_h%delete
CALL u%delete
CALL v%delete
CALL w%delete
CALL rhs%delete
END SUBROUTINE test_mop_zerob
!---------------------------------------------------------------------------------
! Verify exact adjointness <y_f, P x_c> == <R y_f, x_c> across geometric MG level
!---------------------------------------------------------------------------------
SUBROUTINE test_mg_adjoint
TYPE(oft_ml_hcurl_grad_vecspace) :: ml_vs
CLASS(oft_vector), POINTER :: xc,yc,xf,yf,yf_ref
REAL(r8) :: lhs,rhs,rel_err,yf_err
REAL(r8), POINTER :: vals(:) => NULL()
INTEGER(i4) :: i
ml_vs%ML_FE_rep=>ML_hcurl_grad
!---Create coarse and fine test vectors across the top geometric level
CALL ML_hcurl_grad%set_level(ML_hcurl_grad%nlevels-1,propogate=.TRUE.)
CALL ML_hcurl_grad%vec_create(xc)
CALL ML_hcurl_grad%vec_create(yc)
CALL xc%get_local(vals)
DO i=1,SIZE(vals)
  vals(i)=SIN(REAL(i,r8)*0.37d0 + 0.1d0)
END DO
CALL xc%restore_local(vals,add=.TRUE.)
DEALLOCATE(vals)
CALL hcurl_grad_gzerop%apply(xc)
!---Fine level vectors
CALL ML_hcurl_grad%set_level(ML_hcurl_grad%nlevels,propogate=.TRUE.)
CALL ML_hcurl_grad%vec_create(xf)
CALL ML_hcurl_grad%vec_create(yf)
CALL ML_hcurl_grad%vec_create(yf_ref)
CALL yf%get_local(vals)
DO i=1,SIZE(vals)
  vals(i)=COS(REAL(i,r8)*0.53d0 + 0.2d0)
END DO
CALL yf%restore_local(vals,add=.TRUE.)
DEALLOCATE(vals)
CALL yf_ref%copy(yf)
!---Compute <y_f, interp(x_c)>
CALL ML_hcurl_grad%set_level(ML_hcurl_grad%nlevels-1,propogate=.TRUE.)
CALL ml_vs%interp(xc,xf)
lhs=yf%dot(xf)
!---Compute <inject(y_f), x_c> and verify y_f is unchanged
CALL ml_vs%inject(yf,yc)
rhs=yc%dot(xc)
rel_err=ABS(lhs-rhs)/MAX(ABS(lhs),ABS(rhs),1.d-30)
CALL yf_ref%add(1.d0,-1.d0,yf)
yf_err=yf_ref%norm(2)/MAX(yf%norm(2),1.d-30)
rel_err=MAX(rel_err,yf_err)
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='hcurl_grad.results')
  WRITE(io_unit,*)rel_err
  CLOSE(io_unit)
END IF
CALL xc%delete
CALL yc%delete
CALL xf%delete
CALL yf%delete
CALL yf_ref%delete
DEALLOCATE(xc,yc,xf,yf,yf_ref)
END SUBROUTINE test_mg_adjoint
end program test_hcurl_grad
