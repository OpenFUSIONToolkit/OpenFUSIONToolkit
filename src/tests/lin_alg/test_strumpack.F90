!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file test_strumpack.F90
!
!> Regression tests for direct solves of sparse matrices using STRUMPACK.
!!
!! A sparse, non-symmetric matrix from a 2D convection-diffusion stencil is
!! factored and solved, then modified and refactored three times and the
!! solutions are compared against known exact solutions:
!! - Same non-zero pattern (`update(.FALSE.)`)
!! - Full factorization (`update(.TRUE.)`)
!! - Changed non-zero pattern with `update(.FALSE.)`, which must be detected
!!   and trigger a full factorization
!!
!! Solver options are set directly or read from XML (`use_xml`), in which case
!! the parsed options are also written to the results file.
!!
!! @authors Chris Hansen
!! @date October 2026
!! @ingroup testing
!------------------------------------------------------------------------------
PROGRAM test_strumpack
USE oft_base
USE oft_la_base, ONLY: oft_vector
USE oft_native_la, ONLY: oft_native_vector, oft_native_matrix
USE oft_lu, ONLY: oft_lusolver
IMPLICIT NONE
INTEGER(i4) :: io_unit,ierr,nrows
INTEGER(i4) :: nsize = 60
INTEGER(i4) :: comp_type = -1
INTEGER(i4) :: leaf_size = -1
INTEGER(i4) :: min_sep_size = -1
REAL(r8) :: comp_rel_tol = -1.d0
LOGICAL :: use_xml = .FALSE.
NAMELIST/test_strumpack_options/nsize,comp_type,leaf_size,min_sep_size, &
  comp_rel_tol,use_xml
!---Initialize enviroment
CALL oft_init
!---Read in options
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,test_strumpack_options,IOSTAT=ierr)
CLOSE(io_unit)
!---Check available solvers
#if !defined( HAVE_STRUMPACK )
WRITE(*,*)'SKIP TEST'
CALL oft_finalize
#endif
CALL test_solve()
!---Finalize enviroment
CALL oft_finalize
CONTAINS
!---------------------------------------------------------------------------------
!> Factor and solve a system, then refactor after modifying the matrix
!! and output the relative error of each solution.
!---------------------------------------------------------------------------------
SUBROUTINE test_solve()
INTEGER(i4) :: i
REAL(r8) :: err(4)
REAL(r8), ALLOCATABLE, DIMENSION(:) :: x_exact
CLASS(oft_vector), POINTER :: u,g
TYPE(oft_native_matrix), TARGET :: mat
TYPE(oft_lusolver) :: linv
TYPE(xml_node) :: solver_node
!---Build matrix
CALL build_sparse(mat,5)
nrows=mat%nr
CALL create_local_vec(u)
CALL create_local_vec(g)
!---Setup solver
IF(use_xml)THEN
  CALL xml_get_element(oft_env%xml,"solver",solver_node,ierr)
  IF(ierr/=0)CALL oft_abort('Could not find XML node.','test_solve',__FILE__)
  CALL linv%setup_from_xml(solver_node)
ELSE
  linv%package='strumpack'
#ifdef HAVE_STRUMPACK
  linv%strumpack_struct%comp_type=comp_type
  linv%strumpack_struct%leaf_size=leaf_size
  linv%strumpack_struct%min_sep_size=min_sep_size
  linv%strumpack_struct%comp_rel_tol=comp_rel_tol
#endif
END IF
linv%A=>mat
!---Initial solve
ALLOCATE(x_exact(nrows))
DO i=1,nrows
  x_exact(i)=SIN(3.d0*REAL(i,8)/REAL(nrows,8))
END DO
err(1)=check_solve(linv,mat,u,g,x_exact)
!---Modify matrix and refactor (same pattern)
CALL shift_diagonal(mat,0.5d0)
CALL linv%update(.FALSE.)
err(2)=check_solve(linv,mat,u,g,x_exact)
!---Modify matrix and refactor (full factorization)
CALL shift_diagonal(mat,0.5d0)
CALL linv%update(.TRUE.)
err(3)=check_solve(linv,mat,u,g,x_exact)
!---Change non-zero pattern but request same-pattern refactor (must be detected)
CALL mat%delete()
CALL build_sparse(mat,9)
CALL linv%update(.FALSE.)
err(4)=check_solve(linv,mat,u,g,x_exact)
!---Cleanup
CALL linv%delete()
DEALLOCATE(x_exact)
CALL mat%delete()
CALL u%delete()
CALL g%delete()
DEALLOCATE(u,g)
!---Output results
WRITE(*,*)'Relative errors =',err
IF(oft_env%head_proc)THEN
  OPEN(NEWUNIT=io_unit,FILE='strumpack.results')
  DO i=1,4
    WRITE(io_unit,*)err(i)
  END DO
#ifdef HAVE_STRUMPACK
  IF(use_xml)THEN
    WRITE(io_unit,*)linv%package
    WRITE(io_unit,*)linv%strumpack_struct%comp_type
    WRITE(io_unit,*)linv%strumpack_struct%comp_rel_tol
    WRITE(io_unit,*)linv%strumpack_struct%leaf_size
    WRITE(io_unit,*)linv%strumpack_struct%min_sep_size
  END IF
#endif
  CLOSE(io_unit)
END IF
END SUBROUTINE test_solve
!---------------------------------------------------------------------------------
!> Build sparse matrix from a 5-point convection-diffusion stencil on a
!! `nsize` x `nsize` grid, optionally with diagonal neighbors (9-point stencil)
!---------------------------------------------------------------------------------
SUBROUTINE build_sparse(mat,npts)
TYPE(oft_native_matrix), INTENT(inout) :: mat
INTEGER(i4), INTENT(in) :: npts
INTEGER(i4) :: i,j,k,row,nnz
INTEGER(i4) :: di(9) = [0,-1,1,0,0,-1,1,-1,1], dj(9) = [0,0,0,-1,1,-1,-1,1,1]
REAL(r8) :: stencil(9) = [4.d0,-1.2d0,-0.8d0,-1.1d0,-0.9d0,-0.1d0,-0.05d0,-0.05d0,-0.1d0]
mat%nr=nsize**2; mat%nc=nsize**2
mat%nrg=mat%nr; mat%ncg=mat%nc
ALLOCATE(mat%kr(mat%nr+1),mat%lc(npts*mat%nr),mat%M(npts*mat%nr))
nnz=0
mat%kr(1)=1
DO j=1,nsize
  DO i=1,nsize
    row=(j-1)*nsize+i
    DO k=1,npts
      IF((i+di(k)<1).OR.(i+di(k)>nsize).OR.(j+dj(k)<1).OR.(j+dj(k)>nsize))CYCLE
      nnz=nnz+1
      mat%lc(nnz)=row+di(k)+dj(k)*nsize
      mat%M(nnz)=stencil(k)
      IF((k==1).AND.(npts>5))mat%M(nnz)=mat%M(nnz)+1.d0 ! Keep diagonally dominant
    END DO
    mat%kr(row+1)=nnz+1
  END DO
END DO
mat%nnz=nnz
END SUBROUTINE build_sparse
!---------------------------------------------------------------------------------
!> Add `shift` to the diagonal of the matrix
!---------------------------------------------------------------------------------
SUBROUTINE shift_diagonal(mat,shift)
TYPE(oft_native_matrix), INTENT(inout) :: mat
REAL(r8), INTENT(in) :: shift
INTEGER(i4) :: i,j
DO i=1,mat%nr
  DO j=mat%kr(i),mat%kr(i+1)-1
    IF(mat%lc(j)==i)mat%M(j)=mat%M(j)+shift
  END DO
END DO
END SUBROUTINE shift_diagonal
!---------------------------------------------------------------------------------
!> Solve system with RHS corresponding to `x_exact` and return relative error
!---------------------------------------------------------------------------------
FUNCTION check_solve(linv,mat,u,g,x_exact) RESULT(error)
TYPE(oft_lusolver), INTENT(inout) :: linv
TYPE(oft_native_matrix), INTENT(inout) :: mat
CLASS(oft_vector), INTENT(inout) :: u,g
REAL(r8), INTENT(in) :: x_exact(:)
REAL(r8) :: error
INTEGER(i4) :: i,j
REAL(r8), POINTER, DIMENSION(:) :: vals
!---Compute RHS (b = A*x_exact)
ALLOCATE(vals(nrows))
vals=0.d0
DO i=1,mat%nr
  DO j=mat%kr(i),mat%kr(i+1)-1
    vals(i)=vals(i)+mat%M(j)*x_exact(mat%lc(j))
  END DO
END DO
CALL g%restore_local(vals)
!---Solve and compare to exact solution
CALL u%set(0.d0)
CALL linv%apply(u,g)
CALL u%get_local(vals)
error=SQRT(SUM((vals-x_exact)**2))/SQRT(SUM(x_exact**2))
DEALLOCATE(vals)
END FUNCTION check_solve
!---------------------------------------------------------------------------------
!> Create a local (unstitched) native vector of size `nrows`
!---------------------------------------------------------------------------------
SUBROUTINE create_local_vec(vec)
CLASS(oft_vector), POINTER, INTENT(out) :: vec
ALLOCATE(oft_native_vector::vec)
SELECT TYPE(this=>vec)
  CLASS IS(oft_native_vector)
    this%n=nrows; this%ng=nrows
    this%nslice=this%n
    ALLOCATE(this%v(this%n))
    ALLOCATE(this%stitch_info)
    ALLOCATE(this%stitch_info%be(this%n))
    this%stitch_info%full=.TRUE.
    this%stitch_info%nbe=0
    this%stitch_info%be=.FALSE.
END SELECT
END SUBROUTINE create_local_vec
END PROGRAM test_strumpack
