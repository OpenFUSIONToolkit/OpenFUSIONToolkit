!---------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!---------------------------------------------------------------------------
!> @file thin_wall_hodlr.F90
!
!> Hierarchical Off-Diagonal Low Rank matrix approximation functionality for ThinCurr
!!  - SVD compression
!!  - Adaptive Cross Approximation
!!
!! @authors Chris Hansen
!! @date Feb 2024
!! @ingroup doxy_oft_physics
!---------------------------------------------------------------------------
MODULE thin_wall_hodlr
USE thin_wall
IMPLICIT NONE
#include "local.h"
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type :: oft_tw_block
  INTEGER(4) :: np = 0
  INTEGER(4) :: nelems = 0
  INTEGER(4) :: ncells = 0
  INTEGER(4) :: parent = -1
  REAL(8) :: center(3) = 0.d0
  REAL(8) :: extent = 0.d0
  INTEGER(4), POINTER :: ielem(:) => NULL()
  INTEGER(4), POINTER :: inv_map(:) => NULL()
  INTEGER(4), POINTER :: ipts(:) => NULL()
  INTEGER(4), POINTER :: icell(:) => NULL()
end type oft_tw_block
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type oft_tw_level
  INTEGER(4) :: nblocks = 0
  INTEGER(4), POINTER, DIMENSION(:,:) :: mat_mask => NULL()
  TYPE(oft_tw_block), POINTER, DIMENSION(:) :: blocks => NULL()
end type oft_tw_level
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
type, extends(oft_noop_matrix) :: oft_tw_hodlr_op
  INTEGER(4) :: nlevels = 0 !< Number of levels in block heirarchy
  INTEGER(4) :: nblocks = 0 !< Number of blocks on the lowest level
  INTEGER(4) :: ndense = 0 !< Number of diagonal interactions
  INTEGER(4) :: nsparse = 0 !< Number of compressed off-diagonal interactions
  INTEGER(4) :: leaf_target = 1500 !< Target size for leaves on lowest level
  INTEGER(4) :: aca_min_its = 20 !< Minimum number of ACA+ iterations
  INTEGER(4) :: min_rank = 10 !< Minimum rank of compressed off-diagonal matrices
  REAL(8) :: L_svd_tol = -1.d0 !< SVD tolerance for inductance matrix
  REAL(8) :: L_aca_tol = -1.d0 !< ACA+ tolerance for inductance matrix
  REAL(8) :: B_svd_tol = -1.d0 !< SVD tolerance for B-field operators
  REAL(8) :: B_aca_tol = -1.d0 !< ACA+ tolerance for B-field operators
  INTEGER(4), POINTER, DIMENSION(:,:) :: dense_blocks => NULL() !< Indices of diagonal interactions
  INTEGER(4), POINTER, DIMENSION(:,:) :: sparse_blocks => NULL() !< Indices of compressed off-diagonal interactions
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:) :: aca_U_mats => NULL() !< U matrices for compressed off-diagonal interactions (L)
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:) :: aca_V_mats => NULL() !< V matrices for compressed off-diagonal interactions (L)
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:) :: aca_dense => NULL() !< Fallback dense matrices for compressed off-diagonal interactions (L)
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:) :: dense_mats => NULL() !< Dense matrices for diagonal interactions (L)
  TYPE(oft_native_dense_matrix) :: hole_Vcoil_mat !< Dense coupling matrix to holes and Vcoils (L)
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:,:) :: aca_BU_mats => NULL() !< U matrices for compressed off-diagonal interactions (B-field)
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:,:) :: aca_BV_mats => NULL() !< V matrices for compressed off-diagonal interactions (B-field)
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:,:) :: aca_B_dense => NULL() !< Fallback dense matrices for compressed off-diagonal interactions (B-field)
  TYPE(oft_native_dense_matrix), POINTER, DIMENSION(:,:) :: dense_B_mats => NULL() !< Dense matrices for diagonal interactions (B-field)
  REAL(8), POINTER, DIMENSION(:,:,:) :: hole_Vcoil_Bmat => NULL() !< Dense coupling matrix to holes and Vcoils (B-field)
  REAL(8), POINTER, DIMENSION(:,:,:) :: Icoil_Bmat => NULL() !< Dense coupling matrix to Icoils (B-fieldß)
  TYPE(oft_tw_block), POINTER, DIMENSION(:) :: blocks => NULL() !< REMOVE
  TYPE(oft_tw_level), POINTER, DIMENSION(:) :: levels => NULL() !< Block heirarchy 
  type(tw_type), pointer :: tw_obj => NULL()
contains
  !> Setup HODLR by performing partitioning and tagging block-block interactions
  procedure :: setup => tw_hodlr_setup
  !> Compute compressed inductance operator
  procedure :: compute_L => tw_hodlr_Lcompute
  !> Compute compressed B-field operators
  procedure :: compute_B => tw_hodlr_Bcompute
  !> Apply the matrix (L)
  procedure :: apply_real => tw_hodlr_Lapply
  !> Apply B-field operator
  procedure :: apply_bop => tw_hodlr_Bapply
  !> Assemble matrix (L)
  procedure :: assemble => tw_hodlr_Lassemble
end type oft_tw_hodlr_op
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
TYPE :: cmat_container
  COMPLEX(c8), POINTER, DIMENSION(:,:) :: M => NULL()
END TYPE
!---------------------------------------------------------------------------
!> Complex block-Jacobi preconditioner for ThinCurr HODLR matrices
!---------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_csolver) :: oft_tw_hodlr_bjpre
  INTEGER(4) :: max_block_size = 0
  COMPLEX(c8) :: alpha = (1.d0,0.d0)
  COMPLEX(c8) :: beta =  (1.d0,0.d0)
  TYPE(oft_tw_hodlr_op), POINTER :: mf_obj => NULL()
  TYPE(oft_native_matrix), POINTER :: Rmat => NULL()
  TYPE(cmat_container), POINTER, DIMENSION(:) :: inverse_mats => NULL()
CONTAINS
  !> Solve system
  PROCEDURE :: apply => bjprecond_apply
  !> Clean-up internal storage
  PROCEDURE :: delete => bjprecond_delete
END TYPE oft_tw_hodlr_bjpre
!---------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------
TYPE :: rmat_container
  REAL(r8), POINTER, DIMENSION(:,:) :: M => NULL()
END TYPE
!---------------------------------------------------------------------------
!> Real block-Jacobi preconditioner for ThinCurr HODLR matrices
!---------------------------------------------------------------------------
TYPE, PUBLIC, EXTENDS(oft_solver) :: oft_tw_hodlr_rbjpre
  INTEGER(4) :: max_block_size = 0
  REAL(r8) :: alpha = 1.d0
  REAL(r8) :: beta = 1.d0
  TYPE(oft_tw_hodlr_op), POINTER :: mf_obj => NULL()
  TYPE(oft_native_matrix), POINTER :: Rmat => NULL()
  TYPE(rmat_container), POINTER, DIMENSION(:) :: inverse_mats => NULL()
CONTAINS
  !> Solve system
  PROCEDURE :: apply => rbjprecond_apply
  !> Clean-up internal storage
  PROCEDURE :: delete => rbjprecond_delete
END TYPE oft_tw_hodlr_rbjpre
CONTAINS
!------------------------------------------------------------------------------
!> Compute mutual inductance matrix between two thin-wall models
!------------------------------------------------------------------------------
SUBROUTINE tw_compute_LmatHole(row_obj,col_obj,Lmat)
TYPE(tw_type), INTENT(in) :: row_obj !< Thin-wall model object for rows
TYPE(tw_type), INTENT(in) :: col_obj !< Thin-wall model object for columns
REAL(8), CONTIGUOUS, POINTER, INTENT(inout) :: Lmat(:,:) !< Mutual inductance matrix
INTEGER(4) :: i,ii,j,jj,ik,jk,imin,jmax,io_unit,ierr,file_counts(2),iquad
INTEGER(8) :: counts(4)
REAL(8) :: tmp,dl_min,dl_max,f(3),elapsed_time
REAL(8) :: pt_i(3),rgop(3,3),area_i,norm_i(3),evec_i(3,3),pts_i(3,3)
REAL(8) :: pt_j(3),cgop(3,3),area_j,norm_j(3),evec_j(3,3),pts_j(3,3)
LOGICAL :: vvclose_flag, close_flag,exists
CLASS(oft_bmesh), POINTER :: rmesh,cmesh
TYPE(oft_quad_type), ALLOCATABLE :: quads(:)
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
! ALLOCATE(Lmat(col_obj%nelems,row_obj%nelems))
!
! WRITE(*,*)'Building element<->element inductance matrix'
Lmat=0.d0
! CALL mytimer%tick
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
    IF(row_obj%kfh(i+1)-row_obj%kfh(i)==0)CYCLE
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
    DO ii=row_obj%kfh(i),row_obj%kfh(i+1)-1
        ik=ABS(row_obj%lfh(1,ii))
        DO jj=1,3
        jk=col_obj%pmap(cmesh%lc(jj,j))
        IF(jk<=0)CYCLE
        dl_min = SIGN(1,row_obj%lfh(1,ii))*DOT_PRODUCT(evec_i(:,row_obj%lfh(2,ii)),evec_j(:,jj))*tmp
        !$omp atomic
        Lmat(jk,ik) = Lmat(jk,ik) + dl_min
        END DO
        DO jj=col_obj%kfh(j),col_obj%kfh(j+1)-1
        jk=ABS(col_obj%lfh(1,jj))+col_obj%np_active
        dl_min = SIGN(1,row_obj%lfh(1,ii))*SIGN(1,col_obj%lfh(1,jj))* &
            DOT_PRODUCT(evec_i(:,row_obj%lfh(2,ii)),evec_j(:,col_obj%lfh(2,jj)))*tmp
        !$omp atomic
        Lmat(jk,ik) = Lmat(jk,ik) + dl_min
        END DO
    END DO
    END DO
END DO
!$omp end do nowait
!---Add passive coils to model
!$omp do
DO i=1,row_obj%np_active+row_obj%nholes
    DO j=1,col_obj%n_vcoils
    jj=col_obj%nholes+j
    Lmat(i,jj) = row_obj%Ael2coil(i,j)
    END DO
END DO
!$omp end do nowait
!$omp do
DO i=1,row_obj%nholes
    DO j=1,col_obj%n_vcoils
    jj=row_obj%np_active+col_obj%nholes+j
    Lmat(jj,i) = row_obj%Ael2coil(row_obj%np_active+i,j)
    END DO
END DO
!$omp end do nowait
!$omp do
DO i=1,row_obj%n_vcoils
    ii=row_obj%nholes+i
    DO j=1,col_obj%n_vcoils
    jj=col_obj%np_active+col_obj%nholes+j
    Lmat(jj,ii) = row_obj%Acoil2coil(i,j)
    END DO
END DO
!$omp end parallel
Lmat = Lmat/(4.d0*pi)
DO i=1,18
    CALL quads(i)%delete()
END DO
DEALLOCATE(quads)
! elapsed_time=mytimer%tock()
! WRITE(*,*)'  Time = ',elapsed_time
DEBUG_STACK_POP
END SUBROUTINE tw_compute_LmatHole
!------------------------------------------------------------------------------
!> Compute mutual inductance matrix between two thin-wall models
!------------------------------------------------------------------------------
SUBROUTINE tw_compute_Lmatblock(row_obj,col_obj,Lmat,row_block,col_block)
TYPE(tw_type), INTENT(in) :: row_obj !< Thin-wall model object for rows
TYPE(tw_type), INTENT(in) :: col_obj !< Thin-wall model object for columns
TYPE(oft_tw_block), INTENT(in) :: row_block !< Block of rows to compute
TYPE(oft_tw_block), INTENT(in) :: col_block !< Block of columns to compute
REAL(8), INTENT(out) :: Lmat(:,:) !< Mutual inductance matrix
INTEGER(4) :: i,ii,j,jj,ik,jk,imin,jmax,io_unit,ierr,file_counts(2),iquad,irow,jcol
INTEGER(4) :: rowi,coli
INTEGER(8) :: counts(4)
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: rpmap,cpmap
REAL(8) :: tmp,dl_min,dl_max,f(3),elapsed_time
REAL(8) :: pt_i(3),rgop(3,3),area_i,norm_i(3),evec_i(3,3),pts_i(3,3)
REAL(8) :: pt_j(3),cgop(3,3),area_j,norm_j(3),evec_j(3,3),pts_j(3,3)
LOGICAL :: vvclose_flag, close_flag,exists
CLASS(oft_bmesh), POINTER :: rmesh,cmesh
TYPE(oft_quad_type), ALLOCATABLE :: quads(:)
type(oft_timer) :: mytimer
DEBUG_STACK_PUSH
!
! WRITE(*,*)'Building element<->element inductance matrix block'
Lmat=0.d0
! CALL mytimer%tick
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
DO irow=1,row_block%ncells
    i = row_block%icell(irow)
    ! CALL rmesh%jacobian(i,f,rgop,area_i)
    ! CALL rmesh%norm(i,f,norm_i)
    area_i=rmesh%ca(i)
    DO ii=1,3
    pts_i(:,ii)=rmesh%r(:,rmesh%lc(ii,i))
    ! evec_i(:,ii)=cross_product(rgop(:,ii),norm_i)
    evec_i(:,ii)=row_obj%qbasis(:,ii,i)
    END DO
    !---Compute inter-edge inductances
    DO jcol=1,col_block%ncells
    j = col_block%icell(jcol)
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
        ! !$omp simd collapse(1) private(pt_i,pt_j) reduction(+:tmp)
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
        ik=row_block%inv_map(rmesh%lc(ii,i))
        IF(ik==0)CYCLE
        DO jj=1,3
        jk=col_block%inv_map(cmesh%lc(jj,j))
        IF(jk==0)CYCLE
        dl_min = DOT_PRODUCT(evec_i(:,ii),evec_j(:,jj))*tmp
        ! !$omp atomic
        Lmat(jk,ik) = Lmat(jk,ik) + dl_min
        END DO
    END DO
    END DO
END DO
Lmat = Lmat/(4.d0*pi)
DO i=1,18
    CALL quads(i)%delete()
END DO
DEALLOCATE(quads) !,rpmap,cpmap)
! elapsed_time=mytimer%tock()
! WRITE(*,*)'  Time = ',elapsed_time
DEBUG_STACK_POP
END SUBROUTINE tw_compute_Lmatblock
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE tw_compute_Bops_hole(self,Bop,Bop_dr)
TYPE(tw_type), INTENT(inout) :: self
REAL(8), INTENT(out) :: Bop(:,:,:) !< Magnetic field evaluation matrix
REAL(8), INTENT(out) :: Bop_dr(:,:,:) !< Magnetic field evaluation matrix
REAL(r8) :: evec_i(3,3),evec_j(3),pts_i(3,3),pt_i(3),pt_j(3),diffvec(3),ecc(3)
REAL(r8) :: r1,z1,rmag,cvec(3),cpt(3),tmp,area_i,dl_min,dl_max,norm_j(3),f(3),pot_tmp,pot_last
REAL(r8), ALLOCATABLE :: atmp(:,:,:)
REAL(8), PARAMETER :: B_dx = 1.d-6
INTEGER(4) :: i,ii,j,jj,ik,jk,k,kk,iquad
LOGICAL :: is_neighbor
CLASS(oft_bmesh), POINTER :: bmesh
TYPE(oft_quad_type), ALLOCATABLE :: quads(:)
bmesh=>self%mesh
ALLOCATE(quads(18))
DO i=1,18
    CALL set_quad_2d(quads(i),i)
END DO
Bop=0.d0
f=1.d0/3.d0
! WRITE(*,*)'Building element->element magnetic reconstruction operator'
!$omp parallel private(ii,j,jj,ik,pts_i,tmp,pt_i,pt_j,evec_i, &
!$omp atmp,i,area_i,dl_min,dl_max,norm_j,diffvec,is_neighbor,iquad)
ALLOCATE(atmp(3,3,bmesh%np))
!$omp do schedule(dynamic,100)
DO i=1,bmesh%nc
  IF(self%kfh(i+1)-self%kfh(i)==0)CYCLE
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
  DO ii=self%kfh(i),self%kfh(i+1)-1
    ik=ABS(self%lfh(1,ii))!+self%np_active
    DO j=1,bmesh%np
      DO jj=1,3
        tmp=SIGN(1,self%lfh(1,ii))*atmp(jj,self%lfh(2,ii),j)
        !$omp atomic
        Bop(j,ik,jj) = Bop(j,ik,jj) + tmp
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
! WRITE(*,*)'Building vcoil->element magnetic reconstruction operator'
IF(self%n_vcoils>0)THEN
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
        Bop(i,self%nholes+j,jj) = &
          Bop(i,self%nholes+j,jj) + ecc(jj)
      END DO
    END DO
  END DO
END IF
Bop=Bop/(4.d0*pi)
!
Bop_dr=0.d0
! WRITE(*,*)'Building icoil->element magnetic reconstruction operator'
IF(self%n_icoils>0)THEN
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
        Bop_dr(i,j,jj) = Bop_dr(i,j,jj) + ecc(jj)
      END DO
    END DO
  END DO
END IF
Bop_dr=Bop_dr*mu0/(4.d0*pi)
END SUBROUTINE tw_compute_Bops_hole
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE tw_compute_Bops_block(self,Bop,row_block,col_block,dir)
TYPE(tw_type), INTENT(inout) :: self
REAL(8), INTENT(out) :: Bop(:,:) !< Magnetic field evaluation matrix
TYPE(oft_tw_block), INTENT(in) :: row_block !< Block of rows to compute
TYPE(oft_tw_block), INTENT(in) :: col_block !< Block of columns to compute
INTEGER(4), INTENT(in) :: dir
REAL(r8) :: evec_i(3,3),evec_j(3),pts_i(3,3),pt_i(3),pt_j(3),diffvec(3),ecc(3)
REAL(r8) :: r1,z1,rmag,cvec(3),cpt(3),tmp,area_i,dl_min,dl_max,norm_j(3),f(3),pot_tmp,pot_last
REAL(r8), ALLOCATABLE :: atmp(:,:,:)
REAL(8), PARAMETER :: B_dx = 1.d-6
INTEGER(4) :: i,ii,j,jj,ik,jk,k,kk,iquad,irow,jcol
LOGICAL :: is_neighbor
CLASS(oft_bmesh), POINTER :: bmesh
TYPE(oft_quad_type), ALLOCATABLE :: quads(:)
bmesh=>self%mesh
ALLOCATE(quads(18))
DO i=1,18
    CALL set_quad_2d(quads(i),i)
END DO
Bop=0.d0
f=1.d0/3.d0
!WRITE(*,*)'Building element->element magnetic reconstruction operator'
ALLOCATE(atmp(3,3,col_block%np))
DO irow=1,row_block%ncells
    i = row_block%icell(irow)
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
    DO jcol=1,col_block%np
    j=col_block%ipts(jcol)
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
            atmp(1,ik,jcol) = diffvec(2)*evec_i(3,ik) - diffvec(3)*evec_i(2,ik)
            atmp(2,ik,jcol) = diffvec(3)*evec_i(1,ik) - diffvec(1)*evec_i(3,ik)
            atmp(3,ik,jcol) = diffvec(1)*evec_i(2,ik) - diffvec(2)*evec_i(1,ik)
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
            atmp(:,ik,jcol) = diffvec*area_i
        END DO
    END IF
    END DO
    DO ii=1,3
        ik=row_block%inv_map(bmesh%lc(ii,i)) !self%pmap(bmesh%lc(ii,i))
        IF(ik==0)CYCLE
        DO j=1,col_block%np
            Bop(j,ik) = Bop(j,ik) + atmp(dir,ii,j)
        END DO
    END DO
END DO
Bop=Bop/(4.d0*pi)
DEALLOCATE(atmp)
DO i=1,18
    CALL quads(i)%delete()
END DO
DEALLOCATE(quads)
END SUBROUTINE tw_compute_Bops_block
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE tw_hodlr_setup(self,required)
class(oft_tw_hodlr_op), intent(inout) :: self
LOGICAL, INTENT(in) :: required
INTEGER(4) :: i,j,k,l,nlast,branch_count,pind_i,pind_j,max_pts,top_level,level,io_unit,ierr
INTEGER(4) :: size_out,nblocks_progress(9),nblocks_complete,counts(4),iblock,jblock,sparse_count(2)
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: kc_tmp,kr_tmp
INTEGER(4), ALLOCATABLE, DIMENSION(:,:) :: block_interaction,row_count
INTEGER(8) :: mat_size(3),compressed_size
LOGICAL :: is_masked,is_neighbor
LOGICAL, ALLOCATABLE, DIMENSION(:) :: near_closure,cell_mark
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: dense_blocks
REAL(8) :: corners_i(3,8),corners_j(3,8),rmin(3),elapsed_time,avg_size
REAL(8), ALLOCATABLE, DIMENSION(:) :: point_block
TYPE(oft_native_dense_matrix), POINTER :: mat_tmp => NULL()
type(oft_timer) :: mytimer
!---Tuning options
INTEGER(4) :: target_size = 1500
INTEGER(4) :: aca_min_its = 20
INTEGER(4) :: min_rank = 10
REAL(8) :: L_svd_tol = -1.d0
REAL(8) :: L_aca_rel_tol = -1.d0
REAL(8) :: B_svd_tol = -1.d0
REAL(8) :: B_aca_rel_tol = -1.d0
NAMELIST/thincurr_hodlr_options/target_size,min_rank,aca_min_its,L_svd_tol,L_aca_rel_tol, &
  B_svd_tol,B_aca_rel_tol
!
IF(self%tw_obj%kpmap_inv(self%tw_obj%np_active+1)/=self%tw_obj%np_active+1)THEN
  WRITE(*,*)self%tw_obj%kpmap_inv(self%tw_obj%np_active+1),self%tw_obj%np_active+1
  CALL oft_abort( &
    "HODLR not supported with periodic grids","tw_hodlr_setup",__FILE__)
END IF
!---
OPEN(NEWUNIT=io_unit, FILE=oft_env%ifile)
READ(io_unit,thincurr_hodlr_options, IOSTAT=ierr)
CLOSE(io_unit)
if(ierr<0)THEN
  IF(required)THEN
    CALL oft_abort('No "thincurr_hodlr_options" found in input file.', &
      'tw_hodlr_setup',__FILE__)
  ELSE
    WRITE(*,*)'No "thincurr_hodlr_options" found in input file, using full representation'
    RETURN
  END IF
END IF
if(ierr>0)call oft_abort('Error parsing "thincurr_hodlr_options" in input file.', &
  'tw_hodlr_setup',__FILE__)
IF((L_svd_tol>0.d0).AND.(min_rank>aca_min_its))CALL oft_abort('"min_rank" must be <= "aca_min_its"', &
  "tw_hodlr_setup",__FILE__)
self%leaf_target=target_size
self%min_rank=min_rank
self%aca_min_its=aca_min_its
self%L_svd_tol=L_svd_tol
IF(L_aca_rel_tol>0.d0)self%L_aca_tol=L_aca_rel_tol*self%L_svd_tol
self%B_svd_tol=B_svd_tol
IF(B_aca_rel_tol>0.d0)self%B_aca_tol=B_aca_rel_tol*self%B_svd_tol
IF((self%L_svd_tol<0.d0).AND.(self%B_svd_tol<0.d0))RETURN
!---Partition mesh
WRITE(*,*)'Partitioning grid for block low rank compressed operators'
CALL tw_part_mesh(self%tw_obj,self%leaf_target,self%nlevels,self%levels)
sparse_count=0
DO i=1,self%nlevels
  ! WRITE(*,*)i,self%levels(i)%nblocks
  ALLOCATE(self%levels(i)%mat_mask(self%levels(i)%nblocks,self%levels(i)%nblocks))
  self%levels(i)%mat_mask = 0
  IF(i>1)THEN
    DO j=1,self%levels(i)%nblocks
      DO k=1,self%levels(i)%nblocks
        IF(self%levels(i-1)%mat_mask(self%levels(i)%blocks(j)%parent,self%levels(i)%blocks(k)%parent)/=0)THEN
          ! WRITE(*,*)'Masking',i,j,k
          self%levels(i)%mat_mask(j,k)=-1
          self%levels(i)%mat_mask(k,j)=-1
        END IF
      END DO
    END DO
  END IF
  DO j=1,self%levels(i)%nblocks
    DO k=j+1,self%levels(i)%nblocks
      ! WRITE(*,*)i,j,k,magnitude(self%levels(i)%blocks(j)%center-self%levels(i)%blocks(k)%center), &
      ! (self%levels(i)%blocks(j)%extent+self%levels(i)%blocks(k)%extent)
      IF(self%levels(i)%mat_mask(j,k)/=0)CYCLE
      IF(magnitude(self%levels(i)%blocks(j)%center-self%levels(i)%blocks(k)%center)/ &
        (self%levels(i)%blocks(j)%extent+self%levels(i)%blocks(k)%extent) > 1.1d0)THEN
        self%levels(i)%mat_mask(j,k)=2
        self%levels(i)%mat_mask(k,j)=-2
      ELSE IF(i==self%nlevels)THEN
        IF(magnitude(self%levels(i)%blocks(j)%center-self%levels(i)%blocks(k)%center)/ &
          (self%levels(i)%blocks(j)%extent+self%levels(i)%blocks(k)%extent) > 1.1d0)THEN
          self%levels(i)%mat_mask(j,k)=2
          self%levels(i)%mat_mask(k,j)=-2
        ELSE
          self%levels(i)%mat_mask(j,k)=3
          self%levels(i)%mat_mask(k,j)=-3
        END IF
      END IF
    END DO
    IF(i==self%nlevels)self%levels(i)%mat_mask(j,j)=1
  END DO
  counts(1) = SUM(COUNT(self%levels(i)%mat_mask==0,DIM=1))
  counts(2) = SUM(COUNT(self%levels(i)%mat_mask==1,DIM=1))
  counts(3) = SUM(COUNT(self%levels(i)%mat_mask==2,DIM=1))
  counts(4) = SUM(COUNT(self%levels(i)%mat_mask==3,DIM=1))
  ! WRITE(*,*)'  ',counts
  self%ndense=self%ndense+counts(2)
  self%nsparse=self%nsparse+counts(3)+counts(4)
  sparse_count=sparse_count+counts(3:4)
END DO
! WRITE(*,*)self%ndense,self%nsparse
ALLOCATE(self%dense_blocks(3,self%ndense),self%sparse_blocks(3,self%nsparse))
self%ndense=0
self%nsparse=0
DO i=1,self%nlevels
  DO j=1,self%levels(i)%nblocks
    DO k=1,self%levels(i)%nblocks
      IF(self%levels(i)%mat_mask(j,k)==1)THEN
        self%ndense=self%ndense+1
        self%dense_blocks(:,self%ndense)=[i,j,k]
      ELSE IF(self%levels(i)%mat_mask(j,k)>=2)THEN
        self%nsparse=self%nsparse+1
        self%sparse_blocks(:,self%nsparse)=[i,j,k]
      END IF
    END DO
  END DO
END DO
! WRITE(*,*)self%ndense,self%nsparse
! CALL oft_abort("","",__FILE__)
self%nblocks=self%levels(self%nlevels)%nblocks
self%blocks=>self%levels(self%nlevels)%blocks
WRITE(*,*)'  nBlocks =        ',self%nblocks
ALLOCATE(point_block(self%tw_obj%mesh%np))
point_block=0.d0
avg_size=0.d0
DO i=1,self%nblocks
  avg_size=avg_size+self%blocks(i)%nelems
  ! DO j=1,self%blocks(i)%nelems
  !     point_block(self%tw_obj%lpmap_inv(self%blocks(i)%ielem))=i*1.d0
  ! END DO
  DO j=1,self%blocks(i)%np
    point_block(self%blocks(i)%ipts)=i*1.d0
  END DO
END DO
WRITE(*,*)'  Avg block size = ',INT(avg_size/REAL(self%nblocks,8),4)
IF(self%L_aca_tol<=0.d0)THEN
  sparse_count=[0,sparse_count(1)+sparse_count(2)]
END IF
WRITE(*,*)'  # of SVD =       ',sparse_count(2)
WRITE(*,*)'  # of ACA =       ',sparse_count(1)
CALL self%tw_obj%mesh%save_vertex_scalar(point_block, 'ACA_Block')
ALLOCATE(cell_mark(self%tw_obj%mesh%nc))
DO j=1,self%nlevels
  DO i=1,self%levels(j)%nblocks
    IF(oft_debug_print(1))WRITE(*,*)i,self%levels(j)%blocks(i)%nelems
    !
    IF(ASSOCIATED(self%levels(j)%blocks(i)%inv_map))CYCLE
    ALLOCATE(self%levels(j)%blocks(i)%inv_map(self%tw_obj%mesh%np))
    self%levels(j)%blocks(i)%inv_map = 0
    DO l=1,self%levels(j)%blocks(i)%nelems
      k=self%tw_obj%lpmap_inv(self%levels(j)%blocks(i)%ielem(l))
      self%levels(j)%blocks(i)%inv_map(k)=l
    END DO
    !
    cell_mark=.FALSE.
    DO l=1,self%levels(j)%blocks(i)%nelems
      pind_j = self%tw_obj%lpmap_inv(self%levels(j)%blocks(i)%ielem(l))
      DO k=self%tw_obj%mesh%kpc(pind_j),self%tw_obj%mesh%kpc(pind_j+1)-1
        cell_mark(self%tw_obj%mesh%lpc(k))=.TRUE.
      END DO
    END DO
    self%levels(j)%blocks(i)%ncells=COUNT(cell_mark)
    ALLOCATE(self%levels(j)%blocks(i)%icell(self%levels(j)%blocks(i)%ncells))
    self%levels(j)%blocks(i)%ncells=0
    DO l=1,self%tw_obj%mesh%nc
      IF(cell_mark(l))THEN
        self%levels(j)%blocks(i)%ncells=self%levels(j)%blocks(i)%ncells+1
        self%levels(j)%blocks(i)%icell(self%levels(j)%blocks(i)%ncells)=l
      END IF
    END DO
  END DO
END DO
DEALLOCATE(cell_mark)
END SUBROUTINE tw_hodlr_setup
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE tw_hodlr_Lcompute(self,save_file)
class(oft_tw_hodlr_op), intent(inout) :: self
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: save_file
INTEGER(4) :: i,j,k,l,nlast,branch_count,pind_i,pind_j,max_pts,top_level,level
INTEGER(4) :: nblocks_progress(9),nblocks_complete,counts(4),iblock,jblock
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: kc_tmp,kr_tmp
INTEGER(4), ALLOCATABLE, DIMENSION(:,:) :: block_interaction,row_count
INTEGER(8) :: mat_size(3),compressed_size,size_out
LOGICAL :: is_masked,is_neighbor,mat_updated
LOGICAL, ALLOCATABLE, DIMENSION(:) :: near_closure,cell_mark
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: dense_blocks
REAL(8) :: corners_i(3,8),corners_j(3,8),rmin(3),elapsed_time,avg_size
REAL(8), ALLOCATABLE, DIMENSION(:) :: point_block
TYPE(oft_native_dense_matrix), POINTER :: mat_tmp => NULL()
type(oft_timer) :: mytimer
WRITE(*,*)'Building block low rank inductance operator'
IF(self%L_svd_tol<0.d0)CALL oft_abort('"L_svd_tol" must be > 0','tw_Lmat_MF_Lcompute',__FILE__)
!---Build matrix with ACA+ or SVD compression of off-diagonal blocks
CALL mytimer%tick()
ALLOCATE(self%dense_mats(self%ndense))
ALLOCATE(self%aca_U_mats(self%nsparse))
ALLOCATE(self%aca_V_mats(self%nsparse))
ALLOCATE(self%aca_dense(self%nsparse))
WRITE(*,*)'  Building hole and Vcoil columns'
IF(self%tw_obj%nholes+self%tw_obj%n_vcoils>0)THEN
  ALLOCATE(self%hole_Vcoil_mat%M(self%tw_obj%nelems,self%tw_obj%nholes+self%tw_obj%n_vcoils))
  CALL tw_compute_LmatHole(self%tw_obj,self%tw_obj,self%hole_Vcoil_mat%M)
END IF
IF(PRESENT(save_file))CALL load_from_file()
compressed_size=0
mat_updated=.FALSE.
WRITE(*,*)'  Building diagonal blocks'
!$omp parallel private(i,level,j,k,size_out,avg_size) reduction(+:compressed_size) &
!$omp reduction(.OR.:mat_updated)
!$omp single
nblocks_progress = INT(self%ndense*[0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0],4)
nblocks_complete = 0
!$omp end single
!$omp do schedule(static,1)
DO i=1,self%ndense
  level = self%dense_blocks(1,i)
  j = self%dense_blocks(2,i)
  k = self%dense_blocks(3,i)
  !---Build full representation
  IF(.NOT.ASSOCIATED(self%dense_mats(i)%M))THEN
    mat_updated=.TRUE.
    ALLOCATE(self%dense_mats(i)%M(self%levels(level)%blocks(j)%nelems,self%levels(level)%blocks(k)%nelems))
    CALL tw_compute_Lmatblock(self%tw_obj,self%tw_obj,self%dense_mats(i)%M, &
      self%levels(level)%blocks(k),self%levels(level)%blocks(j))
  END IF
  compressed_size = compressed_size + self%levels(level)%blocks(j)%nelems*self%levels(level)%blocks(k)%nelems
  !$omp critical
  nblocks_complete = nblocks_complete + 1
  DO k=1,9
  IF(nblocks_complete==nblocks_progress(k))WRITE(*,'(5X,I2,A1)')k*10,'%'
  END DO
  !$omp end critical
END DO
!$omp single
IF(self%L_aca_tol>0.d0)THEN
  WRITE(*,*)'  Building off-diagonal blocks using ACA+'
ELSE
  WRITE(*,*)'  Building off-diagonal blocks with SVD compression'
END IF
nblocks_progress = INT(self%nsparse*[0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0],4)
nblocks_complete = 0
!$omp end single
!$omp do schedule(dynamic,1)
DO i=1,self%nsparse
  level = self%sparse_blocks(1,i)
  j = self%sparse_blocks(2,i)
  k = self%sparse_blocks(3,i)
  !---Build full representation
  IF(ASSOCIATED(self%aca_dense(i)%M))THEN
    size_out=self%levels(level)%blocks(k)%nelems*self%levels(level)%blocks(j)%nelems
  ELSE IF(ASSOCIATED(self%aca_U_mats(i)%M))THEN
    size_out=SIZE(self%aca_V_mats(i)%M,DIM=1)*(self%levels(level)%blocks(k)%nelems+self%levels(level)%blocks(j)%nelems)
  ELSE
    mat_updated=.TRUE.
    size_out=-1
    IF((self%L_aca_tol>0.d0).AND.self%levels(level)%mat_mask(j,k)==2)THEN
      CALL aca_approx(i,self%L_aca_tol,size_out)
      IF(size_out>0)CALL compress_aca(i,self%L_svd_tol,size_out)
      IF(size_out==-1)WRITE(*,*)'ACA+ failure',self%levels(level)%blocks(j)%nelems, &
        self%levels(level)%blocks(k)%nelems
    END IF
    !---IF ACA+ failed or diagonal compute dense matrix
    IF(size_out<0)THEN
      ALLOCATE(self%aca_dense(i)%M(self%levels(level)%blocks(k)%nelems,self%levels(level)%blocks(j)%nelems))
      CALL tw_compute_Lmatblock(self%tw_obj,self%tw_obj,self%aca_dense(i)%M, &
        self%levels(level)%blocks(j),self%levels(level)%blocks(k))
      CALL compress_block(i,self%L_svd_tol,size_out)
    END IF
  END IF
  compressed_size = compressed_size + size_out
  !$omp critical
  nblocks_complete = nblocks_complete + 1
  DO k=1,9
    IF(nblocks_complete==nblocks_progress(k))WRITE(*,'(5X,I2,A1)')k*10,'%'
  END DO
  !$omp end critical
END DO
!$omp end parallel
elapsed_time=mytimer%tock()
WRITE(*,'(5X,A,F6.1,A,ES9.2,A,ES9.2,A)')'Compression ratio:',compressed_size*1.d2/(REAL(self%tw_obj%np_active,8)**2), &
    "%  (",REAL(compressed_size,8),"/",REAL(self%tw_obj%np_active,8)**2,")"
WRITE(*,'(5X,2A)')'Time = ',time_to_string(elapsed_time)
!
IF(mat_updated.AND.PRESENT(save_file))CALL save_to_file()
CONTAINS
SUBROUTINE save_to_file()
INTEGER(4) :: ierr,io_unit,hash_tmp(6),file_counts(6)
CHARACTER(LEN=5) :: matrix_id
!---Save computed matrix to file
IF(TRIM(save_file)/='none')THEN
  hash_tmp(1) = self%tw_obj%nelems
  hash_tmp(2) = self%tw_obj%mesh%nc
  hash_tmp(3) = self%ndense
  hash_tmp(4) = self%nsparse
  hash_tmp(5) = oft_simple_hash(C_LOC(self%tw_obj%mesh%lc),INT(4*3*self%tw_obj%mesh%nc,8))
  hash_tmp(6) = oft_simple_hash(C_LOC(self%tw_obj%mesh%r),INT(8*3*self%tw_obj%mesh%np,8))
  WRITE(*,*)'  Saving HODLR matrix to file: ',TRIM(save_file)
  CALL hdf5_create_file(TRIM(save_file))
  CALL hdf5_write(hash_tmp,TRIM(save_file),'MODEL_hash')
  !
  hash_tmp(6)=0
  DO i=1,self%ndense
    level = self%dense_blocks(1,i)
    j = self%dense_blocks(2,i)
    k = self%dense_blocks(3,i)
    hash_tmp(1:3) = [level,j,k]
    hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ielem),INT(4*self%levels(level)%blocks(j)%nelems,8))
    hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ielem),INT(4*self%levels(level)%blocks(k)%nelems,8))
    WRITE(matrix_id,'(I5)')i
    CALL hdf5_write(hash_tmp,TRIM(save_file),'DENSE_hash_'//TRIM(ADJUSTL(matrix_id)))
    CALL hdf5_write(self%dense_mats(i)%M,TRIM(save_file),'DENSE_'//TRIM(ADJUSTL(matrix_id)))
  END DO
  !
  DO i=1,self%nsparse
    level = self%sparse_blocks(1,i)
    j = self%sparse_blocks(2,i)
    k = self%sparse_blocks(3,i)
    hash_tmp(1:3) = [level,j,k]
    hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ielem),INT(4*self%levels(level)%blocks(j)%nelems,8))
    hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ielem),INT(4*self%levels(level)%blocks(k)%nelems,8))
    WRITE(matrix_id,'(I5)')i
    IF(ASSOCIATED(self%aca_dense(i)%M))THEN
      hash_tmp(6)=-1
      CALL hdf5_write(hash_tmp,TRIM(save_file),'SPARSE_hash_'//TRIM(ADJUSTL(matrix_id)))
      CALL hdf5_write(self%aca_dense(i)%M,TRIM(save_file),'SPARSE_'//TRIM(ADJUSTL(matrix_id)))
    ELSE
      hash_tmp(6)=SIZE(self%aca_V_mats(i)%M,DIM=1)
      CALL hdf5_write(hash_tmp,TRIM(save_file),'SPARSE_hash_'//TRIM(ADJUSTL(matrix_id)))
      CALL hdf5_write(self%aca_U_mats(i)%M,TRIM(save_file),'SPARSE_U_'//TRIM(ADJUSTL(matrix_id)))
      CALL hdf5_write(self%aca_V_mats(i)%M,TRIM(save_file),'SPARSE_V_'//TRIM(ADJUSTL(matrix_id)))
    END IF
  END DO
  ! TODO: Save create appropriate hash to enable saving Vcoil matrix
  ! WRITE(io_unit)hash_tmp
  ! WRITE(io_unit)self%hole_Vcoil_mat%M
  ! CLOSE(io_unit)
END IF
END SUBROUTINE save_to_file
!
SUBROUTINE load_from_file()
LOGICAL :: exists
INTEGER(4) :: ierr,io_unit,hash_tmp(6),file_counts(6)
CHARACTER(LEN=5) :: matrix_id
!---Try to load from file
IF(TRIM(save_file)/='none')THEN
  INQUIRE(FILE=TRIM(save_file),EXIST=exists)
  IF(exists)THEN
    hash_tmp(1) = self%tw_obj%nelems
    hash_tmp(2) = self%tw_obj%mesh%nc
    hash_tmp(3) = self%ndense
    hash_tmp(4) = self%nsparse
    hash_tmp(5) = oft_simple_hash(C_LOC(self%tw_obj%mesh%lc),INT(4*3*self%tw_obj%mesh%nc,8))
    hash_tmp(6) = oft_simple_hash(C_LOC(self%tw_obj%mesh%r),INT(8*3*self%tw_obj%mesh%np,8))
    WRITE(*,*)'  Reading HODLR matrix from file: ',TRIM(save_file)
    CALL hdf5_read(file_counts,TRIM(save_file),'MODEL_hash',success=exists)
    IF(exists.AND.ALL(file_counts==hash_tmp))THEN
      hash_tmp(6)=0
      DO i=1,self%ndense
        level = self%dense_blocks(1,i)
        j = self%dense_blocks(2,i)
        k = self%dense_blocks(3,i)
        hash_tmp(1:3) = [level,j,k]
        hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ielem),INT(4*self%levels(level)%blocks(j)%nelems,8))
        hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ielem),INT(4*self%levels(level)%blocks(k)%nelems,8))
        WRITE(matrix_id,'(I5)')i
        CALL hdf5_read(file_counts,TRIM(save_file),'DENSE_hash_'//TRIM(ADJUSTL(matrix_id)),success=exists)
        IF(exists.AND.ALL(file_counts==hash_tmp))THEN
          ALLOCATE(self%dense_mats(i)%M(self%levels(level)%blocks(j)%nelems,self%levels(level)%blocks(k)%nelems))
          CALL hdf5_read(self%dense_mats(i)%M,TRIM(save_file),'DENSE_'//TRIM(ADJUSTL(matrix_id)),success=exists)
          IF(.NOT.exists)DEALLOCATE(self%dense_mats(i)%M)
        END IF
      END DO
      !
      DO i=1,self%nsparse
        level = self%sparse_blocks(1,i)
        j = self%sparse_blocks(2,i)
        k = self%sparse_blocks(3,i)
        hash_tmp(1:3) = [level,j,k]
        hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ielem),INT(4*self%levels(level)%blocks(j)%nelems,8))
        hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ielem),INT(4*self%levels(level)%blocks(k)%nelems,8))
        WRITE(matrix_id,'(I5)')i
        CALL hdf5_read(file_counts,TRIM(save_file),'SPARSE_hash_'//TRIM(ADJUSTL(matrix_id)),success=exists)
        IF(exists.AND.ALL(file_counts(1:5)==hash_tmp(1:5)))THEN
          IF(file_counts(6)<0)THEN
            ALLOCATE(self%aca_dense(i)%M(self%levels(level)%blocks(k)%nelems,self%levels(level)%blocks(j)%nelems))
            CALL hdf5_read(self%aca_dense(i)%M,TRIM(save_file),'SPARSE_'//TRIM(ADJUSTL(matrix_id)),success=exists)
            IF(.NOT.exists)DEALLOCATE(self%aca_dense(i)%M)
          ELSE
            ALLOCATE(self%aca_U_mats(i)%M(self%levels(level)%blocks(k)%nelems,file_counts(6)))
            CALL hdf5_read(self%aca_U_mats(i)%M,TRIM(save_file),'SPARSE_U_'//TRIM(ADJUSTL(matrix_id)),success=exists)
            IF(.NOT.exists)THEN
              DEALLOCATE(self%aca_U_mats(i)%M)
            ELSE
              ALLOCATE(self%aca_V_mats(i)%M(file_counts(6),self%levels(level)%blocks(j)%nelems))
              CALL hdf5_read(self%aca_V_mats(i)%M,TRIM(save_file),'SPARSE_V_'//TRIM(ADJUSTL(matrix_id)),success=exists)
              IF(.NOT.exists)DEALLOCATE(self%aca_U_mats(i)%M,self%aca_V_mats(i)%M)
            END IF
          END IF
        END IF
      END DO
    ELSE
      WRITE(*,*)'    Ignoring stored matrix, Model hashes do not match'
    END IF
  END IF
END IF
END SUBROUTINE load_from_file
!
SUBROUTINE compress_block(isparse,tol,size_out)
INTEGER(4), INTENT(in) :: isparse
REAL(8), INTENT(in) :: tol
INTEGER(8), INTENT(out) :: size_out
INTEGER(4) :: i,j,io_unit
REAL(8) :: S_sum,frob_K
!
INTEGER(4) :: ilevel,iblock,jblock
INTEGER(4) :: M,N,LDA,LDU,LDVT,LWORK,INFO,MIN_DIM
INTEGER(8) :: full_size
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: IWORK
REAL(8), ALLOCATABLE, DIMENSION(:) :: S,WORK
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: U,VT,Atmp
!
ilevel=self%sparse_blocks(1,isparse)
iblock=self%sparse_blocks(2,isparse)
jblock=self%sparse_blocks(3,isparse)
!
M = SIZE(self%aca_dense(isparse)%M,DIM=1)
N = SIZE(self%aca_dense(isparse)%M,DIM=2)
full_size = INT(M,8)*INT(N,8)
LDA = M
LDU = M
MIN_DIM = min(M,N)
LDVT = MIN_DIM
ALLOCATE(Atmp(M,N))
Atmp = self%aca_dense(isparse)%M
ALLOCATE(S(MIN_DIM),U(LDU,MIN_DIM),VT(LDVT,N))
!---Get optimal workspace size
LWORK = -1
ALLOCATE(WORK(1),IWORK(8*MIN_DIM))
CALL DGESDD('S', M, N, Atmp, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
LWORK = INT(WORK(1))
DEALLOCATE(WORK)
!---Compute SVD of matrix block
ALLOCATE(WORK(LWORK))
CALL DGESDD('S', M, N, Atmp, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
!---Truncate at desired accuracy
! S_sum=SUM(S**2)
! frob_K=S_sum
! DO i=1,MIN_DIM
!     ! IF(SQRT(frob_K/S_sum)<tol)EXIT
!     IF(SQRT(frob_K)<tol)EXIT
!     frob_K = frob_K - S(i)**2
! END DO
frob_K=0.d0
DO i=MIN_DIM,1,-1
  frob_K = frob_K + S(i)**2
  IF(SQRT(frob_K)>tol)EXIT
END DO
i=MAX(MIN(self%min_rank,MIN_DIM),i)
IF(i<full_size/INT(M+N,8))THEN
  ALLOCATE(self%aca_U_mats(isparse)%M(M,i))
  self%aca_U_mats(isparse)%M=U(:,1:i)
  ALLOCATE(self%aca_V_mats(isparse)%M(i,N))
  DO j=1,i
  self%aca_V_mats(isparse)%M(j,:)=VT(j,:)*S(j)
  END DO
  DEALLOCATE(self%aca_dense(isparse)%M)
  ! WRITE(*,*)'Savings:',iblock,jblock,i,M,N
  size_out=INT(i,8)*INT(M+N,8)
ELSE
  ! WRITE(*,*)'No savings:',iblock,jblock
  size_out=full_size
END IF
!
IF( INFO.NE.0 ) THEN
  CALL oft_abort("The algorithm computing SVD failed to converge.","tw_Lmat_MF_setup::compress_block",__FILE__)
END IF
! CALL oft_abort("","",__FILE__)
DEALLOCATE(WORK,IWORK,S,U,VT,Atmp)
END SUBROUTINE compress_block
!
FUNCTION max_masked(vals,mask) RESULT(iloc)
REAL(8), INTENT(in) :: vals(:)
LOGICAL, INTENT(in) :: mask(:)
REAL(8) :: max_tmp
INTEGER(4) :: i,iloc
max_tmp = -1.d99
DO i=1,SIZE(vals)
  IF((vals(i)>max_tmp).AND.mask(i))THEN
  max_tmp=vals(i)
  iloc=i
  END IF
END DO
END FUNCTION max_masked
!
SUBROUTINE aca_approx(isparse,tol,size_out)
INTEGER(4), INTENT(in) :: isparse
REAL(8), INTENT(in) :: tol
INTEGER(8), INTENT(out) :: size_out
!
INTEGER(4) :: ilevel,iblock,jblock
INTEGER(4) :: M,N,MIN_DIM,max_iter,Iref,Jref,Istar,Jstar,i,k,nterms,k1,k2
INTEGER(8) :: full_size
REAL(8) :: tmp,step_size,Istar_val,Jstar_val
REAL(8), ALLOCATABLE :: us(:,:),vs(:,:),RIstar(:,:),RJstar(:,:),RIref(:,:),RJref(:,:)
LOGICAL, ALLOCATABLE :: prevIstar(:),prevJstar(:)
TYPE(oft_tw_block) :: tmp_block
!
ilevel=self%sparse_blocks(1,isparse)
iblock=self%sparse_blocks(2,isparse)
jblock=self%sparse_blocks(3,isparse)
size_out=-2
M = self%levels(ilevel)%blocks(jblock)%nelems
N = self%levels(ilevel)%blocks(iblock)%nelems
MIN_DIM = min(M,N)
IF(MIN_DIM<=self%aca_min_its)RETURN
full_size = INT(M,8)*INT(N,8)

!---Create temporary single element block
tmp_block%nelems=1
ALLOCATE(tmp_block%ielem(tmp_block%nelems))
tmp_block%ielem=1
ALLOCATE(tmp_block%icell(200))
tmp_block%icell=0
ALLOCATE(tmp_block%inv_map(self%tw_obj%mesh%np))
tmp_block%inv_map=0

! If we haven't converged before running for max_iter, we'll stop anyway.
! We want the sum of the two matrices to always be smaller than the dense one
max_iter = MIN(self%leaf_target*2,INT(full_size/INT(M+N,8),4)-1)

!
size_out=-1
ALLOCATE(US(N,max_iter),VS(max_iter,M))
us = 0.d0  ! The left vectors of the approximation
vs = 0.d0  ! The right vectors of the approximation
nterms = 0
ALLOCATE(prevIstar(N),prevJstar(M))
prevIstar = .TRUE.  ! Previously used i^* pivots
prevJstar = .TRUE.  ! Previously used j^* pivots

! Create a buffer for storing the R_{i^*,j} and R_{i, j^*}
ALLOCATE(RIstar(M,1),RJstar(N,1),RIref(M,1),RJref(N,1))
RIstar = 0.d0
RJstar = 0.d0
RIref = 0.d0
RJref = 0.d0

! Choose our starting random reference row and column.
CALL random_number(tmp)
Iref = INT(tmp*(N-1),4)+1
CALL random_number(tmp)
Jref = INT(tmp*(M-1),4)+1

! And collect the corresponding blocks of rows
! DO Iref=1,N
!     IF(prevIstar(Iref))EXIT
! END DO
tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=0
tmp_block%ielem=self%levels(ilevel)%blocks(iblock)%ielem(Iref)
tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=1
k1=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))
k2=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1))+1)-1
tmp_block%ncells=k2-k1+1
tmp_block%icell(1:tmp_block%ncells)=self%tw_obj%mesh%lpc(k1:k2)
CALL tw_compute_Lmatblock(self%tw_obj,self%tw_obj,RIref, &
  tmp_block,self%levels(ilevel)%blocks(jblock))
! DO i=1,nterms
!   RIref(:,1) = RIref(:,1) - us(Iref,i)*vs(i,:)
! END DO
!... and columns
! DO Jref=1,M
!     IF(prevJstar(Jref))EXIT
! END DO
tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=0
tmp_block%ielem=self%levels(ilevel)%blocks(jblock)%ielem(Jref)
tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=1
k1=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))
k2=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1))+1)-1
tmp_block%ncells=k2-k1+1
tmp_block%icell(1:tmp_block%ncells)=self%tw_obj%mesh%lpc(k1:k2)
CALL tw_compute_Lmatblock(self%tw_obj,self%tw_obj,RJref, &
  tmp_block,self%levels(ilevel)%blocks(iblock))
! DO i=1,nterms
!   RJref(:,1) = RJref(:,1) - us(:,i)*vs(i,Jref)
! END DO

step_size = -1.d0
DO k=1,max_iter
  ! Find the column in RIref with the largest entry (step 1 above).
  Jstar = max_masked(ABS(RIref(:,1)), prevJstar) !argmax_not_in_list(maxabsRIref, prevJstar)

  ! Find the row in RJref with the largest entry (step 1 above).
  Istar = max_masked(ABS(RJref(:,1)), prevIstar) !argmax_not_in_list(maxabsRJref, prevIstar)

  ! Check if we should pivot first based on row or based on column (step 2 above)
  Jstar_val = ABS(RIref(Jstar,1)) !maxabsRIref[Jstar]
  Istar_val = ABS(RJref(Istar,1)) !maxabsRJref[Istar]
  IF(Istar_val > Jstar_val)THEN
    ! If we pivot first on the row, then calculate the corresponding row
    ! of the residual matrix.
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=0
    tmp_block%ielem=self%levels(ilevel)%blocks(iblock)%ielem(Istar)
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=1
    k1=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))
    k2=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1))+1)-1
    tmp_block%ncells=k2-k1+1
    tmp_block%icell(1:tmp_block%ncells)=self%tw_obj%mesh%lpc(k1:k2)
    CALL tw_compute_Lmatblock(self%tw_obj,self%tw_obj,RIstar, &
      tmp_block,self%levels(ilevel)%blocks(jblock))
    DO i=1,nterms
      RIstar(:,1) = RIstar(:,1) - us(Istar,i)*vs(i,:)
    END DO

    ! Then find the largest entry in that row vector to identify which
    ! column to pivot on. (See step 3 above)
    Jstar = max_masked(ABS(RIstar(:,1)), prevJstar) !argmax_not_in_list(np.abs(RIstar), prevJstar)

    ! Calculate the corresponding residual column!
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=0
    tmp_block%ielem=self%levels(ilevel)%blocks(jblock)%ielem(Jstar)
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=1
    k1=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))
    k2=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1))+1)-1
    tmp_block%ncells=k2-k1+1
    tmp_block%icell(1:tmp_block%ncells)=self%tw_obj%mesh%lpc(k1:k2)
    CALL tw_compute_Lmatblock(self%tw_obj,self%tw_obj,RJstar, &
      tmp_block,self%levels(ilevel)%blocks(iblock))
    DO i=1,nterms
      RJstar(:,1) = RJstar(:,1) - us(:,i)*vs(i,Jstar)
    END DO
  ELSE
    ! If we pivot first on the column, then calculate the corresponding column
    ! of the residual matrix.
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=0
    tmp_block%ielem=self%levels(ilevel)%blocks(jblock)%ielem(Jstar)
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=1
    k1=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))
    k2=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1))+1)-1
    tmp_block%ncells=k2-k1+1
    tmp_block%icell(1:tmp_block%ncells)=self%tw_obj%mesh%lpc(k1:k2)
    CALL tw_compute_Lmatblock(self%tw_obj,self%tw_obj,RJstar, &
      tmp_block,self%levels(ilevel)%blocks(iblock))
    DO i=1,nterms
      RJstar(:,1) = RJstar(:,1) - us(:,i)*vs(i,Jstar)
    END DO

    ! Then find the largest entry in that row vector to identify which
    ! column to pivot on.  (See step 3 above)
    Istar = max_masked(ABS(RJstar(:,1)), prevIstar) !argmax_not_in_list(np.abs(RJstar), prevIstar)

    ! Calculate the corresponding residual row!
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=0
    tmp_block%ielem=self%levels(ilevel)%blocks(iblock)%ielem(Istar)
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=1
    k1=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))
    k2=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1))+1)-1
    tmp_block%ncells=k2-k1+1
    tmp_block%icell(1:tmp_block%ncells)=self%tw_obj%mesh%lpc(k1:k2)
    CALL tw_compute_Lmatblock(self%tw_obj,self%tw_obj,RIstar, &
      tmp_block,self%levels(ilevel)%blocks(jblock))
    DO i=1,nterms
      RIstar(:,1) = RIstar(:,1) - us(Istar,i)*vs(i,:)
    END DO
  END IF
  ! Record the pivot row and column so that we don't re-use them.
  prevIstar(Istar) = .FALSE. !prevIstar.append(Istar)
  prevJstar(Jstar) = .FALSE. !prevJstar.append(Jstar)

  ! Add the new rank-1 outer product to the approximation (see step 4 above)
  nterms = nterms + 1
  vs(nterms,:) = RIstar(:,1)/RIstar(Jstar,1) !vs.append(RIstar / RIstar[Jstar])
  us(:,nterms) = RJstar(:,1) !us.append(RJstar.copy())

  ! How "large" was this update to the approximation?
  step_size = SQRT(SUM(us(:,nterms)**2)*SUM(vs(nterms,:)**2))
  IF(oft_debug_print(2))WRITE(*,*)k,Istar,Jstar,step_size,tol

  ! The convergence criteria will simply be whether the Frobenius norm of the
  ! step is smaller than the user provided tolerance.
  IF((k>=self%aca_min_its).AND.(step_size < tol))EXIT

  ! We also break here if this is the last iteration to avoid wasting effort
  ! updating the reference row/column
  IF(k == max_iter - 1)THEN
    step_size = -1.d0
    EXIT !None, None
  END IF

  !! If we didn't converge, let's prep the reference residual row and column for the next iteration:

  ! If we pivoted on the reference row, then choose a new reference row.
  IF(Iref==Istar)THEN
    DO Iref=1,N
      IF(prevIstar(Iref))EXIT
    END DO
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=0
    tmp_block%ielem=self%levels(ilevel)%blocks(iblock)%ielem(Iref)
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=1
    k1=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))
    k2=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1))+1)-1
    tmp_block%ncells=k2-k1+1
    tmp_block%icell(1:tmp_block%ncells)=self%tw_obj%mesh%lpc(k1:k2)
    CALL tw_compute_Lmatblock(self%tw_obj,self%tw_obj,RIref, &
      tmp_block,self%levels(ilevel)%blocks(jblock))
    DO i=1,nterms
      RIref(:,1) = RIref(:,1) - us(Iref,i)*vs(i,:)
    END DO
  ELSE
    ! If we didn't change the reference row of the residual matrix "R",
    ! update the row to account for the new components of the approximation.
    RIref(:,1) = RIref(:,1) - us(Iref,nterms)*vs(nterms,:) !RIref -= us[-1][Iref : Iref + 1][:, None] * vs[-1][None, :]
  END IF

  ! If we pivoted on the reference column, then choose a new reference column.
  if(Jref==Jstar)THEN
  DO Jref=1,M
    IF(prevJstar(Jref))EXIT
  END DO
  tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=0
  tmp_block%ielem=self%levels(ilevel)%blocks(jblock)%ielem(Jref)
  tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=1
  k1=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))
  k2=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1))+1)-1
  tmp_block%ncells=k2-k1+1
  tmp_block%icell(1:tmp_block%ncells)=self%tw_obj%mesh%lpc(k1:k2)
  CALL tw_compute_Lmatblock(self%tw_obj,self%tw_obj,RJref, &
    tmp_block,self%levels(ilevel)%blocks(iblock))
  DO i=1,nterms
    RJref(:,1) = RJref(:,1) - us(:,i)*vs(i,Jref)
  END DO
  else
  ! If we didn't change the reference column of the residual matrix "R",
  ! update the column to account for the new components of the approximation.
  RJref(:,1) = RJref(:,1) - us(:,nterms)*vs(nterms,Jref) !RJref -= vs[-1][Jref : Jref + 1][None, :] * us[-1][:, None]
  END IF
END DO
IF(step_size > 0.d0)THEN
  ! WRITE(*,*)'Savings:   ',iblock,jblock,nterms,MIN(M,N)
  ! Return the left and right approximation matrices.
  ! The approximate is such that:
  ! M ~ U_ACA.dot(V_ACA)
  ALLOCATE(self%aca_U_mats(isparse)%M(M,nterms))
  DO i=1,nterms
    DO k=1,M
      self%aca_U_mats(isparse)%M(k,i)=vs(i,k)
    END DO
  END DO
  ALLOCATE(self%aca_V_mats(isparse)%M(nterms,N))
  DO i=1,nterms
    DO k=1,N
      self%aca_V_mats(isparse)%M(i,k)=us(k,i)
    END DO
  END DO
  ! DEALLOCATE(self%sparse_dense(isparse)%M)
  size_out=INT(nterms,8)*INT(M+N,8)
END IF
DEALLOCATE(RIstar,RJstar,RIref,RJref)
DEALLOCATE(us,vs,prevIstar,prevJstar)
DEALLOCATE(tmp_block%ielem,tmp_block%icell,tmp_block%inv_map)
! return U_ACA, V_ACA
END SUBROUTINE aca_approx
!
SUBROUTINE compress_aca(isparse,tol,size_out)
INTEGER(4), INTENT(in) :: isparse
REAL(8), INTENT(in) :: tol
INTEGER(8), INTENT(out) :: size_out
INTEGER(4) :: i,j,io_unit
REAL(8) :: S_sum,frob_K
!
INTEGER(4) :: ilevel,iblock,jblock
INTEGER(4) :: M,N,K,LDA,LDU,LDVT,LWORK,INFO,MIN_DIM
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: IWORK
REAL(8), ALLOCATABLE, DIMENSION(:) :: S,WORK
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: U,VT,Atmp,QU,QV,RU,RV
!
ilevel=self%sparse_blocks(1,isparse)
iblock=self%sparse_blocks(2,isparse)
jblock=self%sparse_blocks(3,isparse)
! Get QR factorization of U and V matrices
! WRITE(*,*)'QR U'
CALL get_qr(self%aca_U_mats(isparse)%M,QU,RU)
! WRITE(*,*)'QR VT'
ALLOCATE(Atmp(SIZE(self%aca_V_mats(isparse)%M,2),SIZE(self%aca_V_mats(isparse)%M,1)))
Atmp = TRANSPOSE(self%aca_V_mats(isparse)%M)
CALL get_qr(Atmp,QV,RV)
DEALLOCATE(Atmp)
M=SIZE(QU,DIM=1)
N=SIZE(QV,DIM=1)
K=SIZE(RU,DIM=1)
! WRITE(*,*)K
ALLOCATE(Atmp(K,K))
CALL DGEMM('N', 'T', K, K, K, 1.d0, RU, K, RV, K, 0.d0, Atmp, K)
!
LDA = K
LDU = K
MIN_DIM = K !MIN(M,N)
LDVT = MIN_DIM
ALLOCATE(S(MIN_DIM),U(LDU,MIN_DIM),VT(LDVT,K))
!---Get optimal workspace size
LWORK = -1
ALLOCATE(WORK(1),IWORK(8*MIN_DIM))
CALL DGESDD('S', K, K, Atmp, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
IF( INFO.NE.0 ) THEN
  CALL oft_abort("The algorithm computing SVD failed to converge.","tw_Lmat_MF_setup::compress_aca",__FILE__)
END IF
LWORK = INT(WORK(1))
DEALLOCATE(WORK)
!---Compute SVD of matrix block
ALLOCATE(WORK(LWORK))
CALL DGESDD('S', K, K, Atmp, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
IF( INFO.NE.0 ) THEN
  CALL oft_abort("The algorithm computing SVD failed to converge.","tw_Lmat_MF_setup::compress_aca",__FILE__)
END IF
!---Truncate at desired accuracy
frob_K=0.d0
DO i=MIN_DIM,1,-1
  frob_K = frob_K + S(i)**2
  IF(SQRT(frob_K)>tol)EXIT
END DO
i=MAX(MIN(self%min_rank,MIN_DIM),i)
IF((i==MIN_DIM).AND.(MIN_DIM>self%aca_min_its))THEN
  i=MIN_DIM
  WRITE(*,*)'SVD recompression failed'
END IF
!i=MIN(i,MIN_DIM)
!---Replace U and V matrices
DEALLOCATE(self%aca_U_mats(isparse)%M,self%aca_V_mats(isparse)%M)
ALLOCATE(self%aca_U_mats(isparse)%M(self%levels(ilevel)%blocks(jblock)%nelems,i))
ALLOCATE(self%aca_V_mats(isparse)%M(i,self%levels(ilevel)%blocks(iblock)%nelems))
! U = UQ.dot(U[:, :r] * S[:r])
! V = VT[:r, :].dot(VQ.T)
DO j=1,i
  U(:,j) = U(:,j)*S(j)
END DO
CALL DGEMM('N', 'N', M, i, K, 1.d0, QU, M, U, K, 0.d0, self%aca_U_mats(isparse)%M, M)
!self%aca_U_mats(isparse)%M = MATMUL(QU,U(:,1:i))
CALL DGEMM('N', 'T', i, N, K, 1.d0, VT, K, QV, N, 0.d0, self%aca_V_mats(isparse)%M, i)
!self%aca_V_mats(isparse)%M = MATMUL(VT(1:i,:),TRANSPOSE(QV))
size_out=INT(i,8)*INT(M+N,8)
! WRITE(*,*)'Compressed:',iblock,jblock,i,MIN(M,N)
DEALLOCATE(WORK,IWORK,S,U,VT,Atmp,QU,QV,RU,RV)
END SUBROUTINE compress_aca
END SUBROUTINE tw_hodlr_Lcompute
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE tw_hodlr_Bcompute(self,save_file)
class(oft_tw_hodlr_op), intent(inout) :: self
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: save_file
INTEGER(4) :: i,ii,j,k,l,nlast,branch_count,pind_i,pind_j,max_pts,top_level,level,dim,flip
INTEGER(4) :: nblocks_progress(9),nblocks_complete,counts(4),iblock,jblock
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: kc_tmp,kr_tmp
INTEGER(4), ALLOCATABLE, DIMENSION(:,:) :: block_interaction,row_count
INTEGER(8) :: mat_size(3),compressed_size,size_out
LOGICAL :: is_masked,is_neighbor,mat_updated
LOGICAL, ALLOCATABLE, DIMENSION(:) :: near_closure,cell_mark
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: dense_blocks
REAL(8) :: corners_i(3,8),corners_j(3,8),rmin(3),elapsed_time,avg_size
REAL(8), ALLOCATABLE, DIMENSION(:) :: point_block
TYPE(oft_native_dense_matrix), POINTER :: mat_tmp => NULL()
type(oft_timer) :: mytimer
WRITE(*,*)'Building block low rank magnetic field operator'
IF(self%B_svd_tol<0.d0)CALL oft_abort('"B_svd_tol" must be > 0','tw_Lmat_MF_Bcompute',__FILE__)
!---Build matrix with SVD compression of off-diagonal blocks
! approximation with ACA+ planned
CALL mytimer%tick()
ALLOCATE(self%dense_B_mats(self%ndense,3))
ALLOCATE(self%aca_BU_mats(2*self%nsparse,3))
ALLOCATE(self%aca_BV_mats(2*self%nsparse,3))
ALLOCATE(self%aca_B_dense(2*self%nsparse,3))
!
WRITE(*,*)'  Building hole and Vcoil columns'
IF(MAX(self%tw_obj%n_icoils,self%tw_obj%nholes+self%tw_obj%n_vcoils)>0)THEN
  ALLOCATE(self%hole_Vcoil_Bmat(self%tw_obj%mesh%np,MAX(1,self%tw_obj%nholes+self%tw_obj%n_vcoils),3))
  ALLOCATE(self%Icoil_Bmat(self%tw_obj%mesh%np,MAX(1,self%tw_obj%n_icoils),3))
  CALL tw_compute_Bops_hole(self%tw_obj,self%hole_Vcoil_Bmat,self%Icoil_Bmat)
  self%tw_obj%Bdr=>self%Icoil_Bmat
END IF
!
IF(PRESENT(save_file))CALL read_from_file()
WRITE(*,*)'  Building diagonal blocks'
compressed_size=0
mat_updated=.FALSE.
!$omp parallel private(i,ii,level,j,k,size_out,avg_size,dim,flip) reduction(+:compressed_size) &
!$omp reduction(.OR.:mat_updated)
!$omp single
nblocks_progress = INT(self%ndense*[0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0],4)
nblocks_complete = 0
!$omp end single
!$omp do schedule(static,1)
DO i=1,self%ndense
  level = self%dense_blocks(1,i)
  j = self%dense_blocks(2,i)
  k = self%dense_blocks(3,i)
  !---Build full representation
  DO dim=1,3
    IF(.NOT.ASSOCIATED(self%dense_B_mats(i,dim)%M))THEN
      mat_updated=.TRUE.
      ALLOCATE(self%dense_B_mats(i,dim)%M(self%levels(level)%blocks(j)%np,self%levels(level)%blocks(k)%nelems))
      CALL tw_compute_Bops_block(self%tw_obj,self%dense_B_mats(i,dim)%M, &
        self%levels(level)%blocks(k),self%levels(level)%blocks(j),dim)
    END IF
  END DO
  compressed_size = compressed_size + 3*self%levels(level)%blocks(j)%np*self%levels(level)%blocks(k)%nelems
  !$omp critical
  nblocks_complete = nblocks_complete + 1
  DO k=1,9
    IF(nblocks_complete==nblocks_progress(k))WRITE(*,'(5X,I2,A1)')k*10,'%'
  END DO
  !$omp end critical
END DO
!$omp single
IF(self%B_aca_tol>0.d0)THEN
  WRITE(*,*)'  Building off-diagonal blocks using ACA+'
ELSE
  WRITE(*,*)'  Building off-diagonal blocks with SVD compression'
END IF
nblocks_progress = INT(self%nsparse*[0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0],4)
nblocks_complete = 0
!$omp end single
!$omp do schedule(dynamic,1)
DO i=1,self%nsparse
  level = self%sparse_blocks(1,i)
  j = self%sparse_blocks(2,i)
  k = self%sparse_blocks(3,i)
  DO flip=1,2
    ii=(i-1)*2+flip
    IF(flip==2)THEN
      k = self%sparse_blocks(2,i)
      j = self%sparse_blocks(3,i)
    END IF
    DO dim=1,3
      IF(ASSOCIATED(self%aca_B_dense(ii,dim)%M))THEN
        size_out=self%levels(level)%blocks(k)%nelems*self%levels(level)%blocks(j)%nelems
      ELSE IF(ASSOCIATED(self%aca_BU_mats(ii,dim)%M))THEN
        size_out=SIZE(self%aca_BV_mats(ii,dim)%M,DIM=1)*(self%levels(level)%blocks(k)%nelems+self%levels(level)%blocks(j)%nelems)
      ELSE
        mat_updated=.TRUE.
        !---Build full representation
        size_out=-1
        IF((self%B_aca_tol>0.d0).AND.ABS(self%levels(level)%mat_mask(j,k))==2)THEN
        CALL aca_approx(ii,level,j,k,dim,self%B_aca_tol,size_out)
        IF(size_out>0)CALL compress_aca(ii,level,j,k,dim,self%B_svd_tol,size_out)
        IF(size_out==-1)WRITE(*,*)'ACA+ failure',self%levels(level)%blocks(j)%np, &
          self%levels(level)%blocks(k)%nelems
        END IF
        !---IF ACA+ failed or diagonal compute dense matrix
        IF(size_out<0)THEN
          ALLOCATE(self%aca_B_dense(ii,dim)%M(self%levels(level)%blocks(k)%np,self%levels(level)%blocks(j)%nelems))
          CALL tw_compute_Bops_block(self%tw_obj,self%aca_B_dense(ii,dim)%M, &
            self%levels(level)%blocks(j),self%levels(level)%blocks(k),dim)
          CALL compress_block(ii,j,k,dim,self%B_svd_tol,size_out)
          ! size_out = self%levels(level)%blocks(j)%np*self%levels(level)%blocks(k)%nelems
        END IF
      END IF
      compressed_size = compressed_size + size_out
    END DO
  END DO
  !$omp critical
  nblocks_complete = nblocks_complete + 1
  DO k=1,9
    IF(nblocks_complete==nblocks_progress(k))WRITE(*,'(5X,I2,A1)')k*10,'%'
  END DO
  !$omp end critical
END DO
!$omp end parallel
elapsed_time=mytimer%tock()
WRITE(*,'(5X,A,F6.1,A,ES9.2,A,ES9.2,A)')'Compression ratio:', &
  compressed_size*1.d2/(3*REAL(self%tw_obj%np_active,8)*REAL(self%tw_obj%mesh%np,8)), &
  "%  (",REAL(compressed_size,8),"/",3*REAL(self%tw_obj%np_active,8)*REAL(self%tw_obj%mesh%np,8),")"
WRITE(*,'(5X,2A)')'Time = ',time_to_string(elapsed_time)
!
IF(mat_updated.AND.PRESENT(save_file))CALL save_to_file()
CONTAINS
SUBROUTINE save_to_file()
INTEGER(4) :: ierr,io_unit,hash_tmp(6),file_counts(6)
CHARACTER(LEN=1) :: dim_id
CHARACTER(LEN=5) :: matrix_id
!---Save computed matrix to file
IF(TRIM(save_file)/='none')THEN
  hash_tmp(1) = self%tw_obj%nelems
  hash_tmp(2) = self%tw_obj%mesh%nc
  hash_tmp(3) = self%ndense
  hash_tmp(4) = self%nsparse
  hash_tmp(5) = oft_simple_hash(C_LOC(self%tw_obj%mesh%lc),INT(4*3*self%tw_obj%mesh%nc,8))
  hash_tmp(6) = oft_simple_hash(C_LOC(self%tw_obj%mesh%r),INT(8*3*self%tw_obj%mesh%np,8))
  WRITE(*,*)'  Saving HODLR matrix from file: ',TRIM(save_file)
  CALL hdf5_create_file(TRIM(save_file))
  CALL hdf5_write(hash_tmp,TRIM(save_file),'MODEL_hash')
  !
  hash_tmp(6)=0
  DO i=1,self%ndense
    level = self%dense_blocks(1,i)
    j = self%dense_blocks(2,i)
    k = self%dense_blocks(3,i)
    hash_tmp(1:3) = [level,j,k]
    hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ipts),INT(4*self%levels(level)%blocks(j)%np,8))
    hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ielem),INT(4*self%levels(level)%blocks(k)%nelems,8))
    WRITE(matrix_id,'(I5)')i
    CALL hdf5_write(hash_tmp,TRIM(save_file),'DENSE_hash_'//TRIM(ADJUSTL(matrix_id)))
    DO dim=1,3
      WRITE(dim_id,'(I1)')dim
      CALL hdf5_write(self%dense_B_mats(i,dim)%M,TRIM(save_file),'DENSE_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id)
    END DO
  END DO
  !
  DO i=1,self%nsparse
    level = self%sparse_blocks(1,i)
    j = self%sparse_blocks(2,i)
    k = self%sparse_blocks(3,i)
    DO flip=1,2
      ii=(i-1)*2+flip
      IF(flip==2)THEN
        k = self%sparse_blocks(2,i)
        j = self%sparse_blocks(3,i)
      END IF
      WRITE(matrix_id,'(I5)')ii
      DO dim=1,3
        hash_tmp(1:3) = [level,j,k]
        hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ielem),INT(4*self%levels(level)%blocks(j)%nelems,8))
        hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ipts),INT(4*self%levels(level)%blocks(k)%np,8))
        WRITE(dim_id,'(I1)')dim
        IF(ASSOCIATED(self%aca_B_dense(ii,dim)%M))THEN
          hash_tmp(6)=-1
          CALL hdf5_write(hash_tmp,TRIM(save_file),'SPARSE_hash_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id)
          CALL hdf5_write(self%aca_B_dense(ii,dim)%M,TRIM(save_file),'SPARSE_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id)
        ELSE
          hash_tmp(6)=SIZE(self%aca_BV_mats(ii,dim)%M,DIM=1)
          CALL hdf5_write(hash_tmp,TRIM(save_file),'SPARSE_hash_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id)
          CALL hdf5_write(self%aca_BU_mats(ii,dim)%M,TRIM(save_file),'SPARSE_U_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id)
          CALL hdf5_write(self%aca_BV_mats(ii,dim)%M,TRIM(save_file),'SPARSE_V_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id)
        END IF
      END DO
    END DO
  END DO
  ! TODO: Save create appropriate hash to enable saving coil matrices
  ! WRITE(io_unit)hash_tmp
  ! WRITE(io_unit)self%hole_Vcoil_mat%M
  ! CLOSE(io_unit)
END IF
END SUBROUTINE save_to_file
SUBROUTINE read_from_file()
INTEGER(4) :: ierr,io_unit,hash_tmp(6),file_counts(6)
LOGICAL :: exists
CHARACTER(LEN=1) :: dim_id
CHARACTER(LEN=5) :: matrix_id
!---Try to load from file
IF(TRIM(save_file)/='none')THEN
  INQUIRE(FILE=TRIM(save_file),EXIST=exists)
  IF(exists)THEN
    hash_tmp(1) = self%tw_obj%nelems
    hash_tmp(2) = self%tw_obj%mesh%nc
    hash_tmp(3) = self%ndense
    hash_tmp(4) = self%nsparse
    hash_tmp(5) = oft_simple_hash(C_LOC(self%tw_obj%mesh%lc),INT(4*3*self%tw_obj%mesh%nc,8))
    hash_tmp(6) = oft_simple_hash(C_LOC(self%tw_obj%mesh%r),INT(8*3*self%tw_obj%mesh%np,8))
    WRITE(*,*)'  Reading HODLR matrix from file: ',TRIM(save_file)
    CALL hdf5_read(file_counts,TRIM(save_file),'MODEL_hash',success=exists)
    ! OPEN(NEWUNIT=io_unit,FILE=TRIM(save_file),FORM='UNFORMATTED')
    ! READ(io_unit, IOSTAT=ierr)file_counts
    IF(exists.AND.ALL(file_counts==hash_tmp))THEN
      hash_tmp(6)=0
      DO i=1,self%ndense
        level = self%dense_blocks(1,i)
        j = self%dense_blocks(2,i)
        k = self%dense_blocks(3,i)
        hash_tmp(1:3) = [level,j,k]
        hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ipts),INT(4*self%levels(level)%blocks(j)%np,8))
        hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ielem),INT(4*self%levels(level)%blocks(k)%nelems,8))
        WRITE(matrix_id,'(I5)')i
        CALL hdf5_read(file_counts,TRIM(save_file),'DENSE_hash_'//TRIM(ADJUSTL(matrix_id)),success=exists)
        ! READ(io_unit, IOSTAT=ierr)file_counts
        IF(exists.AND.ALL(file_counts==hash_tmp))THEN
          DO dim=1,3
            ALLOCATE(self%dense_B_mats(i,dim)%M(self%levels(level)%blocks(j)%np,self%levels(level)%blocks(k)%nelems))
            WRITE(dim_id,'(I1)')dim
            CALL hdf5_read(self%dense_B_mats(i,dim)%M,TRIM(save_file),'DENSE_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id,success=exists)
            ! READ(io_unit, IOSTAT=ierr)self%dense_B_mats(i,dim)%M
            IF(.NOT.exists)DEALLOCATE(self%dense_B_mats(i,dim)%M)
          END DO
        END IF
      END DO
      !
      DO i=1,self%nsparse
        level = self%sparse_blocks(1,i)
        j = self%sparse_blocks(2,i)
        k = self%sparse_blocks(3,i)
        DO flip=1,2
          ii=(i-1)*2+flip
          IF(flip==2)THEN
            k = self%sparse_blocks(2,i)
            j = self%sparse_blocks(3,i)
          END IF
          WRITE(matrix_id,'(I5)')ii
          DO dim=1,3
            WRITE(dim_id,'(I1)')dim
            hash_tmp(1:3) = [level,j,k]
            hash_tmp(4) = oft_simple_hash(C_LOC(self%levels(level)%blocks(j)%ielem),INT(4*self%levels(level)%blocks(j)%nelems,8))
            hash_tmp(5) = oft_simple_hash(C_LOC(self%levels(level)%blocks(k)%ipts),INT(4*self%levels(level)%blocks(k)%np,8))
            CALL hdf5_read(file_counts,TRIM(save_file),'SPARSE_hash_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id,success=exists)
            ! READ(io_unit, IOSTAT=ierr)file_counts
            IF(exists.AND.ALL(file_counts(1:5)==hash_tmp(1:5)))THEN
              IF(file_counts(6)<0)THEN
                ALLOCATE(self%aca_B_dense(ii,dim)%M(self%levels(level)%blocks(k)%np,self%levels(level)%blocks(j)%nelems))
                CALL hdf5_read(self%aca_B_dense(ii,dim)%M,TRIM(save_file),'SPARSE_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id,success=exists)
                ! READ(io_unit, IOSTAT=ierr)self%aca_B_dense(ii,dim)%M
                IF(.NOT.exists)DEALLOCATE(self%aca_B_dense(ii,dim)%M)
              ELSE
                ALLOCATE(self%aca_BU_mats(ii,dim)%M(self%levels(level)%blocks(k)%np,file_counts(6)))
                CALL hdf5_read(self%aca_BU_mats(ii,dim)%M,TRIM(save_file),'SPARSE_U_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id,success=exists)
                ! READ(io_unit, IOSTAT=ierr)self%aca_BU_mats(ii,dim)%M
                IF(.NOT.exists)THEN
                  DEALLOCATE(self%aca_BU_mats(ii,dim)%M)
                ELSE
                  ALLOCATE(self%aca_BV_mats(ii,dim)%M(file_counts(6),self%levels(level)%blocks(j)%nelems))
                  CALL hdf5_read(self%aca_BV_mats(ii,dim)%M,TRIM(save_file),'SPARSE_V_'//TRIM(ADJUSTL(matrix_id))//'_'//dim_id,success=exists)
                  ! READ(io_unit, IOSTAT=ierr)self%aca_BV_mats(ii,dim)%M
                  IF(.NOT.exists)DEALLOCATE(self%aca_BU_mats(ii,dim)%M,self%aca_BV_mats(ii,dim)%M)
                END IF
              END IF
            END IF
          END DO
        END DO
      END DO
      ! CLOSE(io_unit)
    ELSE
      WRITE(*,*)'  Ignoring stored matrix: Model hashes do not match'
      ! CLOSE(io_unit)
    END IF
  END IF
END IF
END SUBROUTINE read_from_file
!
SUBROUTINE compress_block(isparse,iblock,jblock,dim,tol,size_out)
INTEGER(4), INTENT(in) :: isparse,iblock,jblock,dim
REAL(8), INTENT(in) :: tol
INTEGER(8), INTENT(out) :: size_out
INTEGER(4) :: i,j,io_unit
REAL(8) :: S_sum,frob_K
INTEGER(4) :: M,N,LDA,LDU,LDVT,LWORK,INFO,MIN_DIM
INTEGER(8) :: full_size
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: IWORK
REAL(8), ALLOCATABLE, DIMENSION(:) :: S,WORK
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: U,VT,Atmp
!
M = SIZE(self%aca_B_dense(isparse,dim)%M,DIM=1)
N = SIZE(self%aca_B_dense(isparse,dim)%M,DIM=2)
full_size = INT(M,8)*INT(N,8)
LDA = M
LDU = M
MIN_DIM = min(M,N)
LDVT = MIN_DIM
ALLOCATE(Atmp(M,N))
Atmp = self%aca_B_dense(isparse,dim)%M
ALLOCATE(S(MIN_DIM),U(LDU,MIN_DIM),VT(LDVT,N))
!---Get optimal workspace size
LWORK = -1
ALLOCATE(WORK(1),IWORK(8*MIN_DIM))
CALL DGESDD('S', M, N, Atmp, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
LWORK = INT(WORK(1))
DEALLOCATE(WORK)
!---Compute SVD of matrix block
ALLOCATE(WORK(LWORK))
CALL DGESDD('S', M, N, Atmp, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
!---Truncate at desired accuracy
frob_K=0.d0
DO i=MIN_DIM,1,-1
  frob_K = frob_K + S(i)**2
  IF(SQRT(frob_K)>tol)EXIT
END DO
i=MAX(MIN(self%min_rank,MIN_DIM),i)
IF(i<full_size/INT(M+N,8))THEN
  ALLOCATE(self%aca_BU_mats(isparse,dim)%M(M,i))
  self%aca_BU_mats(isparse,dim)%M=U(:,1:i)
  ALLOCATE(self%aca_BV_mats(isparse,dim)%M(i,N))
  DO j=1,i
      self%aca_BV_mats(isparse,dim)%M(j,:)=VT(j,:)*S(j)
  END DO
  DEALLOCATE(self%aca_B_dense(isparse,dim)%M)
  ! WRITE(*,*)'Savings SVD:',isparse,dim,i,M,N
  size_out=INT(i,8)*INT(M+N,8)
ELSE
  ! WRITE(*,*)'No savings:',iblock,jblock
  size_out=full_size
END IF
!
IF( INFO.NE.0 ) THEN
  CALL oft_abort("The algorithm computing SVD failed to converge.","tw_Lmat_MF_Bsetup::compress_block",__FILE__)
END IF
DEALLOCATE(WORK,IWORK,S,U,VT,Atmp)
END SUBROUTINE compress_block
!
FUNCTION max_masked(vals,mask) RESULT(iloc)
REAL(8), INTENT(in) :: vals(:)
LOGICAL, INTENT(in) :: mask(:)
REAL(8) :: max_tmp
INTEGER(4) :: i,iloc
max_tmp = -1.d99
DO i=1,SIZE(vals)
  IF((vals(i)>max_tmp).AND.mask(i))THEN
    max_tmp=vals(i)
    iloc=i
  END IF
END DO
END FUNCTION max_masked
!
SUBROUTINE aca_approx(isparse,ilevel,iblock,jblock,dim,tol,size_out)
INTEGER(4), INTENT(in) :: isparse,ilevel,iblock,jblock,dim
REAL(8), INTENT(in) :: tol
INTEGER(8), INTENT(out) :: size_out
INTEGER(4) :: M,N,MIN_DIM,max_iter,Iref,Jref,Istar,Jstar,i,k,nterms,k1,k2
INTEGER(8) :: full_size
REAL(8) :: tmp,step_size,Istar_val,Jstar_val
REAL(8), ALLOCATABLE :: us(:,:),vs(:,:),RIstar(:,:),RJstar(:,:),RIref(:,:),RJref(:,:)
LOGICAL, ALLOCATABLE :: prevIstar(:),prevJstar(:)
TYPE(oft_tw_block) :: tmp_block
!
size_out=-2
M = self%levels(ilevel)%blocks(jblock)%np
N = self%levels(ilevel)%blocks(iblock)%nelems
MIN_DIM = min(M,N)
IF(MIN_DIM<=self%aca_min_its)RETURN
full_size = INT(M,8)*INT(N,8)

!---Create temporary single element block
tmp_block%np=1
ALLOCATE(tmp_block%ipts(tmp_block%np))
tmp_block%ipts=1
!
tmp_block%nelems=1
ALLOCATE(tmp_block%ielem(tmp_block%nelems))
tmp_block%ielem=1
ALLOCATE(tmp_block%icell(200))
tmp_block%icell=0
ALLOCATE(tmp_block%inv_map(self%tw_obj%mesh%np))
tmp_block%inv_map=0

! If we haven't converged before running for max_iter, we'll stop anyway.
! We want the sum of the two matrices to always be smaller than the dense one
max_iter = MIN(self%leaf_target*2,INT(full_size/INT(M+N,8),4)-1)

!
size_out=-1
ALLOCATE(US(N,max_iter),VS(max_iter,M))
us = 0.d0  ! The left vectors of the approximation
vs = 0.d0  ! The right vectors of the approximation
nterms = 0
ALLOCATE(prevIstar(N),prevJstar(M))
prevIstar = .TRUE.  ! Previously used i^* pivots
prevJstar = .TRUE.  ! Previously used j^* pivots

! Create a buffer for storing the R_{i^*,j} and R_{i, j^*}
ALLOCATE(RIstar(M,1),RJstar(1,N),RIref(M,1),RJref(1,N))
RIstar = 0.d0
RJstar = 0.d0
RIref = 0.d0
RJref = 0.d0

! Choose our starting random reference row and column.
CALL random_number(tmp)
Iref = INT(tmp*(N-1),4)+1
CALL random_number(tmp)
Jref = INT(tmp*(M-1),4)+1

! And collect the corresponding blocks of rows
! DO Iref=1,N
!     IF(prevIstar(Iref))EXIT
! END DO
tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=0
tmp_block%ielem=self%levels(ilevel)%blocks(iblock)%ielem(Iref)
tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=1
k1=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))
k2=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1))+1)-1
tmp_block%ncells=k2-k1+1
tmp_block%icell(1:tmp_block%ncells)=self%tw_obj%mesh%lpc(k1:k2)
CALL tw_compute_Bops_block(self%tw_obj,RIref, &
  tmp_block,self%levels(ilevel)%blocks(jblock),dim)
! DO i=1,nterms
!   RIref(:,1) = RIref(:,1) - us(Iref,i)*vs(i,:)
! END DO
!... and columns
! DO Jref=1,M
!     IF(prevJstar(Jref))EXIT
! END DO
tmp_block%ipts=self%levels(ilevel)%blocks(jblock)%ipts(Jref)
CALL tw_compute_Bops_block(self%tw_obj,RJref, &
  self%levels(ilevel)%blocks(iblock),tmp_block,dim)
! DO i=1,nterms
!   RJref(:,1) = RJref(:,1) - us(:,i)*vs(i,Jref)
! END DO

step_size = -1.d0
DO k=1,max_iter
  ! Find the column in RIref with the largest entry (step 1 above).
  Jstar = max_masked(ABS(RIref(:,1)), prevJstar) !argmax_not_in_list(maxabsRIref, prevJstar)

  ! Find the row in RJref with the largest entry (step 1 above).
  Istar = max_masked(ABS(RJref(1,:)), prevIstar) !argmax_not_in_list(maxabsRJref, prevIstar)

  ! Check if we should pivot first based on row or based on column (step 2 above)
  Jstar_val = ABS(RIref(Jstar,1)) !maxabsRIref[Jstar]
  Istar_val = ABS(RJref(1,Istar)) !maxabsRJref[Istar]
  IF(Istar_val > Jstar_val)THEN
    ! If we pivot first on the row, then calculate the corresponding row
    ! of the residual matrix.
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=0
    tmp_block%ielem=self%levels(ilevel)%blocks(iblock)%ielem(Istar)
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=1
    k1=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))
    k2=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1))+1)-1
    tmp_block%ncells=k2-k1+1
    tmp_block%icell(1:tmp_block%ncells)=self%tw_obj%mesh%lpc(k1:k2)
    CALL tw_compute_Bops_block(self%tw_obj,RIstar, &
      tmp_block,self%levels(ilevel)%blocks(jblock),dim)
    DO i=1,nterms
      RIstar(:,1) = RIstar(:,1) - us(Istar,i)*vs(i,:)
    END DO

    ! Then find the largest entry in that row vector to identify which
    ! column to pivot on. (See step 3 above)
    Jstar = max_masked(ABS(RIstar(:,1)), prevJstar) !argmax_not_in_list(np.abs(RIstar), prevJstar)

    ! Calculate the corresponding residual column!
    tmp_block%ipts=self%levels(ilevel)%blocks(jblock)%ipts(Jstar)
    CALL tw_compute_Bops_block(self%tw_obj,RJstar, &
      self%levels(ilevel)%blocks(iblock),tmp_block,dim)
    DO i=1,nterms
      RJstar(1,:) = RJstar(1,:) - us(:,i)*vs(i,Jstar)
    END DO
  ELSE
    ! If we pivot first on the column, then calculate the corresponding column
    ! of the residual matrix.
    tmp_block%ipts=self%levels(ilevel)%blocks(jblock)%ipts(Jstar)
    CALL tw_compute_Bops_block(self%tw_obj,RJstar, &
      self%levels(ilevel)%blocks(iblock),tmp_block,dim)
    DO i=1,nterms
      RJstar(1,:) = RJstar(1,:) - us(:,i)*vs(i,Jstar)
    END DO

    ! Then find the largest entry in that row vector to identify which
    ! column to pivot on.  (See step 3 above)
    Istar = max_masked(ABS(RJstar(1,:)), prevIstar) !argmax_not_in_list(np.abs(RJstar), prevIstar)

    ! Calculate the corresponding residual row!
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=0
    tmp_block%ielem=self%levels(ilevel)%blocks(iblock)%ielem(Istar)
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=1
    k1=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))
    k2=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1))+1)-1
    tmp_block%ncells=k2-k1+1
    tmp_block%icell(1:tmp_block%ncells)=self%tw_obj%mesh%lpc(k1:k2)
    CALL tw_compute_Bops_block(self%tw_obj,RIstar, &
      tmp_block,self%levels(ilevel)%blocks(jblock),dim)
    DO i=1,nterms
      RIstar(:,1) = RIstar(:,1) - us(Istar,i)*vs(i,:)
    END DO
  END IF
  ! Record the pivot row and column so that we don't re-use them.
  prevIstar(Istar) = .FALSE. !prevIstar.append(Istar)
  prevJstar(Jstar) = .FALSE. !prevJstar.append(Jstar)

  ! Add the new rank-1 outer product to the approximation (see step 4 above)
  nterms = nterms + 1
  vs(nterms,:) = RIstar(:,1)/RIstar(Jstar,1) !vs.append(RIstar / RIstar[Jstar])
  us(:,nterms) = RJstar(1,:) !us.append(RJstar.copy())

  ! How "large" was this update to the approximation?
  step_size = SQRT(SUM(us(:,nterms)**2)*SUM(vs(nterms,:)**2))
  IF(oft_debug_print(2))WRITE(*,*)k,Istar,Jstar,step_size,tol

  ! The convergence criteria will simply be whether the Frobenius norm of the
  ! step is smaller than the user provided tolerance.
  IF((k>=self%aca_min_its).AND.(step_size < tol))EXIT

  ! We also break here if this is the last iteration to avoid wasting effort
  ! updating the reference row/column
  IF(k == max_iter - 1)THEN
    step_size = -1.d0
    EXIT !None, None
  END IF

  !! If we didn't converge, let's prep the reference residual row and column for the next iteration:

  ! If we pivoted on the reference row, then choose a new reference row.
  IF(Iref==Istar)THEN
    DO Iref=1,N
      IF(prevIstar(Iref))EXIT
    END DO
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=0
    tmp_block%ielem=self%levels(ilevel)%blocks(iblock)%ielem(Iref)
    tmp_block%inv_map(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))=1
    k1=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1)))
    k2=self%tw_obj%mesh%kpc(self%tw_obj%lpmap_inv(tmp_block%ielem(1))+1)-1
    tmp_block%ncells=k2-k1+1
    tmp_block%icell(1:tmp_block%ncells)=self%tw_obj%mesh%lpc(k1:k2)
    CALL tw_compute_Bops_block(self%tw_obj,RIref, &
      tmp_block,self%levels(ilevel)%blocks(jblock),dim)
    DO i=1,nterms
      RIref(:,1) = RIref(:,1) - us(Iref,i)*vs(i,:)
    END DO
  ELSE
    ! If we didn't change the reference row of the residual matrix "R",
    ! update the row to account for the new components of the approximation.
    RIref(:,1) = RIref(:,1) - us(Iref,nterms)*vs(nterms,:) !RIref -= us[-1][Iref : Iref + 1][:, None] * vs[-1][None, :]
  END IF

  ! If we pivoted on the reference column, then choose a new reference column.
  IF(Jref==Jstar)THEN
    DO Jref=1,M
        IF(prevJstar(Jref))EXIT
    END DO
    tmp_block%ipts=self%levels(ilevel)%blocks(jblock)%ipts(Jref)
    CALL tw_compute_Bops_block(self%tw_obj,RJref, &
      self%levels(ilevel)%blocks(iblock),tmp_block,dim)
    DO i=1,nterms
      RJref(1,:) = RJref(1,:) - us(:,i)*vs(i,Jref)
    END DO
  ELSE
    ! If we didn't change the reference column of the residual matrix "R",
    ! update the column to account for the new components of the approximation.
    RJref(1,:) = RJref(1,:) - us(:,nterms)*vs(nterms,Jref) !RJref -= vs[-1][Jref : Jref + 1][None, :] * us[-1][:, None]
  END IF
END DO
IF(step_size > 0.d0)THEN
  ! WRITE(*,*)'Savings:   ',iblock,jblock,nterms,MIN(M,N)
  ! Return the left and right approximation matrices.
  ! The approximate is such that:
  ! M ~ U_ACA.dot(V_ACA)
  ALLOCATE(self%aca_BU_mats(isparse,dim)%M(M,nterms))
  DO i=1,nterms
    DO k=1,M
      self%aca_BU_mats(isparse,dim)%M(k,i)=vs(i,k)
    END DO
  END DO
  ALLOCATE(self%aca_BV_mats(isparse,dim)%M(nterms,N))
  DO i=1,nterms
    DO k=1,N
      self%aca_BV_mats(isparse,dim)%M(i,k)=us(k,i)
    END DO
  END DO
  ! DEALLOCATE(self%sparse_dense(isparse,dim)%M)
  size_out=INT(nterms,8)*INT(M+N,8)
END IF
DEALLOCATE(RIstar,RJstar,RIref,RJref)
DEALLOCATE(us,vs,prevIstar,prevJstar)
DEALLOCATE(tmp_block%ipts,tmp_block%ielem,tmp_block%icell,tmp_block%inv_map)
END SUBROUTINE aca_approx
!
SUBROUTINE compress_aca(isparse,ilevel,iblock,jblock,dim,tol,size_out)
INTEGER(4), INTENT(in) :: isparse,ilevel,iblock,jblock,dim
REAL(8), INTENT(in) :: tol
INTEGER(8), INTENT(out) :: size_out
INTEGER(4) :: i,j,io_unit
REAL(8) :: S_sum,frob_K,tol_loc
INTEGER(4) :: M,N,K,LDA,LDU,LDVT,LWORK,INFO,MIN_DIM
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: IWORK
REAL(8), ALLOCATABLE, DIMENSION(:) :: S,WORK
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: U,VT,Atmp,QU,QV,RU,RV
!
! Get QR factorization of U and V matrices
! WRITE(*,*)'QR U'
CALL get_qr(self%aca_BU_mats(isparse,dim)%M,QU,RU)
! WRITE(*,*)'QR VT'
ALLOCATE(Atmp(SIZE(self%aca_BV_mats(isparse,dim)%M,2),SIZE(self%aca_BV_mats(isparse,dim)%M,1)))
Atmp = TRANSPOSE(self%aca_BV_mats(isparse,dim)%M)
CALL get_qr(Atmp,QV,RV)
DEALLOCATE(Atmp)
M=SIZE(QU,DIM=1)
N=SIZE(QV,DIM=1)
K=SIZE(RU,DIM=1)
! WRITE(*,*)K
ALLOCATE(Atmp(K,K))
CALL DGEMM('N', 'T', K, K, K, 1.d0, RU, K, RV, K, 0.d0, Atmp, K)
!
LDA = K
LDU = K
MIN_DIM = K !MIN(M,N)
LDVT = MIN_DIM
ALLOCATE(S(MIN_DIM),U(LDU,MIN_DIM),VT(LDVT,K))
!---Get optimal workspace size
LWORK = -1
ALLOCATE(WORK(1),IWORK(8*MIN_DIM))
CALL DGESDD('S', K, K, Atmp, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
IF( INFO.NE.0 ) THEN
  CALL oft_abort("The algorithm computing SVD failed to converge.","tw_Lmat_MF_setup::compress_aca",__FILE__)
END IF
LWORK = INT(WORK(1))
DEALLOCATE(WORK)
!---Compute SVD of matrix block
ALLOCATE(WORK(LWORK))
CALL DGESDD('S', K, K, Atmp, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
IF( INFO.NE.0 ) THEN
  CALL oft_abort("The algorithm computing SVD failed to converge.","tw_Lmat_MF_setup::compress_aca",__FILE__)
END IF
!---Truncate at desired accuracy
tol_loc = tol!*MAX(1.d0,S(1))
frob_K=0.d0
DO i=MIN_DIM,1,-1
  frob_K = frob_K + S(i)**2
  IF(SQRT(frob_K)>tol_loc)EXIT
END DO
i=MAX(MIN(self%min_rank,MIN_DIM),i)
IF((i==MIN_DIM).AND.(MIN_DIM>self%min_rank))THEN
  WRITE(*,*)'SVD recompression failed',S(1),S(MIN_DIM),tol_loc
END IF
!i=MIN(i,MIN_DIM)
!---Replace U and V matrices
DEALLOCATE(self%aca_BU_mats(isparse,dim)%M,self%aca_BV_mats(isparse,dim)%M)
ALLOCATE(self%aca_BU_mats(isparse,dim)%M(self%levels(ilevel)%blocks(jblock)%np,i))
ALLOCATE(self%aca_BV_mats(isparse,dim)%M(i,self%levels(ilevel)%blocks(iblock)%nelems))
! U = UQ.dot(U[:, :r] * S[:r])
! V = VT[:r, :].dot(VQ.T)
DO j=1,i
  U(:,j) = U(:,j)*S(j)
END DO
CALL DGEMM('N', 'N', M, i, K, 1.d0, QU, M, U, K, 0.d0, self%aca_BU_mats(isparse,dim)%M, M)
!self%aca_BU_mats(isparse,dim)%M = MATMUL(QU,U(:,1:i))
CALL DGEMM('N', 'T', i, N, K, 1.d0, VT, K, QV, N, 0.d0, self%aca_BV_mats(isparse,dim)%M, i)
!self%aca_BV_mats(isparse,dim)%M = MATMUL(VT(1:i,:),TRANSPOSE(QV))
size_out=INT(i,8)*INT(M+N,8)
! WRITE(*,*)'Compressed:',iblock,jblock,i,MIN(M,N)
DEALLOCATE(WORK,IWORK,S,U,VT,Atmp,QU,QV,RU,RV)
END SUBROUTINE compress_aca
END SUBROUTINE tw_hodlr_Bcompute
!------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in] a Source field
!! @param[out] b Result of matrix product
!------------------------------------------------------------------------------
subroutine tw_hodlr_Lapply(self,a,b)
class(oft_tw_hodlr_op), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a
class(oft_vector), intent(inout) :: b
integer(4) :: i,j,k,iblock,jblock,level
real(8), allocatable :: row_tmp(:),col_tmp(:),int_tmp1(:),int_tmp2(:)
REAL(8), pointer :: avals(:),bvals(:)
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL b%set(0.d0)
CALL b%get_local(bvals)
!$omp parallel private(k,iblock,jblock,level,row_tmp,col_tmp,int_tmp1,int_tmp2)
ALLOCATE(col_tmp(self%tw_obj%np_active),row_tmp(self%tw_obj%np_active))
ALLOCATE(int_tmp1(self%tw_obj%np_active),int_tmp2(self%tw_obj%np_active))
int_tmp2=0.d0
!$omp do schedule(dynamic,1)
DO j=1,self%nsparse
  level = self%sparse_blocks(1,j)
  iblock = self%sparse_blocks(2,j)
  jblock = self%sparse_blocks(3,j)
  IF(ASSOCIATED(self%aca_dense(j)%M))THEN
    !$omp simd
    DO k=1,self%levels(level)%blocks(iblock)%nelems
      col_tmp(k)=avals(self%levels(level)%blocks(iblock)%ielem(k))
    END DO
    CALL dgemv('N',self%levels(level)%blocks(jblock)%nelems,self%levels(level)%blocks(iblock)%nelems,1.d0,self%aca_dense(j)%M, &
      self%levels(level)%blocks(jblock)%nelems,col_tmp,1,0.d0,row_tmp,1)
    !$omp simd
    DO k=1,self%levels(level)%blocks(jblock)%nelems
      int_tmp2(self%levels(level)%blocks(jblock)%ielem(k))=int_tmp2(self%levels(level)%blocks(jblock)%ielem(k))  + row_tmp(k)
    END DO
    !
    !$omp simd
    DO k=1,self%levels(level)%blocks(jblock)%nelems
      col_tmp(k)=avals(self%levels(level)%blocks(jblock)%ielem(k))
    END DO
    CALL dgemv('T',self%levels(level)%blocks(jblock)%nelems,self%levels(level)%blocks(iblock)%nelems,1.d0,self%aca_dense(j)%M, &
      self%levels(level)%blocks(jblock)%nelems,col_tmp,1,0.d0,row_tmp,1)
    !$omp simd
    DO k=1,self%levels(level)%blocks(iblock)%nelems
      int_tmp2(self%levels(level)%blocks(iblock)%ielem(k))=int_tmp2(self%levels(level)%blocks(iblock)%ielem(k)) + row_tmp(k)
    END DO
  ELSE
    !---Upper side
    !$omp simd
    DO k=1,self%levels(level)%blocks(iblock)%nelems
      col_tmp(k)=avals(self%levels(level)%blocks(iblock)%ielem(k))
    END DO
    k=SIZE(self%aca_V_mats(j)%M,DIM=1)
    CALL dgemv('N',k,self%levels(level)%blocks(iblock)%nelems,1.d0,self%aca_V_mats(j)%M, &
      k,col_tmp,1,0.d0,int_tmp1,1)
    CALL dgemv('N',self%levels(level)%blocks(jblock)%nelems,k,1.d0,self%aca_U_mats(j)%M, &
      self%levels(level)%blocks(jblock)%nelems,int_tmp1,1,0.d0,row_tmp,1)
    !$omp simd
    DO k=1,self%levels(level)%blocks(jblock)%nelems
      int_tmp2(self%levels(level)%blocks(jblock)%ielem(k))=int_tmp2(self%levels(level)%blocks(jblock)%ielem(k))+row_tmp(k)
    END DO
    !---Apply other side of diagonal
    !$omp simd
    DO k=1,self%levels(level)%blocks(jblock)%nelems
      col_tmp(k)=avals(self%levels(level)%blocks(jblock)%ielem(k))
    END DO
    k=SIZE(self%aca_V_mats(j)%M,DIM=1)
    CALL dgemv('T',self%levels(level)%blocks(jblock)%nelems,k,1.d0,self%aca_U_mats(j)%M, &
      self%levels(level)%blocks(jblock)%nelems,col_tmp,1,0.d0,int_tmp1,1)
    CALL dgemv('T',k,self%levels(level)%blocks(iblock)%nelems,1.d0,self%aca_V_mats(j)%M, &
      k,int_tmp1,1,0.d0,row_tmp,1)
    !$omp simd
    DO k=1,self%levels(level)%blocks(iblock)%nelems
      int_tmp2(self%levels(level)%blocks(iblock)%ielem(k))=int_tmp2(self%levels(level)%blocks(iblock)%ielem(k))+row_tmp(k)
    END DO
  END IF
END DO
!$omp end do nowait
!$omp do schedule(static,1)
DO j=1,self%ndense
  level = self%dense_blocks(1,j)
  iblock = self%dense_blocks(2,j)
  jblock = self%dense_blocks(3,j)
  !$omp simd
  DO k=1,self%levels(level)%blocks(iblock)%nelems
    col_tmp(k)=avals(self%levels(level)%blocks(iblock)%ielem(k))
  END DO
  CALL dgemv('N',self%levels(level)%blocks(jblock)%nelems,self%levels(level)%blocks(iblock)%nelems,1.d0,self%dense_mats(j)%M, &
    self%levels(level)%blocks(jblock)%nelems,col_tmp,1,0.d0,row_tmp,1)
  !$omp simd
  DO k=1,self%levels(level)%blocks(jblock)%nelems
    int_tmp2(self%levels(level)%blocks(jblock)%ielem(k))=int_tmp2(self%levels(level)%blocks(jblock)%ielem(k))+row_tmp(k)
  END DO
END DO
!$omp critical
bvals(1:self%tw_obj%np_active)=bvals(1:self%tw_obj%np_active)+int_tmp2
!$omp end critical
DEALLOCATE(col_tmp,row_tmp,int_tmp1,int_tmp2)
!$omp end parallel
IF(self%tw_obj%nholes+self%tw_obj%n_vcoils>0)THEN
  CALL dgemv('N',self%tw_obj%nelems,self%tw_obj%nholes+self%tw_obj%n_vcoils,1.d0,self%hole_Vcoil_mat%M, &
    self%tw_obj%nelems,avals(self%tw_obj%np_active+1:self%tw_obj%nelems),1,1.d0,bvals,1)
  avals(self%tw_obj%np_active+1:self%tw_obj%nelems)=0.d0 ! Prevent double diagonal contributions
  CALL dgemv('T',self%tw_obj%nelems,self%tw_obj%nholes+self%tw_obj%n_vcoils,1.d0,self%hole_Vcoil_mat%M, &
    self%tw_obj%nelems,avals,1,1.d0,bvals(self%tw_obj%np_active+1:self%tw_obj%nelems),1)
END IF
CALL b%restore_local(bvals)
DEALLOCATE(avals,bvals)
end subroutine tw_hodlr_Lapply
!------------------------------------------------------------------------------
!> Apply the matrix to a field.
!!
!! b = self * a
!!
!! @note This subroutine is a dummy routine used to specify the interface
!! of the member function and catch errors in uninitialized matrices.
!!
!! @param[in] a Source field
!! @param[out] b Result of matrix product
!------------------------------------------------------------------------------
subroutine tw_hodlr_Bapply(self,a,bx,by,bz)
class(oft_tw_hodlr_op), intent(inout) :: self
class(oft_vector), target, intent(inout) :: a
class(oft_vector), intent(inout) :: bx,by,bz
integer(4) :: i,j,jj,k,dim,iblock,jblock,level,flip
real(8), allocatable :: row_tmp(:),col_tmp(:),int_tmp1(:),int_tmp2(:,:)
REAL(8), pointer :: avals(:),bvals(:,:),btmp(:)
NULLIFY(avals,bvals)
CALL a%get_local(avals)
CALL bx%set(0.d0)
CALL by%set(0.d0)
CALL bz%set(0.d0)
ALLOCATE(bvals(bx%n,3))
btmp=>bvals(:,1)
CALL bx%get_local(btmp)
btmp=>bvals(:,2)
CALL by%get_local(btmp)
btmp=>bvals(:,3)
CALL bz%get_local(btmp)
!$omp parallel private(k,dim,jj,flip,iblock,jblock,level,row_tmp,col_tmp,int_tmp1,int_tmp2)
ALLOCATE(col_tmp(self%tw_obj%np_active),row_tmp(self%tw_obj%mesh%np))
ALLOCATE(int_tmp1(self%tw_obj%mesh%np),int_tmp2(self%tw_obj%mesh%np,3))
int_tmp2=0.d0
!$omp do schedule(dynamic,1)
DO j=1,self%nsparse
  level = self%sparse_blocks(1,j)
  iblock = self%sparse_blocks(2,j)
  jblock = self%sparse_blocks(3,j)
  DO flip=1,2
    jj=(j-1)*2+flip
    IF(flip==2)THEN
      jblock = self%sparse_blocks(2,j)
      iblock = self%sparse_blocks(3,j)
    END IF
    DO dim=1,3
      IF(ASSOCIATED(self%aca_B_dense(jj,dim)%M))THEN
        !$omp simd
        DO k=1,self%levels(level)%blocks(iblock)%nelems
          col_tmp(k)=avals(self%levels(level)%blocks(iblock)%ielem(k))
        END DO
        CALL dgemv('N',self%levels(level)%blocks(jblock)%np,self%levels(level)%blocks(iblock)%nelems,1.d0,self%aca_B_dense(jj,dim)%M, &
          self%levels(level)%blocks(jblock)%np,col_tmp,1,0.d0,row_tmp,1)
        !$omp simd
        DO k=1,self%levels(level)%blocks(jblock)%np
          int_tmp2(self%levels(level)%blocks(jblock)%ipts(k),dim)=int_tmp2(self%levels(level)%blocks(jblock)%ipts(k),dim)  + row_tmp(k)
        END DO
      ELSE
        !$omp simd
        DO k=1,self%levels(level)%blocks(iblock)%nelems
          col_tmp(k)=avals(self%levels(level)%blocks(iblock)%ielem(k))
        END DO
        k=SIZE(self%aca_BV_mats(jj,dim)%M,DIM=1)
        CALL dgemv('N',k,self%levels(level)%blocks(iblock)%nelems,1.d0,self%aca_BV_mats(jj,dim)%M, &
          k,col_tmp,1,0.d0,int_tmp1,1)
        CALL dgemv('N',self%levels(level)%blocks(jblock)%np,k,1.d0,self%aca_BU_mats(jj,dim)%M, &
          self%levels(level)%blocks(jblock)%np,int_tmp1,1,0.d0,row_tmp,1)
        !$omp simd
        DO k=1,self%levels(level)%blocks(jblock)%np
          int_tmp2(self%levels(level)%blocks(jblock)%ipts(k),dim)=int_tmp2(self%levels(level)%blocks(jblock)%ipts(k),dim)+row_tmp(k)
        END DO
      END IF
    END DO
  END DO
END DO
! !$omp end do nowait
!$omp do schedule(static,1)
DO j=1,self%ndense
  level = self%dense_blocks(1,j)
  iblock = self%dense_blocks(2,j)
  jblock = self%dense_blocks(3,j)
  !$omp simd
  DO k=1,self%levels(level)%blocks(iblock)%nelems
    col_tmp(k)=avals(self%levels(level)%blocks(iblock)%ielem(k))
  END DO
  DO dim=1,3
    CALL dgemv('N',self%levels(level)%blocks(jblock)%np,self%levels(level)%blocks(iblock)%nelems,1.d0,self%dense_B_mats(j,dim)%M, &
      self%levels(level)%blocks(jblock)%np,col_tmp,1,0.d0,row_tmp,1)
    !$omp simd
    DO k=1,self%levels(level)%blocks(jblock)%np
      int_tmp2(self%levels(level)%blocks(jblock)%ipts(k),dim)=int_tmp2(self%levels(level)%blocks(jblock)%ipts(k),dim)+row_tmp(k)
    END DO
  END DO
END DO
!$omp critical
DO dim=1,3
  bvals(1:self%tw_obj%mesh%np,dim)=bvals(1:self%tw_obj%mesh%np,dim)+int_tmp2(:,dim)
END DO
!$omp end critical
DEALLOCATE(col_tmp,row_tmp,int_tmp1,int_tmp2)
!$omp end parallel
IF((self%tw_obj%nholes+self%tw_obj%n_vcoils)>0)THEN
  DO dim=1,3
    CALL dgemv('N',self%tw_obj%mesh%np,self%tw_obj%nholes+self%tw_obj%n_vcoils,1.d0,self%hole_Vcoil_Bmat(:,:,dim), &
      self%tw_obj%mesh%np,avals(self%tw_obj%np_active+1:self%tw_obj%nelems),1,1.d0,bvals(:,dim),1)
  END DO
END IF
CALL bx%restore_local(bvals(:,1))
CALL by%restore_local(bvals(:,2))
CALL bz%restore_local(bvals(:,3))
DEALLOCATE(avals,bvals)
end subroutine tw_hodlr_Bapply
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
subroutine tw_hodlr_Lassemble(self,diag)
class(oft_tw_hodlr_op), intent(inout) :: self
class(oft_vector), optional, target, intent(inout) :: diag
integer(4) :: i,j,level,iblock
real(8), POINTER, DIMENSION(:) :: vals
IF(present(diag))THEN
  IF(associated(self%D))CALL self%D%delete
  CALL diag%new(self%D)
  NULLIFY(vals)
  CALL self%D%get_local(vals)
  !---Get hole diagonal values
  IF(self%tw_obj%nholes+self%tw_obj%n_vcoils>0)THEN
    DO i=1,self%tw_obj%nholes+self%tw_obj%n_vcoils
      vals(i+self%tw_obj%np_active)=self%hole_Vcoil_mat%M(i+self%tw_obj%np_active,i)
    END DO
  END IF
  DO i=1,self%ndense
    level=self%dense_blocks(1,i)
    iblock=self%dense_blocks(2,i)
    DO j=1,self%levels(level)%blocks(iblock)%nelems
      vals(self%levels(level)%blocks(iblock)%ielem(j))=self%dense_mats(i)%M(j,j)
    END DO
  END DO
  CALL self%D%restore_local(vals)
  DEALLOCATE(vals)
  CALL diag%add(0.d0,1.d0,self%D)
END IF
self%nr=self%tw_obj%nelems; self%nrg=self%tw_obj%nelems
self%nc=self%tw_obj%nelems; self%ncg=self%tw_obj%nelems
end subroutine tw_hodlr_Lassemble
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
SUBROUTINE tw_part_mesh(self,leaf_target,nlevels,levels)
TYPE(tw_type), INTENT(inout) :: self
INTEGER(4), INTENT(in) :: leaf_target
INTEGER(4), INTENT(out) :: nlevels
TYPE(oft_tw_level), POINTER, INTENT(out) :: levels(:)
INTEGER(4) :: i,j
REAL(8) :: xs
! INTEGER(4), ALLOCATABLE, DIMENSION(:) :: child_inds
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: active_pts

!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
TYPE :: tw_oct_tree
  INTEGER(4) :: npts = 0
  INTEGER(4) :: nchildren = 0
  INTEGER(4) :: depth = -1
  INTEGER(4) :: dind = -1
  INTEGER(4) :: pind = -1
  REAL(8) :: bounds(2,3) = 0.d0
  REAL(8) :: extent = 1.d99
  REAL(8) :: center(3) = 0.d0
  INTEGER(4), POINTER, CONTIGUOUS, DIMENSION(:) :: pts => NULL()
  TYPE(tw_oct_tree), POINTER, DIMENSION(:) :: children => NULL()
  TYPE(tw_oct_tree), POINTER :: parent => NULL()
END TYPE tw_oct_tree
!------------------------------------------------------------------------------
! CLASS tw_oct_tree_list
!------------------------------------------------------------------------------
!> Need docs
!------------------------------------------------------------------------------
TYPE :: tw_oct_tree_level
  INTEGER(4) :: nblocks = 0
  INTEGER(4) :: nc_dense = 0
  INTEGER(4) :: nc_sparse = 0
  INTEGER(4), CONTIGUOUS, POINTER, DIMENSION(:) :: kr_dense => NULL()
  INTEGER(4), CONTIGUOUS, POINTER, DIMENSION(:) :: lc_dense => NULL()
  INTEGER(4), CONTIGUOUS, POINTER, DIMENSION(:) :: kc_dense => NULL()
  INTEGER(4), CONTIGUOUS, POINTER, DIMENSION(:) :: lr_dense => NULL()
  INTEGER(4), CONTIGUOUS, POINTER, DIMENSION(:) :: lrm_dense => NULL()
  INTEGER(4), CONTIGUOUS, POINTER, DIMENSION(:) :: kr_sparse => NULL()
  INTEGER(4), CONTIGUOUS, POINTER, DIMENSION(:) :: lc_sparse => NULL()
  INTEGER(4), POINTER, DIMENSION(:,:) :: dense_mask => NULL()
  TYPE(tw_oct_tree), POINTER, DIMENSION(:) :: block => NULL()
END TYPE tw_oct_tree_level

TYPE(tw_oct_tree) :: mesh_tree

!---Initialize tree
mesh_tree%npts=self%mesh%np !self%np_active
ALLOCATE(mesh_tree%pts(mesh_tree%npts))
ALLOCATE(active_pts(3,mesh_tree%npts))
DO i=1,self%mesh%np
  ! IF(self%pmap(i)==0)CYCLE
  ! mesh_tree%pts(self%pmap(i))=self%pmap(i)
  ! active_pts(:,self%pmap(i))=self%mesh%r(:,i)
  mesh_tree%pts(i)=i
  active_pts(:,i)=self%mesh%r(:,i)
END DO
!---find bounding box around mesh
DO i=1,3
  mesh_tree%bounds(1,i)=MINVAL(active_pts(i,:))
  mesh_tree%bounds(2,i)=MAXVAL(active_pts(i,:))
  IF(mesh_tree%bounds(2,i)-mesh_tree%bounds(1,i)<1.d-3)THEN
  mesh_tree%bounds(1,i)=mesh_tree%bounds(1,i)-5.d-4
  mesh_tree%bounds(2,i)=mesh_tree%bounds(2,i)+5.d-4
  END IF
  xs = mesh_tree%bounds(2,i) - mesh_tree%bounds(1,i)
  mesh_tree%bounds(1,i) = mesh_tree%bounds(1,i) - xs*1.d-2
  mesh_tree%bounds(2,i) = mesh_tree%bounds(2,i) + xs*1.d-2
END DO
! mesh_tree%extents=mesh_tree%bounds
!mesh_tree%center=SUM(mesh_tree%bounds,DIM=1)/2.d0
mesh_tree%center=SUM(active_pts,DIM=2)/REAL(mesh_tree%npts,8)
mesh_tree%depth=1
!---Recursively subdivide
nlevels=1
IF(mesh_tree%npts>leaf_target)THEN
  CALL subdivide_leaf(mesh_tree,2)
ELSE
  nlevels=1
END IF
DEALLOCATE(active_pts)
ALLOCATE(levels(nlevels))
CALL count_levels(mesh_tree,1)
DO i=1,nlevels
  ALLOCATE(levels(i)%blocks(levels(i)%nblocks))
  levels(i)%nblocks=0
END DO
CALL fill_levels(mesh_tree,1,1)
CONTAINS
!
RECURSIVE SUBROUTINE count_levels(leaf,depth)
TYPE(tw_oct_tree), INTENT(inout) :: leaf
INTEGER(4), INTENT(in) :: depth
INTEGER(4) :: j
!---Ensure we are at the bottom of each tree
IF(leaf%nchildren>0)THEN
  DO j=1,leaf%nchildren
    CALL count_levels(leaf%children(j),depth+1)
  END DO
ELSE
  DO j=depth+1,nlevels
    levels(j)%nblocks=levels(j)%nblocks+1
  END DO
END IF
!---Add leaf if non-empty
IF(leaf%npts>0)levels(depth)%nblocks=levels(depth)%nblocks+1
END SUBROUTINE count_levels
!
RECURSIVE SUBROUTINE fill_levels(leaf,depth,parent)
TYPE(tw_oct_tree), INTENT(inout) :: leaf
INTEGER(4), INTENT(in) :: depth,parent
INTEGER(4) :: i,j,block_id
!---Add leaf if non-empty
IF(leaf%npts==0)RETURN
levels(depth)%nblocks=levels(depth)%nblocks+1
!
block_id=levels(depth)%nblocks
levels(depth)%blocks(block_id)%np=leaf%npts
ALLOCATE(levels(depth)%blocks(block_id)%ipts(leaf%npts))
levels(depth)%blocks(block_id)%ipts=leaf%pts
levels(depth)%blocks(block_id)%nelems=0
DO i=1,leaf%npts
  IF(self%pmap(leaf%pts(i))==0)CYCLE
  levels(depth)%blocks(block_id)%nelems=levels(depth)%blocks(block_id)%nelems+1
END DO
ALLOCATE(levels(depth)%blocks(block_id)%ielem(levels(depth)%blocks(block_id)%nelems))
levels(depth)%blocks(block_id)%nelems=0
DO i=1,leaf%npts
  IF(self%pmap(leaf%pts(i))==0)CYCLE
  levels(depth)%blocks(block_id)%nelems=levels(depth)%blocks(block_id)%nelems+1
  levels(depth)%blocks(block_id)%ielem(levels(depth)%blocks(block_id)%nelems)=self%pmap(leaf%pts(i))
END DO
! levels(depth)%blocks(block_id)%ielem=leaf%pts
levels(depth)%blocks(block_id)%center=leaf%center
levels(depth)%blocks(block_id)%extent=leaf%extent
levels(depth)%blocks(block_id)%parent=parent
!---Ensure we are at the bottom of each tree
IF(leaf%nchildren>0)THEN
  DO j=1,leaf%nchildren
    CALL fill_levels(leaf%children(j),depth+1,block_id)
  END DO
ELSE
  DO j=depth+1,nlevels
    levels(j)%nblocks=levels(j)%nblocks+1
    !
    block_id=levels(j)%nblocks
    levels(j)%blocks(block_id)%np=leaf%npts
    ALLOCATE(levels(j)%blocks(block_id)%ipts(leaf%npts))
    levels(j)%blocks(block_id)%ipts=leaf%pts
    levels(j)%blocks(block_id)%nelems=0
    DO i=1,leaf%npts
      IF(self%pmap(leaf%pts(i))==0)CYCLE
      levels(j)%blocks(block_id)%nelems=levels(j)%blocks(block_id)%nelems+1
    END DO
    ALLOCATE(levels(j)%blocks(block_id)%ielem(levels(j)%blocks(block_id)%nelems))
    levels(j)%blocks(block_id)%nelems=0
    DO i=1,leaf%npts
      IF(self%pmap(leaf%pts(i))==0)CYCLE
      levels(j)%blocks(block_id)%nelems=levels(j)%blocks(block_id)%nelems+1
      levels(j)%blocks(block_id)%ielem(levels(j)%blocks(block_id)%nelems)=self%pmap(leaf%pts(i))
    END DO
    ! ALLOCATE(levels(j)%blocks(block_id)%ielem(leaf%npts))
    ! levels(j)%blocks(block_id)%ielem=leaf%pts
    levels(j)%blocks(block_id)%center=leaf%center
    levels(j)%blocks(block_id)%extent=leaf%extent
    levels(j)%blocks(block_id)%parent=levels(j-1)%nblocks
  END DO
END IF
!---Destroy leaf on way back up
leaf%npts=0
DEALLOCATE(leaf%pts)
END SUBROUTINE fill_levels
!---Helper function
RECURSIVE SUBROUTINE subdivide_leaf(leaf,depth)
TYPE(tw_oct_tree), TARGET, INTENT(inout) :: leaf
INTEGER(4), INTENT(in) :: depth
INTEGER(4) :: j,k,kk,i,ii,npts,sep_max(1)
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: child_mark
REAL(8) :: xm,ym,zm,rleft(3),rdenom(3),locs(3),xs,extents(2,3),split_center(3)
REAL(8), ALLOCATABLE, DIMENSION(:) :: vals
REAL(8), PARAMETER :: pt_tol=0.d-6
TYPE(tw_oct_tree) :: child_tmp(2)
! WRITE(*,*)'Subdividing block',leaf%ncells,nblocks
IF(nlevels>20)CALL oft_abort("Maximum recursion", "subdivide_leaf", __FILE__)
nlevels=MAX(nlevels,depth)
ALLOCATE(vals(leaf%npts),child_mark(leaf%npts))
child_mark=0
locs=0.d0
DO j=1,3
  vals=active_pts(j,leaf%pts)
  CALL sort_array(vals,child_mark,leaf%npts)
  split_center(j)=vals(INT(leaf%npts/2,4))
  locs(j)=SUM((vals-split_center(j))**2)
END DO
DEALLOCATE(vals,child_mark)
!---Split cell in half in each direction
xm = split_center(1) !(leaf%bounds(1,1)+leaf%bounds(2,1))/2.d0
ym = split_center(2) !(leaf%bounds(1,2)+leaf%bounds(2,2))/2.d0
zm = split_center(3) !(leaf%bounds(1,3)+leaf%bounds(2,3))/2.d0
DO k=1,leaf%npts
  IF(ABS(active_pts(1,leaf%pts(k))-xm)<1.d-8)xm=xm+(leaf%bounds(2,1)-leaf%bounds(1,1))*1.d-4
  IF(ABS(active_pts(2,leaf%pts(k))-ym)<1.d-8)ym=ym+(leaf%bounds(2,2)-leaf%bounds(1,2))*1.d-4
  IF(ABS(active_pts(3,leaf%pts(k))-zm)<1.d-8)zm=zm+(leaf%bounds(2,3)-leaf%bounds(1,3))*1.d-4
END DO
!---Split cell in half
child_tmp(1)%bounds=leaf%bounds
child_tmp(2)%bounds=leaf%bounds
sep_max=MAXLOC(locs)
j=sep_max(1)
!WRITE(*,*)'Split',j,locs
locs=[xm,ym,zm]
child_tmp(1)%bounds(1,j)=locs(j)
child_tmp(2)%bounds(2,j)=locs(j)
! WRITE(*,*)0
! WRITE(*,*)leaf%bounds(:,1)
! WRITE(*,*)leaf%bounds(:,2)
! WRITE(*,*)leaf%bounds(:,3)
! DO j=1,8
!   WRITE(*,*)j
!   WRITE(*,*)child_tmp(j)%bounds(:,1)
!   WRITE(*,*)child_tmp(j)%bounds(:,2)
!   WRITE(*,*)child_tmp(j)%bounds(:,3)
! END DO
!---
ALLOCATE(child_mark(leaf%npts))
child_mark=0
leaf%nchildren=0
DO j=1,2
  child_tmp(j)%depth=depth
  child_tmp(j)%parent=>leaf
  ! child_mark=0
  rleft = child_tmp(j)%bounds(1,:)
  rdenom = child_tmp(j)%bounds(2,:)-child_tmp(j)%bounds(1,:)
  DO k=1,leaf%npts
    locs = (active_pts(:,leaf%pts(k))-rleft)/rdenom
    IF(ALL(locs>=-pt_tol).AND.ALL(locs<1.d0+pt_tol).AND.child_mark(k)==0)child_mark(k)=j
  END DO
  child_tmp(j)%npts=COUNT(child_mark==j)
  IF(child_tmp(j)%npts==0)CYCLE
  ! IF(child_tmp(j)%npts<200)WRITE(*,*)'Small leaf',child_tmp(j)%npts
  leaf%nchildren=leaf%nchildren+1
  ALLOCATE(child_tmp(j)%pts(child_tmp(j)%npts))
  extents(1,:)=1.d99
  extents(2,:)=-1.d99
  child_tmp(j)%npts=0
  child_tmp(j)%center=0.d0
  DO k=1,leaf%npts
    IF(child_mark(k)/=j)CYCLE
    child_tmp(j)%npts=child_tmp(j)%npts+1
    child_tmp(j)%pts(child_tmp(j)%npts)=leaf%pts(k)
    child_tmp(j)%center=child_tmp(j)%center+active_pts(:,leaf%pts(k))
    DO kk=1,3
      extents(1,kk)=MIN(extents(1,kk),active_pts(kk,leaf%pts(k)))
      extents(2,kk)=MAX(extents(2,kk),active_pts(kk,leaf%pts(k)))
    END DO
  END DO
  DO kk=1,3
    IF(extents(2,kk)-extents(1,kk)<1.d-8)THEN
      extents(1,kk)=extents(1,kk)-5.d-9
      extents(2,kk)=extents(2,kk)+5.d-9
    END IF
    xs = extents(2,kk) - extents(1,kk)
    extents(1,kk) = MAX(child_tmp(j)%bounds(1,kk),extents(1,kk) - xs*1.d-2)
    extents(2,kk) = MIN(child_tmp(j)%bounds(2,kk),extents(2,kk) + xs*1.d-2)
    !---Shrink bounds to match extent
    child_tmp(j)%bounds(1,kk) = MAX(child_tmp(j)%bounds(1,kk),extents(1,kk))
    child_tmp(j)%bounds(2,kk) = MIN(child_tmp(j)%bounds(2,kk),extents(2,kk))
  END DO
  child_tmp(j)%center=child_tmp(j)%center/REAL(child_tmp(j)%npts,8)
  child_tmp(j)%extent=0.d0
  DO k=1,child_tmp(j)%npts
    child_tmp(j)%extent=MAX(child_tmp(j)%extent,magnitude(active_pts(:,child_tmp(j)%pts(k))-child_tmp(j)%center))
  END DO
END DO
IF(ANY(child_mark==0))THEN
  WRITE(*,*)leaf%depth
  WRITE(*,*)leaf%bounds(1,:)
  WRITE(*,*)leaf%bounds(2,:)
  WRITE(*,*)"---"
  DO k=1,leaf%npts
    IF(child_mark(k)==0)THEN
      WRITE(*,*)active_pts(:,leaf%pts(k))
      DO j=1,2!8
        rleft = child_tmp(j)%bounds(1,:)
        rdenom = child_tmp(j)%bounds(2,:)-child_tmp(j)%bounds(1,:)
        locs = (active_pts(:,leaf%pts(k))-rleft)/rdenom
        WRITE(*,*)j,locs
      END DO
    END IF
  END DO
  CALL oft_abort("Lost element in partitioning","",__FILE__)
END IF
DEALLOCATE(child_mark)
!---Subdivide further
IF(leaf%nchildren==0)RETURN
! WRITE(*,*)'Parent',depth,leaf%nchildren,leaf%ncells,branch_target==0
ALLOCATE(leaf%children(leaf%nchildren))
leaf%nchildren=0
DO j=1,2
  IF(child_tmp(j)%npts>0)THEN
    leaf%nchildren=leaf%nchildren+1
    CALL clone_tree(child_tmp(j),leaf%children(leaf%nchildren))
    ! WRITE(*,*)'  Child block',child_tmp(j)%ncells,child_tmp(j)%ncells<=branch_target
    ! DO kk=1,3
    !   WRITE(*,*)'    ',child_tmp(j)%bounds(:,kk)
    ! END DO
    IF(leaf%children(leaf%nchildren)%npts>leaf_target)THEN
      CALL subdivide_leaf(leaf%children(leaf%nchildren),depth+1)
    ! ELSE
    !   nleaves=nleaves+1
    END IF
  END IF
END DO
! WRITE(*,*)'Exit',depth
END SUBROUTINE subdivide_leaf
SUBROUTINE clone_tree(in_tree,out_tree)
TYPE(tw_oct_tree), INTENT(inout) :: in_tree,out_tree
out_tree%dind=in_tree%dind
out_tree%pind=in_tree%pind
out_tree%depth=in_tree%depth
out_tree%bounds=in_tree%bounds
out_tree%extent=in_tree%extent
out_tree%center=in_tree%center
out_tree%npts=in_tree%npts
out_tree%nchildren=in_tree%nchildren
out_tree%pts=>in_tree%pts
out_tree%children=>in_tree%children
out_tree%parent=>in_tree%parent
END SUBROUTINE clone_tree
END SUBROUTINE tw_part_mesh
!
SUBROUTINE get_qr(Mat,Q,R)
REAL(8), INTENT(in) :: Mat(:,:)
REAL(8), ALLOCATABLE, INTENT(out) :: Q(:,:),R(:,:)
!
INTEGER(4) :: i,j,M,N,K,LWORK,INFO
REAL(8), ALLOCATABLE, DIMENSION(:) :: WORK,TAU
M = SIZE(Mat,DIM=1)
N = SIZE(Mat,DIM=2)
K = MIN(M,N)
ALLOCATE(Q(M,N),TAU(MIN(M,N)))
Q = Mat
LWORK = -1
ALLOCATE(WORK(1))
CALL dgeqrf( M, N, Q, M, TAU, WORK, LWORK, INFO )
IF(INFO/=0)CALL oft_abort("QR factorization failed","tw_Lmat_MF_setup::get_qr",__FILE__)
LWORK = INT(WORK(1))
DEALLOCATE(WORK)
ALLOCATE(WORK(LWORK))
CALL dgeqrf( M, N, Q, M, TAU, WORK, LWORK, INFO )
IF(INFO/=0)CALL oft_abort("QR factorization failed","tw_Lmat_MF_setup::get_qr",__FILE__)
DEALLOCATE(WORK)
ALLOCATE(R(K,K))
R=0.d0
DO i=1,M
  DO j=i,N
    R(i,j)=Q(i,j)
  END DO
END DO
ALLOCATE(WORK(1))
LWORK=-1
CALL DORGQR( M, N, N, Q, M, TAU, WORK, LWORK, INFO )
IF(INFO/=0)CALL oft_abort("QR factorization failed","tw_Lmat_MF_setup::get_qr",__FILE__)
LWORK = INT(WORK(1))
DEALLOCATE(WORK)
ALLOCATE(WORK(LWORK))
CALL DORGQR( M, N, N, Q, M, TAU, WORK, LWORK, INFO )
IF(INFO/=0)CALL oft_abort("QR factorization failed","tw_Lmat_MF_setup::get_qr",__FILE__)
DEALLOCATE(WORK,TAU)
END SUBROUTINE get_qr
!---------------------------------------------------------------------------
! SUBROUTINE: bjprecond_apply
!---------------------------------------------------------------------------
!> Precondition a linear system using a Block-Jacobi method
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!---------------------------------------------------------------------------
RECURSIVE SUBROUTINE bjprecond_apply(self,u,g)
CLASS(oft_tw_hodlr_bjpre), INTENT(inout) :: self
CLASS(oft_cvector), INTENT(inout) :: u,g
!---
INTEGER(4) :: i,j,k,l,n,info,level,iblock
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: imap
COMPLEX(8), POINTER, DIMENSION(:) :: utmp,gtmp,uloc,gloc
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%inverse_mats))THEN
  ALLOCATE(imap(self%mf_obj%tw_obj%nelems))
  self%max_block_size = 0
  ALLOCATE(self%inverse_mats(self%mf_obj%ndense+1))
  DO i=1,self%mf_obj%ndense
    level=self%mf_obj%dense_blocks(1,i)
    iblock=self%mf_obj%dense_blocks(2,i)
    n=self%mf_obj%levels(level)%blocks(iblock)%nelems
    self%max_block_size = MAX(self%max_block_size,n)
    ALLOCATE(self%inverse_mats(i)%M(n,n))
    self%inverse_mats(i)%M=self%alpha*self%mf_obj%dense_mats(i)%M
    imap=0
    DO j=1,n
      imap(self%mf_obj%levels(level)%blocks(iblock)%ielem(j))=j
    END DO
    DO j=1,n
      k=self%mf_obj%levels(level)%blocks(iblock)%ielem(j)
      DO l=self%Rmat%kr(k),self%Rmat%kr(k+1)-1
        IF(imap(self%Rmat%lc(l))>0)THEN
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) = &
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) + self%beta*self%Rmat%M(l)
        END IF
      END DO
    END DO
  END DO
  !---Get holes
  n=self%mf_obj%tw_obj%nholes+self%mf_obj%tw_obj%n_vcoils
  IF(n>0)THEN
    i = self%mf_obj%ndense+1
    self%max_block_size = MAX(self%max_block_size,n)
    ALLOCATE(self%inverse_mats(i)%M(n,n))
    imap=0
    DO j=1,n
      DO k=1,n
        self%inverse_mats(i)%M(j,k) = &
          self%alpha*self%mf_obj%hole_Vcoil_mat%M(j+self%mf_obj%tw_obj%np_active,k)
      END DO
      imap(j+self%mf_obj%tw_obj%np_active)=j
    END DO
    DO j=1,n
      k=j+self%mf_obj%tw_obj%np_active
      DO l=self%Rmat%kr(k),self%Rmat%kr(k+1)-1
        IF(imap(self%Rmat%lc(l))>0)THEN
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) = &
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) + self%beta*self%Rmat%M(l)
        END IF
      END DO
    END DO
  END IF
  DEALLOCATE(imap)
  DO i=1,self%mf_obj%ndense+1
    level=self%mf_obj%dense_blocks(1,i)
    iblock=self%mf_obj%dense_blocks(2,i)
    IF(i==self%mf_obj%ndense+1)THEN
      IF(.NOT.ASSOCIATED(self%inverse_mats(i)%M))EXIT
      n = self%mf_obj%tw_obj%nholes+self%mf_obj%tw_obj%n_vcoils
    ELSE
      n = self%mf_obj%levels(level)%blocks(iblock)%nelems
    END IF
    ! oft_env%pm=.TRUE.
    CALL lapack_matinv(n,self%inverse_mats(i)%M,info)
    ! CALL lapack_cholesky(n,self%inverse_mats(i)%M,info)
    IF(info/=0)CALL oft_abort("factorization failed","",__FILE__)
    ! oft_env%pm=.FALSE.
  END DO
END IF
!
NULLIFY(utmp,gtmp)
CALL u%get_local(utmp)
CALL g%get_local(gtmp)
utmp=(0.d0,0.d0)
!$omp parallel private(uloc,gloc,i,level,iblock,n,j)
ALLOCATE(uloc(self%max_block_size),gloc(self%max_block_size))
!$omp do schedule(static,1)
DO i=1,self%mf_obj%ndense
  level=self%mf_obj%dense_blocks(1,i)
  iblock=self%mf_obj%dense_blocks(2,i)
  n = self%mf_obj%levels(level)%blocks(iblock)%nelems
  DO j=1,n
    gloc(j) = gtmp(self%mf_obj%levels(level)%blocks(iblock)%ielem(j))
  END DO
  uloc=0.d0
  CALL zgemv('N',n,n,(1.d0,0.d0),self%inverse_mats(i)%M,n,gloc,1,(0.d0,0.d0),uloc,1)
  DO j=1,n
    utmp(self%mf_obj%levels(level)%blocks(iblock)%ielem(j)) = uloc(j)
  END DO
  ! WRITE(*,*)i,MINVAL(uloc(1:n)),MAXVAL(uloc(1:n))
END DO
!$omp end do nowait
!$omp single
n = self%mf_obj%tw_obj%nholes+self%mf_obj%tw_obj%n_vcoils
IF(n > 0)THEN
  DO j=1,n
    gloc(j) = gtmp(self%mf_obj%tw_obj%np_active+j)
  END DO
  uloc=0.d0
  i = self%mf_obj%ndense+1
  CALL zgemv('N',n,n,(1.d0,0.d0),self%inverse_mats(i)%M,n,gloc,1,(0.d0,0.d0),uloc,1)
  DO j=1,n
    utmp(self%mf_obj%tw_obj%np_active+j) = uloc(j)
  END DO
END IF
!$omp end single
DEALLOCATE(uloc,gloc)
!$omp end parallel
CALL u%restore_local(utmp)
CALL g%set((0.d0,0.d0))
DEALLOCATE(utmp,gtmp)
DEBUG_STACK_POP
END SUBROUTINE bjprecond_apply
!---------------------------------------------------------------------------
!> Destroy Block-Jacobi preconditioner and deallocate all internal storage
!---------------------------------------------------------------------------
subroutine bjprecond_delete(self)
class(oft_tw_hodlr_bjpre), intent(inout) :: self
INTEGER(i4) :: i
DEBUG_STACK_PUSH
!---Destroy local solvers
IF(ASSOCIATED(self%inverse_mats))THEN
  DO i=1,self%mf_obj%nblocks+1
    IF(ASSOCIATED(self%inverse_mats(i)%M))DEALLOCATE(self%inverse_mats(i)%M)
  END DO
  DEALLOCATE(self%inverse_mats)
END IF
!---Reset
self%max_block_size=0
NULLIFY(self%mf_obj)
DEBUG_STACK_POP
end subroutine bjprecond_delete
!---------------------------------------------------------------------------
! SUBROUTINE: bjprecond_apply
!---------------------------------------------------------------------------
!> Precondition a linear system using a Block-Jacobi method
!!
!! @param[in,out] u Guess/Solution field
!! @param[in,out] g RHS/Residual field
!---------------------------------------------------------------------------
RECURSIVE SUBROUTINE rbjprecond_apply(self,u,g)
CLASS(oft_tw_hodlr_rbjpre), INTENT(inout) :: self
CLASS(oft_vector), INTENT(inout) :: u,g
!---
INTEGER(4) :: i,j,k,l,n,info,level,iblock,jblock
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: imap
REAL(8), POINTER, DIMENSION(:) :: utmp,gtmp,uloc,gloc
DEBUG_STACK_PUSH
IF(.NOT.ASSOCIATED(self%inverse_mats))THEN
  ALLOCATE(imap(self%mf_obj%tw_obj%nelems))
  self%max_block_size = 0
  ALLOCATE(self%inverse_mats(self%mf_obj%ndense+1))
  DO i=1,self%mf_obj%ndense
    level=self%mf_obj%dense_blocks(1,i)
    iblock=self%mf_obj%dense_blocks(2,i)
    n=self%mf_obj%levels(level)%blocks(iblock)%nelems
    self%max_block_size = MAX(self%max_block_size,n)
    ALLOCATE(self%inverse_mats(i)%M(n,n))
    self%inverse_mats(i)%M=self%alpha*self%mf_obj%dense_mats(i)%M
    imap=0
    DO j=1,n
      imap(self%mf_obj%levels(level)%blocks(iblock)%ielem(j))=j
    END DO
    DO j=1,n
      k=self%mf_obj%levels(level)%blocks(iblock)%ielem(j)
      DO l=self%Rmat%kr(k),self%Rmat%kr(k+1)-1
        IF(imap(self%Rmat%lc(l))>0)THEN
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) = &
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) + self%beta*self%Rmat%M(l)
        END IF
      END DO
    END DO
  END DO
  !---Get holes
  n=self%mf_obj%tw_obj%nholes+self%mf_obj%tw_obj%n_vcoils
  IF(n>0)THEN
    i = self%mf_obj%ndense+1
    self%max_block_size = MAX(self%max_block_size,n)
    ALLOCATE(self%inverse_mats(i)%M(n,n))
    imap=0
    DO j=1,n
      DO k=1,n
        self%inverse_mats(i)%M(j,k) = &
          self%alpha*self%mf_obj%hole_Vcoil_mat%M(j+self%mf_obj%tw_obj%np_active,k)
      END DO
      imap(j+self%mf_obj%tw_obj%np_active)=j
    END DO
    DO j=1,n
      k=j+self%mf_obj%tw_obj%np_active
      DO l=self%Rmat%kr(k),self%Rmat%kr(k+1)-1
        IF(imap(self%Rmat%lc(l))>0)THEN
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) = &
          self%inverse_mats(i)%M(j,imap(self%Rmat%lc(l))) + self%beta*self%Rmat%M(l)
        END IF
      END DO
    END DO
  END IF
  DEALLOCATE(imap)
  DO i=1,self%mf_obj%ndense+1
    level=self%mf_obj%dense_blocks(1,i)
    iblock=self%mf_obj%dense_blocks(2,i)
    IF(i==self%mf_obj%ndense+1)THEN
      IF(.NOT.ASSOCIATED(self%inverse_mats(i)%M))EXIT
      n = self%mf_obj%tw_obj%nholes+self%mf_obj%tw_obj%n_vcoils
    ELSE
      n = self%mf_obj%levels(level)%blocks(iblock)%nelems
    END IF
    ! oft_env%pm=.TRUE.
    ! CALL lapack_matinv(n,self%inverse_mats(i)%M,info)
    CALL lapack_cholesky(n,self%inverse_mats(i)%M,info)
    IF(info/=0)CALL oft_abort("factorization failed","",__FILE__)
    ! oft_env%pm=.FALSE.
  END DO
END IF
!
NULLIFY(utmp,gtmp)
CALL u%get_local(utmp)
CALL g%get_local(gtmp)
utmp=0.d0
!$omp parallel private(uloc,gloc,i,level,iblock,n,j)
ALLOCATE(uloc(self%max_block_size),gloc(self%max_block_size))
!$omp do schedule(static,1)
DO i=1,self%mf_obj%ndense
  level=self%mf_obj%dense_blocks(1,i)
  iblock=self%mf_obj%dense_blocks(2,i)
  n = self%mf_obj%levels(level)%blocks(iblock)%nelems
  DO j=1,n
    gloc(j) = gtmp(self%mf_obj%levels(level)%blocks(iblock)%ielem(j))
  END DO
  uloc=0.d0
  CALL dgemv('N',n,n,1.d0,self%inverse_mats(i)%M,n,gloc,1,0.d0,uloc,1)
  DO j=1,n
    utmp(self%mf_obj%levels(level)%blocks(iblock)%ielem(j)) = uloc(j)
  END DO
  ! WRITE(*,*)i,MINVAL(uloc(1:n)),MAXVAL(uloc(1:n))
END DO
!$omp end do nowait
!$omp single
n = self%mf_obj%tw_obj%nholes+self%mf_obj%tw_obj%n_vcoils
IF(n > 0)THEN
  DO j=1,n
    gloc(j) = gtmp(self%mf_obj%tw_obj%np_active+j)
  END DO
  uloc=0.d0
  i = self%mf_obj%ndense+1
  CALL dgemv('N',n,n,1.d0,self%inverse_mats(i)%M,n,gloc,1,0.d0,uloc,1)
  DO j=1,n
    utmp(self%mf_obj%tw_obj%np_active+j) = uloc(j)
  END DO
END IF
!$omp end single
DEALLOCATE(uloc,gloc)
!$omp end parallel
CALL u%restore_local(utmp)
CALL g%set(0.d0)
DEALLOCATE(utmp,gtmp)
DEBUG_STACK_POP
END SUBROUTINE rbjprecond_apply
!---------------------------------------------------------------------------
!> Destroy Block-Jacobi preconditioner and deallocate all internal storage
!---------------------------------------------------------------------------
subroutine rbjprecond_delete(self)
class(oft_tw_hodlr_rbjpre), intent(inout) :: self
INTEGER(i4) :: i
DEBUG_STACK_PUSH
!---Destroy local solvers
IF(ASSOCIATED(self%inverse_mats))THEN
  DO i=1,self%mf_obj%nblocks+1
    IF(ASSOCIATED(self%inverse_mats(i)%M))DEALLOCATE(self%inverse_mats(i)%M)
  END DO
  DEALLOCATE(self%inverse_mats)
END IF
!---Reset
self%max_block_size=0
NULLIFY(self%mf_obj)
DEBUG_STACK_POP
end subroutine rbjprecond_delete
END MODULE thin_wall_hodlr
