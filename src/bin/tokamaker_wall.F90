!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!> Experiment specific components.
!!
!! @authors Chris Hansen
!! @date March 2015
!! @ingroup doxy_tokamaker
!------------------------------------------------------------------------------
MODULE exp_geom
USE oft_base
USE oft_gs, ONLY: gs_eq
IMPLICIT NONE
PRIVATE
INTEGER(4) :: ltx_vv_reg = 6
REAL(8) :: ltx_vv_amp = 2.d1
PUBLIC exp_setup
CONTAINS
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE exp_setup(gs)
CLASS(gs_eq), INTENT(inout) :: gs
INTEGER(4) :: ierr
CALL ltx_setup(gs,ierr)
END SUBROUTINE exp_setup
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
SUBROUTINE ltx_setup(gs,ierr)
CLASS(gs_eq), INTENT(inout) :: gs
INTEGER(4), INTENT(out) :: ierr
INTEGER(4) :: io_unit
NAMELIST/ltx_options/ltx_vv_reg,ltx_vv_amp
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,ltx_options,IOSTAT=ierr)
IF(ierr>0)CALL oft_abort('Error parsing "ltx_options" in input file.','ltx_setup',__FILE__)
CLOSE(io_unit)
IF(ierr<0)RETURN
WRITE(*,*)
WRITE(*,'(A)')'*** Setting options for the Lithium Tokamak eXperiment ***'
gs%set_eta=>ltx_eta_set
END SUBROUTINE ltx_setup
!---------------------------------------------------------------------------------
!> Needs docs
!---------------------------------------------------------------------------------
FUNCTION ltx_eta_set(rc,id) RESULT(eta)
real(8), intent(in) :: rc(2)
integer(4), intent(in) :: id
real(8) :: eta
eta=1.d0
!---Vacuum vessel effective resistivity
IF(id==ltx_vv_reg)THEN
  IF((rc(1)>.21d0.AND.rc(1)<0.59d0).OR.ABS(rc(2))<0.36d0)eta=ltx_vv_amp
END IF
END FUNCTION ltx_eta_set
END MODULE exp_geom  
!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> Handle non-axisymmetric walls
!!
!! @authors Chris Hansen
!! @date July 2017
!! @ingroup doxy_tokamaker
!------------------------------------------------------------------------------
MODULE nonax_wall
USE oft_base
USE axi_green, ONLY: green
IMPLICIT NONE
CONTAINS
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
function biot_savart_elem(rc,thetas,ri_xyz) result(bout)
real(8), intent(in) :: rc(2,2)
real(8), intent(in) :: thetas(2)
real(8), intent(in) :: ri_xyz(3)
real(8) :: bout(3),rc_xyz(3,2),dl(3),rsep(3)
real(8) :: rmag,theta,theta_p
bout=0.d0
!---Forward toroidal
theta_p=thetas(1)
theta=thetas(2)
rc_xyz(:,1)=(/COS(theta_p)*rc(1,1),SIN(theta_p)*rc(1,1),rc(2,1)/)
rc_xyz(:,2)=(/COS(theta)*rc(1,1),SIN(theta)*rc(1,1),rc(2,1)/)
rsep = ri_xyz - SUM(rc_xyz,DIM=2)/2.d0
dl = rc_xyz(:,2) - rc_xyz(:,1)
rmag = SQRT(SUM(rsep**2))
bout = bout + cross_product(dl,rsep)/(rmag**3)
!---Middle step
theta=thetas(2)
rc_xyz(:,1)=(/COS(theta)*rc(1,1),SIN(theta)*rc(1,1),rc(2,1)/)
rc_xyz(:,2)=(/COS(theta)*rc(1,2),SIN(theta)*rc(1,2),rc(2,2)/)
rsep = ri_xyz - SUM(rc_xyz,DIM=2)/2.d0
dl = rc_xyz(:,2) - rc_xyz(:,1)
rmag = SQRT(SUM(rsep**2))
bout = bout + cross_product(dl,rsep)/(rmag**3)
!---Backward toroidal
theta_p=thetas(2)
theta=thetas(1)
rc_xyz(:,1)=(/COS(theta_p)*rc(1,2),SIN(theta_p)*rc(1,2),rc(2,2)/)
rc_xyz(:,2)=(/COS(theta)*rc(1,2),SIN(theta)*rc(1,2),rc(2,2)/)
rsep = ri_xyz - SUM(rc_xyz,DIM=2)/2.d0
dl = rc_xyz(:,2) - rc_xyz(:,1)
rmag = SQRT(SUM(rsep**2))
bout = bout + cross_product(dl,rsep)/(rmag**3)
!---Home step
theta=thetas(1)
rc_xyz(:,1)=(/COS(theta)*rc(1,2),SIN(theta)*rc(1,2),rc(2,2)/)
rc_xyz(:,2)=(/COS(theta)*rc(1,1),SIN(theta)*rc(1,1),rc(2,1)/)
rsep = ri_xyz - SUM(rc_xyz,DIM=2)/2.d0
dl = rc_xyz(:,2) - rc_xyz(:,1)
rmag = SQRT(SUM(rsep**2))
bout = bout + cross_product(dl,rsep)/(rmag**3)
!
bout=bout/(4.d0*pi)
end function biot_savart_elem
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE nonax_rescouple(ncoils,rc,extent,result,drivers,ndriver,correction, &
                            sensors,nsensor,inc_mirror,eta,grid_rz,grid_correction)
INTEGER(4), INTENT(in) :: ncoils
REAL(8), INTENT(in) :: rc(2,ncoils),extent(2)
REAL(8), INTENT(out) :: result(ncoils),correction(nsensor)
REAL(8), INTENT(in) :: drivers(2,ndriver),sensors(6,nsensor)
INTEGER(4), INTENT(in) :: ndriver,nsensor
LOGICAL, INTENT(in) :: inc_mirror
REAL(8), OPTIONAL, INTENT(in) :: eta(ncoils)
REAL(8), OPTIONAL, INTENT(in) :: grid_rz(:,:)
REAL(8), OPTIONAL, INTENT(inout) :: grid_correction(:,:,:)
!---
INTEGER(4) :: i,j,k,l,info,N,LWORK,nphi_grid
integer(4), allocatable :: IPIV(:)
INTEGER(4), PARAMETER :: ntheta=3,nphi=40
REAL(8) :: phi,phi0,phi1,dx(3),dy(3),btmp(3)
REAL(8) :: relem_row(3,4),thetas(2),rcen(3),dym,dxm
REAL(8), ALLOCATABLE, DIMENSION(:) :: etatmp,WORK,vert_couple
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: Mmat,dltmp,corr_mat
WRITE(*,*)'Computing non-axisymmetric coupling: resistive limit'
!--- Create coupling matrix
phi0=extent(1)
phi1=extent(2)
ALLOCATE(Mmat((ncoils-1)*nphi,(ncoils-1)*nphi))
ALLOCATE(etatmp(ncoils),vert_couple((ncoils-1)*nphi))
ALLOCATE(dltmp(4,(ncoils-1)*nphi))
etatmp=1.d0
IF(PRESENT(eta))etatmp=eta
IF(PRESENT(grid_correction))THEN
  IF(.NOT.PRESENT(grid_rz))CALL oft_abort('Grid correction requires "grid_rz"','nonax_rescouple',__FILE__)
  nphi_grid=SIZE(grid_correction,DIM=2)
ELSE
  nphi_grid=0
END IF
Mmat=0.d0
dltmp=0.d0
vert_couple=0.d0
!
ALLOCATE(corr_mat(nsensor,(ncoils-1)*nphi))
corr_mat=0.d0
!$omp parallel do private(j,phi,relem_row,dx,dy,dxm,dym)
DO i=1,ncoils-1
  DO j=1,nphi
    IF(j==1)THEN
      phi = (j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
      relem_row(:,1)=(/COS(phi)*rc(1,i),SIN(phi)*rc(1,i),rc(2,i)/)
      phi = (j)*(phi1-phi0)/REAL(nphi,8)+phi0
      relem_row(:,2)=(/COS(phi)*rc(1,i),SIN(phi)*rc(1,i),rc(2,i)/)
      relem_row(:,3)=(/COS(phi)*rc(1,i+1),SIN(phi)*rc(1,i+1),rc(2,i+1)/)
      phi = (j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
      relem_row(:,4)=(/COS(phi)*rc(1,i+1),SIN(phi)*rc(1,i+1),rc(2,i+1)/)
      dx = (relem_row(:,2)+relem_row(:,3))/2.d0-(relem_row(:,1)+relem_row(:,4))/2.d0
      dy = (relem_row(:,1)+relem_row(:,2))/2.d0-(relem_row(:,3)+relem_row(:,4))/2.d0
      dxm = magnitude(dx)/2.d0
      dym = magnitude(dy)/2.d0
    END IF
    dltmp(1,(i-1)*nphi+j)=dltmp(1,(i-1)*nphi+j)+dym/etatmp(i)
    dltmp(2,(i-1)*nphi+j)=dltmp(2,(i-1)*nphi+j)+dxm
    dltmp(3,(i-1)*nphi+j)=dltmp(3,(i-1)*nphi+j)+dym/etatmp(i+1)
    dltmp(4,(i-1)*nphi+j)=dltmp(4,(i-1)*nphi+j)+dxm
    IF(i>1)dltmp(3,(i-2)*nphi+j)=dltmp(3,(i-2)*nphi+j)+dym/etatmp(i)
    IF(i<ncoils-1)dltmp(1,(i)*nphi+j)=dltmp(1,(i)*nphi+j)+dym/etatmp(i+1)
    IF(j>1)dltmp(2,(i-1)*nphi+j-1)=dltmp(2,(i-1)*nphi+j-1)+dxm
    IF(j<nphi)dltmp(4,(i-1)*nphi+j+1)=dltmp(4,(i-1)*nphi+j+1)+dxm
  END DO
END DO
!$omp parallel do private(j,phi,relem_row,thetas,k,rcen,btmp)
DO i=1,ncoils-1
  DO j=1,nphi
    phi = (j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
    relem_row(:,1)=(/COS(phi)*rc(1,i),SIN(phi)*rc(1,i),rc(2,i)/)
    phi = (j)*(phi1-phi0)/REAL(nphi,8)+phi0
    relem_row(:,2)=(/COS(phi)*rc(1,i),SIN(phi)*rc(1,i),rc(2,i)/)
    relem_row(:,3)=(/COS(phi)*rc(1,i+1),SIN(phi)*rc(1,i+1),rc(2,i+1)/)
    phi = (j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
    relem_row(:,4)=(/COS(phi)*rc(1,i+1),SIN(phi)*rc(1,i+1),rc(2,i+1)/)
    !
    DO k=1,ndriver
      vert_couple((i-1)*nphi+j) = vert_couple((i-1)*nphi+j) &
          + green(rc(1,i),rc(2,i),drivers(1,k),drivers(2,k))*rc(1,i)*(phi1-phi0)/REAL(nphi,8)
      vert_couple((i-1)*nphi+j) = vert_couple((i-1)*nphi+j) &
          - green(rc(1,i+1),rc(2,i+1),drivers(1,k),drivers(2,k))*rc(1,i+1)*(phi1-phi0)/REAL(nphi,8)
    END DO
    !---
    Mmat((i-1)*nphi+j,(i-1)*nphi+j)= magnitude(relem_row(:,2)-relem_row(:,1))/dltmp(1,(i-1)*nphi+j) &
                                    + magnitude(relem_row(:,3)-relem_row(:,2))/dltmp(2,(i-1)*nphi+j) &
                                    + magnitude(relem_row(:,4)-relem_row(:,3))/dltmp(3,(i-1)*nphi+j) &
                                    + magnitude(relem_row(:,1)-relem_row(:,4))/dltmp(4,(i-1)*nphi+j)
    IF(i>1)THEN
      Mmat((i-2)*nphi+j,(i-1)*nphi+j)=Mmat((i-2)*nphi+j,(i-1)*nphi+j) &
      - magnitude(relem_row(:,2)-relem_row(:,1))/dltmp(1,(i-1)*nphi+j)
    END IF
    IF(i<ncoils-1)THEN
      Mmat((i)*nphi+j,(i-1)*nphi+j)=Mmat((i)*nphi+j,(i-1)*nphi+j) &
      - magnitude(relem_row(:,4)-relem_row(:,3))/dltmp(3,(i-1)*nphi+j)
    END IF
    !
    IF(j>1)THEN
      Mmat((i-1)*nphi+j-1,(i-1)*nphi+j)=Mmat((i-1)*nphi+j-1,(i-1)*nphi+j) &
      - magnitude(relem_row(:,3)-relem_row(:,2))/dltmp(2,(i-1)*nphi+j)
    END IF
    IF(j<nphi)THEN
      Mmat((i-1)*nphi+j+1,(i-1)*nphi+j)=Mmat((i-1)*nphi+j+1,(i-1)*nphi+j) &
      - magnitude(relem_row(:,1)-relem_row(:,4))/dltmp(4,(i-1)*nphi+j)
    END IF
    thetas(1)=(j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
    thetas(2)=(j)*(phi1-phi0)/REAL(nphi,8)+phi0
    DO k=1,nsensor
      phi=sensors(2,k)
      rcen=(/COS(phi)*sensors(1,k),SIN(phi)*sensors(1,k),sensors(3,k)/)
      btmp=biot_savart_elem(rc(:,i:i+1),thetas,rcen)
      corr_mat(k,(i-1)*nphi+j) = sensors(4,k)*(COS(phi)*btmp(1) + SIN(phi)*btmp(2)) &
                                + sensors(5,k)*(-SIN(phi)*btmp(1) + COS(phi)*btmp(2)) &
                                + sensors(6,k)*btmp(3)
      !
      IF(inc_mirror)THEN
        phi=phi+pi !phi1+phi0
        rcen=(/COS(phi)*sensors(1,k),SIN(phi)*sensors(1,k),sensors(3,k)/)
        btmp=biot_savart_elem(rc(:,i:i+1),thetas,rcen)
        corr_mat(k,(i-1)*nphi+j) = corr_mat(k,(i-1)*nphi+j) &
                                  + sensors(4,k)*(COS(phi)*btmp(1) + SIN(phi)*btmp(2)) &
                                  + sensors(5,k)*(-SIN(phi)*btmp(1) + COS(phi)*btmp(2)) &
                                  + sensors(6,k)*btmp(3)
      END IF
    END DO
  END DO
END DO
DEALLOCATE(dltmp)
!---Invert Mmat
N = (ncoils-1)*nphi
LWORK = MAX(1,2*N)
ALLOCATE(IPIV(N))
CALL dgetrf(N,N,Mmat,N,IPIV,INFO)
IF(INFO/=0)WRITE(*,*)INFO
LWORK=-1
ALLOCATE(WORK(1))
CALL dgetri(N,Mmat,N,IPIV,WORK,LWORK,INFO)
LWORK=INT(WORK(1),4)
DEALLOCATE(WORK)
ALLOCATE(WORK(LWORK))
CALL dgetri(N,Mmat,N,IPIV,WORK,LWORK,INFO)
IF(INFO/=0)WRITE(*,*)INFO
DEALLOCATE(IPIV,WORK)
vert_couple=MATMUL(Mmat,vert_couple)
DEALLOCATE(Mmat)
!---Toroidally average
result=0.d0
DO j=1,nphi
  result(1)=result(1)+vert_couple(j)/REAL(nphi,8)
  result(ncoils)=result(ncoils)-vert_couple((ncoils-2)*nphi+j)/REAL(nphi,8)
  DO i=2,ncoils-1
    result(i)=result(i)+(vert_couple((i-1)*nphi+j)-vert_couple((i-2)*nphi+j))/REAL(nphi,8)
  END DO
END DO
dxm=SQRT(DOT_PRODUCT(result,result))
result=result/dxm
vert_couple=vert_couple/dxm
correction=MATMUL(corr_mat,vert_couple)
IF(nphi_grid>0)THEN
  WRITE(*,*)'  Computing grid correction'
  grid_correction=0.d0
  !$omp parallel do private(j,phi,thetas,k,rcen,btmp,l)
  DO i=1,ncoils-1
    DO j=1,nphi
      thetas(1)=(j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
      thetas(2)=(j)*(phi1-phi0)/REAL(nphi,8)+phi0
      !
      DO k=1,SIZE(grid_rz,DIM=2)
        DO l=1,nphi_grid
          phi = (l-1)*2.d0*pi/REAL(nphi_grid,8)
          rcen = [COS(phi)*grid_rz(1,k),SIN(phi)*grid_rz(1,k),grid_rz(2,k)]
          btmp=biot_savart_elem(rc(:,i:i+1),thetas,rcen)
          grid_correction(:,l,k) = grid_correction(:,l,k) + &
            [COS(phi)*btmp(1)+SIN(phi)*btmp(2), COS(phi)*btmp(2)-SIN(phi)*btmp(1), btmp(3)]*vert_couple((i-1)*nphi+j)
          IF(inc_mirror)THEN
            phi=phi+pi
            rcen = [COS(phi)*grid_rz(1,k),SIN(phi)*grid_rz(1,k),grid_rz(2,k)]
            btmp=biot_savart_elem(rc(:,i:i+1),thetas,rcen)
            grid_correction(:,l,k) = grid_correction(:,l,k) + &
              + [COS(phi)*btmp(1)+SIN(phi)*btmp(2), COS(phi)*btmp(2)-SIN(phi)*btmp(1), btmp(3)]*vert_couple((i-1)*nphi+j)
          END IF
        END DO
      END DO
    END DO
  END DO
END IF
WRITE(*,*)'  Finished'
DEALLOCATE(vert_couple,etatmp,corr_mat)
END SUBROUTINE nonax_rescouple
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE nonax_indcouple(ncoils,rc,extent,result,drivers,ndriver,correction, &
                            sensors,nsensor,inc_mirror,grid_rz,grid_correction)
INTEGER(4), INTENT(in) :: ncoils
REAL(8), INTENT(in) :: rc(2,ncoils),extent(2)
REAL(8), INTENT(out) :: result(ncoils),correction(nsensor)
REAL(8), INTENT(in) :: drivers(2,ndriver),sensors(6,nsensor)
INTEGER(4), INTENT(in) :: ndriver,nsensor
LOGICAL, INTENT(in) :: inc_mirror
REAL(8), OPTIONAL, INTENT(in) :: grid_rz(:,:)
REAL(8), OPTIONAL, INTENT(inout) :: grid_correction(:,:,:)
!---
INTEGER(4) :: i,j,k,kk,l,info,N,LWORK,nphi_grid
integer(4), allocatable :: IPIV(:)
INTEGER(4), PARAMETER :: ntheta=3,nphi=40
REAL(8) :: norm(3),phi,phi0,phi1,dx(3),dy(3),btmp(3)
REAL(8) :: relem_row(3,4),thetas(2),rcen(3),dym,dxm
REAL(8), ALLOCATABLE, DIMENSION(:) :: WORK,vert_couple
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: Mmat,corr_mat
WRITE(*,*)'Computing non-axisymmetric coupling: inductive limit'
!--- Create coupling matrix
phi0=extent(1)
phi1=extent(2)
ALLOCATE(Mmat((ncoils-1)*nphi,(ncoils-1)*nphi))
ALLOCATE(vert_couple((ncoils-1)*nphi))
Mmat=0.d0
vert_couple=0.d0
IF(PRESENT(grid_correction))THEN
  IF(.NOT.PRESENT(grid_rz))CALL oft_abort('Grid correction requires "grid_rz"','nonax_indcouple',__FILE__)
  nphi_grid=SIZE(grid_correction,DIM=2)
ELSE
  nphi_grid=0
END IF
!
ALLOCATE(corr_mat(nsensor,(ncoils-1)*nphi))
corr_mat=0.d0
!$omp parallel do private(j,k,kk,phi,relem_row,rcen,dx,dy,norm,dxm,dym,thetas,btmp)
DO i=1,ncoils-1
  DO j=1,nphi
    phi = (j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
    relem_row(:,1)=(/COS(phi)*rc(1,i),SIN(phi)*rc(1,i),rc(2,i)/)
    phi = (j)*(phi1-phi0)/REAL(nphi,8)+phi0
    relem_row(:,2)=(/COS(phi)*rc(1,i),SIN(phi)*rc(1,i),rc(2,i)/)
    relem_row(:,3)=(/COS(phi)*rc(1,i+1),SIN(phi)*rc(1,i+1),rc(2,i+1)/)
    phi = (j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
    relem_row(:,4)=(/COS(phi)*rc(1,i+1),SIN(phi)*rc(1,i+1),rc(2,i+1)/)
    rcen=(relem_row(:,1)+relem_row(:,2)+relem_row(:,3)+relem_row(:,4))/4.d0
    dx = (relem_row(:,2)+relem_row(:,3))/2.d0-(relem_row(:,1)+relem_row(:,4))/2.d0
    dy = (relem_row(:,1)+relem_row(:,2))/2.d0-(relem_row(:,3)+relem_row(:,4))/2.d0
    norm = cross_product(dx,dy)
    dxm = magnitude(dx)
    dym = magnitude(dy)
    DO k=1,ncoils-1
      DO kk=1,nphi
        thetas(1)=(kk-1)*(phi1-phi0)/REAL(nphi,8)+phi0
        thetas(2)=(kk)*(phi1-phi0)/REAL(nphi,8)+phi0
        btmp=biot_savart_elem(rc(:,k:k+1),thetas,rcen)
        Mmat((i-1)*nphi+j,(k-1)*nphi+kk) = DOT_PRODUCT(btmp,norm)
      END DO
    END DO
    !
    DO k=1,ndriver
      vert_couple((i-1)*nphi+j) = vert_couple((i-1)*nphi+j) &
          + green(rc(1,i),rc(2,i),drivers(1,k),drivers(2,k))*rc(1,i)*(phi1-phi0)/REAL(nphi,8)
      vert_couple((i-1)*nphi+j) = vert_couple((i-1)*nphi+j) &
          - green(rc(1,i+1),rc(2,i+1),drivers(1,k),drivers(2,k))*rc(1,i+1)*(phi1-phi0)/REAL(nphi,8)
    END DO
    thetas(1)=(j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
    thetas(2)=(j)*(phi1-phi0)/REAL(nphi,8)+phi0
    DO k=1,nsensor
      phi=sensors(2,k)
      rcen=(/COS(phi)*sensors(1,k),SIN(phi)*sensors(1,k),sensors(3,k)/)
      btmp=biot_savart_elem(rc(:,i:i+1),thetas,rcen)
      corr_mat(k,(i-1)*nphi+j) = sensors(4,k)*(COS(phi)*btmp(1) + SIN(phi)*btmp(2)) &
                                + sensors(5,k)*(-SIN(phi)*btmp(1) + COS(phi)*btmp(2)) &
                                + sensors(6,k)*btmp(3)
      !
      IF(inc_mirror)THEN
        phi=phi+pi !phi1+phi0
        rcen=(/COS(phi)*sensors(1,k),SIN(phi)*sensors(1,k),sensors(3,k)/)
        btmp=biot_savart_elem(rc(:,i:i+1),thetas,rcen)
        corr_mat(k,(i-1)*nphi+j) = corr_mat(k,(i-1)*nphi+j) &
                                  + sensors(4,k)*(COS(phi)*btmp(1) + SIN(phi)*btmp(2)) &
                                  + sensors(5,k)*(-SIN(phi)*btmp(1) + COS(phi)*btmp(2)) &
                                  + sensors(6,k)*btmp(3)
      END IF
    END DO
  END DO
END DO
!---Invert Mmat
N = (ncoils-1)*nphi
LWORK = MAX(1,2*N)
ALLOCATE(IPIV(N))
CALL dgetrf(N,N,Mmat,N,IPIV,INFO)
IF(INFO/=0)WRITE(*,*)INFO
LWORK=-1
ALLOCATE(WORK(1))
CALL dgetri(N,Mmat,N,IPIV,WORK,LWORK,INFO)
LWORK=INT(WORK(1),4)
DEALLOCATE(WORK)
ALLOCATE(WORK(LWORK))
CALL dgetri(N,Mmat,N,IPIV,WORK,LWORK,INFO)
IF(INFO/=0)WRITE(*,*)INFO
DEALLOCATE(IPIV,WORK)
vert_couple=MATMUL(Mmat,vert_couple)
DEALLOCATE(Mmat)
!---Toroidally average
result=0.d0
DO j=1,nphi
  result(1)=result(1)+vert_couple(j)/REAL(nphi,8)
  result(ncoils)=result(ncoils)-vert_couple((ncoils-2)*nphi+j)/REAL(nphi,8)
  DO i=2,ncoils-1
    result(i)=result(i)+(vert_couple((i-1)*nphi+j)-vert_couple((i-2)*nphi+j))/REAL(nphi,8)
  END DO
END DO
dxm=SQRT(DOT_PRODUCT(result,result))
result=result/dxm
vert_couple=vert_couple/dxm
correction=MATMUL(corr_mat,vert_couple)
IF(nphi_grid>0)THEN
  WRITE(*,*)'  Computing grid correction'
  grid_correction=0.d0
  !$omp parallel do private(j,phi,thetas,k,rcen,btmp,l) collapse(2)
  DO i=1,ncoils-1
    DO j=1,nphi
      thetas(1)=(j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
      thetas(2)=(j)*(phi1-phi0)/REAL(nphi,8)+phi0
      !
      DO k=1,SIZE(grid_rz,DIM=2)
        DO l=1,nphi_grid
          phi = (l-1)*2.d0*pi/REAL(nphi_grid,8)
          rcen = [COS(phi)*grid_rz(1,k),SIN(phi)*grid_rz(1,k),grid_rz(2,k)]
          btmp=biot_savart_elem(rc(:,i:i+1),thetas,rcen)
          grid_correction(:,l,k) = grid_correction(:,l,k) + &
            [COS(phi)*btmp(1)+SIN(phi)*btmp(2), COS(phi)*btmp(2)-SIN(phi)*btmp(1), btmp(3)]*vert_couple((i-1)*nphi+j)
          IF(inc_mirror)THEN
            phi=phi+pi
            rcen = [COS(phi)*grid_rz(1,k),SIN(phi)*grid_rz(1,k),grid_rz(2,k)]
            btmp=biot_savart_elem(rc(:,i:i+1),thetas,rcen)
            grid_correction(:,l,k) = grid_correction(:,l,k) + &
              + [COS(phi)*btmp(1)+SIN(phi)*btmp(2), COS(phi)*btmp(2)-SIN(phi)*btmp(1), btmp(3)]*vert_couple((i-1)*nphi+j)
          END IF
        END DO
      END DO
    END DO
  END DO
END IF
WRITE(*,*)'  Finished'
DEALLOCATE(vert_couple,corr_mat)
END SUBROUTINE nonax_indcouple
!---------------------------------------------------------------------------------
!> Needs Docs
!---------------------------------------------------------------------------------
SUBROUTINE nonax_eigs(ncoils,rc,extent,result,eig_ind,correction,sensors,nsensor, &
                      inc_mirror,eta,grid_rz,grid_correction)
INTEGER(4), INTENT(in) :: ncoils
REAL(8), INTENT(in) :: rc(2,ncoils),extent(2)
REAL(8), INTENT(out) :: result(ncoils),correction(nsensor)
REAL(8), INTENT(in) :: sensors(6,nsensor)
INTEGER(4), INTENT(in) :: eig_ind,nsensor
LOGICAL, INTENT(in) :: inc_mirror
REAL(8), OPTIONAL, INTENT(in) :: eta(ncoils)
REAL(8), OPTIONAL, INTENT(in) :: grid_rz(:,:)
REAL(8), OPTIONAL, INTENT(inout) :: grid_correction(:,:,:)
!---
INTEGER(4) :: N,LDVL,LDVR,LDA,LWORK
INTEGER(4), ALLOCATABLE, DIMENSION(:) :: IPIV
REAL(8), ALLOCATABLE, DIMENSION(:) :: WR,WI,WORK
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: VL,VR
CHARACTER(LEN=1) :: JOBVL,JOBVR
!---
INTEGER(4) :: i,j,k,kk,l,info,nphi_grid
INTEGER(4), PARAMETER :: ntheta=3,nphi=40
REAL(8) :: norm(3),phi,phi0,phi1,dx(3),dy(3),btmp(3)
REAL(8) :: relem_row(3,4),thetas(2),rcen(3),dym,dxm
REAL(8), ALLOCATABLE, DIMENSION(:) :: etatmp,eig_vec
REAL(8), ALLOCATABLE, DIMENSION(:,:) :: Rmat,Lmat,dltmp,corr_mat
WRITE(*,*)'Computing non-axisymmetric eigenvalues'
!--- Create coupling matrix
phi0=extent(1)
phi1=extent(2)
ALLOCATE(Rmat((ncoils-1)*nphi,(ncoils-1)*nphi),Lmat((ncoils-1)*nphi,(ncoils-1)*nphi))
ALLOCATE(etatmp(ncoils))
ALLOCATE(dltmp(4,(ncoils-1)*nphi))
etatmp=1.d0
IF(PRESENT(eta))etatmp=eta
IF(PRESENT(grid_correction))THEN
  IF(.NOT.PRESENT(grid_rz))CALL oft_abort('Grid correction requires "grid_rz"','nonax_indcouple',__FILE__)
  nphi_grid=SIZE(grid_correction,DIM=2)
ELSE
  nphi_grid=0
END IF
Rmat=0.d0
Lmat=0.d0
dltmp=0.d0
!
ALLOCATE(corr_mat(nsensor,(ncoils-1)*nphi))
corr_mat=0.d0
! Inductance matrix
!$omp parallel do private(j,k,kk,phi,relem_row,rcen,dx,dy,norm,dxm,dym,thetas,btmp)
DO i=1,ncoils-1
  DO j=1,nphi
    phi = (j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
    relem_row(:,1)=(/COS(phi)*rc(1,i),SIN(phi)*rc(1,i),rc(2,i)/)
    phi = (j)*(phi1-phi0)/REAL(nphi,8)+phi0
    relem_row(:,2)=(/COS(phi)*rc(1,i),SIN(phi)*rc(1,i),rc(2,i)/)
    relem_row(:,3)=(/COS(phi)*rc(1,i+1),SIN(phi)*rc(1,i+1),rc(2,i+1)/)
    phi = (j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
    relem_row(:,4)=(/COS(phi)*rc(1,i+1),SIN(phi)*rc(1,i+1),rc(2,i+1)/)
    rcen=(relem_row(:,1)+relem_row(:,2)+relem_row(:,3)+relem_row(:,4))/4.d0
    dx = (relem_row(:,2)+relem_row(:,3))/2.d0-(relem_row(:,1)+relem_row(:,4))/2.d0
    dy = (relem_row(:,1)+relem_row(:,2))/2.d0-(relem_row(:,3)+relem_row(:,4))/2.d0
    norm = cross_product(dx,dy)
    dxm = magnitude(dx)
    dym = magnitude(dy)
    DO k=1,ncoils-1
      DO kk=1,nphi
        thetas(1)=(kk-1)*(phi1-phi0)/REAL(nphi,8)+phi0
        thetas(2)=(kk)*(phi1-phi0)/REAL(nphi,8)+phi0
        btmp=biot_savart_elem(rc(:,k:k+1),thetas,rcen)
        Lmat((i-1)*nphi+j,(k-1)*nphi+kk) = DOT_PRODUCT(btmp,norm)
      END DO
    END DO
    !
    thetas(1)=(j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
    thetas(2)=(j)*(phi1-phi0)/REAL(nphi,8)+phi0
    DO k=1,nsensor
      phi=sensors(2,k)
      rcen=(/COS(phi)*sensors(1,k),SIN(phi)*sensors(1,k),sensors(3,k)/)
      btmp=biot_savart_elem(rc(:,i:i+1),thetas,rcen)
      corr_mat(k,(i-1)*nphi+j) = sensors(4,k)*(COS(phi)*btmp(1) + SIN(phi)*btmp(2)) &
                                + sensors(5,k)*(-SIN(phi)*btmp(1) + COS(phi)*btmp(2)) &
                                + sensors(6,k)*btmp(3)
      !
      IF(inc_mirror)THEN
        phi=phi+pi !phi1+phi0
        rcen=(/COS(phi)*sensors(1,k),SIN(phi)*sensors(1,k),sensors(3,k)/)
        btmp=biot_savart_elem(rc(:,i:i+1),thetas,rcen)
        corr_mat(k,(i-1)*nphi+j) = corr_mat(k,(i-1)*nphi+j) &
                                  + sensors(4,k)*(COS(phi)*btmp(1) + SIN(phi)*btmp(2)) &
                                  + sensors(5,k)*(-SIN(phi)*btmp(1) + COS(phi)*btmp(2)) &
                                  + sensors(6,k)*btmp(3)
      END IF
    END DO
  END DO
END DO
! Resistance matrix
!$omp parallel do private(j,phi,relem_row,dx,dy,dxm,dym)
DO i=1,ncoils-1
  DO j=1,nphi
    IF(j==1)THEN
      phi = (j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
      relem_row(:,1)=(/COS(phi)*rc(1,i),SIN(phi)*rc(1,i),rc(2,i)/)
      phi = (j)*(phi1-phi0)/REAL(nphi,8)+phi0
      relem_row(:,2)=(/COS(phi)*rc(1,i),SIN(phi)*rc(1,i),rc(2,i)/)
      relem_row(:,3)=(/COS(phi)*rc(1,i+1),SIN(phi)*rc(1,i+1),rc(2,i+1)/)
      phi = (j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
      relem_row(:,4)=(/COS(phi)*rc(1,i+1),SIN(phi)*rc(1,i+1),rc(2,i+1)/)
      dx = (relem_row(:,2)+relem_row(:,3))/2.d0-(relem_row(:,1)+relem_row(:,4))/2.d0
      dy = (relem_row(:,1)+relem_row(:,2))/2.d0-(relem_row(:,3)+relem_row(:,4))/2.d0
      dxm = magnitude(dx)/2.d0
      dym = magnitude(dy)/2.d0
    END IF
    dltmp(1,(i-1)*nphi+j)=dltmp(1,(i-1)*nphi+j)+dym/etatmp(i)
    dltmp(2,(i-1)*nphi+j)=dltmp(2,(i-1)*nphi+j)+dxm
    dltmp(3,(i-1)*nphi+j)=dltmp(3,(i-1)*nphi+j)+dym/etatmp(i+1)
    dltmp(4,(i-1)*nphi+j)=dltmp(4,(i-1)*nphi+j)+dxm
    IF(i>1)dltmp(3,(i-2)*nphi+j)=dltmp(3,(i-2)*nphi+j)+dym/etatmp(i)
    IF(i<ncoils-1)dltmp(1,(i)*nphi+j)=dltmp(1,(i)*nphi+j)+dym/etatmp(i+1)
    IF(j>1)dltmp(2,(i-1)*nphi+j-1)=dltmp(2,(i-1)*nphi+j-1)+dxm
    IF(j<nphi)dltmp(4,(i-1)*nphi+j+1)=dltmp(4,(i-1)*nphi+j+1)+dxm
  END DO
END DO
!$omp parallel do private(j,phi,relem_row,thetas,k,rcen,btmp)
DO i=1,ncoils-1
  DO j=1,nphi
    phi = (j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
    relem_row(:,1)=(/COS(phi)*rc(1,i),SIN(phi)*rc(1,i),rc(2,i)/)
    phi = (j)*(phi1-phi0)/REAL(nphi,8)+phi0
    relem_row(:,2)=(/COS(phi)*rc(1,i),SIN(phi)*rc(1,i),rc(2,i)/)
    relem_row(:,3)=(/COS(phi)*rc(1,i+1),SIN(phi)*rc(1,i+1),rc(2,i+1)/)
    phi = (j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
    relem_row(:,4)=(/COS(phi)*rc(1,i+1),SIN(phi)*rc(1,i+1),rc(2,i+1)/)
    !---
    Rmat((i-1)*nphi+j,(i-1)*nphi+j)= magnitude(relem_row(:,2)-relem_row(:,1))/dltmp(1,(i-1)*nphi+j) &
                                    + magnitude(relem_row(:,3)-relem_row(:,2))/dltmp(2,(i-1)*nphi+j) &
                                    + magnitude(relem_row(:,4)-relem_row(:,3))/dltmp(3,(i-1)*nphi+j) &
                                    + magnitude(relem_row(:,1)-relem_row(:,4))/dltmp(4,(i-1)*nphi+j)
    ! Neighbor contributions
    IF(i>1)THEN
      Rmat((i-2)*nphi+j,(i-1)*nphi+j)=Rmat((i-2)*nphi+j,(i-1)*nphi+j) &
      - magnitude(relem_row(:,2)-relem_row(:,1))/dltmp(1,(i-1)*nphi+j)
    END IF
    IF(i<ncoils-1)THEN
      Rmat((i)*nphi+j,(i-1)*nphi+j)=Rmat((i)*nphi+j,(i-1)*nphi+j) &
      - magnitude(relem_row(:,4)-relem_row(:,3))/dltmp(3,(i-1)*nphi+j)
    END IF
    !
    IF(j>1)THEN
      Rmat((i-1)*nphi+j-1,(i-1)*nphi+j)=Rmat((i-1)*nphi+j-1,(i-1)*nphi+j) &
      - magnitude(relem_row(:,3)-relem_row(:,2))/dltmp(2,(i-1)*nphi+j)
    END IF
    IF(j<nphi)THEN
      Rmat((i-1)*nphi+j+1,(i-1)*nphi+j)=Rmat((i-1)*nphi+j+1,(i-1)*nphi+j) &
      - magnitude(relem_row(:,1)-relem_row(:,4))/dltmp(4,(i-1)*nphi+j)
    END IF
  END DO
END DO
DEALLOCATE(dltmp)
!---Invert Rmat
N = (ncoils-1)*nphi
LWORK = MAX(1,2*N)
ALLOCATE(IPIV(N))
CALL dgetrf(N,N,Rmat,N,IPIV,INFO)
IF(INFO/=0)WRITE(*,*)INFO
LWORK=-1
ALLOCATE(WORK(1))
CALL dgetri(N,Rmat,N,IPIV,WORK,LWORK,INFO)
LWORK=INT(WORK(1),4)
DEALLOCATE(WORK)
ALLOCATE(WORK(LWORK))
CALL dgetri(N,Rmat,N,IPIV,WORK,LWORK,INFO)
IF(INFO/=0)WRITE(*,*)INFO
Lmat=MATMUL(Rmat,Lmat)
DEALLOCATE(IPIV,WORK,Rmat)
!--- Compute eigenvalues
JOBVL = 'V'
JOBVR = 'V'
N = (ncoils-1)*nphi
LDA = N
LDVL = N
LDVR = N
ALLOCATE(WR(N),WI(N),VL(LDVL,N),VR(LDVR,N),WORK(1))
LWORK=-1
CALL DGEEV(JOBVL, JOBVR, N, Lmat, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
LWORK=INT(WORK(1),4)
DEALLOCATE(WORK)
ALLOCATE(WORK(LWORK))
CALL DGEEV(JOBVL, JOBVR, N, Lmat, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
WRITE(*,*)-1.d0/WR(1:5)
ALLOCATE(eig_vec((ncoils-1)*nphi))
eig_vec=REAL(VR(:,eig_ind),8)
DEALLOCATE(WI,WR,VL,VR,WORK,Lmat)
!---Toroidally average
result=0.d0
DO j=1,nphi
  result(1)=result(1)+eig_vec(j)/REAL(nphi,8)
  result(ncoils)=result(ncoils)-eig_vec((ncoils-2)*nphi+j)/REAL(nphi,8)
  DO i=2,ncoils-1
    result(i)=result(i)+(eig_vec((i-1)*nphi+j)-eig_vec((i-2)*nphi+j))/REAL(nphi,8)
  END DO
END DO
dxm=SQRT(DOT_PRODUCT(result,result))
result=result/dxm
eig_vec=eig_vec/dxm
correction=MATMUL(corr_mat,eig_vec)
IF(nphi_grid>0)THEN
  WRITE(*,*)'  Computing grid correction'
  grid_correction=0.d0
  !$omp parallel do private(j,phi,thetas,k,rcen,btmp,l)
  DO i=1,ncoils-1
    DO j=1,nphi
      thetas(1)=(j-1)*(phi1-phi0)/REAL(nphi,8)+phi0
      thetas(2)=(j)*(phi1-phi0)/REAL(nphi,8)+phi0
      !
      DO k=1,SIZE(grid_rz,DIM=2)
        DO l=1,nphi_grid
          phi = (l-1)*2.d0*pi/REAL(nphi_grid,8)
          rcen = [COS(phi)*grid_rz(1,k),SIN(phi)*grid_rz(1,k),grid_rz(2,k)]
          btmp=biot_savart_elem(rc(:,i:i+1),thetas,rcen)
          grid_correction(:,l,k) = grid_correction(:,l,k) + &
            [COS(phi)*btmp(1)+SIN(phi)*btmp(2), COS(phi)*btmp(2)-SIN(phi)*btmp(1), btmp(3)]*eig_vec((i-1)*nphi+j)
          IF(inc_mirror)THEN
            phi=phi+pi
            rcen = [COS(phi)*grid_rz(1,k),SIN(phi)*grid_rz(1,k),grid_rz(2,k)]
            btmp=biot_savart_elem(rc(:,i:i+1),thetas,rcen)
            grid_correction(:,l,k) = grid_correction(:,l,k) + &
              + [COS(phi)*btmp(1)+SIN(phi)*btmp(2), COS(phi)*btmp(2)-SIN(phi)*btmp(1), btmp(3)]*eig_vec((i-1)*nphi+j)
          END IF
        END DO
      END DO
    END DO
  END DO
END IF
WRITE(*,*)'  Finished'
DEALLOCATE(eig_vec,etatmp,corr_mat)
END SUBROUTINE nonax_eigs
END MODULE nonax_wall
!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> Driver program for GS equilibria
!!
!! @authors Chris Hansen
!! @date March 2014
!! @ingroup doxy_tokamaker
!------------------------------------------------------------------------------
program tokamaker_wall
USE oft_base
USE oft_sort, ONLY: sort_array
USE oft_io, ONLY: hdf5_create_file, hdf5_write, hdf5_create_group
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct_surf
USE fem_base, ONLY: oft_afem_type, oft_ml_fem_type
USE fem_composite, ONLY: oft_ml_fem_comp_type
USE oft_la_base, ONLY: oft_vector
USE oft_lag_basis, ONLY: oft_lag_setup_bmesh, oft_scalar_bfem, oft_lag_setup
USE oft_gs_profiles, ONLY: zero_flux_func
USE oft_gs, ONLY: gs_eq, gs_setup_walls, gs_cond_source, gs_vacuum_solve
USE oft_gs_fit, ONLY: fit_load, fit_constraint_ptr, gs_active
USE axi_green, ONLY: decay_eigenmodes
USE nonax_wall, ONLY: nonax_rescouple, nonax_indcouple, nonax_eigs
USE exp_geom, ONLY: exp_setup
IMPLICIT NONE
#include "local.h"
INTEGER(4) :: i,ierr,io_unit
TYPE(gs_eq) :: mygs
TYPE(multigrid_mesh) :: mg_mesh
TYPE(oft_ml_fem_type), TARGET :: ML_oft_lagrange,ML_oft_blagrange
TYPE(oft_ml_fem_comp_type), TARGET :: ML_oft_vlagrange
!---GS input options
INTEGER(4) :: order = 1
INTEGER(4) :: maxits = 30
INTEGER(4) :: ninner = 4
INTEGER(4) :: mode = 0
INTEGER(4) :: dcon_npsi = -1
INTEGER(4) :: dcon_ntheta = -1
INTEGER(4) :: eqdsk_nr = -1
INTEGER(4) :: eqdsk_nz = -1
LOGICAL :: pm = .FALSE.
REAL(8) :: urf = .3d0
REAL(8) :: nl_tol = 1.d-6
REAL(8) :: pnorm = 0.d0
REAL(8) :: alam = 0.d0
REAL(8) :: beta_mr = .3d0
REAL(8) :: f_offset = 0.d0
REAL(8) :: itor_target = -1.d0
REAL(8) :: rmin = 0.d0
REAL(8) :: R0_target = -1.d0
REAL(8) :: V0_target = -1.d99
REAL(8) :: rbounds(2) = (/-1.d99,1.d99/)
REAL(8) :: zbounds(2) = (/-1.d99,1.d99/)
REAL(8) :: eqdsk_rbounds(2) = (/-1.d99,1.d99/)
REAL(8) :: eqdsk_zbounds(2) = (/-1.d99,1.d99/)
REAL(8) :: init_r0(2)=[-1.d0,0.d0]
REAL(8) :: init_a=-1.d0
REAL(8) :: init_kappa=1.d0
REAL(8) :: init_delta=0.d0
LOGICAL :: free_boundary = .FALSE.
LOGICAL :: has_plasma = .TRUE.
LOGICAL :: save_mug = .FALSE.
LOGICAL :: fast_boundary = .TRUE.
LOGICAL :: limited_only = .FALSE.
CHARACTER(LEN=OFT_PATH_SLEN) :: coil_file = 'none'
CHARACTER(LEN=OFT_PATH_SLEN) :: limiter_file = 'none'
CHARACTER(LEN=OFT_PATH_SLEN) :: eqdsk_filename = 'gTokaMaker'
CHARACTER(LEN=39) :: eqdsk_run_info = ''
CHARACTER(LEN=OFT_PATH_SLEN) :: eqdsk_limiter_file = 'none'
!---Wall specific options
LOGICAL :: mirror_wall = .FALSE.
LOGICAL :: grid_3d = .FALSE.
NAMELIST/tokamaker_options/order,pm,mode,maxits,ninner,urf,nl_tol,itor_target,pnorm, &
alam,beta_mr,free_boundary,coil_file,limiter_file,f_offset,dcon_npsi,dcon_ntheta, &
has_plasma,rbounds,zbounds,rmin,R0_target,V0_target,save_mug,fast_boundary, &
limited_only,eqdsk_filename,eqdsk_nr,eqdsk_nz,eqdsk_rbounds,eqdsk_zbounds,eqdsk_run_info, &
eqdsk_limiter_file,init_r0,init_a,init_kappa,init_delta
NAMELIST/tokamaker_wall_options/mirror_wall,grid_3d
!------------------------------------------------------------------------------
! Initialize enviroment
!------------------------------------------------------------------------------
CALL oft_init
!------------------------------------------------------------------------------
! Load settings
!------------------------------------------------------------------------------
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,tokamaker_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr<0)CALL oft_abort('No "tokamaker_options" found in input file.', &
  'tokamaker_wall',__FILE__)
IF(ierr>0)CALL oft_abort('Error parsing "tokamaker_options" in input file.', &
  'tokamaker_wall',__FILE__)
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,tokamaker_wall_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr<0)CALL oft_abort('No "tokamaker_wall_options" found in input file.', &
  'tokamaker_wall',__FILE__)
IF(ierr>0)CALL oft_abort('Error parsing "tokamaker_wall_options" in input file.', &
  'tokamaker_wall',__FILE__)
!------------------------------------------------------------------------------
! Setup Mesh
!------------------------------------------------------------------------------
CALL multigrid_construct_surf(mg_mesh)
CALL mygs%xdmf%setup("TokaMaker")
CALL mg_mesh%smesh%setup_io(mygs%xdmf,order)
!------------------------------------------------------------------------------
! Setup Lagrange Elements
!------------------------------------------------------------------------------
CALL oft_lag_setup(mg_mesh,order,ML_oft_lagrange,ML_oft_blagrange,ML_oft_vlagrange)
CALL mygs%setup(ML_oft_blagrange)
!------------------------------------------------------------------------------
! Setup experimental geometry
!------------------------------------------------------------------------------
CALL exp_setup(mygs)
mygs%free=.TRUE.
mygs%pnorm=0.d0
mygs%alam=1.d-7
mygs%boundary_limiter=.FALSE.
mygs%has_plasma=.FALSE.
mygs%compute_chi=.FALSE.
IF(TRIM(coil_file)/='none')THEN
  mygs%coil_file=coil_file
  CALL mygs%load_coils(ignore_inmesh=.TRUE.)
  ! DO i=1,mygs%ncoils_ext
  !   mygs%coils_ext(i)%curr=0.d0 ! Zero all external coils
  ! END DO
  ! DO i=1,mygs%ncoil_regs
  !   mygs%coil_regions(i)%curr=0.d0 ! Zero all internal coil regions
  ! END DO
  DO i=1,mygs%ncond_eigs
    mygs%cond_regions(i)%weights=0.d0 ! Zero all conducting regions
  END DO
END IF
CALL gs_setup_walls(mygs,skip_load=.FALSE.)
CALL mygs%init()
ALLOCATE(zero_flux_func::mygs%I)
mygs%I%f_offset=0.d0
ALLOCATE(zero_flux_func::mygs%P)
CALL compute_eddy(mygs)
!------------------------------------------------------------------------------
! Terminate
!------------------------------------------------------------------------------
CALL oft_finalize
CONTAINS
!------------------------------------------------------------------------------
!> Needs Docs
!------------------------------------------------------------------------------
subroutine compute_eddy(self)
type(gs_eq), target, intent(inout) :: self
TYPE(fit_constraint_ptr), POINTER, DIMENSION(:) :: conlist => NULL()
INTEGER(4) :: ncons
INTEGER(4) :: i,j,k,l,istart,nsensor
INTEGER(4), PARAMETER :: nphi_3d=36
INTEGER(4), ALLOCATABLE :: theta_ind(:)
INTEGER(4), POINTER :: lctmp(:,:)
REAL(8) :: rcenter(2),tmark
REAL(8), ALLOCATABLE :: rc(:,:),eigs(:),eta(:),eig_vec(:,:),theta(:)
REAL(8), ALLOCATABLE :: outtmp(:,:),sensors(:,:),correction(:,:),corr_axi(:,:)
REAL(8), POINTER :: rz_grid(:,:),rz_correction(:,:,:),psi_vals(:)
CHARACTER(LEN=2) :: num_str,num_str2,cond_tag,eig_tag
CLASS(oft_vector), POINTER :: source_temp
!---
IF(oft_debug_print(2))THEN
  WRITE(*,'(2A)')oft_indent,'Setup internal regions:'
  CALL oft_increase_indent
END IF
!
gs_active=>self
CALL fit_load('fit.in',conlist)
ncons=SIZE(conlist)
nsensor=0
DO i=1,ncons
  CALL conlist(i)%con%setup_comp()
  nsensor=nsensor+conlist(i)%con%ncomp
END DO
ALLOCATE(sensors(6,nsensor))
nsensor=0
DO i=1,ncons
  DO j=1,conlist(i)%con%ncomp
    nsensor=nsensor+1
    sensors(1:3,nsensor)=conlist(i)%con%comp_r(:,j)
    sensors(4:6,nsensor)=conlist(i)%con%comp_n(:,j)
  END DO
END DO
IF(grid_3d)THEN
  CALL mg_mesh%smesh%tessellate(rz_grid,lctmp,order)
  ALLOCATE(rz_correction(3,nphi_3d,SIZE(rz_grid,DIM=2)))
  CALL hdf5_write(nphi_3d*1.d0,'wall_eig.rst','ngrid_3d')
END IF
!---
CALL hdf5_create_file('wall_eig.rst')
CALL self%psi%new(source_temp)
self%ncond_eigs=0
DO i=1,self%ncond_regs
  !---
  ALLOCATE(rc(2,self%cond_regions(i)%nc_quad))
  ALLOCATE(eta(self%cond_regions(i)%nc_quad))
  rc=self%cond_regions(i)%rc
  eta=1.d0
  !---
  rcenter=SUM(rc,DIM=2)/REAL(self%cond_regions(i)%nc_quad,8)
  ALLOCATE(theta(self%cond_regions(i)%nc_quad))
  ALLOCATE(theta_ind(self%cond_regions(i)%nc_quad))
  DO j=1,self%cond_regions(i)%nc_quad
    theta(j)=ATAN2(rc(2,j)-rcenter(2),rc(1,j)-rcenter(1))
    IF(theta(j)<0.d0)theta(j)=2.d0*pi+theta(j)
    theta_ind(j)=j
  END DO
  !---Find largest poloidal gap to split wall
  CALL sort_array(theta,theta_ind,self%cond_regions(i)%nc_quad)
  istart=1
  tmark=theta(1)
  DO j=1,self%cond_regions(i)%nc_quad-1
    IF(theta(j+1)-theta(j)>tmark)THEN
      istart=j+1
      tmark=theta(j+1)-theta(j)
    END IF
  END DO
  DO j=1,istart-1
    theta(j)=theta(j)+2*pi
  END DO
  CALL sort_array(theta,theta_ind,self%cond_regions(i)%nc_quad)
  DO j=1,self%cond_regions(i)%nc_quad
    rc(:,j)=self%cond_regions(i)%rc(:,theta_ind(j))
    IF(ASSOCIATED(self%set_eta))eta(j)=self%set_eta(rc(:,j),self%cond_regions(i)%id)
  END DO
  DEALLOCATE(theta)
  !---
  ALLOCATE(eigs(self%cond_regions(i)%nc_quad))
  ALLOCATE(eig_vec(self%cond_regions(i)%nc_quad,self%cond_regions(i)%nc_quad))
  eigs=0.d0
  eig_vec=0.d0
  IF(self%cond_regions(i)%neigs>0)THEN
    WRITE(num_str,'(I2.2)')i
    IF(.NOT.self%cond_regions(i)%continuous)THEN
      ALLOCATE(correction(nsensor,self%cond_regions(i)%neigs))
      ALLOCATE(corr_axi(ncons,self%cond_regions(i)%neigs))
      DO j=1,self%cond_regions(i)%neigs
        WRITE(num_str2,'(I2.2)')j
        k = self%cond_regions(i)%mind(j)
        SELECT CASE(self%cond_regions(i)%mtype(j))
          CASE(1)
            IF(grid_3d)THEN
              CALL nonax_indcouple(self%cond_regions(i)%nc_quad,rc,self%cond_regions(i)%extent, &
              eigs,self%coils_ext(k)%pt,self%coils_ext(k)%ncoils,correction(:,j),sensors,nsensor,mirror_wall, &
              rz_grid,rz_correction)
            ELSE 
              CALL nonax_indcouple(self%cond_regions(i)%nc_quad,rc,self%cond_regions(i)%extent, &
              eigs,self%coils_ext(k)%pt,self%coils_ext(k)%ncoils,correction(:,j),sensors,nsensor,mirror_wall)
            END IF
          CASE(2)
            IF(grid_3d)THEN
              CALL nonax_rescouple(self%cond_regions(i)%nc_quad,rc,self%cond_regions(i)%extent, &
              eigs,self%coils_ext(k)%pt,self%coils_ext(k)%ncoils,correction(:,j),sensors,nsensor,mirror_wall, &
              eta,rz_grid,rz_correction)
            ELSE 
              CALL nonax_rescouple(self%cond_regions(i)%nc_quad,rc,self%cond_regions(i)%extent, &
              eigs,self%coils_ext(k)%pt,self%coils_ext(k)%ncoils,correction(:,j),sensors,nsensor,mirror_wall,eta)
            END IF
          CASE(3)
            IF(grid_3d)THEN
              CALL nonax_eigs(self%cond_regions(i)%nc_quad,rc,self%cond_regions(i)%extent,eigs, &
              k,correction(:,j),sensors,nsensor,mirror_wall,eta,rz_grid,rz_correction)
            ELSE 
              CALL nonax_eigs(self%cond_regions(i)%nc_quad,rc,self%cond_regions(i)%extent,eigs, &
              k,correction(:,j),sensors,nsensor,mirror_wall,eta)
            END IF
          CASE DEFAULT
            CALL oft_abort("Invalid mode type specified","compute_eddy",__FILE__)
        END SELECT
        correction(:,j) = correction(:,j)/self%cond_regions(i)%coverage
        IF(grid_3d)THEN
          rz_correction=rz_correction/self%cond_regions(i)%coverage
          ALLOCATE(outtmp(3,SIZE(rz_grid,DIM=2)))
          outtmp=SQRT(SUM(rz_correction**2,DIM=2)/REAL(SIZE(rz_grid,DIM=2),8))
          CALL mg_mesh%smesh%save_vertex_scalar(outtmp(1,:),self%xdmf, 'Br_corr')
          CALL mg_mesh%smesh%save_vertex_scalar(outtmp(2,:),self%xdmf, 'Bt_corr')
          CALL mg_mesh%smesh%save_vertex_scalar(outtmp(3,:),self%xdmf, 'Bz_corr')
          DEALLOCATE(outtmp)
          CALL hdf5_write(rz_correction(1,:,:), 'wall_eig.rst', 'rz_corr_r'//num_str//'_'//num_str2)
          CALL hdf5_write(rz_correction(2,:,:), 'wall_eig.rst', 'rz_corr_t'//num_str//'_'//num_str2)
          CALL hdf5_write(rz_correction(3,:,:), 'wall_eig.rst', 'rz_corr_z'//num_str//'_'//num_str2)
        END IF
        eig_vec(:,j) = eigs
        !---Get axisymmetric effect
        DO k=1,self%cond_regions(i)%nc_quad
          self%cond_regions(i)%cond_vals(theta_ind(k),j) = eigs(k)/4.d0
        END DO
        self%cond_regions(i)%weights(j)=1.d0
        CALL gs_cond_source(self,i,j,source_temp)
        tmark=source_temp%dot(source_temp)
        WRITE(*,*)tmark
        CALL gs_vacuum_solve(self,self%psi,source_temp)
        tmark=self%psi%dot(self%psi)
        WRITE(*,*)tmark
        NULLIFY(psi_vals)
        CALL self%psi%get_local(psi_vals)
        WRITE(cond_tag,'(I2.2)')i
        WRITE(eig_tag,'(I2.2)')j
        CALL mg_mesh%smesh%save_vertex_scalar(psi_vals,self%xdmf,'Eig_'//cond_tag//'_'//eig_tag)
        DEALLOCATE(psi_vals)
        ! CALL self%solve
        DO k=1,ncons
          IF(conlist(k)%con%ncomp==0)CYCLE
          corr_axi(k,j)=conlist(k)%con%eval(self)
        END DO
        self%cond_regions(i)%weights(j)=0.d0
      END DO
      nsensor=0
      DO k=1,ncons
        IF(conlist(k)%con%ncomp==0)CYCLE
        nsensor=nsensor+1
      END DO
      ALLOCATE(outtmp(self%cond_regions(i)%neigs,nsensor))
      outtmp=0.d0
      nsensor=0
      l=1
      DO k=1,ncons
        IF(conlist(k)%con%ncomp==0)CYCLE
        DO j=1,conlist(k)%con%ncomp
          nsensor=nsensor+1
          outtmp(:,l)=outtmp(:,l)+correction(nsensor,:)
        END DO
        outtmp(:,l)=outtmp(:,l)-corr_axi(k,:)
        l=l+1
      END DO
      CALL hdf5_write(outtmp, 'wall_eig.rst', 'corr_'//num_str)
      DEALLOCATE(correction,corr_axi,outtmp)
    ELSE
      CALL decay_eigenmodes(self%cond_regions(i)%nc_quad,rc,eigs,eig_vec,eta)
      DO j=1,self%cond_regions(i)%neigs
        eig_vec(:,j) = eig_vec(:,self%cond_regions(i)%mind(j))
      END DO
    END IF
    ALLOCATE(outtmp(self%cond_regions(i)%nc_quad,self%cond_regions(i)%neigs))
    DO j=1,self%cond_regions(i)%nc_quad
      outtmp(theta_ind(j),:)=eig_vec(j,1:self%cond_regions(i)%neigs)
    END DO
    CALL hdf5_write(outtmp, 'wall_eig.rst', 'eig_'//num_str)
    DEALLOCATE(outtmp)
  END IF
  !---
  IF(oft_debug_print(2))THEN
    WRITE(*,'(2A,I6,A)')oft_indent,'Found ',self%cond_regions(i)%nc,' conductor cells'
    WRITE(*,'(2A,2E11.3)')oft_indent,'Rcenter = ',REAL(rcenter,4)
  END IF
  DEALLOCATE(rc,eigs,eta,eig_vec,theta_ind)
END DO
IF(grid_3d)DEALLOCATE(rz_grid,rz_correction,lctmp)
IF(oft_debug_print(2))CALL oft_decrease_indent
end subroutine compute_eddy
end program tokamaker_wall
