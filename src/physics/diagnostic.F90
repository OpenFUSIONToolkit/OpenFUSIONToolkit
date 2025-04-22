!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
! MODULE: diagnostic
!
!> @file diagnostic.F90
!
!> Classes and subroutines used for synthetic diagnostics.
!!
!! @authors Chris Hansen
!! @date March 2013
!! @ingroup doxy_oft_physics
!------------------------------------------------------------------------------
MODULE diagnostic
USE oft_base
USE oft_quadrature
USE oft_mesh_type, ONLY: oft_mesh, oft_bmesh, mesh_findcell2, cell_is_curved
USE oft_io, ONLY: oft_bin_file
USE fem_utils, ONLY: fem_interp, bfem_interp
IMPLICIT NONE
#include "local.h"
!---------------------------------------------------------------------------------
!> Synthetic field diagnostic
!!
!! Provides setup and driver for sampling a n-vector at specified locations.
!! Locations are mapped during the setup procedure so that subsequent calls to
!! eval have very low overhead.
!---------------------------------------------------------------------------------
TYPE :: field_probe
  INTEGER(i4) :: dim = 3 !< Dimension of field being sampled
  INTEGER(i4) :: nprobes = 0 !< Number of probe locations
  INTEGER(i4), POINTER, DIMENSION(:) :: cells => NULL() !< Cell containing point
  REAL(r8) :: tol=1.d-8 !< Off-mesh tolerance in logical coordinates
  REAL(r8), POINTER, DIMENSION(:,:) :: pts => NULL() !< Location of points in space
  REAL(r8), POINTER, DIMENSION(:,:) :: pts_log => NULL() !< Location of points in logical space
  LOGICAL :: project = .FALSE. !< Project points to mesh boundary
  CHARACTER(LEN=OFT_PATH_SLEN) :: filename = '' !< History file name
  CLASS(fem_interp), POINTER :: B => NULL() !< Interpolation operator for field
  CLASS(oft_mesh), POINTER :: mesh => NULL() !< Mesh containing field
  TYPE(oft_bin_file) :: hist_file
CONTAINS
  !> Initialize point list and setup ownership
  PROCEDURE :: setup => field_probe_setup
  !> Evalute field at all probe locations
  PROCEDURE :: eval => field_probe_eval
  !> Setup history file for repeated saves
  PROCEDURE :: setup_save => field_probe_setup_save
  !> Sample and save the result to a history file
  PROCEDURE :: save => field_probe_save
END TYPE field_probe
!---------------------------------------------------------------------------------
!> Synthetic flux diagnostic
!!
!! Provides setup and driver for computing the flux of a 3-vector at a specified
!! surface.
!!
!! @Note This method defines the surface using a circular cut plane and requires
!! an internal boundary at the desired surface.
!!
!! @warning Currently this diagnostic does not take into account mesh curvature.
!! This may result in errors when the internal boundary intersects a curved boundary.
!---------------------------------------------------------------------------------
TYPE :: flux_probe
  INTEGER(i4) :: order = 3 !< Order of quadrature required
  INTEGER(i4) :: nf = 0 !< Number of surface faces
  INTEGER(i4), POINTER, DIMENSION(:) :: cells => NULL() !< List of cells containing face
  INTEGER(i4), POINTER, DIMENSION(:) :: find => NULL() !< List of face indices in cell
  REAL(r8) :: tol=1.d-8 !< Tolerance for on surface test
  REAL(r8) :: hcpc(3) = 0.d0 !< Center coordinate for circular cut plane
  REAL(r8) :: hcpv(3) = 0.d0 !< Normal vector for circular cut plane
  REAL(r8), POINTER, DIMENSION(:,:) :: norm => NULL() !< List of face normals
  CLASS(fem_interp), POINTER :: B => NULL() !< Interpolation operator for field
  CLASS(oft_mesh), POINTER :: mesh => NULL() !< Mesh containing field
CONTAINS
  !> Determine surface and setup mappings
  PROCEDURE :: setup => flux_probe_setup
  !> Evalute flux through surface
  PROCEDURE :: eval => flux_probe_eval
END TYPE flux_probe
CONTAINS
!------------------------------------------------------------------------------
!> Initialize point list and setup ownership
!!
!! Sampling locations may be set in the code directly, via `pts`, or loaded
!! from a file, via `filename`. If `filename` is specified the number of points
!! read from the file is returned in `npts`.
!------------------------------------------------------------------------------
SUBROUTINE field_probe_setup(self,npts,pts,filename)
CLASS(field_probe), INTENT(inout) :: self
INTEGER(i4), INTENT(inout) :: npts !< Number of probes
REAL(r8), OPTIONAL, INTENT(in) :: pts(:,:) !< Array of probe locations [3,npts] (optional)
CHARACTER(LEN=*), OPTIONAL, INTENT(in) :: filename !< File containing probe locations (optional)
!---
INTEGER(i4) :: i,cell,ierr,io_unit,minf
INTEGER(i4), ALLOCATABLE :: ptmp(:),pown(:)
REAL(r8) :: f(4),fmin,fmax,ftol
REAL(r8), ALLOCATABLE :: fout(:)
LOGICAL :: file_exist
!---
IF(PRESENT(pts))THEN
  self%nprobes=npts
  ALLOCATE(self%pts(3,self%nprobes))
  DO i=1,self%nprobes
    self%pts(:,i)=pts(:,i)
  END DO
ELSE IF(PRESENT(filename))THEN
  INQUIRE(FILE=TRIM(filename),EXIST=file_exist)
  IF(.NOT.file_exist)CALL oft_abort('Probe file does not exist.','field_probe_setup',__FILE__)
  !---
  OPEN(NEWUNIT=io_unit,FILE=TRIM(filename))
  READ(io_unit,*)npts
  self%nprobes=npts
  ALLOCATE(self%pts(3,self%nprobes))
  DO i=1,self%nprobes
    READ(io_unit,*)self%pts(:,i)
  END DO
  CLOSE(io_unit)
ELSE
  CALL oft_abort('Neither locations nor a probe file were specified.','field_probe_setup',__FILE__)
END IF
IF(self%project)CALL project_points_to_boundary(self%mesh,self%nprobes,self%pts)
!---
ALLOCATE(self%cells(self%nprobes),self%pts_log(4,self%nprobes))
ALLOCATE(ptmp(self%nprobes),pown(self%nprobes),fout(self%nprobes))
ptmp=-1
self%cells=0
DO i=1,self%nprobes
  CALL mesh_findcell2(self%mesh,self%cells(i),self%pts(:,i),4,f)
  minf=self%mesh%in_cell(f, self%tol)
  ! fmin=MINVAL(f)
  ! fmax=MAXVAL(f)
  ! fout(i)=MAX(-fmin,fmax-1.d0)
  !---Disable point if off mesh
  ! IF(.NOT.(( fmax<=1.d0+self%tol ).AND.( fmin>=-self%tol )))THEN
  IF(minf/=0)THEN
    fmin=MINVAL(f(1:3))
    fmax=MAXVAL(f(1:3))
    fout(i)=MAX(-fmin,fmax-1.d0)
    self%cells(i)=0
  ELSE
    self%pts_log(:,i)=f
    ptmp(i)=oft_env%rank
  END IF
END DO
!---Ensure single ownership
pown=oft_mpi_max(ptmp,self%nprobes)
DO i=1,self%nprobes
  !---Report points outside mesh
  IF(pown(i)==-1)THEN
#ifdef HAVE_MPI
    CALL MPI_Reduce(fout(i),ftol,1,OFT_MPI_R8,MPI_MIN,0,oft_env%COMM,ierr)
#else
    ftol=fout(i)
#endif
    IF(oft_env%rank==0)THEN
      WRITE(*,*)'Field probe point ',INT(i,2),' offmesh'
      WRITE(*,*)'  fout = ',REAL(ftol,4),REAL(self%tol,4)
      WRITE(*,*)'  r    = ',REAL(self%pts(:,i),4)
    END IF
  END IF
  IF(self%cells(i)==0)CYCLE
  !---Disable point if not owned
  IF(.NOT.(pown(i)==oft_env%rank))THEN
    self%cells(i)=0
  END IF
END DO
DEALLOCATE(ptmp,pown,fout)
END SUBROUTINE field_probe_setup
!------------------------------------------------------------------------------
!> Evalute field at all probe locations
!------------------------------------------------------------------------------
SUBROUTINE field_probe_eval(self,vals)
CLASS(field_probe), INTENT(inout) :: self
REAL(r8), INTENT(inout) :: vals(:,:) !< Fields at all probe locations [3,npts]
!---
INTEGER(i4) :: i,ierr
REAL(r8) :: goptmp(3,4),v
REAL(r8), ALLOCATABLE :: vtmp(:,:)
!---
ALLOCATE(vtmp(self%dim,self%nprobes))
vtmp=0.d0
DO i=1,self%nprobes
  IF(self%cells(i)>0)THEN
    CALL self%mesh%jacobian(self%cells(i),self%pts_log(:,i),goptmp,v)
    CALL self%B%interp(self%cells(i),self%pts_log(:,i),goptmp,vtmp(:,i))
  END IF
END DO
#ifdef HAVE_MPI
CALL MPI_ALLREDUCE(vtmp,vals,self%dim*self%nprobes,OFT_MPI_R8,MPI_SUM,oft_env%COMM,ierr)
#else
vals=vtmp
#endif
DEALLOCATE(vtmp)
END SUBROUTINE field_probe_eval
!------------------------------------------------------------------------------
!> Setup history file for repeated saves
!------------------------------------------------------------------------------
SUBROUTINE field_probe_setup_save(self,filename)
CLASS(field_probe), INTENT(inout) :: self
CHARACTER(LEN=*), INTENT(in) :: filename !< Filename for history file
!---
INTEGER(i4) :: i,j,io_unit
CHARACTER(LEN=4) :: fid
CHARACTER(LEN=2) :: did
CHARACTER(LEN=OFT_SLEN) :: desc_str
self%filename=filename
!---Setup history file
IF(oft_env%head_proc)THEN
100 FORMAT('B-Field? packed(',I2,',',I4,')')
  WRITE(desc_str,100)self%dim, self%nprobes
  CALL self%hist_file%setup(filename)
  CALL self%hist_file%add_field('time', 'r8', desc="Simulation time [s]")
  CALL self%hist_file%add_field('B', 'r8', desc=desc_str, fsize=self%dim*self%nprobes)
  CALL self%hist_file%write_header
  CALL self%hist_file%open
  CALL self%hist_file%close
END IF
END SUBROUTINE field_probe_setup_save
!------------------------------------------------------------------------------
!> Sample and save the result to the history file
!------------------------------------------------------------------------------
SUBROUTINE field_probe_save(self,time)
CLASS(field_probe), INTENT(inout) :: self
REAL(r8), INTENT(in) :: time !< Time of signal sample
INTEGER(i4) :: nhist,io_unit
REAL(r8), ALLOCATABLE :: vals(:,:),output(:)
!---Sample field
nhist=self%dim*self%nprobes + 1
ALLOCATE(vals(self%dim,self%nprobes),output(nhist))
CALL self%eval(vals)
!---Output to history file
IF(oft_env%head_proc)THEN
  output(1)=time
  output(2:nhist)=RESHAPE(vals,(/nhist-1/))
  CALL self%hist_file%open
  CALL self%hist_file%write(data_r8=output)
  CALL self%hist_file%close
END IF
DEALLOCATE(vals,output)
END SUBROUTINE field_probe_save
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE flux_probe_setup(self)
CLASS(flux_probe), INTENT(inout) :: self
!---
INTEGER(i4) :: i,j,cell
INTEGER(i4), ALLOCATABLE :: fmap(:)
REAL(r8) :: r(3),rcc(3),rdv(3)
CLASS(oft_mesh), POINTER :: mesh
!---
mesh=>self%mesh
allocate(fmap(mesh%nf))
CALL get_inverse_map(mesh%lbf,mesh%nbf,fmap,mesh%nf)
self%nf=0
do i=1,mesh%nf
  rcc=0.d0
  do j=1,3
    r=mesh%r(:,mesh%lf(j,i))-self%hcpc
    rcc=rcc+r/3.d0
    rdv(j)=DOT_PRODUCT(r,self%hcpv)
  end do
  if(ALL(ABS(rdv)<self%tol))then
    if(SUM(cross_product(rcc,self%hcpv)**2)<1.d0)then
      if(fmap(i)>0)then
        if(.NOT.mesh%fstitch%leo(fmap(i)))cycle
      end if
      self%nf=self%nf+1
    end if
  end if
end do
IF(self%nf>0)THEN
  ALLOCATE(self%cells(self%nf),self%find(self%nf))
  ALLOCATE(self%norm(3,self%nf))
  self%nf=0
  do i=1,mesh%nf
    rcc=0.d0
    do j=1,3
      r=mesh%r(:,mesh%lf(j,i))-self%hcpc
      rcc=rcc+r/3.d0
      rdv(j)=DOT_PRODUCT(r,self%hcpv)
    end do
    if(ALL(ABS(rdv)<self%tol))then
      if(SUM(cross_product(rcc,self%hcpv)**2)<1.d0)THEN
        if(fmap(i)>0)then
          if(.NOT.mesh%fstitch%leo(fmap(i)))cycle
        end if
        cell=mesh%lfc(1,i)
        self%nf=self%nf+1
        self%cells(self%nf)=cell
        DO j=1,4
          IF(ABS(mesh%lcf(j,cell))==i)EXIT
        END DO
        self%find(self%nf)=j
        !---
        r=mesh%r(:,mesh%lf(2,i))-mesh%r(:,mesh%lf(1,i))
        rcc=mesh%r(:,mesh%lf(3,i))-mesh%r(:,mesh%lf(1,i))
        rcc=cross_product(r,rcc)/2.d0
        self%norm(:,self%nf)=rcc*SIGN(1.d0,DOT_PRODUCT(rcc,self%hcpv))
      END IF
    end if
  end do
END IF
deallocate(fmap)
END SUBROUTINE flux_probe_setup
!------------------------------------------------------------------------------
!> Needs docs
!------------------------------------------------------------------------------
SUBROUTINE flux_probe_eval(self,tflux)
CLASS(flux_probe), INTENT(inout) :: self
REAL(r8), INTENT(inout) :: tflux
!---
INTEGER(i4) :: i,j,m,fmap(3)
REAL(r8) :: goptmp(3,4),v,bcc(3),f(4),reg
TYPE(oft_quad_type) :: quad
!---Set quadrature order
CALL self%mesh%bmesh%quad_rule(self%order,quad)
!---
reg=0.d0
DO i=1,self%nf
  j=1
  DO m=1,4
    IF(self%find(i)==m)CYCLE
    fmap(j)=m
    j=j+1
  END DO
  f=0.d0
  DO m=1,quad%np
    f(fmap)=quad%pts(:,m)
    CALL self%mesh%jacobian(self%cells(i),f,goptmp,v)
    CALL self%B%interp(self%cells(i),f,goptmp,bcc)
    reg=reg+DOT_PRODUCT(bcc,self%norm(:,i))*quad%wts(m)
  END DO
END DO
!---Reduction over processors
tflux=oft_mpi_sum(reg)
!---Delete quadrature object
CALL quad%delete
END SUBROUTINE flux_probe_eval
!------------------------------------------------------------------------------
!> Project a set of points to the mesh boundary
!!
!! Projection is performed by finding the closest point to a set of known points
!! on the boundary mesh. Boundary points are defined by a given 2D quadrature
!! rule. This provides a relatively evenly spaced set of points over each boundary
!! triangle.
!------------------------------------------------------------------------------
SUBROUTINE project_points_to_boundary(mesh,npts,pts,order)
CLASS(oft_mesh), INTENT(inout) :: mesh
INTEGER(i4), INTENT(in) :: npts !< Number of points
REAL(r8), INTENT(inout) :: pts(:,:) !< List of points [3,npts]
INTEGER(i4), OPTIONAL, INTENT(in) :: order !< Order of 2D quadrature rule used (optional)
!---
INTEGER(i4) :: i,j,k,l,m,ierr,quad_order
REAL(r8) :: dtmp,ptmp(3)
REAL(r8), ALLOCATABLE :: dist(:),distout(:),distin(:),ptstmp(:,:)
TYPE(oft_quad_type) :: quad
#ifdef HAVE_MPI
#ifdef OFT_MPI_F08
TYPE(mpi_status) :: stat
#else
INTEGER(i4) :: stat(MPI_STATUS_SIZE)
#endif
#endif
!---
quad_order=4
IF(PRESENT(order))quad_order=order
CALL mesh%bmesh%quad_rule(quad_order,quad)
ALLOCATE(dist(npts),ptstmp(3,npts))
dist=1.d99
DO i=1,npts
  DO j=1,mesh%bmesh%nc
    DO m=1,quad%np
      ptmp=mesh%bmesh%log2phys(j,quad%pts(:,m))
      dtmp=SUM((ptmp-pts(:,i))**2)
      IF(dtmp<dist(i))THEN
        dist(i)=dtmp
        ptstmp(:,i)=ptmp
      END IF
    END DO
  END DO
END DO
if(mesh%fullmesh)then
  pts=ptstmp
  CALL quad%delete
  DEALLOCATE(ptstmp,dist)
  DEBUG_STACK_POP
  return
endif
ALLOCATE(distout(npts),distin(npts))
distout=dist
!---Initialize counters for ring update
j=oft_env%rank+1
if(j>oft_env%nprocs-1)j=0
k=oft_env%rank-1
if(k<0)k=oft_env%nprocs-1
l=k
CALL oft_mpi_barrier(ierr) ! Wait for all processes
!---Loop over processors to find global igrnd
#ifdef HAVE_MPI
DO WHILE(l.NE.oft_env%rank)
  CALL MPI_SEND(distout,npts,OFT_MPI_R8,j,1,oft_env%COMM,ierr)
  CALL MPI_RECV(distin,npts,OFT_MPI_R8,k,1,oft_env%COMM,stat,ierr)
  distout=distin
  !---Test if igrnd is on local mesh and update
  DO i=1,npts
    IF(distin(i)<dist(i).OR.(l<oft_env%rank.AND.distin(i)==dist(i)))dist(i)=-dist(i)
  END DO
  !---Increment ring comm
  l=l-1
  if(l<0)l=oft_env%nprocs-1
END DO
#endif
!---
DO i=1,npts
  IF(dist(i)<0.d0)ptstmp(:,i)=0.d0
END DO
!---Ensure single ownership
#ifdef HAVE_MPI
CALL MPI_ALLREDUCE(ptstmp,pts,3*npts,OFT_MPI_R8,MPI_SUM,oft_env%COMM,ierr)
#else
pts=ptstmp
#endif
CALL quad%delete
DEALLOCATE(ptstmp,dist,distin,distout)
END SUBROUTINE project_points_to_boundary
!------------------------------------------------------------------------------
!> Evaluate the toroidally averaged toroidal flux of a 3-vector
!!
!! @note This requires your geometry is oriented with one of the principle axes
!! as the axis of toroidal symmetry.
!------------------------------------------------------------------------------
FUNCTION tfluxfun(mesh,field,quad_order,axis) RESULT(tflux)
CLASS(oft_mesh), INTENT(inout) :: mesh
CLASS(fem_interp), INTENT(inout) :: field !< Input field
INTEGER(i4), INTENT(in) :: quad_order !< Desired quadrature order
INTEGER(i4), OPTIONAL, INTENT(in) :: axis !< Index of axis coordinate (optional)
REAL(r8) :: tflux !< Toroidally averaged flux
LOGICAL :: curved
INTEGER(i4) :: i,m,raxis,ptind(2)
REAL(r8) :: goptmp(3,4),vol,det,pt(3),bcc(3),rop(3)
TYPE(oft_quad_type) :: quad
DEBUG_STACK_PUSH
!---
ptind=(/1,2/)
raxis=3
IF(PRESENT(axis))raxis=axis
IF(raxis==1)ptind(1)=3
IF(raxis==2)ptind(2)=3
!---Set quadrature order
CALL mesh%quad_rule(quad_order,quad)
tflux=0.d0
!$omp parallel do default(firstprivate) shared(field,quad,raxis,ptind) private(curved) reduction(+:tflux)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i)
  DO m=1,quad%np
    IF(curved.OR.(m==1))CALL mesh%jacobian(i,quad%pts(:,m),goptmp,vol)
    pt=mesh%log2phys(i,quad%pts(:,m))
    CALL field%interp(i,quad%pts(:,m),goptmp,bcc)
    bcc=cross_product(pt,bcc)
    tflux = tflux + bcc(raxis)/SUM(pt(ptind)**2)*vol*quad%wts(m)
  END DO
END DO
!---Global reduction and cleanup
tflux=oft_mpi_sum(tflux)/(2*pi)
CALL quad%delete
DEBUG_STACK_POP
END FUNCTION tfluxfun
!------------------------------------------------------------------------------
!> Evaluate the volume integral of a scalar
!------------------------------------------------------------------------------
FUNCTION scal_int(mesh,field,quad_order) RESULT(energy)
CLASS(oft_mesh), INTENT(inout) :: mesh
CLASS(fem_interp), INTENT(inout) :: field !< Input field
INTEGER(i4), INTENT(in) :: quad_order !< Desired quadrature order
REAL(r8) :: energy !< \f$ \int u dV \f$
LOGICAL :: curved
INTEGER(i4) :: i,m
REAL(r8) :: goptmp(3,4),vol,bcc(1)
TYPE(oft_quad_type) :: quad
DEBUG_STACK_PUSH
!---Setup
CALL mesh%quad_rule(quad_order,quad)
energy=0.d0
!$omp parallel do default(firstprivate) shared(field,quad) private(curved) reduction(+:energy)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i)
  DO m=1,quad%np
    IF(curved.OR.(m==1))CALL mesh%jacobian(i,quad%pts(:,m),goptmp,vol)
    CALL field%interp(i,quad%pts(:,m),goptmp,bcc)
    energy=energy+bcc(1)*vol*quad%wts(m)
  END DO
END DO
!---Global reduction and cleanup
energy=oft_mpi_sum(energy)
CALL quad%delete
DEBUG_STACK_POP
END FUNCTION scal_int
!------------------------------------------------------------------------------
!> Evaluate the field energy of a scalar
!------------------------------------------------------------------------------
FUNCTION scal_energy(mesh,field,quad_order) RESULT(energy)
CLASS(oft_mesh), INTENT(inout) :: mesh
CLASS(fem_interp), INTENT(inout) :: field !< Input field
INTEGER(i4), INTENT(in) :: quad_order !< Desired quadrature order
REAL(r8) :: energy !< \f$ \int u^2 dV \f$
REAL(r8) :: goptmp(3,4),vol,bcc(1)
INTEGER(i4) :: i,m
LOGICAL :: curved
TYPE(oft_quad_type) :: quad
DEBUG_STACK_PUSH
!---Setup
CALL mesh%quad_rule(quad_order,quad)
energy=0.d0
!$omp parallel do default(firstprivate) shared(field,quad) private(curved) reduction(+:energy)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i)
  DO m=1,quad%np
    IF(curved.OR.(m==1))CALL mesh%jacobian(i,quad%pts(:,m),goptmp,vol)
    CALL field%interp(i,quad%pts(:,m),goptmp,bcc)
    energy = energy + (bcc(1)**2)*vol*quad%wts(m)
  END DO
END DO
!---Global reduction and cleanup
energy=oft_mpi_sum(energy)
CALL quad%delete
DEBUG_STACK_POP
END FUNCTION scal_energy
!------------------------------------------------------------------------------
!> Evaluate the field energy of a 3-vector
!------------------------------------------------------------------------------
FUNCTION vec_energy(mesh,field,quad_order) RESULT(energy)
CLASS(oft_mesh), INTENT(inout) :: mesh
CLASS(fem_interp), INTENT(inout) :: field !< Input field
INTEGER(i4), INTENT(in) :: quad_order !< Desired quadrature order
REAL(r8) :: energy !< \f$ \int \left| \textbf{u} \right|^2 dV \f$
REAL(r8) :: goptmp(3,4),vol,bcc(3)
INTEGER(i4) :: i,m
LOGICAL :: curved
TYPE(oft_quad_type) :: quad
DEBUG_STACK_PUSH
!---Setup
CALL mesh%quad_rule(quad_order,quad)
energy=0.d0
!$omp parallel do default(firstprivate) shared(field,quad) private(curved) reduction(+:energy)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i)
  DO m=1,quad%np
    IF(curved.OR.(m==1))CALL mesh%jacobian(i,quad%pts(:,m),goptmp,vol)
    CALL field%interp(i,quad%pts(:,m),goptmp,bcc)
    energy = energy + SUM(bcc**2)*vol*quad%wts(m)
  END DO
END DO
!---Global reduction and cleanup
energy=oft_mpi_sum(energy)
CALL quad%delete
DEBUG_STACK_POP
END FUNCTION vec_energy
!------------------------------------------------------------------------------
!> Evaluate the field energy of a 3-vector with a scalar weight field
!------------------------------------------------------------------------------
FUNCTION weighted_vec_energy(mesh,field,weight,quad_order) RESULT(energy)
CLASS(oft_mesh), INTENT(inout) :: mesh
CLASS(fem_interp), INTENT(inout) :: field !< Input field \f$ (u) \f$
CLASS(fem_interp), INTENT(inout) :: weight !< Weight field \f$ (\omega) \f$
INTEGER(i4), INTENT(in) :: quad_order !< Desired quadrature order
REAL(r8) :: energy !< \f$ \int \left( \omega * \left| \textbf{u} \right|^2 \right) dV \f$
REAL(r8) :: goptmp(3,4),vol,bcc(3),wcc(1)
INTEGER(i4) :: i,m
LOGICAL :: curved
TYPE(oft_quad_type) :: quad
DEBUG_STACK_PUSH
!---Setup
CALL mesh%quad_rule(quad_order,quad)
energy=0.d0
!$omp parallel do  default(firstprivate) shared(field,weight,quad) private(curved) reduction(+:energy)
DO i=1,mesh%nc
  curved=cell_is_curved(mesh,i)
  DO m=1,quad%np
    IF(curved.OR.(m==1))CALL mesh%jacobian(i,quad%pts(:,m),goptmp,vol)
    CALL field%interp(i,quad%pts(:,m),goptmp,bcc)
    CALL weight%interp(i,quad%pts(:,m),goptmp,wcc)
    energy = energy + wcc(1)*SUM(bcc**2)*vol*quad%wts(m)
  END DO
END DO
!---Global reduction and cleanup
energy=oft_mpi_sum(energy)
CALL quad%delete
DEBUG_STACK_POP
END FUNCTION weighted_vec_energy
!------------------------------------------------------------------------------
!> Evaluate the boundary integral of a scalar field
!------------------------------------------------------------------------------
FUNCTION scal_surf_int(mesh,field,quad_order) RESULT(energy)
CLASS(oft_mesh), INTENT(inout) :: mesh
CLASS(fem_interp), INTENT(inout) :: field !< Input field
INTEGER(i4), INTENT(in) :: quad_order !< Desired quadrature order
REAL(r8) :: energy !< \f$ \int u dS \f$
INTEGER(i4) :: i,m,j,face,cell,ptmap(3)
REAL(r8) :: vol,area,flog(4),etmp(1),sgop(3,3),vgop(3,4)
TYPE(oft_quad_type) :: quad
DEBUG_STACK_PUSH
!---Setup
CALL mesh%bmesh%quad_rule(quad_order,quad)
energy=0.d0
!$omp parallel do default(firstprivate) shared(field,quad) reduction(+:energy)
do i=1,mesh%bmesh%nc
  CALL mesh%get_surf_map(i,cell,ptmap) ! Find parent cell and logical coordinate mapping
  !---Loop over quadrature points
  do m=1,quad%np
    call mesh%bmesh%jacobian(i,quad%pts(:,m),sgop,area)
    !---Evaluate in cell coordinates
    CALL mesh%surf_to_vol(quad%pts(:,m),ptmap,flog)
    call mesh%jacobian(cell,flog,vgop,vol)
    call field%interp(cell,flog,vgop,etmp)
    energy = energy + etmp(1)*area*quad%wts(m)
  end do
end do
!---Global reduction and cleanup
energy=oft_mpi_sum(energy)
CALL quad%delete
DEBUG_STACK_POP
END FUNCTION scal_surf_int
!------------------------------------------------------------------------------
!> Evaluate the boundary integral of a boundary scalar field
!------------------------------------------------------------------------------
FUNCTION bscal_surf_int(mesh,field,quad_order,reg_mask) RESULT(energy)
CLASS(oft_bmesh), INTENT(inout) :: mesh
CLASS(bfem_interp), INTENT(inout) :: field !< Input field
INTEGER(i4), INTENT(in) :: quad_order !< Desired quadrature order
INTEGER(i4), OPTIONAL, INTENT(in) :: reg_mask !< Region to integrate over
REAL(r8) :: energy !< \f$ \int u dS \f$
INTEGER(i4) :: i,m,cell
REAL(r8) :: area,etmp(1),sgop(3,3)
TYPE(oft_quad_type) :: quad
DEBUG_STACK_PUSH
!---Setup
CALL mesh%quad_rule(quad_order,quad)
energy=0.d0
!$omp parallel do default(firstprivate) shared(field,quad,reg_mask) reduction(+:energy)
do i=1,mesh%nc
  IF(PRESENT(reg_mask))THEN
    IF(mesh%reg(i)/=reg_mask)CYCLE
  END IF
  !---Loop over quadrature points
  do m=1,quad%np
    call mesh%jacobian(i,quad%pts(:,m),sgop,area)
    call field%interp(i,quad%pts(:,m),sgop,etmp)
    energy = energy + etmp(1)*area*quad%wts(m)
  end do
end do
!---Global reduction and cleanup
energy=oft_mpi_sum(energy)
CALL quad%delete
DEBUG_STACK_POP
END FUNCTION bscal_surf_int
!------------------------------------------------------------------------------
!> Evaluate the boundary flux of a vector field
!------------------------------------------------------------------------------
FUNCTION vec_surf_int(mesh,field,quad_order) RESULT(energy)
CLASS(oft_mesh), INTENT(inout) :: mesh
CLASS(fem_interp), INTENT(inout) :: field !< Input field
INTEGER(i4), INTENT(in) :: quad_order !< Desired quadrature order
REAL(r8) :: energy !< \f$ \int \textbf{u} \cdot \textbf{dS} \f$
INTEGER(i4) :: i,m,j,face,cell,ptmap(3)
REAL(r8) :: vol,area,flog(4),norm(3),etmp(3),sgop(3,3),vgop(3,4)
TYPE(oft_quad_type) :: quad
DEBUG_STACK_PUSH
!---Setup
CALL mesh%bmesh%quad_rule(quad_order,quad)
energy=0.d0
!$omp parallel do default(firstprivate) shared(field,quad) reduction(+:energy)
do i=1,mesh%bmesh%nc
  CALL mesh%get_surf_map(i,cell,ptmap) ! Find parent cell and logical coordinate mapping
  !---Loop over quadrature points
  do m=1,quad%np
    call mesh%bmesh%jacobian(i,quad%pts(:,m),sgop,area)
    call mesh%bmesh%norm(i,quad%pts(:,m),norm)
    !---Evaluate in cell coordinates
    CALL mesh%surf_to_vol(quad%pts(:,m),ptmap,flog)
    call mesh%jacobian(cell,flog,vgop,vol)
    call field%interp(cell,flog,vgop,etmp)
    energy = energy + DOT_PRODUCT(etmp,norm)*area*quad%wts(m)
  end do
end do
!---Global reduction and cleanup
energy=oft_mpi_sum(energy)
CALL quad%delete
DEBUG_STACK_POP
END FUNCTION vec_surf_int
END MODULE diagnostic
