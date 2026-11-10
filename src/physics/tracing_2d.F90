!---------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!---------------------------------------------------------------------------------
!> @file tracing.F90
!
!> 2D tracing implementation for Open FUSION Toolkit (OFT)
!!
!! @authors Chris Hansen
!! @date Feburary 2012
!! @ingroup doxy_oft_physics
!---------------------------------------------------------------------------------
module tracing_2d
USE oft_base
USE oft_mesh_type, ONLY: oft_bmesh, bmesh_findcell
USE oft_la_base, ONLY: oft_vector
USE fem_utils, ONLY: bfem_interp
USE oft_lag_basis, ONLY: oft_blag_eval, oft_blag_geval
implicit none
!---------------------------------------------------------------------------------
!> Abstract tracer class for 2D grids
!---------------------------------------------------------------------------------
type, abstract :: tracer
  integer :: neq = 2 !< Number of ODE equations to advance
  logical :: initialized = .FALSE. !< Has tracer been initialized
  logical :: inv = .FALSE. !< Perform tracing in inverse coordinates?
  integer(i4) :: maxsteps = 1e6 !< Maximum number of timesteps
  integer(i4) :: nsteps = 0 !< Number of timesteps taken
  integer(i4) :: cell = 0 !< Current cell containing tracer
  integer(i4) :: status = 0 !< Status of tracer
  real(r8) :: tol = 1.d-4 !< Tolerance for ODE solver
  real(r8) :: t = 0.d0 !< Current time
  real(r8) :: dt = 0.d0 !< Timestep (for fixed step methods)
  real(r8) :: tmax = 0.d0 !< Maximum time for integration
  real(r8) :: y(2) = 0.d0 !< Position in real coordinates (eg. X,Y)
  real(r8) :: dy(2) = 0.d0 !< Change in position in real coordinates (eg. X,Y)
  real(r8) :: f(3) = 0.d0 !< Logical position in cell
  real(r8) :: gop(3,3) = 0.d0 !< Gradients of the logical coordinates at the last field evaluation
  real(r8) :: raxis = 0.d0 !< Radial location of axis for inverse coordinates
  real(r8) :: zaxis = 0.d0 !< Vertical location of axis for inverse coordinates
  real(r8), pointer :: v(:) => NULL() !< Current ODE solution vector
  real(r8), pointer :: dv(:) => NULL() !< ODE right-hand side at the last evaluation
  class(bfem_interp), pointer :: B => NULL() !< Interpolation object for field evaluation
  class(oft_bmesh), pointer :: mesh => NULL()
contains
  !> Setup tracer and initialize ODE solver
  procedure(tracer_setup), deferred :: setup
  !> Take one step of the ODE solver
  procedure(tracer_step), deferred :: step
  !> Destroy tracer and deallocate internal storage
  procedure(tracer_step), deferred :: delete
end type tracer
!
abstract interface
!---------------------------------------------------------------------------------
!> Interface definition for \ref tracer::setup
!---------------------------------------------------------------------------------
   subroutine tracer_setup(self,y,cell)
      import tracer
      class(tracer), intent(inout) :: self !< Tracer object
      real(8), intent(in) :: y(2) !< Initial position in physical coordinates
      integer(4), optional, intent(in) :: cell !< Starting guess for cell
   end subroutine tracer_setup
!---------------------------------------------------------------------------------
!> Interface definition for \ref tracer::step and \ref tracer::delete
!---------------------------------------------------------------------------------
   subroutine tracer_step(self)
      import tracer
      class(tracer), intent(inout) :: self !< Tracer object
   end subroutine tracer_step
end interface
!---------------------------------------------------------------------------------
!> Tracer implementation using LSODE as the ODE solver
!---------------------------------------------------------------------------------
type, extends(tracer) :: tracer_lsode
  integer :: lrw = 0 !< Size of LSODE floating point work array
  integer :: liw = 0 !< Size of LSODE integer work array
  integer :: itask = 0 !< Task identifier for LSODE
  integer :: istate = 0 !< Internal state for LSODE
  integer, pointer, dimension(:) :: iwork => NULL() !< LSODE integer work array
  real(r8), pointer, dimension(:) :: rwork => NULL() !< LSODE floating point work array
contains
  procedure :: setup => trace_setup_lsode
  procedure :: step => trace_advance_lsode
  procedure :: delete => trace_delete_lsode
end type tracer_lsode
!---------------------------------------------------------------------------------
!> Tracer implementation using an embedded explicit Runge-Kutta method (Bogacki-Shampine 3(2))
!!
!! The field interpolated from a finite element representation is only piecewise smooth,
!! with jumps in the field or its derivatives at cell boundaries. Steps are therefore
!! limited to end just past the (linearly extrapolated) exit point of the current cell, so
!! that each step integrates a smooth right-hand side. This one-step method keeps no
!! history across cell boundaries, so its error estimate remains reliable and the step
!! size is set by the tolerance instead of by repeated crossings.
!---------------------------------------------------------------------------------
type, extends(tracer) :: tracer_rk
  real(r8) :: h = 0.d0 !< Step size for the next step
  real(r8), pointer, dimension(:,:) :: k => NULL() !< Stage derivatives [neq,4]
contains
  procedure :: setup => trace_setup_rk
  procedure :: step => trace_advance_rk
  procedure :: delete => trace_delete_rk
end type tracer_rk
!---------------------------------------------------------------------------------
!> Tracer implementation using a simple forward Euler method as the ODE solver
!---------------------------------------------------------------------------------
type, extends(tracer) :: tracer_euler
  real(r8), pointer :: dv_last(:) => NULL() !< ODE LHS at last step (for error estimation)
contains
  procedure :: setup => trace_setup_euler
  procedure :: step => trace_advance_euler
  procedure :: delete => trace_delete_euler
end type tracer_euler
!---------------------------------------------------------------------------------
!> Abstract interpolation class for inverse mappings
!---------------------------------------------------------------------------------
type, abstract, extends(bfem_interp) :: cylinv_interp
  real(r8) :: rho = 0.d0 !< Radial position
  real(r8) :: t = 0.d0 !< Poloidal angle
  real(r8) :: eps = 1.d-8 !< Epsilon for 1/rho singularity
end type cylinv_interp
class(tracer), pointer :: active_tracer => NULL() !< Current tracer for thread (OpenMP)
!$omp threadprivate(active_tracer)
! integer(4), parameter, private :: nlocpts = 2e4
contains
!------------------------------------------------------------------------------
!> Allocate tracer for current thread
!!
!! Available types:
!! - `1`: Adams method using LSODE
!! - `2`: Forward Euler
!! - `3`: Bogacki-Shampine 3(2) Runge-Kutta method with steps limited at cell boundaries
!------------------------------------------------------------------------------
subroutine set_tracer(type)
integer(i4), intent(in) :: type !< Desired type (see above)
!$omp parallel
if(associated(active_tracer))then
  call active_tracer%delete
  deallocate(active_tracer)
end if
select case(type)
  case(1)
    allocate(tracer_lsode::active_tracer)
    active_tracer%status=0
  case(2)
    allocate(tracer_euler::active_tracer)
    active_tracer%status=0
  case(3)
    allocate(tracer_rk::active_tracer)
    active_tracer%status=0
  case default
    call oft_abort('Invalid tracer type.','set_tracer',__FILE__)
end select
!$omp end parallel
end subroutine set_tracer
!------------------------------------------------------------------------------
!> Trace field line for one flux surface transit in inverse coordinates
!------------------------------------------------------------------------------
subroutine tracinginv_fs(mesh,pt,ptout,vout)
class(oft_bmesh), target, intent(in) :: mesh !< Mesh for tracing
real(8), intent(in) :: pt(2) !< Starting point [2]
real(8), optional, intent(inout) :: ptout(:,:) !< Points on surface [3,:]
real(8), optional, intent(inout) :: vout(:,:) !< Leading ODE RHS components at the last evaluation of each step [:,maxsteps]
real(8) :: z,tp,yp(2)
integer(4) :: ncross
!---
IF(present(ptout))THEN
  ptout(2:3,1)=pt
  ptout(1,1)=0.d0
END IF
!---
ncross=0
active_tracer%mesh=>mesh
call active_tracer%setup(pt)
z=pt(2)
active_tracer%status=0
do while(active_tracer%nsteps<active_tracer%maxsteps)
    yp=active_tracer%y
    tp=active_tracer%t
    call active_tracer%step
    if(active_tracer%status<0)exit
    if(active_tracer%cell==0)then
      active_tracer%status=-2
      exit
    end if
    IF(present(ptout))THEN
      ptout(2:3,active_tracer%nsteps+1)=active_tracer%y
      ptout(1,active_tracer%nsteps+1)=active_tracer%t
    END IF
    IF(present(vout))vout(:,active_tracer%nsteps)=active_tracer%dv(1:SIZE(vout,1))
    if(active_tracer%t>=2.d0*pi)then
        active_tracer%status=1
        exit
    end if
    if(active_tracer%t-tp<active_tracer%tol*1.d-4)then
        active_tracer%status=-5
        exit
    end if
end do
IF(active_tracer%nsteps>=active_tracer%maxsteps)active_tracer%status=-1
end subroutine tracinginv_fs
!------------------------------------------------------------------------------
!> Trace field line and save path to file
!------------------------------------------------------------------------------
subroutine tracing_line(pt,filename)
real(r8), intent(in) :: pt(2) !< Starting point
character(LEN=*), intent(in) :: filename !< Output filename
integer(i4) :: io_unit
open(NEWUNIT=io_unit,FILE=TRIM(filename))
call active_tracer%setup(pt)
do while(active_tracer%nsteps<active_tracer%maxsteps)
    call active_tracer%step
    if(active_tracer%cell==0)exit
    write(io_unit,*)active_tracer%y
end do
close(io_unit)
end subroutine tracing_line
!------------------------------------------------------------------------------
!> Evaluate B-field for tracing using LSODE (called by LSODE)
!------------------------------------------------------------------------------
subroutine tracing_eval_B(neq,t,y,ydot)
integer, intent(in) :: neq !< Number of total ODE eqns (from LSODE)
real(r8), intent(in) :: t !< Current time (from LSODE)
real(r8), intent(in) :: y(neq) !< State vector at current time (from LSODE)
real(r8), intent(out) :: ydot(neq) !< dy/dt vector
real(r8) :: goptmp(3,3),pttmp(3),v,fmin,fmax,tol=1.d-4
active_tracer%y=y(1:2)
active_tracer%v=y
ydot=active_tracer%dv
if(active_tracer%cell==0)return
pttmp=[active_tracer%y(1),active_tracer%y(2),0.d0]
call bmesh_findcell(active_tracer%mesh,active_tracer%cell,pttmp,active_tracer%f)
if(active_tracer%cell==0)return
fmin=minval(active_tracer%f)
fmax=maxval(active_tracer%f)
call active_tracer%mesh%jacobian(active_tracer%cell,active_tracer%f,goptmp,v)
active_tracer%gop=goptmp
call active_tracer%B%interp(active_tracer%cell,active_tracer%f,goptmp,ydot)
active_tracer%dy=ydot(1:2)
active_tracer%dv=ydot
if(.NOT.(( fmax<=1.d0+tol ).AND.( fmin>=-tol )))active_tracer%cell=0
end subroutine tracing_eval_B
!------------------------------------------------------------------------------
!> Evaluate B-field for tracing in inverse coordinates using LSODE (called by LSODE)
!------------------------------------------------------------------------------
subroutine tracing_eval_Binv(neq,t,y,ydot)
integer, intent(in) :: neq !< Number of total ODE eqns (from LSODE)
real(r8), intent(in) :: t !< Current time (from LSODE)
real(r8), intent(in) :: y(neq) !< State vector at current time (from LSODE)
real(r8), intent(out) :: ydot(neq) !< dy/dt vector
real(r8) :: goptmp(3,3),v,fmin,fmax,tol=1.d-4
real(r8) :: rho,s,c,pttmp(3)
!---
rho=y(1)
s=sin(t)
c=cos(t)
!---
active_tracer%y(1)=active_tracer%raxis+rho*c
active_tracer%y(2)=active_tracer%zaxis+rho*s
active_tracer%v=y
ydot=active_tracer%dv
if(active_tracer%cell==0)return
pttmp=[active_tracer%y(1),active_tracer%y(2),0.d0]
call bmesh_findcell(active_tracer%mesh,active_tracer%cell,pttmp,active_tracer%f)
if(active_tracer%cell==0)return
fmin=minval(active_tracer%f)
fmax=maxval(active_tracer%f)
!---
call active_tracer%mesh%jacobian(active_tracer%cell,active_tracer%f,goptmp,v)
active_tracer%gop=goptmp
!---
SELECT TYPE(this=>active_tracer%B)
CLASS IS(cylinv_interp)
  this%rho=rho
  this%t=t
  call this%interp(active_tracer%cell,active_tracer%f,goptmp,ydot)
CLASS DEFAULT
  CALL oft_abort('Tracer interpolator must be of type "cylinv_interp"', &
    "tracing_eval_Binv",__FILE__)
END SELECT
! active_tracer%dy=ydot(1:2) ! \todo Update to give B in physical coordinates
active_tracer%dv=ydot
if(.NOT.(( fmax<=1.d0+tol ).AND.( fmin>=-tol )))active_tracer%cell=0
end subroutine tracing_eval_Binv
!------------------------------------------------------------------------------
!> Initialize tracer using LSODE
!------------------------------------------------------------------------------
subroutine trace_setup_lsode(self,y,cell)
class(tracer_lsode), intent(inout) :: self !< LSODE tracer object
real(8), intent(in) :: y(2) !< Initial position
integer(4), optional, intent(in) :: cell !< Guess for starting cell
real(8) :: pttmp(3)
if(self%neq<=0)call oft_abort('Invalid number of equations.','trace_setup_lsode',__FILE__)
!---
self%istate=1
self%itask=2
!---
self%lrw=20+16*self%neq
self%liw=20
IF(.NOT.self%initialized)THEN
  allocate(self%rwork(self%lrw))
  allocate(self%iwork(self%liw))
  allocate(self%v(self%neq))
  allocate(self%dv(self%neq))
END IF
self%rwork=0.d0
self%iwork=0
self%v=0.d0
self%dv=0.d0
!---
self%t=0.d0
self%tmax=2.d0*pi
self%nsteps=0
!---Find starting point
self%cell=0
if(present(cell))self%cell=cell
IF(self%inv)THEN
  self%v(1)=y(1)-self%raxis
  self%v(2)=0.d0
  self%rwork(1)=2.d0*pi
  self%itask=5
ELSE
  self%v(1:2)=y
END IF
!---
self%y=y
self%dy=0.d0
pttmp=[self%y(1),self%y(2),0.d0]
call bmesh_findcell(self%mesh,self%cell,pttmp,self%f)
self%initialized=.TRUE.
self%status=0
end subroutine trace_setup_lsode
!------------------------------------------------------------------------------
!> Advance the tracer one step using LSODE
!------------------------------------------------------------------------------
subroutine trace_advance_lsode(self)
class(tracer_lsode), intent(inout) :: self !< LSODE tracer object
IF(self%inv)THEN
  call dlsode(tracing_eval_Binv,self%neq,self%v,self%t,self%tmax,1,self%tol, &
    self%tol,self%itask,self%istate,0,self%rwork,self%lrw,self%iwork,self%liw, &
    tracing_eval_Binv,10)
ELSE
  call dlsode(tracing_eval_B,self%neq,self%v,self%t,self%tmax,1,self%tol, &
    self%tol,self%itask,self%istate,0,self%rwork,self%lrw,self%iwork,self%liw, &
    tracing_eval_B,10)
END IF
self%nsteps=self%nsteps+1
IF(self%istate<0)self%status=-4
end subroutine trace_advance_lsode
!------------------------------------------------------------------------------
!> Destroy LSODE tracer and deallocate internal storage
!------------------------------------------------------------------------------
subroutine trace_delete_lsode(self)
class(tracer_lsode), intent(inout) :: self !< LSODE tracer object
self%nsteps=0
self%initialized=.FALSE.
IF(ASSOCIATED(self%rwork))deallocate(self%rwork)
IF(ASSOCIATED(self%iwork))deallocate(self%iwork)
IF(ASSOCIATED(self%v))deallocate(self%v)
IF(ASSOCIATED(self%dv))deallocate(self%dv)
end subroutine trace_delete_lsode
!------------------------------------------------------------------------------
!> Initialize simple forward Euler tracer
!------------------------------------------------------------------------------
subroutine trace_setup_euler(self,y,cell)
class(tracer_euler), intent(inout) :: self !< Tracer object
real(8), intent(in) :: y(2) !< Initial position
integer(4), optional, intent(in) :: cell !< Guess for starting cell
real(8) :: pttmp(3)
if(self%neq<=0)call oft_abort('Invalid number of equations.','trace_setup_euler',__FILE__)
!---
IF(.NOT.self%initialized)THEN
  allocate(self%v(self%neq))
  allocate(self%dv(self%neq))
  allocate(self%dv_last(self%neq))
END IF
self%v=0.d0
self%dv=0.d0
self%dv_last=0.d0
!---
self%t=0.d0
self%tmax=2.d0*pi
self%dt=self%tmax/1.d4
self%nsteps=0
!---Find starting point
self%cell=0
if(present(cell))self%cell=cell
IF(self%inv)THEN
  self%v(1)=y(1)-self%raxis
  self%v(2)=0.d0
ELSE
  self%v(1:2)=y
END IF
!---
self%y=y
self%dy=0.d0
pttmp=[self%y(1),self%y(2),0.d0]
call bmesh_findcell(self%mesh,self%cell,pttmp,self%f)
self%initialized=.TRUE.
self%status=0
end subroutine trace_setup_euler
!------------------------------------------------------------------------------
!> Advance the tracer one step using a forward Euler method with adaptive timestep control
!------------------------------------------------------------------------------
subroutine trace_advance_euler(self)
class(tracer_euler), intent(inout) :: self !< Tracer object
self%dv_last=self%dv
IF(self%inv)THEN
  CALL tracing_eval_Binv(self%neq,self%t,self%v,self%dv)
ELSE
  CALL tracing_eval_B(self%neq,self%t,self%v,self%dv)
END IF
IF(self%t>0.d0)THEN ! Approximate error over last step and update timestep
  self%dv_last=self%dt*(self%dv - self%dv_last)/2.d0 ! LTE = h^2 * y'' / 2 = h^2 * (yd_1 - yd_0)/(2*h)
  self%dt=MIN(self%tol/SQRT(SUM(self%dv_last**2)),2.d0*self%dt)
  IF(self%t+self%dt>self%tmax)self%dt=self%tmax-self%t ! Don't overshoot target time
END IF
self%v=self%v+self%dt*self%dv
self%t=self%t+self%dt
self%nsteps=self%nsteps+1
end subroutine trace_advance_euler
!------------------------------------------------------------------------------
!> Destroy Euler tracer and deallocate internal storage
!------------------------------------------------------------------------------
subroutine trace_delete_euler(self)
class(tracer_euler), intent(inout) :: self !< Tracer object
self%nsteps=0
self%initialized=.FALSE.
IF(ASSOCIATED(self%v))deallocate(self%v)
IF(ASSOCIATED(self%dv))deallocate(self%dv)
IF(ASSOCIATED(self%dv_last))deallocate(self%dv_last)
end subroutine trace_delete_euler
!------------------------------------------------------------------------------
!> Evaluate the tracing ODE right-hand side for the RK tracer
!!
!! Also updates the tracer's current cell, logical position and field.
!------------------------------------------------------------------------------
SUBROUTINE trace_rk_eval(self,t,y,ydot)
CLASS(tracer_rk), INTENT(inout) :: self !< RK tracer object
REAL(r8), INTENT(in) :: t !< Time (or poloidal angle for inverse coordinates)
REAL(r8), INTENT(in) :: y(:) !< State vector [neq]
REAL(r8), INTENT(out) :: ydot(:) !< ODE right-hand side [neq]
REAL(r8) :: ytmp(self%neq)
ytmp=y ! Copy as the evaluation routines overwrite `self%v`
IF(self%inv)THEN
  CALL tracing_eval_Binv(self%neq,t,ytmp,ydot)
ELSE
  CALL tracing_eval_B(self%neq,t,ytmp,ydot)
END IF
END SUBROUTINE trace_rk_eval
!------------------------------------------------------------------------------
!> Initialize tracer using the Bogacki-Shampine 3(2) Runge-Kutta method
!------------------------------------------------------------------------------
SUBROUTINE trace_setup_rk(self,y,cell)
CLASS(tracer_rk), INTENT(inout) :: self !< RK tracer object
REAL(8), INTENT(in) :: y(2) !< Initial position
INTEGER(4), OPTIONAL, INTENT(in) :: cell !< Guess for starting cell
REAL(r8) :: pttmp(3),d0,d1,d2,h0
REAL(r8), ALLOCATABLE, DIMENSION(:) :: v0,dv0,sc
IF(self%neq<=0)CALL oft_abort('Invalid number of equations.','trace_setup_rk',__FILE__)
!---
IF(.NOT.self%initialized)THEN
  ALLOCATE(self%v(self%neq),self%dv(self%neq),self%k(self%neq,4))
END IF
self%v=0.d0
self%dv=0.d0
self%k=0.d0
!---
self%t=0.d0
self%tmax=2.d0*pi
self%nsteps=0
!---Find starting point
self%cell=0
IF(PRESENT(cell))self%cell=cell
IF(self%inv)THEN
  self%v(1)=y(1)-self%raxis
  self%v(2)=0.d0
ELSE
  self%v(1:2)=y
END IF
!---
self%y=y
self%dy=0.d0
pttmp=[self%y(1),self%y(2),0.d0]
CALL bmesh_findcell(self%mesh,self%cell,pttmp,self%f)
self%initialized=.TRUE.
self%status=0
!---Initial step size (Hairer, Norsett & Wanner, Solving ODEs I, Sec. II.4)
ALLOCATE(v0(self%neq),dv0(self%neq),sc(self%neq))
v0=self%v
CALL trace_rk_eval(self,self%t,v0,dv0)
sc=self%tol*(1.d0+ABS(v0))
d0=SQRT(SUM((v0/sc)**2)/self%neq)
d1=SQRT(SUM((dv0/sc)**2)/self%neq)
IF((d0<1.d-5).OR.(d1<1.d-5))THEN
  h0=1.d-6
ELSE
  h0=1.d-2*d0/d1
END IF
h0=MIN(h0,self%tmax)
CALL trace_rk_eval(self,self%t+h0,v0+h0*dv0,self%k(:,2))
d2=SQRT(SUM(((self%k(:,2)-dv0)/sc)**2)/self%neq)/h0
IF(MAX(d1,d2)<=1.d-15)THEN
  self%h=MAX(1.d-6,h0*1.d-3)
ELSE
  self%h=(1.d-2/MAX(d1,d2))**(1.d0/4.d0)
END IF
self%h=MIN(1.d2*h0,self%h,self%tmax)
!---Restore field, cell and state at the starting point
CALL trace_rk_eval(self,self%t,v0,self%dv)
self%v=v0
DEALLOCATE(v0,dv0,sc)
END SUBROUTINE trace_setup_rk
!------------------------------------------------------------------------------
!> Advance the tracer by one accepted step of the Bogacki-Shampine 3(2) method
!!
!! The step is limited to end just past the point where a straight-line extrapolation
!! of the current direction leaves the current cell. A step that is rejected by the
!! error test, or whose stages leave the mesh, is retried with a smaller step size.
!! If no acceptable step can be found the tracer is marked as having left the mesh
!! (`cell=0`, `status=-4`).
!------------------------------------------------------------------------------
SUBROUTINE trace_advance_rk(self)
CLASS(tracer_rk), INTENT(inout) :: self !< RK tracer object
INTEGER(i4) :: i,cell0,nfail
REAL(r8) :: h,err,fac,facmax,texit,dl,dx(2),f0(3),gop0(3,3),y0(2)
REAL(r8) :: v0(self%neq),dv0(self%neq),ynew(self%neq),ytmp(self%neq),sc(self%neq)
LOGICAL :: last,clipped,outside
REAL(r8), PARAMETER :: c2=0.5d0, c3=0.75d0
REAL(r8), PARAMETER :: a32=0.75d0
REAL(r8), PARAMETER :: b1=2.d0/9.d0, b2=1.d0/3.d0, b3=4.d0/9.d0
REAL(r8), PARAMETER :: e1=-5.d0/72.d0, e2=1.d0/12.d0, e3=1.d0/9.d0, e4=-1.d0/8.d0
REAL(r8), PARAMETER :: clip_rel=1.d-6 !< Overshoot past the estimated cell exit, relative to the distance
REAL(r8) :: clip_min
!---Save state at the start of the step
cell0=self%cell; f0=self%f; gop0=self%gop; y0=self%y
v0=self%v; dv0=self%dv
facmax=5.d0
nfail=0
clip_min=MAX(1.d-6,MIN(1.d-3,self%tol)) ! Overshoot past the estimated cell exit, relative to the nominal step
DO
  h=self%h
  last=.FALSE.
  clipped=.FALSE.
  !---End the step just past the exit of the current cell (linear estimate)
  IF(cell0>0)THEN
    IF(self%inv)THEN
      dx(1)=dv0(1)*COS(self%t)-v0(1)*SIN(self%t)
      dx(2)=dv0(1)*SIN(self%t)+v0(1)*COS(self%t)
    ELSE
      dx=dv0(1:2)
    END IF
    texit=HUGE(1.d0)
    DO i=1,3
      dl=gop0(1,i)*dx(1)+gop0(2,i)*dx(2)
      IF(dl<0.d0)texit=MIN(texit,-f0(i)/dl)
    END DO
    texit=MAX(texit,0.d0)*(1.d0+clip_rel)+clip_min*h
    IF(texit<h)THEN
      h=texit
      clipped=.TRUE.
    END IF
  END IF
  IF(self%inv.AND.(self%t+h>=self%tmax))THEN ! Do not pass one full poloidal transit
    h=self%tmax-self%t
    last=.TRUE.
  END IF
  !---Stages (first same as last)
  self%k(:,1)=dv0
  CALL trace_rk_eval(self,self%t+c2*h,v0+(c2*h)*self%k(:,1),self%k(:,2))
  outside=(self%cell==0)
  IF(.NOT.outside)THEN
    CALL trace_rk_eval(self,self%t+c3*h,v0+(a32*h)*self%k(:,2),self%k(:,3))
    outside=(self%cell==0)
  END IF
  IF(.NOT.outside)THEN
    ynew=v0+h*(b1*self%k(:,1)+b2*self%k(:,2)+b3*self%k(:,3))
    CALL trace_rk_eval(self,self%t+h,ynew,self%k(:,4))
    outside=(self%cell==0)
  END IF
  !---Error test
  IF(outside)THEN
    err=1.d10
  ELSE
    ytmp=h*(e1*self%k(:,1)+e2*self%k(:,2)+e3*self%k(:,3)+e4*self%k(:,4))
    sc=self%tol*(1.d0+MAX(ABS(v0),ABS(ynew)))
    err=SQRT(SUM((ytmp/sc)**2)/self%neq)
  END IF
  IF(err<=1.d0)THEN
    self%v=ynew
    self%dv=self%k(:,4)
    IF(last)THEN
      self%t=self%tmax
    ELSE
      self%t=self%t+h
    END IF
    IF(err<1.d-10)THEN
      fac=facmax
    ELSE
      fac=MIN(facmax,MAX(0.2d0,0.9d0*err**(-1.d0/3.d0)))
    END IF
    IF(clipped.OR.last)THEN ! Truncated steps should not reduce the nominal step size
      self%h=MAX(self%h,h*fac)
    ELSE
      self%h=h*fac
    END IF
    EXIT
  END IF
  !---Reject step and restore state at the start of the step
  nfail=nfail+1
  self%h=h*MAX(0.2d0,0.9d0*err**(-1.d0/3.d0))
  facmax=1.d0
  self%cell=cell0; self%f=f0; self%gop=gop0; self%y=y0
  self%v=v0; self%dv=dv0
  IF((self%t+self%h==self%t).OR.(nfail>50))THEN
    self%cell=0
    self%status=-4
    RETURN
  END IF
END DO
self%nsteps=self%nsteps+1
END SUBROUTINE trace_advance_rk
!------------------------------------------------------------------------------
!> Destroy RK tracer and deallocate internal storage
!------------------------------------------------------------------------------
SUBROUTINE trace_delete_rk(self)
CLASS(tracer_rk), INTENT(inout) :: self !< RK tracer object
self%nsteps=0
self%initialized=.FALSE.
IF(ASSOCIATED(self%v))DEALLOCATE(self%v)
IF(ASSOCIATED(self%dv))DEALLOCATE(self%dv)
IF(ASSOCIATED(self%k))DEALLOCATE(self%k)
END SUBROUTINE trace_delete_rk
end module tracing_2d
