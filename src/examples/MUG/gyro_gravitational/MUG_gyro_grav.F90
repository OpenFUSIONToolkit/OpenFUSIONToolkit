!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!!MUG Example: Gyroviscous stabilization of the gravitational instability    {#doc_mug_gyro_grav_ex}
!!============================
!!
!![TOC]
!!
!! This example uses the \ref xmhd "extended MHD" module to reproduce the low-\f$ \beta \f$ case of
!! [Ferraro and Jardin, Phys. Plasmas 13, 092101 (2006)](https://doi.org/10.1063/1.2236277), where
!! Braginskii's ion gyroviscosity stabilizes the gravitational (Rayleigh-Taylor) instability of a plasma
!! supported against gravity by a magnetic field gradient. Equation numbers refer to that paper.
!!
!! The growth rate is measured for a range of \f$ \beta \f$ and compared with the local dispersion
!! relation (Eq. 15) and the result of Roberts and Taylor (Eq. 11). This exercises the full gyroviscous
!! stress in a non-uniform plasma, together with gravity and density and field gradients.
!!
!!\section doc_mug_gyro_grav_ex_setup Problem setup
!!
!! The paper works in units of a reference density \f$ \bar{n} \f$ and field \f$ \bar{B} \f$, with time
!! in units of \f$ \Omega_i^{-1} = m_i / e \bar{B} \f$ and length in units of \f$ d_i = c / \omega_{pi} \f$.
!! Its geometry is rotated here, since gravity in MUG acts along \f$ -\hat{z} \f$ (see \ref xmhd::g_accel "g_accel"):
!!  - Gravity \f$ -g \hat{z} \f$ and density \f$ n = n_0 e^{z/L_n} \f$, so the heavier plasma is on top
!!  - Magnetic field \f$ B(z) \hat{y} \f$ along the ignorable direction
!!  - Perturbation \f$ \delta n = \epsilon n_0 \cos(k_x x) \cos(\pi z / L_z) \f$
!!
!! The domain is a thin slab centered on the origin, periodic in \f$ x \f$ and \f$ y \f$ (one cell along the
!! field) and bounded by conducting, no-slip walls at \f$ z = \pm L_z / 2 \f$.
!!
!! The equilibrium is isothermal with \f$ p = p_i = n k T \f$ (`te_factor=0`), and all dissipation
!! (resistivity, viscosity, particle diffusion, and thermal conduction) is turned off. The field follows from force
!! balance, \f$ \frac{d}{dz} \left( n k T + B^2 / 2 \mu_0 \right) = -m_i n g \f$, which gives
!! \f[ B^2 = B_0^2 - 2 \mu_0 \left( m_i g L_n + k T \right) (n - n_0) \f]
!! This replaces \f$ \Gamma T \f$ in the paper's Eq. 13c with \f$ T \f$, which is what force balance
!! requires for \f$ p = n k T \f$.
!!
!! \warning The Hall term is off by default (`xmhd_hall=F`), and the theory curves then drop the Hall
!! term (\f$ 1 / \Omega \f$) from Eqs. 11 and 15. With the Hall term enabled (`run_scan.py --hall`),
!! the gravitational mode initially grows at the expected rate, but a faster mode
!! (\f$ \gamma \approx 0.045 \Omega_i \f$) then takes over. It also appears without gravity and does
!! not change with the time step or element order, so it is not the gravitational instability and is
!! still under investigation.
!!
!! \warning Even without the Hall term, a second mode grows at \f$ \gamma \approx 0.007 \Omega_i \f$ for all
!! \f$ \beta \f$ at the default resolution. It grows from round-off: with no perturbation (`eps=0`) and the
!! initial field projection converged (kinetic energy after the first step about \f$ 2 \times 10^{-26} \f$),
!! it reaches the same amplitude as with the default `eps`, so tightening the initial force balance does
!! not remove it. It overtakes the gravitational mode late in runs with \f$ \beta \gtrsim 2.5 \times 10^{-3} \f$,
!! where that mode is slow or stable (at \f$ \beta = 3 \times 10^{-3} \f$, stable by Eq. 15, the fit returns this
!! mode). It needs both gyroviscosity and a non-uniform equilibrium, but not gravity, and it grows faster at
!! higher resolution (\f$ \gamma \approx 0.019 \Omega_i \f$ with twice the cells in \f$ x \f$ and \f$ z \f$),
!! so it is numerical and still under investigation. For \f$ \beta \le 2.5 \times 10^{-3} \f$ the measured
!! growth rate agrees with Eq. 15 to within about 2%. Runs may still fail near the end (at
!! \f$ \beta = 10^{-3} \f$ the time step collapses once the kinetic energy reaches about \f$ 10^{-9} \f$).
!! To avoid this late stage, and the second mode at higher \f$ \beta \f$, `run_scan.py` fits \f$ \ln E_k \f$ only over the exponential growth
!! phase, the longest part of the run where its slope is constant to within 10% (`--fit-tol`), and prints
!! that window. Still check the kinetic energy history in `xmhd.hist` before trusting a measured growth rate.
!!
!!\section doc_mug_gyro_grav_ex_code Code walk through
!!
!!\subsection doc_mug_gyro_grav_ex_code_helper Helper module
!!
!! The helper module defines an implementation of \ref fem_utils::fem_interp "fem_interp" that
!! evaluates the initial density and magnetic field.
! START SOURCE
MODULE gyro_grav_helpers
USE oft_base
USE fem_utils, ONLY: fem_interp
IMPLICIT NONE
!------------------------------------------------------------------------------
!> Initial condition for the gravitational instability (SI units)
!------------------------------------------------------------------------------
TYPE, EXTENDS(fem_interp) :: grav_interp
  CHARACTER(LEN=1) :: field = 'n' !< Field to evaluate ('n' or 'b')
  REAL(r8) :: n0 = 1.d19 !< Density at the midplane [m^-3]
  REAL(r8) :: B0 = 1.d0 !< Field at the midplane [T]
  REAL(r8) :: Ln = 1.d0 !< Density gradient scale length [m]
  REAL(r8) :: bcoef = 0.d0 !< \f$ 2 \mu_0 (m_i g L_n + k T) \f$
  REAL(r8) :: eps = 1.d-6 !< Relative amplitude of the density perturbation
  REAL(r8) :: kx = 1.d0 !< Wavenumber along x [m^-1]
  REAL(r8) :: Lz = 1.d0 !< Distance between the walls [m]
CONTAINS
  !> Evaluate field
  PROCEDURE :: interp => grav_interp_apply
END TYPE grav_interp
CONTAINS
!------------------------------------------------------------------------------
!> Evaluate the density (equilibrium plus perturbation) or the equilibrium field
!------------------------------------------------------------------------------
SUBROUTINE grav_interp_apply(self,cell,f,gop,val)
CLASS(grav_interp), INTENT(inout) :: self !< Interpolation object
INTEGER(i4), INTENT(in) :: cell !< Cell for interpolation
REAL(r8), INTENT(in) :: f(:) !< Position in cell in logical coord [4]
REAL(r8), INTENT(in) :: gop(3,4) !< Logical gradient vectors at f [3,4]
REAL(r8), INTENT(out) :: val(:) !< Field value at the point
REAL(r8) :: pt(3),neq
pt=self%mesh%log2phys(cell,f)
neq=self%n0*EXP(pt(3)/self%Ln)
SELECT CASE(self%field)
  CASE('n')
    val(1)=neq + self%eps*self%n0*COS(self%kx*pt(1))*COS(pi*pt(3)/self%Lz)
  CASE('b')
    val=[0.d0, SQRT(self%B0**2 - self%bcoef*(neq-self%n0)), 0.d0]
  CASE DEFAULT
    CALL oft_abort('Unknown field','grav_interp_apply',__FILE__)
END SELECT
END SUBROUTINE grav_interp_apply
END MODULE gyro_grav_helpers
!!\subsection doc_mug_gyro_grav_ex_code_driver Driver program
!!
!! The driver reads the case parameters in normalized units, builds the FE spaces, projects the
!! initial fields, sets gravity, and runs the nonlinear model with a small perturbation.
!!
!! \note The linearized model (\ref xmhd::xmhd_lin_run "xmhd_lin_run") is not used: it omits the
!! density dependence of \f$ J \times B / n \f$, which provides the buoyancy drive when the
!! equilibrium is supported by a magnetic field gradient.
PROGRAM MUG_gyro_grav
USE oft_base
!---Grid
USE multigrid, ONLY: multigrid_mesh
USE multigrid_build, ONLY: multigrid_construct
!---Linear algebra
USE oft_la_base, ONLY: oft_vector, oft_matrix
USE oft_solver_base, ONLY: oft_solver
USE oft_solver_utils, ONLY: create_cg_solver, create_diag_pre
!---Lagrange FE space
USE oft_lag_basis, ONLY: oft_lag_setup
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_getmop, oft_lag_project
!---H1 FE space (Grad(H^1) subspace)
USE oft_h1_basis, ONLY: oft_h1_setup
USE oft_h1_operators, ONLY: h1_setup_interp
!---Full H(Curl) FE space
USE oft_hcurl_basis, ONLY: oft_hcurl_setup, oft_hcurl_grad_setup
USE oft_hcurl_operators, ONLY: hcurl_setup_interp
USE oft_hcurl_grad_operators, ONLY: hcurl_grad_setup_interp, hcurl_grad_getmop, &
  oft_hcurl_grad_project, oft_hcurl_grad_gzerop
!---Physics
USE mhd_utils, ONLY: mu0, proton_mass, elec_charge
USE xmhd, ONLY: xmhd_run, xmhd_minlev, xmhd_sub_fields, g_accel, xmhd_ML_hcurl, &
  xmhd_ML_H1, xmhd_ML_hcurl_grad, xmhd_ML_H1grad, xmhd_ML_lagrange, xmhd_ML_vlagrange
USE gyro_grav_helpers, ONLY: grav_interp
IMPLICIT NONE
!!\subsubsection doc_mug_gyro_grav_ex_code_vars Local variables and input options
!!
!! The case is defined in the paper's normalized units (`beta`, `g`, `Ln`, `I0`), with the reference
!! density `n_bar` [m^-3] and field `B_bar` [T] setting the physical scales. The ion mass is the proton
!! mass, so `mu_ion=1` must be set in `xmhd_options`. The domain size is taken from the mesh.
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v,b,vel,den,temp
TYPE(xmhd_sub_fields) :: equil_fields
TYPE(grav_interp), TARGET :: init_field
TYPE(oft_hcurl_grad_gzerop), TARGET :: hcurl_grad_gzerop
TYPE(multigrid_mesh) :: mg_mesh
INTEGER(i4) :: io_unit,ierr
REAL(r8) :: d_i,omega_i,g_phys,T0,Lx,zmin,zmax
INTEGER(i4) :: order = 2 !< FE order
INTEGER(i4) :: minlev = 1 !< Lowest FE level for multigrid
REAL(r8) :: n_bar = 1.d19 !< Reference density [m^-3]
REAL(r8) :: B_bar = 1.d-2 !< Reference field [T]
REAL(r8) :: beta = 1.d-3 !< Plasma beta at the midplane
REAL(r8) :: g = 5.d-3 !< Gravity [d_i Omega_i^2]
REAL(r8) :: Ln = 1.d2 !< Density gradient scale length [d_i]
REAL(r8) :: I0 = 1.d2 !< Field at the midplane [B_bar]
REAL(r8) :: eps = 1.d-6 !< Relative amplitude of the density perturbation
NAMELIST/gyro_grav_options/order,minlev,n_bar,B_bar,beta,g,Ln,I0,eps
!!\subsubsection doc_mug_gyro_grav_ex_code_setup Mesh and FE setup
CALL oft_init
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,gyro_grav_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr/=0)CALL oft_abort('Error reading "gyro_grav_options"','MUG_gyro_grav',__FILE__)
!---Ground the H(Curl) gauge at the top wall, away from the periodic boundaries in x
!   (corner vertices on both periodic boundaries have too many periodic copies)
CALL multigrid_construct(mg_mesh,[0.d0,0.d0,1.d0])
!---Lagrange
ALLOCATE(xmhd_ML_lagrange,xmhd_ML_vlagrange)
CALL oft_lag_setup(mg_mesh,order,xmhd_ML_lagrange,ML_vlag_obj=xmhd_ML_vlagrange,minlev=minlev)
CALL lag_setup_interp(xmhd_ML_lagrange)
!---Grad(H^1) subspace
ALLOCATE(xmhd_ML_H1)
CALL oft_h1_setup(mg_mesh,order+1,xmhd_ML_H1,minlev=minlev)
CALL h1_setup_interp(xmhd_ML_H1)
!---H(Curl) subspace
ALLOCATE(xmhd_ML_hcurl)
CALL oft_hcurl_setup(mg_mesh,order,xmhd_ML_hcurl,minlev=minlev)
CALL hcurl_setup_interp(xmhd_ML_hcurl)
!---Full H(Curl) space
ALLOCATE(xmhd_ML_hcurl_grad,xmhd_ML_H1grad)
CALL oft_hcurl_grad_setup(xmhd_ML_hcurl,xmhd_ML_H1,xmhd_ML_hcurl_grad,xmhd_ML_H1grad,minlev)
CALL hcurl_grad_setup_interp(xmhd_ML_hcurl_grad,xmhd_ML_H1)
hcurl_grad_gzerop%ML_hcurl_grad_rep=>xmhd_ML_hcurl_grad
!!\subsubsection doc_mug_gyro_grav_ex_code_params Physical parameters
!!
!! Normalized parameters are converted to SI units, and the equilibrium temperature is set from
!! \f$ \beta = 2 \mu_0 n_0 k T / B_0^2 \f$ at the midplane.
omega_i=elec_charge*B_bar/proton_mass
d_i=SQRT(proton_mass/(mu0*n_bar*elec_charge**2))
g_phys=g*d_i*omega_i**2
Lx=MAXVAL(mg_mesh%mesh%r(1,:))-MINVAL(mg_mesh%mesh%r(1,:))
zmin=-oft_mpi_max(-MINVAL(mg_mesh%mesh%r(3,:)))
zmax=oft_mpi_max(MAXVAL(mg_mesh%mesh%r(3,:)))
Lx=oft_mpi_max(Lx)
init_field%mesh=>mg_mesh%mesh
init_field%n0=n_bar
init_field%B0=I0*B_bar
init_field%Ln=Ln*d_i
init_field%eps=eps
init_field%kx=2.d0*pi/Lx
init_field%Lz=zmax-zmin
T0=beta*init_field%B0**2/(2.d0*mu0*n_bar*elec_charge)
init_field%bcoef=2.d0*mu0*(proton_mass*g_phys*init_field%Ln + elec_charge*T0)
IF(oft_env%head_proc)THEN
  WRITE(*,*)
  WRITE(*,'(A)')'Gravitational instability case'
  WRITE(*,'(A,ES12.4)')'  d_i [m]        = ',d_i
  WRITE(*,'(A,ES12.4)')'  Omega_i [s^-1] = ',omega_i
  WRITE(*,'(A,ES12.4)')'  T [eV]         = ',T0
  WRITE(*,'(A,ES12.4)')'  g [m/s^2]      = ',g_phys
  WRITE(*,'(A,2ES12.4)')'  Lx, Lz [d_i]   = ',Lx/d_i,init_field%Lz/d_i
  WRITE(*,*)
END IF
IF(ABS(zmax+zmin)>1.d-6*init_field%Lz)CALL oft_abort('Mesh must be centered in z', &
  'MUG_gyro_grav',__FILE__)
IF(init_field%B0**2<=init_field%bcoef*n_bar*(EXP(zmax/init_field%Ln)-1.d0)) &
  CALL oft_abort('Field vanishes inside the domain (increase I0)','MUG_gyro_grav',__FILE__)
!!\subsubsection doc_mug_gyro_grav_ex_code_ic Initial conditions
!!
!! The magnetic field is projected onto the H(Curl) space and the density onto the Lagrange space.
!! The temperature is uniform and the plasma starts at rest.
!!
!! The projection of the magnetic field must be converged with a relative tolerance. The field is supported
!! almost entirely by its own pressure (\f$ \beta \sim 10^{-3} \f$, and the gravitational drive is smaller
!! still), so a small relative error in \f$ B \f$ gives a force imbalance much larger than the seeded
!! perturbation. The right-hand side of the projection is small (it scales with the cell size), so the
!! absolute tolerance used for the density (`atol=1.d-12`) would stop the solve well before that level.
NULLIFY(mop)
CALL hcurl_grad_getmop(xmhd_ML_hcurl_grad%current_level,mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=0.d0
minv%rtol=1.d-14
CALL create_diag_pre(minv%pre)
CALL xmhd_ML_hcurl_grad%vec_create(u)
CALL xmhd_ML_hcurl_grad%vec_create(v)
CALL xmhd_ML_hcurl_grad%vec_create(b)
init_field%field='b'
CALL oft_hcurl_grad_project(xmhd_ML_hcurl_grad%current_level,init_field,v)
CALL hcurl_grad_gzerop%apply(v)
CALL b%set(0.d0)
CALL minv%apply(b,v)
CALL u%delete
CALL v%delete
CALL mop%delete
DEALLOCATE(u,v,mop)
CALL minv%pre%delete
CALL minv%delete
!---Density
NULLIFY(mop)
CALL oft_lag_getmop(xmhd_ML_lagrange%current_level,mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=1.d-12
CALL create_diag_pre(minv%pre)
CALL xmhd_ML_lagrange%vec_create(v)
CALL xmhd_ML_lagrange%vec_create(den)
init_field%field='n'
CALL oft_lag_project(xmhd_ML_lagrange%current_level,init_field,v)
CALL den%set(0.d0)
CALL minv%apply(den,v)
CALL v%delete
CALL mop%delete
DEALLOCATE(v,mop)
CALL minv%pre%delete
CALL minv%delete
!---Temperature and velocity
CALL xmhd_ML_lagrange%vec_create(temp)
CALL temp%set(T0)
CALL xmhd_ML_vlagrange%vec_create(vel)
CALL vel%set(0.d0)
!!\subsubsection doc_mug_gyro_grav_ex_code_run Run simulation
!!
!! Gravity points along \f$ -\hat{z} \f$, which corresponds to \f$ g_{accel} > 0 \f$ (see \ref xmhd::g_accel "g_accel").
!! The growth rate is then obtained from the kinetic energy in the history file `xmhd.hist`,
!! \f$ \gamma = \frac{1}{2} \partial_t \ln E_k \f$ (see `run_scan.py`).
g_accel=g_phys
xmhd_minlev=minlev
equil_fields%B=>b
equil_fields%V=>vel
equil_fields%Ne=>den
equil_fields%Ti=>temp
CALL xmhd_run(equil_fields)
CALL oft_finalize
END PROGRAM MUG_gyro_grav
!!\section doc_mug_gyro_grav_ex_input Input files
!!
!! The case is run, and the growth rates compared with theory, by `run_scan.py`, which writes an
!! input file of the following form for each value of \f$ \beta \f$.
!!
! OFT_DOC_INCLUDE: oft.in
!!
!!
!! The preconditioner is set in `oft_in.xml` (block-Jacobi with LU on each block).
!!
! OFT_DOC_INCLUDE: oft_in.xml
