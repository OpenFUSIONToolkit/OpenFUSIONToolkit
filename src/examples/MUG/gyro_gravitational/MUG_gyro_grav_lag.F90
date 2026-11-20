!------------------------------------------------------------------------------
! Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
!
! SPDX-License-Identifier: LGPL-3.0-only
!------------------------------------------------------------------------------
!!MUG Example: Gravitational instability with the Lagrange-only model    {#doc_mug_gyro_grav_lag_ex}
!!============================
!!
!![TOC]
!!
!! This example runs the gyroviscous gravitational instability case of \ref doc_mug_gyro_grav_ex with the
!! \ref xmhd_lag "Lagrange-only" extended MHD module, where the magnetic field is represented using vector
!! Lagrange elements with divergence cleaning (see \ref doc_mug_xmhd_eqns_lag). The problem setup, input file,
!! and post-processing are the same as for that example, and the growth rate scan can be run with
!!
!! \verbatim
!! python run_scan.py --exe ./MUG_gyro_grav_lag --plot
!! \endverbatim
!!
!! The only differences from \ref doc_mug_gyro_grav_ex are the FE spaces that are constructed and the
!! projection of the initial magnetic field onto the vector Lagrange space. That projection must be converged
!! with a relative tolerance (see \ref doc_mug_gyro_grav_lag_ex_code_ic), or the force imbalance of the discrete
!! equilibrium exceeds the seeded perturbation (`eps`) and hides the gravitational mode for
!! \f$ \beta \gtrsim 10^{-3} \f$. With the projection converged, the initial kinetic energy from the
!! imbalance is about \f$ 3 \times 10^{-24} \f$ (about \f$ 2 \times 10^{-23} \f$ with the default `eps`), and
!! the measured growth rates agree with those of \ref doc_mug_gyro_grav_ex and with Eq. 15 (without the Hall
!! term) to about 1%. A larger perturbation can still be set with `run_scan.py --eps`.
!!
!!\section doc_mug_gyro_grav_lag_ex_code Code walk through
!!
!!\subsection doc_mug_gyro_grav_lag_ex_code_helper Helper module
!!
!! The helper module defines an implementation of \ref fem_utils::fem_interp "fem_interp" that
!! evaluates the initial density and magnetic field.
! START SOURCE
MODULE gyro_grav_lag_helpers
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
END MODULE gyro_grav_lag_helpers
!!\subsection doc_mug_gyro_grav_lag_ex_code_driver Driver program
!!
!! The driver reads the case parameters in normalized units, builds the Lagrange FE spaces, projects the
!! initial fields, sets gravity, and runs the nonlinear model with a small perturbation.
PROGRAM MUG_gyro_grav_lag
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
USE oft_lag_operators, ONLY: lag_setup_interp, oft_lag_getmop, oft_lag_project, &
  oft_lag_vgetmop, oft_lag_vproject
!---Physics
USE mhd_utils, ONLY: mu0, proton_mass, elec_charge
USE xmhd_lag, ONLY: xmhd_run, xmhd_minlev, xmhd_sub_fields, g_accel, xmhd_ML_lagrange, &
  xmhd_ML_vlagrange
USE gyro_grav_lag_helpers, ONLY: grav_interp
IMPLICIT NONE
!!\subsubsection doc_mug_gyro_grav_lag_ex_code_vars Local variables and input options
!!
!! The case is defined in the paper's normalized units (`beta`, `g`, `Ln`, `I0`), with the reference
!! density `n_bar` [m^-3] and field `B_bar` [T] setting the physical scales. The ion mass is the proton
!! mass, so `mu_ion=1` must be set in `xmhd_options`. The domain size is taken from the mesh.
CLASS(oft_solver), POINTER :: minv => NULL()
CLASS(oft_matrix), POINTER :: mop => NULL()
CLASS(oft_vector), POINTER :: u,v,b,vel,den,temp
TYPE(xmhd_sub_fields) :: equil_fields
TYPE(grav_interp), TARGET :: init_field
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
!!\subsubsection doc_mug_gyro_grav_lag_ex_code_setup Mesh and FE setup
CALL oft_init
OPEN(NEWUNIT=io_unit,FILE=oft_env%ifile)
READ(io_unit,gyro_grav_options,IOSTAT=ierr)
CLOSE(io_unit)
IF(ierr/=0)CALL oft_abort('Error reading "gyro_grav_options"','MUG_gyro_grav_lag',__FILE__)
!---Place the mesh ground point at the top wall, away from the periodic boundaries in x
!   (corner vertices on both periodic boundaries have too many periodic copies)
CALL multigrid_construct(mg_mesh,[0.d0,0.d0,1.d0])
!---Lagrange
ALLOCATE(xmhd_ML_lagrange,xmhd_ML_vlagrange)
CALL oft_lag_setup(mg_mesh,order,xmhd_ML_lagrange,ML_vlag_obj=xmhd_ML_vlagrange,minlev=minlev)
CALL lag_setup_interp(xmhd_ML_lagrange)
!!\subsubsection doc_mug_gyro_grav_lag_ex_code_params Physical parameters
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
  'MUG_gyro_grav_lag',__FILE__)
IF(init_field%B0**2<=init_field%bcoef*n_bar*(EXP(zmax/init_field%Ln)-1.d0)) &
  CALL oft_abort('Field vanishes inside the domain (increase I0)','MUG_gyro_grav_lag',__FILE__)
!!\subsubsection doc_mug_gyro_grav_lag_ex_code_ic Initial conditions
!!
!! The magnetic field is projected onto the vector Lagrange space and the density onto the Lagrange space.
!! The temperature is uniform and the plasma starts at rest.
!!
!! The projection of the magnetic field must be converged with a relative tolerance. The field is supported
!! almost entirely by its own pressure (\f$ \beta \sim 10^{-3} \f$, and the gravitational drive is smaller
!! still), so a relative error of \f$ 10^{-7} \f$ in \f$ B \f$ already gives a force imbalance that is much
!! larger than the seeded perturbation. The right-hand side of the vector Lagrange projection scales with the
!! cell volume (\f$ \sim 10^{-6} \f$ here), so the absolute tolerance used for the density
!! (`atol=1.d-12`) stops the solve after a few iterations at that level of error.
NULLIFY(mop)
CALL oft_lag_vgetmop(xmhd_ML_vlagrange%current_level,mop,"none")
CALL create_cg_solver(minv)
minv%A=>mop
minv%its=-3
minv%atol=0.d0
minv%rtol=1.d-14
CALL create_diag_pre(minv%pre)
CALL xmhd_ML_vlagrange%vec_create(u)
CALL xmhd_ML_vlagrange%vec_create(v)
CALL xmhd_ML_vlagrange%vec_create(b)
init_field%field='b'
CALL oft_lag_vproject(xmhd_ML_lagrange%current_level,init_field,v)
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
!!\subsubsection doc_mug_gyro_grav_lag_ex_code_run Run simulation
!!
!! Gravity points along \f$ -\hat{z} \f$, which corresponds to \f$ g_{accel} > 0 \f$ (see \ref xmhd_lag::g_accel "g_accel").
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
END PROGRAM MUG_gyro_grav_lag
!!\section doc_mug_gyro_grav_lag_ex_input Input files
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
