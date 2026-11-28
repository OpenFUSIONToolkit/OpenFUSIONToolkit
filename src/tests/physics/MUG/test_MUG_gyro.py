#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: rotation of the polarization of a shear Alfven wave in a triply periodic box (xMHD, H(curl) magnetic
# field) by Braginskii gyroviscosity and by the Hall term in Ohm's law. Each splits the two circular polarizations;
# the velocity after half a period is compared with linear theory, which fixes the sign of each term and the
# electron inertia coefficient. Shared with test_MUG_hall.py, which holds the Hall-only cases
import os
import sys
import numpy as np
import pytest
from scipy.linalg import expm
test_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(test_dir, '..', '..')))
from oft_testing import run_OFT

# Input file template, filled in by the setup function
oft_in_template = """
&runtime_options
 ppn=1
 debug=0
 test_run=T
 use_petsc={7}
/

&mesh_options
 meshname='cube'
 cad_type=92
 nlevels=1
 nbase=1
 grid_order=1
/

&cube_options
 mesh_type=1
 ni=1,2,4
 rscale=1.,1.,2.
 shift=-0.5,-0.5,-1.
 ref_per=T,T,T
/

&test_gyro_options
 order=2
 minlev=2
 delta=1.e-3
 ti0={0}
 linear={6}
/

&xmhd_options
 xmhd_jcb=T
 xmhd_advec=F
 xmhd_hall={4}
 xmhd_adv_temp=F
 xmhd_adv_den=F
 xmhd_diss_centered=T
 bbc='bc'
 vbc='all'
 dt=8.e-7
 eta={1}
 nu_visc={2}
 nu_gyro={3}
 me_factor={5}
 xmhd_lin_euler_start=T
 nsteps=125
 rst_freq=1000
 lin_tol=1.E-12
 nl_tol=1.E-5
 nclean=5000
 rst_ind=0
 maxextrap=2
 ittarget=100
/
"""

# Wave and plasma parameters (must match test_MUG_gyro.F90)
k_wave = np.pi      # Wavenumber (wavelength = 2 m)
v_alf = 1.e4        # Alfven speed
n0 = 1.e19          # Density
mu_ion = 2.0        # Ion mass (amu)
t_final = 1.e-4     # One half wave period (nsteps*dt)
eta = 1.e-3         # Resistivity (diffusivity units)
nu_visc = 1.e-3     # Viscosity (diffusivity units)
ti0 = 2.0           # Ion temperature (eV)
phys_tol = 1.e-2    # Allowed deviation; the gyroviscous and Hall rotations each change V_y by ~0.47
# Physical constants (as in mhd_utils.F90)
mu0 = 4.e-7*np.pi
proton_mass = 1.672621777E-27
elec_mass = 9.10938291E-31
elec_charge = 1.602176565E-19
B0 = v_alf*np.sqrt(mu0*mu_ion*n0*proton_mass)

def alfven_linear(nu_gyro, hall=False, me_factor=100.):
    '''Projections of V_x and V_y onto sin(kz) and cos(kz) at t_final from linear theory

    Linearized momentum and induction equations for a single Fourier mode with
    k || B0 = B0 zhat and normalized perturbations U=V/v_A and b=dB/B0. Braginskii's
    gyroviscous stress gives rho dV_x/dt = -eta_4 k^2 V_y and rho dV_y/dt = eta_4 k^2 V_x,
    with eta_4/rho = T_i/B0 (T_i in eV).

    With `hall` the Hall term J x B/(mu0 e n) in Ohm's law gives db/dt = k^2 v_A d_i (zhat x b), and electron
    inertia multiplies db/dt by (1 + k^2 d_e^2), with d_e^2 = me_factor m_e/(mu0 e^2 n0). The electron pressure
    term has no effect on this wave (uniform density).
    '''
    g = nu_gyro*k_wave**2*ti0/B0
    h = 0.0
    de2 = 0.0
    if hall:
        h = k_wave**2*B0/(mu0*elec_charge*n0)
        de2 = me_factor*elec_mass/(mu0*elec_charge**2*n0)
    ikv = 1j*k_wave*v_alf
    dk2 = nu_visc*k_wave**2
    ek2 = eta*k_wave**2
    # State = [U_x, U_y, b_x, b_y]
    mat = np.array([
        [-dk2, -g, ikv, 0.0],
        [g, -dk2, 0.0, ikv],
        [ikv, 0.0, -ek2, -h],
        [0.0, ikv, h, -ek2]
    ])
    # Electron inertia acts as a mass term on dB/dt
    mass = np.diag([1.0, 1.0, 1.0 + de2*k_wave**2, 1.0 + de2*k_wave**2])
    # Traveling wave initialized as U_x = b_x = sin(kz) = Re[-i exp(ikz)]
    x0 = -1j*np.array([1.0, 0.0, 1.0, 0.0])
    x = expm(np.linalg.solve(mass, mat)*t_final) @ x0
    return np.array([[-x[0].imag, x[0].real], [-x[1].imag, x[1].real]])

def gyro_case(nu_gyro, hall=False, me_factor=100., linear=False, exe='./test_MUG_gyro', petsc=False):
    '''Write the input for an Alfven wave with ion temperature ti0 and run test_MUG_gyro

    `hall` enables the Hall, electron pressure, and electron inertia terms, with electron mass factor `me_factor`.
    `linear` advances the wave with the linearized model (xmhd_lin_run, started with a backward Euler step) about
    the uniform equilibrium instead of the nonlinear model. `exe` selects the driver program (`./test_MUG_gyro_lag`
    for the Lagrange-only module). `petsc` runs with the PETSc linear algebra backend.
    '''
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(ti0, eta, nu_visc, nu_gyro, 'T' if hall else 'F', me_factor,
                                         'T' if linear else 'F', 'T' if petsc else 'F'))
    # Single FE level: block-Jacobi/LU preconditioner inside one GMRES iteration (as in mug_testing.mug_command),
    # also used for PETSc runs (additive Schwarz without overlap)
    return run_OFT("{0} oft.in oft_in_xmhd.xml".format(exe), 1, 600)

def validate_result(nu_gyro, err_exp, hall=False, me_factor=100.):
    '''Compare the final transverse velocity with linear theory

    The deviation must be within `phys_tol` and may exceed its regression value `err_exp` by at most 1%, and no
    compressional velocity may be generated.
    '''
    proj = np.loadtxt('gyro.results')
    proj_exp = alfven_linear(nu_gyro, hall, me_factor)
    err = np.max(np.abs(proj[:2,:] - proj_exp))
    print("Linear theory [[Vx_sin,Vx_cos],[Vy_sin,Vy_cos]] = {0}".format(proj_exp.tolist()))
    print("Simulation    [[Vx_sin,Vx_cos],[Vy_sin,Vy_cos]] = {0}".format(proj[:2,:].tolist()))
    retval = True
    print("Deviation = {0}".format(err))
    if err > phys_tol:
        print("FAILED: Deviation from linear theory too large! ({0})".format(err))
        retval = False
    if err > 1.01*err_exp:
        print("FAILED: Deviation from linear theory increased!")
        print("  Expected = {0}".format(err_exp))
        print("  Actual =   {0}".format(err))
        retval = False
    if abs(proj[2,:]).max() > phys_tol:
        print("FAILED: Compressional velocity generated! ({0})".format(abs(proj[2,:]).max()))
        retval = False
    return retval

#============================================================================
# Each case runs with the nonlinear and the linearized model, which share the linear theory and the deviation
# tolerance but have separate regression values
#============================================================================
# Reference case without gyroviscosity
@pytest.mark.parametrize("linear,petsc_flag,err_exp", ((False, False, 8.439817E-004),
                                                  pytest.param(True, False, 7.539158E-04, marks=pytest.mark.linear),
                                                  (False, True, 8.439817E-004)))
def test_alfven_ref(linear, petsc_flag, err_exp):
    '''Linearly polarized shear Alfven wave without gyroviscosity, which keeps its polarization.'''
    assert gyro_case(-1., linear=linear, petsc=petsc_flag)
    assert validate_result(0.0, err_exp)

# Gyroviscous rotation of the polarization
@pytest.mark.parametrize("linear,petsc_flag,err_exp", ((False, False, 1.848943E-003),
                                                  pytest.param(True, False, 2.152501E-03, marks=pytest.mark.linear),
                                                  (False, True, 1.848943E-003)))
def test_gyro_rotation(linear, petsc_flag, err_exp):
    '''Shear Alfven wave with Braginskii's gyroviscous stress (S. I. Braginskii, Rev. Plasma Phys. 1, 205 (1965)),
    whose polarization rotates from x toward y (B along +z) at half the split between the circular polarizations,
    eta_4 k^2/rho = k^2 T_i/B0. The opposite sign of the stress would rotate it the other way.'''
    assert gyro_case(1., linear=linear, petsc=petsc_flag)
    assert validate_result(1.0, err_exp)

#============================================================================
# Gyroviscosity with the Hall term (Hall-only cases are in test_MUG_hall.py)
@pytest.mark.parametrize("linear,petsc_flag,err_exp", ((False, False, 2.029630E-03),
                                                  pytest.param(True, False, 1.941121E-03, marks=pytest.mark.linear),
                                                  (False, True, 2.029630E-03)))
def test_hall_gyro(linear, petsc_flag, err_exp):
    '''Shear Alfven wave with both gyroviscosity and the Hall term (default electron mass factor of 100). Their
    rotation rates are in the ratio beta_i/2 ~ 0.96 and, with the conventions above, add, which checks the relative
    sign of the two terms.'''
    assert gyro_case(1., hall=True, linear=linear, petsc=petsc_flag)
    assert validate_result(1.0, err_exp, hall=True)
