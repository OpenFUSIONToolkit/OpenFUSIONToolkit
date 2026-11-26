#------------------------------------------------------------------------------
# Flexible Unstructured Simulation Infrastructure with Open Numerics (Open FUSION Toolkit)
#
# SPDX-License-Identifier: LGPL-3.0-only
#------------------------------------------------------------------------------
# MUG tests: viscous damping of a sound wave in a triply periodic box (xMHD, H(curl) magnetic field), with a
# uniform field that sets the direction of the parallel stress but is not advanced. The velocity after one period is
# compared with linear theory for the kin, iso, and parallel (Braginskii) viscous stresses
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
 use_petsc={6}
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
 ni=1,1,4
 rscale=1.,1.,2.
 shift=-0.5,-0.5,-1.
 ref_per=T,T,T
/

&test_visc_options
 order=2
 minlev=2
 delta=1.e-3
 b_dir={0}
/

&xmhd_options
 xmhd_jcb=F
 xmhd_advec=F
 xmhd_adv_den=T
 xmhd_adv_temp=T
 xmhd_hall=F
 xmhd_visc_heat=F
 bbc='bc'
 vbc='all'
 dt=1.e-6
 eta=1.
 visc_type='{1}'
 nu_visc={2}
 nu_par={3}
 d_dens={4}
 nsteps=100
 rst_freq=1000
 lin_tol=1.E-12
 nl_tol=1.E-5
 nclean=500
 rst_ind=0
 maxextrap=2
 kappa_par={5}
 kappa_perp={5}
 ittarget=100
/
"""

# Wave and plasma parameters (must match test_MUG_visc.F90)
k_wave = np.pi      # Wavenumber (wavelength = 2 m)
c_sound = 2.e4      # Sound speed (test_MUG_visc.F90 sets T0 so that c_s = v_sound)
t_final = 1.e-4     # One wave period (nsteps*dt)
gamma = 5.0/3.0     # Adiabatic index
d_dens = 1.0        # Density diffusivity
kappa = 1.e-6       # Thermal diffusivity
nu_visc = 500.0     # Viscosity used for damped cases
b_magic = '0.816496580927726,0.,0.577350269189626' # cos^2 = 1/3 to z, so b.W.b = 0
phys_tol = 1.e-2    # Allowed amplitude deviation; a factor-of-2 coefficient error shifts the amplitude by ~0.2
phase_tol = 2.e-2   # Allowed total deviation, dominated by the numerical phase error of the wave

def sound_linear(nu_long):
    '''Projections of V_z onto sin(kz) and cos(kz) at t_final from linear theory

    `nu_long` is the effective kinematic viscosity acting on the longitudinal velocity. Linearized continuity, momentum, and temperature equations for a single Fourier
    mode with normalized perturbations N=n1/n0, U=v1/c_s, and T=T1/T0, including
    longitudinal viscosity, density diffusion, and thermal conduction.
    '''
    ikc = 1j*k_wave*c_sound
    mat = np.array([
        [-d_dens*k_wave**2, -ikc, 0.0],
        [-ikc/gamma, -nu_long*k_wave**2, -ikc/gamma],
        [0.0, -(gamma-1.0)*ikc, -(gamma-1.0)*kappa*k_wave**2]
    ])
    # Traveling wave initialized as sin(kz) = Re[-i exp(ikz)]
    x0 = -1j*np.array([0.6, 0.6, 0.4])
    u_hat = (expm(mat*t_final) @ x0)[1]/0.6
    return np.array([-u_hat.imag, u_hat.real])

def visc_case(visc_type,nu_visc,nu_par,b_dir,exe='./test_MUG_visc',petsc=False):
    '''Write the input for a damped sound wave and run test_MUG_visc

    `b_dir` is the (unnormalized) field direction as a namelist string. `exe` selects the driver program
    (`./test_MUG_visc_lag` for the Lagrange-only module). `petsc` runs with the PETSc linear algebra backend.
    '''
    os.chdir(test_dir)
    with open('oft.in','w+') as fid:
        fid.write(oft_in_template.format(b_dir,visc_type,nu_visc,nu_par,d_dens,kappa,'T' if petsc else 'F'))
    # Single FE level: block-Jacobi/LU preconditioner inside one GMRES iteration (as in mug_testing.mug_command),
    # also used for PETSc runs (additive Schwarz without overlap)
    return run_OFT("{0} oft.in oft_in_xmhd.xml".format(exe), 1, 600)

def validate_result(nu_long,err_exp):
    '''Compare the final velocity with linear theory

    The amplitude must match linear theory within `phys_tol`, the total deviation (including phase error) must be
    within `phase_tol` and may exceed its regression value `err_exp` by at most 1%, and no transverse velocity may
    be generated.
    '''
    proj = np.loadtxt('visc.results')
    proj_exp = sound_linear(nu_long)
    amp_err = abs(proj[2,0] - proj_exp[0])
    err = np.max(np.abs(proj[2,:] - proj_exp))
    trans = np.max(np.abs(proj[:2,:]))
    print("Linear theory [sin,cos] = {0}".format(proj_exp))
    print("Simulation    [sin,cos] = {0}".format(proj[2,:]))
    retval = True
    print("Deviation = {0}".format(err))
    if amp_err > phys_tol:
        print("FAILED: Amplitude deviation from linear theory too large! ({0})".format(amp_err))
        retval = False
    if err > phase_tol:
        print("FAILED: Phase deviation from linear theory too large! ({0})".format(err))
        retval = False
    if err > 1.01*err_exp:
        print("FAILED: Deviation from linear theory increased!")
        print("  Expected = {0}".format(err_exp))
        print("  Actual =   {0}".format(err))
        retval = False
    if trans > phys_tol:
        print("FAILED: Transverse velocity generated! ({0})".format(trans))
        retval = False
    return retval

#============================================================================
# Isotropic viscosity
@pytest.mark.parametrize("petsc_flag", (False, True))
def test_kin(petsc_flag):
    '''Sound wave damped by the `kin` viscous stress (nu grad V), which acts on the longitudinal velocity with
    coefficient nu.'''
    assert visc_case('kin',nu_visc,-1.,'0.,0.,1.',petsc=petsc_flag)
    assert validate_result(nu_visc,6.246040061832365E-003)

@pytest.mark.parametrize("petsc_flag", (False, True))
def test_iso(petsc_flag):
    '''Sound wave damped by the `iso` viscous stress (nu W, with W the rate of strain tensor), which acts on the
    longitudinal velocity with coefficient 4/3 nu.'''
    assert visc_case('iso',nu_visc,-1.,'0.,0.,1.',petsc=petsc_flag)
    assert validate_result(4.0*nu_visc/3.0,7.595011428105933E-003)

#============================================================================
# Parallel viscosity
@pytest.mark.parametrize("petsc_flag", (False, True))
def test_par_aligned(petsc_flag):
    '''Sound wave along B damped by Braginskii's parallel stress (S. I. Braginskii, Rev. Plasma Phys. 1, 205
    (1965)) with eta_0 = m_i n0 nu_par, which acts on the longitudinal velocity with coefficient 4/3 nu_par, the same
    as `iso`.'''
    assert visc_case('kin',0.,nu_visc,'0.,0.,1.',petsc=petsc_flag)
    assert validate_result(4.0*nu_visc/3.0,7.728128494748115E-003)

@pytest.mark.parametrize("petsc_flag", (False, True))
def test_par_perp(petsc_flag):
    '''Sound wave across B damped by the parallel stress, which acts on the longitudinal velocity with coefficient
    1/3 nu_par, a quarter of the field-aligned case.'''
    assert visc_case('kin',0.,nu_visc,'1.,0.,0.',petsc=petsc_flag)
    assert validate_result(nu_visc/3.0,3.337417275579408E-003)

@pytest.mark.parametrize("petsc_flag", (False, True))
def test_par_magic(petsc_flag):
    '''Sound wave at the angle to B where b.W.b = 0 (cos^2 = 1/3), so the parallel stress vanishes and the wave
    is undamped.'''
    assert visc_case('kin',0.,nu_visc,b_magic,petsc=petsc_flag)
    assert validate_result(0.0,1.452759324069408E-003)
